// web-fastq-human-remover — remove reads matching a reference genome, in the browser.
// Copyright (C) 2026 Guillaume Gautreau — MaIAGE (UR 1404), INRAE
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed WITHOUT ANY WARRANTY; see the GNU General Public
// License for more details: https://www.gnu.org/licenses/
//
// Method derived from Cleanifier (MIT) — see NOTICE for attribution.

//! fastq-clean-local — read decontamination core, compiled to wasm64.
//!
//! Cuckoo filter with two candidate windows over canonical k-mers. Everything
//! lives in the WebAssembly linear memory; the JS caller decides the layout:
//!   [0 .. IO_BASE)        Rust stack and static data
//!   [IO_BASE .. IDX_BASE) input/output buffers
//!   [IDX_BASE .. end)     the index
//!
//! FASTQ/FASTA parsing happens here rather than in JS, so a whole block crosses
//! the boundary once instead of once per read.

#![no_std]
#![feature(simd_wasm64)]

use core::panic::PanicInfo;

#[panic_handler]
fn panic(_: &PanicInfo) -> ! {
    core::arch::wasm64::unreachable()
}

// ————————————————————————————— global state —————————————————————————————

static mut IDX_BASE: u64 = 0;
static mut IDX_MASK: u64 = 0; // nbuckets - 1, nbuckets being a power of two
static mut K: u32 = 25; // weight: number of significant positions
static mut W: u32 = 25; // span: total width of the seed
static mut CONTIGUOUS: bool = true;
static mut SEED_POS: [u8; 64] = [0; 64]; // significant positions, ascending
static mut SEED_MASK_INT: u64 = 0;       // significant positions as a bit mask
/// Coverage bitmap for the current read, 1024 bases max.
static mut COVER: [u64; 16] = [0; 16];
static mut RNG: u64 = 0x2545_f491_4f6c_dd1d;

/// 0 reads seen · 1 reads dropped · 2 bases seen · 3 k-mers queried
/// 4 k-mers found · 5 insertions · 6 failed insertions · 7 occupied slots
static mut STATS: [u64; 8] = [0; 8];

const MAX_KICKS: u32 = 500;
const EMPTY: u64 = 0;

// ————————————————————————————— primitives —————————————————————————————

#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^= x >> 31;
    x
}

/// 16-bit fingerprint, never zero (zero encodes an empty slot).
#[inline(always)]
fn fingerprint(h: u64) -> u64 {
    let f = (h >> 48) & 0xffff;
    if f == 0 { 1 } else { f }
}

#[inline(always)]
unsafe fn bucket_ptr(b: u64) -> *mut u64 {
    (IDX_BASE + b * 8) as *mut u64
}

#[inline(always)]
fn slot(bucket: u64, i: u32) -> u64 {
    (bucket >> (i * 16)) & 0xffff
}

#[inline(always)]
fn set_slot(bucket: u64, i: u32, v: u64) -> u64 {
    (bucket & !(0xffffu64 << (i * 16))) | (v << (i * 16))
}

/// Branchless search for a fingerprint across the four slots of a bucket.
#[inline(always)]
fn has_fp(bucket: u64, fp: u64) -> bool {
    let x = bucket ^ (fp * 0x0001_0001_0001_0001);
    (x.wrapping_sub(0x0001_0001_0001_0001) & !x & 0x8000_8000_8000_8000) != 0
}

#[inline(always)]
fn next_rand() -> u64 {
    unsafe {
        RNG ^= RNG << 13;
        RNG ^= RNG >> 7;
        RNG ^= RNG << 17;
        RNG
    }
}

// ——————————————————————————— k-mer encoding ———————————————————————————

/// Cleanifier's fast encoding (`dnaencode_fast.py`): `(ascii >> 1) & 7`, which
/// yields **A=0, C=1, T=2, G=3** — note T and G, not the usual A/C/G/T order.
/// N lands on 7 and anything else on a value > 3, so one test rejects them all.
///
/// Getting this wrong is silent: codes stay in range and the index simply never
/// matches. It cost a full debugging round — every seed of a human read came
/// back absent from the human index.
#[inline(always)]
fn code(b: u8) -> u8 {
    let c = (b >> 1) & 7;
    if c > 3 { 255 } else { c }
}

/// Gather the significant positions of a spaced seed out of the sliding window
/// and pack them into a contiguous 2-bit code.
///
/// `window` holds the last W bases, most recent one in the low bits, so the
/// base at seed position `p` sits at bit offset `2*(W-1-p)`. W can reach 33,
/// i.e. 66 bits, hence the u128.
#[inline(always)]
unsafe fn gather(window: u128) -> u64 {
    if CONTIGUOUS {
        return (window as u64) & kmask();
    }
    let mut out: u64 = 0;
    let mut i = 0usize;
    while i < K as usize {
        let p = SEED_POS[i] as u32;
        out = (out << 2) | (((window >> (2 * (W - 1 - p))) & 3) as u64);
        i += 1;
    }
    out
}

#[inline(always)]
unsafe fn kmask() -> u64 {
    if K >= 32 { u64::MAX } else { (1u64 << (2 * K)) - 1 }
}

/// Reverse complement of a packed k-mer, by bit twiddling rather than by
/// gathering a second window.
///
/// With the A=0/C=1/T=2/G=3 encoding, complementing is simply XOR 2 per base
/// (A<->T is 0<->2, C<->G is 1<->3), i.e. XOR with the 0b10 pattern — not the
/// `3 - b` of the classic A/C/G/T order.
///
/// The shortcut is only valid because Cleanifier requires **symmetric** seeds
/// (`mask.py: check_mask`): reverse-complementing the gathered code equals
/// gathering from the reverse-complemented window.
#[inline(always)]
unsafe fn revcomp(mut x: u64) -> u64 {
    x = ((x & 0x3333_3333_3333_3333) << 2) | ((x >> 2) & 0x3333_3333_3333_3333);
    x = ((x & 0x0f0f_0f0f_0f0f_0f0f) << 4) | ((x >> 4) & 0x0f0f_0f0f_0f0f_0f0f);
    x = ((x & 0x00ff_00ff_00ff_00ff) << 8) | ((x >> 8) & 0x00ff_00ff_00ff_00ff);
    x = ((x & 0x0000_ffff_0000_ffff) << 16) | ((x >> 16) & 0x0000_ffff_0000_ffff);
    x = (x << 32) | (x >> 32);
    x >>= 64 - 2 * K;
    x ^ comp_mask()
}

/// 0b1010…10 over 2K bits: complements every base at once.
#[inline(always)]
unsafe fn comp_mask() -> u64 {
    let full = 0xAAAA_AAAA_AAAA_AAAAu64;
    if K >= 32 { full } else { full & ((1u64 << (2 * K)) - 1) }
}

/// Canonical form. Cleanifier uses `rcmode="max"`, so we do too — the choice is
/// arbitrary as long as indexing and lookup agree.
#[inline(always)]
unsafe fn canonical(window: u128) -> u64 {
    let fw = gather(window);
    let rc = revcomp(fw);
    if fw >= rc { fw } else { rc }
}

// ————————————————————————————— index API —————————————————————————————

/// Prepare the index. `log2_buckets` sets its size: nbuckets = 2^log2_buckets,
/// each bucket being 8 bytes holding four fingerprints.
///
/// The power of two is required: the alternate location is derived by XOR,
/// which is only an involution over a full mask. With an arbitrary modulus,
/// cuckoo eviction could not find its way back.
/// `seed_bits` describes the seed: bit p set means position p is significant.
/// `w` is the total span. A contiguous k-mer is simply seed_bits = 2^k - 1, w = k.
/// The weight k is derived by popcount, so it is never inconsistent with the seed.
#[no_mangle]
pub unsafe extern "C" fn idx_init(base: u64, log2_buckets: u32, seed_bits: u64, w: u32) -> u64 {
    NATIVE = false;
    IDX_BASE = base;
    IDX_MASK = (1u64 << log2_buckets) - 1;
    W = w;
    K = seed_bits.count_ones();
    CONTIGUOUS = K == W;
    let mut i = 0u32;
    let mut n_pos = 0usize;
    while i < w {
        if (seed_bits >> i) & 1 == 1 {
            SEED_POS[n_pos] = i as u8;
            n_pos += 1;
        }
        i += 1;
    }
    SEED_MASK_INT = seed_bits;
    RNG = 0x2545_f491_4f6c_dd1d ^ (K as u64);
    let n = IDX_MASK + 1;
    let mut i = 0u64;
    while i < n {
        bucket_ptr(i).write(EMPTY);
        i += 1;
    }
    let mut s = 0;
    while s < 8 {
        STATS[s] = 0;
        s += 1;
    }
    n
}

#[inline(always)]
unsafe fn positions(key: u64) -> (u64, u64, u64) {
    let h = mix(key);
    let fp = fingerprint(h);
    let b1 = h & IDX_MASK;
    let b2 = b1 ^ (mix(fp) & IDX_MASK);
    (b1, b2, fp)
}

#[inline(always)]
unsafe fn insert_key(key: u64) -> bool {
    let (b1, b2, fp) = positions(key);
    // already there?
    if has_fp(bucket_ptr(b1).read(), fp) || has_fp(bucket_ptr(b2).read(), fp) {
        return true;
    }
    let mut cur_fp = fp;
    // free slot in either candidate bucket
    for &cand in [b1, b2].iter() {
        let v = bucket_ptr(cand).read();
        let mut i = 0u32;
        while i < 4 {
            if slot(v, i) == EMPTY {
                bucket_ptr(cand).write(set_slot(v, i, cur_fp));
                STATS[5] += 1;
                STATS[7] += 1;
                return true;
            }
            i += 1;
        }
    }
    // cascading eviction
    let mut kicks = 0u32;
    let mut b = if (next_rand() & 1) == 0 { b1 } else { b2 };
    while kicks < MAX_KICKS {
        let i = (next_rand() & 3) as u32;
        let v = bucket_ptr(b).read();
        let victim = slot(v, i);
        bucket_ptr(b).write(set_slot(v, i, cur_fp));
        cur_fp = victim;
        b ^= mix(cur_fp) & IDX_MASK; // involution: jumps to the alternate location
        let v2 = bucket_ptr(b).read();
        let mut j = 0u32;
        while j < 4 {
            if slot(v2, j) == EMPTY {
                bucket_ptr(b).write(set_slot(v2, j, cur_fp));
                STATS[5] += 1;
                STATS[7] += 1;
                return true;
            }
            j += 1;
        }
        kicks += 1;
    }
    STATS[6] += 1;
    false
}

#[inline(always)]
unsafe fn lookup_key(key: u64) -> bool {
    let (b1, b2, fp) = positions(key);
    if has_fp(bucket_ptr(b1).read(), fp) {
        return true;
    }
    has_fp(bucket_ptr(b2).read(), fp)
}

// ———————————————————— indexing a FASTA block ————————————————————

/// Insert every canonical k-mer of a FASTA block.
/// Header lines (`>`) break the k-mer context, as does any non-ACGT character.
/// Returns the number of bytes consumed; the caller carries the remainder over.
#[no_mangle]
pub unsafe extern "C" fn index_fasta_block(ptr: u64, len: u64, last: u32) -> u64 {
    let data = ptr as *const u8;
    let wmask: u128 = if W >= 64 { u128::MAX } else { (1u128 << (2 * W)) - 1 };

    let mut window: u128 = 0;
    let mut valid: u32 = 0;
    let mut i: u64 = 0;
    let mut consumed: u64 = 0;
    let mut in_header = false;

    while i < len {
        let b = data.add(i as usize).read();
        if b == b'\n' {
            in_header = false;
            consumed = i + 1;
            i += 1;
            continue;
        }
        if b == b'>' {
            in_header = true;
            valid = 0;
            i += 1;
            continue;
        }
        if in_header || b == b'\r' {
            i += 1;
            continue;
        }
        let c = code(b);
        if c == 255 {
            valid = 0;
        } else {
            window = ((window << 2) | (c as u128)) & wmask;
            valid += 1;
            if valid >= W {
                insert_key(canonical(window));
            }
        }
        i += 1;
    }
    // at end of stream everything is consumed
    if last != 0 { len } else { consumed }
}

// ———————————————————— filtering a FASTQ block ————————————————————

/// Filter a block of complete FASTQ records.
///
/// Kept records are written to `out_ptr`. A read is dropped as soon as the
/// fraction of its k-mers present in the index reaches `threshold_permille`.
/// `keep_matching` inverts the decision (extract matches instead of removing them).
///
/// Returns the number of input bytes consumed; the amount written is read back
/// through `last_out_len()`. The caller carries unconsumed bytes to the next block.
#[no_mangle]
pub unsafe extern "C" fn filter_fastq_block(
    ptr: u64, len: u64,
    out_ptr: u64, out_cap: u64,
    threshold_permille: u32,
    keep_matching: u32,
    last: u32,
) -> u64 {
    let data = ptr as *const u8;
    let out = out_ptr as *mut u8;
    let mut i: u64 = 0;
    let mut consumed: u64 = 0;
    let mut written: u64 = 0;

    loop {
        // delimit one 4-line record
        let start = i;
        let mut ends = [0u64; 4];
        let mut nl = 0;
        let mut j = i;
        while j < len && nl < 4 {
            if data.add(j as usize).read() == b'\n' {
                ends[nl] = j;
                nl += 1;
            }
            j += 1;
        }
        if nl < 4 {
            // incomplete record: stop here unless this is the very end
            if last != 0 && start < len {
                consumed = len;
            }
            break;
        }
        let seq_start = ends[0] + 1;
        let seq_end = ends[1];
        let contaminated = read_is_contaminated(data, seq_start, seq_end, threshold_permille);
        let keep = if keep_matching != 0 { contaminated } else { !contaminated };

        STATS[0] += 1;
        STATS[2] += seq_end - seq_start;
        if !keep {
            STATS[1] += 1;
        } else {
            let size = ends[3] + 1 - start;
            if written + size > out_cap {
                break; // output buffer full: hand control back
            }
            let mut n = 0u64;
            while n < size {
                out.add((written + n) as usize)
                    .write(data.add((start + n) as usize).read());
                n += 1;
            }
            written += size;
        }
        i = ends[3] + 1;
        consumed = i;
        if i >= len {
            break;
        }
    }
    LAST_OUT = written;
    consumed
}

static mut LAST_OUT: u64 = 0;

#[no_mangle]
pub unsafe extern "C" fn last_out_len() -> u64 {
    LAST_OUT
}

#[inline(always)]
/// Cleanifier's criterion is **base coverage**, not a fraction of seeds
/// (`cleanifier_filter.py: classify_issh`): every matching seed marks its
/// significant positions in a bitmap, and the read is dropped when the covered
/// fraction of its length exceeds the threshold. Counting matching seeds
/// instead under-filters by roughly 8 points on real human reads.
unsafe fn read_is_contaminated(data: *const u8, start: u64, end: u64, threshold_permille: u32) -> bool {
    let seq_len = end - start;
    if seq_len == 0 || seq_len > 1024 {
        return false;
    }
    let mut i = 0usize;
    while i < 16 { COVER[i] = 0; i += 1; }

    let wmask: u128 = if W >= 64 { u128::MAX } else { (1u128 << (2 * W)) - 1 };
    let mut window: u128 = 0;
    let mut valid: u32 = 0;
    let mut total: u64 = 0;
    let mut found: u64 = 0;

    let mut p = start;
    while p < end {
        let c = code(data.add(p as usize).read());
        if c == 255 {
            valid = 0;
        } else {
            window = ((window << 2) | (c as u128)) & wmask;
            valid += 1;
            if valid >= W {
                total += 1;
                if contains(canonical(window)) {
                    found += 1;
                    // the seed starting at this offset covers its significant positions
                    cover(p - start - (W as u64 - 1));
                }
            }
        }
        p += 1;
    }
    STATS[3] += total;
    STATS[4] += found;

    let nbits = cover_popcount(seq_len);
    // cleanifier: `threshold * seq_len < nbits_set`, strict
    nbits * 1000 > threshold_permille as u64 * seq_len
}

#[inline(always)]
unsafe fn cover(pos: u64) {
    let word = (pos >> 6) as usize;
    let off = pos & 63;
    if word >= 16 { return; }
    COVER[word] |= SEED_MASK_INT << off;
    if off > 0 && word + 1 < 16 {
        COVER[word + 1] |= SEED_MASK_INT >> (64 - off);
    }
}

#[inline(always)]
unsafe fn cover_popcount(seq_len: u64) -> u64 {
    let mut n = 0u64;
    let mut bits_left = seq_len;
    let mut i = 0usize;
    while i < 16 && bits_left > 0 {
        let take = if bits_left >= 64 { 64 } else { bits_left };
        let m = if take == 64 { u64::MAX } else { (1u64 << take) - 1 };
        n += (COVER[i] & m).count_ones() as u64;
        bits_left -= take;
        i += 1;
    }
    n
}

/// Fraction of k-mers present, in permille — used to decide a pair without
/// writing anything.
#[no_mangle]
pub unsafe extern "C" fn score_read(ptr: u64, start: u64, end: u64) -> u64 {
    let data = ptr as *const u8;
    let k = K;
    let shift = 2 * (k - 1);
    let kmask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fw, mut rc, mut valid, mut total, mut found) = (0u64, 0u64, 0u32, 0u64, 0u64);
    let mut p = start;
    while p < end {
        let c = code(data.add(p as usize).read());
        if c == 255 {
            valid = 0;
        } else {
            let c = c as u64;
            fw = ((fw << 2) | c) & kmask;
            rc = (rc >> 2) | ((3 - c) << shift);
            valid += 1;
            if valid >= k {
                total += 1;
                if lookup_key(if fw <= rc { fw } else { rc }) {
                    found += 1;
                }
            }
        }
        p += 1;
    }
    if total == 0 { 0 } else { (found * 1000) / total }
}

// ——————————————— Cleanifier's native windowed-cuckoo format ———————————————
//
// Read path for the published `.filter` files. Ported from a Python reference
// that agrees with Cleanifier on 60 000/60 000 keys of the human index,
// false positives included.
//
// All parameters come from the `.info` sidecar; nothing is guessed here.

static mut NATIVE: bool = false;
static mut N_BASE: u64 = 0;            // byte offset of the filter array
static mut N_NSUB: u64 = 1;
static mut N_NWINDOWS: u64 = 0;
static mut N_WINSIZE: u64 = 2;
static mut N_BITS_SLOT: u64 = 16;
static mut N_BITS_WINDOW: u64 = 32;
static mut N_BITS_SUBFILTER: u64 = 0;
static mut N_NFP: u64 = 0;             // 2^fpr_bits - 1
static mut N_Q: u32 = 29;              // half of log2(universe)
static mut N_QBITS: u32 = 58;
static mut N_CODEMASK: u64 = 0;
static mut N_A_SUB: u64 = 0;
static mut N_A_WIN: u64 = 0;
static mut N_A_FP: u64 = 0;
static mut N_A_OFF: u64 = 0;
static mut N_M1: u64 = 0;              // SWAR masks over 2*windowsize slots
static mut N_M2: u64 = 0;
static mut N_FPM: u64 = 0;

/// Configure the native reader. `a_*` are the multipliers of the four
/// `linear<odd>` hash functions, in the order subfilter/window/fingerprint/offset.
#[no_mangle]
pub unsafe extern "C" fn native_init(
    base: u64, universe_bits: u32, nsub: u64, nslots: u64, nwindows: u64,
    windowsize: u64, fpr_bits: u32,
    a_sub: u64, a_win: u64, a_fp: u64, a_off: u64,
) -> u64 {
    N_BASE = base;
    N_QBITS = universe_bits;
    N_Q = universe_bits / 2;
    N_CODEMASK = if universe_bits >= 64 { u64::MAX } else { (1u64 << universe_bits) - 1 };
    N_NSUB = nsub;
    N_NWINDOWS = nwindows;
    N_WINSIZE = windowsize;
    N_NFP = (1u64 << fpr_bits) - 1;

    let window_bits = bits_for(windowsize);
    N_BITS_SLOT = fpr_bits as u64 + window_bits + 1;
    N_BITS_WINDOW = N_BITS_SLOT * windowsize;

    // nslots is the total across subfilters; each one is padded up to 512 bits
    let per_sub = (nslots + nsub - 1) / nsub;
    let mut bps = per_sub * N_BITS_SLOT;
    bps += 512 - (bps % 512);
    N_BITS_SUBFILTER = bps;

    N_A_SUB = a_sub; N_A_WIN = a_win; N_A_FP = a_fp; N_A_OFF = a_off;

    // SWAR masks spanning both windows (2 * windowsize slots)
    let slots = 2 * windowsize;
    let (mut m1, mut m2) = (0u64, 0u64);
    let mut s = 0u64;
    while s < slots {
        m1 |= 1u64 << (s * N_BITS_SLOT);
        m2 |= 1u64 << (s * N_BITS_SLOT + N_BITS_SLOT - 1);
        s += 1;
    }
    N_M1 = m1;
    N_M2 = m2;

    // compute_choice_fp_mask: slot index in the offset bits; the choice bit is
    // set on the low half, which carries window 2.
    let mut fpm = 0u64;
    let fpb = fpr_bits as u64;
    s = 0;
    while s < windowsize {
        fpm |= s << (s * N_BITS_SLOT + fpb);
        fpm |= 1u64 << (s * N_BITS_SLOT + fpb + window_bits);
        s += 1;
    }
    s = 0;
    while s < windowsize {
        fpm |= s << ((windowsize + s) * N_BITS_SLOT + fpb);
        s += 1;
    }
    N_FPM = fpm;

    NATIVE = true;
    N_BITS_SUBFILTER
}

#[inline(always)]
fn bits_for(k: u64) -> u64 {
    if k == 0 { return 0; }
    let mut n = 0u64;
    while (1u64 << n) < k { n += 1; }
    n
}

/// The swap-and-multiply shared by every `linear<odd>` hash.
#[inline(always)]
unsafe fn swap_mul(a: u64, code: u64) -> u64 {
    let swap = ((code << N_Q) ^ (code >> N_Q)) & N_CODEMASK;
    a.wrapping_mul(swap) & N_CODEMASK
}

/// `compile_get_bucket` has three branches depending on nbuckets; reproducing
/// all three matters — the offset hash lands on the "even" one.
#[inline(always)]
unsafe fn get_bucket(a: u64, code: u64, nbuckets: u64) -> u64 {
    let mut swap = swap_mul(a, code);
    if nbuckets & (nbuckets - 1) == 0 {
        let bits = bits_for(nbuckets);
        if bits == 0 { return 0; }
        return (swap >> (N_QBITS as u64 - bits)) & ((1u64 << bits) - 1);
    }
    if nbuckets % 2 == 0 {
        swap ^= swap >> N_Q;
    }
    swap % nbuckets
}

/// Reads `N_BITS_WINDOW` bits starting at bit `pos`. The load is unaligned by
/// design: Cleanifier packs slots without padding, so a window straddles bytes.
#[inline(always)]
unsafe fn window_bits(pos: u64) -> u64 {
    let p = (N_BASE + (pos >> 3)) as *const u64;
    let chunk = p.read_unaligned();
    (chunk >> (pos & 7)) & ((1u64 << N_BITS_WINDOW) - 1)
}

#[inline(always)]
unsafe fn native_lookup(code: u64) -> bool {
    let fp = get_bucket(N_A_FP, code, N_NFP) + 1;
    let w1 = get_bucket(N_A_WIN, code, N_NWINDOWS);
    let off = get_bucket(N_A_OFF, fp, N_NWINDOWS - N_WINSIZE + 1);
    let w2 = (w1 + off) % N_NWINDOWS;
    let sf = get_bucket(N_A_SUB, code, N_NSUB);
    let base = sf * N_BITS_SUBFILTER;

    let wb1 = window_bits(base + w1 * N_BITS_SLOT);
    let wb2 = window_bits(base + w2 * N_BITS_SLOT);
    let combined = (wb1 << N_BITS_WINDOW) | wb2;

    // hasvalue(): one SWAR sweep across both windows
    let x = combined ^ (N_M1.wrapping_mul(fp) | N_FPM);
    (x.wrapping_sub(N_M1) & !x & N_M2) != 0
}

/// Single entry point so the FASTQ path does not care which index is loaded.
#[inline(always)]
unsafe fn contains(code: u64) -> bool {
    if NATIVE { native_lookup(code) } else { lookup_key(code) }
}

/// Set the seed without allocating an index — used when the seed comes from a
/// `.info` file rather than from the UI.
#[no_mangle]
pub unsafe extern "C" fn set_seed(seed_bits: u64, w: u32) {
    SEED_MASK_INT = seed_bits;
    W = w;
    K = seed_bits.count_ones();
    CONTIGUOUS = K == W;
    let mut i = 0u32;
    let mut n = 0usize;
    while i < w {
        if (seed_bits >> i) & 1 == 1 {
            SEED_POS[n] = i as u8;
            n += 1;
        }
        i += 1;
    }
}

/// First canonical code of a raw ACGT buffer — lets a test compare the encoder
/// against the Python reference without going through the index at all.
#[no_mangle]
pub unsafe extern "C" fn first_canonical(ptr: u64, len: u64) -> u64 {
    let data = ptr as *const u8;
    let wmask: u128 = if W >= 64 { u128::MAX } else { (1u128 << (2 * W)) - 1 };
    let mut window: u128 = 0;
    let mut valid: u32 = 0;
    let mut p = 0u64;
    while p < len {
        let c = code(data.add(p as usize).read());
        if c == 255 { valid = 0; } else {
            window = ((window << 2) | (c as u128)) & wmask;
            valid += 1;
            if valid >= W { return canonical(window); }
        }
        p += 1;
    }
    u64::MAX
}

/// Forward gather only, no canonicalisation — isolates revcomp from gather.
#[no_mangle]
pub unsafe extern "C" fn first_forward(ptr: u64, len: u64) -> u64 {
    let data = ptr as *const u8;
    let wmask: u128 = if W >= 64 { u128::MAX } else { (1u128 << (2 * W)) - 1 };
    let mut window: u128 = 0;
    let mut valid: u32 = 0;
    let mut p = 0u64;
    while p < len {
        let c = code(data.add(p as usize).read());
        if c == 255 { valid = 0; } else {
            window = ((window << 2) | (c as u128)) & wmask;
            valid += 1;
            if valid >= W { return gather(window); }
        }
        p += 1;
    }
    u64::MAX
}

/// Probe one key through whichever index is active (debugging aid for the UI).
#[no_mangle]
pub unsafe extern "C" fn probe(code: u64) -> u32 {
    if contains(code) { 1 } else { 0 }
}

// ————————————————————————————— misc —————————————————————————————

#[no_mangle]
pub unsafe extern "C" fn stat(i: u32) -> u64 {
    if i < 8 { STATS[i as usize] } else { 0 }
}

#[no_mangle]
pub unsafe extern "C" fn reset_stats() {
    let mut s = 0;
    while s < 8 {
        STATS[s] = 0;
        s += 1;
    }
}

/// Fill rate in permille, so the UI can warn about an undersized index.
#[no_mangle]
pub unsafe extern "C" fn load_permille() -> u64 {
    let slots = (IDX_MASK + 1) * 4;
    if slots == 0 { 0 } else { (STATS[7] * 1000) / slots }
}

#[no_mangle]
pub extern "C" fn mem_pages() -> u64 {
    core::arch::wasm64::memory_size(0) as u64
}

#[no_mangle]
pub extern "C" fn mem_grow(pages: u64) -> u64 {
    let prev = core::arch::wasm64::memory_grow(0, pages as usize);
    if prev == usize::MAX { u64::MAX } else { prev as u64 }
}
