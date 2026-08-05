#!/usr/bin/env python3
"""Lecteur indépendant du format d'index Cleanifier (windowed cuckoo filter).

Réimplémente le chemin de requête à partir du seul .filter + .info, sans
importer cleanifier. Sert à prouver que le format est correctement décodé
avant de le porter en Rust.
"""
import pickle, sys
import numpy as np

MASK64 = (1 << 64) - 1


def bitsfor(k):
    if k == 0:
        return 0
    n = 0
    while (1 << n) < k:
        n += 1
    return n


class CleanifierIndex:
    def __init__(self, prefix):
        with open(f"{prefix}.info", "rb") as f:
            finfo, optinfo, appinfo = pickle.load(f)
        self.finfo, self.appinfo = finfo, appinfo

        self.universe = int(finfo["universe"])
        self.qbits = bitsfor(self.universe)
        assert 4 ** (self.qbits // 2) == self.universe
        self.q = self.qbits // 2
        self.codemask = (1 << self.qbits) - 1

        self.nsubfilters = int(finfo["nsubfilters"])
        self.windowsize = int(finfo["windowsize"])
        self.nwindows = int(finfo["nwindows_per_subfilter"])
        self.fpr_bits = int(finfo["target_fpr"])
        self.nfingerprints = 2 ** self.fpr_bits - 1

        self.window_bits = bitsfor(self.windowsize)
        self.bits_per_slot = self.fpr_bits + self.window_bits + 1
        self.bits_per_window = self.bits_per_slot * self.windowsize
        self.choice_mask = 1 << (self.fpr_bits + self.window_bits)

        nslots = int(finfo["nslots"])
        nslots_per_sub = -(-nslots // self.nsubfilters)
        bps = nslots_per_sub * self.bits_per_slot
        bps += 512 - (bps % 512)
        self.bits_per_subfilter = bps

        self.k = appinfo["k"]
        self.seed = tuple(appinfo["mask"])
        self.span = max(self.seed) + 1
        self.rcmode = appinfo["rcmode"]

        # les quatre fonctions de hachage, chacune "linear<odd>"
        self.a_sub = int(finfo["subfilter_hashfunc_str"][6:])
        self.a_win = int(finfo["window_hashfunc_str"][6:])
        self.a_fp = int(finfo["fingerprint_hashfunc_str"][6:])
        self.a_off = int(finfo["offset_hashfunc_str"][6:])

        self.arr = np.fromfile(f"{prefix}.filter", dtype=np.uint64)
        self.raw = self.arr.view(np.uint8)

    # ————————————————————— hachage —————————————————————

    def _swap_mul(self, a, code):
        swap = ((code << self.q) ^ (code >> self.q)) & self.codemask
        return (a * swap) & self.codemask

    def _get_bucket(self, a, code, nbuckets):
        """compile_get_bucket : trois variantes selon nbuckets."""
        swap = self._swap_mul(a, code)
        if nbuckets & (nbuckets - 1) == 0:            # puissance de 2
            bits = bitsfor(nbuckets)
            if bits == 0:
                return 0
            return (swap >> (self.qbits - bits)) & ((1 << bits) - 1)
        if nbuckets % 2 == 0:                          # pair
            swap ^= (swap >> self.q)
            return swap % nbuckets
        return swap % nbuckets                         # impair

    def subfilter(self, code):
        return self._get_bucket(self.a_sub, code, self.nsubfilters)

    def fingerprint(self, code):
        return self._get_bucket(self.a_fp, code, self.nfingerprints) + 1

    def window1(self, code):
        return self._get_bucket(self.a_win, code, self.nwindows)

    def offset(self, fp):
        return self._get_bucket(self.a_off, fp, self.nwindows - self.windowsize + 1)

    # ————————————————————— accès binaire —————————————————————

    def get_window_bits(self, pos):
        """compile_load_value : charge un u64 à l'octet pos>>3 puis décale de pos&7."""
        byte = pos >> 3
        chunk = int.from_bytes(self.raw[byte:byte + 16].tobytes(), "little")
        return (chunk >> (pos & 7)) & ((1 << self.bits_per_window) - 1)

    def _masks(self, slots):
        m1 = int(("0" * (self.bits_per_slot - 1) + "1") * slots, 2)
        m2 = int(("1" + "0" * (self.bits_per_slot - 1)) * slots, 2)
        return m1, m2

    def has_fp(self, window, fp):
        """hasvalue() : recherche SWAR de fp dans les slots de la fenêtre.

        Toute l'arithmétique doit rester en uint64 : en Python, `~x` sur un
        entier positif donne un négatif non borné et `x - m1` peut passer sous
        zéro, alors que numba travaille modulo 2^64.
        """
        slots = self.windowsize
        m1, m2 = self._masks(slots)
        fp_mask = 0
        for s in range(1, slots):
            fp_mask |= s << (s * self.bits_per_slot + self.fpr_bits)
        x = (window ^ (((m1 * fp) & MASK64 | fp_mask))) & MASK64
        return (((x - m1) & MASK64) & (~x & MASK64) & m2) != 0

    # ————————————————————— requête —————————————————————

    def _choice_fp_mask(self):
        """compute_choice_fp_mask : slot number in the offset bits, plus the
        choice bit set on the low half (which carries window 2)."""
        slots, bps, fpb = self.windowsize, self.bits_per_slot, self.fpr_bits
        wb = self.window_bits
        m = 0
        for s in range(slots):
            m |= s << (s * bps + fpb)
            m |= 1 << (s * bps + fpb + wb)
        for s in range(slots):
            m |= s << ((slots + s) * bps + fpb)
        return m

    def has_fp_combined(self, window64, fp):
        """compile_lookup(..., choice=True) : one SWAR test over both windows."""
        m1, m2 = self._masks(2 * self.windowsize)
        fpm = self._choice_fp_mask()
        x = (window64 ^ (((m1 * fp) & MASK64 | fpm))) & MASK64
        return (((x - m1) & MASK64) & (~x & MASK64) & m2) != 0

    def lookup(self, code):
        """lookup_in_subfilter_fast_combined: cleanifier wires this variant
        whenever bits_per_slot * windowsize * 2 <= 64, which holds here (64)."""
        fp = self.fingerprint(code)
        w1 = self.window1(code)
        w2 = (w1 + self.offset(fp)) % self.nwindows
        base = self.subfilter(code) * self.bits_per_subfilter
        wb1 = self.get_window_bits(base + w1 * self.bits_per_slot)
        wb2 = self.get_window_bits(base + w2 * self.bits_per_slot)
        combined = ((wb1 << self.bits_per_window) | wb2) & MASK64
        return self.has_fp_combined(combined, fp)

    # ————————————————————— k-mers —————————————————————

    def canonical_codes(self, seq):
        """Codes canoniques des seeds espacées d'une séquence ACGT."""
        tr = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
        span, seed, k = self.span, self.seed, self.k
        vals = [tr.get(c, -1) for c in seq]
        out = []
        for i in range(len(vals) - span + 1):
            win = vals[i:i + span]
            if -1 in win:
                continue
            fw = 0
            for p in seed:
                fw = (fw << 2) | win[p]
            rc = 0
            for p in reversed(seed):
                rc = (rc << 2) | (3 - win[p])
            out.append(max(fw, rc) if self.rcmode == "max" else min(fw, rc))
        return out


if __name__ == "__main__":
    prefix = sys.argv[1]
    fasta = sys.argv[2]
    idx = CleanifierIndex(prefix)
    print(f"k={idx.k} span={idx.span} rcmode={idx.rcmode}")
    print(f"subfilters={idx.nsubfilters} windows={idx.nwindows} windowsize={idx.windowsize}")
    print(f"bits/slot={idx.bits_per_slot} bits/subfilter={idx.bits_per_subfilter} fpr_bits={idx.fpr_bits}")
    print(f"filter array: {idx.arr.size} x uint64 = {idx.arr.nbytes} octets")

    seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">"))
    codes = idx.canonical_codes(seq)
    print(f"\nseeds extraites du FASTA : {len(codes):,}")

    hits = sum(idx.lookup(c) for c in codes[:20000])
    n = min(20000, len(codes))
    print(f"retrouvées dans l'index  : {hits:,}/{n:,}  ({100*hits/n:.2f} %)")

    rng = np.random.default_rng(3)
    rand = [int(x) % idx.universe for x in rng.integers(0, 2**62, 20000)]
    fp = sum(idx.lookup(c) for c in rand)
    print(f"faux positifs (aléatoire): {fp:,}/20,000  ({100*fp/20000:.3f} %)"
          f"  [attendu ~{100*2*idx.windowsize/2**idx.fpr_bits:.3f} %]")
