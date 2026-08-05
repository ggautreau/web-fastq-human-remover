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

import { readInfo } from './pickle.js';

// worker.js — the decontamination pipeline. All heavy work happens here;
// the main thread only drives the UI.
//
// Nothing is ever materialised: FASTQ files stream through
// (read -> gunzip -> WASM -> write). The only resident memory is the index
// plus two buffers of a few MB, whatever the input size.

const PAGE = 65536;
const IN_CAP = 32 << 20;   // input buffer
const OUT_CAP = 40 << 20;  // output: headroom, a whole block may be kept

// IO_BASE must NOT be hardcoded. wasm-ld lays the stack and static data at the
// bottom of linear memory and their exact extent is only known after linking,
// so we read __heap_base, which the module exports for this purpose.
// (A guessed value of 1 MiB landed 112 bytes inside Rust's own area and
// silently corrupted every counter — no crash, no error.)
let IO_BASE = 0, IDX_BASE = 0;

const align = n => Number((BigInt(n) + 65535n) & ~65535n);

// The module imports its memory with an explicit maximum, so every Memory we
// hand it must declare the same ceiling. Nothing is *shared*: a single worker
// owns the index, which means the page needs no COOP/COEP headers and can be
// served from plain static hosting such as GitHub Pages.
const MAX_PAGES = 262144n;   // 16 GB, V8's memory64 ceiling

let mem = null, wasm = null, view = null;
let indexReady = false;
let cancelled = false;

const post = (type, data) => postMessage({ type, ...data });

function human(n) {
  const u = ['B', 'KB', 'MB', 'GB'];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
  return `${n.toFixed(n < 10 && i > 0 ? 1 : 0)} ${u[i]}`;
}

const popcount = b => { let n = 0; while (b) { n += Number(b & 1n); b >>= 1n; } return n; };

// ————————————————————————————— setup —————————————————————————————

async function init({ wasmUrl, log2Buckets, seedBits, span }) {
  const indexBytes = (2 ** log2Buckets) * 8;
  const bytes = await (await fetch(wasmUrl)).arrayBuffer();
  const mod = await WebAssembly.compile(bytes);

  // probe instantiation, solely to learn where Rust's data ends
  const probe = new WebAssembly.Memory({ initial: 32n, maximum: MAX_PAGES, address: 'i64' });
  const heapBase = Number((await WebAssembly.instantiate(mod, { env: { memory: probe } }))
                            .exports.__heap_base.value);
  IO_BASE = align(heapBase);
  IDX_BASE = IO_BASE + IN_CAP + OUT_CAP;

  const total = IDX_BASE + indexBytes;
  const pages = Math.ceil(total / PAGE);

  mem = new WebAssembly.Memory({ initial: BigInt(pages), maximum: MAX_PAGES, address: 'i64' });
  wasm = (await WebAssembly.instantiate(mod, { env: { memory: mem } })).exports;
  view = new Uint8Array(mem.buffer);

  wasm.idx_init(BigInt(IDX_BASE), log2Buckets, BigInt(seedBits), span);
  indexReady = false;

  post('ready', {
    memory: human(total),
    index: human(indexBytes),
    buckets: 2 ** log2Buckets,
    kmerCapacity: (2 ** log2Buckets) * 4,
    weight: popcount(BigInt(seedBits)),
    span
  });
}

// ————————————————— fetching a published index from Zenodo —————————————————
//
// The plain download URL (zenodo.org/records/.../files/...) sends no CORS
// headers, so the browser cannot read it. The **API** URL for the same file
// does send `Access-Control-Allow-Origin: *`, which is what makes this work
// without a proxy — nothing is relayed through a third party.

const ZENODO = (record, name) =>
  `https://zenodo.org/api/records/${record}/files/${name}/content`;

// Reads a POSIX tar header block. Returns null on the end-of-archive marker.
function tarHeader(block) {
  if (block[0] === 0) return null;
  const dec = new TextDecoder();
  const str = (o, n) => dec.decode(block.subarray(o, o + n)).replace(/\0.*$/, '').trim();
  const size = parseInt(str(124, 12), 8) || 0;
  return { name: str(0, 100), size, type: String.fromCharCode(block[156] || 48) };
}

// Streams the tar, routing the .filter straight into linear memory and keeping
// the small .info aside. The archive is never held in JS.
async function fetchIndexFromZenodo({ wasmUrl, record, archive, saveDir }) {
  cancelled = false;

  // Single request, no preliminary Range probe: Zenodo does not advertise
  // accept-ranges, and a range request on a 6 GB file makes the server walk the
  // whole thing — it stalled for 20 s before sending a byte. We open the stream
  // once and read the first tar header out of it instead.
  post('progress', { phase: 'index', fraction: 0, detail: 'contacting Zenodo…' });
  const reader = (await fetch(ZENODO(record, archive))).body.getReader();

  const concat = (a, b) => { const c = new Uint8Array(a.length + b.length); c.set(a); c.set(b, a.length); return c; };
  let pending = new Uint8Array(0);

  // pull until the first 512-byte header is complete
  while (pending.length < 512) {
    const { done, value } = await reader.read();
    if (done) throw new Error('archive ended before its first header');
    pending = concat(pending, value);
  }
  const first = tarHeader(pending.subarray(0, 512));
  if (!first) throw new Error('unreadable tar header');
  // The first header already gives the format, for free: .filter is the
  // windowed cuckoo filter we read, .hash is the exact index (a 3-way bucketed
  // hash table) which this build does not implement.
  if (/\.hash$/i.test(first.name)) {
    await reader.cancel();
    throw new Error('this is the exact index (a hash table), which is not supported — '
      + 'use the probabilistic one instead');
  }
  if (!/\.filter$/i.test(first.name)) {
    await reader.cancel();
    throw new Error(`unexpected first archive entry '${first.name}'`);
  }
  post('progress', { phase: 'index', fraction: 0, detail: `${first.name} — ${human(first.size)}` });

  const bytes = await (await fetch(wasmUrl)).arrayBuffer();
  const mod = await WebAssembly.compile(bytes);
  const probeMem = new WebAssembly.Memory({ initial: 32n, maximum: MAX_PAGES, address: 'i64' });
  const heapBase = Number((await WebAssembly.instantiate(mod, { env: { memory: probeMem } }))
                            .exports.__heap_base.value);
  IO_BASE = align(heapBase);
  IDX_BASE = IO_BASE + IN_CAP + OUT_CAP;
  const pages = Math.ceil((IDX_BASE + first.size + PAGE) / PAGE);
  mem = new WebAssembly.Memory({ initial: BigInt(pages), maximum: MAX_PAGES, address: 'i64' });
  wasm = (await WebAssembly.instantiate(mod, { env: { memory: mem } })).exports;
  view = new Uint8Array(mem.buffer);

  const t0 = performance.now();
  let entry = null, written = 0, infoBytes = null, filterLen = 0;
  let total = pending.length, lastPost = 0;
  let eof = false;
  // Optional: mirror the extracted entries to disk as they stream past, so the
  // index can be re-loaded later without downloading it again. Costs no extra
  // transfer — the same bytes are written out on their way to linear memory.
  let saveHandle = null;
  const opened = [];
  if (saveDir) {
    saveHandle = async (name) => {
      const fh = await saveDir.getFileHandle(name, { create: true });
      const w = await fh.createWritable();
      opened.push({ name, w });
      return w;
    };
  }

  for (;;) {
    // consume everything currently buffered
    for (;;) {
      if (!entry) {
        if (pending.length < 512) break;
        const h = tarHeader(pending.subarray(0, 512));
        pending = pending.subarray(512);
        if (!h) { eof = true; break; }
        entry = { ...h, left: h.size, pad: (512 - (h.size % 512)) % 512 };
        if (/\.filter$/i.test(entry.name)) { entry.kind = 'filter'; written = 0; filterLen = entry.size; }
        else if (/\.info$/i.test(entry.name)) { entry.kind = 'info'; infoBytes = new Uint8Array(entry.size); }
        else entry.kind = 'skip';
        if (saveHandle && entry.kind !== 'skip') entry.sink = await saveHandle(entry.name);
        continue;
      }
      const take = Math.min(entry.left, pending.length);
      if (take > 0) {
        const chunk = pending.subarray(0, take);
        if (entry.kind === 'filter') { view.set(chunk, IDX_BASE + written); written += take; }
        else if (entry.kind === 'info') { infoBytes.set(chunk, entry.size - entry.left); }
        if (entry.sink) await entry.sink.write(chunk.slice());
        entry.left -= take;
        pending = pending.subarray(take);
      }
      if (entry.left > 0) break;
      const p = Math.min(entry.pad, pending.length);
      pending = pending.subarray(p);
      entry.pad -= p;
      if (entry.pad > 0) break;
      entry = null;
    }
    if (eof) break;

    if (cancelled) {
      await reader.cancel();
      for (const o of opened) { try { await o.w.abort(); } catch {} }
      post('cancelled', {}); return;
    }
    const { done, value } = await reader.read();
    if (done) break;
    total += value.length;
    pending = pending.length ? concat(pending, value) : value;

    // throttled: one message every ~80 ms keeps the UI thread free
    const now = performance.now();
    if (now - lastPost > 80) {
      lastPost = now;
      const mbps = total / 1048576 / ((now - t0) / 1000);
      post('progress', {
        phase: 'index',
        fraction: filterLen ? written / filterLen : 0,
        detail: `${human(written)} of ${human(filterLen)} · ${mbps.toFixed(0)} MB/s`
      });
    }
  }

  if (cancelled) { post('cancelled', {}); return; }
  for (const o of opened) await o.w.close();
  if (!infoBytes) throw new Error('no .info entry found in the archive');
  const seconds = (performance.now() - t0) / 1000;
  const meta = readInfo(infoBytes);
  applyNative(meta);
  indexReady = true;
  post('nativeDone', {
    seconds: +seconds.toFixed(1), size: human(written),
    rate: `${Math.round(total / 1048576 / seconds)} MB/s`,
    weight: meta.k, span: meta.span, subfilters: meta.nsubfilters,
    windows: meta.nwindows, fprBits: meta.fprBits, rcmode: meta.rcmode,
    saved: opened.map(o => o.name)
  });
}

/// Pushes the parameters read from a .info into the wasm core.
function applyNative(meta) {
  if (meta.filtertype !== 'windowed_cuckoo')
    throw new Error(`unsupported index type '${meta.filtertype}'`);
  const seedBits = meta.seed.reduce((a, p) => a | (1n << BigInt(p)), 0n);
  wasm.set_seed(seedBits, meta.span);
  wasm.native_init(
    BigInt(IDX_BASE), Math.round(Math.log2(meta.universe)), BigInt(meta.nsubfilters),
    BigInt(meta.nslots), BigInt(meta.nwindows), BigInt(meta.windowsize), meta.fprBits,
    BigInt(parseInt(meta.hash.subfilter.slice(6))), BigInt(parseInt(meta.hash.window.slice(6))),
    BigInt(parseInt(meta.hash.fingerprint.slice(6))), BigInt(parseInt(meta.hash.offset.slice(6))));
}

// ————————————————— loading Cleanifier's native index —————————————————

// Streams the .filter straight into linear memory. The file is never
// materialised in JS: chunks land directly at their final offset.
async function loadNative({ wasmUrl, filter, info }) {
  cancelled = false;
  const meta = readInfo(new Uint8Array(await info.arrayBuffer()));
  if (meta.filtertype !== 'windowed_cuckoo')
    throw new Error(`unsupported index type '${meta.filtertype}'`);

  const bytes = await (await fetch(wasmUrl)).arrayBuffer();
  const mod = await WebAssembly.compile(bytes);
  const probe = new WebAssembly.Memory({ initial: 32n, maximum: MAX_PAGES, address: 'i64' });
  const heapBase = Number((await WebAssembly.instantiate(mod, { env: { memory: probe } }))
                            .exports.__heap_base.value);
  IO_BASE = align(heapBase);
  IDX_BASE = IO_BASE + IN_CAP + OUT_CAP;

  const total = IDX_BASE + filter.size + PAGE;
  const pages = Math.ceil(total / PAGE);
  mem = new WebAssembly.Memory({ initial: BigInt(pages), maximum: MAX_PAGES, address: 'i64' });
  wasm = (await WebAssembly.instantiate(mod, { env: { memory: mem } })).exports;
  view = new Uint8Array(mem.buffer);

  const t0 = performance.now();
  const reader = filter.stream().getReader();
  let off = 0;
  for (;;) {
    if (cancelled) { await reader.cancel(); post('cancelled', {}); return; }
    const { done, value } = await reader.read();
    if (done) break;
    view.set(value, IDX_BASE + off);
    off += value.length;
    post('progress', { phase: 'index', fraction: off / filter.size,
                       detail: `${human(off)} of ${human(filter.size)}` });
  }
  const seconds = (performance.now() - t0) / 1000;

  applyNative(meta);
  indexReady = true;
  post('nativeDone', {
    seconds: +seconds.toFixed(1),
    size: human(off),
    rate: seconds > 0.2 ? `${Math.round(off / 1048576 / seconds)} MB/s` : 'instant',
    weight: meta.k, span: meta.span, subfilters: meta.nsubfilters,
    windows: meta.nwindows, fprBits: meta.fprBits, rcmode: meta.rcmode
  });
}

// ————————————————————————— building the index —————————————————————————

function decodedStream(file) {
  let s = file.stream();
  if (/\.(gz|bgz)$/i.test(file.name)) s = s.pipeThrough(new DecompressionStream('gzip'));
  return s.getReader();
}

async function buildIndex({ files }) {
  cancelled = false;
  let readTotal = 0;
  const totalSize = files.reduce((s, f) => s + f.size, 0);
  const t0 = performance.now();

  for (const f of files) {
    const reader = decodedStream(f);
    let carry = 0;   // bytes carried over at the head of the buffer

    for (;;) {
      if (cancelled) { await reader.cancel(); post('cancelled', {}); return; }
      const { done, value } = await reader.read();
      if (done) break;
      readTotal += value.length;

      let pos = 0;
      while (pos < value.length) {
        const n = Math.min(IN_CAP - carry, value.length - pos);
        view.set(value.subarray(pos, pos + n), IO_BASE + carry);
        pos += n;
        const avail = carry + n;
        const consumed = Number(wasm.index_fasta_block(BigInt(IO_BASE), BigInt(avail), 0));
        carry = avail - consumed;
        if (carry > 0 && consumed > 0) view.copyWithin(IO_BASE, IO_BASE + consumed, IO_BASE + avail);
        if (carry === IN_CAP) carry = 0;  // line longer than the buffer: cut it
      }
      post('progress', {
        phase: 'index',
        fraction: totalSize ? readTotal / totalSize : 0,
        detail: `${human(readTotal)} read`
      });
    }
    if (carry > 0) wasm.index_fasta_block(BigInt(IO_BASE), BigInt(carry), 1);
  }

  indexReady = true;
  const load = Number(wasm.load_permille()) / 10;
  post('indexDone', {
    seconds: +((performance.now() - t0) / 1000).toFixed(1),
    kmers: Number(wasm.stat(5)),
    failed: Number(wasm.stat(6)),
    loadPct: +load.toFixed(1),
    warning: load > 90
      ? 'Index nearly saturated — increase its size, some k-mers were dropped.'
      : null
  });
}

// ————————————————————————— filtering samples —————————————————————————

async function filter({ samples, thresholdPermille, keepMatching, outputDir, gzipOutput }) {
  cancelled = false;
  if (!indexReady) { post('error', { message: 'No index has been built yet.' }); return; }

  for (const s of samples) {
    const t0 = performance.now();
    wasm.reset_stats();
    let paired = null;
    try {
      if (s.r2) paired = await filterPair(s, thresholdPermille, keepMatching, outputDir, gzipOutput);
      else await filterSingle(s, thresholdPermille, keepMatching, outputDir, gzipOutput);
    } catch (e) {
      post('error', { message: `${s.name}: ${e.name} — ${e.message}` });
      continue;
    }
    if (cancelled) { post('cancelled', {}); return; }
    post('sampleDone', {
      name: s.name,
      seconds: +((performance.now() - t0) / 1000).toFixed(1),
      reads: paired ? paired.reads : Number(wasm.stat(0)),
      dropped: paired ? paired.dropped : Number(wasm.stat(1)),
      paired: !!paired
    });
  }
  post('done', {});
}

// Uniform sink whether or not we compress. With gzip the bytes go through a
// CompressionStream piped into the file, so the stream stays a stream — nothing
// extra is buffered, whatever the output size.
async function openOutput(dir, name, gzip) {
  const finalName = gzip && !/\.gz$/i.test(name) ? `${name}.gz` : name;
  const fh = await dir.getFileHandle(finalName, { create: true });
  const file = await fh.createWritable();
  if (!gzip) {
    return {
      write: c => file.write(c),
      close: () => file.close(),
      abort: () => file.abort(),
    };
  }
  const cs = new CompressionStream('gzip');
  const done = cs.readable.pipeTo(file);          // resolves when the file is fully written
  const w = cs.writable.getWriter();
  return {
    write: c => w.write(c),
    close: async () => { await w.close(); await done; },
    abort: async () => { try { await w.abort(); } catch {} try { await file.abort(); } catch {} },
  };
}

async function filterSingle(sample, threshold, keepMatching, dir, gzip) {
  const w = await openOutput(dir, outputName(sample.r1.name, keepMatching), gzip);
  const reader = decodedStream(sample.r1);
  let carry = 0, read = 0;
  try {
    for (;;) {
      if (cancelled) { await reader.cancel(); await w.abort(); return; }
      const { done, value } = await reader.read();
      if (done) break;
      read += value.length;
      let pos = 0;
      while (pos < value.length) {
        const n = Math.min(IN_CAP - carry, value.length - pos);
        view.set(value.subarray(pos, pos + n), IO_BASE + carry);
        pos += n;
        const avail = carry + n;
        const consumed = Number(wasm.filter_fastq_block(
          BigInt(IO_BASE), BigInt(avail),
          BigInt(IO_BASE + IN_CAP), BigInt(OUT_CAP),
          threshold, keepMatching ? 1 : 0, 0));
        const out = Number(wasm.last_out_len());
        // slice(), not subarray(): a view on a SharedArrayBuffer cannot be written out
        if (out > 0) await w.write(view.slice(IO_BASE + IN_CAP, IO_BASE + IN_CAP + out));
        carry = avail - consumed;
        if (carry > 0 && consumed > 0) view.copyWithin(IO_BASE, IO_BASE + consumed, IO_BASE + avail);
        if (consumed === 0 && carry === IN_CAP) throw new Error('record larger than the 32 MiB buffer');
      }
      post('progress', {
        phase: 'filter', name: sample.name,
        fraction: sample.r1.size ? read / sample.r1.size : 0,
        detail: `${Number(wasm.stat(0)).toLocaleString('en-US')} reads`
      });
    }
    if (carry > 0) {
      wasm.filter_fastq_block(BigInt(IO_BASE), BigInt(carry),
        BigInt(IO_BASE + IN_CAP), BigInt(OUT_CAP), threshold, keepMatching ? 1 : 0, 1);
      const out = Number(wasm.last_out_len());
      if (out > 0) await w.write(view.slice(IO_BASE + IN_CAP, IO_BASE + IN_CAP + out));
    }
    await w.close();
  } catch (e) { try { await w.abort(); } catch {} throw e; }
}

// Paired-end: one joint decision. A pair leaves or stays whole — dropping a
// read without its mate would break mate pairing for every downstream tool.
async function filterPair(sample, threshold, keepMatching, dir, gzip) {
  const w1 = await openOutput(dir, outputName(sample.r1.name, keepMatching), gzip);
  const w2 = await openOutput(dir, outputName(sample.r2.name, keepMatching), gzip);
  const it1 = records(sample.r1);
  const it2 = records(sample.r2);
  const SCRATCH = IO_BASE + IN_CAP + OUT_CAP - (1 << 20);
  let n = 0, dropped = 0, read = 0;
  const totalSize = sample.r1.size + sample.r2.size;

  try {
    for (;;) {
      if (cancelled) { await w1.abort(); await w2.abort(); return null; }
      const a = await it1.next(), b = await it2.next();
      if (a.done || b.done) {
        if (a.done !== b.done) throw new Error('R1 and R2 hold a different number of reads');
        break;
      }
      const r1 = a.value, r2 = b.value;
      read += r1.raw.length + r2.raw.length;

      const s1 = score(r1.seq, SCRATCH);
      const s2 = score(r2.seq, SCRATCH);
      const contaminated = Math.max(s1, s2) >= threshold;  // worst mate decides both
      const keep = keepMatching ? contaminated : !contaminated;

      n++;
      if (keep) { await w1.write(r1.raw.slice()); await w2.write(r2.raw.slice()); }
      else dropped++;

      if ((n & 0x3fff) === 0) {
        post('progress', {
          phase: 'filter', name: sample.name,
          fraction: totalSize ? read / totalSize : 0,
          detail: `${n.toLocaleString('en-US')} pairs`
        });
      }
    }
    await w1.close(); await w2.close();
    return { reads: n * 2, dropped: dropped * 2 };
  } catch (e) { try { await w1.abort(); await w2.abort(); } catch {} throw e; }
}

function score(seq, scratch) {
  view.set(seq, scratch);
  return Number(wasm.score_read(BigInt(scratch), 0n, BigInt(seq.length)));
}

// FASTQ records off a stream, never materialising the file.
async function* records(file) {
  const reader = decodedStream(file);
  let buf = new Uint8Array(0);
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    const merged = new Uint8Array(buf.length + value.length);
    merged.set(buf); merged.set(value, buf.length);
    buf = merged;
    let pos = 0;
    for (;;) {
      const ends = [];
      let j = pos;
      while (j < buf.length && ends.length < 4) {
        if (buf[j] === 10) ends.push(j);
        j++;
      }
      if (ends.length < 4) break;
      yield { raw: buf.subarray(pos, ends[3] + 1), seq: buf.subarray(ends[0] + 1, ends[1]) };
      pos = ends[3] + 1;
    }
    buf = buf.slice(pos);
  }
}

// Always yields an uncompressed base name; openOutput appends .gz when asked,
// so a gzipped input does not silently produce a plain file called .fastq.gz.
function outputName(name, keepMatching) {
  const suffix = keepMatching ? '.matched' : '.clean';
  const m = name.match(/^(.*?)(\.f(?:ast)?q)(\.gz|\.bgz)?$/i);
  return m ? `${m[1]}${suffix}${m[2]}` : `${name}${suffix}.fastq`;
}

// ————————————————————————————— dispatch —————————————————————————————

onmessage = async (e) => {
  const { cmd } = e.data;
  try {
    if (cmd === 'init') await init(e.data);
    else if (cmd === 'index') await buildIndex(e.data);
    else if (cmd === 'loadIndex') await loadNative(e.data);
    else if (cmd === 'fetchIndex') await fetchIndexFromZenodo(e.data);
    else if (cmd === 'filter') await filter(e.data);
    else if (cmd === 'cancel') cancelled = true;
  } catch (err) {
    post('error', { message: `${err.name}: ${err.message}` });
  }
};
