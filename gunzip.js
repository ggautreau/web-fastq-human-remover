// web-fastq-human-remover — multi-member gzip reader.
// Copyright (C) 2026 Guillaume Gautreau — MaIAGE (UR 1404), INRAE
// GPL-3.0 — see LICENSE.
//
// `DecompressionStream('gzip')` decodes exactly ONE gzip member, then throws
// "Junk found after end of compressed data". Sequencing FASTQ is very often
// multi-member: every bgzip/BGZF file is, and so is any concatenated .gz.
// A real 409 MB MGI file measured here held 240 members.
//
// Three paths, cheapest first:
//
//   BGZF   the header carries a 'BC' extra subfield with the block length, so
//          members are walked without scanning anything.
//   single try the whole file as one member — the common case, costing nothing
//          more than a plain DecompressionStream.
//   multi  scan ONCE for every member start. Scanning per member would be
//          quadratic: 240 members over 409 MB meant re-reading ~49 GB.
//
// File slices are read lazily, so a 3 GB input is never held in memory.

const GZIP_ID1 = 0x1f, GZIP_ID2 = 0x8b, DEFLATE = 8;
const FEXTRA = 4;

/** BGZF block length at `pos`, or null when this is not a BGZF block. */
async function bgzfBlockSize(file, pos) {
  const head = new Uint8Array(await file.slice(pos, pos + 18).arrayBuffer());
  if (head.length < 18) return null;
  if (head[0] !== GZIP_ID1 || head[1] !== GZIP_ID2 || head[2] !== DEFLATE) return null;
  if (!(head[3] & FEXTRA)) return null;
  const xlen = head[10] | (head[11] << 8);
  const extra = new Uint8Array(await file.slice(pos + 12, pos + 12 + xlen).arrayBuffer());
  let i = 0;
  while (i + 4 <= extra.length) {
    const slen = extra[i + 2] | (extra[i + 3] << 8);
    if (extra[i] === 0x42 && extra[i + 1] === 0x43 && slen === 2) {   // 'B','C'
      return (extra[i + 4] | (extra[i + 5] << 8)) + 1;                // BSIZE is size-1
    }
    i += 4 + slen;
  }
  return null;
}

/** Does a plausible gzip header start at `off`? */
function looksLikeHeader(buf, off) {
  if (off + 10 > buf.length) return false;
  if (buf[off] !== GZIP_ID1 || buf[off + 1] !== GZIP_ID2 || buf[off + 2] !== DEFLATE) return false;
  if (buf[off + 3] & 0xe0) return false;                 // reserved flag bits
  const xfl = buf[off + 8], os = buf[off + 9];
  if (xfl !== 0 && xfl !== 2 && xfl !== 4) return false;
  return os <= 13 || os === 255;
}

function concat(a, b) {
  const c = new Uint8Array(a.length + b.length);
  c.set(a); c.set(b, a.length);
  return c;
}

/** Passes bytes through untouched while reporting how many have gone by. */
function counting(stream, onBytes) {
  let n = 0;
  return stream.pipeThrough(new TransformStream({
    transform(chunk, ctrl) { n += chunk.length; onBytes(n); ctrl.enqueue(chunk); }
  }));
}

async function inflateMember(blob) {
  return new Uint8Array(await new Response(
    blob.stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
}

/**
 * Every member start, in one pass. A signature alone is not proof — 1f 8b 08
 * does occur inside deflate data — so the following header is validated, and
 * the caller still merges any offset whose slice will not inflate.
 */
async function memberStarts(file, onProgress) {
  const CHUNK = 4 << 20;
  const starts = [0];
  let pos = 0;
  let carry = new Uint8Array(0);
  while (pos < file.size) {
    const buf = new Uint8Array(await file.slice(pos, Math.min(pos + CHUNK, file.size)).arrayBuffer());
    const view = carry.length ? concat(carry, buf) : buf;
    const base = pos - carry.length;
    for (let i = 0; i + 10 <= view.length; i++) {
      if (view[i] !== GZIP_ID1) continue;                // cheap reject first
      if (looksLikeHeader(view, i)) {
        const at = base + i;
        if (at > starts[starts.length - 1]) starts.push(at);
      }
    }
    carry = view.slice(Math.max(0, view.length - 32));
    pos += buf.length;
    if (onProgress) onProgress(pos, file.size);
  }
  starts.push(file.size);
  return starts;
}

/**
 * Decompressed bytes of a gzip file, whatever its member layout.
 *
 * `onProgress(compressedBytesRead, fileSize)` reports position in the *input*.
 * Progress must be measured there, not on the output: comparing decompressed
 * bytes to a compressed file size sends the bar past 100 % by the compression
 * ratio.
 */
export function gunzipStream(file, onProgress) {
  return new ReadableStream({
    async start(controller) {
      try {
        // ——— BGZF: lengths live in the headers ———
        if (await bgzfBlockSize(file, 0) !== null) {
          let pos = 0;
          while (pos < file.size) {
            const bsize = await bgzfBlockSize(file, pos);
            if (bsize === null) break;                   // trailing bytes: stop cleanly
            const part = await inflateMember(file.slice(pos, pos + bsize));
            if (part.length) controller.enqueue(part);
            pos += bsize;
            if (onProgress) onProgress(pos, file.size);
          }
          controller.close();
          return;
        }

        // ——— everything else: scan first, THEN decode ———
        //
        // Trying a single DecompressionStream first was tempting, but it emits
        // bytes before it discovers the second member; restarting would then
        // duplicate them. So the member map is built before anything goes out.
        const starts = await memberStarts(file, null);

        // one member: stream it, no slicing, no extra cost
        if (starts.length === 2) {
          const src = onProgress ? counting(file.stream(), n => onProgress(n, file.size))
                                 : file.stream();
          const rd = src.pipeThrough(new DecompressionStream('gzip')).getReader();
          for (;;) {
            const { done, value } = await rd.read();
            if (done) break;
            if (value.length) controller.enqueue(value);
          }
          controller.close();
          return;
        }

        // ——— multi-member: walk the map ———
        let i = 0;
        while (i + 1 < starts.length) {
          let part = null, next = -1;
          for (let j = i + 1; j < starts.length; j++) {
            try {
              part = await inflateMember(file.slice(starts[i], starts[j]));
              next = j;
              break;
            } catch { /* false signature: extend to the following one */ }
          }
          if (part === null) throw new Error(`cannot inflate member at byte ${starts[i]}`);
          if (part.length) controller.enqueue(part);
          i = next;
          if (onProgress) onProgress(starts[i], file.size);
        }
        controller.close();
      } catch (e) {
        controller.error(e);
      }
    }
  });
}

/** True when the file starts with a gzip signature. */
export async function isGzip(file) {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return head.length === 2 && head[0] === GZIP_ID1 && head[1] === GZIP_ID2;
}
