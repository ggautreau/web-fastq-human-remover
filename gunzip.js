// web-fastq-human-remover — multi-member gzip reader.
// Copyright (C) 2026 Guillaume Gautreau — MaIAGE (UR 1404), INRAE
// GPL-3.0 — see LICENSE.
//
// `DecompressionStream('gzip')` decodes exactly ONE gzip member and then throws
// "Junk found after end of compressed data". Sequencing FASTQ files are very
// often several members concatenated — every bgzip/BGZF file is, and so is any
// `cat a.gz b.gz`. Feeding those straight to DecompressionStream fails on the
// second member, which is what happens with real MGI/Illumina data.
//
// This module walks the members and decodes each one on its own. Two cases:
//
//   BGZF  — the header carries an extra 'BC' subfield holding BSIZE, so the
//           block length is known without decoding anything. Free and exact.
//   plain — no length anywhere, so the next member has to be found by scanning
//           for a gzip signature and validating the header that follows.
//
// Slices of a File are read lazily, so a 3 GB file is never held in memory.

const GZIP_ID1 = 0x1f, GZIP_ID2 = 0x8b, DEFLATE = 8;
const FEXTRA = 4;

/** BGZF block size at `pos`, or null when this is not a BGZF block. */
async function bgzfBlockSize(file, pos) {
  const head = new Uint8Array(await file.slice(pos, pos + 18).arrayBuffer());
  if (head.length < 18) return null;
  if (head[0] !== GZIP_ID1 || head[1] !== GZIP_ID2 || head[2] !== DEFLATE) return null;
  if (!(head[3] & FEXTRA)) return null;
  const xlen = head[10] | (head[11] << 8);
  const extra = new Uint8Array(await file.slice(pos + 12, pos + 12 + xlen).arrayBuffer());
  let i = 0;
  while (i + 4 <= extra.length) {
    const si1 = extra[i], si2 = extra[i + 1];
    const slen = extra[i + 2] | (extra[i + 3] << 8);
    if (si1 === 0x42 && si2 === 0x43 && slen === 2) {          // 'B','C'
      const bsize = extra[i + 4] | (extra[i + 5] << 8);
      return bsize + 1;                                        // BSIZE is size-1
    }
    i += 4 + slen;
  }
  return null;
}

/** Does a plausible gzip header start at `off` inside `buf`? */
function looksLikeHeader(buf, off) {
  if (off + 10 > buf.length) return false;
  if (buf[off] !== GZIP_ID1 || buf[off + 1] !== GZIP_ID2 || buf[off + 2] !== DEFLATE) return false;
  const flg = buf[off + 3];
  if (flg & 0xe0) return false;              // reserved bits must be zero
  const xfl = buf[off + 8], os = buf[off + 9];
  if (xfl !== 0 && xfl !== 2 && xfl !== 4) return false;
  return os <= 13 || os === 255;
}

/**
 * Offset of the next member after `from`, or file.size if this is the last one.
 * Only used for plain concatenated gzip; BGZF never gets here.
 */
async function nextMemberStart(file, from) {
  const CHUNK = 1 << 20;
  let pos = from;
  let carry = new Uint8Array(0);
  while (pos < file.size) {
    const buf = new Uint8Array(await file.slice(pos, Math.min(pos + CHUNK, file.size)).arrayBuffer());
    const view = carry.length ? concat(carry, buf) : buf;
    const base = pos - carry.length;
    for (let i = 0; i + 10 <= view.length; i++) {
      if (looksLikeHeader(view, i)) {
        const at = base + i;
        if (at > from) return at;
      }
    }
    carry = view.subarray(Math.max(0, view.length - 32));
    pos += buf.length;
  }
  return file.size;
}

function concat(a, b) {
  const c = new Uint8Array(a.length + b.length);
  c.set(a); c.set(b, a.length);
  return c;
}

async function inflateMember(blob) {
  return new Uint8Array(await new Response(
    blob.stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
}

// Decoding in parallel batches was tried and measured: no gain (29 vs 32 MB/s),
// so members are handled one at a time and the code stays simple. The cost is
// the per-member DecompressionStream setup, not the serialisation.

/**
 * Decompressed bytes of a gzip file, member by member.
 *
 * Single-member files — the common case — take the direct path and cost nothing
 * more than a plain DecompressionStream.
 */
export function gunzipStream(file) {
  return new ReadableStream({
    async start(controller) {
      try {
        // fast path: one member, stream it straight through
        const firstSize = await bgzfBlockSize(file, 0);
        if (firstSize === null) {
          const only = await nextMemberStart(file, 0);
          if (only >= file.size) {
            const rd = file.stream().pipeThrough(new DecompressionStream('gzip')).getReader();
            for (;;) {
              const { done, value } = await rd.read();
              if (done) break;
              if (value.length) controller.enqueue(value);
            }
            controller.close();
            return;
          }
        }

        let pos = 0;
        while (pos < file.size) {
          const bsize = await bgzfBlockSize(file, pos);
          const end = bsize !== null ? pos + bsize : await nextMemberStart(file, pos);
          if (end <= pos) throw new Error(`stalled at byte ${pos} of ${file.name}`);
          const part = await inflateMember(file.slice(pos, end));
          if (part.length) controller.enqueue(part);
          pos = end;
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
