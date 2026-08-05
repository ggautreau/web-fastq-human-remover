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

// pickle.js — minimal reader for the pickle protocol, enough to parse
// Cleanifier's `.info` sidecar without any Python round-trip.
//
// The file is a 3-tuple (filterinfo, optinfo, appinfo) of plain dicts holding
// str / int / tuple / numpy scalars. We implement only the opcodes those
// produce; anything else throws loudly rather than returning a half-decoded
// object, because silently wrong index parameters would mean silently wrong
// filtering.

const OP = {
  MARK: 0x28, STOP: 0x2e, EMPTY_TUPLE: 0x29, EMPTY_DICT: 0x7d, EMPTY_LIST: 0x5d,
  BININT: 0x4a, BININT1: 0x4b, BININT2: 0x4d, NONE: 0x4e,
  BINPUT: 0x71, LONG_BINPUT: 0x72, BINGET: 0x68, LONG_BINGET: 0x6a,
  SETITEM: 0x73, SETITEMS: 0x75, TUPLE: 0x74, APPEND: 0x61, APPENDS: 0x65,
  BINUNICODE: 0x58, SHORT_BINBYTES: 0x43, BINBYTES: 0x42, GLOBAL: 0x63, REDUCE: 0x52,
  PROTO: 0x80, NEWOBJ: 0x81, TUPLE1: 0x85, TUPLE2: 0x86, TUPLE3: 0x87,
  NEWTRUE: 0x88, NEWFALSE: 0x89, LONG1: 0x8a, LONG4: 0x8b,
  SHORT_BINUNICODE: 0x8c, BINUNICODE8: 0x8d, BINBYTES8: 0x8e,
  EMPTY_SET: 0x8f, FROZENSET: 0x91, STACK_GLOBAL: 0x93, MEMOIZE: 0x94, FRAME: 0x95,
};

// Marks the stack so TUPLE/SETITEMS know where the group started.
const MARKER = Symbol('mark');

export function loadPickle(bytes) {
  const u8 = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  const dec = new TextDecoder();
  let i = 0;
  const stack = [];
  const memo = new Map();

  const popMark = () => {
    const k = stack.lastIndexOf(MARKER);
    if (k < 0) throw new Error('pickle: unbalanced MARK');
    return stack.splice(k).slice(1);
  };
  const str = n => { const s = dec.decode(u8.subarray(i, i + n)); i += n; return s; };

  // Numpy scalars arrive as numpy.core.multiarray.scalar(dtype, raw_bytes).
  // We only need their numeric value, so we resolve them to a JS number or
  // BigInt and drop the dtype machinery entirely.
  const numpyScalar = (dtypeName, raw) => {
    const d = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    switch (dtypeName) {
      case 'u8': return d.getBigUint64(0, true);
      case 'i8': return d.getBigInt64(0, true);
      case 'u4': return d.getUint32(0, true);
      case 'i4': return d.getInt32(0, true);
      case 'u2': return d.getUint16(0, true);
      case 'i2': return d.getInt16(0, true);
      case 'u1': return d.getUint8(0);
      case 'i1': return d.getInt8(0);
      case 'f8': return d.getFloat64(0, true);
      case 'f4': return d.getFloat32(0, true);
      default: throw new Error(`pickle: unsupported numpy dtype '${dtypeName}'`);
    }
  };

  while (i < u8.length) {
    const op = u8[i++];
    switch (op) {
      case OP.PROTO: i += 1; break;
      case OP.FRAME: i += 8; break;
      case OP.STOP: return stack.pop();

      case OP.MARK: stack.push(MARKER); break;
      case OP.MEMOIZE: memo.set(memo.size, stack[stack.length - 1]); break;
      case OP.BINPUT: memo.set(u8[i++], stack[stack.length - 1]); break;
      case OP.LONG_BINPUT: memo.set(dv.getUint32(i, true), stack[stack.length - 1]); i += 4; break;
      case OP.BINGET: stack.push(memo.get(u8[i++])); break;
      case OP.LONG_BINGET: stack.push(memo.get(dv.getUint32(i, true))); i += 4; break;

      case OP.NONE: stack.push(null); break;
      case OP.NEWTRUE: stack.push(true); break;
      case OP.NEWFALSE: stack.push(false); break;

      case OP.BININT: stack.push(dv.getInt32(i, true)); i += 4; break;
      case OP.BININT1: stack.push(u8[i++]); break;
      case OP.BININT2: stack.push(dv.getUint16(i, true)); i += 2; break;
      case OP.LONG1: {                       // arbitrary-precision int, little-endian
        const n = u8[i++];
        let v = 0n;
        for (let k = n - 1; k >= 0; k--) v = (v << 8n) | BigInt(u8[i + k]);
        if (n && (u8[i + n - 1] & 0x80)) v -= 1n << BigInt(8 * n);   // two's complement
        i += n;
        stack.push(v);
        break;
      }

      case OP.SHORT_BINUNICODE: stack.push(str(u8[i++])); break;
      case OP.BINUNICODE: { const n = dv.getUint32(i, true); i += 4; stack.push(str(n)); break; }
      case OP.SHORT_BINBYTES: { const n = u8[i++]; stack.push(u8.slice(i, i + n)); i += n; break; }
      case OP.BINBYTES: { const n = dv.getUint32(i, true); i += 4; stack.push(u8.slice(i, i + n)); i += n; break; }

      case OP.EMPTY_TUPLE: stack.push([]); break;
      case OP.TUPLE1: { const a = stack.pop(); stack.push([a]); break; }
      case OP.TUPLE2: { const b = stack.pop(), a = stack.pop(); stack.push([a, b]); break; }
      case OP.TUPLE3: { const c = stack.pop(), b = stack.pop(), a = stack.pop(); stack.push([a, b, c]); break; }
      case OP.TUPLE: stack.push(popMark()); break;
      case OP.EMPTY_LIST: stack.push([]); break;
      case OP.APPEND: { const v = stack.pop(); stack[stack.length - 1].push(v); break; }
      case OP.APPENDS: { const items = popMark(); stack[stack.length - 1].push(...items); break; }

      case OP.EMPTY_DICT: stack.push({}); break;
      case OP.SETITEM: { const v = stack.pop(), k = stack.pop(); stack[stack.length - 1][k] = v; break; }
      case OP.SETITEMS: {
        const items = popMark();
        const d = stack[stack.length - 1];
        for (let k = 0; k < items.length; k += 2) d[items[k]] = items[k + 1];
        break;
      }

      case OP.GLOBAL: {                       // module \n name \n
        let s = '';
        while (u8[i] !== 0x0a) s += String.fromCharCode(u8[i++]);
        i++;
        let n = '';
        while (u8[i] !== 0x0a) n += String.fromCharCode(u8[i++]);
        i++;
        stack.push({ __global__: `${s}.${n}` });
        break;
      }
      case OP.STACK_GLOBAL: {
        const name = stack.pop(), mod = stack.pop();
        stack.push({ __global__: `${mod}.${name}` });
        break;
      }

      case OP.REDUCE: {
        const args = stack.pop();
        const fn = stack.pop();
        const g = fn && fn.__global__;
        if (g === 'numpy.core.multiarray.scalar' || g === 'numpy._core.multiarray.scalar') {
          const [dtype, raw] = args;
          const name = dtype && dtype.__dtype__ ? dtype.__dtype__ : 'u8';
          stack.push(numpyScalar(name, raw));
        } else if (g && /dtype$/.test(g)) {
          stack.push({ __dtype__: String(args[0]).replace(/^[<>|]/, '') });
        } else {
          // unknown callable: keep the arguments so the caller can inspect them
          stack.push({ __reduce__: g, args });
        }
        break;
      }

      case OP.NEWOBJ: { const args = stack.pop(); const cls = stack.pop(); stack.push({ __obj__: cls, args }); break; }

      // dtype objects carry state through BUILD; we ignore it, the name is enough
      case 0x62 /* BUILD */: { const st = stack.pop(); const o = stack[stack.length - 1];
        if (o && o.__dtype__ === undefined && Array.isArray(st)) o.__state__ = st;
        break; }

      default:
        throw new Error(`pickle: unsupported opcode 0x${op.toString(16)} at offset ${i - 1}`);
    }
  }
  throw new Error('pickle: reached end of data without STOP');
}

/** Reads Cleanifier's `.info` and returns the parameters the filter needs. */
export function readInfo(bytes) {
  const obj = loadPickle(bytes);
  if (!Array.isArray(obj) || obj.length < 3) throw new Error('.info: expected a 3-tuple');
  const [filterinfo, , appinfo] = obj;
  const num = v => (typeof v === 'bigint' ? Number(v) : v);
  const seed = (appinfo.mask || []).map(num);
  return {
    filtertype: filterinfo.filtertype,
    universe: num(filterinfo.universe),
    nsubfilters: num(filterinfo.nsubfilters),
    nslots: num(filterinfo.nslots),
    nwindows: num(filterinfo.nwindows_per_subfilter),
    windowsize: num(filterinfo.windowsize),
    fprBits: num(filterinfo.target_fpr),
    maxsteps: num(filterinfo.maxsteps),
    hash: {
      subfilter: filterinfo.subfilter_hashfunc_str,
      window: filterinfo.window_hashfunc_str,
      fingerprint: filterinfo.fingerprint_hashfunc_str,
      offset: filterinfo.offset_hashfunc_str,
    },
    k: num(appinfo.k),
    seed,
    span: seed.length ? Math.max(...seed) + 1 : num(appinfo.k),
    rcmode: appinfo.rcmode,
    checksum: filterinfo.checksum === undefined ? null : String(filterinfo.checksum),
  };
}
