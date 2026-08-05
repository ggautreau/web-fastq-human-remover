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

// app.js — UI only. No heavy work on this thread: parsing 1 GB here freezes
// the tab for about 8 seconds, so everything goes to the worker.

const $ = s => document.querySelector(s);
const fmt = n => n.toLocaleString('en-US');
const human = n => {
  const u = ['B', 'KB', 'MB', 'GB'];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
  return `${n.toFixed(n < 10 && i > 0 ? 1 : 0)} ${u[i]}`;
};

// Cleanifier's published seed for its human index (cleanifier_index.py, k=29, w=33).
// check_mask() enforces symmetry, which is what makes the reverse-complement
// shortcut in the Rust core valid.
const SEEDS = {
  'cleanifier': { mask: '######_#######_###_#######_######', label: "Cleanifier spaced seed (k=29, w=33)" },
  'k25':        { mask: '#'.repeat(25), label: 'Contiguous k-mer, k=25' },
  'k31':        { mask: '#'.repeat(31), label: 'Contiguous k-mer, k=31' },
};

const state = {
  source: 'ref',            // 'ref' = build from FASTA, 'index' = load a prebuilt one
  ref: [], idxFiles: [], samples: [], dir: null,
  indexReady: false, busy: false, worker: null, t0: 0
};

// ————————————————————————————— log —————————————————————————————

function log(msg, cls = '') {
  const el = $('#log');
  const t = new Date().toLocaleTimeString('en-US', { hour12: false });
  el.innerHTML += `\n<span class="${cls}">[${t}] ${msg}</span>`;
  el.scrollTop = el.scrollHeight;
}

// ————————————————————————————— theme —————————————————————————————

const theme0 = localStorage.getItem('fqc-theme')
  || (matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark');
document.documentElement.setAttribute('data-theme', theme0);
$('#theme').onclick = () => {
  const t = document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
  document.documentElement.setAttribute('data-theme', t);
  localStorage.setItem('fqc-theme', t);
};

// ————————————————————————— drop zones —————————————————————————

// ————————————————————— reference source: prebuilt index vs FASTA —————————————————————

function selectSource(src) {
  state.source = src;
  state.indexReady = false;
  $('#srcIndex').setAttribute('aria-pressed', String(src === 'index'));
  $('#srcRef').setAttribute('aria-pressed', String(src === 'ref'));
  $('#paneIndex').hidden = src !== 'index';
  $('#idxOk').hidden = true;
  $('#paneRef').hidden = src !== 'ref';
  $('#indexState').textContent = 'No index in memory.';
  $('#idxState').textContent = 'No index loaded.';
  refreshRun();
}
$('#srcIndex').onclick = () => selectSource('index');
$('#srcRef').onclick = () => selectSource('ref');

function wireDrop(zone, input, onFiles) {
  zone.onclick = () => input.click();
  input.onchange = () => { onFiles([...input.files]); input.value = ''; };
  ['dragenter', 'dragover'].forEach(e => zone.addEventListener(e, ev => {
    ev.preventDefault(); zone.classList.add('over');
  }));
  ['dragleave', 'drop'].forEach(e => zone.addEventListener(e, ev => {
    ev.preventDefault(); zone.classList.remove('over');
  }));
  zone.addEventListener('drop', ev => {
    const fs = [...(ev.dataTransfer?.files || [])];
    if (fs.length) onFiles(fs);
  });
}

// ————————————————————— reference and index sizing —————————————————————

// Rough base count: ~1 byte per base uncompressed, DNA FASTA compressing about
// 3.5x. Deliberately generous — an oversized index only costs memory, an
// undersized one silently drops k-mers.
function estimatedBases(files) {
  return files.reduce((s, f) => s + f.size * (/\.(gz|bgz)$/i.test(f.name) ? 3.5 : 1), 0);
}

function log2BucketsFor(bases) {
  const slots = bases / 0.95;              // target 95% load
  return Math.min(31, Math.max(10, Math.ceil(Math.log2(Math.max(1024, slots / 4)))));
}

function refreshRef() {
  $('#chipsRef').innerHTML = state.ref.map((f, i) =>
    `<span class="chip">${f.name} <span style="color:var(--ink-muted)">${human(f.size)}</span>
     <span class="x" data-i="${i}">✕</span></span>`).join('');
  $('#chipsRef').querySelectorAll('.x').forEach(x => x.onclick = e => {
    state.ref.splice(+e.target.dataset.i, 1); refreshRef();
  });
  $('#btnIndex').disabled = state.ref.length === 0 || state.busy;

  if (!state.ref.length) { $('#noteIndex').textContent = ''; refreshSizes(0); return; }
  const bases = estimatedBases(state.ref);
  refreshSizes(log2BucketsFor(bases));
  $('#noteIndex').textContent = `≈ ${fmt(Math.round(bases / 1e6))} Mb of reference`;
}

function refreshSizes(auto) {
  const sel = $('#size');
  const cur = sel.value;
  const opts = [`<option value="auto">automatic${auto ? ` (2^${auto})` : ''}</option>`];
  for (let lb = 16; lb <= 31; lb++) {
    opts.push(`<option value="${lb}">2^${lb} · ${human((2 ** lb) * 8)} · ~${fmt((2 ** lb) * 4 / 1e6)} M k-mers</option>`);
  }
  sel.innerHTML = opts.join('');
  sel.value = [...sel.options].some(o => o.value === cur) ? cur : 'auto';
}

wireDrop($('#dropRef'), $('#fileRef'), fs => {
  const keep = fs.filter(f => /\.(fa|fasta|fna)(\.gz|\.bgz)?$/i.test(f.name));
  if (keep.length < fs.length) log(`${fs.length - keep.length} file(s) ignored: unrecognised extension`, 'err');
  state.ref.push(...keep);
  refreshRef();
});

// ————————————————————— prebuilt index files —————————————————————

function refreshIdx() {
  $('#chipsIdx').innerHTML = state.idxFiles.map((f, i) =>
    `<span class="chip">${f.name} <span style="color:var(--ink-muted)">${human(f.size)}</span>
     <span class="x" data-i="${i}">✕</span></span>`).join('');
  $('#chipsIdx').querySelectorAll('.x').forEach(x => x.onclick = e => {
    state.idxFiles.splice(+e.target.dataset.i, 1); refreshIdx();
  });
  const hasFilter = state.idxFiles.some(f => /\.filter$/i.test(f.name));
  const hasInfo = state.idxFiles.some(f => /\.info$/i.test(f.name));
  $('#btnLoadIdx').disabled = !(hasFilter && hasInfo) || state.busy;
  if (state.idxFiles.length && !(hasFilter && hasInfo)) {
    $('#idxState').textContent = hasFilter ? 'missing the .info file' : 'missing the .filter file';
  }
}

wireDrop($('#dropIdx'), $('#fileIdx'), fs => {
  const keep = fs.filter(f => /\.(filter|info)$/i.test(f.name));
  if (keep.length < fs.length) log(`${fs.length - keep.length} file(s) ignored: expected .filter and .info`, 'err');
  state.idxFiles = [...state.idxFiles.filter(f => !keep.some(k => k.name === f.name)), ...keep];
  refreshIdx();
});

// Zenodo record holding all three published indexes
const ZENODO_RECORD = '19913374';

document.querySelectorAll('.dlbtn').forEach(b => b.onclick = async () => {
  const archive = b.dataset.archive;

  // The directory picker must be called straight from the click, before any
  // await, or the browser no longer considers it user-initiated.
  let saveDir = null;
  if ($('#saveIdx').checked) {
    if (!window.showDirectoryPicker) {
      log('this browser cannot save to disk (no File System Access) — downloading without saving', 'err');
    } else {
      try { saveDir = await window.showDirectoryPicker({ mode: 'readwrite' }); }
      catch (e) {
        if (e.name === 'AbortError') { log('save cancelled — download not started'); return; }
        log(`save folder: ${e.message}`, 'err'); return;
      }
      log(`the extracted .filter and .info will be written to ${saveDir.name}/`);
    }
  }

  busy(true);
  state.dlStart = performance.now();
  $('#dlProg').hidden = false;
  $('#dlName').textContent = `${archive} · ${b.dataset.size}`;
  $('#dlDetail').textContent = 'contacting Zenodo…';
  $('#dlEta').textContent = '';
  $('#dlBar').style.width = '0%';
  $('#dlBar').classList.add('waiting');
  progress(0, 'contacting Zenodo…');
  log(`downloading ${archive} (${b.dataset.size}) from Zenodo — streamed straight into memory`);
  if (+(b.dataset.bytes || 0) > 1e9) log('large file: Zenodo can take ~20 s before the first byte arrives');
  if (state.worker) { state.worker.terminate(); state.worker = null; }
  worker().postMessage({ cmd: 'fetchIndex', wasmUrl: 'fqclean.wasm',
                         record: ZENODO_RECORD, archive, saveDir });
});

// Cancelling a download stops the stream in the worker and tears the worker
// down, so the half-filled linear memory is released rather than left behind.
$('#btnCancelDl').onclick = () => {
  log('cancelling download…');
  $('#dlDetail').textContent = 'cancelling…';
  if (state.worker) state.worker.postMessage({ cmd: 'cancel' });
  setTimeout(() => {
    if (state.busy && state.worker) {
      state.worker.terminate(); state.worker = null;
      log('download cancelled — memory released', 'err');
      hideDl(); progress(0, 'idle'); idle();
    }
  }, 1500);
};

function hideDl() { $('#dlProg').hidden = true; }

$('#btnLoadIdx').onclick = () => {
  const filter = state.idxFiles.find(f => /\.filter$/i.test(f.name));
  const info = state.idxFiles.find(f => /\.info$/i.test(f.name));
  busy(true);
  progress(0, 'reading index…');
  log(`loading prebuilt index: ${filter.name} (${human(filter.size)})`);
  if (state.worker) { state.worker.terminate(); state.worker = null; }
  worker().postMessage({ cmd: 'loadIndex', wasmUrl: 'fqclean.wasm', filter, info });
};

// ————————————————————— samples and mate pairing —————————————————————

// Groups R1/R2. The _1/_2 and _R1/_R2 patterns cover the vast majority of
// sequencer output; anything else is treated as single-end.
function pairUp(files) {
  const left = [...files];
  const groups = [];
  const key = n => n.replace(/\.(f(ast)?q)(\.gz|\.bgz)?$/i, '').replace(/_(R?)[12](_001)?$/i, '');
  const mate = n => {
    const m = n.replace(/\.(f(ast)?q)(\.gz|\.bgz)?$/i, '').match(/_R?([12])(_001)?$/i);
    return m ? +m[1] : 0;
  };
  while (left.length) {
    const f = left.shift();
    const n = mate(f.name);
    if (n === 0) { groups.push({ name: f.name, r1: f, r2: null }); continue; }
    const k = key(f.name);
    const j = left.findIndex(g => key(g.name) === k && mate(g.name) === (n === 1 ? 2 : 1));
    if (j >= 0) {
      const other = left.splice(j, 1)[0];
      const [r1, r2] = n === 1 ? [f, other] : [other, f];
      groups.push({ name: k, r1, r2 });
    } else groups.push({ name: f.name, r1: f, r2: null });
  }
  return groups;
}

function refreshSamples() {
  $('#chipsFq').innerHTML = state.samples.map((s, i) => {
    const size = s.r1.size + (s.r2?.size || 0);
    return `<span class="chip">${s.name}
      ${s.r2 ? '<span class="paired">R1+R2</span>' : ''}
      <span style="color:var(--ink-muted)">${human(size)}</span>
      <span class="x" data-i="${i}">✕</span></span>`;
  }).join('');
  $('#chipsFq').querySelectorAll('.x').forEach(x => x.onclick = e => {
    state.samples.splice(+e.target.dataset.i, 1); refreshSamples();
  });
  refreshRun();
}

wireDrop($('#dropFq'), $('#fileFq'), fs => {
  const fq = fs.filter(f => /\.(fastq|fq)(\.gz|\.bgz)?$/i.test(f.name));
  if (fq.length < fs.length) log(`${fs.length - fq.length} file(s) ignored: unrecognised extension`, 'err');
  state.samples = pairUp([...state.samples.flatMap(s => [s.r1, s.r2].filter(Boolean)), ...fq]);
  refreshSamples();
});

// ————————————————————— threshold —————————————————————

function refreshThreshold() {
  const pct = +$('#threshold').value;
  $('#thresholdVal').textContent = pct + ' %';
  const seed = SEEDS[$('#seed').value];
  const weight = [...seed.mask].filter(c => c === '#').length;
  const one = 100 / Math.max(1, 150 - seed.mask.length + 1);   // one hit on a 150 bp read
  $('#thresholdNote').innerHTML =
    `Fraction of a read's seeds that must be present in the reference for it to be removed. ` +
    `On a 150 bp read this seed yields ${150 - seed.mask.length + 1} positions, so anything below ` +
    `<b>${one.toFixed(1)} %</b> means a single shared seed is enough. Weight ${weight}, span ${seed.mask.length}.`;
}
$('#threshold').oninput = refreshThreshold;
$('#seed').onchange = refreshThreshold;

// ————————————————————— output directory —————————————————————

$('#btnDir').onclick = async () => {
  if (!window.showDirectoryPicker) {
    log('This browser does not expose the File System Access API — use Chrome or Edge.', 'err');
    return;
  }
  try {
    state.dir = await window.showDirectoryPicker({ mode: 'readwrite' });
    $('#dirState').textContent = `→ ${state.dir.name}`;
    log(`output directory: ${state.dir.name}`, 'ok');
    refreshRun();
  } catch (e) {
    if (e.name !== 'AbortError') log(`directory: ${e.message}`, 'err');
  }
};

function refreshRun() {
  $('#btnRun').disabled = !(state.indexReady && state.samples.length && state.dir && !state.busy);
}

// ————————————————————— worker —————————————————————

function worker() {
  if (state.worker) return state.worker;
  const w = new Worker('worker.js', { type: 'module' });
  w.onmessage = e => handle(e.data);
  w.onerror = e => { log(`worker: ${e.message}`, 'err'); idle(); };
  state.worker = w;
  return w;
}

function progress(fraction, text) {
  $('#prog').style.width = Math.round(Math.min(1, Math.max(0, fraction)) * 100) + '%';
  $('#progTxt').textContent = text;
  $('#progPct').textContent = fraction > 0 ? Math.round(fraction * 100) + ' %' : '';
}

function handle(m) {
  switch (m.type) {
    case 'ready':
      log(`index allocated: ${m.index} · ${fmt(m.buckets)} buckets · capacity ~${fmt(m.kmerCapacity)} seeds (weight ${m.weight}, span ${m.span})`);
      worker().postMessage({ cmd: 'index', files: state.ref });
      break;

    case 'progress':
      progress(m.fraction, m.phase === 'index' ? `indexing — ${m.detail}` : `${m.name} — ${m.detail}`);
      if (!$('#dlProg').hidden && m.phase === 'index') {
        // Zenodo takes ~20 s to start sending a multi-GB file; keep the bar
        // animated until the first byte lands so it does not look frozen
        if (m.fraction > 0) $('#dlBar').classList.remove('waiting');
        const pct = Math.min(1, m.fraction) * 100;
        $('#dlBar').style.width = pct.toFixed(2) + '%';
        $('#dlDetail').textContent = m.detail;
        const el = (performance.now() - state.dlStart) / 1000;
        if (m.fraction > 0.01 && el > 1) {
          const left = el * (1 - m.fraction) / m.fraction;
          const mn = Math.floor(left / 60), s = Math.round(left % 60);
          const pctTxt = m.fraction < 0.1 ? (m.fraction * 100).toFixed(1) : Math.round(m.fraction * 100);
          $('#dlEta').textContent =
            `${pctTxt} % · ${mn > 0 ? mn + ' min ' : ''}${s} s left`;
        }
      }
      break;

    case 'nativeDone': {
      state.indexReady = true;
      hideDl();
      // retour visuel explicite : sans lui, un index chargé ne se distingue pas
      // d'un index absent, et le bouton de dépôt manuel prête à confusion
      $('#idxOk').hidden = false;
      $('#idxOkTitle').textContent = `Index ready — ${m.size}`;
      $('#idxOkDetail').textContent =
        `weight ${m.weight}, span ${m.span}, ${m.subfilters} subfilters. `
        + (m.saved && m.saved.length
            ? `Saved to disk (${m.saved.join(', ')}) — next time, drop those two files instead. `
            : '')
        + `Next: drop your FASTQ files in step 2, pick an output folder in step 3, then Run.`;
      if (m.saved && m.saved.length) log(`saved to disk: ${m.saved.join(', ')}`, 'ok');
      progress(1, 'index ready');
      log(`native index loaded: ${m.size} in ${m.seconds} s (${m.rate}) · `
        + `weight ${m.weight}, span ${m.span}, ${m.subfilters} subfilters, ${m.fprBits}-bit fingerprints, rcmode ${m.rcmode}`, 'ok');
      $('#idxState').textContent = `ready · ${m.size} · k=${m.weight}, w=${m.span}`;
      idle();
      break;
    }

    case 'indexDone':
      state.indexReady = true;
      progress(1, 'index ready');
      log(`index built in ${m.seconds} s · ${fmt(m.kmers)} seeds inserted · load ${m.loadPct} %`, 'ok');
      if (m.failed > 0) log(`${fmt(m.failed)} insertion(s) dropped — index too small`, 'err');
      if (m.warning) log(m.warning, 'err');
      $('#indexState').textContent = `ready · ${fmt(m.kmers)} seeds · load ${m.loadPct} %`;
      $('#idxOk').hidden = false;
      $('#idxOkTitle').textContent = `Index ready — ${fmt(m.kmers)} seeds`;
      $('#idxOkDetail').textContent =
        `load ${m.loadPct} %. Next: drop your FASTQ files in step 2, pick an output folder in step 3, then Run.`;
      idle();
      break;

    case 'sampleDone': {
      const pct = m.reads ? (m.dropped / m.reads * 100) : 0;
      const rate = m.seconds > 0 ? Math.round(m.reads / m.seconds) : 0;
      $('#cardRes').style.display = '';
      const tr = document.createElement('tr');
      tr.innerHTML = `<td>${m.name}${m.paired ? ' <span class="paired">PE</span>' : ''}</td>
        <td>${fmt(m.reads)}</td><td>${fmt(m.dropped)}</td>
        <td class="pct">${pct.toFixed(2)} %</td><td>${m.seconds} s</td><td>${fmt(rate)} reads/s</td>`;
      $('#tblRes tbody').appendChild(tr);
      log(`${m.name}: ${fmt(m.dropped)}/${fmt(m.reads)} reads removed (${pct.toFixed(2)} %) in ${m.seconds} s`, 'ok');
      break;
    }

    case 'done':
      progress(1, 'finished');
      log(`finished in ${((performance.now() - state.t0) / 1000).toFixed(1)} s`, 'ok');
      idle();
      break;

    case 'cancelled':
      hideDl();
      progress(0, 'cancelled');
      log($('#dlProg').hidden ? 'cancelled — partial output files were discarded'
                              : 'download cancelled', 'err');
      idle();
      break;

    case 'error':
      hideDl();
      log(m.message, 'err');
      idle();
      break;
  }
}

function busy(v) {
  state.busy = v;
  $('#btnIndex').disabled = v || !state.ref.length;
  document.querySelectorAll('.dlbtn').forEach(b => b.disabled = v);
  $('#btnCancel').disabled = !v;
  $('#btnDir').disabled = v;
  refreshRun();
}
const idle = () => busy(false);

// ————————————————————— actions —————————————————————

$('#btnIndex').onclick = () => {
  const seed = SEEDS[$('#seed').value];
  const seedBits = [...seed.mask].reduce((a, c, i) => c === '#' ? a | (1n << BigInt(i)) : a, 0n);
  const sel = $('#size').value;
  const lb = sel === 'auto' ? log2BucketsFor(estimatedBases(state.ref)) : +sel;
  const bytes = (2 ** lb) * 8;

  state.indexReady = false;
  busy(true);
  progress(0, 'allocating index…');
  log(`building: seed ${seed.mask} · index 2^${lb} buckets (${human(bytes)})`);
  if (bytes > 4 * 1024 ** 3) log('index above 4 GB: needs memory64 support and enough free RAM', 'err');

  // recreate the worker so we always start from a clean linear memory
  if (state.worker) { state.worker.terminate(); state.worker = null; }
  worker().postMessage({
    cmd: 'init', wasmUrl: 'fqclean.wasm', log2Buckets: lb,
    seedBits: seedBits.toString(), span: seed.mask.length
  });
};

$('#btnRun').onclick = () => {
  busy(true);
  state.t0 = performance.now();
  $('#tblRes tbody').innerHTML = '';
  const threshold = Math.round(+$('#threshold').value * 10);   // % -> permille
  const keepMatching = $('#mode').value === 'extract';
  log(`filtering ${state.samples.length} sample(s) · threshold ${$('#threshold').value} %` +
      (keepMatching ? ' · extract mode' : '') + ($('#gzipOut').checked ? ' · gzip output' : ''));
  worker().postMessage({
    cmd: 'filter', samples: state.samples,
    thresholdPermille: threshold, keepMatching, outputDir: state.dir,
    gzipOutput: $('#gzipOut').checked
  });
};

$('#btnCancel').onclick = () => { worker().postMessage({ cmd: 'cancel' }); log('cancelling…'); };

$('#btnClear').onclick = () => {
  if (state.worker) { state.worker.terminate(); state.worker = null; }
  Object.assign(state, { ref: [], samples: [], indexReady: false, busy: false });
  refreshRef(); refreshSamples();
  $('#indexState').textContent = 'No index in memory.';
  $('#idxOk').hidden = true;
  $('#tblRes tbody').innerHTML = '';
  $('#cardRes').style.display = 'none';
  progress(0, 'idle');
  log('reset');
};

// ————————————————————— environment checks —————————————————————

(function check() {
  const missing = [];
  if (typeof DecompressionStream === 'undefined') missing.push('DecompressionStream');
  if (!window.showDirectoryPicker) missing.push('File System Access');
  try { new WebAssembly.Memory({ initial: 1n, maximum: 262144n, address: 'i64' }); }
  catch { missing.push('WebAssembly memory64'); }

  if (missing.length) {
    log(`unavailable: ${missing.join(', ')}`, 'err');
    log('Needs a recent Chrome or Edge (memory64 and File System Access).', 'err');
  } else {
    log('environment supported — drop a reference FASTA to begin.', 'ok');
  }
  refreshRef(); refreshIdx(); refreshSamples(); refreshThreshold();
  selectSource('ref');
})();
