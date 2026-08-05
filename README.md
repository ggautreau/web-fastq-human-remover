# web-fastq-human-remover

Split a FASTQ against a reference genome — **entirely in your browser**. No server, no upload,
no install.

Each sample yields **two files**: the reads that match the reference, and those that do not.
Which side is contamination and which is signal is your call, not the tool's — so it writes both
and lets you decide.

**→ [ggautreau.github.io/web-fastq-human-remover](https://ggautreau.github.io/web-fastq-human-remover/)**

Rust core compiled to WebAssembly (memory64): a Cuckoo filter over canonical spaced seeds. It reads
Cleanifier's published indexes directly, or builds one locally from a FASTA you supply.

> **Unofficial, independent port.** The method comes from
> [Cleanifier](https://gitlab.com/rahmannlab/cleanifier), released under the MIT licence by Jens
> Zentgraf, Johanna Elena Schmitz and Sven Rahmann (Algorithmic Bioinformatics, Saarland
> University). This port has **no connection to them and carries no endorsement from them**.
> For publishable work, run Cleanifier natively.

## Please cite

> Zentgraf, J., Schmitz, J.E., Rahmann, S.: *Cleanifier: Contamination removal from microbial
> sequences using spaced seeds of a human pangenome index.* Bioinformatics (2025).
> https://doi.org/10.1093/bioinformatics/btaf632

If you use their pre-built index:

> Schmitz, J. E., Zentgraf, J., & Rahmann, S. (2025). *Human index for Cleanifier* (0.1.0). Zenodo.
> https://doi.org/10.5281/zenodo.15639519

## Two ways to get a reference

**Prebuilt index** — the published indexes download straight into the tab from Zenodo and are
untarred on the fly. Nothing hits the disk unless you ask, nothing goes through a proxy.

| Index | Size | Note |
|---|---|---|
| Ribosomal RNA | 282 MB | easiest to try; seed k=25, w=31 |
| Human, probabilistic | 6.36 GB | needs ~7 GB of free RAM; seed k=29, w=33 |

Tick *also save to disk* and the extracted `.filter`/`.info` are written to a folder you pick, as
the same bytes stream past — no second transfer. Next session, drop those two files instead.

Zenodo also hosts an *exact* human index; it is a 3-way bucketed hash table (`.hash`) rather than a
cuckoo filter, which this build does not read. It covers the same genome as the probabilistic one.

**Reference genome** — drop a FASTA and the index is built locally, sized on it: a few MB for PhiX
or a bacterium, gigabytes only for a whole human genome.

## Accuracy

Against the published human index on **HG002** (GIAB, the sample from the paper's supplement),
100 000 Illumina reads of 125 bp, threshold 0.5:

| | Reads matching |
|---|---|
| Cleanifier `--sensitive` | 93 642 (93.64 %) |
| **This port** | **93 637 (93.64 %)** |

Five reads out of 100 000 differ (0.005 %), from N handling. Seed-level agreement is exact:
15 291/18 600 found by both.

The native index reader agrees with Cleanifier on **60 000/60 000 keys**, false positives included.

Note Cleanifier's default mode is `sampling_issh`, which subsamples seeds; `--sensitive` runs the
exhaustive `classify_issh`. This port implements the latter, so compare against `--sensitive`.

## Test data

Ready to use in [`testdata/samples/`](testdata/samples/) — four gzipped paired-end samples, 1.3 MB
total. Drop a pair into the app directly. `testdata/gen_testsets.py` regenerates them. Measured:

| Sample | matched, human index | matched, rRNA index |
|---|---|---|
| `with_human` | **94.1 %** | 0.3 % |
| `without_human` | **0.0 %** | 50.0 % |
| `with_rrna` | 50.0 % | **100.0 %** |
| `without_rrna` | **0.0 %** | **0.0 %** |

The off-diagonal figures are correct, not errors: `with_rrna` is half human 18S, which genuinely
belongs to the human genome; `without_human` is half *E. coli* 16S, which is genuinely ribosomal.

## Running locally

Any static server works — the page needs **no special headers**:

```bash
python3 serve.py     # http://127.0.0.1:8767/
```

## Build

```bash
./build.sh
```

Needs nightly: `wasm64-unknown-unknown` is Tier 3 with no precompiled `core`, and
`core::arch::wasm64` is still unstable (rust-lang#90599).

## Why it works on GitHub Pages

`SharedArrayBuffer` requires COOP/COEP headers, which GitHub Pages cannot set. This build sidesteps
that entirely: the wasm memory is **not shared** (one worker owns the index), so no isolation is
needed. Verified with `crossOriginIsolated: false` and `SharedArrayBuffer` undefined.

The trade-off is no multi-threading. Throughput is ~14 000 reads/s single-threaded against the
human index, which is dominated by random-access TLB misses rather than by the lack of threads.

## Implementation notes

- **The encoding is A=0, C=1, T=2, G=3**, not the usual A/C/G/T — it comes from
  `(ascii >> 1) & 7`, which also maps N to 7. Complementing is XOR 2, not `3 - b`. Get it wrong and
  every seed of a human read comes back absent from the human index, silently.
- **The decision is base coverage**, not a fraction of seeds: each matching seed marks its
  significant positions in a bitmap, and the read goes when covered bases exceed
  `threshold × length`.
- **`compile_get_bucket` has three branches** by bucket count — power of two, even, odd. The offset
  hash uses the *even* one.
- **Never hardcode the wasm heap base**: `wasm-ld` lays stack and static data at the bottom of
  linear memory, and `__heap_base` is the only reliable answer.
- Everything runs in a worker; FASTQ files stream through and are never materialised.
- **`DecompressionStream('gzip')` decodes only the FIRST gzip member** and then throws *"Junk found
  after end of compressed data"*. Sequencing FASTQ is very often multi-member — every bgzip/BGZF
  file is, and so is any concatenated `.gz`. A real 409 MB MGI file tested here held **240 members**.
  `gunzip.js` walks them: BGZF block lengths come free from the `BC` extra subfield, plain
  concatenations get one scan for validated gzip headers.
- **The member map is built before any byte is emitted.** Trying a single stream first and falling
  back on failure looks cheaper, but `DecompressionStream` emits data *before* it discovers the
  second member, so restarting would silently duplicate it.
- **Scan once, not per member.** Re-scanning from each member start is quadratic: 240 members over
  409 MB meant re-reading ~49 GB.

## Limitations

- Reads longer than 32 MiB are not handled (irrelevant for Illumina).
- Multi-member gzip is slower than a single member (~210 MB/s): 76 MB/s on 4 MB members,
  ~30 MB/s on 64 KB BGZF blocks. The cost is the per-member stream setup.
- Cleanifier's exact `.hash` index is not supported.
- Build depends on nightly and one unstable feature (`simd_wasm64`).
- The locally-built index path is validated against synthetic ground truth, not the Cleanifier binary.

## Licence

**GNU General Public License v3.0** — see [LICENSE](LICENSE).

The method derives from [Cleanifier](https://gitlab.com/rahmannlab/cleanifier), released under the
**MIT licence**, copyright 2019-2025 Jens Zentgraf, Johanna Elena Schmitz and Sven Rahmann
(Saarland University). The MIT licence permits this relicensing under GPL-3.0; its notice is
retained in [NOTICE](NOTICE) as it requires.

---

Created using Claude Opus 5.0, prompted by G. Gautreau — MaIAGE unit (UR 1404, Mathématiques et
Informatique Appliquées du Génome à l'Environnement), INRAE.
