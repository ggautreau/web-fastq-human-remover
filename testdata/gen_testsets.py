#!/usr/bin/env python3
"""Builds four paired-end test samples with known ground truth.

Each sample is deliberately unambiguous for one of the two published indexes,
so a wrong result is obvious rather than plausible:

  with_human      100 % real human reads (GIAB HG002)
  without_human   random DNA + E. coli 16S — nothing human
  with_rrna       E. coli 16S + human 18S — all ribosomal
  without_rrna    random DNA only — no ribosomal RNA anywhere

Random DNA is used for the negatives on purpose: a real bacterial genome would
carry rRNA operons and quietly break the "without_rrna" guarantee.
"""
import random, sys, os

random.seed(2026)
READLEN, PAIRS = 150, 5000
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def rc(s):
    return ''.join(COMP.get(c, 'N') for c in reversed(s))


def load_fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def rand_dna(n):
    return ''.join(random.choice('ACGT') for _ in range(n))


def pairs_from(seq, tag, n):
    """Simulate fragments: R1 from the start, R2 as the reverse complement of the end."""
    out = []
    frag_min = 2 * READLEN + 20
    for i in range(n):
        flen = random.randint(frag_min, min(len(seq), frag_min + 300))
        p = random.randint(0, max(0, len(seq) - flen))
        frag = seq[p:p + flen]
        out.append((f'{tag}_{i}', frag[:READLEN], rc(frag[-READLEN:])))
    return out


def pairs_random(tag, n):
    out = []
    for i in range(n):
        frag = rand_dna(2 * READLEN + random.randint(20, 300))
        out.append((f'{tag}_{i}', frag[:READLEN], rc(frag[-READLEN:])))
    return out


def write_pair(name, recs):
    for mate in (1, 2):
        with open(f'{name}_R{mate}.fastq', 'w') as f:
            for rid, r1, r2 in recs:
                s = r1 if mate == 1 else r2
                f.write(f'@{rid}/{mate}\n{s}\n+\n{"I" * len(s)}\n')
    print(f'  {name}_R1.fastq / {name}_R2.fastq — {len(recs)} pairs')


def read_fastq_pairs(p1, p2, n):
    """Take real paired reads straight from the GIAB files."""
    out = []
    with open(p1) as f1, open(p2) as f2:
        while len(out) < n:
            h1, s1, _, _ = (f1.readline(), f1.readline(), f1.readline(), f1.readline())
            h2, s2, _, _ = (f2.readline(), f2.readline(), f2.readline(), f2.readline())
            if not h1 or not h2:
                break
            s1, s2 = s1.strip(), s2.strip()
            if 'N' in s1 or 'N' in s2:      # keep the truth clean
                continue
            out.append((f'HUMAN_{len(out)}', s1, s2))
    return out


if __name__ == '__main__':
    scratch = sys.argv[1] if len(sys.argv) > 1 else '.'
    r16 = load_fasta('/tmp/J01859.1.fa')       # E. coli 16S  — SILVA SSU
    r18 = load_fasta('/tmp/NR_003286.4.fa')    # human 18S    — SILVA SSU

    print('building test sets:')

    human = read_fastq_pairs(f'{scratch}/hg002_R1.fastq', f'{scratch}/hg002_R2.fastq', PAIRS)
    write_pair('sample_with_human', human)

    no_human = pairs_random('RANDOM', PAIRS // 2) + pairs_from(r16, 'ECOLI16S', PAIRS // 2)
    random.shuffle(no_human)
    write_pair('sample_without_human', no_human)

    rrna = pairs_from(r16, 'ECOLI16S', PAIRS // 2) + pairs_from(r18, 'HUMAN18S', PAIRS // 2)
    random.shuffle(rrna)
    write_pair('sample_with_rrna', rrna)

    write_pair('sample_without_rrna', pairs_random('RANDOM', PAIRS))

    print("""
expected results (threshold 0.5):

  sample                 human index      rRNA index
  ---------------------------------------------------------
  with_human             ~all removed     ~none removed
  without_human          ~none removed    half removed (the 16S half)
  with_rrna              18S half removed ~all removed
  without_rrna           ~none removed    ~none removed""")
