#!/usr/bin/env python3
"""Jeu de test à vérité terrain connue.

Référence : génome aléatoire. Reads : moitié extraits de ce génome (donc à
écarter), moitié tirés indépendamment (donc à conserver). Les seconds portent
un identifiant préfixé KEEP_, ce qui permet de compter exactement les erreurs.
"""
import random, sys

random.seed(11)
GENOME = 200_000
N_CONTA = 5_000
N_PROPRE = 5_000
READLEN = 150

genome = ''.join(random.choice('ACGT') for _ in range(GENOME))

with open('ref.fasta', 'w') as f:
    f.write('>ref_synthetique\n')
    for i in range(0, len(genome), 70):
        f.write(genome[i:i + 70] + '\n')

def qual(n):
    return 'I' * n

reads = []
# contaminants : sous-chaînes exactes du génome de référence
for i in range(N_CONTA):
    p = random.randint(0, GENOME - READLEN)
    seq = genome[p:p + READLEN]
    if random.random() < 0.5:  # une moitié sur le brin inverse
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        seq = ''.join(comp[c] for c in reversed(seq))
    reads.append((f'CONTA_{i}', seq))

# propres : tirés indépendamment, aucun lien avec la référence
for i in range(N_PROPRE):
    seq = ''.join(random.choice('ACGT') for _ in range(READLEN))
    reads.append((f'KEEP_{i}', seq))

random.shuffle(reads)
with open('test.fastq', 'w') as f:
    for name, seq in reads:
        f.write(f'@{name}\n{seq}\n+\n{qual(len(seq))}\n')

print(f'ref.fasta   : {GENOME:,} pb')
print(f'test.fastq  : {len(reads):,} reads ({N_CONTA:,} à écarter, {N_PROPRE:,} à garder)')
