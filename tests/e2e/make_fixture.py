#!/usr/bin/env python3
"""Generate the end-to-end fixture: simulated PIP-seq V4 nanopore reads with known isoform truth.

Run once, from the repository root; the outputs are committed and CI never regenerates them:

    python3 tests/e2e/make_fixture.py

Inputs (committed): tests/e2e/fixture/reference.fa (GENCODE v32 isoforms of 12 genes),
tests/e2e/fixture/genes.tsv, tests/whitelist/synthetic_barcodes.csv.
Outputs: tests/e2e/fixture/{reads.fq.gz, truth.tsv, cells.tsv}.

Read layout (forward architecture): noise(4) bc1 ATG bc2 GAG bc3 TCGAG bc4 UMI(12) revcomp(cDNA).
Half the reads are reverse-complemented whole (the polyA-strand architecture), so both regex arms
run; about 5% carry a broken bc2-bc3 linker (GAG -> GAT), so both regexes miss and the seal must
rescue them. Errors are i.i.d. at the measured ONT rates (substitution 0.59%, insertion 0.33%,
deletion 0.78%) over the whole read, barcode included, so some barcodes need edit-1 correction
and some become ambiguous. Standard library only; fixed seeds.
"""
import gzip
import random
from collections import defaultdict
from pathlib import Path

FIX = Path("tests/e2e/fixture")
WL = Path("tests/whitelist/synthetic_barcodes.csv")
N_CELLS = 4
DEPTH = 16          # reads per gene per cell
DOM = 0.8           # dominant-isoform fraction; the dominant isoform alternates by cell
P_SUB, P_INS, P_DEL = 0.0059, 0.0033, 0.0078
P_REVERSE = 0.5
P_BROKEN_LINKER = 0.05
BASES = "ACGT"

rng = random.Random(20260930)
err = random.Random(12345)   # errors draw from their own stream


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rnd(n):
    return "".join(rng.choice(BASES) for _ in range(n))


def add_error(seq):
    out = []
    for b in seq:
        r = err.random()
        if r < P_DEL:
            continue
        if r < P_DEL + P_SUB:
            b = err.choice([x for x in BASES if x != b])
        out.append(b)
        if err.random() < P_INS:
            out.append(err.choice(BASES))
    return "".join(out)


seq, cur = {}, None
for ln in open(FIX / "reference.fa"):
    ln = ln.strip()
    if ln.startswith(">"):
        cur = ln[1:].split()[0]
        seq[cur] = []
    else:
        seq[cur].append(ln)
seq = {t: "".join(s) for t, s in seq.items()}

gene_tx = defaultdict(list)
for ln in list(open(FIX / "genes.tsv"))[1:]:
    t, g, i, _l, _u = ln.rstrip("\n").split("\t")
    gene_tx[g].append((int(i), t))
gene_tx = {g: [t for _, t in sorted(v)] for g, v in gene_tx.items()}

rows = [ln.strip().split(",") for ln in open(WL)][1:]   # columns: bc4, bc3, bc2, bc1
cells, used = [], set()
while len(cells) < N_CELLS:
    b = (rng.choice(rows)[3], rng.choice(rows)[2], rng.choice(rows)[1], rng.choice(rows)[0])
    if b not in used:
        used.add(b)
        cells.append((f"cell{len(cells)}", *b, len(cells) % 2))

reads, truth = [], defaultdict(int)
for cid, b1, b2, b3, b4, state in cells:
    for g, txs in gene_tx.items():
        dom = state if state < len(txs) else 0
        for k in range(DEPTH):
            if rng.random() < DOM or len(txs) == 1:
                t = txs[dom]
            else:
                t = rng.choice([x for x in txs if x != txs[dom]])
            linker2 = "GAT" if rng.random() < P_BROKEN_LINKER else "GAG"
            read = rnd(4) + b1 + "ATG" + b2 + linker2 + b3 + "TCGAG" + b4 + rnd(12) + rc(seq[t])
            if rng.random() < P_REVERSE:
                read = rc(read)
            reads.append((f"{cid}.{g}.{t}.{k}", add_error(read)))
            truth[(cid, t)] += 1

rng.shuffle(reads)
with gzip.GzipFile(FIX / "reads.fq.gz", "wb", mtime=0) as fo:
    for name, s in reads:
        fo.write(f"@{name}\n{s}\n+\n{'I' * len(s)}\n".encode())
with open(FIX / "cells.tsv", "w") as fo:
    fo.write("cell\tbc1\tbc2\tbc3\tbc4\tdominant_isoform_idx\n")
    for c in cells:
        fo.write("\t".join(map(str, c)) + "\n")
with open(FIX / "truth.tsv", "w") as fo:
    fo.write("cell\ttranscript_id\ttrue_reads\n")
    for (cid, t), n in sorted(truth.items()):
        fo.write(f"{cid}\t{t}\t{n}\n")
print(f"{len(cells)} cells, {len(gene_tx)} genes, {len(reads)} reads, {len(truth)} truth rows")
