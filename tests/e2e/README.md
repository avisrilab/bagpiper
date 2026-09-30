# End-to-end test

`tests/e2e.rs` runs the built binary on a small simulated fixture and compares the result with the
recorded output of the same release in `expected/<version>/` (the version comes from
`bagpiper --version`).

- Leg A: `barcode` on `fixture/reads.fq.gz`. Records the passed reads (read, barcode, UMI, cDNA
  length and FNV-1a digest), the failed read names, and the `Total:` line; the assignment-path
  counters, when printed, must sum to `Matched`.
- Leg B: `count --b1` on `fixture/aligned.bam`, a frozen alignment of the 0.1.1 barcode output
  (minimap2 2.31-r1302, `-ax map-ont --for-only -N 200 -p 0.9` against `fixture/reference.fa`,
  then `samtools sort -n`). 0.1.0 and 0.1.1 produce identical matrices on this BAM; the eqclass
  dedup fix of 0.1.1 is not exercised here, so leg A carries the difference between them.
- Leg C: `count --r1 --reference`, only for releases with the internal aligner.

Legs B and C also check the matrix against the simulated truth: all four cells are called, and on
genes with at least 10% unique sequence (`genes.tsv`, column `gene_unique_frac`) the dominant
isoform is right in at least 95% of cell-genes.

## Fixture

Simulated only: no real reads and no real whitelist are in this repository. `make_fixture.py`
(standard library, fixed seeds) writes `reads.fq.gz`, `truth.tsv` and `cells.tsv` from
`reference.fa` (GENCODE v32 isoforms of 12 genes, three from each unique-sequence bin: below 1%,
1-5%, 5-10%, and 10% or more) and `tests/whitelist/synthetic_barcodes.csv`. 4 cells x 12 genes x
16 reads; half the reads reverse-complemented, about 5% with a broken bc2-bc3 linker (seal path),
i.i.d. errors at the measured ONT rates (substitution 0.59%, insertion 0.33%, deletion 0.78%).
`gene_unique_frac` is the unique-31-mer fraction of the gene's isoforms, computed in the
BenchDrop-seq revision analysis; it cannot be recomputed from the fixture. `aligned.cb_cell.tsv` maps each barcode in the frozen BAM to its simulated cell by majority.

## Running

    cargo test --test e2e                                  # this build
    BAGPIPER_BIN=/path/to/bagpiper cargo test --test e2e   # another build, e.g. an older release

A new release records its expected output once, from its own build:

    BAGPIPER_E2E_BLESS=1 cargo test --release --test e2e
