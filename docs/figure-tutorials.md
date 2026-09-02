# Figure tutorials and interpretation guide

Each output is a reusable analysis component. The commands below write a new directory, so runs do not overwrite one another.

## Before you plot

Roary's CSV includes gene-cluster metadata followed by one column per genome; the Rtab is a binary gene-by-genome matrix. Read the [Roary output documentation](https://sanger-pathogens.github.io/Roary/) before changing column names or filtering rows. Keep the original output immutable and record the Roary command, version, GFF annotation source, `-i`, `-cd`, and QC rules.

For the supplied table, the word “core” in the filename is historical. It is not a strict core matrix: only 13 of 196 families are present in all 344 represented genomes.

## Baseline panel figures

Run all baseline figures with:

```bash
python scripts/roary_pangenome_figures.py \
  --input data/core_genes_reps.csv \
  --outdir output/current_table \
  --seed 42
```

### 1. Gene-family prevalence distribution

**Input:** a supported gene-by-genome presence/absence table.

**Read it as:** the percentage of supplied genomes in which each family occurs. A broad distribution means the selected panel includes variable families.

**Do not claim:** total pangenome partition sizes unless the input is the complete, unfiltered Roary matrix. Prevalence can also be affected by assembly fragmentation, annotation, and clustering.

### 2. Presence/absence heatmap

**Input:** a supported table; `--max-genes` and `--max-genomes` set readable display limits.

```bash
python scripts/roary_pangenome_figures.py \
  --input roary_output/gene_presence_absence.csv \
  --outdir output/heatmap_120x80 \
  --max-genes 120 --max-genomes 80 --seed 42
```

**Read it as:** rows are selected because their prevalence is most variable, then sorted by presence within a deterministic genome subset. Blocks are hypotheses to check against a phylogeny and metadata, not inferred clades.

**Reuse:** pair a selected matrix with a vetted Newick tree in [ggtree + ComplexHeatmap Phylogenomics](https://github.com/mbilal-OU/ggtree-complexheatmap-phylogenomics) when evolutionary alignment is needed.

### 3. Per-genome content

**Input:** a supported table.

**Read it as:** each point is the number of supplied families in one genome. Low values may reflect genuine gene loss, lower assembly/annotation completeness, or panel selection.

**Next check:** compare outlying genomes to assembly statistics and annotation counts before concluding biological loss.

### 4. Prevalence and copy number

**Input:** a Roary CSV containing `No. isolates` and `No. sequences`.

**Read it as:** a y-axis value above one means Roary reports more than one sequence per represented isolate, which can flag multi-copy families or clusters that were not split. The supplied panel has a mean of one for all rows, so it is a negative diagnostic.

**Do not claim:** confirmed paralogy without inspecting the genes, annotations, and clustering settings.

### 5. Jaccard gene-content PCoA

**Input:** a supported table. Add metadata only after resolving genome identifiers exactly.

```bash
python scripts/roary_pangenome_figures.py \
  --input roary_output/gene_presence_absence.csv \
  --outdir output/ordination \
  --metadata metadata.tsv --group-column lineage --seed 42
```

`metadata.tsv` needs an identifier column named `genome`, `sample`, `isolate`, or `accession`, and a `lineage` column in this example.

**Read it as:** nearby points have more similar gene-content profiles under Jaccard distance. Axes are descriptive. Colouring does not establish significant group separation; use a pre-specified test and dispersion diagnostics.

## Full-pangenome figures

These are enabled only for an unfiltered `gene_presence_absence.csv` or `gene_presence_absence.Rtab`.

```bash
python scripts/roary_pangenome_figures.py \
  --input roary_output/gene_presence_absence.csv \
  --outdir output/full_roary \
  --full-pangenome --iterations 1000 --seed 42
```

### 6. Permutation accumulation summary

The script samples `--iterations` random genome orders, calculates cumulative observed families and newly observed families at each added genome, then draws median and 95% empirical intervals. The seed and count are saved in `run_manifest.json`.

**Use:** report whether additional sampling continues to discover families within the cohort, together with uncertainty over order.

**Limit:** this is not an open/closed verdict. For model-based openness work and uncertainty, use [PanGenome Openness Estimator](https://github.com/mbilal-OU/PanGenome-Openness-Estimator).

### 7. Threshold sensitivity

This figure counts families meeting 90%, 95%, 99%, and 100% prevalence. It shows how a “core-like” membership count changes with its definition.

**Limit:** changing the threshold changes a reporting definition. It does not rerun or modify Roary's gene clustering.

## Figures that need additional biological inputs

| Figure | Required additional input | Appropriate repository or tool |
|---|---|---|
| Tree-aligned gene-content heatmap | Vetted core-genome tree plus metadata | [ggtree + ComplexHeatmap Phylogenomics](https://github.com/mbilal-OU/ggtree-complexheatmap-phylogenomics) |
| Gene-content association plot | Pre-specified phenotype, covariates, multiple-testing plan | [Scoary](https://github.com/AdmiralenOla/Scoary) or a controlled GWAS workflow |
| Gene neighbourhood or synteny figure | Gene coordinates or GenBank/GFF annotations | [Matplotlib Genomic Figures](https://github.com/mbilal-OU/matplotlib-genomic-figures) |
| Interactive phylogeny, metadata, and PAV browser | Tree, metadata, PAV matrix | [Phandango](https://jameshadfield.github.io/phandango/) |
| Partitioned pangenome graph | Graph-aware pangenome output | [PPanGGOLiN](https://ppanggolin.readthedocs.io/) |

The report is stronger when every figure has the inputs needed to support it than when a template produces unsupported panels.
