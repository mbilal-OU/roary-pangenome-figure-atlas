# Roary Pangenome Figure Atlas

Reproducible, evidence-aware figures for bacterial gene presence/absence data produced by [Roary](https://sanger-pathogens.github.io/Roary/). The repository turns a Roary table into a report with publication-ready PNG and SVG figures, reusable TSV summaries, and a lightweight HTML index.

This is not a gallery of attractive heatmaps. Each figure has a defined input, command, interpretation, and limitation. The current example is deliberately described as a **representative gene-family panel**, not a complete pangenome analysis.

## What the checked-in data supports

`data/core_genes_reps.csv` contains 196 gene families across 344 GCF genomes. It includes 13 families present in every supplied genome; the remaining families occur in 50.0% to 100.0% of genomes. It is therefore useful for descriptive presence/absence, prevalence, copy-number, and gene-content structure figures. It does **not** support a claim about total pangenome size, open/closed status, or core/shell/cloud partition sizes.

The checked-in output was generated with seed `42`. The exact summary is in [run_manifest.json](figures/current_table/run_manifest.json).

| Figure | What it answers | Figure | Tutorial and interpretation |
|---|---|---|---|
| Prevalence distribution | Which supplied families are common or variable? | [PNG](figures/current_table/01_prevalence_distribution.png) | [Read](docs/figure-tutorials.md#1-gene-family-prevalence-distribution) |
| Presence/absence heatmap | Which variable families co-occur across a reproducible display subset? | [PNG](figures/current_table/02_presence_absence_heatmap.png) | [Read](docs/figure-tutorials.md#2-presenceabsence-heatmap) |
| Genome content | How many supplied families occur in each genome? | [PNG](figures/current_table/03_genome_gene_content.png) | [Read](docs/figure-tutorials.md#3-per-genome-content) |
| Copy-number diagnostic | Are prevalent families also represented by more than one sequence per isolate? | [PNG](figures/current_table/04_prevalence_copy_number.png) | [Read](docs/figure-tutorials.md#4-prevalence-and-copy-number) |
| Jaccard PCoA | Do genomes separate in the supplied gene-content space? | [PNG](figures/current_table/05_gene_content_pcoa.png) | [Read](docs/figure-tutorials.md#5-jaccard-gene-content-pcoa) |

The [HTML report](figures/current_table/report.html), [gene-family prevalence table](figures/current_table/gene_family_prevalence.tsv), and [genome-content table](figures/current_table/genome_content.tsv) are reusable project outputs.

## Quick start

```bash
git clone https://github.com/mbilal-OU/roary-pangenome-figure-atlas.git
cd roary-pangenome-figure-atlas
python -m pip install -r requirements.txt

python scripts/roary_pangenome_figures.py \
  --input data/core_genes_reps.csv \
  --outdir output/current_table \
  --seed 42
```

Open `output/current_table/report.html`. Every PNG also has an SVG version for manuscript editing.

## Use your own Roary data

### A representative or filtered table

Use this when your input is like the included file: Roary metadata columns followed by accession-named genome columns. This route produces Figures 1 to 5 only.

```bash
python scripts/roary_pangenome_figures.py \
  --input path/to/core_genes_reps.csv \
  --outdir output/representative_panel \
  --seed 42 --max-genes 75 --max-genomes 60
```

### A complete Roary pangenome matrix

Use the unfiltered `gene_presence_absence.csv` from the Roary output directory. This additionally generates a permutation accumulation summary and a threshold-sensitivity plot. The same route accepts `gene_presence_absence.Rtab`.

```bash
python scripts/roary_pangenome_figures.py \
  --input roary_output/gene_presence_absence.csv \
  --outdir output/full_roary \
  --full-pangenome --iterations 1000 --seed 42
```

The optional full-input figures are documented in [the tutorial](docs/figure-tutorials.md#full-pangenome-figures). The workflow refuses `--full-pangenome` for a filtered table, so it cannot silently overstate an incomplete input.

### Colour PCoA by your metadata

Supply a CSV or TSV with a genome identifier column named `genome`, `sample`, `isolate`, or `accession`, plus the grouping column.

```bash
python scripts/roary_pangenome_figures.py \
  --input roary_output/gene_presence_absence.csv \
  --outdir output/host_coloured \
  --metadata metadata.tsv --group-column host \
  --full-pangenome --iterations 1000 --seed 42
```

## Figure contract

| Output | Input required | Reproducibility controls | Correct interpretation |
|---|---|---|---|
| `01_prevalence_distribution` | Any supported table | Deterministic | Describes supplied families, not the whole pangenome unless the input is complete. |
| `02_presence_absence_heatmap` | Any supported table | Deterministic rank-based subset | A readable panel of variable families. The display subset is not a statistical sample or a phylogeny. |
| `03_genome_gene_content` | Any supported table | Deterministic | Counts supplied families per genome. Differences can reflect biology, annotation, assembly, or input filtering. |
| `04_prevalence_copy_number` | Roary CSV containing sequence/isolate fields | Deterministic | Flags possible multi-copy or clustering behaviour. It does not prove paralogy. |
| `05_gene_content_pcoa` | Any supported table | Deterministic | Descriptive Jaccard structure. Colouring or separation is not evidence of association without a pre-specified test. |
| `06_permutation_accumulation` | Complete Roary matrix plus `--full-pangenome` | `--iterations`, `--seed` | Shows order-sensitive accumulation uncertainty. It does not label a pangenome open or closed. |
| `07_threshold_sensitivity` | Complete Roary matrix plus `--full-pangenome` | Deterministic | Shows how a prevalence definition changes membership counts. It does not alter Roary clustering. |

## Biological reporting guardrails

- Report the Roary version, annotation pipeline, genome quality criteria, identity threshold, `-cd` setting, genome count, and all random seeds.
- A core-like prevalence threshold is a definition. It is not a universal biological boundary.
- Use accumulation curves as descriptive, order-aware summaries. For openness modelling and uncertainty, use the dedicated [PanGenome Openness Estimator](https://github.com/mbilal-OU/PanGenome-Openness-Estimator).
- PCoA axes describe variation in the supplied matrix. They do not establish population structure, host association, or phylogenetic relatedness without supporting analysis.
- If the cohort spans multiple species or a high taxonomic rank, avoid language such as “the species core genome.” State the sampling scope instead.

## Related visualization and analysis repositories

| Repository | Use it with this atlas when you need |
|---|---|
| [ggtree + ComplexHeatmap Phylogenomics](https://github.com/mbilal-OU/ggtree-complexheatmap-phylogenomics) | A tree-aligned presence/absence heatmap with clade annotations. Requires a defensible Newick tree and metadata. |
| [ggplot2 Omics Grammar](https://github.com/mbilal-OU/ggplot2-omics-grammar) | Manuscript-ready layered statistical figures. |
| [Plotly Interactive Omics](https://github.com/mbilal-OU/plotly-interactive-omics) | Interactive exploration of metadata, QC, and ordination. |
| [Matplotlib Genomic Figures](https://github.com/mbilal-OU/matplotlib-genomic-figures) | Custom genomic tracks, synteny, and sequence-aware layouts. |
| [Seaborn Biological Statistics](https://github.com/mbilal-OU/seaborn-biological-statistics) | Distributional and replicate-aware biological statistics figures. |
| [Shiny Omics Explorer](https://github.com/mbilal-OU/shiny-omics-explorer) | A shareable R/Shiny interface for exploratory omics data. |
| [Gnuplot Bioinformatics CLI](https://github.com/mbilal-OU/gnuplot-bioinformatics-cli) | Fast, headless figures on HPC systems. |
| [PanPhyloFlow](https://github.com/mbilal-OU/PanPhyloFlow) | Reproducible pangenome-to-core-phylogeny processing. |
| [PanOrd](https://github.com/mbilal-OU/PanOrd) | Dedicated ordination workbench for pangenome matrices. |

## Repository layout

```text
data/                       checked-in representative gene-family table
figures/current_table/      checked-in reproducible figures and tables
scripts/                    command-line figure generator
docs/                       report and figure tutorials
tests/                      end-to-end figure-pipeline check
notebooks/                  archived exploratory notebook (.ipynb)
```

The archived notebook is kept for provenance. The supported and tested route is `scripts/roary_pangenome_figures.py`.

## Research basis and verification

Roary documents the semantics of `gene_presence_absence.csv` and the binary `gene_presence_absence.Rtab` matrix in its [output guide](https://sanger-pathogens.github.io/Roary/). The figure decisions and limitations are expanded in [the report](docs/report.md), with links to Roary, pangenome-estimation, Phandango, Panaroo, and PPanGGOLiN literature.

GitHub Actions runs the supplied data through the pipeline and checks that the report, manifest, tables, and five baseline figures are created. Run the same check locally with:

```bash
python -m unittest discover -s tests -v
```
