# Report audit: Roary pangenome heatmaps

## Executive assessment

The original repository contained one static heatmap, one interactive HTML file, and an exploratory notebook stored with an `.R` extension even though it is a Jupyter notebook. Its README referenced missing paths and described the input as a core-gene file without checking whether all rows were core. The static heatmap also used a random genome selection and gave no report-level route for using another dataset.

The rebuilt repository has a tested command-line pipeline, five checked-in figures from the supplied panel, SVG exports, TSV summaries, an HTML report, explicit input contracts, and two full-matrix figures that are locked behind `--full-pangenome`. The exploratory notebook is now correctly labelled `.ipynb`; it is archival, while the Python pipeline is the supported route.

## Current data evidence

The checked-in `core_genes_reps.csv` has 196 supplied families across 344 genome accessions. The baseline run finds 13 families in all genomes and a median family prevalence of 63.08%. This confirms that the file is a variable representative panel, not a strict all-genome core matrix. The generated figures make descriptive statements about this panel only.

## Why these figures

Roary identifies gene clusters and writes both a detailed CSV and a binary Rtab matrix. The official output guide defines the sequence-count and isolate-count fields used in the copy-number diagnostic and makes the matrix suitable for downstream analysis. [Roary documentation](https://sanger-pathogens.github.io/Roary/)

The recommended core report starts with prevalence, a readable PAV heatmap, genome content, copy-number diagnostics, and gene-content ordination because these answer different questions from the same matrix. A full matrix can additionally support order-aware accumulation and threshold sensitivity. The protocol by Sitto and Battistuzzi identifies `summary_statistics.txt` and `gene_presence_absence.csv` as the central Roary outputs for interpreting pangenome composition. [Estimating Pangenomes with Roary](https://academic.oup.com/mbe/article/37/3/933/5652084)

Tree-linked PAV displays and interactive investigation are valuable once a defensible phylogeny is available. Phandango was designed to combine bacterial population-genomics outputs in a browser rather than treating a heatmap as an isolated result. [Hadfield et al., 2018](https://academic.oup.com/bioinformatics/article/34/2/292/4212949)

Different pangenome callers and parameter choices can change the gene-presence matrix itself. Panaroo's benchmarking work compares accumulation curves from Roary and Panaroo matrices, so a reported pangenome conclusion should preserve both the input matrix and pipeline settings. [Tonkin-Hill et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9977150/)

For graph-aware partitioning, PPanGGOLiN models persistent, shell, and cloud components and provides statistics and visual summaries. That is a different analytical object from a filtered Roary panel, so the atlas links to it rather than creating a misleading graph figure. [PPanGGOLiN documentation](https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeStat.html)

## Recommended report order

1. Cohort and quality-control table: genome counts, taxonomic scope, assembly and annotation criteria.
2. Pangenome construction table: tool version, identity threshold, core threshold, and command.
3. Prevalence distribution and PAV heatmap.
4. Genome gene-content distribution and PCoA with pre-specified metadata.
5. Tree-aligned PAV heatmap only when the tree is independently justified.
6. Full-matrix accumulation and threshold sensitivity, with seed and iteration count.
7. Trait-association or gene-neighbourhood figures only when their required phenotype or coordinate inputs are available.

## Claims this repository will not make

- “This is a strict core genome” from the current input.
- “The pangenome is open or closed” from one curve, one threshold, or a filtered table.
- “These clusters are phylogenetic lineages” from a heatmap or PCoA alone.
- “A family is a paralogue” from a copy-number diagnostic alone.

Those boundaries are part of the report design, not a limitation of the visual style.
