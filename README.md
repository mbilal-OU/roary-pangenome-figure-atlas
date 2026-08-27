# Prokaryotic Core Pangenome Heatmap Tutorial

This repository demonstrates how to build an **interactive heatmap** and a **static figure**
from a **core-gene file** (`core_genes_reps.csv`) derived from a *prokaryotic pangenome* analysis
(e.g. Roary / Panaroo).

The focus is on:

- **Prokaryotic pangenomes** (bacteria / archaea)
- A **core genes file** summarizing genes present in all genomes
- Using R to create:
  - an **interactive heatmap** with hover (via [`heatmaply`](https://cran.r-project.org/package=heatmaply))
  - a **static PNG heatmap** (via [`pheatmap`](https://cran.r-project.org/package=pheatmap))

---

## Repository Structure

```text
prokaryotic-core-pangenome-heatmap/
├── data/
│   └── core_genes_reps_example.csv       # optional small example (same structure as real core_genes_reps.csv)
├── notebooks/
│   └── core_heatmap_pipeline.R           # full R script for Google Colab / RStudio
├── figs/
│   └── core_genes_heatmap.png            # static PNG heatmap (25 genomes)
├── html/
│   └── core_genes_heatmap.html           # interactive heatmap (heatmaply)
└── README.md
```

> **Note:** The runnable script in this repository is `notebooks/core_heatmap_pipeline.R`. Use this filename when following or adapting the tutorial.
