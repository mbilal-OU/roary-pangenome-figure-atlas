#!/usr/bin/env python3
"""Create evidence-aware figures from Roary presence/absence output.

The script accepts either a Roary gene_presence_absence.csv table or a
representative table with Roary metadata columns and GCF_/GCA_ genome columns.
It deliberately reserves full-pangenome statements for an explicit complete
gene_presence_absence.csv input.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform


METADATA_COLUMNS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
    "absent",
    "present",
    "coverage",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate reusable, evidence-aware Roary pangenome figures."
    )
    parser.add_argument("--input", required=True, type=Path, help="Roary CSV or .Rtab input")
    parser.add_argument("--outdir", type=Path, default=Path("output"))
    parser.add_argument("--metadata", type=Path, help="Optional TSV/CSV with genome and group columns")
    parser.add_argument("--group-column", help="Metadata column used for the optional group figure")
    parser.add_argument("--iterations", type=int, default=1000, help="Permutation count for accumulation")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible sampling")
    parser.add_argument("--max-genes", type=int, default=75, help="Rows shown in the heatmap")
    parser.add_argument("--max-genomes", type=int, default=60, help="Columns shown in the heatmap")
    parser.add_argument("--full-pangenome", action="store_true", help="Enable figures requiring a complete gene_presence_absence.csv")
    return parser.parse_args()


def read_input(path: Path) -> tuple[pd.DataFrame, list[str], pd.DataFrame, bool]:
    if not path.exists():
        raise FileNotFoundError(f"Input does not exist: {path}")
    if path.suffix.lower() == ".rtab":
        raw = pd.read_csv(path, sep="\t")
        if raw.shape[1] < 3:
            raise ValueError("An Rtab file needs a gene column and at least two genome columns.")
        genes = raw.iloc[:, 0].fillna("unnamed_gene").astype(str)
        genomes = list(raw.columns[1:])
        matrix = raw.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int).clip(0, 1)
        matrix.index = make_unique(genes)
        details = pd.DataFrame({"Gene": matrix.index, "Annotation": ""})
        return details, genomes, matrix, True

    raw = pd.read_csv(path, low_memory=False)
    genome_columns = [c for c in raw.columns if str(c).startswith(("GCF_", "GCA_"))]
    if not genome_columns:
        genome_columns = [c for c in raw.columns if c not in METADATA_COLUMNS and not str(c).startswith("Unnamed")]
    if len(genome_columns) < 2:
        raise ValueError("Could not find at least two genome columns. Use Roary CSV or gene_presence_absence.Rtab.")
    gene_column = "Gene" if "Gene" in raw.columns else raw.columns[0]
    matrix = raw[genome_columns].fillna("").astype(str).apply(
        lambda col: (~col.str.strip().isin(["", "0", "NA", "nan"])).astype(int)
    )
    matrix.index = make_unique(raw[gene_column].fillna("unnamed_gene").astype(str))
    details = raw.copy()
    details["Gene"] = matrix.index
    return details, genome_columns, matrix, False


def make_unique(values: pd.Series) -> pd.Index:
    seen: dict[str, int] = {}
    out = []
    for value in values:
        count = seen.get(value, 0)
        seen[value] = count + 1
        out.append(value if count == 0 else f"{value}__{count + 1}")
    return pd.Index(out, name="Gene")


def save(fig: plt.Figure, outdir: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def select_heatmap(matrix: pd.DataFrame, max_genes: int, max_genomes: int) -> pd.DataFrame:
    prevalence = matrix.mean(axis=1)
    variable_rows = (prevalence * (1 - prevalence)).sort_values(ascending=False).index[:max_genes]
    genome_burden = matrix.loc[variable_rows].sum(axis=0).sort_values()
    if len(genome_burden) > max_genomes:
        positions = np.linspace(0, len(genome_burden) - 1, max_genomes).round().astype(int)
        genome_names = genome_burden.index[positions]
    else:
        genome_names = genome_burden.index
    return matrix.loc[variable_rows, genome_names].loc[
        matrix.loc[variable_rows, genome_names].sum(axis=1).sort_values().index
    ]


def prevalence_figure(prevalence: pd.Series, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    sns.histplot(prevalence * 100, binwidth=5, color="#1f6f8b", edgecolor="white", ax=ax)
    ax.axvline(99, color="#b63d32", linestyle="--", linewidth=1.4, label="99% reference")
    ax.set(xlabel="Gene-family prevalence across supplied genomes (%)", ylabel="Number of gene families")
    ax.set_title("Gene-family prevalence distribution")
    ax.legend(frameon=False)
    save(fig, outdir, "01_prevalence_distribution")


def heatmap_figure(matrix: pd.DataFrame, outdir: Path, max_genes: int, max_genomes: int) -> None:
    panel = select_heatmap(matrix, max_genes, max_genomes)
    fig, ax = plt.subplots(figsize=(13, 9))
    sns.heatmap(panel, cmap=sns.color_palette(["#f3f4f6", "#1f6f8b"], as_cmap=True),
                cbar_kws={"label": "Presence (1) / absence (0)"}, yticklabels=panel.index,
                xticklabels=panel.columns, ax=ax)
    ax.set(xlabel="Genome", ylabel="Gene family")
    ax.set_title("Most variable gene families (deterministic display subset)")
    ax.tick_params(axis="x", rotation=90, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    save(fig, outdir, "02_presence_absence_heatmap")


def genome_burden_figure(matrix: pd.DataFrame, outdir: Path) -> None:
    burden = matrix.sum(axis=0).sort_values()
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.plot(range(1, len(burden) + 1), burden.values, color="#1f6f8b", linewidth=1.3)
    ax.fill_between(range(1, len(burden) + 1), burden.values, color="#1f6f8b", alpha=0.16)
    ax.set(xlabel="Genome rank (fewest to most supplied gene families)", ylabel="Gene families present")
    ax.set_title("Per-genome content in the supplied gene-family table")
    save(fig, outdir, "03_genome_gene_content")


def prevalence_copy_figure(details: pd.DataFrame, prevalence: pd.Series, outdir: Path) -> None:
    if "No. sequences" not in details.columns:
        return
    copies = pd.to_numeric(details["No. sequences"], errors="coerce") / pd.to_numeric(
        details.get("No. isolates", np.nan), errors="coerce"
    )
    fig, ax = plt.subplots(figsize=(7.2, 5))
    ax.scatter(prevalence.values * 100, copies, color="#c77d2f", alpha=0.75, s=22, linewidths=0)
    ax.axhline(1, color="#5b6470", linestyle="--", linewidth=1)
    ax.set(xlabel="Gene-family prevalence (%)", ylabel="Mean copies per represented isolate")
    ax.set_title("Prevalence and copy-number diagnostic")
    save(fig, outdir, "04_prevalence_copy_number")


def pcoa_figure(matrix: pd.DataFrame, outdir: Path, metadata: Optional[pd.DataFrame], group_column: Optional[str]) -> None:
    distances = squareform(pdist(matrix.T.values, metric="jaccard"))
    n = distances.shape[0]
    center = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * center @ (distances ** 2) @ center
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = np.maximum(values[:2], 0)
    coords = vectors[:, :2] * np.sqrt(positive)
    total = values[values > 0].sum()
    explained = [100 * value / total if total > 0 else 0 for value in positive]
    fig, ax = plt.subplots(figsize=(7.4, 5.7))
    colors = "#5a8f7b"
    if metadata is not None and group_column and group_column in metadata.columns:
        genome_col = next((c for c in metadata.columns if c.lower() in {"genome", "sample", "isolate", "accession"}), None)
        if genome_col:
            groups = metadata.set_index(genome_col).reindex(matrix.columns)[group_column].fillna("Unassigned")
            for group in groups.unique():
                idx = groups.eq(group).values
                ax.scatter(coords[idx, 0], coords[idx, 1], s=28, alpha=0.8, label=str(group))
            ax.legend(title=group_column, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=26, alpha=0.8, color=colors)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=26, alpha=0.8, color=colors)
    ax.axhline(0, color="#d1d5db", linewidth=0.8)
    ax.axvline(0, color="#d1d5db", linewidth=0.8)
    ax.set(xlabel=f"PCoA 1 ({explained[0]:.1f}% positive variance)", ylabel=f"PCoA 2 ({explained[1]:.1f}% positive variance)")
    ax.set_title("Jaccard PCoA of supplied gene content")
    save(fig, outdir, "05_gene_content_pcoa")


def full_pangenome_figures(matrix: pd.DataFrame, outdir: Path, iterations: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    n_genomes = matrix.shape[1]
    curves = np.empty((iterations, n_genomes), dtype=float)
    new_genes = np.empty_like(curves)
    values = matrix.values.astype(bool)
    for i in range(iterations):
        order = rng.permutation(n_genomes)
        observed = np.zeros(values.shape[0], dtype=bool)
        for j, column in enumerate(order):
            newly_seen = values[:, column] & ~observed
            new_genes[i, j] = newly_seen.sum()
            observed |= values[:, column]
            curves[i, j] = observed.sum()
    x = np.arange(1, n_genomes + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, series, color in [("Pangenome size", curves, "#1f6f8b"), ("New families", new_genes, "#c77d2f")]:
        median = np.median(series, axis=0)
        low, high = np.quantile(series, [0.025, 0.975], axis=0)
        ax.plot(x, median, label=label, color=color, linewidth=2)
        ax.fill_between(x, low, high, color=color, alpha=0.16)
    ax.set(xlabel="Genomes sampled", ylabel="Gene families", title=f"Permutation accumulation summary ({iterations:,} orders; seed {seed})")
    ax.legend(frameon=False)
    save(fig, outdir, "06_permutation_accumulation")

    prevalence = matrix.mean(axis=1)
    thresholds = np.array([0.90, 0.95, 0.99, 1.00])
    counts = np.array([(prevalence >= threshold).sum() for threshold in thresholds])
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(thresholds * 100, counts, marker="o", color="#b63d32", linewidth=2)
    ax.set(xlabel="Prevalence threshold (%)", ylabel="Gene families meeting threshold")
    ax.set_title("Core-like membership sensitivity to threshold")
    for x_value, y_value in zip(thresholds * 100, counts):
        ax.annotate(f"{y_value:,}", (x_value, y_value), xytext=(0, 7), textcoords="offset points", ha="center")
    save(fig, outdir, "07_threshold_sensitivity")


def write_report(outdir: Path, details: pd.DataFrame, genomes: list[str], matrix: pd.DataFrame, full: bool, args: argparse.Namespace) -> None:
    prevalence = matrix.mean(axis=1)
    summary = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "input": str(args.input),
        "input_scope": "complete gene_presence_absence table" if full else "representative/partial gene-family table",
        "gene_families": int(matrix.shape[0]),
        "genomes": len(genomes),
        "strictly_present_in_all_genomes": int((prevalence == 1).sum()),
        "median_prevalence_percent": round(float(prevalence.median() * 100), 2),
        "seed": args.seed,
        "iterations": args.iterations if full else None,
    }
    (outdir / "run_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame({"gene_family": matrix.index, "prevalence": prevalence.values, "prevalence_percent": prevalence.values * 100}).to_csv(
        outdir / "gene_family_prevalence.tsv", sep="\t", index=False
    )
    pd.DataFrame({"genome": matrix.columns, "gene_families_present": matrix.sum(axis=0).values}).to_csv(
        outdir / "genome_content.tsv", sep="\t", index=False
    )
    scope_note = (
        "Full-pangenome input was explicitly enabled; accumulation and threshold figures are descriptive permutation summaries."
        if full
        else "This input is a supplied representative/partial table. It supports descriptive panel-level figures only; it cannot establish pangenome size, openness, or partition sizes."
    )
    figure_list = "\n".join(f"<li><a href='{path.name}'>{path.stem.replace('_', ' ')}</a></li>" for path in sorted(outdir.glob("*.png")))
    (outdir / "report.html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>Roary figure report</title>"
        "<style>body{font-family:system-ui;max-width:900px;margin:2rem auto;line-height:1.5}code{background:#f3f4f6;padding:.15rem .3rem}</style></head><body>"
        "<h1>Roary pangenome figure report</h1>"
        f"<p><strong>Scope:</strong> {scope_note}</p>"
        f"<p><strong>Input:</strong> <code>{args.input}</code><br><strong>Gene families:</strong> {matrix.shape[0]:,}<br><strong>Genomes:</strong> {len(genomes):,}<br><strong>Strict all-genome families:</strong> {(prevalence == 1).sum():,}</p>"
        f"<h2>Figures</h2><ul>{figure_list}</ul>"
        "<p>Interpret figures alongside the repository tutorial and report. Do not use one accumulation trajectory or a threshold count alone to label a pangenome open or closed.</p>"
        "</body></html>\n"
    )


def main() -> None:
    args = arguments()
    if args.iterations < 20:
        raise ValueError("--iterations must be at least 20 for a stable uncertainty ribbon.")
    details, genomes, matrix, is_rtab = read_input(args.input)
    requested_full = args.full_pangenome
    input_named_full = args.input.name == "gene_presence_absence.csv" or is_rtab
    full = requested_full and input_named_full
    if requested_full and not input_named_full:
        raise ValueError("--full-pangenome requires gene_presence_absence.csv or gene_presence_absence.Rtab, not a filtered representative table.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata = None
    if args.metadata:
        separator = "\t" if args.metadata.suffix.lower() in {".tsv", ".tab"} else ","
        metadata = pd.read_csv(args.metadata, sep=separator)
    sns.set_theme(style="whitegrid", context="notebook")
    prevalence = matrix.mean(axis=1)
    prevalence_figure(prevalence, args.outdir)
    heatmap_figure(matrix, args.outdir, args.max_genes, args.max_genomes)
    genome_burden_figure(matrix, args.outdir)
    prevalence_copy_figure(details, prevalence, args.outdir)
    pcoa_figure(matrix, args.outdir, metadata, args.group_column)
    if full:
        full_pangenome_figures(matrix, args.outdir, args.iterations, args.seed)
    write_report(args.outdir, details, genomes, matrix, full, args)
    print(f"Wrote figures, tables, and report to {args.outdir}")


if __name__ == "__main__":
    main()
