#!/usr/bin/env python3
"""Export matched-sample genus relative abundances for MaAsLin3 / ALDEx3."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "integrated" / "results" / "additional_analyses_2026-07-19"
TABLES = OUT / "tables"
GENUS_TSV = (
    ROOT
    / "microbiome_analysis"
    / "zr24558.16S_250813.zymo"
    / "00...AllSamples.Bac16Sv34"
    / "Heatmap"
    / "6.genus"
    / "taxa_abun_table_L6.txt"
)
ASV_CSV = (
    ROOT
    / "microbiome_analysis"
    / "zr24558.16S_250813.zymo"
    / "00...AllSamples.Bac16Sv34"
    / "DADA2_ASV_Distribution"
    / "ASV_Abundance_Table.csv"
)


def clean_label(x: str) -> str:
    s = str(x).strip()
    s = re.sub(r"\.+$", "", s)
    s = re.sub(r"\.\.+", ".", s)
    return s


def genus_name(taxon: str) -> str:
    # Prefer g__Name; fall back to last semicolon field.
    t = str(taxon)
    m = re.search(r"g__([A-Za-z0-9_]+)", t)
    if m:
        return m.group(1)
    parts = [p for p in t.split(";") if p and p not in ("__", "g__")]
    return parts[-1] if parts else "Unknown"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    both = pd.read_csv(TABLES / "matched_scfa_16s_qpcr_both.csv")
    both["customer_label"] = both["customer_label"].map(clean_label)
    both["sample_key"] = (
        both["donor_id"].astype(str)
        + "__"
        + both["aliquot_id"].astype(str)
        + "__"
        + both["carbohydrate"].astype(str)
        + "__"
        + both["well_repeat"].astype(str)
        + "__"
        + both["timepoint"].astype(str)
    )

    g = pd.read_csv(GENUS_TSV, sep="\t")
    tax = g["Taxon"].map(genus_name)
    mat = g.drop(columns=["Taxon"])
    mat.columns = [clean_label(c) for c in mat.columns]
    # Aggregate duplicate genus names
    mat.insert(0, "genus", tax)
    mat = mat.groupby("genus", as_index=True).sum(numeric_only=True)
    # relative within sample (columns are samples)
    rel = mat.div(mat.sum(axis=0).replace(0, pd.NA), axis=1).fillna(0.0)

    # Transpose to samples × genera
    rel_t = rel.T
    rel_t.index.name = "customer_label"
    rel_t = rel_t.reset_index()

    # ASV counts for ALDEx3 (aggregate to genus via taxonomy not available here;
    # use genus relative × depth as pseudo-counts; true counts need tax map.
    asv = pd.read_csv(ASV_CSV)
    sample_col = asv.columns[0]
    depth = asv.drop(columns=[sample_col]).sum(axis=1)
    depth_df = pd.DataFrame(
        {
            "customer_label": [clean_label(x) for x in asv[sample_col]],
            "asv_read_depth": depth.values.astype(float),
        }
    )

    merged = both[["sample_key", "customer_label", "donor_id", "aliquot_id", "carbohydrate", "timepoint", "group", "gene_copies_per_ul"]].merge(
        rel_t, on="customer_label", how="inner"
    )
    merged = merged.merge(depth_df, on="customer_label", how="left")

    # Drop non-feature cols for relative matrix export
    meta_cols = [
        "sample_key",
        "customer_label",
        "donor_id",
        "aliquot_id",
        "carbohydrate",
        "timepoint",
        "group",
        "gene_copies_per_ul",
        "asv_read_depth",
    ]
    feature_cols = [c for c in merged.columns if c not in meta_cols]
    # Integer-ish count approximation for ALDEx3 Dirichlet: round(rel * depth)
    counts = merged.copy()
    for c in feature_cols:
        counts[c] = (merged[c].astype(float) * merged["asv_read_depth"].fillna(0)).round().astype(int)

    out_rel = merged[["sample_key"] + feature_cols]
    out_counts = counts[["sample_key"] + feature_cols]
    out_rel.to_csv(TABLES / "genus_relative_for_da.csv", index=False)
    out_counts.to_csv(TABLES / "genus_counts_matched.csv", index=False)
    merged[meta_cols].to_csv(TABLES / "genus_da_sample_meta.csv", index=False)
    print(
        f"Wrote {out_rel.shape[0]} samples × {len(feature_cols)} genera to genus_relative_for_da.csv"
    )


if __name__ == "__main__":
    main()
