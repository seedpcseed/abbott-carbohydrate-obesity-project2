#!/usr/bin/env python3
"""Build exact well-matched SCFA–16S–qPCR analysis table (Phase A)."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "integrated" / "results" / "additional_analyses_2026-07-19"
AUDIT = OUT / "data_audit"
TABLES = OUT / "tables"

SCFA_CSV = (
    ROOT
    / "scfa_metabolomics"
    / "data"
    / "removed_qcs_quant_results_20250722_PFBBr_PS2720_20250731.csv"
)
QPCR_CSV = (
    ROOT
    / "microbiome_analysis"
    / "zr24558.16S_250813.zymo"
    / "00...AllSamples.Bac16Sv34"
    / "Sample_Information"
    / "absolute.abundance.csv"
)
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
META_XLSX = (
    ROOT
    / "scfa_metabolomics"
    / "metadata"
    / "SEED LAB HMM Submission Sheet 202405.xlsx"
)
CANON_SCFA = ROOT / "integrated" / "metadata" / "canonical_experimental_units_scfa.csv"


def clean_label(x: str) -> str:
    s = str(x).strip()
    s = re.sub(r"\.+$", "", s)
    s = s.replace("..", ".")
    return s


def parse_micro_label(label: str) -> dict | None:
    lab = clean_label(label)
    m = re.match(
        r"^(?P<aliquot>\d+[A-Za-z])\.(?P<cond>R[12]|S[12]|NC)\.*(?P<time>0H|48H)$",
        lab,
        flags=re.I,
    )
    if not m:
        return None
    aliquot = m.group("aliquot").upper()
    cond = m.group("cond").upper()
    timepoint = m.group("time").upper()
    if cond.startswith("R"):
        carb, well = "RDC", int(cond[1])
    elif cond.startswith("S"):
        carb, well = "SDC", int(cond[1])
    else:
        carb, well = "No Added Carb", 1
    donor = re.match(r"^(\d+)", aliquot).group(1)
    return {
        "customer_label": lab,
        "donor_id": donor,
        "aliquot_id": aliquot,
        "carbohydrate": carb,
        "well_repeat": well,
        "timepoint": timepoint,
        "well_id": f"{aliquot}::{carb}::{well}",
    }


def parse_scfa_sample_id(sample_id: str) -> dict | None:
    s = str(sample_id).strip().lower()
    m = re.match(
        r"^(?P<root>\d+)(?P<aliq>[a-z])\s+(?P<carb>r\d|s\d|nc)\s+(?P<time>\d+)h$",
        s,
    )
    if not m:
        return None
    root = m.group("root")
    aliq = m.group("aliq").upper()
    carb_tok = m.group("carb")
    timepoint = f"{int(m.group('time'))}H"
    if carb_tok.startswith("r"):
        carb, well = "RDC", int(carb_tok[1:])
    elif carb_tok.startswith("s"):
        carb, well = "SDC", int(carb_tok[1:])
    else:
        carb, well = "No Added Carb", 1
    aliquot = f"{root}{aliq}"
    return {
        "scfa_sample_id": s,
        "donor_id": root,
        "aliquot_id": aliquot,
        "carbohydrate": carb,
        "well_repeat": well,
        "timepoint": timepoint,
        "well_id": f"{aliquot}::{carb}::{well}",
    }


def load_group_map() -> dict[str, str]:
    canon = pd.read_csv(CANON_SCFA, dtype=str)
    g = (
        canon.drop_duplicates("donor_id")
        .set_index("donor_id")["group"]
        .to_dict()
    )
    return g


def load_scfa_long() -> pd.DataFrame:
    raw = pd.read_csv(SCFA_CSV)
    raw["sample_id"] = raw["sampleid"].astype(str).str.strip()
    raw = raw.loc[
        ~raw["sample_id"].str.contains(
            "QC|Blank|blank|Standard|pool|Pool", case=False, na=False
        )
    ]
    parsed = raw["sample_id"].map(parse_scfa_sample_id)
    raw = raw.loc[parsed.notna()].copy()
    meta = pd.DataFrame(list(raw.loc[parsed.notna(), "sample_id"].map(parse_scfa_sample_id)))
    raw = raw.reset_index(drop=True)
    meta = meta.reset_index(drop=True)
    df = pd.concat([raw, meta], axis=1)
    analytes = {
        "acetate": "Acetate",
        "propionate": "Propionate",
        "butyrate": "Butyrate",
    }
    # Column names in export may vary in case
    cols = {c.lower(): c for c in df.columns}
    long_rows = []
    for key, label in analytes.items():
        src = None
        for cand in (label, label.lower(), key, key.capitalize()):
            if cand.lower() in cols:
                src = cols[cand.lower()]
                break
        if src is None:
            # fuzzy
            matches = [c for c in df.columns if key in c.lower()]
            if matches:
                src = matches[0]
        if src is None:
            raise KeyError(f"Could not find SCFA column for {key}")
        part = df[
            [
                "scfa_sample_id",
                "donor_id",
                "aliquot_id",
                "carbohydrate",
                "well_repeat",
                "timepoint",
                "well_id",
                src,
            ]
        ].copy()
        part = part.rename(columns={src: "concentration"})
        part["analyte"] = key
        long_rows.append(part)
    long = pd.concat(long_rows, ignore_index=True)
    wide = long.pivot_table(
        index=[
            "scfa_sample_id",
            "donor_id",
            "aliquot_id",
            "carbohydrate",
            "well_repeat",
            "timepoint",
            "well_id",
        ],
        columns="analyte",
        values="concentration",
        aggfunc="mean",
    ).reset_index()
    wide.columns.name = None
    for a in ("acetate", "propionate", "butyrate"):
        if a not in wide.columns:
            wide[a] = np.nan
    wide["total_scfa"] = wide[["acetate", "propionate", "butyrate"]].sum(
        axis=1, min_count=1
    )
    return wide


def load_qpcr() -> pd.DataFrame:
    q = pd.read_csv(QPCR_CSV)
    q.columns = [c.strip() for c in q.columns]
    gene_col = [c for c in q.columns if c.lower().startswith("gene_copies")][0]
    genome_col = [c for c in q.columns if "genome_copies" in c.lower()][0]
    dna_col = [c for c in q.columns if "dna" in c.lower()][0]
    q["customer_label"] = q["customer_label"].map(clean_label)
    q = q.rename(
        columns={
            gene_col: "gene_copies_per_ul",
            genome_col: "genome_copies_per_ul_star",
            dna_col: "dna_ng_per_ul",
        }
    )
    parsed = q["customer_label"].map(parse_micro_label)
    keep = parsed.notna()
    out = pd.DataFrame([x for x in parsed if x is not None])
    meta = q.loc[keep, ["customer_label", "Ct", "gene_copies_per_ul", "genome_copies_per_ul_star", "dna_ng_per_ul", "sample_id"]].reset_index(drop=True)
    out = out.reset_index(drop=True)
    return pd.concat([out, meta.drop(columns=["customer_label"])], axis=1)


def load_genus_rel() -> pd.DataFrame:
    g = pd.read_csv(GENUS_TSV, sep="\t")
    fusi = g.loc[g["Taxon"].astype(str).str.contains("Fusicatenibacter", case=False, na=False)]
    if fusi.empty:
        raise ValueError("Fusicatenibacter not found in genus table")
    rel = fusi.drop(columns=["Taxon"]).sum(axis=0)
    bif = g.loc[g["Taxon"].astype(str).str.contains(r"^g__Bifidobacterium$", case=False, na=False, regex=True)]
    bif_rel = bif.drop(columns=["Taxon"]).sum(axis=0) if len(bif) else pd.Series(0.0, index=rel.index)
    out = pd.DataFrame(
        {
            "customer_label": [clean_label(x) for x in rel.index],
            "fusicatenibacter_rel": rel.values,
            "bifidobacterium_rel": bif_rel.reindex(rel.index).fillna(0).values,
        }
    )
    return out


def load_asv_depth() -> pd.DataFrame:
    # ASV table: first column is sample label (`seqs`), remaining columns are ASVs.
    asv = pd.read_csv(ASV_CSV)
    sample_col = asv.columns[0]
    depth = asv.drop(columns=[sample_col]).sum(axis=1)
    return pd.DataFrame(
        {
            "customer_label": [clean_label(x) for x in asv[sample_col]],
            "asv_read_depth": depth.values.astype(float),
        }
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    group_map = load_group_map()
    scfa = load_scfa_long()
    qpcr = load_qpcr()
    genus = load_genus_rel()
    depth = load_asv_depth()

    micro = qpcr.merge(genus, on="customer_label", how="outer", indicator="qpcr_genus")
    micro = micro.merge(depth, on="customer_label", how="left")
    # re-parse where qpcr missing but genus present
    for idx, row in micro.loc[micro["donor_id"].isna()].iterrows():
        p = parse_micro_label(row["customer_label"])
        if p:
            for k, v in p.items():
                micro.at[idx, k] = v

    key = ["donor_id", "aliquot_id", "carbohydrate", "well_repeat", "timepoint"]
    matched = scfa.merge(
        micro,
        on=key,
        how="outer",
        suffixes=("_scfa", "_micro"),
        indicator="scfa_micro",
    )
    # unify well_id
    matched["well_id"] = matched["well_id_scfa"].fillna(matched.get("well_id_micro"))
    if "well_id_micro" in matched.columns:
        matched = matched.drop(columns=["well_id_micro"], errors="ignore")
    matched = matched.drop(columns=["well_id_scfa"], errors="ignore")

    matched["group"] = matched["donor_id"].map(group_map)
    matched["group_label"] = matched["group"].map(
        {"control": "healthy-weight", "case": "obesity"}
    )
    matched["fusicatenibacter_qpcr_scaled"] = (
        matched["fusicatenibacter_rel"] * matched["gene_copies_per_ul"]
    )
    matched["bifidobacterium_qpcr_scaled"] = (
        matched["bifidobacterium_rel"] * matched["gene_copies_per_ul"]
    )
    matched["genome_over_gene_ratio"] = (
        matched["genome_copies_per_ul_star"] / matched["gene_copies_per_ul"]
    )

    # diagnostics
    join_counts = (
        matched["scfa_micro"]
        .value_counts()
        .rename_axis("join")
        .reset_index(name="n")
    )
    join_counts.to_csv(AUDIT / "join_counts.csv", index=False)

    unmatched_scfa = matched.loc[matched["scfa_micro"] == "left_only"].copy()
    unmatched_micro = matched.loc[matched["scfa_micro"] == "right_only"].copy()
    unmatched_scfa.to_csv(AUDIT / "unmatched_scfa_rows.csv", index=False)
    unmatched_micro.to_csv(AUDIT / "unmatched_micro_rows.csv", index=False)

    both = matched.loc[matched["scfa_micro"] == "both"].copy()
    dup = both.duplicated(key + ["scfa_sample_id", "customer_label"], keep=False)
    both.loc[dup].to_csv(AUDIT / "duplicate_join_keys.csv", index=False)

    # trajectory completeness for both assays
    traj = (
        both.groupby(
            ["donor_id", "aliquot_id", "carbohydrate", "well_repeat"], as_index=False
        )
        .agg(
            has_0=("timepoint", lambda s: "0H" in set(s)),
            has_48=("timepoint", lambda s: "48H" in set(s)),
            n_rows=("timepoint", "size"),
            n_qpcr=("gene_copies_per_ul", lambda s: s.notna().sum()),
            n_fusi=("fusicatenibacter_qpcr_scaled", lambda s: s.notna().sum()),
        )
    )
    traj["paired_0_48"] = traj["has_0"] & traj["has_48"]
    traj.to_csv(AUDIT / "trajectory_completeness.csv", index=False)

    # missing planned from Phase 0
    planned_missing = pd.read_csv(
        ROOT / "integrated" / "metadata" / "phase0_scfa_missing_planned_donors.csv"
    )
    missing_design = pd.read_csv(
        ROOT / "integrated" / "metadata" / "phase0_scfa_missing_design_samples.csv"
    )
    inventory = pd.DataFrame(
        [
            {"metric": "scfa_rows", "value": len(scfa)},
            {"metric": "qpcr_study_rows", "value": int(qpcr["donor_id"].notna().sum())},
            {"metric": "matched_both", "value": int((matched["scfa_micro"] == "both").sum())},
            {"metric": "unmatched_scfa", "value": int((matched["scfa_micro"] == "left_only").sum())},
            {"metric": "unmatched_micro", "value": int((matched["scfa_micro"] == "right_only").sum())},
            {"metric": "donors_matched", "value": both["donor_id"].nunique()},
            {"metric": "aliquots_matched", "value": both["aliquot_id"].nunique()},
            {
                "metric": "paired_trajectories",
                "value": int(traj["paired_0_48"].sum()),
            },
            {
                "metric": "missing_planned_donors",
                "value": len(planned_missing),
            },
            {
                "metric": "missing_scfa_design_samples",
                "value": len(missing_design),
            },
            {
                "metric": "genome_gene_ratio_median",
                "value": float(both["genome_over_gene_ratio"].median(skipna=True)),
            },
            {
                "metric": "qpcr_qc_documented",
                "value": 0,
            },
            {
                "metric": "qpcr_assay_block_note",
                "value": "LOD/LOQ/primer/efficiency not in repo; provisional thresholds only",
            },
        ]
    )
    inventory.to_csv(AUDIT / "phaseA_inventory.csv", index=False)
    planned_missing.to_csv(AUDIT / "missing_planned_donors.csv", index=False)
    missing_design.to_csv(AUDIT / "missing_scfa_design_samples.csv", index=False)

    # analysis table: matched both preferred; keep all with keys for audit
    analysis_cols = [
        "donor_id",
        "aliquot_id",
        "carbohydrate",
        "well_repeat",
        "timepoint",
        "well_id",
        "group",
        "group_label",
        "scfa_sample_id",
        "customer_label",
        "acetate",
        "propionate",
        "butyrate",
        "total_scfa",
        "asv_read_depth",
        "fusicatenibacter_rel",
        "bifidobacterium_rel",
        "Ct",
        "gene_copies_per_ul",
        "genome_copies_per_ul_star",
        "dna_ng_per_ul",
        "fusicatenibacter_qpcr_scaled",
        "bifidobacterium_qpcr_scaled",
        "scfa_micro",
    ]
    for c in analysis_cols:
        if c not in matched.columns:
            matched[c] = np.nan
    analysis = matched[analysis_cols].sort_values(
        ["donor_id", "aliquot_id", "carbohydrate", "well_repeat", "timepoint"]
    )
    analysis.to_csv(TABLES / "matched_scfa_16s_qpcr.csv", index=False)
    both[analysis_cols].to_csv(TABLES / "matched_scfa_16s_qpcr_both.csv", index=False)

    print(inventory.to_string(index=False))
    print(f"Wrote {TABLES / 'matched_scfa_16s_qpcr.csv'}")
    print(f"Wrote audit under {AUDIT}")


if __name__ == "__main__":
    main()
