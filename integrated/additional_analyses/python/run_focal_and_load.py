#!/usr/bin/env python3
"""Total bacterial load models + focal Fusicatenibacter–SCFA analyses."""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "integrated" / "results" / "additional_analyses_2026-07-19"
TABLES = OUT / "tables"
FOCAL = OUT / "focal_fusicatenibacter"
LOAD = OUT / "total_load"
DIAG = OUT / "diagnostics"

DETECTION_FLOOR_REL = 1e-6  # relative abundance floor before scaling
BOOT_N = 2000
RNG = np.random.default_rng(20260719)


def log10_floor(x: pd.Series | np.ndarray, floor: float) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)
    ok = np.isfinite(arr) & (arr > 0)
    floored = np.maximum(arr[ok], floor)
    out[ok] = np.log10(floored)
    # zeros: assign floor
    zero = np.isfinite(arr) & (arr <= 0)
    out[zero] = np.log10(floor)
    return out


def aggregate_changes(matched: pd.DataFrame) -> pd.DataFrame:
    """Exact-well changes, then mean across wells within aliquot, then A/B within donor."""
    cols = [
        "acetate",
        "propionate",
        "butyrate",
        "total_scfa",
        "gene_copies_per_ul",
        "fusicatenibacter_rel",
        "fusicatenibacter_qpcr_scaled",
        "bifidobacterium_qpcr_scaled",
    ]
    both = matched.loc[matched["scfa_micro"] == "both"].copy()
    # detection floor on qPCR-scaled from rel floor * gene copies median per sample? Use absolute floor.
    # Prefer: floor = DETECTION_FLOOR_REL * gene_copies when available, else small absolute.
    both["fusi_floor"] = np.where(
        both["gene_copies_per_ul"].notna(),
        DETECTION_FLOOR_REL * both["gene_copies_per_ul"],
        1.0,
    )
    both["log10_fusi_qpcr"] = [
        log10_floor([v], max(f, 1e-12))[0]
        for v, f in zip(both["fusicatenibacter_qpcr_scaled"], both["fusi_floor"])
    ]
    both["log10_fusi_rel"] = log10_floor(both["fusicatenibacter_rel"], DETECTION_FLOOR_REL)
    both["log10_load"] = log10_floor(both["gene_copies_per_ul"], 1.0)

    well_id = ["donor_id", "aliquot_id", "carbohydrate", "well_repeat", "group", "group_label"]
    wide_parts = []
    value_cols = {
        "butyrate": "butyrate",
        "acetate": "acetate",
        "propionate": "propionate",
        "total_scfa": "total_scfa",
        "log10_fusi_qpcr": "log10_fusi_qpcr",
        "log10_fusi_rel": "log10_fusi_rel",
        "log10_load": "log10_load",
        "fusicatenibacter_qpcr_scaled": "fusi_qpcr",
        "fusicatenibacter_rel": "fusi_rel",
    }
    for tp, suffix in (("0H", "0"), ("48H", "48")):
        sub = both.loc[both["timepoint"] == tp, well_id + list(value_cols.keys())].copy()
        ren = {k: f"{v}_{suffix}" for k, v in value_cols.items()}
        sub = sub.rename(columns=ren)
        wide_parts.append(sub)
    well = wide_parts[0].merge(wide_parts[1], on=well_id, how="inner")
    for base in (
        "butyrate",
        "acetate",
        "propionate",
        "total_scfa",
        "log10_fusi_qpcr",
        "log10_fusi_rel",
        "log10_load",
        "fusi_qpcr",
        "fusi_rel",
    ):
        well[f"delta_{base}"] = well[f"{base}_48"] - well[f"{base}_0"]

    # aliquot mean across wells
    aliq = (
        well.groupby(
            ["donor_id", "aliquot_id", "carbohydrate", "group", "group_label"],
            as_index=False,
        )
        .mean(numeric_only=True)
    )
    # donor mean across A/B
    donor = (
        aliq.groupby(["donor_id", "carbohydrate", "group", "group_label"], as_index=False)
        .mean(numeric_only=True)
    )
    return well, aliq, donor


def net_against_nc(donor: pd.DataFrame, value: str) -> pd.DataFrame:
    piv = donor.pivot(index=["donor_id", "group", "group_label"], columns="carbohydrate", values=value)
    out = piv.reset_index()
    for carb in ("RDC", "SDC"):
        if carb in out.columns and "No Added Carb" in out.columns:
            out[f"net_{carb}"] = out[carb] - out["No Added Carb"]
        else:
            out[f"net_{carb}"] = np.nan
    return out


def spearman_perm(x, y, n_perm=10000, rng=RNG):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 4:
        return np.nan, np.nan, n
    rho, _ = stats.spearmanr(x, y)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stats.spearmanr(x, rng.permutation(y)).correlation
    p = (np.sum(np.abs(null) >= abs(rho)) + 1) / (n_perm + 1)
    return float(rho), float(p), n


def donor_bootstrap_spearman(x, y, n_boot=BOOT_N, rng=RNG):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 4:
        return np.nan, np.nan, np.nan
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r = stats.spearmanr(x[idx], y[idx]).correlation
        if np.isfinite(r):
            boots.append(r)
    boots = np.asarray(boots)
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)), len(boots)


def leave_one_out(x, y, labels):
    rows = []
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    labels = np.asarray(labels)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y, labels = x[mask], y[mask], labels[mask]
    for i, lab in enumerate(labels):
        keep = np.ones(len(x), dtype=bool)
        keep[i] = False
        if keep.sum() < 4:
            continue
        rho = stats.spearmanr(x[keep], y[keep]).correlation
        rows.append({"left_out_donor": lab, "spearman_rho": rho, "n": int(keep.sum())})
    return pd.DataFrame(rows)


def fit_load_models(matched: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    both = matched.loc[
        (matched["scfa_micro"] == "both") & matched["gene_copies_per_ul"].notna()
    ].copy()
    both["log10_load"] = log10_floor(both["gene_copies_per_ul"], 1.0)
    both["time_48"] = (both["timepoint"] == "48H").astype(int)
    both["carb"] = pd.Categorical(
        both["carbohydrate"], categories=["No Added Carb", "RDC", "SDC"]
    )
    both["group_f"] = pd.Categorical(both["group"], categories=["control", "case"])

    # Aggregate to aliquot×carb×time to stabilize nesting (mean across wells)
    agg = (
        both.groupby(
            [
                "donor_id",
                "aliquot_id",
                "carbohydrate",
                "timepoint",
                "group",
            ],
            as_index=False,
            observed=True,
        )
        .agg(log10_load=("log10_load", "mean"), time_48=("time_48", "max"))
    )
    agg["carb"] = pd.Categorical(
        agg["carbohydrate"], categories=["No Added Carb", "RDC", "SDC"]
    )
    agg["group_f"] = pd.Categorical(agg["group"], categories=["control", "case"])

    # donor-level change model
    d0 = agg.loc[agg["timepoint"] == "0H"]
    d48 = agg.loc[agg["timepoint"] == "48H"]
    ch = d0.merge(
        d48,
        on=["donor_id", "aliquot_id", "carbohydrate", "group"],
        suffixes=("_0", "_48"),
    )
    ch["delta_log10_load"] = ch["log10_load_48"] - ch["log10_load_0"]
    ch["carb"] = pd.Categorical(
        ch["carbohydrate"], categories=["No Added Carb", "RDC", "SDC"]
    )
    ch["group_f"] = pd.Categorical(ch["group"], categories=["control", "case"])
    donor_ch = (
        ch.groupby(["donor_id", "carbohydrate", "group"], as_index=False, observed=True)
        .agg(delta_log10_load=("delta_log10_load", "mean"))
    )
    donor_ch["carb"] = pd.Categorical(
        donor_ch["carbohydrate"], categories=["No Added Carb", "RDC", "SDC"]
    )
    donor_ch["group_f"] = pd.Categorical(
        donor_ch["group"], categories=["control", "case"]
    )

    # Mixed model on aliquot-level deltas
    model_rows = []
    try:
        md = MixedLM.from_formula(
            "delta_log10_load ~ C(carb, Treatment(reference='No Added Carb')) * C(group_f, Treatment(reference='control'))",
            groups="donor_id",
            data=ch,
        )
        fit = md.fit(reml=True, method="lbfgs")
        for name, est, se, p in zip(
            fit.fe_params.index, fit.fe_params.values, fit.bse_fe.values, fit.pvalues.values
        ):
            model_rows.append(
                {
                    "model": "aliquot_delta_MixedLM",
                    "term": name,
                    "estimate": est,
                    "se": se,
                    "pvalue": p,
                    "n_obs": len(ch),
                    "n_donors": ch["donor_id"].nunique(),
                    "converged": bool(fit.converged),
                }
            )
    except Exception as e:
        model_rows.append(
            {
                "model": "aliquot_delta_MixedLM",
                "term": "ERROR",
                "estimate": np.nan,
                "se": np.nan,
                "pvalue": np.nan,
                "n_obs": len(ch),
                "n_donors": ch["donor_id"].nunique(),
                "converged": False,
                "note": str(e),
            }
        )

    # Simple donor-level pairwise contrasts vs NC
    contrast_rows = []
    for carb in ("RDC", "SDC"):
        sub = donor_ch.loc[donor_ch["carbohydrate"].isin([carb, "No Added Carb"])]
        piv = sub.pivot(index="donor_id", columns="carbohydrate", values="delta_log10_load")
        if carb in piv.columns and "No Added Carb" in piv.columns:
            d = (piv[carb] - piv["No Added Carb"]).dropna()
            if len(d) >= 3:
                t = stats.ttest_1samp(d, 0.0)
                contrast_rows.append(
                    {
                        "contrast": f"{carb}_vs_NC_delta_log10_load",
                        "mean_diff": float(d.mean()),
                        "sd": float(d.std(ddof=1)),
                        "n_donors": int(d.shape[0]),
                        "t_stat": float(t.statistic),
                        "pvalue": float(t.pvalue),
                        "ci95_low": float(d.mean() - 1.96 * d.sem()),
                        "ci95_high": float(d.mean() + 1.96 * d.sem()),
                    }
                )
    # SDC vs RDC
    piv = donor_ch.pivot(index="donor_id", columns="carbohydrate", values="delta_log10_load")
    if {"RDC", "SDC"}.issubset(piv.columns):
        d = (piv["SDC"] - piv["RDC"]).dropna()
        t = stats.ttest_1samp(d, 0.0)
        contrast_rows.append(
            {
                "contrast": "SDC_vs_RDC_delta_log10_load",
                "mean_diff": float(d.mean()),
                "sd": float(d.std(ddof=1)),
                "n_donors": int(d.shape[0]),
                "t_stat": float(t.statistic),
                "pvalue": float(t.pvalue),
                "ci95_low": float(d.mean() - 1.96 * d.sem()),
                "ci95_high": float(d.mean() + 1.96 * d.sem()),
            }
        )

    return pd.DataFrame(model_rows), pd.DataFrame(contrast_rows), donor_ch, ch


def run_focal(donor: pd.DataFrame, aliq: pd.DataFrame) -> None:
    # Primary: net SDC Δlog10 fusi qpcr vs net SDC Δbutyrate
    net_but = net_against_nc(donor, "delta_butyrate").rename(
        columns={"net_SDC": "net_delta_butyrate_SDC", "net_RDC": "net_delta_butyrate_RDC"}
    )
    net_fusi = net_against_nc(donor, "delta_log10_fusi_qpcr").rename(
        columns={
            "net_SDC": "net_delta_log10_fusi_qpcr_SDC",
            "net_RDC": "net_delta_log10_fusi_qpcr_RDC",
        }
    )
    net_fusi_rel = net_against_nc(donor, "delta_log10_fusi_rel").rename(
        columns={
            "net_SDC": "net_delta_log10_fusi_rel_SDC",
            "net_RDC": "net_delta_log10_fusi_rel_RDC",
        }
    )
    net_load = net_against_nc(donor, "delta_log10_load").rename(
        columns={"net_SDC": "net_delta_log10_load_SDC", "net_RDC": "net_delta_log10_load_RDC"}
    )
    net_ace = net_against_nc(donor, "delta_acetate")
    net_pro = net_against_nc(donor, "delta_propionate")
    net_tot = net_against_nc(donor, "delta_total_scfa")

    primary = (
        net_but[["donor_id", "group", "group_label", "net_delta_butyrate_SDC", "net_delta_butyrate_RDC"]]
        .merge(
            net_fusi[
                [
                    "donor_id",
                    "net_delta_log10_fusi_qpcr_SDC",
                    "net_delta_log10_fusi_qpcr_RDC",
                ]
            ],
            on="donor_id",
        )
        .merge(
            net_fusi_rel[
                [
                    "donor_id",
                    "net_delta_log10_fusi_rel_SDC",
                    "net_delta_log10_fusi_rel_RDC",
                ]
            ],
            on="donor_id",
        )
        .merge(
            net_load[["donor_id", "net_delta_log10_load_SDC", "net_delta_log10_load_RDC"]],
            on="donor_id",
        )
    )
    primary.to_csv(FOCAL / "donor_net_change_table.csv", index=False)

    # Primary test
    x = primary["net_delta_log10_fusi_qpcr_SDC"]
    y = primary["net_delta_butyrate_SDC"]
    rho, p, n = spearman_perm(x, y)
    lo, hi, nboot = donor_bootstrap_spearman(x, y)
    # also Pearson for slope interpretability at donor level after OLS
    mask = x.notna() & y.notna()
    slope = intercept = r2 = np.nan
    if mask.sum() >= 4:
        lr = stats.linregress(x[mask], y[mask])
        slope, intercept, r2 = lr.slope, lr.intercept, lr.rvalue**2

    primary_res = pd.DataFrame(
        [
            {
                "analysis": "primary_donor_net_SDC",
                "exposure": "net_delta_log10_fusi_qpcr_SDC",
                "outcome": "net_delta_butyrate_SDC",
                "spearman_rho": rho,
                "permutation_p": p,
                "bootstrap_ci95_low": lo,
                "bootstrap_ci95_high": hi,
                "n_donors": n,
                "ols_slope_per_log10": slope,
                "ols_intercept": intercept,
                "ols_r2": r2,
                "n_bootstrap": nboot,
            }
        ]
    )
    primary_res.to_csv(FOCAL / "primary_fusicatenibacter_butyrate.csv", index=False)
    loo = leave_one_out(x, y, primary["donor_id"])
    loo.to_csv(FOCAL / "primary_leave_one_donor_out.csv", index=False)

    # Specificity family
    specs = []
    pairs = [
        ("net_delta_log10_fusi_qpcr_SDC", "net_delta_butyrate_SDC", "primary_SDC_butyrate"),
        ("net_delta_log10_fusi_qpcr_RDC", "net_delta_butyrate_RDC", "RDC_butyrate"),
        ("net_delta_log10_fusi_rel_SDC", "net_delta_butyrate_SDC", "SDC_rel_fusi_butyrate"),
        ("net_delta_log10_load_SDC", "net_delta_butyrate_SDC", "SDC_total_load_butyrate"),
    ]
    # acetate/propionate/total need merge
    primary2 = primary.merge(
        net_ace[["donor_id", "net_SDC"]].rename(columns={"net_SDC": "net_delta_acetate_SDC"}),
        on="donor_id",
    ).merge(
        net_pro[["donor_id", "net_SDC"]].rename(columns={"net_SDC": "net_delta_propionate_SDC"}),
        on="donor_id",
    ).merge(
        net_tot[["donor_id", "net_SDC"]].rename(columns={"net_SDC": "net_delta_total_scfa_SDC"}),
        on="donor_id",
    )
    pairs.extend(
        [
            ("net_delta_log10_fusi_qpcr_SDC", "net_delta_acetate_SDC", "SDC_acetate"),
            ("net_delta_log10_fusi_qpcr_SDC", "net_delta_propionate_SDC", "SDC_propionate"),
            ("net_delta_log10_fusi_qpcr_SDC", "net_delta_total_scfa_SDC", "SDC_total_scfa"),
        ]
    )
    for ex, ey, name in pairs:
        rho, p, n = spearman_perm(primary2[ex], primary2[ey])
        lo, hi, _ = donor_bootstrap_spearman(primary2[ex], primary2[ey])
        specs.append(
            {
                "analysis": name,
                "exposure": ex,
                "outcome": ey,
                "spearman_rho": rho,
                "permutation_p": p,
                "bootstrap_ci95_low": lo,
                "bootstrap_ci95_high": hi,
                "n_donors": n,
            }
        )
    spec_df = pd.DataFrame(specs)
    # BH within secondary family (exclude primary row if present as first)
    secondary = spec_df.loc[spec_df["analysis"] != "primary_SDC_butyrate"].copy()
    if len(secondary):
        from statsmodels.stats.multitest import multipletests

        reject, q, _, _ = multipletests(secondary["permutation_p"].fillna(1), method="fdr_bh")
        secondary["bh_q"] = q
        secondary["bh_reject_0.1"] = reject
    secondary.to_csv(FOCAL / "specificity_family.csv", index=False)

    # Aliquot-level MixedLM sensitivity
    # Build aliquot net SDC
    def aliq_net(df, value):
        piv = df.pivot_table(
            index=["donor_id", "aliquot_id", "group"],
            columns="carbohydrate",
            values=value,
            aggfunc="mean",
        ).reset_index()
        if "SDC" in piv.columns and "No Added Carb" in piv.columns:
            piv["net"] = piv["SDC"] - piv["No Added Carb"]
        else:
            piv["net"] = np.nan
        return piv[["donor_id", "aliquot_id", "group", "net"]]

    a_but = aliq_net(aliq, "delta_butyrate").rename(columns={"net": "net_butyrate"})
    a_fusi = aliq_net(aliq, "delta_log10_fusi_qpcr").rename(columns={"net": "net_fusi"})
    a = a_but.merge(a_fusi, on=["donor_id", "aliquot_id", "group"])
    a = a.dropna(subset=["net_butyrate", "net_fusi"])
    a.to_csv(FOCAL / "aliquot_net_change_table.csv", index=False)
    mix_rows = []
    try:
        md = MixedLM.from_formula(
            "net_butyrate ~ net_fusi + C(group)",
            groups="donor_id",
            data=a,
        )
        fit = md.fit(reml=True, method="lbfgs")
        for name, est, se, p in zip(
            fit.fe_params.index, fit.fe_params.values, fit.bse_fe.values, fit.pvalues.values
        ):
            mix_rows.append(
                {
                    "term": name,
                    "estimate": est,
                    "se": se,
                    "pvalue": p,
                    "n_obs": len(a),
                    "n_donors": a["donor_id"].nunique(),
                    "converged": bool(fit.converged),
                }
            )
    except Exception as e:
        mix_rows.append(
            {
                "term": "ERROR",
                "estimate": np.nan,
                "se": np.nan,
                "pvalue": np.nan,
                "n_obs": len(a),
                "n_donors": a["donor_id"].nunique(),
                "converged": False,
                "note": str(e),
            }
        )
    pd.DataFrame(mix_rows).to_csv(FOCAL / "aliquot_mixedlm_sensitivity.csv", index=False)

    # Raw SDC change (not net) sensitivity
    raw = donor.loc[donor["carbohydrate"] == "SDC"].copy()
    rho, p, n = spearman_perm(raw["delta_log10_fusi_qpcr"], raw["delta_butyrate"])
    lo, hi, _ = donor_bootstrap_spearman(raw["delta_log10_fusi_qpcr"], raw["delta_butyrate"])
    pd.DataFrame(
        [
            {
                "analysis": "raw_SDC_change_not_net",
                "spearman_rho": rho,
                "permutation_p": p,
                "bootstrap_ci95_low": lo,
                "bootstrap_ci95_high": hi,
                "n_donors": n,
            }
        ]
    ).to_csv(FOCAL / "sensitivity_raw_sdc_change.csv", index=False)


def variance_components(matched: pd.DataFrame) -> pd.DataFrame:
    """Rough ICC from nested ANOVA-style variance of deltas at well level."""
    both = matched.loc[matched["scfa_micro"] == "both"].copy()
    # Well-level 0-48 deltas for butyrate and load
    rows = []
    for analyte_col, label in (
        ("butyrate", "butyrate"),
        ("gene_copies_per_ul", "gene_copies_per_ul"),
    ):
        w0 = both.loc[both["timepoint"] == "0H"]
        w48 = both.loc[both["timepoint"] == "48H"]
        m = w0.merge(
            w48,
            on=["donor_id", "aliquot_id", "carbohydrate", "well_repeat"],
            suffixes=("_0", "_48"),
        )
        m["delta"] = m[f"{analyte_col}_48"] - m[f"{analyte_col}_0"]
        if label == "gene_copies_per_ul":
            m["delta"] = np.log10(np.maximum(m[f"{analyte_col}_48"], 1)) - np.log10(
                np.maximum(m[f"{analyte_col}_0"], 1)
            )
        for carb, sub in m.groupby("carbohydrate"):
            # variance partition using group means
            donor_mean = sub.groupby("donor_id")["delta"].transform("mean")
            aliq_mean = sub.groupby(["donor_id", "aliquot_id"])["delta"].transform("mean")
            var_total = sub["delta"].var(ddof=1)
            var_donor = donor_mean.var(ddof=1)
            var_aliq = (aliq_mean - donor_mean).var(ddof=1)
            var_well = (sub["delta"] - aliq_mean).var(ddof=1)
            rows.append(
                {
                    "outcome": f"delta_{label}",
                    "carbohydrate": carb,
                    "var_total": var_total,
                    "var_donor": var_donor,
                    "var_aliquot_within_donor": var_aliq,
                    "var_well_within_aliquot": var_well,
                    "icc_donor": var_donor / var_total if var_total else np.nan,
                    "n_wells": len(sub),
                    "n_donors": sub["donor_id"].nunique(),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    warnings.filterwarnings("ignore")
    FOCAL.mkdir(parents=True, exist_ok=True)
    LOAD.mkdir(parents=True, exist_ok=True)
    DIAG.mkdir(parents=True, exist_ok=True)

    matched = pd.read_csv(TABLES / "matched_scfa_16s_qpcr.csv", dtype={"donor_id": str})
    well, aliq, donor = aggregate_changes(matched)
    well.to_csv(TABLES / "well_level_deltas.csv", index=False)
    aliq.to_csv(TABLES / "aliquot_level_deltas.csv", index=False)
    donor.to_csv(TABLES / "donor_level_deltas.csv", index=False)

    run_focal(donor, aliq)

    mix, contrasts, donor_ch, aliq_ch = fit_load_models(matched)
    mix.to_csv(LOAD / "mixedlm_terms.csv", index=False)
    contrasts.to_csv(LOAD / "donor_contrasts.csv", index=False)
    donor_ch.to_csv(LOAD / "donor_delta_log10_load.csv", index=False)
    aliq_ch.to_csv(LOAD / "aliquot_delta_log10_load.csv", index=False)

    vc = variance_components(matched)
    vc.to_csv(DIAG / "variance_components_deltas.csv", index=False)

    print("Primary:")
    print(pd.read_csv(FOCAL / "primary_fusicatenibacter_butyrate.csv").to_string(index=False))
    print("\nLoad contrasts:")
    print(contrasts.to_string(index=False))
    print(f"\nWrote focal outputs under {FOCAL}")


if __name__ == "__main__":
    main()
