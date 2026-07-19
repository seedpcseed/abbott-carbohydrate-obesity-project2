"""Obesity-group SCFA trajectory contrasts and equivalence evaluation.

The primary estimand is the obesity-minus-healthy-weight difference in the
0-to-48-hour concentration change, separately for RDC and SDC, for acetate,
propionate, and butyrate. A/B is retained as part of the subject identifier;
R1/R2 and S1/S2 identify biological culture wells.

Formal TOST decisions are made only when an externally justified absolute
margin is populated in ``equivalence_margins.csv``. When no margin is
available, the script reports 90% and 95% confidence intervals and the
post-hoc minimum bound that would contain the 90% interval, without calling
the groups equivalent.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import het_breuschpagan


ROOT = Path(__file__).resolve().parents[1]
SCFA_DIR = ROOT / "scfa_metabolomics"
DATA_FILE = (
    SCFA_DIR
    / "data"
    / "removed_qcs_quant_results_20250722_PFBBr_PS2720_20250731.csv"
)
METADATA_FILES = (
    SCFA_DIR / "metadata" / "SEED LAB HMM Submission Sheet 202405.xlsx",
    SCFA_DIR / "metadata" / "PS2720_xtra_samples_HMM Submission Sheet.xlsx",
)
MARGIN_FILE = SCFA_DIR / "equivalence_margins.csv"
RESULTS_DIR = SCFA_DIR / "results"
FIGURES_DIR = SCFA_DIR / "figures"

ANALYTES = (
    "acetate",
    "propionate",
    "butyrate",
    "5aminovalerate",
    "succinate",
)
PRIMARY_ANALYTES = ("acetate", "propionate", "butyrate")
CARBOHYDRATES = ("No added carbohydrate", "RDC", "SDC")
PRIMARY_CARBOHYDRATES = ("RDC", "SDC")
GROUPS = ("control", "case")
BOOTSTRAP_REPLICATES = 10_000
RANDOM_SEED = 2720

FORMULA = (
    "concentration ~ "
    "C(group, Treatment(reference='control')) * "
    "C(carbohydrate, Treatment(reference='No added carbohydrate')) * "
    "C(timepoint_hr, Treatment(reference=0))"
)


@dataclass
class ModelFit:
    result: Any
    model_type: str
    optimizer: str
    warning_text: str


def load_metadata() -> pd.DataFrame:
    """Read, validate, and deduplicate the two facility submission workbooks."""
    frames: list[pd.DataFrame] = []
    for path in METADATA_FILES:
        frame = pd.read_excel(path, skiprows=9).iloc[3:, :2].copy()
        frame.columns = ["sampleid", "group"]
        frame["source_file"] = path.name
        frames.append(frame)

    metadata = pd.concat(frames, ignore_index=True).dropna(
        subset=["sampleid", "group"]
    )
    metadata["sampleid"] = metadata["sampleid"].astype(str).str.strip()
    metadata["group"] = metadata["group"].astype(str).str.lower().str.strip()

    inconsistent_groups = (
        metadata.groupby("sampleid", observed=True)["group"].nunique().gt(1)
    )
    if inconsistent_groups.any():
        bad = inconsistent_groups[inconsistent_groups].index.tolist()
        raise ValueError(f"Sample IDs map to multiple groups: {bad[:10]}")

    metadata = metadata.drop_duplicates("sampleid", keep="first").copy()
    metadata["subject"] = (
        metadata["sampleid"]
        .str.extract(r"^(\d+[A-Za-z])", expand=False)
        .str.upper()
    )
    metadata["numeric_root"] = metadata["sampleid"].str.extract(
        r"^(\d+)", expand=False
    )
    metadata["carb_repeat"] = (
        metadata["sampleid"]
        .str.extract(r"(?<= )([rs]\d|nc)(?= )", expand=False)
        .str.lower()
    )
    metadata["carbohydrate"] = np.select(
        (
            metadata["carb_repeat"].str.startswith("r", na=False),
            metadata["carb_repeat"].str.startswith("s", na=False),
            metadata["carb_repeat"].eq("nc"),
        ),
        ("RDC", "SDC", "No added carbohydrate"),
        default=None,
    )
    metadata["well_repeat"] = pd.to_numeric(
        metadata["carb_repeat"].str.extract(r"(\d+)", expand=False),
        errors="coerce",
    ).fillna(1).astype(int)
    metadata["timepoint_hr"] = pd.to_numeric(
        metadata["sampleid"].str.extract(r" (\d+)h$", expand=False),
        errors="coerce",
    )
    metadata["well_id"] = (
        metadata["subject"].astype(str)
        + "::"
        + metadata["carbohydrate"].astype(str)
        + "::"
        + metadata["well_repeat"].astype(str)
    )

    required = [
        "subject",
        "carbohydrate",
        "well_repeat",
        "timepoint_hr",
    ]
    incomplete = metadata[required].isna().any(axis=1)
    metadata = metadata.loc[
        ~incomplete & metadata["group"].isin(GROUPS)
    ].copy()
    return metadata


def load_long_data() -> pd.DataFrame:
    """Join concentrations to metadata while retaining biological wells."""
    metadata = load_metadata()
    concentrations = pd.read_csv(DATA_FILE)
    concentrations["sampleid"] = concentrations["sampleid"].astype(str).str.strip()

    sample_data = metadata.merge(
        concentrations,
        on="sampleid",
        how="inner",
        validate="one_to_one",
    )
    long_data = sample_data.melt(
        id_vars=[
            "sampleid",
            "subject",
            "numeric_root",
            "group",
            "carbohydrate",
            "well_repeat",
            "well_id",
            "timepoint_hr",
        ],
        value_vars=ANALYTES,
        var_name="analyte",
        value_name="concentration",
    ).dropna(subset=["concentration"])

    long_data["carbohydrate"] = pd.Categorical(
        long_data["carbohydrate"],
        categories=CARBOHYDRATES,
        ordered=True,
    )
    long_data["group"] = pd.Categorical(
        long_data["group"], categories=GROUPS, ordered=True
    )
    return long_data


def load_margins() -> pd.DataFrame:
    """Load the frozen external margins; blank margins remain unavailable."""
    margins = pd.read_csv(MARGIN_FILE)
    required = {
        "analyte",
        "margin_lower_uM",
        "margin_upper_uM",
        "status",
        "version",
        "source",
        "rationale",
    }
    missing = required.difference(margins.columns)
    if missing:
        raise ValueError(f"Margin file is missing columns: {sorted(missing)}")
    margins["margin_lower_uM"] = pd.to_numeric(
        margins["margin_lower_uM"], errors="coerce"
    )
    margins["margin_upper_uM"] = pd.to_numeric(
        margins["margin_upper_uM"], errors="coerce"
    )
    return margins


def fit_trajectory_model(analyte_data: pd.DataFrame) -> ModelFit:
    """Fit the nested-well model, with documented fallbacks."""
    attempts = (
        ("mixed_subject_well", "lbfgs"),
        ("mixed_subject_well", "powell"),
        ("mixed_subject", "lbfgs"),
        ("mixed_subject", "powell"),
    )
    messages: list[str] = []

    for model_type, optimizer in attempts:
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                if model_type == "mixed_subject_well":
                    model = smf.mixedlm(
                        FORMULA,
                        analyte_data,
                        groups=analyte_data["subject"],
                        re_formula="1",
                        vc_formula={"well": "0 + C(well_id)"},
                    )
                else:
                    model = smf.mixedlm(
                        FORMULA,
                        analyte_data,
                        groups=analyte_data["subject"],
                        re_formula="1",
                    )
                result = model.fit(
                    reml=True,
                    method=optimizer,
                    maxiter=2_000,
                    disp=False,
                )
                warning_text = " | ".join(
                    sorted({str(item.message) for item in caught})
                )
            if bool(getattr(result, "converged", False)):
                return ModelFit(result, model_type, optimizer, warning_text)
            messages.append(f"{model_type}/{optimizer}: did not converge")
        except Exception as exc:  # pragma: no cover - exercised only on fit failure
            messages.append(f"{model_type}/{optimizer}: {exc}")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = smf.ols(FORMULA, analyte_data)
        result = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": analyte_data["subject"]},
        )
        warning_text = " | ".join(
            sorted({str(item.message) for item in caught})
        )
    messages.append("Fallback to donor-clustered OLS")
    return ModelFit(
        result,
        "ols_cluster_fallback",
        "cluster_robust",
        " | ".join(messages + [warning_text]).strip(" |"),
    )


def fixed_parameters(model_fit: ModelFit) -> tuple[pd.Series, pd.DataFrame]:
    """Return fixed-effect estimates and their covariance matrix."""
    result = model_fit.result
    if model_fit.model_type.startswith("mixed"):
        params = result.fe_params
    else:
        params = result.params
    covariance = result.cov_params().loc[params.index, params.index]
    return params, covariance


def design_row(model_fit: ModelFit, group: str, carbohydrate: str, time: int) -> np.ndarray:
    """Construct one fixed-effect design row from the fitted model schema."""
    design_info = model_fit.result.model.data.design_info
    cell = pd.DataFrame(
        {
            "group": [group],
            "carbohydrate": [carbohydrate],
            "timepoint_hr": [time],
        }
    )
    matrix = patsy.build_design_matrices(
        [design_info], cell, return_type="dataframe"
    )[0]
    params, _ = fixed_parameters(model_fit)
    return matrix.loc[:, params.index].to_numpy()[0]


def difference_in_differences(
    model_fit: ModelFit, carbohydrate: str
) -> dict[str, float]:
    """Estimate (Obesity48-Obesity0) - (Healthy48-Healthy0)."""
    params, covariance = fixed_parameters(model_fit)
    contrast = (
        design_row(model_fit, "case", carbohydrate, 48)
        - design_row(model_fit, "case", carbohydrate, 0)
        - design_row(model_fit, "control", carbohydrate, 48)
        + design_row(model_fit, "control", carbohydrate, 0)
    )
    estimate = float(contrast @ params.to_numpy())
    standard_error = float(
        np.sqrt(contrast @ covariance.to_numpy() @ contrast)
    )
    z_value = estimate / standard_error
    ci90_half = stats.norm.ppf(0.95) * standard_error
    ci95_half = stats.norm.ppf(0.975) * standard_error
    return {
        "estimate_uM": estimate,
        "standard_error_uM": standard_error,
        "degrees_freedom": np.inf,
        "inference_distribution": "asymptotic_normal_wald",
        "wald_z": z_value,
        "superiority_p_two_sided": float(2 * stats.norm.sf(abs(z_value))),
        "ci90_lower_uM": estimate - ci90_half,
        "ci90_upper_uM": estimate + ci90_half,
        "ci95_lower_uM": estimate - ci95_half,
        "ci95_upper_uM": estimate + ci95_half,
    }


def donor_deltas(long_data: pd.DataFrame) -> pd.DataFrame:
    """Average wells within donor-cell-time and calculate 48h minus 0h."""
    averaged = (
        long_data.groupby(
            ["subject", "group", "carbohydrate", "timepoint_hr", "analyte"],
            observed=True,
            as_index=False,
        )["concentration"]
        .mean()
    )
    wide = averaged.pivot(
        index=["subject", "group", "carbohydrate", "analyte"],
        columns="timepoint_hr",
        values="concentration",
    ).reset_index()
    wide["delta_48h_minus_0h_uM"] = wide.get(48) - wide.get(0)
    return wide.dropna(subset=["delta_48h_minus_0h_uM"]).copy()


def bootstrap_delta_contrast(
    deltas: pd.DataFrame,
    analyte: str,
    carbohydrate: str,
    seed: int,
) -> dict[str, float | int]:
    """Donor-level nonparametric bootstrap sensitivity interval."""
    subset = deltas.loc[
        deltas["analyte"].eq(analyte)
        & deltas["carbohydrate"].eq(carbohydrate)
    ]
    control = subset.loc[
        subset["group"].eq("control"), "delta_48h_minus_0h_uM"
    ].to_numpy()
    case = subset.loc[
        subset["group"].eq("case"), "delta_48h_minus_0h_uM"
    ].to_numpy()
    rng = np.random.default_rng(seed)
    control_means = rng.choice(
        control,
        size=(BOOTSTRAP_REPLICATES, len(control)),
        replace=True,
    ).mean(axis=1)
    case_means = rng.choice(
        case,
        size=(BOOTSTRAP_REPLICATES, len(case)),
        replace=True,
    ).mean(axis=1)
    differences = case_means - control_means
    return {
        "n_healthy_weight": len(control),
        "n_obesity": len(case),
        "healthy_weight_mean_delta_uM": float(control.mean()),
        "obesity_mean_delta_uM": float(case.mean()),
        "bootstrap_ci90_lower_uM": float(np.quantile(differences, 0.05)),
        "bootstrap_ci90_upper_uM": float(np.quantile(differences, 0.95)),
        "bootstrap_ci95_lower_uM": float(np.quantile(differences, 0.025)),
        "bootstrap_ci95_upper_uM": float(np.quantile(differences, 0.975)),
    }


def leave_one_donor_influence(
    deltas: pd.DataFrame, analyte: str, carbohydrate: str
) -> dict[str, float | str]:
    """Summarize the largest leave-one-donor change in the raw delta contrast."""
    subset = deltas.loc[
        deltas["analyte"].eq(analyte)
        & deltas["carbohydrate"].eq(carbohydrate)
    ].copy()
    baseline = (
        subset.loc[
            subset["group"].eq("case"), "delta_48h_minus_0h_uM"
        ].mean()
        - subset.loc[
            subset["group"].eq("control"), "delta_48h_minus_0h_uM"
        ].mean()
    )
    changes: list[tuple[str, float]] = []
    for subject in subset["subject"].unique():
        reduced = subset.loc[~subset["subject"].eq(subject)]
        estimate = (
            reduced.loc[
                reduced["group"].eq("case"), "delta_48h_minus_0h_uM"
            ].mean()
            - reduced.loc[
                reduced["group"].eq("control"), "delta_48h_minus_0h_uM"
            ].mean()
        )
        changes.append((str(subject), float(abs(estimate - baseline))))
    influential_subject, maximum_change = max(changes, key=lambda item: item[1])
    return {
        "most_influential_subject": influential_subject,
        "maximum_leave_one_out_change_uM": maximum_change,
    }


def tost_result(
    estimate: float,
    standard_error: float,
    lower_margin: float,
    upper_margin: float,
) -> dict[str, float | bool | str]:
    """Run TOST against fixed lower and upper absolute margins."""
    lower_z = (estimate - lower_margin) / standard_error
    upper_z = (estimate - upper_margin) / standard_error
    lower_p = float(stats.norm.sf(lower_z))
    upper_p = float(stats.norm.cdf(upper_z))
    tost_p = max(lower_p, upper_p)
    return {
        "tost_lower_p": lower_p,
        "tost_upper_p": upper_p,
        "tost_p": tost_p,
        "equivalent_alpha_0_05": bool(tost_p < 0.05),
        "equivalence_verdict": (
            "equivalent_within_prespecified_margin"
            if tost_p < 0.05
            else "equivalence_not_established"
        ),
    }


def unavailable_tost() -> dict[str, float | bool | str]:
    """Represent a deliberately unevaluable equivalence test."""
    return {
        "tost_lower_p": np.nan,
        "tost_upper_p": np.nan,
        "tost_p": np.nan,
        "equivalent_alpha_0_05": np.nan,
        "equivalence_verdict": "not_evaluable_margin_unavailable",
    }


def model_diagnostics(
    analyte: str, analyte_data: pd.DataFrame, model_fit: ModelFit
) -> dict[str, Any]:
    """Extract convergence, variance, and residual diagnostics."""
    result = model_fit.result
    try:
        residuals = np.asarray(result.resid, dtype=float)
    except Exception:
        residuals = np.asarray(result.model.endog - result.fittedvalues, dtype=float)

    exog = np.asarray(result.model.exog, dtype=float)
    try:
        bp_lm, bp_p, _, _ = het_breuschpagan(residuals, exog)
    except Exception:
        bp_lm, bp_p = np.nan, np.nan
    shapiro_p = (
        float(stats.shapiro(residuals).pvalue)
        if 3 <= len(residuals) <= 5_000
        else np.nan
    )

    subject_variance = np.nan
    well_variance = np.nan
    residual_variance = float(getattr(result, "scale", np.nan))
    singular = False
    if model_fit.model_type.startswith("mixed"):
        try:
            subject_variance = float(result.cov_re.iloc[0, 0])
        except Exception:
            pass
        if model_fit.model_type == "mixed_subject_well":
            try:
                well_variance = float(np.asarray(result.vcomp)[0])
            except Exception:
                pass
        variance_values = np.asarray(
            [subject_variance, well_variance], dtype=float
        )
        finite_variances = variance_values[np.isfinite(variance_values)]
        total_variance = residual_variance + finite_variances.sum()
        singular = bool(
            finite_variances.size
            and (
                np.any(
                    finite_variances
                    / max(total_variance, np.finfo(float).eps)
                    < 1e-3
                )
                or "boundary" in model_fit.warning_text.lower()
                or "singular" in model_fit.warning_text.lower()
            )
        )

    return {
        "analyte": analyte,
        "model_type": model_fit.model_type,
        "optimizer": model_fit.optimizer,
        "converged": bool(getattr(result, "converged", True)),
        "singular_or_boundary": singular,
        "n_observations": int(len(analyte_data)),
        "n_subjects": int(analyte_data["subject"].nunique()),
        "n_wells": int(analyte_data["well_id"].nunique()),
        "subject_random_intercept_variance": subject_variance,
        "well_random_intercept_variance": well_variance,
        "residual_variance": residual_variance,
        "residual_skewness": float(stats.skew(residuals, bias=False)),
        "residual_excess_kurtosis": float(
            stats.kurtosis(residuals, fisher=True, bias=False)
        ),
        "shapiro_wilk_p": shapiro_p,
        "breusch_pagan_lm": bp_lm,
        "breusch_pagan_p": bp_p,
        "fit_warnings": model_fit.warning_text,
    }


def sample_counts(long_data: pd.DataFrame) -> pd.DataFrame:
    """Count subjects, wells, and analyte observations by design cell."""
    counts = (
        long_data.groupby(
            ["group", "carbohydrate", "timepoint_hr"],
            observed=True,
            as_index=False,
        )
        .agg(
            n_subjects=("subject", "nunique"),
            n_wells=("well_id", "nunique"),
            n_analyte_observations=("concentration", "size"),
        )
        .sort_values(["group", "carbohydrate", "timepoint_hr"])
    )
    counts["group_label"] = counts["group"].map(
        {"control": "Healthy-weight", "case": "Obesity"}
    )
    counts["n_analyte_observations"] = (
        counts["n_analyte_observations"] // len(ANALYTES)
    )
    return counts[
        [
            "group_label",
            "carbohydrate",
            "timepoint_hr",
            "n_subjects",
            "n_wells",
            "n_analyte_observations",
        ]
    ]


def make_forest_plot(contrasts: pd.DataFrame, path: Path) -> None:
    """Plot model CIs; shade equivalence regions only when margins are locked."""
    primary = contrasts.loc[
        contrasts["primary_contrast"]
    ].copy()
    fig, axes = plt.subplots(
        1, len(PRIMARY_ANALYTES), figsize=(13.5, 4.7), sharey=True
    )
    for axis, analyte in zip(axes, PRIMARY_ANALYTES, strict=True):
        panel = primary.loc[primary["analyte"].eq(analyte)].copy()
        panel["y"] = panel["carbohydrate"].map({"RDC": 1, "SDC": 0})

        for _, row in panel.iterrows():
            if pd.notna(row["margin_lower_uM"]) and pd.notna(
                row["margin_upper_uM"]
            ):
                axis.axvspan(
                    row["margin_lower_uM"],
                    row["margin_upper_uM"],
                    color="#d9ead3",
                    alpha=0.22,
                    zorder=0,
                )
            axis.plot(
                [row["ci95_lower_uM"], row["ci95_upper_uM"]],
                [row["y"], row["y"]],
                color="#8c8c8c",
                linewidth=2,
                zorder=1,
            )
            axis.plot(
                [row["ci90_lower_uM"], row["ci90_upper_uM"]],
                [row["y"], row["y"]],
                color="#1f4e79",
                linewidth=5,
                solid_capstyle="round",
                zorder=2,
            )
            axis.scatter(
                row["estimate_uM"],
                row["y"],
                color="#b03a2e",
                s=38,
                zorder=3,
            )
            annotation_y = row["y"] - 0.14 if row["y"] > 0.5 else row["y"] + 0.14
            annotation_va = "top" if row["y"] > 0.5 else "bottom"
            axis.text(
                0.98,
                annotation_y,
                f"frontier ±{row['minimum_margin_to_contain_90_ci_uM']:.2f}",
                transform=axis.get_yaxis_transform(),
                ha="right",
                va=annotation_va,
                fontsize=8,
                color="#555555",
            )

        axis.axvline(0, color="black", linewidth=1, linestyle=":")
        axis.set_title(analyte.capitalize())
        axis.set_xlabel("Difference in 0–48 h change (µM)", fontsize=9)
        axis.set_yticks([0, 1], labels=["SDC", "RDC"])
        axis.grid(axis="x", color="#e6e6e6", linewidth=0.7)

    fig.suptitle(
        "SCFA trajectory contrasts: 90% and 95% confidence intervals",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.01,
        "Thick line: 90% CI; thin line: 95% CI. "
        "No equivalence region is shown because no external margin is locked. "
        "Frontier values are post hoc precision summaries, not margins.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    long_data = load_long_data()
    margins = load_margins()
    deltas = donor_deltas(long_data)
    margin_lookup = margins.set_index("analyte")

    contrast_rows: list[dict[str, Any]] = []
    diagnostic_rows: list[dict[str, Any]] = []
    fitted_models: dict[str, ModelFit] = {}

    for analyte_index, analyte in enumerate(ANALYTES):
        analyte_data = long_data.loc[long_data["analyte"].eq(analyte)].copy()
        model_fit = fit_trajectory_model(analyte_data)
        fitted_models[analyte] = model_fit
        diagnostic_rows.append(
            model_diagnostics(analyte, analyte_data, model_fit)
        )

        margin_row = margin_lookup.loc[analyte]
        for carb_index, carbohydrate in enumerate(CARBOHYDRATES):
            estimates = difference_in_differences(model_fit, carbohydrate)
            bootstrap = bootstrap_delta_contrast(
                deltas,
                analyte,
                carbohydrate,
                RANDOM_SEED + analyte_index * 10 + carb_index,
            )
            influence = leave_one_donor_influence(
                deltas, analyte, carbohydrate
            )
            required_margin = max(
                abs(estimates["ci90_lower_uM"]),
                abs(estimates["ci90_upper_uM"]),
            )

            lower_margin = margin_row["margin_lower_uM"]
            upper_margin = margin_row["margin_upper_uM"]
            valid_margin = (
                pd.notna(lower_margin)
                and pd.notna(upper_margin)
                and lower_margin < 0 < upper_margin
            )
            if valid_margin:
                tost = tost_result(
                    estimates["estimate_uM"],
                    estimates["standard_error_uM"],
                    float(lower_margin),
                    float(upper_margin),
                )
            else:
                tost = unavailable_tost()

            contrast_rows.append(
                {
                    "analyte": analyte,
                    "carbohydrate": carbohydrate,
                    "estimand": (
                        "(Obesity 48h - Obesity 0h) - "
                        "(Healthy-weight 48h - Healthy-weight 0h)"
                    ),
                    "primary_contrast": (
                        analyte in PRIMARY_ANALYTES
                        and carbohydrate in PRIMARY_CARBOHYDRATES
                    ),
                    "unit": "uM_pending_facility_confirmation",
                    **bootstrap,
                    **estimates,
                    "margin_lower_uM": lower_margin,
                    "margin_upper_uM": upper_margin,
                    "margin_status": margin_row["status"],
                    "margin_version": margin_row["version"],
                    "margin_source": margin_row["source"],
                    "margin_rationale": margin_row["rationale"],
                    "minimum_margin_to_contain_90_ci_uM": required_margin,
                    "minimum_margin_is_post_hoc": True,
                    **tost,
                    **influence,
                    "model_type": model_fit.model_type,
                    "optimizer": model_fit.optimizer,
                }
            )

    contrasts = pd.DataFrame(contrast_rows)
    diagnostics = pd.DataFrame(diagnostic_rows)
    counts = sample_counts(long_data)

    primary = contrasts.loc[contrasts["primary_contrast"]]
    if primary["equivalence_verdict"].eq(
        "equivalent_within_prespecified_margin"
    ).all():
        global_verdict = "all_six_primary_contrasts_equivalent"
    elif primary["equivalence_verdict"].eq(
        "not_evaluable_margin_unavailable"
    ).any():
        global_verdict = "not_evaluable_external_margins_unavailable"
    else:
        global_verdict = "global_equivalence_not_established"

    analysis_status = pd.DataFrame(
        [
            {
                "analysis_date": pd.Timestamp.now().date().isoformat(),
                "planned_participants": 40,
                "analyzed_subjects": int(long_data["subject"].nunique()),
                "healthy_weight_subjects": int(
                    long_data.loc[
                        long_data["group"].eq("control"), "subject"
                    ].nunique()
                ),
                "obesity_subjects": int(
                    long_data.loc[
                        long_data["group"].eq("case"), "subject"
                    ].nunique()
                ),
                "subject_definition": "numeric identifier plus A/B suffix",
                "well_definition": (
                    "R1/R2 and S1/S2 biological culture wells; NC well 1"
                ),
                "primary_estimands": 6,
                "global_equivalence_verdict": global_verdict,
                "margin_specification": str(MARGIN_FILE.relative_to(ROOT)),
                "formal_tost_performed": bool(
                    primary["tost_p"].notna().all()
                ),
                "unit_status": "uM label pending facility confirmation",
            }
        ]
    )

    outputs = {
        "obesity_group_scfa_contrasts.csv": contrasts,
        "obesity_group_scfa_model_diagnostics.csv": diagnostics,
        "obesity_group_scfa_sample_counts.csv": counts,
        "obesity_group_scfa_equivalence_status.csv": analysis_status,
    }
    for filename, frame in outputs.items():
        output_path = RESULTS_DIR / filename
        frame.to_csv(output_path, index=False, float_format="%.8g")
        print(f"Wrote {len(frame):>2} rows: {output_path.relative_to(ROOT)}")

    figure_path = FIGURES_DIR / "obesity_group_scfa_equivalence_forest.png"
    make_forest_plot(contrasts, figure_path)
    print(f"Wrote figure: {figure_path.relative_to(ROOT)}")

    manifest = {
        "script": str(Path(__file__).relative_to(ROOT)),
        "data": str(DATA_FILE.relative_to(ROOT)),
        "margin_file": str(MARGIN_FILE.relative_to(ROOT)),
        "outputs": [
            str((RESULTS_DIR / name).relative_to(ROOT)) for name in outputs
        ]
        + [str(figure_path.relative_to(ROOT))],
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "random_seed": RANDOM_SEED,
    }
    manifest_path = RESULTS_DIR / "obesity_group_scfa_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote manifest: {manifest_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
