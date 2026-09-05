"""
Diagnostic experiment for T-zone enrichment mechanisms.


Experiment : activation-localization experiment
    - activation regimes:
        non_local: beta_a_R=1.0, beta_a_T=1.2
        medium:    beta_a_R=2.5, beta_a_T=3.0
        local:     beta_a_R=5.0, beta_a_T=3.0
      with lambda_a_T=0.5 and lambda_a_c=0.005 fixed in all regimes.
    - in each activation regime, sweep lambda_d and the three ERK regimes.
    - each simulation keeps the original burn-in 500 deaths and analyzes the next 1000 deaths.

Statistical analysis:
    - pool the independent replicates within each pre-specified parameter setting;
    - perform the one-sided exact binomial test
          H0: pi_T = p0 = |T|/|W|
          H1: pi_T > p0;
    - report the raw binomial p-value and p < 0.05 indicator;

CSV binomial-test outputs are written to:
    diagnostic_outputs_Tzone_mechanisms/
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List
import importlib.util
import sys

import numpy as np
import pandas as pd

try:
    from scipy.stats import binomtest
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This diagnostic script requires scipy: binomtest.") from exc


# ============================================================
# 0. User-adjustable settings
# ============================================================

# For a quick smoke test use 3 or 5. For final results use 20 if runtime allows.
N_REPLICATES = 10
SEED0 = 7001

BURN_IN_DEATHS = 500
TARGET_DEATHS = 1000
MAX_PROPOSALS = 120_000_000

LAMBDA_D_VALUES = [0.25, 0.5, 1.0, 2.0, 4.0]

ERK_REGIMES = [
    {
        "erk_regime": "strong_ERK_protection",
        "beta_d_R": 1.0,
        "beta_d_T": 0.4,
        "meaning": "large ERK radius, long ERK lifetime",
    },
    {
        "erk_regime": "reference_ERK",
        "beta_d_R": 2.0,
        "beta_d_T": 0.8,
        "meaning": "reference ERK protection",
    },
    {
        "erk_regime": "weak_ERK_protection",
        "beta_d_R": 4.0,
        "beta_d_T": 1.6,
        "meaning": "small ERK radius, short ERK lifetime",
    },
]

ACTIVATION_REGIMES = [
    {
        "activation_regime": "non_local_activation",
        "beta_a_R": 1.0,
        "beta_a_T": 1.2,
        "meaning": "large and long-lived active disks",
    },
    {
        "activation_regime": "medium_activation",
        "beta_a_R": 2.5,
        "beta_a_T": 3.0,
        "meaning": "current visible reference regime",
    },
    {
        "activation_regime": "local_activation",
        "beta_a_R": 5.0,
        "beta_a_T": 3.0,
        "meaning": "small and short-lived active disks",
    },
]

# Fixed T-zone activation contrast for these diagnostic experiments.
FIXED_LAMBDA_A_T = 0.5
FIXED_LAMBDA_A_C = 0.005

OUTPUT_DIR = Path("diagnostic_outputs_Tzone_mechanisms")

# The script will try to import the first existing file from this list.
SIMULATION_FILE_CANDIDATES = [
    "step1code.py",
    "step2code.py",
]


# ============================================================
# 1. Load the user's simulation module
# ============================================================


def load_simulation_module() -> Any:
    here = Path(__file__).resolve().parent
    candidates = []
    for name in SIMULATION_FILE_CANDIDATES:
        candidates.append(here / name)
    # Also search the current working directory, in case the script is run elsewhere.
    cwd = Path.cwd().resolve()
    if cwd != here:
        for name in SIMULATION_FILE_CANDIDATES:
            candidates.append(cwd / name)

    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location("user_simulation_module", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules["user_simulation_module"] = module
            spec.loader.exec_module(module)
            print(f"Loaded simulation module: {path}")
            return module

    tried = "\n".join(str(p) for p in candidates)
    raise FileNotFoundError(
        "Could not find a simulation file. Put this diagnostic script next to one of:\n"
        + "\n".join(SIMULATION_FILE_CANDIDATES)
        + "\n\nTried:\n"
        + tried
    )


# ============================================================
# 2. Shared utilities
# ============================================================


def pstr(x: float) -> str:
    return str(x).replace(".", "p").replace("-", "m")


def seed_for(rep_idx: int, offset: int = 0) -> int:
    return SEED0 + offset + rep_idx


def make_base_parameters(sim: Any) -> Any:
    return sim.Parameters(
        burn_in_deaths=BURN_IN_DEATHS,
        target_deaths=TARGET_DEATHS,
        max_proposals=MAX_PROPOSALS,
    )


def make_parameters(
    sim: Any,
    base: Any,
    *,
    seed: int,
    lambda_d: float,
    beta_d_R: float,
    beta_d_T: float,
    beta_a_R: float = 2.5,
    beta_a_T: float = 3.0,
    lambda_a_T: float = FIXED_LAMBDA_A_T,
    lambda_a_c: float = FIXED_LAMBDA_A_C,
) -> Any:
    return replace(
        base,
        seed=seed,
        lambda_a_T=lambda_a_T,
        lambda_a_c=lambda_a_c,
        beta_a_R=beta_a_R,
        beta_a_T=beta_a_T,
        lambda_d=lambda_d,
        beta_d_R=beta_d_R,
        beta_d_T=beta_d_T,
        max_proposals=MAX_PROPOSALS,
    )


def t_zone_stats_from_arrays(sim: Any, x: np.ndarray, y: np.ndarray, p: Any) -> Dict[str, float | int]:
    n = int(len(x))
    if n == 0:
        area_T = float(sim.t_zone_area(p))
        area_W = float(p.area_W)
        p0 = area_T / area_W if area_W > 0 else np.nan
        return {
            "n_observed_deaths_analyzed": 0,
            "t_zone_binomial_test_k": 0,
            "t_zone_binomial_test_n": 0,
            "t_zone_area_fraction": p0,
            "t_zone_observed_fraction": np.nan,
            "t_zone_density_ratio": np.nan,
            "t_zone_binomial_p_value_greater": np.nan,
        }

    in_T = np.array([sim.is_inside_T_zone(float(xx), float(yy), p) for xx, yy in zip(x, y)], dtype=bool)
    k = int(np.sum(in_T))
    n_out = n - k

    area_T = float(sim.t_zone_area(p))
    area_W = float(p.area_W)
    area_out = area_W - area_T
    p0 = area_T / area_W if area_W > 0 else np.nan

    frac = k / n
    density_T = k / area_T if area_T > 0 else np.nan
    density_out = n_out / area_out if area_out > 0 else np.nan
    ratio = density_T / density_out if density_out > 0 else np.nan
    pval = binomtest(k=k, n=n, p=p0, alternative="greater").pvalue if n > 0 else np.nan

    return {
        "n_observed_deaths_analyzed": n,
        "t_zone_binomial_test_k": k,
        "t_zone_binomial_test_n": n,
        "t_zone_area_fraction": p0,
        "t_zone_observed_fraction": frac,
        "t_zone_density_ratio": ratio,
        "t_zone_binomial_p_value_greater": pval,
    }


def activations_during_window(result: Any, t0: float, t1: float) -> int:
    act_t = np.asarray(result.activation_t, dtype=float)
    if act_t.size == 0:
        return 0
    return int(np.sum((act_t > t0) & (act_t <= t1)))


def row_from_fixed_death_result(
    sim: Any,
    result: Any,
    *,
    experiment: str,
    base_scenario: str,
    replicate: int,
    activation_regime: str,
    erk_regime: str,
    meaning: str = "",
) -> Dict[str, Any]:
    p = result.params
    stats = t_zone_stats_from_arrays(sim, np.asarray(result.death_x), np.asarray(result.death_y), p)
    analysis_start = float(result.burn_in_time)
    analysis_end = float(result.raw_final_time if hasattr(result, "raw_final_time") else result.final_time)
    # result has final_time; raw_final_time is only a CSV column.
    analysis_end = float(result.final_time)
    return {
        "experiment": experiment,
        "base_scenario": base_scenario,
        "replicate": replicate,
        "seed": p.seed,
        "activation_regime": activation_regime,
        "erk_regime": erk_regime,
        "regime_meaning": meaning,
        "lambda_a_T": p.lambda_a_T,
        "lambda_a_c": p.lambda_a_c,
        "beta_a_R": p.beta_a_R,
        "beta_a_T": p.beta_a_T,
        "lambda_d": p.lambda_d,
        "beta_d_R": p.beta_d_R,
        "beta_d_T": p.beta_d_T,
        "burn_in_deaths": p.burn_in_deaths,
        "target_deaths": p.target_deaths,
        "burn_in_time": float(result.burn_in_time),
        "analysis_time_span": float(result.analysis_time_span),
        "observed_death_rate_after_burnin": (
            float(result.n_observed_deaths_analyzed) / float(result.analysis_time_span)
            if float(result.analysis_time_span) > 0
            else np.nan
        ),
        "accepted_activations_total": int(result.n_accepted_activations),
        "accepted_activations_during_analysis": activations_during_window(result, analysis_start, analysis_end),
        "death_candidates": int(result.n_death_candidates),
        "death_rejected": int(result.n_death_rejected),
        "activation_candidates": int(result.n_activation_candidates),
        "activation_rejected": int(result.n_activation_rejected),
        "active_center_expirations": int(result.n_active_expirations),
        "erk_zone_expirations": int(result.n_erk_expirations),
        "all_proposals": int(result.n_proposals),
        "stopped_by_target": bool(result.stopped_by_target),
        **stats,
    }


# ============================================================
# 3. Pooled exact binomial test 
# ============================================================


def pooled_binomial(df: pd.DataFrame, group_cols: List[str], alpha: float = 0.05) -> pd.DataFrame:
    """Pool replicates by parameter setting and run a raw one-sided exact binomial test."""
    rows: List[Dict[str, Any]] = []

    for keys, g in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)

        row = {col: val for col, val in zip(group_cols, keys)}
        pooled_k = int(g["t_zone_binomial_test_k"].sum())
        pooled_n = int(g["t_zone_binomial_test_n"].sum())
        p0 = float(g["t_zone_area_fraction"].iloc[0])

        pooled_fraction = pooled_k / pooled_n if pooled_n > 0 else np.nan
        n_out = pooled_n - pooled_k
        density_ratio = (
            (pooled_k / p0) / (n_out / (1.0 - p0))
            if pooled_n > 0 and n_out > 0 and 0.0 < p0 < 1.0
            else np.nan
        )
        pval = (
            binomtest(k=pooled_k, n=pooled_n, p=p0, alternative="greater").pvalue
            if pooled_n > 0
            else np.nan
        )

        row.update(
            {
                "n_replicates": int(len(g)),
                "pooled_t_zone_deaths_k": pooled_k,
                "pooled_total_deaths_n": pooled_n,
                "t_zone_area_fraction_p0": p0,
                "pooled_t_zone_fraction": pooled_fraction,
                "pooled_t_zone_density_ratio": density_ratio,
                "pooled_binomial_p_value_greater": pval,
                "significant_nominal_0p05": bool(pval < alpha) if np.isfinite(pval) else False,
            }
        )
        rows.append(row)

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(
            ["significant_nominal_0p05", "pooled_binomial_p_value_greater"],
            ascending=[False, True],
        ).reset_index(drop=True)
    return out


# ============================================================
# 4. Experiment 2: Expérience de localisation de l'activation
# ============================================================


def run_experiment_2(sim: Any, out_dir: Path) -> pd.DataFrame:
    """Compare non-local, medium and local activation under the same death/ERK grid."""
    base = make_base_parameters(sim)
    rows: List[Dict[str, Any]] = []

    print("\n=== Activation localization x death/ERK sweep ===")
    for ar_idx, ar in enumerate(ACTIVATION_REGIMES):
        for er_idx, er in enumerate(ERK_REGIMES):
            for ld_idx, lambda_d in enumerate(LAMBDA_D_VALUES):
                base_scenario = (
                    f"{ar['activation_regime']}"
                    f"_{er['erk_regime']}"
                    f"_ld{pstr(lambda_d)}"
                )
                for rep in range(1, N_REPLICATES + 1):
                    p = make_parameters(
                        sim,
                        base,
                        seed=seed_for(
                            rep,
                            offset=(
                                30_000
                                + 5000 * ar_idx
                                + 1000 * er_idx
                                + 100 * ld_idx
                            ),
                        ),
                        lambda_d=lambda_d,
                        beta_d_R=er["beta_d_R"],
                        beta_d_T=er["beta_d_T"],
                        beta_a_R=ar["beta_a_R"],
                        beta_a_T=ar["beta_a_T"],
                    )
                    scenario = f"{base_scenario}__rep{rep:02d}"
                    print(scenario)
                    result = sim.simulate_one(scenario, p)
                    rows.append(
                        row_from_fixed_death_result(
                            sim,
                            result,
                            experiment="activation_localization_x_deathERK",
                            base_scenario=base_scenario,
                            replicate=rep,
                            activation_regime=ar["activation_regime"],
                            erk_regime=er["erk_regime"],
                            meaning=f"activation: {ar['meaning']}; ERK: {er['meaning']}",
                        )
                    )

    df = pd.DataFrame(rows)
    group_cols = [
        "activation_regime",
        "beta_a_R",
        "beta_a_T",
        "erk_regime",
        "beta_d_R",
        "beta_d_T",
        "lambda_d",
    ]
    pooled = pooled_binomial(df, group_cols)
    pooled.to_csv(
        out_dir / "activation_localization_deathERK_pooled_binomial.csv",
        index=False,
    )

    # Compact table for the Section 2.5.3 conclusion (0/15, 15/15, ...).
    robust = pooled.groupby(
        ["activation_regime", "beta_a_R", "beta_a_T"],
        dropna=False,
    ).agg(
        n_settings=("pooled_t_zone_fraction", "size"),
        n_significant_p005=("significant_nominal_0p05", "sum"),
        min_fraction=("pooled_t_zone_fraction", "min"),
        mean_fraction=("pooled_t_zone_fraction", "mean"),
        max_fraction=("pooled_t_zone_fraction", "max"),
        min_density_ratio=("pooled_t_zone_density_ratio", "min"),
        mean_density_ratio=("pooled_t_zone_density_ratio", "mean"),
        max_density_ratio=("pooled_t_zone_density_ratio", "max"),
        min_p_value=("pooled_binomial_p_value_greater", "min"),
        max_p_value=("pooled_binomial_p_value_greater", "max"),
    ).reset_index()
    robust.to_csv(
        out_dir / "activation_regime_robustness_binomial.csv",
        index=False,
    )

    return pooled


# ============================================================
# 5. Main
# ============================================================


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sim = load_simulation_module()

    pooled_2 = run_experiment_2(sim, OUTPUT_DIR)

    print("\nDone.")
    print(f"Outputs are in: {OUTPUT_DIR.resolve()}")
    print("Raw one-sided exact binomial tests")
    print("Output files:")
    print("  1. experiment_activation_localization_deathERK_pooled_binomial.csv")
    print("  2. experiment_activation_regime_robustness_binomial.csv")
    print(f"Experiment  parameter settings: {len(pooled_2)}")


if __name__ == "__main__":
    main()
