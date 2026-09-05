import pandas as pd
import numpy as np
from scipy.stats import binomtest

# ============
# Input file
# ============
input_csv = "step1resultsummary.csv"

# ============
# Output file
# ============
output_csv = "step1pooled_binomial_summary.csv"

alpha = 0.05

df = pd.read_csv(input_csv)

required_cols = [
    "base_scenario",
    "t_zone_binomial_test_k",
    "t_zone_binomial_test_n",
    "t_zone_binomial_null_probability_area_fraction",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in CSV: {missing}")

rows = []

for base_scenario, g in df.groupby("base_scenario", sort=False):
    pooled_k = int(g["t_zone_binomial_test_k"].sum())
    pooled_n = int(g["t_zone_binomial_test_n"].sum())
    p0 = float(g["t_zone_binomial_null_probability_area_fraction"].iloc[0])

    pooled_fraction = pooled_k / pooled_n

    pooled_p_value = binomtest(
        k=pooled_k,
        n=pooled_n,
        p=p0,
        alternative="greater",
    ).pvalue

    # pooled density ratio:
    # density in T-zone / density outside T-zone
    n_out = pooled_n - pooled_k
    if n_out > 0 and 0 < p0 < 1:
        pooled_density_ratio = (pooled_k / p0) / (n_out / (1 - p0))
    else:
        pooled_density_ratio = np.nan

    rows.append({
        "base_scenario": base_scenario,
        "n_replicates": int(len(g)),
        "pooled_t_zone_deaths_k": pooled_k,
        "pooled_total_deaths_n": pooled_n,
        "t_zone_area_fraction_p0": p0,
        "pooled_t_zone_fraction": pooled_fraction,
        "pooled_t_zone_density_ratio": pooled_density_ratio,
        "pooled_binomial_p_value_greater": pooled_p_value,
    })

pooled = pd.DataFrame(rows)

# Significant at the ordinary 5% level (no Bonferroni correction)
pooled["significant_0p05"] = (
    pooled["pooled_binomial_p_value_greater"] < alpha
)

# Sort: significant results first, then smallest p-value
pooled = pooled.sort_values(
    by=["significant_0p05", "pooled_binomial_p_value_greater"],
    ascending=[False, True],
)

pooled.to_csv(output_csv, index=False)

print(f"Saved: {output_csv}")
print(f"Number of parameter settings: {len(pooled)}")
print(
    "Significant at alpha = 0.05:",
    int(pooled["significant_0p05"].sum()),
)
