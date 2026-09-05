import pandas as pd
import numpy as np
from scipy.stats import binomtest

# ============================
# Input / output
# ============================

input_csv = (
    "step2result.csv"
)

output_csv = (
    "step2pooled_binomial_summary.csv"
)

# ============================
# Read replicate-level results
# ============================

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

# ============================
# Pool the replicates
# ============================

rows = []

for base_scenario, g in df.groupby("base_scenario", sort=False):

    pooled_k = int(g["t_zone_binomial_test_k"].sum())
    pooled_n = int(g["t_zone_binomial_test_n"].sum())

    p0 = float(
        g["t_zone_binomial_null_probability_area_fraction"].iloc[0]
    )

    pooled_fraction = pooled_k / pooled_n

    # One-sided exact binomial test:
    # H0: pi = p0
    # HA: pi > p0
    pooled_p_value = binomtest(
        k=pooled_k,
        n=pooled_n,
        p=p0,
        alternative="greater",
    ).pvalue

    n_out = pooled_n - pooled_k

    if n_out > 0 and 0 < p0 < 1:
        pooled_density_ratio = (
            (pooled_k / p0)
            /
            (n_out / (1.0 - p0))
        )
    else:
        pooled_density_ratio = np.nan

    rows.append({
        "base_scenario": base_scenario,
        "n_replicates": len(g),
        "t_zone_deaths_k": pooled_k,
        "total_deaths_n": pooled_n,
        "t_zone_area_fraction_p0": p0,
        "t_zone_observed_fraction": pooled_fraction,
        "t_zone_density_ratio": pooled_density_ratio,
        "binomial_p_value_greater": pooled_p_value,
        "significant_0p05": pooled_p_value < 0.05,
    })

# ============================
# Save CSV only
# ============================

result = pd.DataFrame(rows)

result = result.sort_values(
    "binomial_p_value_greater"
)

result.to_csv(output_csv, index=False)

print(f"Saved: {output_csv}")
print(f"Number of parameter settings: {len(result)}")