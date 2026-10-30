from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT = Path("results/step6")

df = pd.read_csv(
    OUT / "pixel_intensity_comparison_full.csv"
)

required = {
    "roi_row", "reference_row", "compartment",
    "marker", "qupath_mean", "python_mean"
}
assert required.issubset(df.columns)
assert len(df) == 31440
assert df["roi_row"].nunique() == 1310
assert not df[["qupath_mean", "python_mean"]].isna().any().any()

rows = []

for (comp, marker), g in df.groupby(
    ["compartment", "marker"], sort=False
):
    a = g["qupath_mean"].to_numpy(float)
    b = g["python_mean"].to_numpy(float)
    err = b - a

    p05, p95 = np.percentile(a, [5, 95])
    dynamic_range = p95 - p05

    rho = spearmanr(a, b).statistic

    rows.append({
        "compartment": comp,
        "marker": marker,
        "n": len(g),
        "reference_p05": p05,
        "reference_p95": p95,
        "MAE": np.mean(np.abs(err)),
        "P95_abs_error": np.percentile(np.abs(err), 95),
        "max_abs_error": np.max(np.abs(err)),
        "bias": np.mean(err),
        "normalized_MAE_pct": (
            100 * np.mean(np.abs(err)) / dynamic_range
            if dynamic_range > 0 else np.nan
        ),
        "spearman_rho": rho
    })

summary = pd.DataFrame(rows)
summary = summary.sort_values(
    "normalized_MAE_pct", ascending=False
)

summary.to_csv(
    OUT / "marker_robustness_full.csv",
    index=False
)

# Identify the most discrepant cell-marker combinations
df = df.copy()
df["abs_error"] = (
    df["python_mean"] - df["qupath_mean"]
).abs()

outliers = df.nlargest(30, "abs_error")
outliers.to_csv(
    OUT / "top30_intensity_discrepancies.csv",
    index=False
)

print("\n========== MARKER ROBUSTNESS ==========")
print(summary.round(5).to_string(index=False))

print("\n========== OVERALL QC ==========")
print("Cells:", df["roi_row"].nunique())
print("Measurements:", len(df))
print("Marker-compartment combinations:", len(summary))
print("Global MAE:", round(df["abs_error"].mean(), 6))
print("Global P95:", round(df["abs_error"].quantile(.95), 6))

print("\n========== RANKING STABILITY ==========")
print("Minimum Spearman rho:",
      round(summary["spearman_rho"].min(), 6))
print("Median Spearman rho:",
      round(summary["spearman_rho"].median(), 6))
print("Combinations with rho >= 0.99:",
      int((summary["spearman_rho"] >= 0.99).sum()),
      "/", len(summary))

print("\nSTAGE_6D_MARKER_ROBUSTNESS: COMPLETE")
