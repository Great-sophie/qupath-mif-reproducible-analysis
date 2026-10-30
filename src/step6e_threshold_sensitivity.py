from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("results/step6")
df = pd.read_csv(
    OUT / "pixel_intensity_comparison_full.csv"
)

assert len(df) == 31440
assert not df[
    ["qupath_mean", "python_mean"]
].isna().any().any()

targets = [
    ("Cell", "CK (Opal 690)"),
    ("Cell", "CD8 (Opal 540)"),
    ("Cell", "CD68 (Opal 620)"),
    ("Nucleus", "FoxP3 (Opal 570)"),
    ("Cell", "PD1 (Opal 650)"),
    ("Cell", "PDL1 (Opal 520)")
]

rows = []
changes = []

for comp, marker in targets:
    g = df[
        (df["compartment"] == comp) &
        (df["marker"] == marker)
    ].copy()

    assert len(g) == 1310
    assert g["reference_row"].is_unique

    a = g["qupath_mean"].to_numpy(float)
    b = g["python_mean"].to_numpy(float)

    for q in [0.50, 0.75, 0.90]:
        threshold = float(np.quantile(a, q))

        ref_positive = a >= threshold
        py_positive = b >= threshold

        flipped = ref_positive != py_positive

        rows.append({
            "compartment": comp,
            "marker": marker,
            "reference_quantile": q,
            "threshold": threshold,
            "reference_positive": int(ref_positive.sum()),
            "python_positive": int(py_positive.sum()),
            "flipped_cells": int(flipped.sum()),
            "flip_rate_pct": 100 * float(flipped.mean()),
            "positive_to_negative": int(
                (ref_positive & ~py_positive).sum()
            ),
            "negative_to_positive": int(
                (~ref_positive & py_positive).sum()
            )
        })

        for pos in np.where(flipped)[0]:
            changes.append({
                "reference_row": int(
                    g.iloc[pos]["reference_row"]
                ),
                "marker": marker,
                "compartment": comp,
                "reference_quantile": q,
                "threshold": threshold,
                "qupath_mean": a[pos],
                "python_mean": b[pos],
                "direction": (
                    "positive_to_negative"
                    if ref_positive[pos]
                    else "negative_to_positive"
                )
            })

summary = pd.DataFrame(rows)
summary.to_csv(
    OUT / "phenotyping_threshold_sensitivity.csv",
    index=False
)

pd.DataFrame(
    changes,
    columns=[
        "reference_row", "marker", "compartment",
        "reference_quantile", "threshold",
        "qupath_mean", "python_mean", "direction"
    ]
).to_csv(
    OUT / "threshold_discordant_cells.csv",
    index=False
)

print("\n========== STAGE 6E THRESHOLD SENSITIVITY ==========")
print(
    summary[
        ["marker", "compartment", "reference_quantile",
         "flipped_cells", "flip_rate_pct",
         "positive_to_negative", "negative_to_positive"]
    ].round(3).to_string(index=False)
)

print("\n========== SUMMARY ==========")
print("Marker-compartment targets:", len(targets))
print("Diagnostic thresholds:", len(summary))
print("Worst flip rate:",
      round(summary["flip_rate_pct"].max(), 3), "%")
print("Median flip rate:",
      round(summary["flip_rate_pct"].median(), 3), "%")
print("Threshold tests with zero flips:",
      int((summary["flipped_cells"] == 0).sum()),
      "/", len(summary))

print("\nSTAGE_6E_THRESHOLD_SENSITIVITY: COMPLETE")
print("BIOLOGICAL_PHENOTYPING_VALIDATION: PENDING")
