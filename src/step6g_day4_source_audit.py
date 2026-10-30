from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

root = Path.cwd()
out = root / "results/step6"
out.mkdir(parents=True, exist_ok=True)

ref = pd.read_csv(
    root / "LuCa7_roi01_cells.csv.txt",
    sep=None, engine="python"
)

day4 = pd.read_csv(
    root / "results/day4/LuCa7_full_AF_1310_20261030_190648.csv"
)

ref_xy = ref[
    ["Centroid X µm", "Centroid Y µm"]
].to_numpy(float)

day_xy = day4[
    ["x_um", "y_um"]
].to_numpy(float)

tree = cKDTree(ref_xy)
dist, idx = tree.query(day_xy, k=1)

matched = (
    len(ref) == len(day4) == 1310
    and np.isfinite(dist).all()
    and (dist <= 0.5).all()
    and len(np.unique(idx)) == 1310
)

print("========== DAY4 CELL ALIGNMENT ==========")
print("Reference cells:", len(ref))
print("Day4 cells:", len(day4))
print("Unique matches:", len(np.unique(idx)))
print("Maximum distance um:", float(np.max(dist)))
print("Alignment:", "PASS" if matched else "REVIEW_REQUIRED")

if not matched:
    raise SystemExit("Stop: Day4 alignment not validated")

mapping = pd.DataFrame({
    "day4_row": np.arange(len(day4)),
    "reference_row": idx,
    "distance_um": dist,
    "cell_id": day4["cell_id"].astype(str)
})

mapping.to_csv(
    out / "day4_reference_mapping.csv",
    index=False
)

markers = {
    "PDL1": "PDL1 (Opal 520)",
    "CD8": "CD8 (Opal 540)",
    "FoxP3": "FoxP3 (Opal 570)",
    "CD68": "CD68 (Opal 620)",
    "PD1": "PD1 (Opal 650)",
    "CK": "CK (Opal 690)",
    "DAPI": "DAPI",
    "Autofluorescence": "Autofluorescence"
}

rows = []

for short, full in markers.items():
    observed = pd.to_numeric(
        day4[short], errors="coerce"
    ).to_numpy(float)

    for comp in ["Cell", "Nucleus", "Cytoplasm"]:
        original = pd.to_numeric(
            ref.iloc[idx][f"{comp}: {full} mean"],
            errors="coerce"
        ).to_numpy(float)

        valid = (
            np.isfinite(observed) &
            np.isfinite(original)
        )

        if not valid.any():
            continue

        difference = observed[valid] - original[valid]

        rows.append({
            "marker": short,
            "compartment": comp,
            "n_valid": int(valid.sum()),
            "MAE": float(np.mean(np.abs(difference))),
            "max_abs_error": float(np.max(np.abs(difference))),
            "bias": float(np.mean(difference)),
            "fraction_within_0_001": float(
                np.mean(np.abs(difference) <= 0.001)
            )
        })

audit = pd.DataFrame(rows)
audit.to_csv(
    out / "day4_compartment_source_audit.csv",
    index=False
)

print("\n========== DAY4 COMPARTMENT COMPARISON ==========")
print(audit.round(6).to_string(index=False))

print("\n========== CLOSEST REFERENCE BY MARKER ==========")
best = (
    audit.sort_values("MAE")
    .drop_duplicates("marker")
)
print(best.round(6).to_string(index=False))

print("\n========== DAY4 CK/CD8 EXPLORATORY GATE ==========")
ck = pd.to_numeric(day4["CK"], errors="coerce")
cd8 = pd.to_numeric(day4["CD8"], errors="coerce")
double_high = (ck >= 3.5) & (cd8 >= 0.5)

print("CK+/CD8+ count in Day4:", int(double_high.sum()))

print("\n========== THREE DISCORDANT CELLS ==========")
discrepant = [54, 182, 1109]

for ref_row in discrepant:
    hits = mapping[
        mapping["reference_row"] == ref_row
    ]

    if len(hits) != 1:
        print("Unresolved reference row:", ref_row)
        continue

    k = int(hits.iloc[0]["day4_row"])
    r = day4.iloc[k]

    print({
        "reference_row": ref_row,
        "day4_cell_id": str(r["cell_id"]),
        "x_um": float(r["x_um"]),
        "y_um": float(r["y_um"]),
        "day4_CK": float(r["CK"]),
        "day4_CD8": float(r["CD8"])
    })

print("\nSTAGE_6G_SOURCE_AUDIT: COMPLETE")
print("BIOLOGICAL_PHENOTYPING_VALIDATION: PENDING")
