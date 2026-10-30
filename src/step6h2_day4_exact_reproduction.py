from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

root = Path.cwd()
out = root / "results/step6"

day = pd.read_csv(
    root / "results/day4/LuCa7_full_AF_1310_20261030_190648.csv"
)
old = pd.read_csv(
    root / "results/day4/radius_sensitivity/radius_30/ck_cell_neighborhoods.csv"
)
mapping = pd.read_csv(
    out / "day4_reference_mapping.csv"
).sort_values("day4_row").reset_index(drop=True)
cmp = pd.read_csv(
    out / "pixel_intensity_comparison_full.csv"
)

assert len(day) == len(mapping) == 1310
assert day["cell_id"].is_unique
assert old["cell_id"].is_unique

ref_rows = mapping["reference_row"].to_numpy(int)
assert np.array_equal(
    mapping["day4_row"].to_numpy(int),
    np.arange(1310)
)

coords = day[["x_um", "y_um"]].to_numpy(float)
assert np.isfinite(coords).all()

# Reconstruct independent CK/CD8 measurements
p = cmp[
    (cmp["compartment"] == "Cell") &
    cmp["marker"].isin(["CK (Opal 690)", "CD8 (Opal 540)"])
].pivot(
    index="reference_row",
    columns="marker",
    values="python_mean"
).reindex(ref_rows)

assert not p.isna().any().any()

q_ck = day["CK"].to_numpy(float) >= 3.5
q_cd8 = day["CD8"].to_numpy(float) >= 0.5

p_ck = p["CK (Opal 690)"].to_numpy(float) >= 3.5
p_cd8 = p["CD8 (Opal 540)"].to_numpy(float) >= 0.5

# Check the original Day 4 center-cell definition
q_centers = np.flatnonzero(q_ck & ~q_cd8)
p_centers = np.flatnonzero(p_ck & ~p_cd8)

ids = day["cell_id"].astype(str).to_numpy()
old_ids = old["cell_id"].astype(str).to_numpy()

print("========== ORIGINAL DAY4 CENTER QC ==========")
print("Original file rows:", len(old))
print("CK+/CD8- expected:", len(q_centers))
print("Python CK+/CD8-:", len(p_centers))

same_set = set(old_ids) == set(ids[q_centers])
print("Exact cell ID set:", same_set)

if not same_set:
    print("Missing from original:",
          len(set(ids[q_centers]) - set(old_ids)))
    print("Extra in original:",
          len(set(old_ids) - set(ids[q_centers])))
    raise SystemExit(
        "STOP: Original Day4 center definition is different"
    )

# Align all comparisons to original Day 4 row order
id_to_index = {cid: i for i, cid in enumerate(ids)}
center_indices = np.array(
    [id_to_index[cid] for cid in old_ids], dtype=int
)

def calculate(radius, ck, cd8):
    target_ids = np.flatnonzero(~ck & cd8)

    if len(target_ids) == 0:
        raise ValueError("No CK-/CD8+ cells")

    tree = cKDTree(coords[target_ids])
    center_xy = coords[center_indices]

    nearest, _ = tree.query(center_xy, k=1)
    groups = tree.query_ball_point(center_xy, r=radius)

    counts = np.array([len(x) for x in groups])
    edges = set()

    for center_idx, group in zip(center_indices, groups):
        for j in group:
            edges.add((int(center_idx), int(target_ids[j])))

    return counts, nearest, edges

q_count, q_nearest, q_edges = calculate(30, q_ck, q_cd8)
p_count, p_nearest, p_edges = calculate(30, p_ck, p_cd8)

old_count = pd.to_numeric(
    old["n_CD8_only_within_radius"],
    errors="coerce"
).to_numpy(float)

old_has = (
    old["has_CD8_only_within_radius"]
    .astype(str).str.strip().str.lower()
    .isin(["true", "1", "1.0", "yes"])
    .to_numpy()
)

old_nearest = pd.to_numeric(
    old["nearest_CD8_only_um"],
    errors="coerce"
).to_numpy(float)

count_match = np.isfinite(old_count) & (old_count == q_count)
has_match = old_has == (q_count > 0)

# Compare two possible original nearest-distance conventions:
# nearest anywhere, or nearest only if within the radius.
nearest_global_match = np.isclose(
    old_nearest, q_nearest,
    atol=0.01, rtol=0, equal_nan=True
)

nearest_within = np.where(q_nearest <= 30, q_nearest, np.nan)
nearest_within_match = np.isclose(
    old_nearest, nearest_within,
    atol=0.01, rtol=0, equal_nan=True
)

print("\n========== DAY4 ORIGINAL REPRODUCTION ==========")
print("Neighbor count matches:", int(count_match.sum()), "/", len(old))
print("Neighbor status matches:", int(has_match.sum()), "/", len(old))
print("Nearest distance matches (global):",
      int(nearest_global_match.sum()), "/", len(old))
print("Nearest distance matches (within radius):",
      int(nearest_within_match.sum()), "/", len(old))

if not count_match.all():
    mismatch = old.loc[~count_match, ["cell_id"]].copy()
    mismatch["old_count"] = old_count[~count_match]
    mismatch["recomputed_count"] = q_count[~count_match]
    mismatch.to_csv(
        out / "day4_neighbor_count_mismatches.csv",
        index=False
    )

original_pass = (
    count_match.all()
    and has_match.all()
    and (
        nearest_global_match.all()
        or nearest_within_match.all()
    )
)

print("DAY4_ORIGINAL_REPRODUCTION:",
      "PASS" if original_pass else "REVIEW_REQUIRED")

# Compare QuPath-based vs independently remeasured phenotypes
rows = []

for radius in [15, 20, 30, 50, 75]:
    qc, _, qe = calculate(radius, q_ck, q_cd8)
    pc, _, pe = calculate(radius, p_ck, p_cd8)

    rows.append({
        "radius_um": radius,
        "centers": len(center_indices),
        "qupath_edges": len(qe),
        "python_edges": len(pe),
        "lost_edges": len(qe - pe),
        "gained_edges": len(pe - qe),
        "qupath_covered": int((qc > 0).sum()),
        "python_covered": int((pc > 0).sum()),
        "coverage_flips": int(
            np.sum((qc > 0) != (pc > 0))
        ),
        "cells_with_neighbor_count_change": int(
            np.sum(qc != pc)
        )
    })

summary = pd.DataFrame(rows)
summary.to_csv(
    out / "day4_exact_spatial_robustness.csv",
    index=False
)

cell_audit = pd.DataFrame({
    "cell_id": old_ids,
    "reference_row": ref_rows[center_indices],
    "original_neighbor_count": old_count,
    "qupath_recomputed_count": q_count,
    "python_recomputed_count": p_count,
    "qupath_has_neighbor": q_count > 0,
    "python_has_neighbor": p_count > 0
})
cell_audit.to_csv(
    out / "day4_exact_30um_cell_audit.csv",
    index=False
)

print("\n========== EXACT 439-CELL SPATIAL ROBUSTNESS ==========")
print(summary.to_string(index=False))

print("\n========== PRIMARY 30 UM ==========")
r = summary[summary["radius_um"] == 30].iloc[0]
print("CK+/CD8- centers:", int(r["centers"]))
print("Edges:", int(r["qupath_edges"]),
      "->", int(r["python_edges"]))
print("Covered cells:", int(r["qupath_covered"]),
      "->", int(r["python_covered"]))
print("Coverage status flips:", int(r["coverage_flips"]))
print("Neighbor count changed cells:",
      int(r["cells_with_neighbor_count_change"]))

print("\nSTAGE_6H2: COMPLETE")
