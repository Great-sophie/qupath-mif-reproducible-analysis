from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path.cwd()
OUT = ROOT / "results/step6"
OUT.mkdir(parents=True, exist_ok=True)

day4 = pd.read_csv(
    ROOT / "results/day4/LuCa7_full_AF_1310_20261030_190648.csv"
)
mapping = pd.read_csv(
    OUT / "day4_reference_mapping.csv"
).sort_values("day4_row").reset_index(drop=True)

cmp = pd.read_csv(
    OUT / "pixel_intensity_comparison_full.csv"
)

assert len(day4) == len(mapping) == 1310
assert np.array_equal(
    mapping["day4_row"].to_numpy(int),
    np.arange(1310)
)
assert mapping["reference_row"].is_unique

ref_rows = mapping["reference_row"].to_numpy(int)
coords = day4[["x_um", "y_um"]].to_numpy(float)
assert np.isfinite(coords).all()

markers = ["CK (Opal 690)", "CD8 (Opal 540)"]

g = cmp[
    (cmp["compartment"] == "Cell") &
    cmp["marker"].isin(markers)
]

qtable = g.pivot(
    index="reference_row",
    columns="marker",
    values="qupath_mean"
).reindex(ref_rows)

ptable = g.pivot(
    index="reference_row",
    columns="marker",
    values="python_mean"
).reindex(ref_rows)

assert not qtable.isna().any().any()
assert not ptable.isna().any().any()

# Verify that Day 4 corresponds to the original QuPath values
for short, full in [
    ("CK", "CK (Opal 690)"),
    ("CD8", "CD8 (Opal 540)")
]:
    delta = np.abs(
        day4[short].to_numpy(float) -
        qtable[full].to_numpy(float)
    )
    assert np.max(delta) <= 0.0001, (
        f"Day4 source mismatch: {short}"
    )

# Frozen exploratory thresholds
CK_THRESHOLD = 3.5
CD8_THRESHOLD = 0.5

q_ck = day4["CK"].to_numpy(float) >= CK_THRESHOLD
q_cd8 = day4["CD8"].to_numpy(float) >= CD8_THRESHOLD

p_ck = (
    ptable["CK (Opal 690)"].to_numpy(float)
    >= CK_THRESHOLD
)
p_cd8 = (
    ptable["CD8 (Opal 540)"].to_numpy(float)
    >= CD8_THRESHOLD
)

def phenotype(ck, cd8):
    return np.select(
        [ck & cd8, ck & ~cd8, ~ck & cd8],
        ["CK+/CD8+", "CK+/CD8-", "CK-/CD8+"],
        default="CK-/CD8-"
    )

qa = phenotype(q_ck, q_cd8)
pa = phenotype(p_ck, p_cd8)

changed = qa != pa

changes = pd.DataFrame({
    "reference_row": ref_rows,
    "cell_id": day4["cell_id"].astype(str),
    "x_um": coords[:, 0],
    "y_um": coords[:, 1],
    "qupath_phenotype": qa,
    "python_phenotype": pa,
    "changed": changed
})

changes[changes["changed"]].to_csv(
    OUT / "spatial_phenotype_discordances.csv",
    index=False
)

print("========== PHENOTYPE QC ==========")
print("Changed cells:", int(changed.sum()))
print("Changed reference rows:",
      sorted(ref_rows[changed].tolist()))

def spatial_edges(radius, ck, cd8):
    # Define immune neighbors as CK-negative, CD8-positive
    ck_ids = np.flatnonzero(ck)
    immune_ids = np.flatnonzero(~ck & cd8)

    if len(ck_ids) == 0 or len(immune_ids) == 0:
        return set()

    tree = cKDTree(coords[immune_ids])
    neighbor_lists = tree.query_ball_point(
        coords[ck_ids], r=radius
    )

    edges = set()
    for ck_i, neighbors in zip(ck_ids, neighbor_lists):
        for j in neighbors:
            edges.add((int(ck_i), int(immune_ids[j])))

    return edges

rows = []
changed_edges = []

for radius in [15, 20, 30, 50, 75]:
    qe = spatial_edges(radius, q_ck, q_cd8)
    pe = spatial_edges(radius, p_ck, p_cd8)

    lost = qe - pe
    gained = pe - qe

    q_covered = {a for a, b in qe}
    p_covered = {a for a, b in pe}

    rows.append({
        "radius_um": radius,
        "qupath_CK_positive": int(q_ck.sum()),
        "python_CK_positive": int(p_ck.sum()),
        "qupath_immune_CD8_positive": int((~q_ck & q_cd8).sum()),
        "python_immune_CD8_positive": int((~p_ck & p_cd8).sum()),
        "qupath_edges": len(qe),
        "python_edges": len(pe),
        "lost_edges": len(lost),
        "gained_edges": len(gained),
        "qupath_CK_with_CD8_neighbor": len(q_covered),
        "python_CK_with_CD8_neighbor": len(p_covered),
        "CK_neighbor_status_flips": len(
            q_covered.symmetric_difference(p_covered)
        ),
        "CK_coverage_delta": len(p_covered) - len(q_covered)
    })

    for status, edge_set in [
        ("lost", lost), ("gained", gained)
    ]:
        for ck_i, cd8_i in sorted(edge_set):
            changed_edges.append({
                "radius_um": radius,
                "status": status,
                "CK_reference_row": int(ref_rows[ck_i]),
                "CD8_reference_row": int(ref_rows[cd8_i]),
                "CK_cell_id": str(day4.iloc[ck_i]["cell_id"]),
                "CD8_cell_id": str(day4.iloc[cd8_i]["cell_id"]),
                "distance_um": float(
                    np.linalg.norm(
                        coords[ck_i] - coords[cd8_i]
                    )
                )
            })

summary = pd.DataFrame(rows)
summary.to_csv(
    OUT / "ck_cd8_spatial_robustness.csv",
    index=False
)

pd.DataFrame(
    changed_edges,
    columns=[
        "radius_um", "status",
        "CK_reference_row", "CD8_reference_row",
        "CK_cell_id", "CD8_cell_id", "distance_um"
    ]
).to_csv(
    OUT / "ck_cd8_changed_spatial_edges.csv",
    index=False
)

print("\n========== SPATIAL ROBUSTNESS ==========")
print(summary.to_string(index=False))

print("\n========== PRIMARY 30 UM AUDIT ==========")
r30 = summary[summary["radius_um"] == 30].iloc[0]
for col in [
    "qupath_edges", "python_edges",
    "lost_edges", "gained_edges",
    "qupath_CK_with_CD8_neighbor",
    "python_CK_with_CD8_neighbor",
    "CK_neighbor_status_flips"
]:
    print(f"{col}: {int(r30[col])}")

print("\n========== EXISTING DAY4 FILE ==========")
existing = ROOT / (
    "results/day4/radius_sensitivity/"
    "radius_30/ck_cell_neighborhoods.csv"
)
if existing.exists():
    e = pd.read_csv(existing)
    print("Rows:", len(e))
    print("Columns:", list(e.columns))

print("\nSTAGE_6H_SPATIAL_AUDIT: COMPLETE")
print("ORIGINAL_SPATIAL_PIPELINE_REPRODUCTION: PENDING")
