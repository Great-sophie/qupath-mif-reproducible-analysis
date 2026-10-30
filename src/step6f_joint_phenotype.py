from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("results/step6")
d = pd.read_csv(
    OUT / "pixel_intensity_comparison_full.csv"
)

markers = ["CK (Opal 690)", "CD8 (Opal 540)"]
g = d[
    (d["compartment"] == "Cell") &
    d["marker"].isin(markers)
]

ref = g.pivot(
    index="reference_row",
    columns="marker",
    values="qupath_mean"
)

py = g.pivot(
    index="reference_row",
    columns="marker",
    values="python_mean"
)

assert len(ref) == len(py) == 1310
assert ref.index.equals(py.index)
assert not ref.isna().any().any()
assert not py.isna().any().any()

def classify(table):
    ck = table["CK (Opal 690)"] >= 3.5
    cd8 = table["CD8 (Opal 540)"] >= 0.5

    return np.select(
        [
            ck & cd8,
            ck & ~cd8,
            ~ck & cd8
        ],
        [
            "CK+/CD8+",
            "CK+/CD8-",
            "CK-/CD8+"
        ],
        default="CK-/CD8-"
    )

a = classify(ref)
b = classify(py)
changed = a != b

result = pd.DataFrame({
    "reference_row": ref.index.to_numpy(),
    "qupath_phenotype": a,
    "python_phenotype": b,
    "changed": changed,
    "qupath_CK": ref[markers[0]].to_numpy(),
    "python_CK": py[markers[0]].to_numpy(),
    "qupath_CD8": ref[markers[1]].to_numpy(),
    "python_CD8": py[markers[1]].to_numpy()
})

result.to_csv(
    OUT / "ck_cd8_joint_phenotype_sensitivity.csv",
    index=False
)

labels = [
    "CK+/CD8+",
    "CK+/CD8-",
    "CK-/CD8+",
    "CK-/CD8-"
]

confusion = pd.crosstab(
    pd.Categorical(a, categories=labels),
    pd.Categorical(b, categories=labels),
    rownames=["QuPath"],
    colnames=["Python"],
    dropna=False
)

confusion.to_csv(
    OUT / "ck_cd8_joint_phenotype_confusion.csv"
)

print("========== STAGE 6F JOINT PHENOTYPE ==========")
print("Cells:", len(result))
print("Agreement:", int((~changed).sum()))
print("Changed phenotypes:", int(changed.sum()))
print("Agreement rate:", round(100 * (~changed).mean(), 4))
print("Flip rate:", round(100 * changed.mean(), 4))

print("\n========== CONFUSION MATRIX ==========")
print(confusion.to_string())

print("\n========== CHANGED CELLS ==========")
print(
    result[result["changed"]]
    .round(5).head(25).to_string(index=False)
)

print("\n========== DAY 4 DATA INSPECTION ==========")
for name in [
    "results/day4/LuCa7_full_AF_1310_20261030_190648.csv",
    "results/day4/CK_CD8_double_high_review.csv",
    "results/day3/marker_qc/marker_summary.csv"
]:
    path = Path(name)
    if path.exists():
        x = pd.read_csv(path)
        print("\nFILE:", name)
        print("Rows:", len(x))
        print("Columns:", list(x.columns))
    else:
        print("\nNOT FOUND:", name)

print("\nSTAGE_6F_EXPLORATORY_GATE: COMPLETE")
print("BIOLOGICAL_PHENOTYPING_VALIDATION: PENDING")
