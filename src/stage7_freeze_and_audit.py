#!/usr/bin/env python3
"""Stage 7: non-destructive reproducibility freeze and biology review queue.
Run from QuPath mIF Cell Quantification & Spatial Analysis project root.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = Path.cwd()
OUT = ROOT / "results/stage7"
OUT.mkdir(parents=True, exist_ok=True)

REQUIRED = {
    "raw_image": "data/raw/LuCa-7color.tif",
    "original_qupath_objects": "LuCa-7color.tif - resolution #1.qpdata",
    "qupath_reference": "LuCa7_roi01_cells.csv.txt",
    "original_command_history": "qupath_mif_step5_reproducible/provenance/qupath_0.7_original_command_history.txt",
    "roi_export": "results/step6/LuCa7_step6_original_cell_rois.tsv",
    "roi_reference_mapping": "results/step6/cell_roi_reference_mapping.csv",
    "independent_intensity": "results/step6/pixel_intensity_comparison_full.csv",
    "marker_robustness": "results/step6/marker_robustness_full.csv",
    "day4_intensities": "results/day4/LuCa7_full_AF_1310_20261030_190648.csv",
    "day4_mapping": "results/step6/day4_reference_mapping.csv",
    "original_30um_neighbors": "results/day4/radius_sensitivity/radius_30/ck_cell_neighborhoods.csv",
    "reproduced_30um_cells": "results/step6/day4_exact_30um_cell_audit.csv",
    "spatial_sensitivity": "results/step6/day4_exact_spatial_robustness.csv",
}
OPTIONAL = {
    "alternative_reference": "LuCa7_roi01_cells.txt",
    "project_copy": "qupath_project_repro/project.qpproj",
    "project_object_copy": "qupath_project_repro/data/1/data.qpdata",
    "original_roi_geojson": "results/step5/original_roi_reference.geojson",
    "phenotype_discordances": "results/step6/ck_cd8_joint_phenotype_sensitivity.csv",
    "marker_threshold_sensitivity": "results/step6/phenotyping_threshold_sensitivity.csv",
    "mask_count_discrepancy": "results/step6/nucleus_pixel_count_discrepancy_full.csv",
}
missing = [f"{role}: {rel}" for role, rel in REQUIRED.items() if not (ROOT / rel).is_file()]
print("========== STAGE 7 INPUT AUDIT ==========")
if missing:
    print("MISSING REQUIRED FILES:")
    for item in missing:
        print(" -", item)
    raise SystemExit("STAGE7_INPUT_AUDIT: BLOCKED; no PASS claimed")
print("Required files:", len(REQUIRED), "/", len(REQUIRED))

# Selected source scripts are frozen as copies; no raw image or data copy.
script_files = sorted(ROOT.glob("step6*.py"))
self_path = ROOT / "stage7_freeze_and_audit.py"
if self_path.is_file():
    script_files.append(self_path)
for base in [ROOT, ROOT / "scripts", ROOT / "qupath_mif_step5_reproducible"]:
    if base.is_dir():
        script_files.extend(base.glob("*.groovy"))
script_files = sorted({p.resolve() for p in script_files if p.is_file()})
for src in script_files:
    relative = src.relative_to(ROOT)
    dest = OUT / "source_snapshots" / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)

# Hash only specified scientific inputs and source scripts, not an entire project.
def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

manifest = []
for required, sources in [(True, REQUIRED), (False, OPTIONAL)]:
    for role, rel in sources.items():
        path = ROOT / rel
        manifest.append({
            "role": role,
            "required": required,
            "relative_path": rel,
            "exists": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else None,
            "sha256": sha256(path) if path.is_file() else None,
        })
for path in script_files:
    manifest.append({
        "role": "source_code", "required": False,
        "relative_path": str(path.relative_to(ROOT)),
        "exists": True, "bytes": path.stat().st_size,
        "sha256": sha256(path),
    })
pd.DataFrame(manifest).to_csv(OUT / "manifest_sha256.csv", index=False)

def pkg(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None

env = {
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "python": sys.version,
    "os": platform.platform(),
    "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
    "packages": {name: pkg(name) for name in
                 ["numpy", "pandas", "scipy", "Shapely", "tifffile", "scikit-image"]},
}
(OUT / "environment.json").write_text(
    json.dumps(env, indent=2, ensure_ascii=False), encoding="utf-8"
)

ref = pd.read_csv(ROOT / REQUIRED["qupath_reference"], sep=None, engine="python")
day = pd.read_csv(ROOT / REQUIRED["day4_intensities"])
map4 = pd.read_csv(ROOT / REQUIRED["day4_mapping"]).sort_values("day4_row")
old = pd.read_csv(ROOT / REQUIRED["original_30um_neighbors"])
audit = pd.read_csv(ROOT / REQUIRED["reproduced_30um_cells"])
spatial = pd.read_csv(ROOT / REQUIRED["spatial_sensitivity"])
intensity = pd.read_csv(ROOT / REQUIRED["independent_intensity"])
robust = pd.read_csv(ROOT / REQUIRED["marker_robustness"])

assert len(ref) == len(day) == len(map4) == 1310, "Cell count mismatch"
assert len(intensity) == 31440, "Expected 31,440 intensity comparisons"
assert len(robust) == 24, "Expected 24 marker-compartment comparisons"
assert map4["reference_row"].is_unique and map4["day4_row"].is_unique
assert np.array_equal(map4["day4_row"].to_numpy(int), np.arange(1310))
assert np.isfinite(map4["distance_um"]).all()
assert (map4["distance_um"] <= 0.5).all(), "Day4 mapping tolerance exceeded"
assert np.isfinite(intensity[["qupath_mean", "python_mean"]].to_numpy(float)).all()

markers = {
    "PDL1": "Cell: PDL1 (Opal 520) mean",
    "CD8": "Cell: CD8 (Opal 540) mean",
    "FoxP3": "Nucleus: FoxP3 (Opal 570) mean",
    "CD68": "Cytoplasm: CD68 (Opal 620) mean",
    "PD1": "Cell: PD1 (Opal 650) mean",
    "CK": "Cell: CK (Opal 690) mean",
    "DAPI": "Nucleus: DAPI mean",
    "Autofluorescence": "Cell: Autofluorescence mean",
}
indices = map4["reference_row"].to_numpy(int)
source_rows = []
for name, col in markers.items():
    values_day = pd.to_numeric(day[name], errors="coerce").to_numpy(float)
    values_original = pd.to_numeric(ref.iloc[indices][col], errors="coerce").to_numpy(float)
    assert np.isfinite(values_day).all() and np.isfinite(values_original).all()
    difference = values_day - values_original
    source_rows.append({
        "marker": name,
        "qupath_reference_column": col,
        "max_abs_difference": float(np.max(np.abs(difference))),
        "mae": float(np.mean(np.abs(difference))),
        "matches_within_1e-4": bool(np.all(np.abs(difference) <= 1e-4)),
    })
source_df = pd.DataFrame(source_rows)
source_df.to_csv(OUT / "day4_marker_provenance.csv", index=False)
assert source_df["matches_within_1e-4"].all(), "Marker provenance requires review"

original_expected = (day["CK"] >= 3.5) & (day["CD8"] < 0.5)
assert len(old) == 439 and original_expected.sum() == 439
assert day["cell_id"].is_unique and old["cell_id"].is_unique
assert set(old["cell_id"].astype(str)) == set(day.loc[original_expected, "cell_id"].astype(str))
assert len(audit) == 439 and audit["cell_id"].is_unique
original_check = old[["cell_id", "n_CD8_only_within_radius"]].merge(
    audit[["cell_id", "qupath_recomputed_count"]], on="cell_id", validate="one_to_one"
)
assert len(original_check) == 439
assert np.array_equal(
    original_check["n_CD8_only_within_radius"].to_numpy(float),
    original_check["qupath_recomputed_count"].to_numpy(float)
), "Original neighbor counts not reproduced"

r30 = spatial.loc[spatial["radius_um"] == 30]
assert len(r30) == 1
r30 = r30.iloc[0]

# Create a deterministic review queue, WITHOUT generating truth labels.
# Include 3 discordant cells, double-high cases, CD8 threshold-adjacent cases,
# large per-marker measurement differences, plus a reproducible random subset.
selected = {}
def add(reference_row, reason):
    selected.setdefault(int(reference_row), set()).add(reason)
if (ROOT / OPTIONAL["phenotype_discordances"]).is_file():
    phen = pd.read_csv(ROOT / OPTIONAL["phenotype_discordances"])
    for rr in phen.loc[phen["changed"].astype(str).str.lower() == "true", "reference_row"]:
        add(rr, "CK_CD8_discordant")
for rr in indices[(day["CK"].to_numpy() >= 3.5) & (day["CD8"].to_numpy() >= 0.5)]:
    add(rr, "CK_CD8_double_high")
for i in np.argsort(np.abs(day["CD8"].to_numpy(float) - 0.5))[:25]:
    add(indices[i], "near_CD8_0.5")
for comp, marker in [("Nucleus", "FoxP3 (Opal 570)"),
                     ("Cytoplasm", "CD68 (Opal 620)"),
                     ("Cell", "PDL1 (Opal 520)"),
                     ("Cell", "PD1 (Opal 650)")]:
    subset = intensity.loc[(intensity["compartment"] == comp) &
                           (intensity["marker"] == marker)].copy()
    subset["recomputed_abs_error"] = np.abs(subset["python_mean"] - subset["qupath_mean"])
    for rr in subset.nlargest(8, "recomputed_abs_error")["reference_row"]:
        add(rr, "large_remeasurement_delta_" + marker.split(" (")[0])
if len(selected) < 100:
    rng = np.random.default_rng(20261030)
    candidates = np.array(sorted(set(indices) - set(selected)), dtype=int)
    for rr in rng.choice(candidates, size=min(100 - len(selected), len(candidates)), replace=False):
        add(rr, "random_background")
# If overlapping categories yield >100 mandatory cases, retain all to avoid dropping rare cases.
row_by_reference = {int(rr): i for i, rr in enumerate(indices)}
review = []
for rr, reasons in sorted(selected.items()):
    i = row_by_reference[rr]
    v = day.iloc[i]
    review.append({
        "reference_row": rr,
        "cell_id": str(v["cell_id"]),
        "x_um": float(v["x_um"]),
        "y_um": float(v["y_um"]),
        "selection_reason": ";".join(sorted(reasons)),
        "CK": float(v["CK"]),
        "CD8": float(v["CD8"]),
        "FoxP3": float(v["FoxP3"]),
        "CD68": float(v["CD68"]),
    })
review_df = pd.DataFrame(review)
review_df.to_csv(OUT / "biology_review_key.csv", index=False)
blind = review_df[["cell_id", "x_um", "y_um"]].copy()
for col in ["reviewer", "CK_manual", "CD8_manual", "FoxP3_manual", "CD68_manual",
            "segmentation_quality", "ambiguous", "notes"]:
    blind[col] = ""
blind.to_csv(OUT / "biology_review_blinded_template.csv", index=False)

summary = {
    "required_files": len(REQUIRED),
    "source_scripts_snapshotted": len(script_files),
    "reference_cells": len(ref),
    "independent_intensity_measurements": len(intensity),
    "marker_compartment_combinations": len(robust),
    "day4_source_columns_verified": int(source_df["matches_within_1e-4"].sum()),
    "original_day4_neighbor_counts_reproduced": len(original_check),
    "primary_radius_um": 30,
    "original_30um_edges": int(r30["qupath_edges"]),
    "remeasured_30um_edges": int(r30["python_edges"]),
    "original_30um_covered_cells": int(r30["qupath_covered"]),
    "remeasured_30um_covered_cells": int(r30["python_covered"]),
    "30um_coverage_status_flips": int(r30["coverage_flips"]),
    "manual_review_candidates": len(review_df),
    "biological_validation_status": "PENDING",
    "AF_correction_status": "NOT_DEMONSTRATED_IN_DAY4_MARKER_VALUES",
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

readme = f"""# QuPath mIF reproducibility package — Stage 7

Generated: {env['generated_utc']} (UTC)

## Scope and evidence
- Original reference: 1,310 matched cells; 8 fluorescent channels.
- Original Day 4 neighborhood counts: {len(original_check)}/{len(original_check)} reproduced at 30 µm for CK+/CD8− centers.
- Independent TIFF + exported WKT ROI remeasurement: {len(intensity):,} intensity comparisons, **not exact pixel-level numerical reproduction**.
- Primary 30 µm: CK-to-CK−/CD8+ edges {int(r30['qupath_edges'])} → {int(r30['python_edges'])}; covered CK+/CD8− cells {int(r30['qupath_covered'])} → {int(r30['python_covered'])}; coverage flips {int(r30['coverage_flips'])}.
- Day 4 source columns match their original QuPath measurements to within 1e-4: {int(source_df['matches_within_1e-4'].sum())}/8.

## Interpretation boundaries
1. `LuCa7_full_AF...csv` has marker values nearly identical to the original QuPath measurements. **No numerical AF correction is demonstrated by this file**, despite the filename; separate upstream processing/compensation provenance must be checked.
2. The CK ≥ 3.5 and CD8 ≥ 0.5 gates are **exploratory**, not independently validated biological thresholds.
3. Original Watershed label pixels and stored smoothed polygon ROIs can differ; do not claim pixel-exact reproduction.
4. Same-slide robustness does not establish biological phenotyping accuracy, cross-slide generalization, or biological significance.

## Package files
- `manifest_sha256.csv`: curated input/output/source code hashes (paths relative to project root).
- `environment.json`: runtime/library versions.
- `source_snapshots/`: frozen copies of discovered Stage 6 Python and Groovy scripts only.
- `day4_marker_provenance.csv`: source compartment and numeric match for each Day 4 marker.
- `biology_review_blinded_template.csv`: candidate IDs/coordinates for blinded manual review; **empty annotations are intentional**.
- `biology_review_key.csv`: **unblinded** selection reasons/intensity reference, keep hidden from reviewers.
- `summary.json`: machine-readable engineering results.

## Reproducibility procedure
1. Verify file integrity against `manifest_sha256.csv` before rerunning analyses.
2. Record exact QuPath application and extension versions plus detection script and input ROI selection.
3. Use **original QuPath mean measurements** for reproducing original Day 4 results; use independent Python ROI means only for measurement sensitivity tests.
4. Reproduce the original 439-center / 30 µm neighborhood analysis from preserved `day4_exact_30um_cell_audit.csv` and the frozen script.
5. Independently review phenotype threshold provenance and any AF/spectral-unmixing or background compensation before biological claims.

## Biological validation (not yet completed)
- Verify original channel acquisition, spectral unmixing, controls, and compensation provenance.
- Independently assess cell morphology/segmentation and marker localization for the blinded review queue.
- Establish marker thresholds using documented biological controls or a prespecified validated gating procedure.
- Recompute multi-marker phenotypes and spatial endpoints under the validated pipeline; expand to independent slides/cases.

## Data policy
No raw TIFF, QuPath object file, or patient data are copied into this package. Do not upload these files or derived cell-level data publicly without confirming rights, consent, and de-identification.
"""
(OUT / "README_reproducibility.md").write_text(readme, encoding="utf-8")

print("\n========== STAGE 7 FREEZE ==========")
print("Required inputs verified:", len(REQUIRED))
print("Source scripts snapshotted:", len(script_files))
print("Manifest entries:", len(manifest))
print("Day 4 marker sources verified:", len(source_df), "/ 8")
print("Original 30um neighbor counts reproduced:", len(original_check), "/ 439")
print("30um edges:", summary["original_30um_edges"], "->", summary["remeasured_30um_edges"])
print("30um coverage:", summary["original_30um_covered_cells"], "->", summary["remeasured_30um_covered_cells"])
print("30um coverage flips:", summary["30um_coverage_status_flips"])
print("Blinded biology review candidates:", len(review_df))
print("AF correction: NOT ESTABLISHED")
print("Biological validation: PENDING")
print("Output:", OUT)
print("STAGE_7A_ENGINEERING_FREEZE: PASS")
print("STAGE_7B_BIOLOGICAL_VALIDATION: PENDING")
