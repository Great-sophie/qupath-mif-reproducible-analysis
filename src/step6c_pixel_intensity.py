from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import tifffile
from shapely import wkt

try:
    from shapely import contains_xy
    def inside(geometry, xx, yy):
        return contains_xy(geometry, xx, yy)
except ImportError:
    from shapely.geometry import Point
    from shapely.prepared import prep
    def inside(geometry, xx, yy):
        prepared = prep(geometry)
        result = [
            prepared.contains(Point(float(x), float(y)))
            for x, y in zip(xx.ravel(), yy.ravel())
        ]
        return np.asarray(result).reshape(xx.shape)

parser = argparse.ArgumentParser()
parser.add_argument("--limit", type=int, default=None)
args = parser.parse_args()

ROOT = Path.cwd()
OUT = ROOT / "results/step6"
OUT.mkdir(parents=True, exist_ok=True)

ref = pd.read_csv(
    ROOT / "LuCa7_roi01_cells.csv.txt",
    sep=None,
    engine="python"
)

rois = pd.read_csv(
    OUT / "LuCa7_step6_original_cell_rois.tsv",
    sep="\t"
)

mapping = pd.read_csv(
    OUT / "cell_roi_reference_mapping.csv"
)

channels = [
    "PDL1 (Opal 520)",
    "CD8 (Opal 540)",
    "FoxP3 (Opal 570)",
    "CD68 (Opal 620)",
    "PD1 (Opal 650)",
    "CK (Opal 690)",
    "DAPI",
    "Autofluorescence"
]

with tifffile.TiffFile(ROOT / "data/raw/LuCa-7color.tif") as tif:
    series = tif.series[0]
    assert series.axes == "CYX", series.axes
    assert series.shape[0] == 8, series.shape
    image = series.asarray()

assert image.dtype == np.float32
assert len(ref) == len(rois) == len(mapping) == 1310

assert mapping["roi_row"].is_unique
assert mapping["reference_row"].is_unique

mapping = mapping.sort_values("roi_row").reset_index(drop=True)

assert np.array_equal(
    mapping["roi_row"].to_numpy(dtype=int),
    np.arange(len(rois))
)

assert set(mapping["reference_row"].astype(int)) == set(range(1310))

required = [
    f"{comp}: {marker} mean"
    for comp in ["Nucleus", "Cell", "Cytoplasm"]
    for marker in channels
]
assert all(c in ref.columns for c in required)

n = len(rois) if args.limit is None else min(args.limit, len(rois))
assert n > 0

H, W = image.shape[1:]
records = []
geometry_qc = []

for i in range(n):
    row = rois.iloc[i]
    ref_idx = int(mapping.iloc[i]["reference_row"])

    cell = wkt.loads(row["cell_wkt"])
    nucleus = wkt.loads(row["nucleus_wkt"])

    if (
        cell.is_empty or nucleus.is_empty
        or not cell.is_valid or not nucleus.is_valid
        or cell.area <= 0 or nucleus.area <= 0
    ):
        raise ValueError(f"Invalid geometry at ROI row {i}")

    x0 = max(0, int(np.floor(cell.bounds[0] - 2)))
    x1 = min(W, int(np.ceil(cell.bounds[2] + 2)))
    y0 = max(0, int(np.floor(cell.bounds[1] - 2)))
    y1 = min(H, int(np.ceil(cell.bounds[3] + 2)))

    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"Empty image window at ROI row {i}")

    xs = np.arange(x0, x1) + 0.5
    ys = np.arange(y0, y1) + 0.5
    xx, yy = np.meshgrid(xs, ys)

    cmask = inside(cell, xx, yy)
    nmask = inside(nucleus, xx, yy)
    ymask = cmask & ~nmask

    masks = {
        "Nucleus": nmask,
        "Cell": cmask,
        "Cytoplasm": ymask
    }

    geometry_qc.append({
        "roi_row": i,
        "reference_row": ref_idx,
        "cell_pixels": int(cmask.sum()),
        "nucleus_pixels": int(nmask.sum()),
        "cytoplasm_pixels": int(ymask.sum()),
        "nucleus_pixels_outside_cell": int(
            (nmask & ~cmask).sum()
        )
    })

    patch = image[:, y0:y1, x0:x1]

    for compartment, mask in masks.items():
        for ch, marker in enumerate(channels):
            values = patch[ch][mask]
            values = values[np.isfinite(values)]

            reproduced = (
                float(np.mean(values, dtype=np.float64))
                if values.size else np.nan
            )

            column = f"{compartment}: {marker} mean"
            original = float(ref.iloc[ref_idx][column])

            records.append({
                "roi_row": i,
                "reference_row": ref_idx,
                "compartment": compartment,
                "marker": marker,
                "n_pixels": int(values.size),
                "qupath_mean": original,
                "python_mean": reproduced,
                "signed_error": reproduced - original,
                "abs_error": abs(reproduced - original)
            })

    if (i + 1) % 100 == 0:
        print(f"Processed: {i + 1}/{n}", flush=True)

result = pd.DataFrame(records)
gqc = pd.DataFrame(geometry_qc)

suffix = "_pilot" if args.limit is not None else "_full"

result.to_csv(
    OUT / f"pixel_intensity_comparison{suffix}.csv",
    index=False
)
gqc.to_csv(
    OUT / f"roi_pixel_geometry_qc{suffix}.csv",
    index=False
)

def summarize(g):
    a = g["qupath_mean"].to_numpy(float)
    b = g["python_mean"].to_numpy(float)
    valid = np.isfinite(a) & np.isfinite(b)

    if not valid.any():
        return pd.Series({
            "n": 0, "MAE": np.nan,
            "RMSE": np.nan, "max_abs": np.nan,
            "bias": np.nan, "pearson_r": np.nan
        })

    a, b = a[valid], b[valid]
    d = b - a

    corr = (
        float(np.corrcoef(a, b)[0, 1])
        if len(a) > 1 and np.std(a) > 0
        and np.std(b) > 0 else np.nan
    )

    return pd.Series({
        "n": len(a),
        "MAE": float(np.mean(np.abs(d))),
        "RMSE": float(np.sqrt(np.mean(d*d))),
        "max_abs": float(np.max(np.abs(d))),
        "bias": float(np.mean(d)),
        "pearson_r": corr
    })

summary = (
    result.groupby(
        ["compartment", "marker"], sort=False
    )
    .apply(summarize)
    .reset_index()
)

summary.to_csv(
    OUT / f"pixel_intensity_summary{suffix}.csv",
    index=False
)

print("\n========== STAGE 6C PIXEL AUDIT ==========")
print("Cells:", n)
print("Measurements:", len(result))
print("Valid measurements:",
      int(np.isfinite(result["python_mean"]).sum()))
print("Empty compartments:",
      int((gqc[
          ["cell_pixels", "nucleus_pixels", "cytoplasm_pixels"]
      ] == 0).any(axis=1).sum()))
print("Nucleus pixels outside cell:",
      int(gqc["nucleus_pixels_outside_cell"].sum()))

print("\n========== ERROR BY COMPARTMENT ==========")
print(
    result.groupby("compartment", sort=False)
    .apply(summarize)
    .round(6)
    .to_string()
)

print("\n========== MARKER-LEVEL SUMMARY ==========")
print(summary.round(6).to_string(index=False))

print("\n========== LARGEST DISCREPANCIES ==========")
print(
    result.sort_values("abs_error", ascending=False)
    .head(10)
    .round(6)
    .to_string(index=False)
)

print("\nSaved output:", OUT)
print("STAGE_6C_COMPARISON: COMPLETE")
print("QUANTITATIVE_REPRODUCTION: REQUIRES_REVIEW")
