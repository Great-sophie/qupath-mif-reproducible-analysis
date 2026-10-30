from pathlib import Path
import pandas as pd
import numpy as np
import re

ROOT = Path.cwd()
OUT = ROOT / "results/step6"
OUT.mkdir(parents=True, exist_ok=True)

REF = ROOT / "LuCa7_roi01_cells.csv.txt"
ALT = ROOT / "LuCa7_roi01_cells.txt"
IMG = ROOT / "data/raw/LuCa-7color.tif"

df = pd.read_csv(REF, sep=None, engine="python")

print("\n========== STAGE 6A: REFERENCE AUDIT ==========")
print("Reference rows:", len(df))
print("Reference columns:", len(df.columns))
print("Duplicate rows:", int(df.duplicated().sum()))

if ALT.exists():
    other = pd.read_csv(ALT, sep=None, engine="python")
    print("Two reference exports identical:", df.equals(other))

print("\n========== IDENTIFIER / COORDINATE COLUMNS ==========")
for c in df.columns:
    if any(x in c.lower() for x in [
        "id", "centroid", "position", "classification"
    ]):
        print(c)

print("\n========== INTENSITY MEASUREMENTS ==========")

metrics = ["mean", "sum", "std dev", "max", "min", "range"]
measurements = []

for col in df.columns:
    m = re.match(
        r"^([^:]+): (.+) (mean|sum|std dev|max|min|range)$",
        str(col)
    )
    if not m:
        continue

    compartment, marker, metric = m.groups()
    if compartment not in (
        "Nucleus", "Cell", "Cytoplasm", "Membrane"
    ):
        continue

    measurements.append({
        "column": col,
        "compartment": compartment,
        "marker": marker,
        "metric": metric
    })

meta = pd.DataFrame(measurements)

if meta.empty:
    print("No recognized intensity columns.")
else:
    coverage = (
        meta.groupby(["compartment", "marker"])["metric"]
        .apply(lambda x: ",".join(sorted(x)))
        .reset_index()
    )
    print(coverage.to_string(index=False))

    coverage.to_csv(
        OUT / "intensity_column_coverage.csv",
        index=False
    )

    print("\n========== MARKER MEAN QC ==========")

    results = []

    for item in measurements:
        if item["metric"] != "mean":
            continue

        values = pd.to_numeric(
            df[item["column"]], errors="coerce"
        ).to_numpy(dtype=float)

        finite = values[np.isfinite(values)]

        results.append({
            "compartment": item["compartment"],
            "marker": item["marker"],
            "n_total": len(values),
            "n_valid": len(finite),
            "n_missing": len(values) - len(finite),
            "n_negative": int((finite < 0).sum()),
            "median": np.median(finite) if len(finite) else np.nan,
            "p05": np.percentile(finite, 5) if len(finite) else np.nan,
            "p95": np.percentile(finite, 95) if len(finite) else np.nan,
        })

    qc = pd.DataFrame(results)
    qc.to_csv(
        OUT / "marker_intensity_reference_qc.csv",
        index=False
    )
    print(qc.to_string(index=False))

print("\n========== IMAGE METADATA ==========")

try:
    import tifffile

    with tifffile.TiffFile(IMG) as tif:
        for i, series in enumerate(tif.series[:5]):
            print(
                "Series:", i,
                "Shape:", series.shape,
                "Axes:", series.axes,
                "Dtype:", series.dtype,
                "Levels:", len(series.levels)
            )

except Exception as e:
    print("Image metadata unavailable:", e)

print("\n========== ROI / MASK AVAILABILITY ==========")

extensions = {".geojson", ".qpdata", ".qpproj", ".npy", ".npz"}

found = [
    p for p in ROOT.rglob("*")
    if p.is_file() and p.suffix.lower() in extensions
]

for p in sorted(found)[:40]:
    print(p.relative_to(ROOT))

print("Candidate ROI/mask/project files:", len(found))

print("\nSTAGE 6A AUDIT COMPLETE")
print("Output:", OUT)
