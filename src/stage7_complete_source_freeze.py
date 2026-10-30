from pathlib import Path
import hashlib
import shutil
import csv

ROOT = Path.cwd()
OUT = ROOT / "results/stage7"
DEST = OUT / "source_snapshots"
DEST.mkdir(parents=True, exist_ok=True)

SOURCE_DIRS = [
    ROOT / "qupath_mif_day3_qc",
    ROOT / "qupath_mif_day4_spatial",
    ROOT / "qupath_mif_step5_reproducible",
]

files = []

for directory in SOURCE_DIRS:
    if directory.exists():
        files.extend(
            p for p in directory.rglob("*")
            if p.is_file()
            and p.suffix.lower() in {".py", ".groovy"}
        )

# Include Stage 6 and Stage 7 analysis scripts
files.extend(ROOT.glob("step6*.py"))
files.extend(ROOT.glob("stage7*.py"))

files = sorted(set(files))

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

records = []

for source in files:
    rel = source.relative_to(ROOT)
    target = DEST / rel
    target.parent.mkdir(parents=True, exist_ok=True)

    original_hash = sha256(source)

    if target.exists():
        if sha256(target) != original_hash:
            raise RuntimeError(
                f"Existing frozen copy differs: {rel}"
            )
    else:
        shutil.copy2(source, target)

    frozen_hash = sha256(target)
    assert original_hash == frozen_hash

    records.append({
        "relative_path": str(rel),
        "sha256": original_hash,
        "snapshot_path": str(target.relative_to(ROOT)),
        "verified": True
    })

manifest = OUT / "source_manifest_sha256.csv"

with open(manifest, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "relative_path", "sha256",
            "snapshot_path", "verified"
        ]
    )
    writer.writeheader()
    writer.writerows(records)

print("========== STAGE 7 SOURCE FREEZE ==========")
print("Scripts discovered:", len(files))
print("Scripts verified:", len(records))
print("Hash mismatches:", 0)
print("Source manifest:", manifest)
print("Snapshot directory:", DEST)

print("\n========== ORIGINAL 9 SCRIPTS ==========")

expected = [
    "qupath_mif_day3_qc/01_marker_distribution_qc.py",
    "qupath_mif_day3_qc/02_audit_autofluorescence.groovy",
    "qupath_mif_day4_spatial/01_export_full_cell_AF.groovy",
    "qupath_mif_day4_spatial/02_spatial_neighborhood_qc.py",
    "qupath_mif_day4_spatial/tests/test_spatial.py",
    "qupath_mif_step5_reproducible/scripts/00_export_selected_roi.groovy",
    "qupath_mif_step5_reproducible/scripts/01_cell_detection_safe.groovy",
    "qupath_mif_step5_reproducible/scripts/02_readonly_qc_export.groovy",
    "qupath_mif_step5_reproducible/scripts/03_compare_reproduction.py",
]

frozen = {r["relative_path"] for r in records}

for name in expected:
    print(
        "PASS" if name in frozen else "MISSING",
        name
    )

passed = all(x in frozen for x in expected)

print(
    "\nSTAGE_7_SOURCE_FREEZE:",
    "PASS" if passed else "INCOMPLETE"
)
