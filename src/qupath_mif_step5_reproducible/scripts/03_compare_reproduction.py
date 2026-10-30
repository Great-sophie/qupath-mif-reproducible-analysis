#!/usr/bin/env python3
"""Compare cell centroid coordinates in original vs reproduced QuPath CSV.

Example:
 python 03_compare_reproduction.py --reference original.csv --reproduced repro.csv
A matching COUNT alone is not sufficient. Accepts x_um/y_um and some QuPath export names.
"""
import argparse
import csv
import math
from pathlib import Path


def read_points(p):
    with Path(p).open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        names = reader.fieldnames or []
        options = [("x_um", "y_um"), ("Centroid X µm", "Centroid Y µm"),
                   ("Centroid X um", "Centroid Y um")]
        coords = next(((x, y) for x, y in options if x in names and y in names), None)
        if coords is None:
            raise ValueError(f"Coordinates not found in {p}; columns: {names[:20]}")
        points = []
        for row in reader:
            try:
                x, y = float(row[coords[0]]), float(row[coords[1]])
                if math.isfinite(x) and math.isfinite(y):
                    points.append((x, y))
            except (ValueError, TypeError):
                continue
    return points


def compare(reference, reproduced, tolerance):
    import numpy as np
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise SystemExit("Install scipy in your environment: conda install scipy") from exc
    ref = np.asarray(reference, dtype=float).reshape(-1, 2)
    rep = np.asarray(reproduced, dtype=float).reshape(-1, 2)
    if len(ref) == 0 or len(rep) == 0:
        raise ValueError("Reference or reproduction has zero valid centroids")
    # Mutual nearest-neighbor matching prevents double-counting in crowded regions.
    d_ref, i_ref = cKDTree(rep).query(ref)
    d_rep, i_rep = cKDTree(ref).query(rep)
    mutual = [(j, int(i_ref[j]), float(d_ref[j])) for j in range(len(ref))
              if d_ref[j] <= tolerance and i_rep[i_ref[j]] == j
              and d_rep[i_ref[j]] <= tolerance]
    errors = [distance for _, _, distance in mutual]
    return {
        "reference_n": len(ref), "reproduced_n": len(rep),
        "count_match": len(ref) == len(rep),
        "mutual_matches": len(mutual),
        "match_fraction_of_reference": len(mutual) / len(ref),
        "median_distance_um": float(np.median(errors)) if errors else None,
        "maximum_distance_um": max(errors) if errors else None,
        "tolerance_um": tolerance,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reproduced", required=True)
    parser.add_argument("--tolerance-um", type=float, default=0.5)
    args = parser.parse_args()
    if args.tolerance_um <= 0:
        parser.error("--tolerance-um must be positive")
    stats = compare(read_points(args.reference), read_points(args.reproduced), args.tolerance_um)
    for key, value in stats.items():
        print(f"{key}: {value}")
    if stats["count_match"] and stats["match_fraction_of_reference"] >= 0.99:
        print("CENTROID_REPRODUCTION_CHECK: PASS (within specified tolerance)")
    else:
        print("CENTROID_REPRODUCTION_CHECK: REVIEW REQUIRED")


if __name__ == "__main__":
    main()
