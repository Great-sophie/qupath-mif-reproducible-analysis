#!/usr/bin/env python3
"""QuPath LuCa7 marker-intensity QC before phenotype thresholding.

Only computes descriptive statistics and diagnostic plots. Does NOT classify
cells positive/negative or claim validated biomarkers.
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path, help='QuPath per-cell CSV')
    parser.add_argument('--out', required=True, type=Path, help='New output folder')
    args = parser.parse_args()

    try:
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(f'Missing dependency: {exc}. Use conda environment with numpy pandas matplotlib.')

    required = ['cell_id', 'x_um', 'y_um', 'PDL1', 'CD8', 'FoxP3', 'CD68', 'PD1', 'CK', 'DAPI']
    if not args.input.is_file():
        raise SystemExit(f'File does not exist: {args.input}')

    df = pd.read_csv(args.input)
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit('CSV missing fields: ' + ', '.join(missing))
    for k in required[1:]:
        df[k] = pd.to_numeric(df[k], errors='coerce')

    args.out.mkdir(parents=True, exist_ok=True)
    marker_cols = required[3:]
    rows = []
    for marker in marker_cols:
        vals = df[marker].to_numpy(dtype=float)
        valid = vals[np.isfinite(vals)]
        if len(valid):
            qs = np.quantile(valid, [0.05, 0.25, 0.5, 0.75, 0.95])
        else:
            qs = [float('nan')] * 5
        rows.append(dict(marker=marker, n_total=len(vals), n_valid=len(valid),
                         n_nonfinite=int((~np.isfinite(vals)).sum()),
                         n_negative=int((valid < 0).sum()),
                         P05=float(qs[0]), P25=float(qs[1]),
                         P50=float(qs[2]), P75=float(qs[3]), P95=float(qs[4])))
        if len(valid):
            fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
            ax.hist(valid, bins=55)
            ax.axvline(qs[2], linestyle='--', label=f'P50={qs[2]:.3g}')
            ax.axvline(qs[4], linestyle=':', label=f'P95={qs[4]:.3g}')
            ax.set(xlabel=f'{marker}: measured mean fluorescence intensity',
                   ylabel='Cells', title=f'{marker} — descriptive distribution (NOT positivity calls)')
            ax.legend()
            fig.savefig(args.out / f'hist_{marker}.png', dpi=180)
            plt.close(fig)

    pd.DataFrame(rows).to_csv(args.out / 'marker_summary.csv', index=False)

    finite_xy = np.isfinite(df[['x_um', 'y_um']].to_numpy(dtype=float)).all(axis=1)
    if finite_xy.any():
        xy = df.loc[finite_xy, ['x_um', 'y_um']].to_numpy(dtype=float)
        for marker in ['CK', 'CD8', 'FoxP3', 'CD68']:
            x = pd.to_numeric(df.loc[finite_xy, marker], errors='coerce').to_numpy(dtype=float)
            ok = np.isfinite(x)
            if not ok.any():
                continue
            fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
            pts = ax.scatter(xy[ok, 0], xy[ok, 1], c=np.log1p(np.clip(x[ok], 0, None)),
                             s=9, alpha=0.85, linewidths=0)
            fig.colorbar(pts, ax=ax, label=f'log(1 + {marker} intensity)')
            ax.invert_yaxis()  # Image coordinates increase downward
            ax.set(xlabel='Image x (µm)', ylabel='Image y (µm)',
                   title=f'{marker} continuous-intensity map (NOT classified)')
            ax.set_aspect('equal', adjustable='box')
            fig.savefig(args.out / f'spatial_{marker}.png', dpi=180)
            plt.close(fig)

    for m1, m2 in [('CK', 'CD8'), ('CD8', 'PD1'), ('CD68', 'PDL1')]:
        a = df[m1].to_numpy(dtype=float)
        b = df[m2].to_numpy(dtype=float)
        valid = np.isfinite(a) & np.isfinite(b)
        if not valid.any():
            continue
        fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
        ax.scatter(np.log1p(np.clip(a[valid], 0, None)),
                   np.log1p(np.clip(b[valid], 0, None)), s=11, alpha=0.7)
        ax.set(xlabel=f'log(1 + {m1} intensity)',
               ylabel=f'log(1 + {m2} intensity)',
               title=f'{m1} vs {m2} — exploratory co-expression')
        fig.savefig(args.out / f'scatter_{m1}_{m2}.png', dpi=180)
        plt.close(fig)

    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    report = dict(input_file=str(args.input), input_sha256=digest,
                  n_cells=int(len(df)), n_duplicate_ids=int(df.cell_id.duplicated().sum()),
                  n_invalid_xy=int((~finite_xy).sum()),
                  note='NO positivity thresholds estimated or applied. No autofluorescence column in this source CSV.'
                       if 'AF' not in df.columns else
                       'NO positivity thresholds estimated or applied. AF column exists.')
    (args.out / 'audit.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'INPUT: {args.input}')
    print(f'CELLS: {len(df)}')
    print(f'DUPLICATE CELL IDS: {report["n_duplicate_ids"]}')
    print(f'INVALID COORDINATES: {report["n_invalid_xy"]}')
    print('OUTPUT:', args.out)
    print('Created: marker_summary.csv, audit.json, 7 marker histograms, 4 continuous-intensity maps, 3 scatterplots')
    print('NOTE: Distribution quantiles are NOT positive/negative thresholds.')

if __name__ == '__main__':
    main()
