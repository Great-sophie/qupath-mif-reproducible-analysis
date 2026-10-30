#!/usr/bin/env python3
"""Exploratory 30-um CK/CD8 neighborhood QC for a single QuPath ROI.

No marker thresholds are considered validated. No cell identity is inferred.
Preserves spatial positions and randomizes CD8-high labels among CK-low cells.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CK_T = 3.5
CD8_T = 0.5


def group_masks(df, ck_t=CK_T, cd8_t=CD8_T):
    ck = df['CK'].to_numpy(dtype=float)
    cd8 = df['CD8'].to_numpy(dtype=float)
    finite = np.isfinite(ck) & np.isfinite(cd8)
    return {
        'valid': finite,
        'ck_only': finite & (ck >= ck_t) & (cd8 < cd8_t),
        'cd8_only': finite & (ck < ck_t) & (cd8 >= cd8_t),
        'dual': finite & (ck >= ck_t) & (cd8 >= cd8_t),
        'low_low': finite & (ck < ck_t) & (cd8 < cd8_t),
        'ck_low_eligible': finite & (ck < ck_t),
    }


def local_bin(x, y, eligible, n_bins=4):
    """Divide ROI bounding rectangle into bins to preserve coarse spatial intensity."""
    if eligible.sum() == 0:
        return np.zeros(len(x), dtype=int)
    xx = x[eligible]
    yy = y[eligible]
    xmin, xmax = float(xx.min()), float(xx.max())
    ymin, ymax = float(yy.min()), float(yy.max())
    dx = max(xmax - xmin, 1e-9)
    dy = max(ymax - ymin, 1e-9)
    bx = np.clip(((x - xmin) / dx * n_bins).astype(int), 0, n_bins - 1)
    by = np.clip(((y - ymin) / dy * n_bins).astype(int), 0, n_bins - 1)
    return by * n_bins + bx


def summarize(df, radius=30., permutations=999, seed=20261030, ck_t=CK_T, cd8_t=CD8_T):
    needed = ['cell_id', 'x_um', 'y_um', 'CK', 'CD8', 'Autofluorescence']
    absent = [k for k in needed if k not in df.columns]
    if absent:
        raise ValueError('Missing required columns: ' + ', '.join(absent))
    if df['cell_id'].duplicated().any():
        raise ValueError('cell_id must be unique; found duplicate IDs')
    if (not np.isfinite(df[['x_um', 'y_um']].to_numpy(dtype=float)).all()):
        raise ValueError('Some cell coordinates are not finite')
    if radius <= 0 or permutations < 0:
        raise ValueError('radius must be >0 and permutations must be >=0')

    m = group_masks(df, ck_t, cd8_t)
    xy = df[['x_um', 'y_um']].to_numpy(dtype=float)
    ck_idx = np.flatnonzero(m['ck_only'])
    cd8_idx = np.flatnonzero(m['cd8_only'])
    eligible_idx = np.flatnonzero(m['ck_low_eligible'])
    counts = {k: int(v.sum()) for k, v in m.items()}

    if ck_idx.size == 0 or cd8_idx.size == 0:
        raise ValueError('Cannot compute CK/CD8 distances: one candidate group is empty')

    # All fixed CK-only cells against all CK-low eligible positions
    d = np.sqrt(np.square(xy[ck_idx, None, :] - xy[None, eligible_idx, :]).sum(axis=2))
    distance_by_global_id = {int(idx): j for j, idx in enumerate(eligible_idx)}
    obs_columns = np.array([distance_by_global_id[int(idx)] for idx in cd8_idx], dtype=int)
    obs_distances = d[:, obs_columns]
    near_obs = obs_distances <= radius
    n_neighbors = near_obs.sum(axis=1)
    nearest = obs_distances.min(axis=1)
    hit_fraction = float(np.mean(n_neighbors > 0))
    mean_count = float(n_neighbors.mean())

    rng = np.random.default_rng(seed)
    null_global = np.empty(permutations, dtype=float)
    null_local = np.empty(permutations, dtype=float)
    bins = local_bin(xy[:, 0], xy[:, 1], m['ck_low_eligible'])
    all_groups = {}
    for b in np.unique(bins[eligible_idx]):
        positions = np.flatnonzero(bins[eligible_idx] == b)
        n_marked = int(np.isin(eligible_idx[positions], cd8_idx).sum())
        all_groups[int(b)] = (positions, n_marked)

    adjacency = d <= radius
    for i in range(permutations):
        shuffled = rng.choice(eligible_idx.size, size=cd8_idx.size, replace=False)
        null_global[i] = float(np.any(adjacency[:, shuffled], axis=1).mean())
        picks = []
        for positions, n_marked in all_groups.values():
            if n_marked:
                picks.append(rng.choice(positions, size=n_marked, replace=False))
        selected = np.concatenate(picks) if picks else np.empty(0, dtype=int)
        null_local[i] = float(np.any(adjacency[:, selected], axis=1).mean())

    def null_report(arr):
        if not len(arr):
            return {'n': 0, 'mean': None, 'q025': None, 'q975': None,
                    'p_more_neighbors': None, 'p_fewer_neighbors': None}
        return {
            'n': int(len(arr)), 'mean': float(arr.mean()),
            'q025': float(np.quantile(arr, 0.025)),
            'q975': float(np.quantile(arr, 0.975)),
            'p_more_neighbors': float((1 + np.count_nonzero(arr >= hit_fraction)) / (len(arr) + 1)),
            'p_fewer_neighbors': float((1 + np.count_nonzero(arr <= hit_fraction)) / (len(arr) + 1)),
        }

    # AF is a flagging measure, not a deterministic correction/subtraction.
    af = df['Autofluorescence'].to_numpy(dtype=float)
    valid_af = af[np.isfinite(af)]
    af_p95 = float(np.quantile(valid_af, 0.95)) if valid_af.size else None
    af_summary = {'n_valid': int(valid_af.size), 'p95': af_p95}
    for k in ('ck_only', 'cd8_only', 'dual', 'low_low'):
        subset = af[m[k] & np.isfinite(af)]
        af_summary[k] = {
            'n': int(subset.size),
            'median': float(np.median(subset)) if subset.size else None,
            'n_above_global_p95': int((subset >= af_p95).sum()) if subset.size and af_p95 is not None else None
        }

    per_ck = pd.DataFrame({
        'cell_id': df.iloc[ck_idx]['cell_id'].to_numpy(),
        'x_um': xy[ck_idx, 0],
        'y_um': xy[ck_idx, 1],
        'CK': df.iloc[ck_idx]['CK'].to_numpy(dtype=float),
        'CD8': df.iloc[ck_idx]['CD8'].to_numpy(dtype=float),
        'n_CD8_only_within_radius': n_neighbors.astype(int),
        'nearest_CD8_only_um': nearest,
        'has_CD8_only_within_radius': (n_neighbors > 0).astype(int),
    })

    report = {
        'scope': 'one ROI; exploratory only; no independent biological replicates',
        'unit': 'um, centroid-to-centroid Euclidean distance',
        'radius_um': float(radius),
        'thresholds_unvalidated': {'CK': float(ck_t), 'CD8': float(cd8_t)},
        'n_cells_total': int(len(df)),
        'counts': counts,
        'observed': {
            'ck_only_with_any_cd8_only_within_radius_fraction': hit_fraction,
            'ck_only_mean_cd8_only_neighbors_within_radius': mean_count,
            'ck_only_median_nearest_cd8_only_um': float(np.median(nearest)),
        },
        'null_global_shuffle_among_ck_low_cells': null_report(null_global),
        'null_shuffle_within_spatial_4x4_bins': null_report(null_local),
        'af_audit': af_summary,
        'cautions': [
            'CK/CD8 thresholds are exploratory, not validated positivity calls.',
            'The ten double-high candidates are excluded from CK-only and CD8-only comparisons, not deleted.',
            'Random labeling among CK-low cells is only a provisional null and ignores tissue morphology.',
            'Spatially constrained shuffle helps probe coarse spatial heterogeneity but is not a validated tissue null.',
            'One ROI cannot justify general tumor biology or patient-level p-values.',
            'Check segmentation and AF/crosstalk before making immune-phenotype claims.'
        ]
    }
    return report, per_ck, null_global, null_local, m


def plot_outputs(df, output, report, gnull, lnull, masks):
    # Agg works on headless systems and does not require GUI backends.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 7))
    for label, marker in [('low_low', 'o'), ('ck_only', 's'), ('cd8_only', '^'), ('dual', 'x')]:
        sub = df.loc[masks[label]]
        ax.scatter(sub['x_um'], sub['y_um'], s=13, marker=marker, alpha=.75, label=f'{label} (n={len(sub)})')
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)')
    ax.set_title('Exploratory CK/CD8 candidate groups; not validated phenotypes')
    ax.legend(markerscale=1.5, fontsize=8)
    fig.tight_layout()
    fig.savefig(output/'spatial_candidate_groups.png', dpi=180)
    plt.close(fig)

    if len(gnull):
        fig, ax = plt.subplots(figsize=(7.5, 5))
        ax.hist(gnull, alpha=.55, bins=25, label='Global CK-low label shuffle')
        ax.hist(lnull, alpha=.55, bins=25, label='Local 4x4-bin label shuffle')
        ax.axvline(report['observed']['ck_only_with_any_cd8_only_within_radius_fraction'],
                   linestyle='--', linewidth=2, label='Observed')
        ax.set_xlabel('Fraction of CK-only cells with CD8-only candidate within radius')
        ax.set_ylabel('Permutation count')
        ax.set_title('Exploratory spatial random-labeling controls')
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(output/'spatial_permutation_null.png', dpi=180)
        plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--radius', type=float, default=30)
    p.add_argument('--permutations', type=int, default=999)
    p.add_argument('--seed', type=int, default=20261030)
    p.add_argument('--ck-threshold', type=float, default=CK_T)
    p.add_argument('--cd8-threshold', type=float, default=CD8_T)
    a = p.parse_args()
    df = pd.read_csv(a.input)
    report, per_ck, gnull, lnull, masks = summarize(
        df, a.radius, a.permutations, a.seed, a.ck_threshold, a.cd8_threshold)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out/'spatial_qc_summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    per_ck.to_csv(a.out/'ck_cell_neighborhoods.csv', index=False)
    plot_outputs(df, a.out, report, gnull, lnull, masks)
    print('=== EXPLORATORY 30-UM SPATIAL QC ===')
    print('Total cells:', report['n_cells_total'])
    for key, value in report['counts'].items():
        print(key, value)
    print('Observed:', report['observed'])
    print('Global null:', report['null_global_shuffle_among_ck_low_cells'])
    print('Local null:', report['null_shuffle_within_spatial_4x4_bins'])
    print('Autofluorescence audit:', report['af_audit'])
    print('Outputs:', a.out.resolve())
    print('CAUTION: Single-ROI exploratory results, NOT validated cell phenotypes or biological inference.')


if __name__ == '__main__':
    main()
