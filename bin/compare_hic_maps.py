#!/usr/bin/env python3
"""
Quantitative and visual comparison of:
  1. Fair-Trimmed Non-Multi nf-core/hic v2.1.0 (first 5 nt and last 40 nt removed -> 106 nt; 3.20M contacts)
  2. Multi-Contact hic_multi (first 5 nt removed -> 146 nt, iterative 5'+3' rescue + AB/AC/BC all-pairs; 4.45M contacts)
  3. (Optional reference) Untrimmed Non-Multi nf-core/hic v2.1.0 (151 nt; 229k contacts)
Supports any chromosome (--chrom I, II, III, IV, V, X).
"""

import argparse
import json
from pathlib import Path

import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.stats import pearsonr, spearmanr


def resolve_chrom(clr: cooler.Cooler, base_chrom: str) -> str:
    """Resolve chromosome name across UCSC ('chrX') and Ensembl ('X') naming conventions."""
    if base_chrom in clr.chromnames:
        return base_chrom
    alt = f"chr{base_chrom}" if not base_chrom.startswith("chr") else base_chrom[3:]
    if alt in clr.chromnames:
        return alt
    return clr.chromnames[0]


def open_cooler_res(mcool_path: Path, res: int) -> cooler.Cooler:
    """Open resolution from .mcool if present, otherwise fallback to sibling *.<res>_balanced.cool."""
    uri = f"{mcool_path}::/resolutions/{res}"
    try:
        return cooler.Cooler(uri)
    except Exception:
        candidates = list(mcool_path.parent.glob(f"*.{res}_balanced.cool")) + list(
            mcool_path.parent.glob(f"*.{res}.cool")
        )
        if candidates:
            return cooler.Cooler(str(candidates[0]))
        raise


def compute_ps_curve(clr: cooler.Cooler, chroms=None):
    """Compute mean contact frequency P(s) as a function of genomic distance s (in bins)."""
    binsize = clr.binsize
    if chroms is None:
        chroms = [c for c in clr.chromnames if c not in ("MtDNA", "chrM", "M")]
    else:
        chroms = [resolve_chrom(clr, c) for c in chroms]

    max_bins = max(int(clr.chromsizes[c] // binsize) + 2 for c in chroms)
    sum_counts = np.zeros(max_bins, dtype=np.float64)
    n_possible = np.zeros(max_bins, dtype=np.float64)

    for chrom in chroms:
        mat = clr.matrix(balance=False).fetch(chrom)
        n = mat.shape[0]
        for d in range(1, n):
            diag = np.diagonal(mat, offset=d)
            sum_counts[d] += np.nansum(diag)
            n_possible[d] += len(diag)

    valid = n_possible > 0
    distances_bp = np.arange(max_bins)[valid] * binsize
    mean_raw = sum_counts[valid] / n_possible[valid]
    total_cis = np.sum(sum_counts[valid])
    norm_prob = mean_raw / total_cis if total_cis > 0 else mean_raw
    return distances_bp, mean_raw, norm_prob


def matrix_density_stats(clr_nonmulti: cooler.Cooler, clr_multi: cooler.Cooler, base_chrom: str):
    """Compare non-zero pixel density, total cis contacts, and correlation on a chromosome."""
    c1 = resolve_chrom(clr_nonmulti, base_chrom)
    c2 = resolve_chrom(clr_multi, base_chrom)
    m1_raw = clr_nonmulti.matrix(balance=False).fetch(c1).astype(np.float64)
    m2_raw = clr_multi.matrix(balance=False).fetch(c2).astype(np.float64)

    triu_idx = np.triu_indices_from(m1_raw, k=1)
    v1 = m1_raw[triu_idx]
    v2 = m2_raw[triu_idx]

    nz1 = np.sum(v1 > 0)
    nz2 = np.sum(v2 > 0)
    tot_pixels = len(v1)

    r_pearson, _ = pearsonr(np.log1p(v1), np.log1p(v2))
    r_spearman, _ = spearmanr(v1, v2)

    # Distance-stratified Pearson correlation (diagonals 1..100 bins)
    strat_corrs = []
    strat_weights = []
    for d in range(1, min(100, m1_raw.shape[0] - 1)):
        d1 = np.diagonal(m1_raw, offset=d)
        d2 = np.diagonal(m2_raw, offset=d)
        if np.std(d1) > 0 and np.std(d2) > 0:
            c, _ = pearsonr(d1, d2)
            if not np.isnan(c):
                strat_corrs.append(c)
                strat_weights.append(len(d1))
    scc_approx = float(np.average(strat_corrs, weights=strat_weights)) if strat_corrs else 0.0

    tot_cis_1 = int(np.sum(np.triu(m1_raw, k=0)))
    tot_cis_2 = int(np.sum(np.triu(m2_raw, k=0)))

    return {
        "binsize": int(clr_multi.binsize),
        "chrom": base_chrom,
        "total_upper_triangle_pixels": int(tot_pixels),
        "nonmulti_nonzero_pixels": int(nz1),
        "nonmulti_nonzero_pct": round(100.0 * nz1 / tot_pixels, 2),
        "multi_nonzero_pixels": int(nz2),
        "multi_nonzero_pct": round(100.0 * nz2 / tot_pixels, 2),
        "nonzero_pixel_gain_factor": round(nz2 / max(1, nz1), 2),
        "nonmulti_total_cis_contacts": tot_cis_1,
        "multi_total_cis_contacts": tot_cis_2,
        "contact_gain_factor": round(float(tot_cis_2) / max(1.0, float(tot_cis_1)), 3),
        "log1p_pearson_r": round(float(r_pearson), 4),
        "spearman_r": round(float(r_spearman), 4),
        "distance_stratified_corr_1_to_100bins": round(scc_approx, 4),
    }


def make_split_triangle_matrix(m_upper: np.ndarray, m_lower: np.ndarray) -> np.ndarray:
    """Combine m_upper in the upper triangle and m_lower in the lower triangle."""
    combined = np.zeros_like(m_upper, dtype=np.float64)
    iu = np.triu_indices_from(m_upper, k=0)
    il = np.tril_indices_from(m_lower, k=-1)
    combined[iu] = m_upper[iu]
    combined[il] = m_lower[il]
    return combined


def main():
    parser = argparse.ArgumentParser(description="Compare Trimmed Non-Multi vs Multi-Contact Hi-C .mcool files.")
    parser.add_argument("--nonmulti-mcool", required=True, type=Path, help="Trimmed nf-core/hic .mcool (5'-5nt, 3'-40nt)")
    parser.add_argument("--multi-mcool", required=True, type=Path, help="Multi-contact hic_multi .mcool")
    parser.add_argument("--untrimmed-mcool", type=Path, default=None, help="Optional untrimmed nf-core/hic .mcool")
    parser.add_argument("--chrom", default="I", help="Target chromosome to plot (default: I, e.g. X)")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    target_chrom = args.chrom.replace("chr", "")

    # 1. Quantitative metrics across shared resolutions (1kb, 2kb, 4kb, 10kb, 20kb, 200kb)
    shared_resolutions = [1000, 2000, 4000, 10000, 20000, 200000]
    metrics = {"chrom": target_chrom, "resolutions": {}}

    for res in shared_resolutions:
        try:
            c_nm = open_cooler_res(args.nonmulti_mcool, res)
            c_mc = open_cooler_res(args.multi_mcool, res)
            metrics["resolutions"][str(res)] = matrix_density_stats(c_nm, c_mc, target_chrom)
        except Exception as e:
            print(f"Skipping resolution {res}: {e}")

    json_path = args.outdir / f"hic_map_comparison_metrics_chr{target_chrom}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    # 2. Generate Visual Comparison Figure (4 panels)
    fig = plt.figure(figsize=(18, 15.5), dpi=180)
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.25)

    # Panel A: Whole-Chromosome Split-Diagonal Map at 20 kb (Raw Contacts)
    c_nm_20k = open_cooler_res(args.nonmulti_mcool, 20000)
    c_mc_20k = open_cooler_res(args.multi_mcool, 20000)
    ch_nm = resolve_chrom(c_nm_20k, target_chrom)
    ch_mc = resolve_chrom(c_mc_20k, target_chrom)
    m_nm_raw = c_nm_20k.matrix(balance=False).fetch(ch_nm).astype(np.float64)
    m_mc_raw = c_mc_20k.matrix(balance=False).fetch(ch_mc).astype(np.float64)
    chr_len_mb = c_mc_20k.chromsizes[ch_mc] / 1e6

    ax1 = fig.add_subplot(gs[0, 0])
    split_raw = make_split_triangle_matrix(m_mc_raw, m_nm_raw)
    vmax_raw = np.nanpercentile(m_mc_raw[m_mc_raw > 0], 99) if np.any(m_mc_raw > 0) else 100
    im1 = ax1.imshow(
        split_raw + 0.5,
        cmap="YlOrRd",
        norm=LogNorm(vmin=0.5, vmax=max(10.0, vmax_raw)),
        extent=[0, chr_len_mb, chr_len_mb, 0],
    )
    ax1.plot([0, chr_len_mb], [0, chr_len_mb], color="black", lw=0.8)
    r20k = metrics["resolutions"]["20000"]
    ax1.set_title(
        f"A. Chr {target_chrom} Raw Contact Counts (20 kb bins, r={r20k['log1p_pearson_r']})\n"
        f"Upper-Right: Multi ({r20k['multi_total_cis_contacts']:,} cis) | Lower-Left: Trimmed nf-core ({r20k['nonmulti_total_cis_contacts']:,} cis)",
        fontsize=10.5,
        fontweight="bold",
    )
    ax1.set_xlabel(f"Chromosome {target_chrom} Position (Mb)")
    ax1.set_ylabel(f"Chromosome {target_chrom} Position (Mb)")
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label="Raw Contacts / 20 kb bin pair")

    # Panel B: Whole-Chromosome ICE-Normalized Split-Diagonal Map at 20 kb
    m_nm_bal = np.nan_to_num(c_nm_20k.matrix(balance=True).fetch(ch_nm), nan=0.0)
    m_mc_bal = np.nan_to_num(c_mc_20k.matrix(balance=True).fetch(ch_mc), nan=0.0)
    ax2 = fig.add_subplot(gs[0, 1])
    split_bal = make_split_triangle_matrix(m_mc_bal, m_nm_bal)
    nz_bal = m_mc_bal[m_mc_bal > 0]
    vmin_bal = np.percentile(nz_bal, 5) if len(nz_bal) > 0 else 1e-5
    vmax_bal = np.percentile(nz_bal, 98.5) if len(nz_bal) > 0 else 1e-2
    im2 = ax2.imshow(
        np.clip(split_bal, vmin_bal * 0.5, None),
        cmap="YlOrRd",
        norm=LogNorm(vmin=max(1e-6, vmin_bal), vmax=max(1e-4, vmax_bal)),
        extent=[0, chr_len_mb, chr_len_mb, 0],
    )
    ax2.plot([0, chr_len_mb], [0, chr_len_mb], color="black", lw=0.8)
    ax2.set_title(
        f"B. Chr {target_chrom} ICE-Balanced Contact Map (20 kb bins)\n"
        f"Upper-Right: Multi ({r20k['multi_nonzero_pct']}% filled) | Lower-Left: Trimmed nf-core ({r20k['nonmulti_nonzero_pct']}% filled)",
        fontsize=10.5,
        fontweight="bold",
    )
    ax2.set_xlabel(f"Chromosome {target_chrom} Position (Mb)")
    ax2.set_ylabel(f"Chromosome {target_chrom} Position (Mb)")
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04, label="ICE-Balanced Contact Frequency")

    # Panel C: High-Resolution 4 kb Zoom-In on Central TAD Region (6.0 - 12.0 Mb for X, 5.0 - 8.0 Mb otherwise)
    if target_chrom == "X":
        z_start_mb, z_end_mb = 6.0, 12.0
        zoom_res = 4000
    else:
        z_start_mb, z_end_mb = 5.0, 8.0
        zoom_res = 2000

    c_nm_z = open_cooler_res(args.nonmulti_mcool, zoom_res)
    c_mc_z = open_cooler_res(args.multi_mcool, zoom_res)
    z_start_bp = int(z_start_mb * 1e6)
    z_end_bp = int(z_end_mb * 1e6)
    z_nm_bal = np.nan_to_num(
        c_nm_z.matrix(balance=True).fetch(f"{resolve_chrom(c_nm_z, target_chrom)}:{z_start_bp:,}-{z_end_bp:,}"),
        nan=0.0,
    )
    z_mc_bal = np.nan_to_num(
        c_mc_z.matrix(balance=True).fetch(f"{resolve_chrom(c_mc_z, target_chrom)}:{z_start_bp:,}-{z_end_bp:,}"),
        nan=0.0,
    )
    ax3 = fig.add_subplot(gs[1, 0])
    split_zoom = make_split_triangle_matrix(z_mc_bal, z_nm_bal)
    nz_z = z_mc_bal[z_mc_bal > 0]
    vmin_z = np.percentile(nz_z, 8) if len(nz_z) > 0 else 1e-5
    vmax_z = np.percentile(nz_z, 98.5) if len(nz_z) > 0 else 1e-2
    im3 = ax3.imshow(
        np.clip(split_zoom, vmin_z * 0.5, None),
        cmap="YlOrRd",
        norm=LogNorm(vmin=max(1e-6, vmin_z), vmax=max(1e-4, vmax_z)),
        extent=[z_start_mb, z_end_mb, z_end_mb, z_start_mb],
    )
    ax3.plot([z_start_mb, z_end_mb], [z_start_mb, z_end_mb], color="black", lw=0.8)
    rz_stats = metrics["resolutions"][str(zoom_res)]
    ax3.set_title(
        f"C. High-Resolution ICE-Balanced Zoom-In: Chr {target_chrom}:{z_start_mb:.1f}-{z_end_mb:.1f} Mb ({zoom_res//1000} kb bins)\n"
        f"Upper-Right: Multi (+{round((rz_stats['contact_gain_factor']-1)*100, 1)}% cis contacts) | Lower-Left: Trimmed nf-core",
        fontsize=10.5,
        fontweight="bold",
    )
    ax3.set_xlabel(f"Chromosome {target_chrom} Position (Mb)")
    ax3.set_ylabel(f"Chromosome {target_chrom} Position (Mb)")
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04, label=f"ICE-Balanced Contact Frequency ({zoom_res//1000} kb)")

    # Panel D: Contact Distance Decay Curve P(s) on target_chrom at 2 kb Resolution
    c_nm_2k = open_cooler_res(args.nonmulti_mcool, 2000)
    c_mc_2k = open_cooler_res(args.multi_mcool, 2000)
    dist_nm, raw_nm, prob_nm = compute_ps_curve(c_nm_2k, chroms=[target_chrom])
    dist_mc, raw_mc, prob_mc = compute_ps_curve(c_mc_2k, chroms=[target_chrom])

    ax4 = fig.add_subplot(gs[1, 1])
    mask_nm = (dist_nm >= 2000) & (raw_nm > 0)
    mask_mc = (dist_mc >= 2000) & (raw_mc > 0)
    ax4.loglog(
        dist_mc[mask_mc] / 1e3,
        raw_mc[mask_mc],
        label=f"Multi-Contact Chr {target_chrom} ({r20k['multi_total_cis_contacts']:,} cis contacts)",
        color="#d95f02",
        lw=2.2,
    )
    ax4.loglog(
        dist_nm[mask_nm] / 1e3,
        raw_nm[mask_nm],
        label=f"Trimmed nf-core/hic Chr {target_chrom} ({r20k['nonmulti_total_cis_contacts']:,} cis contacts)",
        color="#7570b3",
        lw=2.0,
    )

    if args.untrimmed_mcool and args.untrimmed_mcool.exists():
        c_un_2k = open_cooler_res(args.untrimmed_mcool, 2000)
        dist_un, raw_un, prob_un = compute_ps_curve(c_un_2k, chroms=[target_chrom])
        mask_un = (dist_un >= 2000) & (raw_un > 0)
        ax4.loglog(
            dist_un[mask_un] / 1e3,
            raw_un[mask_un],
            label=f"Untrimmed nf-core/hic Chr {target_chrom}",
            color="#1b9e77",
            lw=1.5,
            alpha=0.8,
        )

    # Inset for normalized P(s) shape overlay
    ax4_in = ax4.inset_axes([0.48, 0.52, 0.48, 0.42])
    ax4_in.loglog(dist_mc[mask_mc] / 1e3, prob_mc[mask_mc], color="#d95f02", lw=1.8, label="Multi P(s)")
    ax4_in.loglog(dist_nm[mask_nm] / 1e3, prob_nm[mask_nm], color="#7570b3", lw=1.5, ls="--", label="Trimmed nf-core P(s)")
    ax4_in.set_title(f"Chr {target_chrom} Normalized P(s) Overlay", fontsize=8.5, fontweight="bold")
    ax4_in.set_xlabel("Distance s (kb)", fontsize=8)
    ax4_in.set_ylabel("P(s)", fontsize=8)
    ax4_in.tick_params(labelsize=7)
    ax4_in.grid(True, which="both", alpha=0.25)
    ax4_in.legend(fontsize=7, loc="upper right")

    ax4.set_title(f"D. Chr {target_chrom} Contact Distance Decay P(s) (2 kb bins)", fontsize=10.5, fontweight="bold")
    ax4.set_xlabel("Genomic Separation s (kb)")
    ax4.set_ylabel("Mean Cis Contacts per 2 kb Bin Pair")
    ax4.grid(True, which="both", alpha=0.3)
    ax4.legend(loc="lower left", fontsize=8.5)

    fig.suptitle(
        f"C. elegans Chromosome {target_chrom} Arima Hi-C+ (HN2_2): Multi-Contact ({r20k['multi_total_cis_contacts']:,} cis) vs. Trimmed nf-core/hic ({r20k['nonmulti_total_cis_contacts']:,} cis)",
        fontsize=13,
        fontweight="bold",
        y=0.98,
    )
    out_png = args.outdir / f"HN2_2_chr{target_chrom}_trimmed_nfcore_vs_multicontact.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved comparison figure to {out_png}")


if __name__ == "__main__":
    main()
