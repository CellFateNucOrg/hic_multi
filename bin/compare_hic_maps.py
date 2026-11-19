#!/usr/bin/env python3
"""
Quantitative and visual comparison of Non-Multi (standard nf-core/hic HiC-Pro)
vs. Multi-Contact (hicpro-multicontact-nf iterative 5'+3' rescue + all-pairs AB/AC/BC)
Hi-C contact maps (.mcool).
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
    """Resolve chromosome name across UCSC ('chrI') and Ensembl ('I') naming conventions."""
    if base_chrom in clr.chromnames:
        return base_chrom
    alt = f"chr{base_chrom}" if not base_chrom.startswith("chr") else base_chrom[3:]
    if alt in clr.chromnames:
        return alt
    return clr.chromnames[0]


def compute_ps_curve(clr: cooler.Cooler, chroms=None):
    """Compute mean contact frequency P(s) as a function of genomic distance s (in bins)."""
    binsize = clr.binsize
    if chroms is None:
        chroms = [c for c in clr.chromnames if c not in ("MtDNA", "chrM", "M")]

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
    """Compare non-zero pixel density, median counts per non-zero pixel, and correlation on a chromosome."""
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

    return {
        "binsize": int(clr_multi.binsize),
        "chrom": chrom,
        "total_upper_triangle_pixels": int(tot_pixels),
        "nonmulti_nonzero_pixels": int(nz1),
        "nonmulti_nonzero_pct": round(100.0 * nz1 / tot_pixels, 2),
        "multi_nonzero_pixels": int(nz2),
        "multi_nonzero_pct": round(100.0 * nz2 / tot_pixels, 2),
        "nonzero_pixel_gain_factor": round(nz2 / max(1, nz1), 2),
        "nonmulti_total_cis_contacts": int(np.sum(m1_raw)),
        "multi_total_cis_contacts": int(np.sum(m2_raw)),
        "contact_gain_factor": round(float(np.sum(m2_raw)) / max(1.0, float(np.sum(m1_raw))), 2),
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
    parser = argparse.ArgumentParser(description="Compare Non-Multi vs Multi-Contact Hi-C .mcool files.")
    parser.add_argument("--nonmulti-mcool", required=True, type=Path)
    parser.add_argument("--multi-mcool", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # 1. Quantitative metrics across shared resolutions (2kb, 4kb, 10kb, 20kb, 200kb)
    shared_resolutions = [2000, 4000, 10000, 20000, 200000]
    metrics = {"resolutions": {}}

    for res in shared_resolutions:
        c_nm = cooler.Cooler(f"{args.nonmulti_mcool}::/resolutions/{res}")
        c_mc = cooler.Cooler(f"{args.multi_mcool}::/resolutions/{res}")
        metrics["resolutions"][str(res)] = matrix_density_stats(c_nm, c_mc, "I")

    with open(args.outdir / "hic_map_comparison_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    # 2. Generate Visual Comparison Figure (4 panels)
    fig = plt.figure(figsize=(18, 15), dpi=180)
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.25)

    # Panel A: Chromosome I Whole-Chromosome Split-Diagonal Map at 20 kb (Raw Contacts)
    c_nm_20k = cooler.Cooler(f"{args.nonmulti_mcool}::/resolutions/20000")
    c_mc_20k = cooler.Cooler(f"{args.multi_mcool}::/resolutions/20000")
    ch_nm = resolve_chrom(c_nm_20k, "I")
    ch_mc = resolve_chrom(c_mc_20k, "I")
    m_nm_I_raw = c_nm_20k.matrix(balance=False).fetch(ch_nm).astype(np.float64)
    m_mc_I_raw = c_mc_20k.matrix(balance=False).fetch(ch_mc).astype(np.float64)

    ax1 = fig.add_subplot(gs[0, 0])
    split_raw = make_split_triangle_matrix(m_mc_I_raw, m_nm_I_raw)
    vmax_raw = np.nanpercentile(m_mc_I_raw[m_mc_I_raw > 0], 99) if np.any(m_mc_I_raw > 0) else 100
    im1 = ax1.imshow(
        split_raw + 0.5,
        cmap="YlOrRd",
        norm=LogNorm(vmin=0.5, vmax=max(10.0, vmax_raw)),
        extent=[0, c_mc_20k.chromsizes[ch_mc] / 1e6, c_mc_20k.chromsizes[ch_mc] / 1e6, 0],
    )
    ax1.plot([0, c_mc_20k.chromsizes[ch_mc] / 1e6], [0, c_mc_20k.chromsizes[ch_mc] / 1e6], color="black", lw=0.8)
    ax1.set_title("A. Chr I Raw Contact Counts (20 kb bins)\nUpper-Right: Multi-Contact (4.45M) | Lower-Left: Non-Multi (229k)", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Chromosome I Position (Mb)")
    ax1.set_ylabel("Chromosome I Position (Mb)")
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label="Raw Contacts / 20 kb bin pair")

    # Panel B: Chromosome I ICE-Normalized Split-Diagonal Map at 20 kb
    m_nm_I_bal = np.nan_to_num(c_nm_20k.matrix(balance=True).fetch(ch_nm), nan=0.0)
    m_mc_I_bal = np.nan_to_num(c_mc_20k.matrix(balance=True).fetch(ch_mc), nan=0.0)
    ax2 = fig.add_subplot(gs[0, 1])
    split_bal = make_split_triangle_matrix(m_mc_I_bal, m_nm_I_bal)
    nz_bal = m_mc_I_bal[m_mc_I_bal > 0]
    vmin_bal = np.percentile(nz_bal, 5) if len(nz_bal) > 0 else 1e-5
    vmax_bal = np.percentile(nz_bal, 98) if len(nz_bal) > 0 else 1e-2
    im2 = ax2.imshow(
        np.clip(split_bal, vmin_bal * 0.5, None),
        cmap="YlOrRd",
        norm=LogNorm(vmin=max(1e-6, vmin_bal), vmax=max(1e-4, vmax_bal)),
        extent=[0, c_mc_20k.chromsizes[ch_mc] / 1e6, c_mc_20k.chromsizes[ch_mc] / 1e6, 0],
    )
    ax2.plot([0, c_mc_20k.chromsizes[ch_mc] / 1e6], [0, c_mc_20k.chromsizes[ch_mc] / 1e6], color="black", lw=0.8)
    ax2.set_title("B. Chr I ICE-Balanced Contact Map (20 kb bins)\nUpper-Right: Multi-Contact | Lower-Left: Non-Multi", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Chromosome I Position (Mb)")
    ax2.set_ylabel("Chromosome I Position (Mb)")
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04, label="ICE-Balanced Contact Frequency")

    # Panel C: High-Resolution 4 kb Zoom-In on Chr I: 4.0 - 9.0 Mb (Split Diagonal)
    c_nm_4k = cooler.Cooler(f"{args.nonmulti_mcool}::/resolutions/4000")
    c_mc_4k = cooler.Cooler(f"{args.multi_mcool}::/resolutions/4000")
    z_nm_raw = c_nm_4k.matrix(balance=False).fetch(f"{resolve_chrom(c_nm_4k, 'I')}:4,000,000-9,000,000").astype(np.float64)
    z_mc_raw = c_mc_4k.matrix(balance=False).fetch(f"{resolve_chrom(c_mc_4k, 'I')}:4,000,000-9,000,000").astype(np.float64)
    ax3 = fig.add_subplot(gs[1, 0])
    split_zoom = make_split_triangle_matrix(z_mc_raw, z_nm_raw)
    vmax_zoom = np.nanpercentile(z_mc_raw[z_mc_raw > 0], 98) if np.any(z_mc_raw > 0) else 20
    im3 = ax3.imshow(
        split_zoom + 0.2,
        cmap="YlOrRd",
        norm=LogNorm(vmin=0.2, vmax=max(5.0, vmax_zoom)),
        extent=[4.0, 9.0, 9.0, 4.0],
    )
    ax3.plot([4.0, 9.0], [4.0, 9.0], color="black", lw=0.8)
    r4k_stats = metrics["resolutions"]["4000"]
    ax3.set_title(
        f"C. High-Resolution Zoom-In: Chr I:4.0-9.0 Mb (4 kb bins)\n"
        f"Upper-Right: Multi ({r4k_stats['multi_nonzero_pct']}% filled) | Lower-Left: Non-Multi ({r4k_stats['nonmulti_nonzero_pct']}% filled)",
        fontsize=11,
        fontweight="bold",
    )
    ax3.set_xlabel("Chromosome I Position (Mb)")
    ax3.set_ylabel("Chromosome I Position (Mb)")
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04, label="Raw Contacts / 4 kb bin pair")

    # Panel D: Contact Decay Curve P(s) at 2 kb Resolution (Raw Mean Contacts & Relative Probability)
    c_nm_2k = cooler.Cooler(f"{args.nonmulti_mcool}::/resolutions/2000")
    c_mc_2k = cooler.Cooler(f"{args.multi_mcool}::/resolutions/2000")
    dist_nm, raw_nm, prob_nm = compute_ps_curve(c_nm_2k)
    dist_mc, raw_mc, prob_mc = compute_ps_curve(c_mc_2k)

    ax4 = fig.add_subplot(gs[1, 1])
    mask_nm = (dist_nm >= 2000) & (raw_nm > 0)
    mask_mc = (dist_mc >= 2000) & (raw_mc > 0)
    ax4.loglog(dist_mc[mask_mc] / 1e3, raw_mc[mask_mc], label="Multi-Contact (Raw mean contacts/bin)", color="#d95f02", lw=2.2)
    ax4.loglog(dist_nm[mask_nm] / 1e3, raw_nm[mask_nm], label="Non-Multi nf-core/hic (Raw mean contacts/bin)", color="#1b9e77", lw=2.2)

    # Inset for normalized P(s) shape overlay
    ax4_in = ax4.inset_axes([0.48, 0.52, 0.48, 0.42])
    ax4_in.loglog(dist_mc[mask_mc] / 1e3, prob_mc[mask_mc], color="#d95f02", lw=1.8, label="Multi P(s)")
    ax4_in.loglog(dist_nm[mask_nm] / 1e3, prob_nm[mask_nm], color="#1b9e77", lw=1.5, ls="--", label="Non-Multi P(s)")
    ax4_in.set_title("Normalized P(s) Shape Overlay", fontsize=8.5, fontweight="bold")
    ax4_in.set_xlabel("Distance s (kb)", fontsize=8)
    ax4_in.set_ylabel("P(s)", fontsize=8)
    ax4_in.tick_params(labelsize=7)
    ax4_in.grid(True, which="both", alpha=0.25)
    ax4_in.legend(fontsize=7, loc="upper right")

    ax4.set_title("D. Contact Distance Decay P(s) Across C. elegans Autosomes & X (2 kb bins)", fontsize=11, fontweight="bold")
    ax4.set_xlabel("Genomic Separation s (kb)")
    ax4.set_ylabel("Mean Cis Contacts per 2 kb Bin Pair")
    ax4.grid(True, which="both", alpha=0.3)
    ax4.legend(loc="lower left", fontsize=9.5)

    fig.suptitle(
        "C. elegans Arima Hi-C+ (HN2_2 Shallow, 5.82M Read Pairs): Multi-Contact (4.45M Contacts) vs. Non-Multi nf-core/hic (229k Contacts)",
        fontsize=13.5,
        fontweight="bold",
        y=0.98,
    )
    out_png = args.outdir / "HN2_2_nonmulti_vs_multicontact_comparison.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved comparison figure to {out_png}")


if __name__ == "__main__":
    main()
