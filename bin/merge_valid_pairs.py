#!/usr/bin/env python3
"""
Merge chunked .validPairs files for a sample, remove PCR duplicates,
and generate HiC-Pro compatible .allValidPairs, .mergestat, and multi-contact QC report.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge and deduplicate HiC-Pro .validPairs chunks and aggregate QC metrics."
    )
    parser.add_argument("-s", "--sample", required=True, help="Sample ID")
    parser.add_argument(
        "-i",
        "--valid-pairs",
        nargs="+",
        required=True,
        type=Path,
        help="Input .validPairs chunk files",
    )
    parser.add_argument(
        "--map-stats",
        nargs="*",
        default=[],
        type=Path,
        help="Per-chunk iterative mapping stats JSON files",
    )
    parser.add_argument(
        "--pair-stats",
        nargs="*",
        default=[],
        type=Path,
        help="Per-chunk multi-contact pairing stats JSON files",
    )
    parser.add_argument(
        "-o",
        "--output-all-valid-pairs",
        required=True,
        type=Path,
        help="Output deduplicated .allValidPairs file",
    )
    parser.add_argument(
        "--output-mergestat",
        required=True,
        type=Path,
        help="Output HiC-Pro compatible .mergestat file",
    )
    parser.add_argument(
        "--output-summary-json",
        required=True,
        type=Path,
        help="Output comprehensive JSON summary report",
    )
    parser.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Do not remove PCR duplicate pairs",
    )
    parser.add_argument(
        "-p", "--threads", type=int, default=4, help="Threads for GNU sort (default: 4)"
    )
    parser.add_argument(
        "-m", "--sort-mem", default="4G", help="Memory limit for GNU sort (default: 4G)"
    )
    return parser.parse_args()


def aggregate_jsons(paths: List[Path]) -> Dict[str, int]:
    totals: Dict[str, int] = {}
    for p in paths:
        if not p.exists():
            continue
        with open(p, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        for k, v in data.items():
            if isinstance(v, int):
                totals[k] = totals.get(k, 0) + v
            elif k == "rounds" and isinstance(v, list):
                for r_obj in v:
                    r_num = r_obj.get("round", 1)
                    for rk, rv in r_obj.items():
                        if rk != "round" and isinstance(rv, int):
                            key = f"round{r_num}_{rk}"
                            totals[key] = totals.get(key, 0) + rv
    return totals


def main():
    args = parse_args()
    args.output_all_valid_pairs.parent.mkdir(parents=True, exist_ok=True)

    # Sort all validPairs chunks by (chr1, pos1, chr2, pos2, strand1, strand2) and deduplicate
    total_before_dedup = 0
    total_after_dedup = 0
    cis_after_dedup = 0
    cis_short_lt_1kb = 0
    cis_1kb_to_20kb = 0
    cis_gt_20kb = 0
    trans_after_dedup = 0
    intra_mate_after_dedup = 0
    from_3frag_after_dedup = 0
    from_4plus_after_dedup = 0

    with tempfile.TemporaryDirectory(prefix="hicpro_merge_") as tmpdir_str:
        sorted_tmp = Path(tmpdir_str) / "sorted_pairs.txt"

        # Use external GNU sort if available for massive multi-million pair files, fallback to Python sort
        try:
            sort_cmd = [
                "sort",
                f"--parallel={args.threads}",
                f"-S{args.sort_mem}",
                "-T",
                tmpdir_str,
                "-k2,2V",
                "-k3,3n",
                "-k5,5V",
                "-k6,6n",
                "-k4,4",
                "-k7,7",
                "-o",
                str(sorted_tmp),
            ] + [str(p) for p in args.valid_pairs]
            res = subprocess.run(sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0:
                raise RuntimeError(res.stderr)
        except Exception:
            # Pure-Python fallback
            rows = []
            for p in args.valid_pairs:
                with open(p, "r", encoding="utf-8") as fin:
                    for line in fin:
                        if not line.strip():
                            continue
                        cols = line.rstrip("\r\n").split("\t")
                        rows.append((cols[1], int(cols[2]), cols[4], int(cols[5]), cols[3], cols[6], line))
            rows.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4], r[5]))
            with open(sorted_tmp, "w", encoding="utf-8") as fout:
                for r in rows:
                    fout.write(r[6])

        prev_coord_key = None
        with open(sorted_tmp, "r", encoding="utf-8") as fin, open(
            args.output_all_valid_pairs, "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                if not line.strip():
                    continue
                total_before_dedup += 1
                cols = line.rstrip("\r\n").split("\t")
                # cols: 0=id, 1=chr1, 2=pos1, 3=strand1, 4=chr2, 5=pos2, 6=strand2, 7=dist, ...
                coord_key = (cols[1], cols[2], cols[3], cols[4], cols[5], cols[6])
                if not args.keep_duplicates and coord_key == prev_coord_key:
                    continue
                prev_coord_key = coord_key
                fout.write(line)
                total_after_dedup += 1

                pair_id = cols[0]
                if "_K3_" in pair_id:
                    from_3frag_after_dedup += 1
                elif "_K" in pair_id and not ("_K2_" in pair_id or "_K3_" in pair_id):
                    from_4plus_after_dedup += 1
                if pair_id.endswith("_M11") or pair_id.endswith("_M22"):
                    intra_mate_after_dedup += 1

                if cols[1] == cols[4]:
                    cis_after_dedup += 1
                    dist = abs(int(cols[5]) - int(cols[2]))
                    if dist < 1000:
                        cis_short_lt_1kb += 1
                    elif dist <= 20000:
                        cis_1kb_to_20kb += 1
                    else:
                        cis_gt_20kb += 1
                else:
                    trans_after_dedup += 1

    map_totals = aggregate_jsons(args.map_stats)
    pair_totals = aggregate_jsons(args.pair_stats)

    # Write HiC-Pro compatible .mergestat
    args.output_mergestat.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_mergestat, "w", encoding="utf-8") as ms_out:
        ms_out.write(f"valid_interaction\t{total_before_dedup}\n")
        ms_out.write(f"valid_interaction_rmdup\t{total_after_dedup}\n")
        ms_out.write(f"trans_interaction\t{trans_after_dedup}\n")
        ms_out.write(f"cis_interaction\t{cis_after_dedup}\n")
        ms_out.write(f"cis_shortRange\t{cis_short_lt_1kb + cis_1kb_to_20kb}\n")
        ms_out.write(f"cis_longRange\t{cis_gt_20kb}\n")
        ms_out.write(f"dangling_end_pairs\t{pair_totals.get('dangling_end_pairs', 0)}\n")
        ms_out.write(f"self_circle_pairs\t{pair_totals.get('self_circle_pairs', 0)}\n")
        ms_out.write(f"religation_adjacent_pairs\t{pair_totals.get('religation_adjacent_pairs', 0)}\n")
        ms_out.write(f"valid_pairs_from_3frag_molecules_rmdup\t{from_3frag_after_dedup}\n")
        ms_out.write(f"valid_pairs_from_4plus_frag_molecules_rmdup\t{from_4plus_after_dedup}\n")
        ms_out.write(f"valid_pairs_intra_mate_rescued_rmdup\t{intra_mate_after_dedup}\n")

    summary = {
        "sample": args.sample,
        "valid_pairs_before_dedup": total_before_dedup,
        "valid_pairs_after_dedup": total_after_dedup,
        "pcr_duplicate_rate": round(
            (total_before_dedup - total_after_dedup) / total_before_dedup, 4
        )
        if total_before_dedup > 0
        else 0.0,
        "cis_after_dedup": cis_after_dedup,
        "cis_short_lt_1kb": cis_short_lt_1kb,
        "cis_1kb_to_20kb": cis_1kb_to_20kb,
        "cis_gt_20kb": cis_gt_20kb,
        "trans_after_dedup": trans_after_dedup,
        "valid_pairs_from_3frag_molecules_rmdup": from_3frag_after_dedup,
        "valid_pairs_from_4plus_frag_molecules_rmdup": from_4plus_after_dedup,
        "valid_pairs_intra_mate_rescued_rmdup": intra_mate_after_dedup,
        "mapping_totals": map_totals,
        "pairing_totals": pair_totals,
    }

    args.output_summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_summary_json, "w", encoding="utf-8") as sj_out:
        json.dump(summary, sj_out, indent=2)

    sys.stderr.write(
        f"[merge_valid_pairs] Sample {args.sample}: {total_before_dedup} valid pairs -> "
        f"{total_after_dedup} unique valid pairs (3-frag: {from_3frag_after_dedup}, "
        f"intra-mate: {intra_mate_after_dedup}).\n"
    )


if __name__ == "__main__":
    main()
