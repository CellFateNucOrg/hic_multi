#!/usr/bin/env python3
"""
Iterative HiC-Pro style Bowtie2 mapper with multi-fragment 5' + 3' restriction rescue.

Algorithm per mate (R1 or R2):
  1. Round 0 (Full-length mapping):
     Align full-length reads with Bowtie2 (--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end).
     Keep alignments with MAPQ >= min_mapq (default 30).
  2. Iterative Restriction Cutting & Rescue (Rounds 1 .. max_rounds):
     For every unmapped / low-MAPQ sequence (initially the unmapped full read; in subsequent
     rounds the unmapped 3' remainder from the previous round):
       a) Locate the 5'-most ligation junction motif.
       b) Split at the exact blunt-end ligation point into:
            - 5' fragment (mapped in this round if length >= min_frag_len)
            - 3' remainder fragment (also mapped in this round if length >= min_frag_len)
       c) Keep any 5' or 3' fragment that maps with MAPQ >= min_mapq.
       d) Any 3' remainder fragment that fails to map (MAPQ < min_mapq) is carried into
          the next round to be cut at the next restriction ligation junction and mapped again.
  3. Output a queryname-sorted BAM containing all mapped fragments per read with SAM tags:
       - ZM:i:<mate>      Mate number (1 or 2)
       - ZI:i:<seg_idx>   0-based 5'->3' segment index along the read
       - ZR:i:<round>     Mapping round (0 = full-length, 1+ = iterative cut round)
       - ZS:i:<start>     0-based start offset in the original read
       - ZE:i:<end>       0-based end offset in the original read
       - ZK:Z:<kind>      FULL, K5 (5' cut piece), or K3 (3' remainder piece)
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam


complement_trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(complement_trans)[::-1]


IUPAC_REGEX: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGTN]",
}


def build_junction_regexes(spec: str) -> List[Tuple[re.Pattern, int, str]]:
    """
    Build compiled regex patterns and cut offsets for ligation junctions.

    Preset 'arima' includes:
      - Biochemical fill-in junctions for DpnII (^GATC) + HinfI (G^ANTC):
          GATCGATC  (cut +4)
          GATCANTC  (cut +4)
          GANTCGATC (cut +5)
          GANTCANTC (cut +5)
      - 8-bp symmetric notation junctions:
          GANTGATC  (cut +4)
          GATCGANTC (cut +4)
          GANTANTC  (cut +4)

    Custom spec format: comma-separated 'MOTIF:OFFSET' or 'LEFT_PART^RIGHT_PART',
    e.g. 'GATCGATC:4,GANTGATC:4,GATCGANTC:4,GANTANTC:4' or 'GATC^GATC,GANT^GATC'.
    """
    if spec.strip().lower() == "arima":
        raw_items = [
            ("GATCGATC", 4),
            ("GATCANTC", 4),
            ("GANTCGATC", 5),
            ("GANTCANTC", 5),
            ("GANTGATC", 4),
            ("GATCGANTC", 4),
            ("GANTANTC", 4),
        ]
    elif spec.strip().lower() == "arima_strict_8bp":
        raw_items = [
            ("GATCGATC", 4),
            ("GANTGATC", 4),
            ("GATCGANTC", 4),
            ("GANTANTC", 4),
        ]
    else:
        raw_items = []
        for token in spec.split(","):
            token = token.strip().upper()
            if not token:
                continue
            if "^" in token:
                cut_offset = token.index("^")
                motif = token.replace("^", "")
                raw_items.append((motif, cut_offset))
            elif ":" in token:
                motif, off_str = token.split(":", 1)
                raw_items.append((motif.strip(), int(off_str.strip())))
            else:
                # Default to midpoint split
                raw_items.append((token, len(token) // 2))

    compiled = []
    for motif, offset in raw_items:
        pattern_str = "".join(IUPAC_REGEX.get(ch, re.escape(ch)) for ch in motif)
        compiled.append((re.compile(pattern_str), offset, motif))
    return compiled


def find_first_ligation_cut(
    seq: str, junctions: List[Tuple[re.Pattern, int, str]]
) -> Optional[Tuple[int, str]]:
    """
    Find the 5'-most ligation junction cut position in `seq`.
    Returns (cut_pos, motif_name) or None if no junction is found.
    """
    best_match_start: Optional[int] = None
    best_cut_pos: Optional[int] = None
    best_motif: Optional[str] = None

    seq_upper = seq.upper()
    for regex, offset, motif in junctions:
        match = regex.search(seq_upper)
        if match is not None:
            m_start = match.start()
            cut_pos = m_start + offset
            if 0 < cut_pos < len(seq):
                if best_match_start is None or m_start < best_match_start:
                    best_match_start = m_start
                    best_cut_pos = cut_pos
                    best_motif = motif

    if best_cut_pos is None:
        return None
    return best_cut_pos, best_motif


def run_bowtie2(
    bt2_index: str,
    fastq_path: Path,
    sam_out_path: Path,
    threads: int,
    seed_len: int,
    score_min: str,
    extra_args: str = "",
) -> None:
    cmd = [
        "bowtie2",
        "--very-sensitive",
        "-L",
        str(seed_len),
        "--score-min",
        score_min,
        "--end-to-end",
        "--reorder",
        "-p",
        str(threads),
        "-x",
        bt2_index,
        "-U",
        str(fastq_path),
        "-S",
        str(sam_out_path),
    ]
    if extra_args:
        cmd.extend(extra_args.split())

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"[iterative_cut_and_map] Bowtie2 failed:\n{proc.stderr}\n")
        raise RuntimeError(f"Bowtie2 exited with code {proc.returncode}")


def get_original_seq_and_qual(aln: pysam.AlignedSegment) -> Tuple[str, str]:
    """
    Extract the original 5'->3' read sequence and ASCII phred33 quality string
    from a pysam AlignedSegment (reversing if Bowtie2 reverse-complemented it).
    """
    seq = aln.query_sequence or ""
    quals = aln.query_qualities
    if quals is not None:
        qual_str = "".join(chr(q + 33) for q in quals)
    else:
        qual_str = "I" * len(seq)

    if aln.is_reverse and not aln.is_unmapped:
        seq = revcomp(seq)
        qual_str = qual_str[::-1]
    return seq, qual_str


def apply_5prime_trim(
    in_fastq: Path, out_fastq: Path, trim_5p: int
) -> None:
    """Optional 5' fixed base trimming before Round 0 if requested."""
    open_fn = gzip.open if str(in_fastq).endswith(".gz") else open
    with open_fn(in_fastq, "rt", encoding="utf-8") as fin, open(
        out_fastq, "w", encoding="utf-8"
    ) as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().rstrip("\r\n")
            plus = fin.readline()
            qual = fin.readline().rstrip("\r\n")
            if len(seq) > trim_5p:
                seq = seq[trim_5p:]
                qual = qual[trim_5p:]
            fout.write(f"{header}{seq}\n{plus}{qual}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="HiC-Pro style Bowtie2 iterative mapper rescuing both 5' and 3' restriction fragments."
    )
    parser.add_argument("-i", "--fastq", required=True, type=Path, help="Input FASTQ (.fq or .fq.gz)")
    parser.add_argument("-x", "--bt2-index", required=True, help="Bowtie2 index prefix")
    parser.add_argument("-m", "--mate", required=True, type=int, choices=[1, 2], help="Mate number (1 or 2)")
    parser.add_argument("-o", "--output-bam", required=True, type=Path, help="Output queryname-sorted BAM")
    parser.add_argument("--stats-json", required=True, type=Path, help="Output mapping stats JSON")
    parser.add_argument("-q", "--min-mapq", type=int, default=30, help="Minimum MAPQ threshold (default: 30)")
    parser.add_argument(
        "--ligation-motifs",
        default="arima",
        help="Ligation junction preset ('arima', 'arima_strict_8bp') or comma-separated MOTIF:CUT_OFFSET list",
    )
    parser.add_argument(
        "--min-frag-len",
        type=int,
        default=20,
        help="Minimum fragment length (bp) after restriction cutting to attempt mapping (default: 20)",
    )
    parser.add_argument(
        "--max-rounds",
        type=int,
        default=4,
        help="Maximum number of iterative restriction-cutting rounds after full-length mapping (default: 4)",
    )
    parser.add_argument(
        "--trim-5p",
        type=int,
        default=0,
        help="Optional number of bases to trim from 5' end of raw reads before Round 0 (default: 0)",
    )
    parser.add_argument(
        "--seed-len-full",
        type=int,
        default=30,
        help="Bowtie2 -L seed length for Round 0 full-length mapping (default: 30)",
    )
    parser.add_argument(
        "--seed-len-cut",
        type=int,
        default=20,
        help="Bowtie2 -L seed length for cut restriction fragments (default: 20)",
    )
    parser.add_argument(
        "--score-min",
        default="L,-0.6,-0.2",
        help="Bowtie2 --score-min function (default: L,-0.6,-0.2 as in HiC-Pro)",
    )
    parser.add_argument("-p", "--threads", type=int, default=4, help="Number of CPU threads (default: 4)")
    return parser.parse_args()


def main():
    args = parse_args()
    junctions = build_junction_regexes(args.ligation_motifs)

    stats = {
        "mate": args.mate,
        "min_mapq": args.min_mapq,
        "min_frag_len": args.min_frag_len,
        "max_rounds": args.max_rounds,
        "total_input_reads": 0,
        "round0_full_mapped_q30": 0,
        "round0_unmapped_or_lowq": 0,
        "rounds": [],
        "total_mapped_segments_q30": 0,
        "reads_with_0_segments": 0,
        "reads_with_1_segment": 0,
        "reads_with_2_segments": 0,
        "reads_with_3plus_segments": 0,
    }

    segments_per_read: Dict[str, int] = {}

    with tempfile.TemporaryDirectory(prefix="hicpro_iter_") as tmpdir_str:
        tmpdir = Path(tmpdir_str)

        # Optional 5' trimming
        if args.trim_5p > 0:
            round0_fastq = tmpdir / "round0_trimmed.fastq"
            apply_5prime_trim(args.fastq, round0_fastq, args.trim_5p)
        else:
            round0_fastq = args.fastq

        round0_sam = tmpdir / "round0.sam"
        run_bowtie2(
            bt2_index=args.bt2_index,
            fastq_path=round0_fastq,
            sam_out_path=round0_sam,
            threads=args.threads,
            seed_len=args.seed_len_full,
            score_min=args.score_min,
        )

        unsorted_bam_path = tmpdir / "accepted_unsorted.bam"
        # Queue element: (orig_read_name, seq, qual, offset_in_orig_read, next_seg_idx)
        unmapped_queue: List[Tuple[str, str, str, int, int]] = []

        with pysam.AlignmentFile(round0_sam, "r") as sam_in:
            header = sam_in.header.to_dict()
            with pysam.AlignmentFile(unsorted_bam_path, "wb", header=header) as bam_out:
                for aln in sam_in:
                    if aln.is_secondary or aln.is_supplementary:
                        continue
                    qname = aln.query_name.split()[0]
                    stats["total_input_reads"] += 1
                    segments_per_read.setdefault(qname, 0)

                    if not aln.is_unmapped and aln.mapping_quality >= args.min_mapq:
                        seq_len = aln.query_length or len(aln.query_sequence or "")
                        aln.query_name = qname
                        aln.set_tag("ZM", args.mate, value_type="i")
                        aln.set_tag("ZI", 0, value_type="i")
                        aln.set_tag("ZR", 0, value_type="i")
                        aln.set_tag("ZS", 0, value_type="i")
                        aln.set_tag("ZE", seq_len, value_type="i")
                        aln.set_tag("ZK", "FULL", value_type="Z")
                        bam_out.write(aln)
                        stats["round0_full_mapped_q30"] += 1
                        stats["total_mapped_segments_q30"] += 1
                        segments_per_read[qname] += 1
                    else:
                        stats["round0_unmapped_or_lowq"] += 1
                        seq, qual = get_original_seq_and_qual(aln)
                        unmapped_queue.append((qname, seq, qual, 0, 0))

                # Iterative cutting & mapping rounds (1 .. max_rounds)
                for r in range(1, args.max_rounds + 1):
                    if not unmapped_queue:
                        break

                    round_fastq = tmpdir / f"round_{r}.fastq"
                    round_sam = tmpdir / f"round_{r}.sam"

                    round_info = {
                        "round": r,
                        "input_unmapped_sequences": len(unmapped_queue),
                        "sequences_with_junction": 0,
                        "sequences_without_junction": 0,
                        "emitted_5p_candidates": 0,
                        "emitted_3p_candidates": 0,
                        "mapped_5p_q30": 0,
                        "mapped_3p_q30": 0,
                        "unmapped_3p_passed_to_next_round": 0,
                    }

                    candidates_written = 0
                    with open(round_fastq, "w", encoding="utf-8") as fq_out:
                        for qname, seq, qual, offset_start, seg_idx in unmapped_queue:
                            cut_res = find_first_ligation_cut(seq, junctions)
                            if cut_res is None:
                                round_info["sequences_without_junction"] += 1
                                continue

                            round_info["sequences_with_junction"] += 1
                            cut_pos, _motif = cut_res

                            # 5' fragment
                            seq_5p = seq[:cut_pos]
                            qual_5p = qual[:cut_pos]
                            span5_start = offset_start
                            span5_end = offset_start + cut_pos

                            # 3' remainder fragment
                            seq_3p = seq[cut_pos:]
                            qual_3p = qual[cut_pos:]
                            span3_start = offset_start + cut_pos
                            span3_end = offset_start + len(seq)

                            next_idx = seg_idx
                            if len(seq_5p) >= args.min_frag_len:
                                fq_header_5p = f"{qname}|K5|S{next_idx}|B{span5_start}|E{span5_end}"
                                fq_out.write(f"@{fq_header_5p}\n{seq_5p}\n+\n{qual_5p}\n")
                                round_info["emitted_5p_candidates"] += 1
                                candidates_written += 1
                                next_idx += 1

                            if len(seq_3p) >= args.min_frag_len:
                                fq_header_3p = f"{qname}|K3|S{next_idx}|B{span3_start}|E{span3_end}"
                                fq_out.write(f"@{fq_header_3p}\n{seq_3p}\n+\n{qual_3p}\n")
                                round_info["emitted_3p_candidates"] += 1
                                candidates_written += 1

                    if candidates_written == 0:
                        stats["rounds"].append(round_info)
                        break

                    run_bowtie2(
                        bt2_index=args.bt2_index,
                        fastq_path=round_fastq,
                        sam_out_path=round_sam,
                        threads=args.threads,
                        seed_len=args.seed_len_cut,
                        score_min=args.score_min,
                    )

                    next_unmapped_queue: List[Tuple[str, str, str, int, int]] = []
                    with pysam.AlignmentFile(round_sam, "r") as r_sam_in:
                        for aln in r_sam_in:
                            if aln.is_secondary or aln.is_supplementary:
                                continue
                            parts = aln.query_name.split("|")
                            orig_qname = parts[0]
                            kind = parts[1]  # K5 or K3
                            seg_idx_val = int(parts[2][1:])
                            span_start = int(parts[3][1:])
                            span_end = int(parts[4][1:])

                            if not aln.is_unmapped and aln.mapping_quality >= args.min_mapq:
                                aln.query_name = orig_qname
                                aln.set_tag("ZM", args.mate, value_type="i")
                                aln.set_tag("ZI", seg_idx_val, value_type="i")
                                aln.set_tag("ZR", r, value_type="i")
                                aln.set_tag("ZS", span_start, value_type="i")
                                aln.set_tag("ZE", span_end, value_type="i")
                                aln.set_tag("ZK", kind, value_type="Z")
                                bam_out.write(aln)
                                stats["total_mapped_segments_q30"] += 1
                                segments_per_read[orig_qname] = segments_per_read.get(orig_qname, 0) + 1
                                if kind == "K5":
                                    round_info["mapped_5p_q30"] += 1
                                else:
                                    round_info["mapped_3p_q30"] += 1
                            else:
                                # If the 3' remainder did not map, pass it to the next round
                                # so it can be cut again at the next restriction site
                                if kind == "K3":
                                    seq_3p, qual_3p = get_original_seq_and_qual(aln)
                                    next_unmapped_queue.append(
                                        (orig_qname, seq_3p, qual_3p, span_start, seg_idx_val)
                                    )
                                    round_info["unmapped_3p_passed_to_next_round"] += 1

                    stats["rounds"].append(round_info)
                    unmapped_queue = next_unmapped_queue

        # Summarize segment counts per read
        for count in segments_per_read.values():
            if count == 0:
                stats["reads_with_0_segments"] += 1
            elif count == 1:
                stats["reads_with_1_segment"] += 1
            elif count == 2:
                stats["reads_with_2_segments"] += 1
            else:
                stats["reads_with_3plus_segments"] += 1

        # Queryname sort output BAM so all segments of a read are adjacent
        args.output_bam.parent.mkdir(parents=True, exist_ok=True)
        pysam.sort("-n", "-@", str(args.threads), "-o", str(args.output_bam), str(unsorted_bam_path))

    args.stats_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.stats_json, "w", encoding="utf-8") as f_json:
        json.dump(stats, f_json, indent=2)

    sys.stderr.write(
        f"[iterative_cut_and_map] Mate {args.mate}: {stats['total_input_reads']} reads -> "
        f"{stats['total_mapped_segments_q30']} mapped segments (Q>={args.min_mapq}) across "
        f"{1 + len(stats['rounds'])} rounds.\n"
    )


if __name__ == "__main__":
    main()
