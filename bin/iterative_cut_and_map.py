#!/usr/bin/env python3
"""
Iterative HiC-Pro style Bowtie2 mapper with multi-fragment 5' + 3' restriction rescue.
Streaming O(1)-memory architecture (no intermediate .sam files on disk, no full-dataset
Python lists/dicts in RAM) suitable for deep 200M+ read libraries.

Algorithm per mate (R1 or R2):
  1. Round 0 (Full-length mapping, optionally 5'-trimmed by --trim-5p):
     Streams reads into Bowtie2 (--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder)
     and parses SAM alignments directly from Bowtie2 stdout:
       - Keep alignments with MAPQ >= min_mapq (default 30) in accepted_unsorted.bam.
       - For unmapped / low-MAPQ reads, immediately locate the 5'-most ligation junction,
         split into 5' (K5) and 3' (K3) candidates, and stream them to round_1.fastq.
  2. Iterative Restriction Cutting & Rescue (Rounds 1 .. max_rounds):
     Streams round_r.fastq through Bowtie2 (-L 20) and parses stdout:
       - Keep any 5' (K5) or 3' (K3) fragment that maps with MAPQ >= min_mapq.
       - Any 3' remainder fragment (K3) that fails to map (MAPQ < min_mapq) is immediately
         checked for the next restriction ligation junction; if found, its 5' and 3' sub-pieces
         are streamed directly to round_{r+1}.fastq.
       - round_r.fastq is unlinked immediately after Round r completes.
  3. Queryname-sorts accepted_unsorted.bam into output_bam and computes exact per-read
     segment multiplicity stats (0, 1, 2, 3+ segments) in an O(1)-memory streaming pass.
"""

import argparse
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple

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


def split_and_emit_candidates(
    qname: str,
    seq: str,
    qual: str,
    offset_start: int,
    seg_idx: int,
    min_frag_len: int,
    junctions: List[Tuple[re.Pattern, int, str]],
    fq_out: Optional[TextIO],
    round_info: Dict[str, int],
) -> int:
    """
    Cut `seq` at its 5'-most ligation junction and write eligible 5' (K5) and 3' (K3)
    candidates to `fq_out`. Returns number of candidates written (0, 1, or 2).
    """
    round_info["input_unmapped_sequences"] += 1
    cut_res = find_first_ligation_cut(seq, junctions)
    if cut_res is None:
        round_info["sequences_without_junction"] += 1
        return 0

    round_info["sequences_with_junction"] += 1
    if fq_out is None:
        return 0

    cut_pos, _motif = cut_res
    seq_5p = seq[:cut_pos]
    qual_5p = qual[:cut_pos]
    span5_start = offset_start
    span5_end = offset_start + cut_pos

    seq_3p = seq[cut_pos:]
    qual_3p = qual[cut_pos:]
    span3_start = offset_start + cut_pos
    span3_end = offset_start + len(seq)

    written = 0
    next_idx = seg_idx
    if len(seq_5p) >= min_frag_len:
        fq_header_5p = f"{qname}|K5|S{next_idx}|B{span5_start}|E{span5_end}"
        fq_out.write(f"@{fq_header_5p}\n{seq_5p}\n+\n{qual_5p}\n")
        round_info["emitted_5p_candidates"] += 1
        written += 1
        next_idx += 1

    if len(seq_3p) >= min_frag_len:
        fq_header_3p = f"{qname}|K3|S{next_idx}|B{span3_start}|E{span3_end}"
        fq_out.write(f"@{fq_header_3p}\n{seq_3p}\n+\n{qual_3p}\n")
        round_info["emitted_3p_candidates"] += 1
        written += 1

    return written


def make_round_info(r: int) -> Dict[str, int]:
    return {
        "round": r,
        "input_unmapped_sequences": 0,
        "sequences_with_junction": 0,
        "sequences_without_junction": 0,
        "emitted_5p_candidates": 0,
        "emitted_3p_candidates": 0,
        "mapped_5p_q30": 0,
        "mapped_3p_q30": 0,
        "unmapped_3p_passed_to_next_round": 0,
    }


def start_bowtie2_stream(
    bt2_index: str,
    fastq_path: Path,
    threads: int,
    seed_len: int,
    score_min: str,
    trim_5p: int = 0,
) -> Tuple[subprocess.Popen, Optional[subprocess.Popen]]:
    """
    Launch Bowtie2 streaming SAM to stdout.
    If trim_5p > 0, streams 5'-trimmed reads via a fast awk pipe into Bowtie2 stdin (-U -).
    """
    bt2_cmd = [
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
    ]

    if trim_5p > 0:
        start_idx = trim_5p + 1
        cat_cmd = "gzip -dc" if str(fastq_path).endswith(".gz") else "cat"
        trim_shell = f"{cat_cmd} '{fastq_path}' | awk 'NR%2==0 {{print substr($0, {start_idx}); next}} {{print}}'"
        trim_proc = subprocess.Popen(
            ["/bin/bash", "-c", trim_shell],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        bt2_cmd.extend(["-U", "-"])
        bt2_proc = subprocess.Popen(
            bt2_cmd,
            stdin=trim_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if trim_proc.stdout is not None:
            trim_proc.stdout.close()
        return bt2_proc, trim_proc
    else:
        bt2_cmd.extend(["-U", str(fastq_path)])
        bt2_proc = subprocess.Popen(
            bt2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return bt2_proc, None


def finish_bowtie2_stream(
    bt2_proc: subprocess.Popen, trim_proc: Optional[subprocess.Popen]
) -> None:
    bt2_err = bt2_proc.stderr.read().decode("utf-8", errors="replace") if bt2_proc.stderr else ""
    bt2_rc = bt2_proc.wait()
    if trim_proc is not None:
        trim_err = trim_proc.stderr.read().decode("utf-8", errors="replace") if trim_proc.stderr else ""
        trim_rc = trim_proc.wait()
        if trim_rc != 0:
            sys.stderr.write(f"[iterative_cut_and_map] 5'-trimming pipe failed:\n{trim_err}\n")
            raise RuntimeError(f"5'-trimming pipe exited with code {trim_rc}")
    if bt2_rc != 0:
        sys.stderr.write(f"[iterative_cut_and_map] Bowtie2 failed:\n{bt2_err}\n")
        raise RuntimeError(f"Bowtie2 exited with code {bt2_rc}")


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

    with tempfile.TemporaryDirectory(prefix="hicpro_iter_") as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        unsorted_bam_path = tmpdir / "accepted_unsorted.bam"
        round1_fastq = tmpdir / "round_1.fastq"
        round1_info = make_round_info(1)
        round1_candidates = 0

        # Round 0: Full-length mapping streamed directly from Bowtie2 stdout
        bt2_proc, trim_proc = start_bowtie2_stream(
            bt2_index=args.bt2_index,
            fastq_path=args.fastq,
            threads=args.threads,
            seed_len=args.seed_len_full,
            score_min=args.score_min,
            trim_5p=args.trim_5p,
        )

        with pysam.AlignmentFile(bt2_proc.stdout, "r") as sam_in:
            header = sam_in.header.to_dict()
            with pysam.AlignmentFile(unsorted_bam_path, "wb", header=header) as bam_out:
                with open(round1_fastq, "w", encoding="utf-8") as fq1_out:
                    for aln in sam_in:
                        if aln.is_secondary or aln.is_supplementary:
                            continue
                        qname = aln.query_name.split()[0]
                        stats["total_input_reads"] += 1

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
                        else:
                            stats["round0_unmapped_or_lowq"] += 1
                            if args.max_rounds >= 1:
                                seq, qual = get_original_seq_and_qual(aln)
                                round1_candidates += split_and_emit_candidates(
                                    qname=qname,
                                    seq=seq,
                                    qual=qual,
                                    offset_start=0,
                                    seg_idx=0,
                                    min_frag_len=args.min_frag_len,
                                    junctions=junctions,
                                    fq_out=fq1_out,
                                    round_info=round1_info,
                                )

                finish_bowtie2_stream(bt2_proc, trim_proc)

                # Iterative cutting & mapping rounds (1 .. max_rounds)
                curr_fastq = round1_fastq
                curr_info = round1_info
                curr_candidates = round1_candidates

                for r in range(1, args.max_rounds + 1):
                    if curr_candidates == 0:
                        if curr_info["input_unmapped_sequences"] > 0:
                            stats["rounds"].append(curr_info)
                        curr_fastq.unlink(missing_ok=True)
                        break

                    next_fastq = tmpdir / f"round_{r + 1}.fastq"
                    next_info = make_round_info(r + 1)
                    next_candidates = 0

                    r_bt2_proc, _ = start_bowtie2_stream(
                        bt2_index=args.bt2_index,
                        fastq_path=curr_fastq,
                        threads=args.threads,
                        seed_len=args.seed_len_cut,
                        score_min=args.score_min,
                        trim_5p=0,
                    )

                    with pysam.AlignmentFile(r_bt2_proc.stdout, "r") as r_sam_in, open(
                        next_fastq, "w", encoding="utf-8"
                    ) as next_fq_out:
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
                                if kind == "K5":
                                    curr_info["mapped_5p_q30"] += 1
                                else:
                                    curr_info["mapped_3p_q30"] += 1
                            else:
                                if kind == "K3":
                                    curr_info["unmapped_3p_passed_to_next_round"] += 1
                                    if r < args.max_rounds:
                                        seq_3p, qual_3p = get_original_seq_and_qual(aln)
                                        next_candidates += split_and_emit_candidates(
                                            qname=orig_qname,
                                            seq=seq_3p,
                                            qual=qual_3p,
                                            offset_start=span_start,
                                            seg_idx=seg_idx_val,
                                            min_frag_len=args.min_frag_len,
                                            junctions=junctions,
                                            fq_out=next_fq_out,
                                            round_info=next_info,
                                        )

                    finish_bowtie2_stream(r_bt2_proc, None)
                    curr_fastq.unlink(missing_ok=True)
                    stats["rounds"].append(curr_info)

                    curr_fastq = next_fastq
                    curr_info = next_info
                    curr_candidates = next_candidates

        # Queryname sort output BAM so all segments of a read are adjacent
        args.output_bam.parent.mkdir(parents=True, exist_ok=True)
        pysam.sort("-n", "-@", str(args.threads), "-o", str(args.output_bam), str(unsorted_bam_path))

    # Streaming O(1)-memory pass over queryname-sorted output BAM to count segments per read
    mapped_reads_seen = 0
    with pysam.AlignmentFile(args.output_bam, "rb") as sorted_bam:
        for _qname, group in itertools.groupby(sorted_bam, key=lambda a: a.query_name):
            mapped_reads_seen += 1
            seg_count = sum(1 for _ in group)
            if seg_count == 1:
                stats["reads_with_1_segment"] += 1
            elif seg_count == 2:
                stats["reads_with_2_segments"] += 1
            else:
                stats["reads_with_3plus_segments"] += 1

    stats["reads_with_0_segments"] = max(0, stats["total_input_reads"] - mapped_reads_seen)

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
