#!/usr/bin/env python3
"""
Multi-contact pairing and HiC-Pro restriction fragment filtering.

Given queryname-sorted BAMs for Mate 1 (R1) and Mate 2 (R2) produced by
`iterative_cut_and_map.py` and a restriction fragment BED file from `digest_genome.py`:
  1. Groups all uniquely mapped (MAPQ >= min_mapq) segments per sequenced DNA molecule (QNAME).
  2. Assigns each segment to its restriction fragment in the reference genome using binary search.
  3. Optionally collapses redundant read-through observations when R1 and R2 sequence across
     the exact same restriction fragment in a multi-fragment chimeric read.
  4. Generates all pairwise contacts C(k, 2) among the k mapped fragments of the molecule
     (e.g., if R1 has fragments A & B and R2 has fragment C -> emits 3 contacts: AB, AC, BC).
  5. Classifies each contact pair following HiC-Pro rules:
       - Valid cis / Valid trans
       - Dangling End (same restriction fragment, inward +/-)
       - Self-Circle (same restriction fragment, outward -/+)
       - Same-Fragment same-strand (dumped)
       - Adjacent restriction fragment religation (tracked in stats; optionally filtered)
  6. Writes standard 12-column HiC-Pro `.validPairs` and JSON summary statistics.
"""

import argparse
import bisect
import itertools
import json
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pysam


@dataclass
class RestrictionBedIndex:
    chrom_order: Dict[str, int]
    starts_by_chrom: Dict[str, List[int]]
    ends_by_chrom: Dict[str, List[int]]
    names_by_chrom: Dict[str, List[str]]
    nums_by_chrom: Dict[str, List[int]]

    @classmethod
    def from_bed(cls, bed_path: Path) -> "RestrictionBedIndex":
        chrom_order: Dict[str, int] = {}
        starts_by_chrom: Dict[str, List[int]] = {}
        ends_by_chrom: Dict[str, List[int]] = {}
        names_by_chrom: Dict[str, List[str]] = {}
        nums_by_chrom: Dict[str, List[int]] = {}

        with open(bed_path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                fields = line.split("\t")
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                frag_name = fields[3] if len(fields) > 3 else f"HIC_{chrom}_{ len(starts_by_chrom.get(chrom, [])) + 1 }"

                if chrom not in chrom_order:
                    chrom_order[chrom] = len(chrom_order)
                    starts_by_chrom[chrom] = []
                    ends_by_chrom[chrom] = []
                    names_by_chrom[chrom] = []
                    nums_by_chrom[chrom] = []

                idx_num = len(starts_by_chrom[chrom]) + 1
                starts_by_chrom[chrom].append(start)
                ends_by_chrom[chrom].append(end)
                names_by_chrom[chrom].append(frag_name)
                nums_by_chrom[chrom].append(idx_num)

        return cls(
            chrom_order=chrom_order,
            starts_by_chrom=starts_by_chrom,
            ends_by_chrom=ends_by_chrom,
            names_by_chrom=names_by_chrom,
            nums_by_chrom=nums_by_chrom,
        )

    def lookup(self, chrom: str, pos_0based: int) -> Tuple[str, int, int]:
        """
        Return (frag_name, frag_num, frag_length) for a 0-based genomic position.
        """
        starts = self.starts_by_chrom.get(chrom)
        if not starts:
            return f"HIC_{chrom}_0", 0, 0
        idx = bisect.bisect_right(starts, pos_0based) - 1
        if idx < 0:
            idx = 0
        elif idx >= len(starts):
            idx = len(starts) - 1
        frag_len = self.ends_by_chrom[chrom][idx] - starts[idx]
        return self.names_by_chrom[chrom][idx], self.nums_by_chrom[chrom][idx], frag_len


@dataclass
class MappedSegment:
    qname: str
    mate: int
    seg_idx: int
    round_num: int
    span_start: int
    span_end: int
    kind: str
    chrom: str
    chrom_rank: int
    pos_5p_1based: int
    midpoint_0based: int
    strand: str
    mapq: int
    aln_len: int
    rfrag_name: str
    rfrag_num: int
    rfrag_len: int


def extract_segment(aln: pysam.AlignedSegment, rbed: RestrictionBedIndex) -> Optional[MappedSegment]:
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        return None
    chrom = aln.reference_name
    if chrom is None:
        return None
    ref_start = aln.reference_start
    ref_end = aln.reference_end or (ref_start + (aln.query_length or 1))
    strand = "-" if aln.is_reverse else "+"
    pos_5p_1based = ref_end if strand == "-" else (ref_start + 1)
    midpoint_0based = (ref_start + ref_end) // 2

    rfrag_name, rfrag_num, rfrag_len = rbed.lookup(chrom, midpoint_0based)
    chrom_rank = rbed.chrom_order.get(chrom, 999999)

    mate = int(aln.get_tag("ZM")) if aln.has_tag("ZM") else 1
    seg_idx = int(aln.get_tag("ZI")) if aln.has_tag("ZI") else 0
    round_num = int(aln.get_tag("ZR")) if aln.has_tag("ZR") else 0
    span_start = int(aln.get_tag("ZS")) if aln.has_tag("ZS") else 0
    span_end = int(aln.get_tag("ZE")) if aln.has_tag("ZE") else (aln.query_length or 0)
    kind = str(aln.get_tag("ZK")) if aln.has_tag("ZK") else "FULL"

    return MappedSegment(
        qname=aln.query_name,
        mate=mate,
        seg_idx=seg_idx,
        round_num=round_num,
        span_start=span_start,
        span_end=span_end,
        kind=kind,
        chrom=chrom,
        chrom_rank=chrom_rank,
        pos_5p_1based=pos_5p_1based,
        midpoint_0based=midpoint_0based,
        strand=strand,
        mapq=aln.mapping_quality,
        aln_len=ref_end - ref_start,
        rfrag_name=rfrag_name,
        rfrag_num=rfrag_num,
        rfrag_len=rfrag_len,
    )


def deduplicate_readthrough_segments(segments: List[MappedSegment]) -> Tuple[List[MappedSegment], int]:
    """
    In a 2x150bp library, R1 and R2 can overlap in the middle of the insert and both sequence
    the same central restriction fragment from opposite strands.
    When a molecule has >= 3 mapped segments across R1 and R2 and an R2 segment maps to the exact
    same restriction fragment (rfrag_name) as an R1 segment, collapse the duplicate observation
    so a 2-fragment chimeric read doesn't generate duplicate contacts or a spurious self-pair,
    and a 3-fragment read where R1 & R2 overlap on the middle fragment yields 3 contacts (AB, AC, BC).

    Note: If a molecule has only 2 segments (1 on R1, 1 on R2) on the same restriction fragment,
    we do NOT collapse them so standard Dangling Ends and Self-Circles are properly classified.
    """
    if len(segments) < 3:
        return segments, 0

    r1_segs = [s for s in segments if s.mate == 1]
    r2_segs = [s for s in segments if s.mate == 2]
    if not r1_segs or not r2_segs:
        return segments, 0

    r1_rfrags = {s.rfrag_name for s in r1_segs}
    collapsed_r2: List[MappedSegment] = []
    collapsed_count = 0
    for s2 in r2_segs:
        if s2.rfrag_name in r1_rfrags:
            collapsed_count += 1
        else:
            collapsed_r2.append(s2)

    return r1_segs + collapsed_r2, collapsed_count


def order_pair(s1: MappedSegment, s2: MappedSegment) -> Tuple[MappedSegment, MappedSegment]:
    key1 = (s1.chrom_rank, s1.chrom, s1.pos_5p_1based)
    key2 = (s2.chrom_rank, s2.chrom, s2.pos_5p_1based)
    if key1 <= key2:
        return s1, s2
    return s2, s1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pair multi-fragment Hi-C alignments into all C(k,2) contacts (AB, AC, BC, ...)."
    )
    parser.add_argument("--bam-r1", required=True, type=Path, help="Queryname-sorted R1 BAM from iterative_cut_and_map.py")
    parser.add_argument("--bam-r2", required=True, type=Path, help="Queryname-sorted R2 BAM from iterative_cut_and_map.py")
    parser.add_argument("-r", "--restriction-bed", required=True, type=Path, help="Restriction fragment BED file")
    parser.add_argument("-o", "--output-valid-pairs", required=True, type=Path, help="Output .validPairs file")
    parser.add_argument("--output-dumped-pairs", type=Path, default=None, help="Optional output file for DE/SC/dumped pairs")
    parser.add_argument("--stats-json", required=True, type=Path, help="Output pairing & multi-contact stats JSON")
    parser.add_argument(
        "--no-dedup-overlapping-mates",
        action="store_true",
        help="Disable collapsing overlapping R1/R2 read-through segments on the same restriction fragment in >=3-segment molecules",
    )
    parser.add_argument(
        "--discard-religation",
        action="store_true",
        help="Filter out inward pairs on immediately adjacent restriction fragments (religation) instead of keeping them as valid pairs",
    )
    parser.add_argument("-p", "--threads", type=int, default=4, help="Threads for BAM merge (default: 4)")
    return parser.parse_args()


def main():
    args = parse_args()
    rbed = RestrictionBedIndex.from_bed(args.restriction_bed)

    stats = {
        "molecules_with_mapped_segments": 0,
        "molecules_1_segment_singleton": 0,
        "molecules_2_segments": 0,
        "molecules_3_segments": 0,
        "molecules_4plus_segments": 0,
        "collapsed_overlapping_mate_segments": 0,
        "total_candidate_pairs_generated": 0,
        "valid_pairs_total": 0,
        "valid_pairs_cis": 0,
        "valid_pairs_cis_short_lt_1kb": 0,
        "valid_pairs_cis_1kb_to_20kb": 0,
        "valid_pairs_cis_gt_20kb": 0,
        "valid_pairs_trans": 0,
        "valid_pairs_from_2frag_molecules": 0,
        "valid_pairs_from_3frag_molecules": 0,
        "valid_pairs_from_4plus_frag_molecules": 0,
        "valid_pairs_intra_mate_rescued": 0,
        "valid_pairs_inter_mate": 0,
        "dangling_end_pairs": 0,
        "self_circle_pairs": 0,
        "same_fragment_other_pairs": 0,
        "religation_adjacent_pairs": 0,
    }

    args.output_valid_pairs.parent.mkdir(parents=True, exist_ok=True)
    dumped_handle = (
        open(args.output_dumped_pairs, "w", encoding="utf-8")
        if args.output_dumped_pairs is not None
        else None
    )

    with tempfile.TemporaryDirectory(prefix="hicpro_pair_") as tmpdir_str:
        merged_bam_path = Path(tmpdir_str) / "merged_namesorted.bam"
        # Merge already queryname-sorted R1 and R2 BAMs in O(N) streaming pass
        pysam.merge(
            "-n",
            "-@",
            str(args.threads),
            "-f",
            str(merged_bam_path),
            str(args.bam_r1),
            str(args.bam_r2),
        )

        with pysam.AlignmentFile(merged_bam_path, "rb") as bam_in, open(
            args.output_valid_pairs, "w", encoding="utf-8"
        ) as vp_out:
            for qname, aln_iter in itertools.groupby(bam_in, key=lambda a: a.query_name):
                raw_segments: List[MappedSegment] = []
                for aln in aln_iter:
                    seg = extract_segment(aln, rbed)
                    if seg is not None:
                        raw_segments.append(seg)

                if not raw_segments:
                    continue

                # Sort segments deterministically: Mate 1 (5'->3') then Mate 2 (5'->3')
                raw_segments.sort(key=lambda s: (s.mate, s.seg_idx, s.span_start))

                if not args.no_dedup_overlapping_mates:
                    segments, collapsed_n = deduplicate_readthrough_segments(raw_segments)
                    stats["collapsed_overlapping_mate_segments"] += collapsed_n
                else:
                    segments = raw_segments

                k = len(segments)
                stats["molecules_with_mapped_segments"] += 1
                if k == 1:
                    stats["molecules_1_segment_singleton"] += 1
                    continue
                elif k == 2:
                    stats["molecules_2_segments"] += 1
                elif k == 3:
                    stats["molecules_3_segments"] += 1
                else:
                    stats["molecules_4plus_segments"] += 1

                # Generate all C(k, 2) pairwise combinations (e.g. AB, AC, BC for k=3)
                for pair_idx, (seg_a, seg_b) in enumerate(itertools.combinations(segments, 2), start=1):
                    stats["total_candidate_pairs_generated"] += 1
                    s1, s2 = order_pair(seg_a, seg_b)
                    pair_id = f"{qname}__C{pair_idx}_K{k}_M{seg_a.mate}{seg_b.mate}"

                    is_same_chrom = s1.chrom == s2.chrom
                    is_same_rfrag = is_same_chrom and (s1.rfrag_name == s2.rfrag_name)

                    if is_same_rfrag:
                        dist = abs(s2.pos_5p_1based - s1.pos_5p_1based)
                        if s1.strand == "+" and s2.strand == "-":
                            stats["dangling_end_pairs"] += 1
                            pair_type = "DE"
                        elif s1.strand == "-" and s2.strand == "+":
                            stats["self_circle_pairs"] += 1
                            pair_type = "SC"
                        else:
                            stats["same_fragment_other_pairs"] += 1
                            pair_type = "SAME_FRAG"

                        if dumped_handle is not None:
                            dumped_handle.write(
                                f"{pair_id}\t{s1.chrom}\t{s1.pos_5p_1based}\t{s1.strand}\t"
                                f"{s2.chrom}\t{s2.pos_5p_1based}\t{s2.strand}\t{dist}\t"
                                f"{s1.rfrag_name}\t{s2.rfrag_name}\t{s1.mapq}\t{s2.mapq}\t{pair_type}\n"
                            )
                        continue

                    # Check adjacent restriction fragment religation
                    is_adjacent_religation = (
                        is_same_chrom
                        and abs(s1.rfrag_num - s2.rfrag_num) == 1
                        and s1.strand == "+"
                        and s2.strand == "-"
                    )
                    if is_adjacent_religation:
                        stats["religation_adjacent_pairs"] += 1
                        if args.discard_religation:
                            if dumped_handle is not None:
                                dist = abs(s2.pos_5p_1based - s1.pos_5p_1based)
                                dumped_handle.write(
                                    f"{pair_id}\t{s1.chrom}\t{s1.pos_5p_1based}\t{s1.strand}\t"
                                    f"{s2.chrom}\t{s2.pos_5p_1based}\t{s2.strand}\t{dist}\t"
                                    f"{s1.rfrag_name}\t{s2.rfrag_name}\t{s1.mapq}\t{s2.mapq}\tRELIG\n"
                                )
                            continue

                    # Valid Hi-C contact!
                    stats["valid_pairs_total"] += 1
                    if seg_a.mate == seg_b.mate:
                        stats["valid_pairs_intra_mate_rescued"] += 1
                    else:
                        stats["valid_pairs_inter_mate"] += 1

                    if k == 2:
                        stats["valid_pairs_from_2frag_molecules"] += 1
                    elif k == 3:
                        stats["valid_pairs_from_3frag_molecules"] += 1
                    else:
                        stats["valid_pairs_from_4plus_frag_molecules"] += 1

                    if is_same_chrom:
                        dist = abs(s2.pos_5p_1based - s1.pos_5p_1based)
                        stats["valid_pairs_cis"] += 1
                        if dist < 1000:
                            stats["valid_pairs_cis_short_lt_1kb"] += 1
                        elif dist <= 20000:
                            stats["valid_pairs_cis_1kb_to_20kb"] += 1
                        else:
                            stats["valid_pairs_cis_gt_20kb"] += 1
                    else:
                        dist = -1
                        stats["valid_pairs_trans"] += 1

                    # HiC-Pro 12-column validPairs format:
                    # read_name chr1 pos1 strand1 chr2 pos2 strand2 insert_size rfrag1 rfrag2 mapq1 mapq2
                    vp_out.write(
                        f"{pair_id}\t{s1.chrom}\t{s1.pos_5p_1based}\t{s1.strand}\t"
                        f"{s2.chrom}\t{s2.pos_5p_1based}\t{s2.strand}\t{dist}\t"
                        f"{s1.rfrag_name}\t{s2.rfrag_name}\t{s1.mapq}\t{s2.mapq}\n"
                    )

    if dumped_handle is not None:
        dumped_handle.close()

    args.stats_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.stats_json, "w", encoding="utf-8") as f_json:
        json.dump(stats, f_json, indent=2)

    sys.stderr.write(
        f"[pair_multicontacts] Processed {stats['molecules_with_mapped_segments']} molecules -> "
        f"{stats['valid_pairs_total']} valid pairs "
        f"(2-frag: {stats['valid_pairs_from_2frag_molecules']}, "
        f"3-frag [AB,AC,BC]: {stats['valid_pairs_from_3frag_molecules']}, "
        f"4+frag: {stats['valid_pairs_from_4plus_frag_molecules']}).\n"
    )


if __name__ == "__main__":
    main()
