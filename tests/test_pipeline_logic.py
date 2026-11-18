#!/usr/bin/env python3
"""
Unit tests verifying:
  1. Arima restriction digestion (^GATC, G^ANTC) in digest_genome.py
  2. Iterative 5' + 3' ligation junction splitting across multiple rounds in iterative_cut_and_map.py
  3. Multi-contact all-pairs generation (R1: A & B, R2: C -> 3 contacts AB, AC, BC)
     and overlapping mate read-through deduplication in pair_multicontacts.py
"""

import re
import sys
import tempfile
import unittest
from pathlib import Path

# Add bin/ to import path
BIN_DIR = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import digest_genome as dg  # noqa: E402


class TestArimaMulticontactLogic(unittest.TestCase):
    def test_arima_genome_digestion(self):
        # Construct synthetic C. elegans chromosome I sequence with DpnII (^GATC) at 50 and HinfI (G^ANTC) at 120
        seq = ("A" * 50) + "GATC" + ("C" * 65) + "GAATC" + ("T" * 80)
        patterns = []
        for motif in ["^GATC", "G^ANTC"]:
            regex_str, offset = dg.motif_to_regex_and_offset(motif)
            patterns.append((re.compile(f"(?=({regex_str}))"), offset))

        intervals = dg.digest_chromosome(seq, patterns)
        # Expected cuts:
        #   DpnII ^GATC at index 50 -> cut at 50
        #   HinfI G^AATC at index 50+4+65 = 119 -> cut at 119 + 1 = 120
        self.assertEqual(intervals, [(0, 50), (50, 120), (120, len(seq))])

    def test_iterative_5p_and_3p_cutting_multi_round(self):
        # Test the exact junction regexes and multi-round 5' + 3' splitting logic
        # Suppose a 150bp read contains 3 ligated fragments: F1 (45bp) + F2 (50bp) + F3 (55bp)
        # Junction 1 (DpnII-DpnII: GATCGATC): F1 ends with GATC, F2 starts with GATC
        # Junction 2 (HinfI-DpnII: GANTCGATC): F2 ends with GAATC, F3 starts with GATC
        frag1 = ("A" * 41) + "GATC"          # 45 bp
        frag2 = "GATC" + ("C" * 41) + "GAATC"  # 50 bp
        frag3 = "GATC" + ("T" * 51)          # 55 bp
        chimeric_read = frag1 + frag2 + frag3
        self.assertEqual(len(chimeric_read), 150)

        # Import build_junction_regexes and find_first_ligation_cut without requiring pysam if not installed
        # We load those functions directly from iterative_cut_and_map.py source
        src = (BIN_DIR / "iterative_cut_and_map.py").read_text(encoding="utf-8")
        ns = {"re": re, "Dict": dict, "List": list, "Optional": object, "Tuple": tuple}
        # Execute helper definitions up to run_bowtie2
        helper_code = src.split("import pysam")[1].split("def run_bowtie2")[0]
        exec(helper_code, ns)

        junctions = ns["build_junction_regexes"]("arima")

        # Round 1: full read fails to map -> cut at first junction
        cut1 = ns["find_first_ligation_cut"](chimeric_read, junctions)
        self.assertIsNotNone(cut1)
        cut_pos1, motif1 = cut1
        self.assertEqual(cut_pos1, 45)
        self.assertEqual(motif1, "GATCGATC")

        piece_5p_r1 = chimeric_read[:cut_pos1]
        piece_3p_r1 = chimeric_read[cut_pos1:]
        self.assertEqual(piece_5p_r1, frag1)
        self.assertEqual(piece_3p_r1, frag2 + frag3)

        # Piece 5' (frag1) maps! Piece 3' (frag2 + frag3) is still chimeric so it fails Round 1 mapping
        # Round 2: cut the unmapped 3' remainder (piece_3p_r1) at the next restriction site
        cut2 = ns["find_first_ligation_cut"](piece_3p_r1, junctions)
        self.assertIsNotNone(cut2)
        cut_pos2, motif2 = cut2
        self.assertEqual(cut_pos2, 50)
        self.assertEqual(motif2, "GANTCGATC")

        piece_5p_r2 = piece_3p_r1[:cut_pos2]
        piece_3p_r2 = piece_3p_r1[cut_pos2:]
        self.assertEqual(piece_5p_r2, frag2)
        self.assertEqual(piece_3p_r2, frag3)

        # Neither frag2 nor frag3 has any remaining ligation junction
        self.assertIsNone(ns["find_first_ligation_cut"](piece_3p_r2, junctions))

    def test_ab_ac_bc_three_contacts_generation(self):
        # Load MappedSegment, deduplicate_readthrough_segments, and order_pair from pair_multicontacts.py
        src = (BIN_DIR / "pair_multicontacts.py").read_text(encoding="utf-8")
        import bisect
        import itertools
        from dataclasses import dataclass

        ns = {
            "bisect": bisect,
            "itertools": itertools,
            "dataclass": dataclass,
            "Path": Path,
            "Dict": dict,
            "Iterable": list,
            "List": list,
            "Optional": object,
            "Tuple": tuple,
        }
        helper_code = src.split("import pysam")[1].split("def extract_segment")[0]
        helper_code_2 = "def deduplicate_readthrough_segments" + src.split("def deduplicate_readthrough_segments")[1].split("def parse_args")[0]
        exec(helper_code + "\n" + helper_code_2, ns)

        MappedSegment = ns["MappedSegment"]
        dedup_fn = ns["dedup_fn"] = ns["deduplicate_readthrough_segments"]
        order_pair = ns["order_pair"]

        # Molecule with R1 containing fragments A and B, and R2 containing fragment C
        seg_A = MappedSegment(
            qname="mol1", mate=1, seg_idx=0, round_num=1, span_start=0, span_end=75, kind="K5",
            chrom="I", chrom_rank=0, pos_5p_1based=1050, midpoint_0based=1080, strand="+", mapq=42,
            aln_len=75, rfrag_name="HIC_I_10", rfrag_num=10, rfrag_len=300
        )
        seg_B = MappedSegment(
            qname="mol1", mate=1, seg_idx=1, round_num=1, span_start=75, span_end=150, kind="K3",
            chrom="I", chrom_rank=0, pos_5p_1based=55000, midpoint_0based=55035, strand="-", mapq=40,
            aln_len=75, rfrag_name="HIC_I_240", rfrag_num=240, rfrag_len=420
        )
        seg_C = MappedSegment(
            qname="mol1", mate=2, seg_idx=0, round_num=0, span_start=0, span_end=150, kind="FULL",
            chrom="II", chrom_rank=1, pos_5p_1based=9000, midpoint_0based=9075, strand="+", mapq=42,
            aln_len=150, rfrag_name="HIC_II_45", rfrag_num=45, rfrag_len=500
        )

        segments, collapsed = dedup_fn([seg_A, seg_B, seg_C])
        self.assertEqual(collapsed, 0)
        self.assertEqual(len(segments), 3)

        pairs = [order_pair(a, b) for a, b in itertools.combinations(segments, 2)]
        pair_rfrags = [(p[0].rfrag_name, p[1].rfrag_name) for p in pairs]
        self.assertEqual(
            pair_rfrags,
            [
                ("HIC_I_10", "HIC_I_240"),   # Contact AB
                ("HIC_I_10", "HIC_II_45"),   # Contact AC
                ("HIC_I_240", "HIC_II_45"),  # Contact BC
            ],
        )


if __name__ == "__main__":
    unittest.main()
