#!/usr/bin/env python3
"""
Digest a reference genome FASTA into restriction fragments (BED6 format)
for Arima Hi-C+ or custom restriction enzyme combinations, compatible with HiC-Pro.

Default Arima Hi-C+ enzymes:
  - DpnII / MboI: ^GATC
  - HinfI:        G^ANTC
"""

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


IUPAC_MAP: Dict[str, str] = {
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
    "N": "[ACGT]",
}

PRESET_ENZYMES: Dict[str, List[str]] = {
    "arima": ["^GATC", "G^ANTC"],
    "arima_4e": ["^GATC", "G^ANTC", "C^TNAG", "T^TAA"],
    "dpnii": ["^GATC"],
    "mboi": ["^GATC"],
    "hindiii": ["A^AGCT_T"],
}


def motif_to_regex_and_offset(cut_motif: str) -> Tuple[str, int]:
    """
    Convert a restriction site motif with '^' indicating the 5' cleavage position
    (e.g. '^GATC' or 'G^ANTC') into a regex string and 0-based cleavage offset.
    """
    cleaned = cut_motif.strip().upper().replace("_", "")
    if "^" not in cleaned:
        raise ValueError(
            f"Restriction motif '{cut_motif}' must contain '^' to mark the cut position (e.g., '^GATC', 'G^ANTC')."
        )
    cut_offset = cleaned.index("^")
    raw_bases = cleaned.replace("^", "")
    regex_parts = []
    for base in raw_bases:
        if base not in IUPAC_MAP:
            raise ValueError(f"Unsupported IUPAC base '{base}' in motif '{cut_motif}'.")
        regex_parts.append(IUPAC_MAP[base])
    return "".join(regex_parts), cut_offset


def open_fasta(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def iter_fasta(path: Path):
    chrom = None
    seq_chunks: List[str] = []
    with open_fasta(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(seq_chunks).upper()
                chrom = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if chrom is not None:
            yield chrom, "".join(seq_chunks).upper()


def digest_chromosome(seq: str, patterns: List[Tuple[re.Pattern, int]]) -> List[Tuple[int, int]]:
    """
    Find all cleavage coordinates along `seq` and return 0-based half-open
    intervals [(start, end), ...] covering the full chromosome length.
    """
    cut_positions = {0, len(seq)}
    for regex, offset in patterns:
        # Use lookahead so overlapping restriction sites are both detected
        for match in regex.finditer(seq):
            cut_pos = match.start() + offset
            if 0 < cut_pos < len(seq):
                cut_positions.add(cut_pos)
    sorted_cuts = sorted(cut_positions)
    return [(sorted_cuts[i], sorted_cuts[i + 1]) for i in range(len(sorted_cuts) - 1)]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Digest reference genome FASTA into HiC-Pro compatible restriction fragment BED."
    )
    parser.add_argument(
        "-f", "--fasta", required=True, type=Path, help="Input reference genome FASTA (.fa or .fa.gz)"
    )
    parser.add_argument(
        "-r",
        "--restriction-sites",
        default="^GATC,G^ANTC",
        help="Comma-separated restriction cut sites with '^' (default: '^GATC,G^ANTC' for Arima Hi-C+) or preset name ('arima').",
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path, help="Output BED file of restriction fragments"
    )
    parser.add_argument(
        "--chrom-sizes",
        type=Path,
        default=None,
        help="Optional output path to also write chromosome sizes (chrom<TAB>size).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    raw_spec = args.restriction_sites.strip()
    if raw_spec.lower() in PRESET_ENZYMES:
        motifs = PRESET_ENZYMES[raw_spec.lower()]
    else:
        motifs = [m.strip() for m in raw_spec.split(",") if m.strip()]

    compiled_patterns: List[Tuple[re.Pattern, int]] = []
    for motif in motifs:
        regex_str, offset = motif_to_regex_and_offset(motif)
        # Positive lookahead allows overlapping motifs
        compiled_patterns.append((re.compile(f"(?=({regex_str}))"), offset))

    total_frags = 0
    chrom_sizes: List[Tuple[str, int]] = []

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out_bed:
        for chrom, seq in iter_fasta(args.fasta):
            chrom_sizes.append((chrom, len(seq)))
            intervals = digest_chromosome(seq, compiled_patterns)
            for idx, (start, end) in enumerate(intervals, start=1):
                frag_name = f"HIC_{chrom}_{idx}"
                out_bed.write(f"{chrom}\t{start}\t{end}\t{frag_name}\t0\t+\n")
                total_frags += 1

    if args.chrom_sizes:
        args.chrom_sizes.parent.mkdir(parents=True, exist_ok=True)
        with open(args.chrom_sizes, "w", encoding="utf-8") as cs_out:
            for chrom, size in chrom_sizes:
                cs_out.write(f"{chrom}\t{size}\n")

    sys.stderr.write(
        f"[digest_genome] Digested {len(chrom_sizes)} chromosomes with {motifs} -> {total_frags} fragments written to {args.output}\n"
    )


if __name__ == "__main__":
    main()
