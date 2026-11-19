# hic_multi

Nextflow DSL2 Hi-C mapping and multi-contact extraction pipeline (`hic_multi`) modeled on [HiC-Pro](https://github.com/nservant/HiC-Pro) and [nf-core/hic](https://nf-co.re/hic/2.1.0/), tailored for **_C. elegans_ Arima Hi-C+ / KAPA 2x150 bp paired-end libraries** running as SLURM jobs on **`izblisbon`**.

---

## Key Algorithm: Iterative 5' + 3' Restriction Rescue & Multi-Contact Pairing

Unlike standard HiC-Pro—which cuts unmapped chimeric reads at the first restriction ligation junction and discards the 3' portion—this pipeline rescues **both the 5' and 3' restriction fragments iteratively across multiple rounds** and forms **all pairwise contacts $\binom{k}{2}$** per sequenced DNA molecule:

1. **Round 0 (Full-length mapping)**:
   - Aligns full-length R1 and R2 reads independently using `bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder`.
   - Keeps alignments with `MAPQ >= 30` (`--min_mapq 30`).
   - Passes unmapped or `MAPQ < 30` reads to Round 1.
2. **Round $r = 1 \dots R_{\max}$ (Iterative 5' + 3' restriction junction splitting & mapping)**:
   - Scans each unmapped sequence from 5' to 3' for the earliest Arima ligation junction (`^GATC` + `G^ANTC` blunt-end fill-in motifs: `GATCGATC`, `GATCANTC`, `GANTCGATC`, `GANTCANTC`, `GANTGATC`, `GATCGANTC`, `GANTANTC`).
   - Splits the read at the exact ligation point into:
     - **5' fragment** (if length $\ge 20$ bp) $\rightarrow$ mapped with `bowtie2` (`-L 20 --score-min L,-0.6,-0.2 --end-to-end`), kept if `MAPQ >= 30`.
     - **3' remainder fragment** (if length $\ge 20$ bp) $\rightarrow$ **also mapped full-length in the same round** with `bowtie2`, kept if `MAPQ >= 30`.
   - If the **3' remainder fragment fails to map (`MAPQ < 30`)**, it is passed into Round $r+1$, split at the *next* restriction ligation junction into a 5' sub-fragment and 3' sub-remainder, and mapped again (up to `--max_rounds`, default `4`).
3. **Multi-Contact All-Pairs Generation (`AB`, `AC`, `BC`)**:
   - All uniquely mapped (`MAPQ >= 30`) fragments from R1 and R2 belonging to the same read pair (`QNAME`) are assigned to their _C. elegans_ restriction fragment (`digest_genome.py`).
   - Redundant read-through observations (where 2x150 bp mates overlap in the middle of a short insert and sequence the exact same restriction fragment from opposite ends) are collapsed so they do not produce duplicate contacts or artificial self-pairs.
   - For a molecule with $k \ge 2$ mapped restriction fragments (e.g., fragments `A` and `B` on R1, and fragment `C` on R2), all $\binom{k}{2}$ pairwise combinations are generated:
     - **`AB`** (intra-mate R1–R1 contact)
     - **`AC`** (inter-mate R1–R2 contact)
     - **`BC`** (inter-mate R1–R2 contact)
   - Each pair is filtered using HiC-Pro's restriction fragment rules (removing Dangling Ends `DE` and Self-Circles `SC` on the same restriction fragment), deduplicated for PCR duplicates, and exported as:
     - Standard 12-column HiC-Pro **`.allValidPairs`** and **`.mergestat`**
     - Multi-contact QC JSON report (**`.multicontact_qc.json`**)
     - Cooler **`.cool`** (`500 bp` base bin) and ICE-balanced multi-resolution **`.mcool`** (`500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 25 kb, 50 kb, 100 kb`).

---

## Quick Start on `izblisbon` (SLURM)

1. Prepare `assets/samplesheet.csv`:
   ```csv
   sample,fastq_1,fastq_2
   N2_ArimaHiC_rep1,/path/to/N2_rep1_R1.fastq.gz,/path/to/N2_rep1_R2.fastq.gz
   N2_ArimaHiC_rep2,/path/to/N2_rep2_R1.fastq.gz,/path/to/N2_rep2_R2.fastq.gz
   ```

2. Edit `GENOME_FASTA` in `submit_izblisbon.sh` and submit to SLURM:
   ```bash
   sbatch submit_izblisbon.sh
   ```
   Or run directly:
   ```bash
   nextflow run main.nf \
       -profile izblisbon,conda \
       --input assets/samplesheet.csv \
       --fasta /path/to/WBcel235/ce11.fa \
       --outdir results_celegans
   ```
   *(If `izblisbon` is configured as a compute node name rather than a SLURM partition name on your cluster, pass `--slurm_partition '' --slurm_cluster_options '--nodelist=izblisbon'`.)*
