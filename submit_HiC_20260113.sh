#!/usr/bin/env bash
#SBATCH --job-name=hic_multi_20260113
#SBATCH --partition=all
#SBATCH --ntasks-per-core=2
#SBATCH --hint=multithread
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=96:00:00
#SBATCH --output=/mnt/meister.data/FischleLab_KarthikEswara/HiC_20260113/hic_multi/hic_multi_%j.log
#SBATCH --error=/mnt/meister.data/FischleLab_KarthikEswara/HiC_20260113/hic_multi/hic_multi_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ==============================================================================
# Launcher script for hic_multi on all 14 libraries in:
#   /mnt/meister.data/FischleLab_KarthikEswara/HiC_20260113/fastq
# Outputs stored in:
#   /mnt/meister.data/FischleLab_KarthikEswara/HiC_20260113/hic_multi
# ==============================================================================

source /home/pmeister/miniconda3/bin/activate env_nf2

export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"
export NXF_SYNTAX_PARSER=v1

PIPELINE_DIR="/home/pmeister/hic_multi"
chmod +x "${PIPELINE_DIR}/bin/"*.py

PROJECT_DIR="/mnt/meister.data/FischleLab_KarthikEswara/HiC_20260113"
RUN_BASE="${PROJECT_DIR}/hic_multi"
mkdir -p "${RUN_BASE}"
cd "${RUN_BASE}"

# Build samplesheet pointing to full-length 150 bp reads in fastq/ (from samples.csv)
ORIG_SAMPLESHEET="${PROJECT_DIR}/samples.csv"
FASTQ_SAMPLESHEET="${RUN_BASE}/samples_fastq.csv"
sed 's|/clipped_reads/|/fastq/|g' "${ORIG_SAMPLESHEET}" > "${FASTQ_SAMPLESHEET}"

echo "[Samplesheet] Generated ${FASTQ_SAMPLESHEET} pointing to ${PROJECT_DIR}/fastq:"
cat "${FASTQ_SAMPLESHEET}"

# Query live cluster usage on partition 'all' to enforce <= 50% of Available CPUs
TOTAL_UP=$(sinfo -h -N -p all -o "%N %C %T" 2>/dev/null | awk '$3 !~ /down|drain|drng|maint|fail/ {split($2,a,"/"); s += a[4]} END {print s+0}')
OTHER_USED=$(squeue -h -t R -p all -o "%u %C" 2>/dev/null | awk -v me="$(whoami)" '$1 != me {s += $2} END {print s+0}')
OWN_RUNNING=$(squeue -h -t R -p all -u "$(whoami)" -o "%C" 2>/dev/null | awk '{s += $1} END {print s+0}')
AVAIL_CPUS=$(( TOTAL_UP - OTHER_USED ))
[[ "${AVAIL_CPUS}" -lt 16 ]] && AVAIL_CPUS=32
MAX_OWN_CPUS=$(( AVAIL_CPUS / 2 ))
REMAINING_OWN_CPUS=$(( MAX_OWN_CPUS - OWN_RUNNING ))
# Ensure at least 4 concurrent 8-CPU tasks (32 CPUs) can run as existing short array tasks rotate
[[ "${REMAINING_OWN_CPUS}" -lt 32 ]] && REMAINING_OWN_CPUS=32

TASK_CPUS=8
MAX_CONCURRENT_TASKS=$(( REMAINING_OWN_CPUS / TASK_CPUS ))
[[ "${MAX_CONCURRENT_TASKS}" -lt 4 ]] && MAX_CONCURRENT_TASKS=4

echo "[SLURM Policy] Total online CPUs: ${TOTAL_UP}, Other users allocated: ${OTHER_USED}, Available: ${AVAIL_CPUS}"
echo "[SLURM Policy] Max own CPUs (50% available): ${MAX_OWN_CPUS}, Own currently running: ${OWN_RUNNING} -> Nextflow queueSize=${MAX_CONCURRENT_TASKS} (${TASK_CPUS}-CPU tasks)"

GENOME_FASTA="/mnt/meister.data/pmeister/WBcel235/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
BOWTIE2_INDEX="/mnt/meister.data/pmeister/WBcel235/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/genome"
OUTDIR="${RUN_BASE}"
WORKDIR="${RUN_BASE}/work"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile izblisbon \
    -w "${WORKDIR}" \
    --slurm_partition 'all' \
    --slurm_cluster_options '--ntasks-per-core=2 --hint=multithread' \
    --slurm_queue_size "${MAX_CONCURRENT_TASKS}" \
    --max_cpus "${TASK_CPUS}" \
    --input "${FASTQ_SAMPLESHEET}" \
    --fasta "${GENOME_FASTA}" \
    --bowtie2_index "${BOWTIE2_INDEX}" \
    --outdir "${OUTDIR}" \
    --restriction_sites '^GATC,G^ANTC' \
    --ligation_motifs 'arima' \
    --trim_5p 5 \
    --min_mapq 30 \
    --min_frag_len 20 \
    --max_rounds 4 \
    --fastq_chunk_size 0 \
    --cool_base_bin 500 \
    --cool_resolutions '500,1000,2000,4000,5000,10000,20000,25000,50000,100000,200000' \
    -resume
