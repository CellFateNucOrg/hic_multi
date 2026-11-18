#!/usr/bin/env bash
#SBATCH --job-name=nf_hic_multicontact
#SBATCH --partition=all
#SBATCH --ntasks-per-core=2
#SBATCH --hint=multithread
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=nf_hic_multicontact_%j.log
#SBATCH --error=nf_hic_multicontact_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ==============================================================================
# Launcher script for hicpro-multicontact-nf on izblisbon.izb.unibe.ch (SLURM)
# Dispatches all compute tasks to SLURM compute nodes (--partition=all) and
# caps total concurrent logical CPUs across own jobs to <= 50% of Available CPUs.
#
# Usage (from izblisbon login node):
#   sbatch submit_izblisbon.sh
# ==============================================================================

chmod +x bin/*.py

# Query live cluster usage on partition 'all' to enforce <= 50% of Available CPUs
TOTAL_UP=$(sinfo -h -N -p all -o "%N %C %T" 2>/dev/null | awk '$3 !~ /down|drain|drng|maint|fail/ {split($2,a,"/"); s += a[4]} END {print s+0}')
OTHER_USED=$(squeue -h -t R -p all -o "%u %C" 2>/dev/null | awk -v me="$(whoami)" '$1 != me {s += $2} END {print s+0}')
AVAIL_CPUS=$(( TOTAL_UP - OTHER_USED ))
[[ "${AVAIL_CPUS}" -lt 4 ]] && AVAIL_CPUS=16
MAX_OWN_CPUS=$(( AVAIL_CPUS / 2 ))
[[ "${MAX_OWN_CPUS}" -lt 2 ]] && MAX_OWN_CPUS=2

# Each ITERATIVE_BOWTIE2_MAP task uses 8 logical CPUs by default (or MAX_OWN_CPUS if smaller)
TASK_CPUS=8
[[ "${MAX_OWN_CPUS}" -lt "${TASK_CPUS}" ]] && TASK_CPUS="${MAX_OWN_CPUS}"
MAX_CONCURRENT_TASKS=$(( MAX_OWN_CPUS / TASK_CPUS ))
[[ "${MAX_CONCURRENT_TASKS}" -lt 1 ]] && MAX_CONCURRENT_TASKS=1

echo "[SLURM Policy] Total online CPUs: ${TOTAL_UP}, Other users allocated: ${OTHER_USED}, Available: ${AVAIL_CPUS}"
echo "[SLURM Policy] Capping own usage to <= 50% of available (${MAX_OWN_CPUS} CPUs -> max ${MAX_CONCURRENT_TASKS} concurrent ${TASK_CPUS}-CPU tasks)"

# Configure paths for your C. elegans reference (WBcel235 / ce11) and samplesheet
SAMPLESHEET="assets/samplesheet.csv"
GENOME_FASTA="/path/to/genomes/Caenorhabditis_elegans/WBcel235/ce11.fa"
OUTDIR="results_celegans_arima"

export NXF_OPTS="-Xms1G -Xmx4G"

nextflow run main.nf \
    -profile izblisbon,conda \
    --slurm_partition 'all' \
    --slurm_cluster_options '--ntasks-per-core=2 --hint=multithread' \
    --slurm_queue_size "${MAX_CONCURRENT_TASKS}" \
    --max_cpus "${TASK_CPUS}" \
    --input "${SAMPLESHEET}" \
    --fasta "${GENOME_FASTA}" \
    --outdir "${OUTDIR}" \
    --restriction_sites '^GATC,G^ANTC' \
    --ligation_motifs 'arima' \
    --min_mapq 30 \
    --min_frag_len 20 \
    --max_rounds 4 \
    --fastq_chunk_size 2500000 \
    --cool_base_bin 500 \
    --cool_resolutions '500,1000,2000,5000,10000,25000,50000,100000' \
    -resume
