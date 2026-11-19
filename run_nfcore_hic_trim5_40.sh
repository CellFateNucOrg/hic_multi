#!/usr/bin/env bash
#SBATCH --job-name=nfcore_hic_trim5_40
#SBATCH --partition=all
#SBATCH --ntasks-per-core=2
#SBATCH --hint=multithread
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/meister.data/pmeister/20260904_arima_HN2_2_shallow_nfcore_trim5_40/nfcore_hic_trim5_40_%j.log
#SBATCH --error=/mnt/meister.data/pmeister/20260904_arima_HN2_2_shallow_nfcore_trim5_40/nfcore_hic_trim5_40_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

WORK_DIR="/mnt/meister.data/pmeister/20260904_arima_HN2_2_shallow_nfcore_trim5_40"
mkdir -p "${WORK_DIR}/fastq"
cd "${WORK_DIR}"

cat << 'EOS' > "${WORK_DIR}/samplesheet.csv"
sample,fastq_1,fastq_2
HN2_2_trim5_40,/mnt/meister.data/pmeister/20260904_arima_HN2_2_shallow_nfcore_trim5_40/fastq/HN2_2_trim5_40_R1.fastq.gz,/mnt/meister.data/pmeister/20260904_arima_HN2_2_shallow_nfcore_trim5_40/fastq/HN2_2_trim5_40_R2.fastq.gz
EOS

RAW_R1="/mnt/meister.data/mdas/20260904_arima_HN2_2_shallow/fastq/HN2_2_S1_L001_R1_001.fastq.gz"
RAW_R2="/mnt/meister.data/mdas/20260904_arima_HN2_2_shallow/fastq/HN2_2_S1_L001_R2_001.fastq.gz"
TRIM_R1="${WORK_DIR}/fastq/HN2_2_trim5_40_R1.fastq.gz"
TRIM_R2="${WORK_DIR}/fastq/HN2_2_trim5_40_R2.fastq.gz"

# 1. Trim first 5 nt and last 40 nt from R1 and R2 in parallel on the compute node
if [[ ! -s "${TRIM_R1}" || ! -s "${TRIM_R2}" ]]; then
    echo "[Trim] Trimming first 5 nt and last 40 nt from R1 and R2..."
    gzip -dc "${RAW_R1}" | awk 'NR%2==0 {print substr($0, 6, length($0)-45); next} {print}' | gzip -1 > "${TRIM_R1}.tmp" &
    PID1=$!
    gzip -dc "${RAW_R2}" | awk 'NR%2==0 {print substr($0, 6, length($0)-45); next} {print}' | gzip -1 > "${TRIM_R2}.tmp" &
    PID2=$!
    wait "${PID1}" "${PID2}"
    mv "${TRIM_R1}.tmp" "${TRIM_R1}"
    mv "${TRIM_R2}.tmp" "${TRIM_R2}"
    echo "[Trim] Done."
fi

# 2. Enforce SLURM cluster policy: <= 50% of Available CPUs
TOTAL_UP=$(sinfo -h -N -p all -o "%N %C %T" 2>/dev/null | awk '$3 !~ /down|drain|drng|maint|fail/ {split($2,a,"/"); s += a[4]} END {print s+0}')
OTHER_USED=$(squeue -h -t R -p all -o "%u %C" 2>/dev/null | awk -v me="$(whoami)" '$1 != me {s += $2} END {print s+0}')
OWN_RUNNING=$(squeue -h -t R -p all -u "$(whoami)" -o "%C" 2>/dev/null | awk '{s += $1} END {print s+0}')
AVAIL_CPUS=$(( TOTAL_UP - OTHER_USED ))
MAX_OWN_CPUS=$(( AVAIL_CPUS / 2 ))
REMAINING_CPUS=$(( MAX_OWN_CPUS - OWN_RUNNING ))
[[ "${REMAINING_CPUS}" -lt 8 ]] && REMAINING_CPUS=8
QUEUE_SIZE=$(( REMAINING_CPUS / 8 ))
[[ "${QUEUE_SIZE}" -lt 1 ]] && QUEUE_SIZE=1
[[ "${QUEUE_SIZE}" -gt 10 ]] && QUEUE_SIZE=10

echo "[SLURM Policy] Total: ${TOTAL_UP}, Other used: ${OTHER_USED}, Avail: ${AVAIL_CPUS}, Max own: ${MAX_OWN_CPUS}, Own running: ${OWN_RUNNING} -> Nextflow queueSize=${QUEUE_SIZE}"

HICPRO_PATH="/home/pmeister/HiC-Pro_install/HiC-Pro_3.1.0/HiC-Pro/scripts:/home/pmeister/HiC-Pro_install/HiC-Pro_3.1.0/bin/utils:/home/pmeister/miniconda3/envs/HiC-Pro_v3.1.0/bin:/home/pmeister/miniconda3/envs/coolns/bin:/home/pmeister/miniconda3/envs/hicexplorer/bin:/home/pmeister/miniconda3/envs/deeptools/bin"
COOL_PATH="/home/pmeister/miniconda3/envs/coolns/bin:${HICPRO_PATH}"

cat << EOC > "${WORK_DIR}/policy_override.config"
singularity {
    enabled = false
}
executor {
    name = 'slurm'
    queueSize = ${QUEUE_SIZE}
}
process {
    executor = 'slurm'
    queue = 'all'
    maxRetries = 2
    scratch = "/scratch/\${System.getenv('USER')}"
    clusterOptions = "--ntasks-per-core=2 --hint=multithread"
    beforeScript = "mkdir -p /scratch/\\\$USER && export PYTHONNOUSERSITE=1 PATH=${HICPRO_PATH}:\\\$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"

    withName: '.*COOLER.*|.*COOLTOOLS.*|.*DUMPMATRIX.*' {
        beforeScript = "mkdir -p /scratch/\\\$USER && export PYTHONNOUSERSITE=1 PATH=${COOL_PATH}:\\\$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"
    }
}
EOC

source /home/pmeister/miniconda3/bin/activate env_nf2
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"
export NXF_SYNTAX_PARSER=v1

GENOME_FILE="/mnt/meister.data/pmeister/WBcel235/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
BWT2_DIR="/mnt/meister.data/pmeister/WBcel235/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index"

nextflow run nf-core/hic \
    -r 2.1.0 \
    -c "${WORK_DIR}/policy_override.config" \
    -w "${WORK_DIR}/work" \
    --input "${WORK_DIR}/samplesheet.csv" \
    --fasta "${GENOME_FILE}" \
    --bwt2_index "${BWT2_DIR}" \
    --digestion "arima" \
    --bin_size "2000,4000,10000,20000,200000" \
    --skip_fastqc \
    --skip_multiqc \
    --skip_compartments \
    --skip_tads \
    --outdir "${WORK_DIR}/results" \
    -resume
