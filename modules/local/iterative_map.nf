process ITERATIVE_BOWTIE2_MAP {
    tag "${meta.id}_chunk${chunk_id}_R${mate}"
    label 'process_medium'
    publishDir "${params.outdir}/mapping/bams/${meta.id}", mode: params.publish_dir_mode, enabled: params.save_intermediate_bams, pattern: "*.bam"
    publishDir "${params.outdir}/mapping/stats/${meta.id}", mode: params.publish_dir_mode, pattern: "*.json"

    input:
    tuple val(meta), val(chunk_id), val(mate), path(fastq_chunk)
    tuple val(index_name), path(bt2_index_files)

    output:
    tuple val(meta), val(chunk_id), val(mate), path("${meta.id}_chunk${chunk_id}_R${mate}.namesorted.bam"), emit: bam
    tuple val(meta), path("${meta.id}_chunk${chunk_id}_R${mate}.mapstats.json")                           , emit: stats

    script:
    def prefix = "${meta.id}_chunk${chunk_id}_R${mate}"
    def bt2_file = bt2_index_files.find { it.toString().endsWith('.1.bt2') || it.toString().endsWith('.1.bt2l') } ?: bt2_index_files[0]
    def bt2_prefix = bt2_file.toString().replaceAll(/(\.rev)?\.[0-9]+\.bt2l?$/, '')
    """
    iterative_cut_and_map.py \\
        --fastq ${fastq_chunk} \\
        --bt2-index ${bt2_prefix} \\
        --mate ${mate} \\
        --output-bam ${prefix}.namesorted.bam \\
        --stats-json ${prefix}.mapstats.json \\
        --min-mapq ${params.min_mapq} \\
        --ligation-motifs "${params.ligation_motifs}" \\
        --min-frag-len ${params.min_frag_len} \\
        --max-rounds ${params.max_rounds} \\
        --trim-5p ${params.trim_5p} \\
        --seed-len-full ${params.seed_len_full} \\
        --seed-len-cut ${params.seed_len_cut} \\
        --score-min "${params.score_min}" \\
        --threads ${task.cpus}
    """
}
