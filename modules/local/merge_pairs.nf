process MERGE_DEDUP_PAIRS {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/hic_results/data/${meta.id}", mode: params.publish_dir_mode, pattern: "*.allValidPairs"
    publishDir "${params.outdir}/hic_results/stats/${meta.id}", mode: params.publish_dir_mode, pattern: "*.{mergestat,json}"

    input:
    tuple val(meta), path(valid_pairs_chunks), path(map_stats_jsons), path(pair_stats_jsons)

    output:
    tuple val(meta), path("${meta.id}.allValidPairs")         , emit: all_valid_pairs
    tuple val(meta), path("${meta.id}_allValidPairs.mergestat"), emit: mergestat
    tuple val(meta), path("${meta.id}.multicontact_qc.json")  , emit: qc_json

    script:
    def keep_dup_flag = params.keep_duplicates ? "--keep-duplicates" : ""
    def sort_mem = task.memory ? "${(task.memory.toGiga() * 0.7).intValue() ?: 2}G" : "4G"
    """
    merge_valid_pairs.py \\
        --sample ${meta.id} \\
        --valid-pairs ${valid_pairs_chunks} \\
        --map-stats ${map_stats_jsons} \\
        --pair-stats ${pair_stats_jsons} \\
        --output-all-valid-pairs ${meta.id}.allValidPairs \\
        --output-mergestat ${meta.id}_allValidPairs.mergestat \\
        --output-summary-json ${meta.id}.multicontact_qc.json \\
        --threads ${task.cpus} \\
        --sort-mem ${sort_mem} \\
        ${keep_dup_flag}
    """
}
