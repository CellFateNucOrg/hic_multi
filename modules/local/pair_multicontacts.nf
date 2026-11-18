process PAIR_MULTICONTACTS {
    tag "${meta.id}_chunk${chunk_id}"
    label 'process_medium'
    publishDir "${params.outdir}/pairing/chunks/${meta.id}", mode: params.publish_dir_mode, pattern: "*.json"
    publishDir "${params.outdir}/pairing/dumped/${meta.id}", mode: params.publish_dir_mode, enabled: params.save_dumped_pairs, pattern: "*.dumpedPairs"

    input:
    tuple val(meta), val(chunk_id), path(bam_r1), path(bam_r2)
    path restriction_bed

    output:
    tuple val(meta), path("${meta.id}_chunk${chunk_id}.validPairs")       , emit: valid_pairs
    tuple val(meta), path("${meta.id}_chunk${chunk_id}.pairstats.json")   , emit: stats
    tuple val(meta), path("${meta.id}_chunk${chunk_id}.dumpedPairs")      , optional: true, emit: dumped_pairs

    script:
    def prefix = "${meta.id}_chunk${chunk_id}"
    def dedup_flag = params.dedup_overlapping_mates ? "" : "--no-dedup-overlapping-mates"
    def relig_flag = params.discard_religation ? "--discard-religation" : ""
    def dumped_arg = params.save_dumped_pairs ? "--output-dumped-pairs ${prefix}.dumpedPairs" : ""
    """
    pair_multicontacts.py \\
        --bam-r1 ${bam_r1} \\
        --bam-r2 ${bam_r2} \\
        --restriction-bed ${restriction_bed} \\
        --output-valid-pairs ${prefix}.validPairs \\
        --stats-json ${prefix}.pairstats.json \\
        --threads ${task.cpus} \\
        ${dedup_flag} \\
        ${relig_flag} \\
        ${dumped_arg}
    """
}
