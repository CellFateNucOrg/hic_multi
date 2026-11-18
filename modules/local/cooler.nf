process BUILD_COOLER {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/cooler/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(all_valid_pairs)
    path chrom_sizes

    output:
    tuple val(meta), path("${meta.id}.${params.cool_base_bin}.cool"), emit: cool
    tuple val(meta), path("${meta.id}.mcool")                       , emit: mcool

    script:
    def balance_flag = params.balance_cooler ? "--balance" : ""
    """
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 5 -p2 6 \\
        ${chrom_sizes}:${params.cool_base_bin} \\
        ${all_valid_pairs} \\
        ${meta.id}.${params.cool_base_bin}.cool

    cooler zoomify \\
        --nproc ${task.cpus} \\
        --resolutions ${params.cool_resolutions} \\
        ${balance_flag} \\
        -o ${meta.id}.mcool \\
        ${meta.id}.${params.cool_base_bin}.cool
    """
}
