process DIGEST_GENOME {
    tag "$fasta.simpleName"
    label 'process_low'
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    path fasta
    val restriction_sites

    output:
    path "${fasta.simpleName}_restriction_fragments.bed", emit: restriction_bed
    path "${fasta.simpleName}.chrom.sizes"              , emit: chrom_sizes

    script:
    """
    digest_genome.py \\
        --fasta ${fasta} \\
        --restriction-sites "${restriction_sites}" \\
        --output ${fasta.simpleName}_restriction_fragments.bed \\
        --chrom-sizes ${fasta.simpleName}.chrom.sizes
    """
}

process BOWTIE2_BUILD {
    tag "$fasta.simpleName"
    label 'process_medium'
    publishDir "${params.outdir}/genome/bowtie2_index", mode: params.publish_dir_mode, enabled: params.save_reference

    input:
    path fasta

    output:
    tuple val(fasta.simpleName), path("bowtie2_index/*"), emit: index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/${fasta.simpleName}
    """
}
