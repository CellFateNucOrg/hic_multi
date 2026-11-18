#!/usr/bin/env nextflow
/*
 * =========================================================================================
 *  hicpro-multicontact-nf: HiC-Pro Iterative 5'+3' Multi-Fragment Rescue & All-Pairs Pipeline
 * =========================================================================================
 *  Organism target : C. elegans (WBcel235 / ce11)
 *  Chemistry       : Arima Hi-C+ (^GATC, G^ANTC) + KAPA Library Prep (2x150 bp PE)
 *  Cluster target  : izblisbon (SLURM executor, nf-core/hic style DSL2 architecture)
 * =========================================================================================
 */

nextflow.enable.dsl = 2

include { DIGEST_GENOME; BOWTIE2_BUILD } from './modules/local/prepare_genome'
include { ITERATIVE_BOWTIE2_MAP        } from './modules/local/iterative_map'
include { PAIR_MULTICONTACTS           } from './modules/local/pair_multicontacts'
include { MERGE_DEDUP_PAIRS            } from './modules/local/merge_pairs'
include { BUILD_COOLER                 } from './modules/local/cooler'

workflow {

    if (!params.input) {
        error "Please provide an input samplesheet via --input <samplesheet.csv>"
    }
    if (!params.fasta) {
        error "Please provide a reference genome FASTA via --fasta <ce11.fa>"
    }

    // 1. Prepare Reference Genome & Arima Restriction Fragment BED
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).first()

    if (params.restriction_bed && params.chrom_sizes) {
        ch_restriction_bed = Channel.fromPath(params.restriction_bed, checkIfExists: true).first()
        ch_chrom_sizes     = Channel.fromPath(params.chrom_sizes, checkIfExists: true).first()
    } else {
        DIGEST_GENOME(ch_fasta, params.restriction_sites)
        ch_restriction_bed = DIGEST_GENOME.out.restriction_bed.first()
        ch_chrom_sizes     = DIGEST_GENOME.out.chrom_sizes.first()
    }

    if (params.bowtie2_index) {
        ch_bt2_index = Channel
            .fromPath("${params.bowtie2_index}*.bt2*", checkIfExists: true)
            .collect()
            .map { files -> [ file(params.bowtie2_index).name, files ] }
            .first()
    } else {
        BOWTIE2_BUILD(ch_fasta)
        ch_bt2_index = BOWTIE2_BUILD.out.index.first()
    }

    // 2. Read Samplesheet and Split Paired-End FASTQs into Parallel SLURM Chunks (like nf-core/hic)
    ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [ id: row.sample ]
            def fq1  = file(row.fastq_1, checkIfExists: true)
            def fq2  = file(row.fastq_2, checkIfExists: true)
            return [ meta, fq1, fq2 ]
        }

    if (params.fastq_chunk_size && params.fastq_chunk_size > 0) {
        ch_chunked_pairs = ch_samplesheet
            .splitFastq(by: params.fastq_chunk_size, pe: true, file: true)
            .map { meta, fq1_chunk, fq2_chunk ->
                // Extract deterministic chunk identifier from splitFastq filename
                def chunk_tag = fq1_chunk.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                return [ meta, chunk_tag, fq1_chunk, fq2_chunk ]
            }
    } else {
        ch_chunked_pairs = ch_samplesheet
            .map { meta, fq1, fq2 ->
                def chunk_tag = fq1.simpleName
                return [ meta, chunk_tag, fq1, fq2 ]
            }
    }

    // Fan out R1 (mate 1) and R2 (mate 2) for independent single-end iterative Bowtie2 mapping
    ch_mates_for_mapping = ch_chunked_pairs
        .flatMap { meta, chunk_id, fq1, fq2 ->
            [
                [ meta, chunk_id, 1, fq1 ],
                [ meta, chunk_id, 2, fq2 ]
            ]
        }

    // 3. Iterative Full-Length -> 5' & 3' Restriction Cutting & Mapping (Q30)
    ITERATIVE_BOWTIE2_MAP(
        ch_mates_for_mapping,
        ch_bt2_index
    )

    // Re-pair R1 and R2 BAMs for each (meta, chunk_id)
    ch_bams_by_chunk = ITERATIVE_BOWTIE2_MAP.out.bam
        .map { meta, chunk_id, mate, bam ->
            [ [ meta, chunk_id ], mate, bam ]
        }
        .groupTuple(by: 0, size: 2)
        .map { key, mates, bams ->
            def meta     = key[0]
            def chunk_id = key[1]
            def idx_r1   = mates.indexOf(1)
            def idx_r2   = mates.indexOf(2)
            return [ meta, chunk_id, bams[idx_r1], bams[idx_r2] ]
        }

    // 4. Multi-Contact Pairing (AB, AC, BC ...) & HiC-Pro Fragment Filtering
    PAIR_MULTICONTACTS(
        ch_bams_by_chunk,
        ch_restriction_bed
    )

    // 5. Merge Chunks per Biological Sample & Remove PCR Duplicates
    ch_valid_pairs_grouped = PAIR_MULTICONTACTS.out.valid_pairs.groupTuple(by: 0)
    ch_map_stats_grouped   = ITERATIVE_BOWTIE2_MAP.out.stats.groupTuple(by: 0)
    ch_pair_stats_grouped  = PAIR_MULTICONTACTS.out.stats.groupTuple(by: 0)

    ch_for_merge = ch_valid_pairs_grouped
        .join(ch_map_stats_grouped, by: 0)
        .join(ch_pair_stats_grouped, by: 0)

    MERGE_DEDUP_PAIRS(ch_for_merge)

    // 6. Build Cooler (.cool) & Balanced Multi-Resolution (.mcool) Contact Matrices
    if (!params.skip_cooler) {
        BUILD_COOLER(
            MERGE_DEDUP_PAIRS.out.all_valid_pairs,
            ch_chrom_sizes
        )
    }
}
