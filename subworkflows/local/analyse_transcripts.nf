//
// ANALYSE TRANSCRITPS
//

 include { STRINGTIE_STRINGTIE         } from '../../modules/nf-core/stringtie/stringtie/main'
 include { GFFCOMPARE                  } from '../../modules/nf-core/gffcompare/main'
 include { GENERATE_COUNTS_DROP        } from '../../modules/local/generate_gene_counts4drop'
 include { GENERATE_ANNOTATION_DROP    } from '../../modules/local/generate_sample_annot'
 include { GENERATE_CONFIG_RUN_AE_DROP } from '../../modules/local/generate_config_runAEdrop'

workflow ANALYSE_TRANSCRIPTS {
    take:
        ch_bam               // channel (mandatory): [ val(meta), [ path(bam) ] ]
        ch_gtf               // channel (mandatory): [ path(gtf) ]
        ch_fasta_fai_meta    // channel (mandatory): [ val(meta), [ path(fasta), path(fai) ]
        gene_counts          // channel [val(meta), path(tsv)] 
        reference_count_file // channel [ path(tsv) ]

    main:
        ch_versions = Channel.empty()

        // DROP
        // Generates count files for samples and merges them with reference count file
        star_count = gene_counts.map{ meta, cnt_file -> cnt_file }.collect()
        star_samp  = gene_counts.map{ meta, cnt_file -> meta }.collect()
        GENERATE_COUNTS_DROP(star_count, star_samp, ch_gtf, reference_count_file)

        // Annotation file
        GENERATE_ANNOTATION_DROP(GENERATE_COUNTS_DROP.out.processed_gene_counts, ch_gtf)
       
        // Config file
        GENERATE_CONFIG_RUN_AE_DROP(
            ch_fasta_fai_meta, 
            ch_gtf, 
            GENERATE_ANNOTATION_DROP.out.sample_annotation_drop,
            GENERATE_COUNTS_DROP.out.processed_gene_counts
        )

        // Stringtie
        STRINGTIE_STRINGTIE(
            ch_bam,
            ch_gtf
        )

        // Compare stringtie results to reference
        GFFCOMPARE(
            STRINGTIE_STRINGTIE.out.transcript_gtf,
            ch_fasta_fai_meta,
            ch_gtf.map{ gtf -> [ [id:gtf.simpleName], gtf ] }
        )

        ch_versions = ch_versions.mix(GENERATE_COUNTS_DROP.out.versions)
        ch_versions = ch_versions.mix(GENERATE_ANNOTATION_DROP.out.versions)
        ch_versions = ch_versions.mix(GENERATE_CONFIG_RUN_AE_DROP.out.versions)
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())

    emit:
        transcript_gtf        = STRINGTIE_STRINGTIE.out.transcript_gtf              // channel: [ val(meta), [ path(transctript_gtf)] ]
        abundance             = STRINGTIE_STRINGTIE.out.abundance                   // channel: [ val(meta), [ path(abundance) ] ]
        coverage_gtf          = STRINGTIE_STRINGTIE.out.coverage_gtf                // channel: [ val(meta), [ path(coverage_gtf) ] ]
        annotated_gtf         = GFFCOMPARE.out.annotated_gtf                        // channel: [ val(meta), [ path(annotated_gtf) ] ]
        stats_gtf             = GFFCOMPARE.out.stats                                // channel: [ val(meta), [ path(stats) ] ]
        processed_gene_counts = GENERATE_COUNTS_DROP.out.processed_gene_counts      // channel: [ path(tsv) ]
        annotation_drop       = GENERATE_ANNOTATION_DROP.out.sample_annotation_drop // channel: [ path(tsv) ]
        config_drop           = GENERATE_CONFIG_RUN_AE_DROP.out.config_drop         // channel: [ path(confg_file.yml) ]
        drop_ae_out           = GENERATE_CONFIG_RUN_AE_DROP.out.drop_ae_out
        versions              = ch_versions                                         // channel: [ path(versions.yml) ]
}
