//
// ANALYSE TRANSCRITPS
//

 include { STRINGTIE_STRINGTIE } from '../../modules/nf-core/stringtie/stringtie/main'
 include { GFFCOMPARE          } from '../../modules/nf-core/gffcompare/main'
 include { DROP_COUNTS         } from '../../modules/local/drop_counts'
 include { DROP_SAMPLE_ANNOT   } from '../../modules/local/drop_sample_annot'
 include { DROP_CONFIG_RUN_AE  } from '../../modules/local/drop_config_runAE'
 include { DROP_CONFIG_RUN_AS  } from '../../modules/local/drop_config_runAS'
 include { DROP_FILTER_RESULTS } from '../../modules/local/drop_filter_results'

workflow ANALYSE_TRANSCRIPTS {
    take:
        ch_bam_bai                    // channel (mandatory): [ val(meta), [ path(bam) ],[ path(bai) ] ]
        ch_gtf                        // channel (mandatory): [ path(gtf) ]
        ch_fasta_fai_meta             // channel (mandatory): [ val(meta), [ path(fasta), path(fai) ]
        gene_counts                   // channel [val(meta), path(tsv)] 
        ch_ref_drop_count_file        // channel [ path(tsv) ]
        ch_ref_drop_annot_file        // channel [ path(tsv) ]
        ch_ref_drop_splice_folder     // channel [ path(folder) ]
        genome                        // channel [val(genome)]
        drop_padjcutoff_ae            // channel [val(drop_padjcutoff_ae)]
        drop_padjcutoff_as            // channel [val(drop_padjcutoff_as)]
        drop_zscorecutoff             // channel [val(drop_zscorecutoff)]
        ch_gene_panel_clinical_filter // channel [ path(tsv) ]

    main:
        ch_versions = Channel.empty()

        // DROP
        // Generates count files for samples and merges them with reference count file
        star_count = gene_counts.map{ meta, cnt_file -> cnt_file }.collect()
        star_samples  = gene_counts.map{ meta, cnt_file -> meta }.collect()
        DROP_COUNTS(
            star_count,
            star_samples,
            ch_gtf,
            ch_ref_drop_count_file
        )
        ch_drop_counts = DROP_COUNTS.out.processed_gene_counts.collect()
        
        // Generates sample annotation
        ch_bam_files = ch_bam_bai.collect{it[1]}
        DROP_SAMPLE_ANNOT(
            ch_bam_files,
            star_samples,
            ch_drop_counts,
            ch_ref_drop_annot_file,
            ch_gtf
        )

        // Generates  config file and runs Aberrant expression module
        DROP_CONFIG_RUN_AE(
            ch_fasta_fai_meta, 
            ch_gtf, 
            DROP_SAMPLE_ANNOT.out.drop_annot,
            ch_drop_counts,
            genome,
            drop_padjcutoff_ae,
            drop_zscorecutoff
        )
        
        // Generates  config file and runs Aberrant splicing module
        ch_bai_files = ch_bam_bai.collect{it[2]}.toList()
        ch_bam_bai_files=ch_bam_files.toList().combine(ch_bai_files)
        DROP_CONFIG_RUN_AS(
            ch_fasta_fai_meta, 
            ch_gtf, 
            DROP_SAMPLE_ANNOT.out.drop_annot,
            ch_bam_bai_files,
            ch_ref_drop_splice_folder,
            genome,
            drop_padjcutoff_as
        )

        ch_out_drop_ae_rds    = DROP_CONFIG_RUN_AE.out.drop_ae_rds    ? DROP_CONFIG_RUN_AE.out.drop_ae_rds.collect()
                                                                        : Channel.empty()
        ch_out_drop_gene_name = DROP_CONFIG_RUN_AE.out.drop_gene_name ? DROP_CONFIG_RUN_AE.out.drop_gene_name.collect()
                                                                        : Channel.empty()
        ch_out_drop_as_tsv    = DROP_CONFIG_RUN_AS.out.drop_as_tsv    ? DROP_CONFIG_RUN_AS.out.drop_as_tsv.collect()
                                                                        : Channel.empty()

        DROP_FILTER_RESULTS(
            star_samples,
            ch_gene_panel_clinical_filter,
            ch_out_drop_ae_rds,
            ch_out_drop_gene_name,
            ch_out_drop_as_tsv
        )

        // Stringtie
        ch_bam = ch_bam_bai.map{ meta, bam, bai -> [meta, [bam]]}
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

        ch_versions = ch_versions.mix(DROP_COUNTS.out.versions)
        ch_versions = ch_versions.mix(DROP_SAMPLE_ANNOT.out.versions)
        ch_versions = ch_versions.mix(DROP_CONFIG_RUN_AE.out.versions)
        ch_versions = ch_versions.mix(DROP_CONFIG_RUN_AS.out.versions)
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())

    emit:
        transcript_gtf        = STRINGTIE_STRINGTIE.out.transcript_gtf     // channel: [ val(meta), [ path(transctript_gtf)] ]
        abundance             = STRINGTIE_STRINGTIE.out.abundance          // channel: [ val(meta), [ path(abundance) ] ]
        coverage_gtf          = STRINGTIE_STRINGTIE.out.coverage_gtf       // channel: [ val(meta), [ path(coverage_gtf) ] ]
        annotated_gtf         = GFFCOMPARE.out.annotated_gtf               // channel: [ val(meta), [ path(annotated_gtf) ] ]
        stats_gtf             = GFFCOMPARE.out.stats                       // channel: [ val(meta), [ path(stats) ] ]
        processed_gene_counts = DROP_COUNTS.out.processed_gene_counts      // channel: [ path(tsv) ]
        annotation_drop       = DROP_SAMPLE_ANNOT.out.drop_annot           // channel: [ path(sample_annotation.tsv) ]
        config_drop_ae        = DROP_CONFIG_RUN_AE.out.config_drop         // channel: [ path(confg_file.yml) ]
        drop_ae_out           = DROP_CONFIG_RUN_AE.out.drop_ae_out
        config_drop_as        = DROP_CONFIG_RUN_AS.out.config_drop         // channel: [ path(confg_file.yml) ]
        drop_as_out           = DROP_CONFIG_RUN_AS.out.drop_as_out
        drop_filter_ae_res    = DROP_FILTER_RESULTS.out.ae_out_filtered
        drop_unfilter_ae_res  = DROP_FILTER_RESULTS.out.ae_out_unfiltered
        drop_filter_as_res    = DROP_FILTER_RESULTS.out.as_out_filtered
        drop_unfilter_as_res  = DROP_FILTER_RESULTS.out.as_out_unfiltered
        versions              = ch_versions                                // channel: [ path(versions.yml) ]
}
