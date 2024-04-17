//
// Annotating SNVs
//

include { ENSEMBLVEP_VEP       } from '../../modules/nf-core/ensemblvep/vep/main'
include { RENAME_FILES         } from '../../modules/local/rename_files'
include { TABIX_BGZIPTABIX     } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX          } from '../../modules/nf-core/tabix/tabix/main'
include { GAWK                 } from '../../modules/nf-core/gawk/main'
include { CREATE_HGNCIDS_FILE  } from '../../modules/local/create_hgncids_file.nf'
include { ENSEMBLVEP_FILTERVEP } from '../../modules/nf-core/ensemblvep/filtervep/main'



workflow ANNOTATE_SNV {
    take:
    vcf                            // channel:   [mandatory] [ val(meta), path(vcf), path(tbi) ]
    val_vep_genome                 // parameter: [mandatory] 'GRCh37' or 'GRCh38'
    val_vep_cache_version          // parameter: [mandatory] default: 110
    ch_vep_cache                   // channel:   [mandatory] [ path(cache) ]
    ch_fasta                       // channel:   [mandatory] [ val(meta), path(fasta) ]
    ch_vep_extra_files             // channel:   [mandatory] [ path(files) ]
    ch_gene_panel_clinical_filter  // channel:   [optional]  [ path(file) ]

    main:
        ch_versions = Channel.empty()

        // Annotate with VEP
        ENSEMBLVEP_VEP(
            vcf.map{ meta, vcf -> [ meta, vcf, [] ] },
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_fasta,
            ch_vep_extra_files
        )

        ENSEMBLVEP_VEP.out.vcf
            .multiMap { meta, vcf ->
                clinical: [ meta + [ set: "clinical" ], vcf ]
                research: [ meta + [ set: "research" ], vcf ]
            }
            .set { ch_clin_research_vcf }

        RENAME_FILES( ch_clin_research_vcf.research )

        TABIX_TABIX( RENAME_FILES.out.output )
        ch_vcf_research = RENAME_FILES.out.output.join(TABIX_TABIX.out.tbi)

        GAWK( ch_gene_panel_clinical_filter.map{it -> [[id:'hgnc'], it]}.collect(),
            [] )

        // Generate Clinical filter
        CREATE_HGNCIDS_FILE( GAWK.out.output )
            .txt
            .set {ch_hgnc_ids}

        //Filter results
        ENSEMBLVEP_FILTERVEP(
            ch_clin_research_vcf.clinical,
            ch_hgnc_ids
        )
        .output
        .set { ch_filtervep_out }

        TABIX_BGZIPTABIX( ch_filtervep_out )
        ch_vcf_clin = TABIX_BGZIPTABIX.out.gz_tbi

        ch_versions = ch_versions.mix( ENSEMBLVEP_VEP.out.versions.first() )
        ch_versions = ch_versions.mix( GAWK.out.versions )
        ch_versions = ch_versions.mix( ENSEMBLVEP_FILTERVEP.out.versions )
        ch_versions = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )
        ch_versions = ch_versions.mix( TABIX_TABIX.out.versions )

    emit:
        //vcf_gz   = ENSEMBLVEP_VEP.out.vcf    // channel: [ val(meta), path(vcf.gz) ]
        //tbi_gz   = TABIX_VEP.out.tbi         // channel: [ val(meta), path(tbi) ]
        tab_gz   = ENSEMBLVEP_VEP.out.tab    // channel: [ val(meta), path(tab.gz) ]
        json_gz  = ENSEMBLVEP_VEP.out.json   // channel: [ val(meta), path(json.gz) ]
        report   = ENSEMBLVEP_VEP.out.report // channel: [ path(html) ]
        ch_vcf_clin = ch_vcf_clin // channel: [ val(meta), path(vcf.gz) path(tbi)]
        ch_vcf_research = ch_vcf_research // channel: [ val(meta), path(vcf.gz) path(tbi)]
        versions = ch_versions               // channel: [ path(versions.yml) ]
}

