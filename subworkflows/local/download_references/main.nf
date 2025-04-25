//
// Download reference files
//

include { GENCODE_DOWNLOAD as FASTA_DOWNLOAD } from '../../../modules/local/download/gencode_download'
include { GENCODE_DOWNLOAD as GTF_DOWNLOAD   } from '../../../modules/local/download/gencode_download'
include { WGET_DOWNLOAD                      } from '../../../modules/local/download/wget_download'
include { VEP_GNOMAD_DOWNLOAD                } from '../../../modules/local/download/vep_gnomad_download'
include { BUILD_VEP_CACHE                    } from '../../../modules/local/build_vep_cache'

workflow DOWNLOAD_REFERENCES {
    take:
        ch_genome                        // channel: [mandatory]  val(genome)
        ch_gencode_annotation_version    // channel: [mandatory]  val(gencode_annotation_version)
        ch_vep_refs_download_unprocessed // channel: [optional]   val(path_to_csv)
        ch_vep_cache_version             // channel: [optional]   val(vep_cache_version)

    main:
        ch_versions = Channel.empty()

        // Download fasta if not provided
        ch_downloaded_fasta = Channel.empty()
        if ( !params.fasta ) {
            FASTA_DOWNLOAD(ch_genome, ch_gencode_annotation_version, "fasta")
            ch_downloaded_fasta = ch_downloaded_fasta.mix( FASTA_DOWNLOAD.out.fasta ).collect()
            ch_versions = ch_versions.mix(FASTA_DOWNLOAD.out.versions)
        }

        // Download gtf if not provided
        ch_downloaded_gtf = Channel.empty()
        if ( !params.gtf ) {
            GTF_DOWNLOAD(ch_genome, ch_gencode_annotation_version, "gtf")
            ch_downloaded_gtf = ch_downloaded_gtf.mix( GTF_DOWNLOAD.out.gtf ).collect()
            ch_versions = ch_versions.mix(GTF_DOWNLOAD.out.versions)
        }

        // Read and store paths in vep_refs_download_unprocessed
        // Download files
        if ( !params.skip_download_vep && !params.vep_cache ) { 
            ch_vep_refs_download_unprocessed.splitCsv(header: true)
                .map { row -> return tuple(row.name, row.path_for_wget) }
                .set { ch_vep_refs_download }
            WGET_DOWNLOAD(ch_vep_refs_download.filter{ it != null })
            ch_versions = ch_versions.mix(WGET_DOWNLOAD.out.versions)
        }

        if ( !params.skip_download_gnomad ) {
            VEP_GNOMAD_DOWNLOAD(ch_genome, ch_vep_cache_version)
            ch_versions = ch_versions.mix(VEP_GNOMAD_DOWNLOAD.out.versions)
        }

        ch_built_vep_cache = Channel.empty()
        ch_built_vep_plugin_file = Channel.empty()
        if ( !params.skip_download_vep && !params.vep_cache ){
            BUILD_VEP_CACHE(WGET_DOWNLOAD.out.downloaded_file.collect(), VEP_GNOMAD_DOWNLOAD.out.gnomad_vcf_tbi.flatten().collect())
            ch_built_vep_cache = ch_built_vep_cache.mix( BUILD_VEP_CACHE.out.vep_cache )
            ch_built_vep_plugin_file = ch_built_vep_plugin_file.mix( BUILD_VEP_CACHE.out.vep_plugin_file )
            ch_versions = ch_versions.mix(BUILD_VEP_CACHE.out.versions)
        }


    emit:
        fasta      = ch_downloaded_fasta      // channel: [ path(fasta) ]
        gtf        = ch_downloaded_gtf        // channel: [ path(gtf) ]
        vep_cache  = ch_built_vep_cache       // channel: [ path(vep_cache) ]
        vep_plugin = ch_built_vep_plugin_file // channel: [ path(vep_plugin) ]
        versions   = ch_versions              // channel: [ path(versions.yml) ]
}
