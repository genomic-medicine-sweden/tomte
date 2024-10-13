//
// Download reference files
//

include { GENCODE_DOWNLOAD as FASTA_DOWNLOAD } from '../../modules/local/gencode_download'
include { GENCODE_DOWNLOAD as GTF_DOWNLOAD   } from '../../modules/local/gencode_download'
include { WGET_DOWNLOAD as WGET_DOWNLOAD     } from '../../modules/local/wget_download'
include { VEP_GNOMAD_DOWNLOAD                } from '../../modules/local/vep_gnomad_download'
include { BUILD_VEP_CACHE                    } from '../../modules/local/build_vep_cache'

workflow DOWNLOAD_REFERENCES {
    take:
        ch_genome                        // channel: [mandatory]  val(genome)
        ch_gencode_annotation_version    // channel: [mandatory]  val(gencode_annotation_version)
        ch_vep_refs_download_unprocessed // channel: [optional]   val(path_to_csv)
        ch_vep_cache_version             // channel: [optional]   val(vep_cache_version)

    main:
        ch_versions = Channel.empty()

        // Download fasta if not provided
        FASTA_DOWNLOAD(ch_genome, ch_gencode_annotation_version, "fasta")

        // Download gtf if not provided
        GTF_DOWNLOAD(ch_genome, ch_gencode_annotation_version, "gtf")

        // Read and store paths in vep_refs_download_unprocessed
        ch_vep_refs_download_unprocessed.splitCsv(header: true)
            .map { row -> return tuple(row.name, row.path_for_wget) }
            .set { ch_vep_refs_download }

        // Download files
        WGET_DOWNLOAD(ch_vep_refs_download.filter{ it != null })
        VEP_GNOMAD_DOWNLOAD(ch_genome, ch_vep_cache_version)

        BUILD_VEP_CACHE(WGET_DOWNLOAD.out.downloaded_file.collect(), VEP_GNOMAD_DOWNLOAD.out.gnomad_vcf_tbi.flatten().collect())

        ch_versions = ch_versions.mix(FASTA_DOWNLOAD.out.versions)
        ch_versions = ch_versions.mix(GTF_DOWNLOAD.out.versions)
        ch_versions = ch_versions.mix(WGET_DOWNLOAD.out.versions)
        ch_versions = ch_versions.mix(VEP_GNOMAD_DOWNLOAD.out.versions)


    emit:
        fasta      = FASTA_DOWNLOAD.out.fasta            // channel: [ path(fasta) ]
        gtf        = GTF_DOWNLOAD.out.gtf                // channel: [ path(gtf) ]
        vep_cache  = BUILD_VEP_CACHE.out.vep_cache       // channel: [ path(vep_cache) ]
        vep_plugin = BUILD_VEP_CACHE.out.vep_plugin_file // channel: [ path(vep_plugin) ]
        versions   = ch_versions                         // channel: [ path(versions.yml) ]
}
