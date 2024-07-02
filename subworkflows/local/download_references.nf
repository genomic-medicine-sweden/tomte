//
// Download reference files
//

include { GENCODE_DOWNLOAD as FASTA_DOWNLOAD } from '../../modules/local/gencode_download'
include { GENCODE_DOWNLOAD as GTF_DOWNLOAD   } from '../../modules/local/gencode_download'
include { VEP_DOWNLOAD                       } from '../../modules/local/vep_annotation_download'

workflow DOWNLOAD_REFERENCES {
    take:
        ch_genome                 // channel: [mandatory]   val(genome)
        ch_genome_version         // channel: [mandatory]   val(genome_version)
        ch_vep_cache_version      // channel: [optional]    val(vep_cache_version)

    main:
        ch_versions = Channel.empty()

        // Download fasta if not provided
        FASTA_DOWNLOAD(ch_genome, ch_genome_version, "fasta")

        // Download gtf if not provided
        GTF_DOWNLOAD(ch_genome, ch_genome_version, "gtf")

        // Download vep references if skip_download_vep = False & params.vep_cache is not provided
        VEP_DOWNLOAD(ch_genome, ch_vep_cache_version)

        ch_versions = ch_versions.mix(FASTA_DOWNLOAD.out.versions)
        ch_versions = ch_versions.mix(GTF_DOWNLOAD.out.versions)
        ch_versions = ch_versions.mix(VEP_DOWNLOAD.out.versions)

    emit:
        fasta      = FASTA_DOWNLOAD.out.fasta     // channel: [ path(fasta) ]
        gtf        = GTF_DOWNLOAD.out.gtf         // channel: [ path(gtf) ]
        vep_cache  = VEP_DOWNLOAD.out.vep_cache   // channel: [ path(vep_cache) ]
        vep_plugin = VEP_DOWNLOAD.out.plugin_file // channel: [ path(vep_plugin) ]
        versions   = ch_versions                  // channel: [ path(versions.yml) ]
}