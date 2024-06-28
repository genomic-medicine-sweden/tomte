//
// Download reference files
//

include { GENCODE_DOWNLOAD as DOWNLOAD_FASTA          } from '../../modules/local/gencode_download'
include { GENCODE_DOWNLOAD as DOWNLOAD_GTF            } from '../../modules/local/gencode_download'

workflow DOWNLOAD_REFERENCES {
    take:
        ch_genome                 // channel: [mandatory]   val(genome)
        ch_genome_version         // channel: [mandatory]   val(genome_version)

    main:
        ch_versions = Channel.empty()

        // Download fasta if not provided
        DOWNLOAD_FASTA(ch_genome, ch_genome_version, "fasta")

        // Download gtf if not provided
        DOWNLOAD_GTF(ch_genome, ch_genome_version, "gtf")

        ch_versions = ch_versions.mix(DOWNLOAD_FASTA.out.versions)
        ch_versions = ch_versions.mix(DOWNLOAD_GTF.out.versions)

    emit:
        fasta    = DOWNLOAD_FASTA.out.fasta // channel: [ path(fasta) ]
        gtf      = DOWNLOAD_GTF.out.gtf     // channel: [ path(gtf) ]
        versions = ch_versions              // channel: [ path(versions.yml) ]
}