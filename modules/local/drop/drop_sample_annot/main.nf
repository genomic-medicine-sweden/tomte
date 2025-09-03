process DROP_SAMPLE_ANNOT {
    tag "DROP_annot_file"
    label 'process_low'

    container "docker.io/clinicalgenomics/drop:1.4.0"

    input:
    tuple val(meta2), path(gtf)
    tuple val(ids), val(single_ends), val(strandednesses), val(sex), path(bam), path(bai)
    path(ref_gene_counts)
    path(ref_annot)
    val(drop_group_samples_ae)
    val(drop_group_samples_as)

    output:
    path('sample_annotation.tsv'), emit: drop_annot
    path "versions.yml"   , emit: versions

    when:
    // Exit if running this module with -profile conda / -profile mamba
    workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0 ||
        { log.error("Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."); return false }

    script:
    def id = ids.join(' ')
    def single_end = single_ends.join(' ')
    def sex_drop = sex.collect { it.replace("1","M").replace("2","F").replace("0","NA").replace("other","NA") }.join(' ')
    def strandedness = strandednesses.join(' ')
    def drop_group = "${drop_group_samples_ae},${drop_group_samples_as}".replace(" ","").replace("[","").replace("]","")
    def reference_count_file = ref_gene_counts ? "--ref_count_file ${ref_gene_counts}" : ''
    def reference_annotation = ref_annot ? "--ref_annot ${ref_annot}" : ''
    """
    # Single_end values are only provided if you start from fastq
    SINGLE_ENDS=(${single_end})
    BAMS=(${bam.join(' ')})

    # Check if single_end values are provided
    updated_single_ends=()
    for ((i=0; i<\${#SINGLE_ENDS[@]}; i++)); do
        if [[ "\${SINGLE_ENDS[i]}" == "null" ]]; then
            result=\$(samtools view -c -f 1 "\${BAMS[i]}" | awk '{print \$1 == 0 ? "true" : "false"}')
            updated_single_ends+=("\$result")
        else
            updated_single_ends+=("\${SINGLE_ENDS[i]}")
        fi
    done

    # Convert updated_single_ends array to space-separated string and save to file
    echo "\${updated_single_ends[*]}" > updated_single_ends.txt

    drop_sample_annot.py \\
        --bam ${bam.join(' ')} \\
        --samples $id \\
        --strandedness $strandedness \\
        --single_end \$(cat updated_single_ends.txt) \\
        --sex $sex_drop \\
        $reference_count_file \\
        $reference_annotation \\
        --drop_group_sample $drop_group \\
        --gtf ${gtf} \\
        --output sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: \$(drop_sample_annot.py --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: \$(drop_sample_annot.py --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
