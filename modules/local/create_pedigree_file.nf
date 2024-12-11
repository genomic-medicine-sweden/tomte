process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    input:
    val meta

    output:
    path ("*.ped"), emit: ped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    outfile_text += "\\n" + [meta.case, meta.sample, '0', '0', meta.sex, '0'].join('\\t')
    """
    echo -e "${outfile_text}" > ${meta.sample}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
    END_VERSIONS
    """
}
