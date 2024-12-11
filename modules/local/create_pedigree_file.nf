process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    input:
    val samples

    output:
    path ("peddy.ped"), emit: ped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    samples.each { sample ->
        def sex_int = sample.sex == "M" ? "1" : sample.sex == "F" ? "2" : "0"
        outfile_text += "\\n" + [sample.case, sample.sample, '0', '0', sex_int, '0'].join('\\t')
    }
    """
    echo -e "${outfile_text}" > peddy.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch peddy.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
    END_VERSIONS
    """
}
