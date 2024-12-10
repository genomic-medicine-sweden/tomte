process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'biocontainers/python:3.8.3'}"

    input:
    val samples

    output:
    path ("*.ped"), emit: ped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def case_name = samples[0].case
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    samples
        .unique { it.sample }
        .each { sample ->
            outfile_text += "\\n" + [sample.case, sample.sample, 0, 0, sample.sex, "affected"].join('\\t')
        }
    """
    echo -e "${outfile_text}" >${case_name}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def case_name = samples[0].case
    """
    touch ${case_name}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
