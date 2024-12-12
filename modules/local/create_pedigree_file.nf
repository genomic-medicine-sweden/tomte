process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    val(samples)

    output:
    path("*.ped"),       emit: ped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def case_name = samples[0].case
    def out   = new File(case_name + ".ped")
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        def sample_name =  samples[i].sample
        if (!samples_list.contains(sample_name)) {
            def sex_int = sample.sex == "M" ? "1" : sample.sex == "F" ? "2" : "0"
            outfile_text += "\\n" + [samples[i].case, sample_name, samples[i].paternal, samples[i].maternal, sex_int, samples[i].phenotype].join('\\t')
            samples_list.add(sample_name)
        }
    }
    """
    echo -e "$outfile_text" >${case_name}.ped

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


// process CREATE_PEDIGREE_FILE {
//     tag "pedigree"
//     label 'process_single'

//     input:
//     val samples

//     output:
//     path ("peddy.ped"), emit: ped
//     path "versions.yml", emit: versions

//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
//     samples.each { sample ->
//         def sex_int = sample.sex == "M" ? "1" : sample.sex == "F" ? "2" : "0"
//         outfile_text += "\\n" + [sample.case, sample.sample, '0', '0', sex_int, '0'].join('\\t')
//     }
//     """
//     echo -e "${outfile_text}" > peddy.ped

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         create_pedigree_file: v1.0
//     END_VERSIONS
//     """

//     stub:
//     """
//     touch peddy.ped

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         create_pedigree_file: v1.0
//     END_VERSIONS
//     """
// }
