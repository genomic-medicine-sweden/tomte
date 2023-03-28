process RNA_DOWNSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(bam)
    path downsample_bed

    output:
    tuple val(meta), path("*_downsmapled.bam")                                , emit: bam
    tuple val(meta), path("*_downsmapled.bam.bai")                            , emit: bai
    tuple val(meta), path("*_downsmapled.bam"), path("*_downsmapled.bam.bai") , emit: bam_bai
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seed_frac       = (args.contains('-s')) ? '' : '-s 0.001'
    def num_reads       = (args.contains('-num_reads')) ? '' : '80000000'

    """
    if [[ \$mode == *"subsample_region"* ]]; then
        samtools view -@ $task.cpus -b -U non_select.bam -L ${downsample_bed} ${bam} | samtools view $seed_frac -@ $task.cpus -b -o select.bam
        if [[ \$mode == *"downsample"* ]]; then
            samtools index non_select.bam
            fraction=\$(samtools idxstats non_select.bam | cut -f3 | awk -v ct="$num_reads" 'BEGIN {total=0} {total += \$1} END {print ct/total}')
            percent=\$(samtools idxstats non_select.bam | cut -f3 | awk -v ct="$num_reads"00 'BEGIN {total=0} {total += \$1} END {print ct/total}')
            percent=\$(echo \$percent | cut -d. -f1)

            if (( \$(echo "\$percent < 100") )); then
                samtools merge -u non_select.bam select.bam -o /dev/stdout | samtools view -@ $task.cpus -b -s \$fraction -o ${prefix}_downsmapled.bam
            else
                samtools merge -u non_select.bam select.bam -o ${prefix}_downsmapled.bam
            fi
        fi
    else
        if [[ \$mode == *"downsample"* ]]; then
            fraction=\$(samtools idxstats ${bam} | cut -f3 | awk -v ct="$num_reads" 'BEGIN {total=0} {total += \$1} END {print ct/total}')
            percent=\$(samtools idxstats ${bam} | cut -f3 | awk -v ct="$num_reads"00 'BEGIN {total=0} {total += \$1} END {print ct/total}')
            percent=\$(echo \$percent | cut -d. -f1)
            
            if (( \$(echo "\$percent < 100") )); then
                samtools view -@ $task.cpus -b ${bam} -s \$fraction -o ${prefix}_downsmapled.bam
            else
                mv ${bam} ${prefix}_downsmapled.bam
            fi
        else
            mv ${bam} ${prefix}_downsmapled.bam
        fi
    fi
    
    samtools index ${prefix}_downsmapled.bam

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
