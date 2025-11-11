// Minimal example module (DSL2)

process HELLO_PROCESS {
    tag { sampleName }
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy'

    input:
    val sampleName

    output:
    path "${sampleName}_greeting.txt"

    script:
    """
    echo "Hello from Nextflow, ${sampleName}!" > ${sampleName}_greeting.txt
    cat ${sampleName}_greeting.txt
    """
}

// module wrapper to accept a channel
workflow HELLO (sample_ch) {
    sample_ch
        .map { it -> it }
        .set { names }

    names.into { in_ch }

    HELLO_PROCESS(in_ch)
}