nextflow.enable.dsl=2

// include a minimal example module
include { HELLO } from './modules/hello/hello'

// Default parameters
params.samples = "samples.csv"
params.outputDir = "results"
params.help = false

workflow {
    if (params.help) {
        println "Usage: nextflow run main.nf --samples samples.csv"
        exit 0
    }

    // simple samples channel: read CSV with header `name` or fallback to built-in example
    def samples_ch = Channel
        .fromPath(params.samples)
        .ifEmpty { Channel.of([['name':'example']]) }
        .splitCsv(header:true)
        .map { row -> row.name }

    // call the HELLO module for each sample name
    hello_out = HELLO(samples_ch)

    // publish greetings
    hello_out.view()
}