process GCBIAS {

    maxForks 3
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2
    tag "$id"
    label "process_deeptools"

    publishDir "GCBias_DeepTools" ,            mode: "symlink", overwrite: true, pattern: "*.png"
    publishDir "GCBias_DeepTools" ,            mode: "symlink", overwrite: true, pattern: "*.txt"


    input:

        tuple val(id), path(dedup_bam)
        tuple val(id), path(dedup_bai)
        path twoBits
        val egsize

    output:

        tuple val(id), path("*.png")                 , emit: "gc_png"
        tuple val(id), path("*.txt")                 , emit: "gc_txt"
        path "deeptools_versions.yml"                , emit: "versions"

    script:

        """
            mkdir -p tmp/
            export MPLCONFIGDIR="tmp/"

            computeGCBias -b ${dedup_bam} \\
            --effectiveGenomeSize ${egsize} \\
            -g ${twoBits} -l 200 \\
            --GCbiasFrequenciesFile ${id}_gcBias_freq.txt \\
            --biasPlot ${id}_gc.png

            cat <<-END_VERSIONS > deeptools_versions.yml
            "GCBIAS":
                deeptools: \$(computeGCBias --version 2>&1 | head -1 | sed 's/computeGCBias //')
            END_VERSIONS

        """



}
