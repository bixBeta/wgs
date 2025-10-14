
process QUALIMAP {

    maxForks 8
    tag "$id"
    label "process_qualimap"

   
    publishDir "STATS/QUALIMAP_RES",             mode: "symlink", overwrite: true, pattern: "*bamqc"
    

    input:
        tuple val(id), path(dup_marked_bam)

    output:
        path("*bamqc")                   , emit: bamqc_out
        
    script:

        """
            qualimap bamqc -bam ${dup_marked_bam} \\
            -nt 16 -c \\
            --java-mem-size=6G \\
            -outdir ${id}.bamqc

        """

}