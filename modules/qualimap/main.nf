
process QUALIMAP {

    maxForks 2
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
            --java-mem-size=40G \\
            -outdir ${id}.bamqc
        
        sed -i "s/bam file = ${id}.DEDUP.bam/bam file = ${id}.bam/g" ${id}.bamqc/genome_results.txt

        """

}
