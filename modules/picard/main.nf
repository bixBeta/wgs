
process MARKDUPS {

    maxForks 8
    tag "$id"
    label "process_high"

   
    publishDir "STATS/SAMTOOLS",             mode: "symlink", overwrite: true, pattern: "*stat*"
    publishDir "STATS/PICARD/",              mode: "symlink", overwrite: true, pattern: "*MarkDuplicates.metrics.txt"
    publishDir "DEDUP_BAMS",                 mode: "symlink", overwrite: true, pattern: "*DEDUP.bam*"

    input:
        tuple val(id), path(primary_bam)
        tuple val(id), path(primary_bai)

    output:

        tuple val(id), path("*.dupMarked.bam")                      , emit: "dupmarked_bam"
        tuple val(id), path("*.dupMarked.bam.bai")                  , emit: "dupmarked_bai"
        path("*.dupMarked.flagstat")                                , emit: "dupmarked_flagstat"
        path("*.dupMarked.idxstats")                                , emit: "dupmarked_idxstats"        

        path("*.MarkDuplicates.metrics.txt")                        , emit: "dup_stats"

        tuple val(id), path("*DEDUP.bam")                            , emit: "dedup_bam"
        tuple val(id), path("*DEDUP.bam.bai")                        , emit: "dedup_bai"
        path("*.DEDUP.flagstat")                                     , emit: "dedup_flagstat"
        path("*.DEDUP.idxstats")                                     , emit: "dedup_idxstats"

    script:

        """
            java -jar /myBin/picard.jar \\
                    MarkDuplicates \\
                    INPUT=${primary_bam} \\
                    OUTPUT=${id}.dupMarked.bam \\
                    ASSUME_SORTED=true \\
                    REMOVE_DUPLICATES=false \\
                    METRICS_FILE=${id}.MarkDuplicates.metrics.txt \\
                    VALIDATION_STRINGENCY=LENIENT \\
                    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
                    TMP_DIR=tmp

            
            

            samtools index ${id}.dupMarked.bam
            samtools flagstat ${id}.dupMarked.bam> ${id}.dupMarked.flagstat
            samtools idxstats ${id}.dupMarked.bam > ${id}.dupMarked.idxstats


            samtools view -b -h -F 0X400 ${id}.dupMarked.bam > ${id}.DEDUP.bam

            samtools index ${id}.DEDUP.bam
            samtools flagstat ${id}.DEDUP.bam > ${id}.DEDUP.flagstat
            samtools idxstats ${id}.DEDUP.bam > ${id}.DEDUP.idxstats

        """

}