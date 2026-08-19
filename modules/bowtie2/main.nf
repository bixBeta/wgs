
process BOWTIE2 {

    maxForks 3
    tag "$id, $genomePrefix"
    label "process_high"

    publishDir "primary_BAMS",          mode: "symlink", overwrite: true, pattern: "*.bam*"
    publishDir "STATS/SAMTOOLS",        mode: "symlink", overwrite: true, pattern: "*stat*"
    publishDir "STATS/BT2",             mode: "symlink", overwrite: true, pattern: "*log"

    input:
        tuple val(id), path(trimmed)
        val genome
        val genomePrefix


    output:
        tuple val(id), path("*.primary.sorted.bam")             , emit: "primary_sorted_bam"
        tuple val(id), path("*.primary.sorted.bam.bai")         , emit: "primary_sorted_bai"
        path("*.primary.log")                                   , emit: "primary_log"
        path("*.primary.flagstat")                              , emit: "primary_flagstat"
        path("*.primary.idxstats")                              , emit: "primary_idxstats"
        path "bowtie2_versions.yml"                             , emit: "versions"

    script:

        """
            (bowtie2 \\
            --no-unal \\
              -x  ${genome}${genomePrefix} \\
              -1 ${trimmed[0]} -2 ${trimmed[1]} \\
              --threads 24 \\
              --rg-id ${id} \\
              --rg SM:${id} \\
              --rg PL:ILLUMINA \\
              --rg LB:${id} \\
              -S - | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${id}.primary.bam)2>${id}.primary.log


        samtools sort ${id}.primary.bam > ${id}.primary.sorted.bam
        samtools index ${id}.primary.sorted.bam 

        samtools flagstat ${id}.primary.sorted.bam > ${id}.primary.flagstat
        samtools idxstats ${id}.primary.sorted.bam > ${id}.primary.idxstats

        cat <<-END_VERSIONS > bowtie2_versions.yml
        "BOWTIE2":
            bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/.*version //')
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS

        """


}
