nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"
params.help             = false
params.fastqs           = false
params.id               = "TREx_ID"
params.listGenomes      = false

// Module Params:
params.bowtie2          = false
params.fastp            = false
params.mode             = "PE"
params.genome           = null
params.bamqc            = false


// Base Channels:

ch_pin                  = channel.value(params.id)
ch_sheet                = channel.fromPath(params.sheet)
ch_mqc_conf             = channel.fromPath("${projectDir}/multiqc_config.yaml")
ch_mqc_logo             = channel.fromPath("${projectDir}/img/trex-mini.png")


if( params.help ) {

log.info """
W  G  S  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/wgs -r main < args ... >

Args:
    * --help           : Prints this help documentation
    * --listGenomes    : Get extended list of genomes available for this pipeline
    * --id             : TREx Project ID 
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -------------------------------------------
        Sample Sheet Example:    
        label   fastq1          fastq2
        SS1     SS1_R1.fastq.gz SS1_R2.fastq.gz
        SS2     SS2_R1.fastq.gz SS2_R2.fastq.gz  
        .
        .
        . etc.
        -------------------------------------------
    * --bwa         : deprecated
    * --bowtie2     : invokes bowtie2 alignment step
    * --fastp       : invokes fastp trimming
    * --mode        : PE (default)
    * --genome      : reference genome ( Available genomes: hg38, mm10, fc9 & dm6 : Use --listGenomes for more details )

"""

    exit 0
}



// BT2 Indices MAP

genomeDir = [

mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/bowtie2index",
hg38            :"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/BT2.ENSEMBL_INDEX",
dm6             :"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Bowtie2.Index",
canFam4         :"/workdir/genomes/Canis_familiaris/canFam4/NCBI/bowtie2",
fc9             :"/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/bowtie2"

]

bt2Prefix = [

mm10            :"/mm10",
hg38            :"/hg38",
dm6             :"/dm6",
canFam4         :"/cf4",
fc9             :"/fc9"

]

// 2Bit MAP  effective genome Size for deeptools

twoBits = [

mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.dna.2bit",
hg38            :"nA",
dm6             :"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.32.dna.2bit",
canFam4         :"/workdir/genomes/Canis_familiaris/canFam4/NCBI/canFam4.2bit",
fc9             :"/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.dna.2bit"

]


gSize = [

mm10            :2652783500,
hg38            :2913022398,
dm6             :142573017,
canFam4         :2482000080,
fc9             :2521863845

]


if( params.listGenomes) {
    
    println("")
    log.info """
    Available BT2 Indices
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    genomeDir.each(printMap)

    exit 0
}


// Load all modules 

include {    FASTP     } from './modules/fastp'
include {    BOWTIE2   } from './modules/bowtie2'
include {    MARKDUPS  } from './modules/picard'
include {    QUALIMAP  } from './modules/qualimap'
include {    MQC       } from './modules/multiqc'
include {    GCBIAS    } from './modules/deeptools'

workflow BTPAIRED {

    // Trimming 

    if( params.fastqs ){
        meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("fastqs/" + row.fastq1), file("fastqs/" + row.fastq2)]] }
                |  view 
    } else {

        meta_ch = ch_sheet
            |  splitCsv( header:true )
            |  map { row -> [row.label, [file(row.fastq1), file(row.fastq2)]] }
            |  view
    }

    if( params.fastp ){

        FASTP(meta_ch)
            .set{fastp_out}
    } else {

        fastp_out = meta_ch
    }

    // Alignment 

    genome = genomeDir[params.genome]
    ch_genome   = channel.value(genome)

    genome_prefix = bt2Prefix[params.genome]
    ch_genome_prefix = channel.value(genome_prefix)

    if( params.bowtie2 ){
        
        if (params.fastp){
            BOWTIE2(fastp_out.trimmed_fqs, ch_genome, ch_genome_prefix)
        } else {
            BOWTIE2(fastp_out, ch_genome, ch_genome_prefix)
        }    
        
        bt2_ch  = BOWTIE2.out.primary_sorted_bam
        bt2i_ch = BOWTIE2.out.primary_sorted_bai

        MARKDUPS(bt2_ch, bt2i_ch)
        // qc_ch = MARKDUPS.out.dupmarked_bam
        qc_ch   = MARKDUPS.out.dedup_bam
        qc_ch_i = MARKDUPS.out.dedup_bai

        QUALIMAP(qc_ch)

        ch_gs   = channel.value(gSize[genome])
        ch_2bit = channel.value(twoBits[genome])


        GCBIAS(qc_ch, qc_ch_i, ch_2bit, ch_gs)

        mqc_ch = BOWTIE2.out.primary_log
                    .concat(
                        BOWTIE2.out.primary_flagstat, 
                        BOWTIE2.out.primary_idxstats,
                        QUALIMAP.out.bamqc_out,
                        MARKDUPS.out.dupmarked_flagstat,
                        MARKDUPS.out.dupmarked_idxstats,
                        MARKDUPS.out.dedup_flagstat,
                        MARKDUPS.out.dedup_idxstats,
                        MARKDUPS.out.dup_stats)
                    .collect()
                    .view()

        MQC(mqc_ch, ch_mqc_conf, ch_mqc_logo)

    }
}



workflow{
    if( params.mode == "PE") {
        BTPAIRED()
    } else {

        exit 0
    }
}
