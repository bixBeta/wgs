mqcgenome =  params.genome 

process MQC {

    label 'process_mqc'
    
    publishDir "Reports/${mqcgenome}", mode: "symlink", overwrite: true
    
    input:

        path "*"
        path(conf)
        path(logo)             

    output:
        path "*html"                    , emit: mqc_out  

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome} 
       multiqc -n ${params.id}_bt2_TESTTTTt --config ${conf} --cl-config "custom_logo: ${logo}" .

    """

}