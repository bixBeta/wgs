process GC_BIAS_REPORT {

    label "process_quarto"
    publishDir "Reports", mode: 'copy'

    input:
    path pngs
    path qmd

    output:
    path "gc_bias_report.html"

    script:
    """
    quarto render ${qmd} --output gc_bias_report.html
    """
}
