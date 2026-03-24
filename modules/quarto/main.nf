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
     export PATH=\$PATH:~/bin
     export PATH=\$PATH:/usr/bin
     export R_HOME=/usr/bin/R
     export R_LIBS=/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library
     export PATH=/usr/local/bin:/usr/bin:/home/fa286/bin:\$PATH
    ~/bin/quarto render ${qmd} --output gc_bias_report.html
    
    """
}
