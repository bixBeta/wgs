process SOFTWAREVERSIONS {

    label 'process_mqc'

    input:
        path versions

    output:
        path "software_versions_mqc.yml"    , emit: mqc_yml

    script:
    """
    dump_software_versions.py
    """

}
