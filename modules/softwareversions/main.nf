process SOFTWAREVERSIONS {

    label 'process_mqc'

    input:
        path versions

    output:
        path "software_versions_mqc.yml"    , emit: mqc_yml

    script:
    """
    #!/usr/bin/env python3

    import yaml
    import glob

    # Collect and deduplicate versions from all input yml files
    versions = {}
    for fname in glob.glob("*.yml"):
        try:
            with open(fname) as fh:
                data = yaml.safe_load(fh)
                if data and isinstance(data, dict):
                    for process, tools in data.items():
                        if process not in versions and isinstance(tools, dict):
                            versions[process] = tools
        except Exception:
            pass

    # Build MultiQC custom content HTML table
    rows = ""
    for process in sorted(versions.keys()):
        tools = versions[process]
        if isinstance(tools, dict):
            for tool, version in sorted(tools.items()):
                rows += f"        <dt>{tool}</dt><dd><samp>{version}</samp></dd>\\n"

    mqc_yaml = f"""id: 'software_versions'
section_name: 'Software Versions'
section_href: 'https://github.com/bixBeta/wgs'
plot_type: 'html'
description: 'Collected at run time from software output.'
data: |
    <dl class="dl-horizontal">
{rows}    </dl>
"""

    with open("software_versions_mqc.yml", "w") as fh:
        fh.write(mqc_yaml)
    """

}
