#!/usr/bin/env python3

import yaml
import glob

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

rows = ""
for process in sorted(versions.keys()):
    tools = versions[process]
    if isinstance(tools, dict):
        for tool, version in sorted(tools.items()):
            rows += "        <dt>{}</dt><dd><samp>{}</samp></dd>\n".format(tool, version)

content = (
    "id: 'software_versions'\n"
    "section_name: 'Software Versions'\n"
    "section_href: 'https://github.com/bixBeta/wgs'\n"
    "plot_type: 'html'\n"
    "description: 'Collected at run time from software output.'\n"
    "data: |\n"
    "    <dl class=\"dl-horizontal\">\n"
    + rows +
    "    </dl>\n"
)

with open("software_versions_mqc.yml", "w") as fh:
    fh.write(content)
