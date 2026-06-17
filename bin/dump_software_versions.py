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
                    if isinstance(tools, dict):
                        for tool, version in tools.items():
                            versions[tool] = version  # deduplicate by tool name
    except Exception:
        pass

# Build HTML table rows
rows = ""
for tool in sorted(versions.keys()):
    version = versions[tool]
    rows += (
        "            <tr>\n"
        "                <td><samp>{}</samp></td>\n"
        "                <td><samp>{}</samp></td>\n"
        "            </tr>\n"
    ).format(tool, version)

table = (
    "    <table class=\"table\" style=\"width:100%\">\n"
    "        <thead>\n"
    "            <tr>\n"
    "                <th>Software</th>\n"
    "                <th>Version</th>\n"
    "            </tr>\n"
    "        </thead>\n"
    "        <tbody>\n"
    + rows +
    "        </tbody>\n"
    "    </table>\n"
)

content = (
    "id: 'software_versions'\n"
    "section_name: 'Software Versions'\n"
    "section_href: 'https://github.com/bixBeta/wgs'\n"
    "plot_type: 'html'\n"
    "description: 'Versions of software tools used in this pipeline run.'\n"
    "data: |\n"
    + table
)

with open("software_versions_mqc.yml", "w") as fh:
    fh.write(content)
