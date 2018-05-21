configfile: "config.json"

rule all:
    input: ""
    params: sge_opts=""

rule some_rule:
    input: ""
    output: ""
    params: sge_opts=""
    shell: ""
