import pandas as pd

configfile: "config.yaml"

sampleTable = pd.read_csv("samples.tsv", index_col="sample", sep="\t")
SAMPLES = sampleTable.index.values



def get_genome(wildcards):
    return sampleTable.at[wildcards.sample, 'genome_file']

def get_taxonomy(wildcards):
    return sampleTable.at[wildcards.sample, 'taxonomy']

def get_biomass(wildcards):
    return sampleTable.at[wildcards.sample, 'biomass']

def get_medium(wildcards):
    return sampleTable.at[wildcards.sample, 'medium']


localrules:
    gapseq_medium,
    install_gapseq,
    all


rule all:
    input:
        expand("models/{sample}/{sample}.RDS", sample=SAMPLES),
    output:
        touch("finished_recon")


rule install_gapseq:
    output:
        testlog="gapseq_test.txt",
        gapseq_bin="gapseq/gapseq"
    params:
        repo=config.get("gapseq_repo", 1),
        branch=config.get("gapseq_branch", 1),
        seqDB=config.get("gapseq_seqDB", 1),
        seqDBversion=config.get("gapseq_seqDBversion", 1)
    threads: 1
    log:
        "logs/install_gapseq.log"
    shell:
        """
        git clone --depth 1 -b {params.branch} {params.repo} > {log}
        cd gapseq
        ./gapseq update-sequences -t Bacteria -D {params.seqDB} -Z {params.seqDBversion} >> ../{log}
        ./gapseq update-sequences -t Archaea -D {params.seqDB} -Z {params.seqDBversion} >> ../{log}
        ./gapseq find -p all -t Bacteria -x -D {params.seqDB} toy/myb71.faa.gz >> ../{log}
        ./gapseq find -p all -t Archaea -x -D {params.seqDB} toy/myb71.faa.gz >> ../{log}
        ./gapseq test > ../{output.testlog}
        """


rule gapseq_find:
    input:
        genome=get_genome,
        testlog=rules.install_gapseq.output.testlog
    params:
        b=config.get("find_b", 1),
        taxonomy=get_taxonomy,
        aligner=config.get("aligner", 1),
        seqDB=config.get("gapseq_seqDB", 1),
        seqDBversion=config.get("gapseq_seqDBversion", 1)
    output:
        rxn="models/{sample}/{sample}-all-Reactions.tbl.gz",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz",
        trsp="models/{sample}/{sample}-Transporter.tbl.gz"
    threads: config.get("find_threads", 1)
    resources:
        mem_mb=config.get("find_mem", 1) * 1000,
        runtime=config.get("find_time", 1)
    log:
        "logs/find/{sample}.log"
    shell:
        """
        # Reactions / Pathways
        gapseq/./gapseq find -p all -b {params.b} -t {params.taxonomy} -m {params.taxonomy} -K {threads} -O -A {params.aligner} -f models/{wildcards.sample} -D {params.seqDB} {input.genome} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-all-Reactions.tbl
        gzip -f models/{wildcards.sample}/{wildcards.sample}-all-Pathways.tbl
        
        # Transporters
        gapseq/./gapseq find-transport -b {params.b} -K {threads} -A {params.aligner} -f models/{wildcards.sample} {input.genome} >> {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-Transporter.tbl
        """

        
rule gapseq_draft:
    input:
        rxn="models/{sample}/{sample}-all-Reactions.tbl.gz",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz",
        trsp="models/{sample}/{sample}-Transporter.tbl.gz"
    params:
        biomass=get_biomass,
        u=config.get("draft_u", 1),
        l=config.get("draft_l", 1)
    output:
        draft="models/{sample}/{sample}-draft.RDS",
        xml="models/{sample}/{sample}-draft.xml.gz"
    threads: config.get("draft_threads", 1)
    resources:
        mem_mb=config.get("draft_mem", 1) * 1000,
        runtime=config.get("draft_time", 1)
    log:
        "logs/draft/{sample}.log"
    shell:
        """
        gapseq/./gapseq draft -r {input.rxn} -t {input.trsp} -b {params.biomass} -p {input.pwy} -u {params.u} -l {params.l} -f models/{wildcards.sample} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-draft.xml
        """
        
rule gapseq_medium:
    input:
        model="models/{sample}/{sample}-draft.RDS",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz"
    params:
        c=config.get("medium_c", "")
    output:
        "models/{sample}/{sample}-medium.csv"
    log:
        "logs/medium/{sample}.log"
    shell:
        """
        par_c="{params.c}"
        if [ -n "$par_c" ]; then
            gapseq/./gapseq medium -m {input.model} -p {input.pwy} -c {params.c} -f models/{wildcards.sample} > {log}
        else
            gapseq/./gapseq medium -m {input.model} -p {input.pwy} -f models/{wildcards.sample} > {log}
        fi
        """
        
rule gapseq_fill:
    input:
        draft="models/{sample}/{sample}-draft.RDS",
        medium=get_medium
    params:
        b=config.get("fill_b", 1),
        mingr=config.get("fill_mingr", 1)
    output:
        model="models/{sample}/{sample}.RDS",
        xml="models/{sample}/{sample}.xml.gz",
        cs="models/{sample}/{sample}-cs.tbl",
        ferm="models/{sample}/{sample}-ferm.tbl"
    threads: config.get("fill_threads", 1)
    resources:
        mem_mb=config.get("fill_mem", 1) * 1000,
        runtime=config.get("fill_time", 1)
    log:
        "logs/fill/{sample}.log"
    shell:
        """
        if grep -q cpd11640 "{input.medium}"; then
            gapseq/./gapseq fill -m {input.draft} -n {input.medium} -b {params.b} -e highH2 -f models/{wildcards.sample} -k {params.mingr} -w > {log}
        else
             gapseq/./gapseq fill -m {input.draft} -n {input.medium} -b {params.b} -f models/{wildcards.sample} -k {params.mingr} -w > {log}
        fi
        
        gzip -f models/{wildcards.sample}/{wildcards.sample}.xml
        """

