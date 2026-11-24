# HLA class I typing, genome-free: no genome index, no HLA reads lost to pseudogenes,
# alt haplotypes or patches of full GENCODE/Ensembl genomes (APM_locus/README.md).
#   OptiType 1.3.5 [Szolek 2014]: razers3 fishing vs its IMGT reference -> ILP.
#   arcasHLA 0.6.0 [Orenbuch 2020]: genotype from FASTQ (kallisto 0.44 = its internal engine).
#   reconcile_hla.py: read support on fished reads, site check (MAFFT), Cellosaurus [Bairoch 2018].
#   select_hla_alleles.py -> IPD-IMGT/HLA CDS [Robinson 2019] -> build_personal_MANE_utr.py
#   (MANE UTRs kept, HLA-A/B/C CDS exchanged) [Aguiar 2019; Morales 2022].
# Typing uses the raw reads of every run; the group (cell line / individual) takes the calls of
# its run with most HLA reads; all runs are compared in <group>_hla_run_concordance.tsv.


rule optitype_fish:
    input: stage_reads("raw")
    output:
        r1=f"{OUT}/runs/{{run}}/hla/optitype/{{run}}_1.fished.fq",
        r2=f"{OUT}/runs/{{run}}/hla/optitype/{{run}}_2.fished.fq",
    threads: config["threads"]
    conda: "../../envs/optitype.yaml"
    shell:
        """
        REF=$(dirname $(readlink -f $(which OptiTypePipeline.py)))/data/hla_reference_rna.fasta
        for i in 1 2; do
            if [ $i = 1 ]; then f={input[0]}; o={output.r1}; else f={input[1]}; o={output.r2}; fi
            zcat $f > $o.tmp.fq
            razers3 -i 95 -m 1 -dr 0 -tc {threads} -o $o.bam $REF $o.tmp.fq
            samtools bam2fq $o.bam > $o; rm $o.tmp.fq $o.bam
        done
        """


rule optitype:
    input: rules.optitype_fish.output
    output: f"{OUT}/runs/{{run}}/hla/optitype/{{run}}_result.tsv"
    params: d=f"{OUT}/runs/{{run}}/hla/optitype/run"
    conda: "../../envs/optitype.yaml"
    log: f"{OUT}/logs/{{run}}_optitype.log"
    shell:
        "export PYTHONNOUSERSITE=1; rm -rf {params.d}; "
        "OptiTypePipeline.py -i {input[0]} {input[1]} --rna -v -o {params.d} > {log} 2>&1; "
        "cp {params.d}/*/*_result.tsv {output}"


rule arcashla:
    input: stage_reads("raw")
    output: f"{OUT}/runs/{{run}}/hla/arcashla/{{run}}.genotype.json"
    params: d=f"{OUT}/runs/{{run}}/hla/arcashla"
    threads: config["threads"]
    conda: "../../envs/arcashla.yaml"
    log: f"{OUT}/logs/{{run}}_arcashla.log"
    shell:
        # arcasHLA reference: built once by server/helpers/download_references.sh (IMGT 3.63.1);
        # default --min_count (75); the old min_count 5 gave junk calls
        "export PYTHONNOUSERSITE=1; "
        "arcasHLA genotype {input} -g A,B,C -t {threads} -o {params.d} --temp {params.d}/tmp > {log} 2>&1; "
        "mv {params.d}/{wildcards.run}_raw_1.genotype.json {output}"


rule reconcile_hla:
    input:
        opt=rules.optitype.output,
        arc=rules.arcashla.output,
        reads=rules.optitype_fish.output,
        nuc=config["imgt_hla_nuc"],
    output:
        calls=f"{OUT}/runs/{{run}}/hla/{{run}}_hla_calls.tsv",
        support=f"{OUT}/runs/{{run}}/hla/{{run}}_hla_support.tsv",
    params:
        cel=lambda wc: f"--cellosaurus {SAMPLES.loc[wc.run, 'cellosaurus']}" if SAMPLES.loc[wc.run, "cellosaurus"] else "",
        over=lambda wc: f"--override {SAMPLES.loc[wc.run, 'hla_override']}" if SAMPLES.loc[wc.run, "hla_override"] else "",
        d=f"{OUT}/runs/{{run}}/hla",
    threads: 4
    conda: "../../envs/hla_reconcile.yaml"
    log: f"{OUT}/logs/{{run}}_reconcile_hla.log"
    shell:
        "python3 {PIPE}/scripts/reconcile_hla.py --sample {wildcards.run} --optitype {input.opt} "
        "--arcas {input.arc} --reads {input.reads} --hla-nuc {input.nuc} {params.cel} {params.over} "
        "--threads {threads} --outdir {params.d} > {log} 2>&1"


rule group_hla:
    input: lambda wc: expand(f"{OUT}/runs/{{run}}/hla/{{run}}_hla_calls.tsv", run=runs_of(wc.group))
    output:
        calls=f"{OUT}/groups/{{group}}/hla/{{group}}_hla_calls.tsv",
        conc=f"{OUT}/groups/{{group}}/hla/{{group}}_hla_run_concordance.tsv",
    conda: "../../envs/py.yaml"
    shell:
        "python3 {PIPE}/scripts/group_hla_calls.py --group {wildcards.group} --calls {input} "
        "--out {output.calls} --concordance {output.conc}"


rule personal_mane:
    input:
        calls=rules.group_hla.output.calls,
        nuc=config["imgt_hla_nuc"],
        mane=config["mane_rna"], gff=config["mane_gff"], prot=config["mane_protein"],
        saf=f"{config['refs']}/MANE/cds_mane.saf",
    output:
        alleles=f"{OUT}/groups/{{group}}/ref/{{group}}_alleles.fa",
        report=f"{OUT}/groups/{{group}}/ref/{{group}}_alleles.tsv",
        fa=f"{OUT}/groups/{{group}}/ref/{{group}}_personal_MANE.fa",
        saf=f"{OUT}/groups/{{group}}/ref/{{group}}_cds_mane.saf",
        faa=f"{OUT}/groups/{{group}}/ref/{{group}}_personal_MANE_protein.faa",
    conda: "../../envs/py.yaml"
    log: f"{OUT}/logs/{{group}}_personal_mane.log"
    shell:
        "python3 {PIPE}/scripts/select_hla_alleles.py --calls {input.calls} --hla-nuc {input.nuc} "
        "  --out {output.alleles} --report {output.report} > {log} 2>&1; "
        "python3 {PIPE}/scripts/build_personal_MANE_utr.py --sample {wildcards.group} --mane {input.mane} "
        "  --gff {input.gff} --saf {input.saf} --alleles {output.alleles} --out {output.fa} "
        "  --out_saf {output.saf} >> {log} 2>&1; "
        "python3 {PIPE}/scripts/build_personal_MANE_protein.py --personal_mrna {output.fa} "
        "  --ref_protein {input.prot} --out {output.faa} >> {log} 2>&1"
