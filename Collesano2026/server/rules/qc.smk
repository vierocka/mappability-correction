# Reads: SRA download [Leinonen 2011] or local FASTQ; FastQC text parsed for mean read length
# (selects the funique table) and mean Phred (≥ 30 expected); Trimmomatic [Bolger 2014] optional.


def fastq_source(wc):
    r = SAMPLES.loc[wc.run]
    return [r.fastq_1, r.fastq_2] if r.source == "fastq" else []


rule fetch:
    input: fastq_source
    output:
        r1=f"{OUT}/runs/{{run}}/fastq/{{run}}_raw_1.fq.gz",
        r2=f"{OUT}/runs/{{run}}/fastq/{{run}}_raw_2.fq.gz",
    params:
        d=f"{OUT}/runs/{{run}}/fastq",
        src=lambda wc: SAMPLES.loc[wc.run, "source"],
        f1=lambda wc: SAMPLES.loc[wc.run, "fastq_1"],
        f2=lambda wc: SAMPLES.loc[wc.run, "fastq_2"],
    threads: 6
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/{{run}}_fetch.log"
    shell:
        """
        if [ {params.src} = fastq ]; then
            ln -sf $(readlink -f {params.f1}) {output.r1}; ln -sf $(readlink -f {params.f2}) {output.r2}
        else
            prefetch --max-size u {wildcards.run} -O {params.d} > {log} 2>&1
            fasterq-dump --split-files -e {threads} -O {params.d} {params.d}/{wildcards.run} >> {log} 2>&1
            # R1 and R2 must differ (the 2026 HPC runs passed R1 twice -> failed HLA typing)
            cmp -s <(head -2 {params.d}/{wildcards.run}_1.fastq | tail -1) <(head -2 {params.d}/{wildcards.run}_2.fastq | tail -1) \
                && {{ echo "R1 == R2: not paired-end?" >> {log}; exit 1; }}
            gzip -c {params.d}/{wildcards.run}_1.fastq > {output.r1}; gzip -c {params.d}/{wildcards.run}_2.fastq > {output.r2}
            rm -rf {params.d}/{wildcards.run}_[12].fastq {params.d}/{wildcards.run}
        fi
        """


rule trimmomatic:
    input: rules.fetch.output
    output:
        r1=f"{OUT}/runs/{{run}}/fastq/{{run}}_trimmed_1.fq.gz",
        r2=f"{OUT}/runs/{{run}}/fastq/{{run}}_trimmed_2.fq.gz",
        u1=temp(f"{OUT}/runs/{{run}}/fastq/{{run}}_unpaired_1.fq.gz"),
        u2=temp(f"{OUT}/runs/{{run}}/fastq/{{run}}_unpaired_2.fq.gz"),
        summary=f"{OUT}/runs/{{run}}/qc/{{run}}_trimmomatic_summary.txt",
    params:
        steps=config["trimmomatic_steps"],
        adapters=config["trimmomatic_adapters"],
    threads: 8
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/{{run}}_trimmomatic.log"
    shell:
        "trimmomatic PE -threads {threads} -phred33 -summary {output.summary} {input} "
        "{output.r1} {output.u1} {output.r2} {output.u2} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 {params.steps} > {log} 2>&1"


def stage_reads(stage):
    def f(wc):
        return [f"{OUT}/runs/{wc.run}/fastq/{wc.run}_{stage}_{i}.fq.gz" for i in (1, 2)]
    return f


rule fastqc:
    input: lambda wc: stage_reads(wc.stage)(wc)
    output: f"{OUT}/runs/{{run}}/qc/{{run}}_qc_{{stage}}.tsv"
    params: d=f"{OUT}/runs/{{run}}/qc/fastqc_{{stage}}", minq=config["min_phred"]
    wildcard_constraints: stage="raw|trimmed"
    threads: 2
    conda: "../../envs/qc.yaml"
    shell:
        "mkdir -p {params.d}; fastqc -q -t {threads} --extract -o {params.d} {input}; "
        "python3 {PIPE}/scripts/parse_fastqc.py --sample {wildcards.run} --stage {wildcards.stage} "
        "--min-phred {params.minq} --fastqc {params.d}/*_fastqc/fastqc_data.txt --out {output}"


def aligned_reads(wc):
    return stage_reads("trimmed" if config["align_trimmed"] else "raw")(wc)
