# Reports: filter-step losses and composition shift per group; one-table summaries.


rule filter_steps:
    input:
        lambda wc: [f"{OUT}/correction/{wc.group}/corrected_counts/{r}_mane_fc_corrected.tsv" for r in runs_of(wc.group)],
        rules.correct.output.cfg,
    output: f"{OUT}/reports/{{group}}_steps.tsv"
    params:
        files=lambda wc: [f"{OUT}/correction/{wc.group}/corrected_counts/{r}_mane_fc_corrected.tsv" for r in runs_of(wc.group)],
        labels=lambda wc: [f"{wc.group}:{r}" for r in runs_of(wc.group)],
    conda: "../../envs/py.yaml"
    shell:
        "python3 {PIPE}/scripts/report_filter_steps.py --samples {params.files} --labels {params.labels} "
        "--resources {PIPE}/resources --min-tpm 1 --out-prefix {OUT}/reports/{wildcards.group}"


rule summary:
    input:
        hla=expand(f"{OUT}/groups/{{group}}/hla/{{group}}_hla_calls.tsv", group=GROUPS),
        conc=expand(f"{OUT}/groups/{{group}}/hla/{{group}}_hla_run_concordance.tsv", group=GROUPS),
        qc=expand(f"{OUT}/runs/{{run}}/qc/{{run}}_qc_{{stage}}.tsv", run=RUNS,
                  stage=["raw", "trimmed"]),
        rl=expand(f"{OUT}/runs/{{run}}/qc/{{run}}_read_length.txt", run=RUNS),
        star=expand(f"{OUT}/runs/{{r.run}}/align/{{r.run}}_{{r.group}}_MANEpersUTR_Log.final.out",
                    r=SAMPLES.itertuples()),
        corr=rules.correct.output,
    output:
        hla=f"{OUT}/summary/hla_types.tsv",
        qc=f"{OUT}/summary/qc_reads.tsv",
        tpm=f"{OUT}/summary/corrected_TPM_matrix.tsv",
        png=f"{OUT}/summary/overview.png",
    params: samples=config["samples"]
    conda: "../../envs/py.yaml"
    shell:
        "python3 {PIPE}/scripts/summarize_server.py --samples {params.samples} --outdir {OUT} "
        "--min-phred {config[min_phred]}"
