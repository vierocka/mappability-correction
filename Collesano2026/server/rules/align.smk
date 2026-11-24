# Personal MANE transcriptome used as a pseudo-genome [Morales 2022], STAR 2.7.10b [Dobin 2013],
# featureCounts (subread 2.1.1) [Liao 2014]: fc = unique pairs, fcM = multimappers split by fraction
# [Deschamps-Francoeur 2020]. Flags as in Laura_MANE_personalized_UTR.sh (original_shell/).


rule mane_saf:
    input: config["mane_rna"]
    output: f"{config['refs']}/MANE/cds_mane.saf"
    conda: "../../envs/py.yaml"
    shell: "python3 {PIPE}/scripts/make_mane_saf.py {input} {output}"


def index_fasta(wc):
    if wc.ref == "basic":
        return config["mane_rna"]
    return f"{OUT}/groups/{wc.ref}/ref/{wc.ref}_personal_MANE.fa"


rule star_index_mane:
    input: index_fasta
    output: directory(f"{OUT}/index/STAR_{{ref}}")
    threads: config["threads"]
    resources: mem_mb=8000
    conda: "../../envs/align.yaml"
    log: f"{OUT}/logs/STAR_index_{{ref}}.log"
    shell:
        "mkdir -p {output}; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} "
        "--genomeFastaFiles {input} --genomeSAindexNbases 12 --genomeChrBinNbits 11 "
        "--outFileNamePrefix {output}/ > {log} 2>&1"


def align_inputs(wc):
    ref = group_of(wc.run) if wc.kind == "MANEpersUTR" else "basic"
    saf = (f"{OUT}/groups/{ref}/ref/{ref}_cds_mane.saf" if wc.kind == "MANEpersUTR"
           else f"{config['refs']}/MANE/cds_mane.saf")
    return {"reads": aligned_reads(wc), "idx": f"{OUT}/index/STAR_{ref}", "saf": saf}


rule star_featurecounts:
    input: unpack(align_inputs)
    output:
        fc=f"{OUT}/correction/{{group}}/raw_counts/{{run}}_{{group}}_{{kind}}_fc.counts.csv",
        fcm=f"{OUT}/correction/{{group}}/raw_counts/{{run}}_{{group}}_{{kind}}_fcM.counts.csv",
        log_final=f"{OUT}/runs/{{run}}/align/{{run}}_{{group}}_{{kind}}_Log.final.out",
    wildcard_constraints: kind="MANEpersUTR|mane"
    params: pre=f"{OUT}/runs/{{run}}/align/{{run}}_{{group}}_{{kind}}_"
    threads: config["threads"]
    resources: mem_mb=16000
    conda: "../../envs/align.yaml"
    log: f"{OUT}/logs/{{run}}_{{group}}_{{kind}}_align.log"
    shell:
        """
        STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} \
            --readFilesIn {input.reads} --readFilesCommand zcat \
            --alignIntronMax 1 --alignIntronMin 2 --outFilterMultimapNmax 20 \
            --outSAMattributes NH HI AS NM MD --runRNGseed 0 \
            --outSAMtype BAM Unsorted --outFileNamePrefix {params.pre} > {log} 2>&1
        featureCounts -T {threads} -F SAF -p -B -C --countReadPairs -s 0 -a {input.saf} \
            -o {output.fc} {params.pre}Aligned.out.bam >> {log} 2>&1
        featureCounts -T {threads} -F SAF -p -B -C --countReadPairs -s 0 -M --fraction -a {input.saf} \
            -o {output.fcm} {params.pre}Aligned.out.bam >> {log} 2>&1
        mv {params.pre}Log.final.out {output.log_final}; rm {params.pre}Aligned.out.bam
        """
