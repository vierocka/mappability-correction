# funique table per read length [Derrien 2012; Storvall 2013].
# Mean read length per mate (FastQC of the aligned reads): 73-77 bp -> L75, 98-102 bp -> L100
# (resources/funique, HPC 2026-08-12). Other lengths: simulated once with the original worker
# (compute_uniqueness.py) and cached in resources/funique/.
# Method: MANE exons ± 50 bp tiled at 1-bp step, pseudo-pairs (R2 = revcomp R1), STAR 2.7.11b
# [Dobin 2013] vs Ensembl GRCh38 primary assembly (194 seqs), funique = MAPQ-255 fraction of R1.
# NOTE: the pseudogene-evidence table (low_uniqueness_MANE_transcripts.tsv) covers L75/L100 only;
# at a new length, transcripts with f < 0.9 that are missing from it are not corrected (ROADMAP).

FUNIQUE_DIR = f"{PIPE}/resources/funique"


checkpoint read_length:
    input:
        f"{OUT}/runs/{{run}}/qc/{{run}}_qc_{'trimmed' if config['align_trimmed'] else 'raw'}.tsv",
    output:
        f"{OUT}/runs/{{run}}/qc/{{run}}_read_length.txt",
    conda: "../../envs/py.yaml"
    shell:
        "python3 {PIPE}/scripts/choose_read_length.py --qc {input} --out {output}"


def group_read_length(group):
    Ls = {int(open(checkpoints.read_length.get(run=r).output[0]).read().split()[0]) for r in runs_of(group)}
    if len(Ls) != 1:
        raise ValueError(f"{group}: runs need different funique tables {sorted(Ls)}; split the group")
    return Ls.pop()


def funique_table(wc):
    return f"{FUNIQUE_DIR}/uniqueness_correction_L{group_read_length(wc.group)}bp.tsv"


rule funique_genome_index:
    input: config["genome_fa"]
    output: directory(f"{config['refs']}/STAR_index_primary_assembly")
    threads: config["threads"]
    resources: mem_mb=40000
    conda: "../../envs/genome.yaml"
    shell:
        "n=$(grep -c '>' {input}); [ $n -eq 194 ] || {{ echo 'genome_fa: expected 194-seq primary assembly'; exit 1; }}; "
        "mkdir -p {output}; STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output} --genomeFastaFiles {input}"


rule funique_simulate:
    input:
        gff=config["mane_gff"],
        fa=config["genome_fa"],
        idx=rules.funique_genome_index.output,
    output: f"{FUNIQUE_DIR}/transcript_uniqueness_factors_L{{L}}bp.tsv"
    params: run=lambda wc: f"{config['refs']}/funique_L{wc.L}"
    threads: config["threads"]
    resources: mem_mb=40000
    conda: "../../envs/genome.yaml"
    log: f"{OUT}/logs/funique_L{{L}}.log"
    shell:
        """
        mkdir -p {params.run}/bam; ulimit -n $(ulimit -Hn)
        python3 {PIPE}/scripts/compute_uniqueness.py simulate-base --gff {input.gff} --genome-fa {input.fa} \
            --read-len {wildcards.L} --flank 50 --outdir {params.run} > {log} 2>&1
        STAR --runThreadN {threads} --genomeDir {input.idx} \
            --readFilesIn {params.run}/simreads/sim_R1.fastq {params.run}/simreads/sim_R2.fastq \
            --outSAMtype BAM Unsorted --outSAMattributes NH HI AS NM --outSAMmultNmax 1 \
            --outFilterMultimapNmax 40 --outFilterMismatchNmax 1 \
            --outFileNamePrefix {params.run}/bam/sim_ >> {log} 2>&1
        python3 {PIPE}/scripts/compute_uniqueness.py score --outdir {params.run} \
            --bam {params.run}/bam/sim_Aligned.out.bam --out-name $(basename {output}) >> {log} 2>&1
        cp {params.run}/$(basename {output}) {output}; rm -rf {params.run}/simreads {params.run}/bam
        """


rule funique_table:
    input: f"{FUNIQUE_DIR}/transcript_uniqueness_factors_L{{L}}bp.tsv"
    output: f"{FUNIQUE_DIR}/uniqueness_correction_L{{L}}bp.tsv"
    params: mane_dir=lambda wc: str(Path(config["mane_rna"]).parent)
    conda: "../../envs/py.yaml"
    shell:
        "PIPE_MANE_DIR={params.mane_dir} PIPE_FACTORS_DIR={FUNIQUE_DIR} PIPE_OUT_DIR={FUNIQUE_DIR} "
        "PIPE_READ_LENS={wildcards.L} python3 {PIPE}/scripts/build_correction_tables.py"
