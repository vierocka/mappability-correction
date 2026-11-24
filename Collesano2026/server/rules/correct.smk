# Correction with the original script (apply_uniqueness_and_average.py, PIPE_BASE/PIPE_CONFIG hooks):
# exclusions (excluded_transcripts/), alt-copy fold, ×funique only if f < 0.9 + pseudogene evidence
# + shared < 0.10 [Robert & Watson 2015], HLA total = allele fcM sum split by unique-read ratio,
# TPM, replicate mean per condition. Then make_hla_tap_tables.py (HLA-A/B/C, TAP1/2 table).
import json


def correct_inputs(wc):
    kinds = ["MANEpersUTR"] + (["mane"] if config["hla_total_source"] == "basic_mane_fc" else [])
    files = [f"{OUT}/correction/{group_of(r)}/raw_counts/{r}_{group_of(r)}_{k}_{c}.counts.csv"
             for r in RUNS for k in kinds for c in ("fc", "fcM")]
    tables = sorted({f"{FUNIQUE_DIR}/uniqueness_correction_L{group_read_length(g)}bp.tsv" for g in GROUPS})
    return files + tables


def drb_ids(names):
    ids = {"HLA-DRB3": "ENST00000383126.8", "HLA-DRB4": "ENST00000457451.7"}
    return {ids[n]: n for n in names.split(";") if n in ids}


rule correct:
    input: correct_inputs
    output:
        expand(f"{OUT}/correction/{{r.group}}/corrected_counts/{{r.run}}_mane_fc_corrected.tsv",
               r=SAMPLES.itertuples()),
        expand(f"{OUT}/correction/{{r.group}}/corrected_counts/{{r.group}}_corrected_mean_{{r.condition}}_mane_fc_corrected.tsv",
               r=SAMPLES.drop_duplicates(["group", "condition"]).itertuples()),
        cfg=f"{OUT}/correction/pipe_config.json",
        hla=f"{OUT}/correction/HLA_corrected_mean_values.tsv",
    params: base=f"{OUT}/correction"
    log: f"{OUT}/logs/correct.log"
    run:
        import shutil
        b = Path(params.base)
        (b / "uniqueness_factors").mkdir(parents=True, exist_ok=True)
        for f in ("non_primary_MANE_transcripts.tsv", "non_coding_MANE_transcripts_no_protein.tsv"):
            shutil.copy(f"{PIPE}/resources/{f}", b / f)
        shutil.copy(f"{FUNIQUE_DIR}/low_uniqueness_MANE_transcripts.tsv", b / "uniqueness_factors")
        cell_cfg, include, conds = {}, {}, []
        for g in GROUPS:
            L = group_read_length(g)
            shutil.copy(f"{FUNIQUE_DIR}/uniqueness_correction_L{L}bp.tsv", b / "uniqueness_factors")
            sub = SAMPLES[SAMPLES.group == g]
            c = {k: list(v.run) for k, v in sub.groupby("condition", sort=False)}
            cell_cfg[g] = {"factors_table": f"uniqueness_correction_L{L}bp.tsv", "conditions": c}
            conds += [(g, k, v) for k, v in c.items()]
            names = ";".join(x for x in sub.include_drb if x)
            if names:
                include[g] = drb_ids(names)
        cfg = {"CELL_CONFIG": cell_cfg, "CELL_SPECIFIC_INCLUDE": include, "CONDITIONS": conds,
               "hla_total_source": config["hla_total_source"]}
        json.dump(cfg, open(output.cfg, "w"), indent=1)
        shell(f"PIPE_BASE={b} PIPE_CONFIG={output.cfg} python3 {PIPE}/scripts/apply_uniqueness_and_average.py > {log} 2>&1")
        shell(f"PIPE_BASE={b} PIPE_CONFIG={output.cfg} python3 {PIPE}/scripts/make_hla_tap_tables.py >> {log} 2>&1")
