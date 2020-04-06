configfile: "config.yaml"
rule all:
    input: f"{config['run_id']}_predictions.txt"

rule generate_feature_matrix:
    input:
        control=config["coordinates"]["control"],
        disease=config["coordinates"]["disease"]
    output:
        f"data/{config['run_id']}_combined_matrix.txt"
    script:
        "gfm_script.py"

rule rf_model:
    input:
        f"data/{config['run_id']}_combined_matrix.txt"
    output:
        f"{config['run_id']}_predictions.txt"
    script:
        "rfm_script.py"
