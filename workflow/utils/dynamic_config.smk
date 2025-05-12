"""
Rules for dynamic Auspice configuration generation
"""

rule generate_dynamic_auspice_config:
    input:
        base_config = lambda wildcards: pathogen_config.get("auspice_config_path", "config/auspice_config.json"),
        metadata = "data/metadata/{pathogen}_metadata.tsv" if segment_mode == "single" else "data/metadata/{pathogen}_metadata_{segment}.tsv",
        master_config = "config/master_config.yaml"
    output:
        dynamic_config = "results/configs/{pathogen}_auspice_config_dynamic.json" if segment_mode == "single" else "results/configs/{pathogen}_auspice_config_{segment}_dynamic.json"
    params:
        pathogen = lambda wildcards: wildcards.pathogen
    log:
        "logs/generate_dynamic_auspice_config_{pathogen}.log" if segment_mode == "single" else "logs/generate_dynamic_auspice_config_{pathogen}_{segment}.log"
    shell:
        """
        python scripts/generate_auspice_config.py \
            --base-config {input.base_config} \
            --metadata {input.metadata} \
            --master-config {input.master_config} \
            --pathogen {params.pathogen} \
            --output {output.dynamic_config} \
            > {log} 2>&1
        """