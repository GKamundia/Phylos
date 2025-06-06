"""
Rules for exporting data to Auspice JSON format
"""

# Dynamically get auspice_config based on pathogen
rule generate_dynamic_auspice_config:
    input:
        base_config = lambda wildcards: pathogen_config.get("auspice_config_path", "config/auspice_config.json"),
        metadata = "data/metadata/{pathogen}_metadata.tsv" if segment_mode == "single" else "data/metadata/{pathogen}_metadata_{segment}.tsv"
    output:
        dynamic_config = "results/configs/{pathogen}_auspice_config_dynamic.json" if segment_mode == "single" else "results/configs/{pathogen}_auspice_config_{segment}_dynamic.json"
    log:
        "logs/generate_dynamic_config_{pathogen}.log" if segment_mode == "single" else "logs/generate_dynamic_config_{pathogen}_{segment}.log"
    run:
        import os
        import shutil
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output.dynamic_config), exist_ok=True)
        
        # Copy the base config to the dynamic config location
        shutil.copy2(input.base_config, output.dynamic_config)
        
        # Log the operation
        with open(log[0], "w") as log_file:
            log_file.write(f"Copied {input.base_config} to {output.dynamic_config}\n")

rule export:
    input:
        # Get input files based on segment mode
        tree = "results/tree/{pathogen}_refined.nwk" if segment_mode == "single" else "results/segments/{segment}/tree/{pathogen}_{segment}_refined.nwk",
        metadata = "data/metadata/{pathogen}_metadata.tsv" if segment_mode == "single" else "data/metadata/{pathogen}_metadata.tsv",
        node_data = expand(
            "results/node_data/{pathogen}_{node_data}.json" if segment_mode == "single" else "results/segments/{segment}/node_data/{pathogen}_{segment}_{node_data}.json",
            node_data=config["export"]["node_data"],
            allow_missing=True
        ),
        dynamic_config = rules.generate_dynamic_auspice_config.output.dynamic_config
    output:
        auspice_json = "results/auspice/{pathogen}.json" if segment_mode == "single" else "results/segments/{segment}/auspice/{pathogen}_{segment}.json"
    params:
        node_data_files = lambda wildcards, input: " ".join([f"--node-data {file}" for file in input.node_data]),
        lat_longs = pathogen_config.get("lat_longs_path", "config/lat_longs.tsv"),
        title = lambda wildcards: f"{config['pathogen_name']} Genomic Surveillance" if segment_mode == "single" else f"{config['pathogen_name']} {wildcards.segment} Segment Surveillance"
    log:
        "logs/export_{pathogen}.log" if segment_mode == "single" else "logs/export_{pathogen}_{segment}.log"
    benchmark:
        "benchmarks/export_{pathogen}.txt" if segment_mode == "single" else "benchmarks/export_{pathogen}_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("export", {}).get("mem_mb", 2000)
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            {params.node_data_files} \
            --auspice-config {input.dynamic_config} \
            --lat-longs {params.lat_longs} \
            --title "{params.title}" \
            --output {output.auspice_json} \
            > {log} 2>&1
        """

# Rule for combining segment-specific JSONs if in multi-segment mode
if segment_mode == "multi" and config.get("workflow", {}).get("create_combined_view", True):
    rule combine_jsons:
        input:
            jsons = expand("results/segments/{segment}/auspice/{pathogen}_{segment}.json", segment=segments, pathogen=output_prefix)
        output:
            combined = f"results/auspice/{output_prefix}_combined.json"
        log:
            f"logs/combine_jsons_{output_prefix}.log"
        benchmark:
            f"benchmarks/combine_jsons_{output_prefix}.txt"
        resources:
            mem_mb = config["resources"].get("combine", {}).get("mem_mb", 2000)
        shell:
            """
            python scripts/combine_segment_jsons.py \
                --input-jsons {input.jsons} \
                --output {output.combined} \
                > {log} 2>&1
            """