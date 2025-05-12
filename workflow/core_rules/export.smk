"""
Rules for exporting data to Auspice JSON format with dynamic configuration
"""

rule export:
    input:
        tree = "results/tree/{pathogen}_tree.nwk" if segment_mode == "single" else "results/segments/{segment}/tree/{pathogen}_{segment}_tree.nwk",
        node_data = expand("results/node_data/{pathogen}_{node_data}.json", 
                           node_data=config["export"]["node_data"], 
                           pathogen=output_prefix) if segment_mode == "single" else 
                    expand("results/segments/{segment}/node_data/{pathogen}_{segment}_{node_data}.json", 
                           node_data=config["export"]["node_data"], 
                           pathogen=output_prefix, 
                           segment="{segment}"),
        dynamic_config = "results/configs/{pathogen}_auspice_config_dynamic.json" if segment_mode == "single" else "results/configs/{pathogen}_auspice_config_{segment}_dynamic.json"
    output:
        auspice_json = "results/auspice/{pathogen}.json" if segment_mode == "single" else "results/segments/{segment}/auspice/{pathogen}_{segment}.json"
    params:
        auspice_config = lambda wildcards, input: input.dynamic_config,
        lat_longs = pathogen_config.get("lat_longs_path", "config/lat_longs.tsv"),
        title = lambda wildcards: f"{config['pathogen_name']} Genomic Surveillance" if segment_mode == "single" else f"{config['pathogen_name']} {wildcards.segment} Segment Surveillance",
        colors = lambda w: f"--colors {config['export']['colors']}" if config.get('export', {}).get('colors') else "",
        extra_node_data = lambda w: " ".join([f"--node-data {data}" for data in config.get("export", {}).get("extra_node_data", [])])
    log:
        f"logs/export_{output_prefix}.log"
    benchmark:
        f"benchmarks/export_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("export", {}).get("mem_mb", 2000)
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.auspice_json})
        
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.node_data} \
            {params.extra_node_data} \
            --auspice-config {params.auspice_config} \
            --lat-longs {params.lat_longs} \
            {params.colors} \
            --title "{params.title}" \
            --output {output.auspice_json} \
            > {log} 2>&1
        """