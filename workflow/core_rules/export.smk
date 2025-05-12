"""
Rules for exporting results to Auspice JSON format
"""

# Export to auspice
rule export:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        metadata = get_input_metadata,
        branch_lengths = f"results/node_data/{output_prefix}_branch_lengths.json",
        traits = f"results/node_data/{output_prefix}_traits.json",
        nt_muts = f"results/node_data/{output_prefix}_nt_muts.json"
    output:
        auspice_json = f"results/auspice/{output_prefix}.json"
    params:
        auspice_config = pathogen_config.get("auspice_config_path", "config/auspice_config.json"),
        lat_longs = pathogen_config.get("lat_longs_path", "config/lat_longs.tsv"),
        title = f"{config['pathogen_name']} Genomic Surveillance",
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
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
            {params.extra_node_data} \
            --auspice-config {params.auspice_config} \
            --lat-longs {params.lat_longs} \
            {params.colors} \
            --title "{params.title}" \
            --output {output.auspice_json} \
            > {log} 2>&1
        """