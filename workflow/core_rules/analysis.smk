"""
Rules for phylogenetic analysis and tree building
"""

# Get alignment input based on whether masking is enabled
def get_alignment_input(wildcards):
    if config.get("mask", {}).get("sites", []) or \
       config.get("mask", {}).get("from_beginning") or \
       config.get("mask", {}).get("from_end"):
        return f"results/masked/{output_prefix}_masked.fasta"
    else:
        return f"results/aligned/{output_prefix}_aligned.fasta"

# Build phylogenetic tree
rule tree:
    input:
        alignment = get_alignment_input
    output:
        tree = f"results/tree/{output_prefix}_tree.nwk"
    params:
        method = config.get("tree", {}).get("method", "iqtree"),
        iqtree_args = config.get("tree", {}).get("iqtree_args", "-ninit 2 -n 2"),
        substitution_model = config.get("tree", {}).get("substitution_model", "GTR")
    threads: config["resources"].get("tree", {}).get("threads", 4)
    log:
        f"logs/tree_{output_prefix}.log"
    benchmark:
        f"benchmarks/tree_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.tree})
        
        # Run tree building
        if [ "{params.method}" = "iqtree" ]; then
            augur tree \
                --alignment {input.alignment} \
                --output {output.tree} \
                --method iqtree \
                --substitution-model {params.substitution_model} \
                --nthreads {threads} \
                {params.iqtree_args} \
                > {log} 2>&1
        else
            # Default to FastTree
            augur tree \
                --alignment {input.alignment} \
                --output {output.tree} \
                --method fasttree \
                --nthreads {threads} \
                > {log} 2>&1
        fi
        """

# Refine tree with time and metadata
rule refine:
    input:
        tree = f"results/tree/{output_prefix}_tree.nwk",
        alignment = get_alignment_input,
        metadata = get_input_metadata
    output:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        node_data = f"results/node_data/{output_prefix}_branch_lengths.json"
    params:
        coalescent = config["refine"].get("coalescent", "opt"),
        date_inference = config["refine"].get("date_inference", "marginal"),
        clock_filter_iqd = config["refine"].get("clock_filter_iqd", 4),
        clock_rate = lambda w: f"--clock-rate {config['refine']['clock_rate']}" if config["refine"].get("clock_rate") else "",
        clock_std_dev = lambda w: f"--clock-std-dev {config['refine']['clock_std_dev']}" if config["refine"].get("clock_std_dev") else "",
        root = lambda w: f"--root {config['refine']['root']}" if config["refine"].get("root") else ""
    threads: config["resources"].get("refine", {}).get("threads", 2)
    log:
        f"logs/refine_{output_prefix}.log"
    benchmark:
        f"benchmarks/refine_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("refine", {}).get("mem_mb", 4000)
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            {params.clock_rate} \
            {params.clock_std_dev} \
            {params.root} \
            --date-confidence \
            > {log} 2>&1
        """

# Reconstruct ancestral sequences and mutations
rule ancestral:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        alignment = get_alignment_input
    output:
        node_data = f"results/node_data/{output_prefix}_nt_muts.json"
    params:
        inference = config.get("ancestral", {}).get("inference", "joint")
    threads: config["resources"].get("ancestral", {}).get("threads", 1)
    log:
        f"logs/ancestral_{output_prefix}.log"
    benchmark:
        f"benchmarks/ancestral_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("ancestral", {}).get("mem_mb", 2000)
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            > {log} 2>&1
        """

# Reconstruct ancestral traits (e.g., country, host)
rule traits:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk",
        metadata = get_input_metadata
    output:
        node_data = f"results/node_data/{output_prefix}_traits.json"
    params:
        columns = lambda w: ",".join(config.get("traits", {}).get("columns", ["country"]))
    log:
        f"logs/traits_{output_prefix}.log"
    benchmark:
        f"benchmarks/traits_{output_prefix}.txt"
    resources:
        mem_mb = config["resources"].get("traits", {}).get("mem_mb", 2000)
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            > {log} 2>&1
        """