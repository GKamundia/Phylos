"""
Rules for phylogenetic analysis and tree building
"""

# Build phylogenetic tree (Step 1: ML search/model selection)
rule ml_tree:
    input:
        alignment = get_alignment_input
    output:
        prefix = f"results/tree/{output_prefix}" if segment_mode == "single" else f"results/segments/{{segment}}/tree/rvf_{{segment}}"
    params:
        n = 10,
        ninit = 100,
        seed = 12345,
    threads: config["resources"].get("tree", {}).get("threads", 0) or os.cpu_count()
    log:
        f"logs/ml_tree_{output_prefix}.log" if segment_mode == "single" else "logs/ml_tree_rvf_{segment}.log"
    benchmark:
        f"benchmarks/ml_tree_{output_prefix}.txt" if segment_mode == "single" else "benchmarks/ml_tree_rvf_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)
    shell:
        """
        python scripts/run_tree.py --alignment {input.alignment} \
            --output {output.prefix}.treefile \
            --log {log} \
            --ml-prefix {output.prefix} \
            --n {params.n} --ninit {params.ninit} --seed {params.seed}
        """

# Step 2: Branch support estimation (publication-ready)
rule support_tree:
    input:
        alignment = get_alignment_input,
        ml_prefix = rules.ml_tree.output.prefix
    output:
        prefix = f"results/tree/{output_prefix}_support" if segment_mode == "single" else f"results/segments/{{segment}}/tree/rvf_{{segment}}_support"
    params:
        seed = 12345,
        bb = 1000,
        alrt = 1000,
    threads: config["resources"].get("tree", {}).get("threads", 0) or os.cpu_count()
    log:
        f"logs/support_tree_{output_prefix}.log" if segment_mode == "single" else "logs/support_tree_rvf_{segment}.log"
    benchmark:
        f"benchmarks/support_tree_{output_prefix}.txt" if segment_mode == "single" else "benchmarks/support_tree_rvf_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)
    shell:
        """
        python scripts/run_tree.py --alignment {input.alignment} \
            --output {output.prefix}.treefile \
            --log {log} \
            --ml-prefix {input.ml_prefix} \
            --support-prefix {output.prefix} \
            --seed {params.seed} --bb {params.bb} --alrt {params.alrt}
        """

# Refine tree with time and metadata
rule refine:
    input:
        tree = f"results/tree/{output_prefix}_tree.nwk" if segment_mode == "single" else "results/segments/{{segment}}/tree/rvf_{{segment}}_tree.nwk",
        alignment = get_alignment_input,
        metadata = get_input_metadata
    output:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else "results/segments/{{segment}}/tree/rvf_{{segment}}_refined.nwk",
        node_data = f"results/node_data/{output_prefix}_branch_lengths.json" if segment_mode == "single" else "results/segments/{{segment}}/node_data/rvf_{{segment}}_branch_lengths.json"
    params:
        coalescent = config["refine"].get("coalescent", "opt"),
        date_inference = config["refine"].get("date_inference", "marginal"),
        clock_filter_iqd = config["refine"].get("clock_filter_iqd", 4)
    threads: 
        config["resources"].get("refine", {}).get("threads", 2)
    log:
        f"logs/refine_{output_prefix}.log" if segment_mode == "single" else "logs/refine_rvf_{segment}.log"
    benchmark:
        f"benchmarks/refine_{output_prefix}.txt" if segment_mode == "single" else "benchmarks/refine_rvf_{segment}.txt"
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
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --date-format "%Y-%m-%d" \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --date-confidence \
            --keep-root \
            > {log} 2>&1
        """

# Reconstruct ancestral sequences and mutations
rule ancestral:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else "results/segments/{segment}/tree/rvf_{segment}_refined.nwk",
        alignment = get_alignment_input
    output:
        node_data = f"results/node_data/{output_prefix}_nt_muts.json" if segment_mode == "single" else "results/segments/{segment}/node_data/rvf_{segment}_nt_muts.json"
    params:
        inference = config.get("ancestral", {}).get("inference", "joint")
    threads: config["resources"].get("ancestral", {}).get("threads", 1)
    log:
        f"logs/ancestral_{output_prefix}.log" if segment_mode == "single" else "logs/ancestral_rvf_{segment}.log"
    benchmark:
        f"benchmarks/ancestral_{output_prefix}.txt" if segment_mode == "single" else "benchmarks/benchmarks_rvf_{segment}.txt"
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

# Reconstruct ancestral traits (e.g., country, division)
rule traits:
    input:
        tree = f"results/tree/{output_prefix}_refined.nwk" if segment_mode == "single" else "results/segments/{segment}/tree/rvf_{segment}_refined.nwk",
        metadata = get_input_metadata
    output:
        node_data = f"results/node_data/{output_prefix}_traits.json" if segment_mode == "single" else "results/segments/{segment}/node_data/rvf_{segment}_traits.json"
    params:
        columns = lambda w: ",".join(config.get("traits", {}).get("columns", ["country"]))
    threads: 
        config["resources"].get("traits", {}).get("threads", 1)
    log:
        f"logs/traits_{output_prefix}.log" if segment_mode == "single" else "logs/traits_rvf_{segment}.log"    
    benchmark:
        f"benchmarks/traits_{output_prefix}.txt" if segment_mode == "single" else "benchmarks/traits_rvf_{segment}.txt"
    resources:
        mem_mb = config["resources"].get("traits", {}).get("mem_mb", 2000)
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession \
            --output {output.node_data} \
            --columns {params.columns} \
            > {log} 2>&1
        """