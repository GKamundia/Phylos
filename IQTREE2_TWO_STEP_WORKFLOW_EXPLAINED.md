# IQ-TREE2 Two-Step Phylogenetic Analysis Workflow: Parameters, Outputs, and Best Practices

## Overview: Why a Two-Step Workflow?

The two-step IQ-TREE2 workflow is the gold standard for robust, publication-ready phylogenetic inference. It separates the maximum likelihood (ML) tree search/model selection (Step 1) from branch support estimation (Step 2), ensuring reproducibility and statistical rigor ([IQ-TREE2 Manual](http://www.iqtree.org/doc/iqtree2-manual), [Best Practices](https://doi.org/10.1093/molbev/msy096)).

**Key update:** This workflow uses `--keep-ident` in both steps, ensuring that all identical sequences are retained in the analysis. This is essential for:

- Accurate sample representation in the final tree
- Downstream compatibility (e.g., with Nextstrain/Auspice, TreeTime)
- Full reproducibility and transparency for publication

_Keeping all identical sequences ensures compatibility with tools like Nextstrain (Auspice JSON) and TreeTime, which assume that each sequence maps to a unique metadata entry._

---

## Step 1: ML Tree Search and Model Selection (with All Identical Sequences Retained)

**Command used:**

```
iqtree2 -s results/segments/L/aligned/rvf_L_aligned.fasta -m MFP -nt AUTO -n 10 --ninit 100 --safe --seed 12345 --keep-ident
```

**Explanation:**

- `-s results/segments/L/aligned/rvf_L_aligned.fasta`: Input alignment (all L segment sequences, including identicals)
- `-m MFP`: ModelFinder Plus; automatically selects the best-fit evolutionary model
- `-nt AUTO`: Use all available CPU threads
- `-n 10`: 10 independent ML tree searches (avoids local optima)
- `--ninit 100`: 100 initial parsimony trees per search
- `--safe`: Safe mode for reproducibility
- `--seed 12345`: Set random seed for reproducibility
- `--keep-ident`: **Retain all identical sequences** (prevents default deduplication)

**What happens:**

- IQ-TREE2 performs multiple ML tree searches, each starting from many initial trees, to find the best tree and model for your data.
- All identical sequences are kept in the analysis, ensuring the output tree and statistics represent the full sample set ([IQ-TREE2 Manual: --keep-ident](http://www.iqtree.org/doc/iqtree2-output#cmdoption-keep-ident)).

**How to extract the best-fit model for Step 2:**

```bash
grep "Best-fit model:" results/segments/L/aligned/rvf_L_aligned.fasta.iqtree
```

---

## Step 2: Branch Support Estimation (with All Identical Sequences Retained)

**Command used:**

```
iqtree2 -s results/segments/L/aligned/rvf_L_aligned.fasta -m GTR+F+I+G4 -nt AUTO --keep-ident --safe --seed 12345 -bb 1000 -alrt 1000 -pre results/segments/L/aligned/rvf_L_aligned_support --redo
```

_Replace `GTR+F+I+G4` with the actual best-fit model from Step 1 (see above)._

**Explanation:**

- `-s results/segments/L/aligned/rvf_L_aligned.fasta`: Input alignment (must match Step 1)
- `-m GTR+F+I+G4`: Use the best-fit model found in Step 1
- `-nt AUTO`: Use all available CPU threads
- `--keep-ident`: **Retain all identical sequences**
- `--safe`: Safe mode for reproducibility
- `--seed 12345`: Set random seed for reproducibility
- `-bb 1000`: 1000 ultrafast bootstrap replicates (branch support)
- `-alrt 1000`: 1000 SH-aLRT support replicates (alternative support metric)
- `-pre results/segments/L/aligned/rvf_L_aligned_support`: Prefix for all output files from this step
- `--redo`: Overwrite previous results if they exist

**What happens:**

- IQ-TREE2 estimates branch support values for the ML tree using two robust methods.
- All identical sequences are retained, so support values correspond to the full sample set.
- _Note:_ When using `--keep-ident`, support estimation (`-bb`, `-alrt`) requires full model re-optimization, so `-ft` (fixed topology) is not used. This ensures that support values are calculated for the ML tree under the specified model and full dataset ([IQ-TREE2 Manual: --keep-ident](http://www.iqtree.org/doc/iqtree2-output#cmdoption-keep-ident)).
- Using `-pre` keeps Step 2 outputs clean and separate from Step 1 files, improving reproducibility and preventing overwriting.
- If you add `-ft results/segments/L/aligned/rvf_L_aligned.fasta.treefile`, IQ-TREE2 will attempt to fix the topology, but support estimation (`-bb`, `-alrt`) requires full ML re-optimization and is not compatible with a fixed tree when using `--keep-ident`. If you try to use both, you may see errors or unexpected behavior. Always omit `-ft` for support estimation with `--keep-ident` ([IQ-TREE2 Manual: --keep-ident](http://www.iqtree.org/doc/iqtree2-output#cmdoption-keep-ident)).
- IQ-TREE2 does not allow branch support estimation (-bb, -alrt) on a fixed tree (-ft) when the model contains empirical frequencies (like +F) or when the model was selected using -m MFP.

---

## Why Use `--keep-ident`?

- **Publication and reproducibility:** Ensures all samples are represented in the final tree, as required for transparent reporting.
- **Downstream compatibility:** Some tools (e.g., Nextstrain/Auspice, TreeTime) require all input sequences to be present in the tree.
- **Avoids errors:** Prevents mismatches between alignment and tree files that can occur if deduplication is used.
- **Best practice:** Recommended for all analyses intended for publication or public sharing ([IQ-TREE2 Manual](http://www.iqtree.org/doc/iqtree2-output#cmdoption-keep-ident)).

---

## Output Files in `results/segments/L/aligned/`

| File                             | Use in Visualization?    | Example Tool           | Description                                                                             | Reference                                                  |
| -------------------------------- | ------------------------ | ---------------------- | --------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| `rvf_L_aligned.fasta`            |                          |                        | Input alignment (all L segment sequences, including identicals)                         | [Manual: Input](http://www.iqtree.org/doc/iqtree2-input)   |
| `rvf_L_aligned.fasta.iqtree`     | ✅ (for model reporting) | Methods section        | Run summary: best-fit model, ML tree stats, search details, and model selection results | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned.fasta.treefile`   | ✅ (raw ML tree)         | FigTree, iTOL          | ML tree in Newick format, with all samples (including identicals) as tips               | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned.fasta.log`        |                          |                        | Full log of the run: command, progress, warnings, and errors                            | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned.fasta.model.gz`   |                          |                        | Compressed file with detailed model selection results (all models tested and scores)    | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned_support.iqtree`   | ✅ (for model reporting) | Methods section        | Support estimation summary (from step 2)                                                | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned_support.treefile` | ✅ (with support)        | FigTree, iTOL, Auspice | ML tree with branch support values (from step 2)                                        | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |
| `rvf_L_aligned_support.log`      |                          |                        | Full log of the support estimation run                                                  | [Manual: Output](http://www.iqtree.org/doc/iqtree2-output) |

**How to interpret these files:**

- `.fasta`: The exact alignment used for all analyses; must be identical for both steps.
- `.iqtree`: Contains the best-fit model (e.g., `GTR+F+I+G4`), log-likelihoods, and search statistics. Use this file to report your model and tree search details.
- `.treefile`: The ML tree, with all samples as tips. This is the tree to use for visualization and downstream analysis.
- `.log`: The full run log, useful for troubleshooting and reproducibility.
- `.model.gz`: All tested models and their scores; useful for supplementary materials or further analysis.
- `_support.*`: Output files from the support estimation step, including the tree with support values and run summary.

---

## Best Practices for Publication-Quality Phylogenetics

- **Always use `--keep-ident`** for publication and downstream compatibility.
- **Set random seeds** for reproducibility.
- **Use sufficient replicates** (`-bb 1000`, `-alrt 1000`) for robust support.
- **Document all parameters and outputs** in your methods and supplementary materials.
- **Reference authoritative resources** in your documentation and publications.

---

## References

- [IQ-TREE2 Manual](http://www.iqtree.org/doc/iqtree2-manual)
- [IQ-TREE2 Output Files](http://www.iqtree.org/doc/iqtree2-output)
- [Best Practices in Phylogenetic Analysis](https://doi.org/10.1093/molbev/msy096)

---

_Prepared for robust, reproducible phylogenetic analysis using IQ-TREE2, with all identical sequences retained for publication and downstream use._

Rule refine:
augur refine --tree results/segments/L/aligned/rvf_L_aligned_support.treefile --alignment results/segments/L/aligned/rvf_L_aligned.fasta --metadata results/segments/L/filtered/rvf_L_metadata.tsv --output-tree results/segments/L/tree/rvf_L_refined.nwk --output-node-data results/segments/L/node_data/rvf_L_branch_lengths.json --coalescent opt --date-inference marginal --date-format "%Y-%m-%d"
--clock-filter-iqd 4 --date-confidence --keep-root

## STEPS THAT WORKED

1. Installed conda install -c conda-forge -c bioconda augur

1. tree building:
   (nextstrain-py39) root@Anarchy:/mnt/c/Users/Anarchy/Documents/Data_Science/NextStrain/rvf-nextstrain# augur tree --method iqtree --alignment results/segments/L/aligned/rvf_L_aligned.fasta --substitution-model auto --output results/segments/L/aligned/rvf_L_aligned.tree --nthreads auto --tree-builder-args="--ninit 100 -n 10 --nstop 500 -bb 1000 --alrt 1000" --override-default-args

1. refine
   (nextstrain-py39) root@Anarchy:/mnt/c/Users/Anarchy/Documents/Data_Science/NextStrain/rvf-nextstrain# augur refine --alignment results/segments/L/aligned/rvf_L_aligned.fasta --tree results/segments/L/aligned/rvf_L_aligned_support.treefile --metadata results/segments/L/filtered/rvf_L_metadata_minimal.tsv --output-tree results/segments/L/refined/rvf_L_refined.tree --output-node-data results/segments/L/refined/rvf_L_refined_node_data.json --date-confidence --coalescent opt --clock-filter-iqd 4 --keep-root --date-inference marginal --timetree

1. Ancestral
   (nextstrain-py39) root@Anarchy:/mnt/c/Users/Anarchy/Documents/Data_Science/NextStrain/rvf-nextstrain# augur ancestral --tree results/segments/L/refined/rvf_L_refined.tree --alignment results/segments/L/aligned/rvf_L_aligned.fasta --output-node-data results/segments/L/node_data/rvf_L_nt_muts.json

1. Traits
   (nextstrain-py39) root@Anarchy:/mnt/c/Users/Anarchy/Documents/Data_Science/NextStrain/rvf-nextstrain# augur traits --tree results/segments/L/refined/rvf_L_refined.tree --metadata results/segments/L/filtered/rvf_L_metadata_minimal.tsv --output results/segments/L/traits/rvf_L_traits.json --columns country

1. Export
   (nextstrain-py39) root@Anarchy:/mnt/c/Users/Anarchy/Documents/Data_Science/NextStrain/rvf-nextstrain# augur export v2 --tree results/segments/L/refined/rvf_L_refined.tree --node-data results/segments/L/refined/rvf_L_refined_node_data.json --node-data results/segments/L/node_data/rvf_L_nt_muts.json --node-data results/segments/L/traits/rvf_L_traits.json --metadata results/segments/L/filtered/rvf_L_metadata_minimal.tsv --output auspice/rvf_L.json
