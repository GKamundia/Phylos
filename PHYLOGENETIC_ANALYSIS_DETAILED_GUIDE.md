# Comprehensive Guide to Phylogenetic Analysis Rules (analysis.smk)

## 🔬 Overview

This file contains **4 Snakemake rules** that form the core of phylogenetic analysis in the Nextstrain pipeline. Think of it as a **recipe book** for building evolutionary family trees from viral genome sequences.

## 🎯 What This File Does (Simple Version)

```
Raw sequences → Aligned sequences → Tree → Time-calibrated tree → Mutations → Geographic traits
                      ↑              ↑            ↑                ↑              ↑
                   (input)      [RULE 1: tree]  [RULE 2: refine]  [RULE 3: ancestral]  [RULE 4: traits]
```

---

## 📊 Pipeline Context & Data Flow

### **Input Requirements:**

- **Aligned sequences**: DNA sequences lined up position by position
- **Metadata**: Information about when/where samples were collected
- **Configuration**: Settings that control how analysis runs

### **Output Products:**

- **Phylogenetic trees**: Family trees showing evolutionary relationships
- **Time estimates**: When mutations and lineage splits occurred
- **Ancestral sequences**: What ancient viruses likely looked like
- **Geographic patterns**: How viruses spread across locations

---

# 🌳 RULE 1: `tree` - Building the Basic Family Tree

## 🎯 **Simple Explanation**

This rule takes aligned viral sequences and figures out how they're related to each other, like building a family tree for viruses.

## 📥 **Inputs**

```python
input:
    alignment = get_alignment_input  # Function that returns aligned sequence file
```

**What `get_alignment_input` returns:**

- **Single segment mode**: `results/aligned.fasta`
- **Multi-segment mode**: `results/segments/{L,M,S}/aligned/rvf_{L,M,S}_aligned.fasta`

### **File Format Details:**

```fasta
>sequence1
ATGCGATCGATCG...
>sequence2
ATGCGATCGATCG...
>sequence3
ATGCGATCGATCG...
```

## 📤 **Outputs**

```python
output:
    tree = f"results/tree/{output_prefix}_tree.nwk" if segment_mode == "single"
           else "results/segments/{segment}/tree/rvf_{segment}_tree.nwk"
```

**Output Examples:**

- **Single**: `results/tree/rvf_tree.nwk`
- **Multi-segment**:
  - `results/segments/L/tree/rvf_L_tree.nwk`
  - `results/segments/M/tree/rvf_M_tree.nwk`
  - `results/segments/S/tree/rvf_S_tree.nwk`

### **Newick Format (.nwk):**

```
((sequence1:0.01,sequence2:0.02):0.005,sequence3:0.03);
```

- Numbers after `:` = branch lengths (evolutionary distance)
- Parentheses = groupings (who's related to whom)

## ⚙️ **Configuration Parameters - Detailed Analysis**

### **method** (default: "iqtree") - Tree Building Algorithm Selection

```yaml
# In config file:
tree:
  method: "iqtree" # Primary choice: "iqtree", "fasttree", "raxml"
```

#### **🔬 Detailed Method Comparison:**

| Method       | Algorithm          | Accuracy | Speed | Memory | Best Use Case               |
| ------------ | ------------------ | -------- | ----- | ------ | --------------------------- |
| **iqtree**   | Maximum Likelihood | ★★★★★    | ★★★☆☆ | ★★★☆☆  | Publication, final analysis |
| **fasttree** | Approximate ML     | ★★★☆☆    | ★★★★★ | ★★★★★  | Large datasets, screening   |
| **raxml**    | Maximum Likelihood | ★★★★☆    | ★★☆☆☆ | ★★☆☆☆  | Traditional alternative     |

#### **🎯 When to Use Each Method:**

**IQ-TREE (Recommended):**

```yaml
tree:
  method: "iqtree"
  iqtree_args: "-ninit 10 -n 10 -bb 1000"
```

- **Best for:** Final analysis, publication-quality trees
- **Strengths:**
  - Automatic model selection (ModelFinder)
  - Ultra-fast bootstrap (UFBoot)
  - Excellent likelihood optimization
  - Active development and support
- **When to use:** <1000 sequences, high-quality analysis needed
- **Runtime:** 10-60 minutes for typical datasets

**FastTree (Speed Priority):**

```yaml
tree:
  method: "fasttree"
  # Note: iqtree_args ignored for fasttree
```

- **Best for:** Large datasets (>1000 sequences), preliminary analysis
- **Strengths:**
  - Very fast execution
  - Low memory usage
  - Good for screening many datasets
- **When to use:** >1000 sequences, exploratory analysis
- **Runtime:** 1-10 minutes for large datasets

**RAxML (Traditional Choice):**

```yaml
tree:
  method: "raxml"
```

- **Best for:** When you need RAxML-specific features
- **Strengths:** Well-established, proven reliability
- **Limitations:** Slower than IQ-TREE, less automated

### **iqtree_args** (default: "-ninit 2 -n 2") - Fine-Tuning IQ-TREE Performance

```yaml
tree:
  iqtree_args: "-ninit 2 -n 2" # Conservative/fast settings
  # Alternative configurations shown below
```

#### **🔧 Comprehensive Parameter Breakdown:**

**Core Search Parameters:**

- **`-ninit X`**: Number of initial parsimony trees to start optimization

  - Default: `2` (very conservative)
  - Recommended: `10` (good balance)
  - Thorough: `100` (maximum accuracy)

- **`-n X`**: Number of independent tree searches
  - Default: `2` (minimal)
  - Recommended: `10` (publication quality)
  - Thorough: `50` (exhaustive search)

**Bootstrap Support Parameters:**

- **`-bb X`**: Ultra-fast bootstrap replicates

  - Standard: `1000` (recommended for publications)
  - Quick: `100` (acceptable for exploratory)
  - Thorough: `10000` (maximum confidence)

- **`-alrt X`**: SH-like approximate likelihood ratio test
  - Standard: `1000` (alternative to bootstrap)
  - Combines well with bootstrap: `-bb 1000 -alrt 1000`

#### **🎯 What Are Bootstrap Replicates? (Detailed Explanation)**

**📖 Simple Explanation:**
Bootstrap replicates measure **how confident we can be** in each branch of our phylogenetic tree. Think of it like asking "If I built this tree 1000 different ways, how often would I get the same result?"

**🔬 Technical Process:**

```
Original Dataset: 100 sequences, 1000 nucleotide positions
                 ↓
Bootstrap Round 1: Randomly sample 1000 positions WITH REPLACEMENT
                   → Some positions appear multiple times, others not at all
                   → Build tree from this "resampled" dataset
                 ↓
Bootstrap Round 2: Different random sample of 1000 positions
                   → Build another tree
                 ↓
... Repeat 1000 times ...
                 ↓
Compare Results: Count how often each branch appears across all 1000 trees
```

**🎲 Bootstrap Sampling Example:**

```
Original alignment positions: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Bootstrap replicate 1:        [3, 7, 7, 1, 9, 2, 5, 1, 4, 8]  ← Note duplicates!
Bootstrap replicate 2:        [2, 1, 6, 6, 3, 10, 4, 9, 7, 2]
Bootstrap replicate 3:        [8, 3, 1, 5, 5, 7, 9, 2, 6, 4]
... continue for 1000 rounds ...
```

**📊 Interpreting Bootstrap Values:**

| Bootstrap Value | Confidence Level | Interpretation                     | Publication Standard |
| --------------- | ---------------- | ---------------------------------- | -------------------- |
| **95-100%**     | Very High        | Strong support, trust this branch  | ✅ Reliable          |
| **80-94%**      | High             | Good support, likely correct       | ✅ Acceptable        |
| **70-79%**      | Moderate         | Some support, interpret cautiously | ⚠️ Questionable      |
| **50-69%**      | Low              | Weak support, unreliable           | ❌ Not supported     |
| **<50%**        | Very Low         | No meaningful support              | ❌ Discard           |

**🌳 Bootstrap Values on Trees:**

```
Example Newick tree with bootstrap values:
((species_A:0.01,species_B:0.02)95:0.005,(species_C:0.03,species_D:0.01)78:0.008);
                                   ↑                                      ↑
                              95% support                              78% support
```

**🔍 Why 1000 Replicates is Standard:**

**Statistical Reasoning:**

- **100 replicates**: Rough estimate, adequate for exploratory analysis
- **1000 replicates**: Industry standard, provides stable confidence estimates
- **10000 replicates**: Diminishing returns, minimal improvement over 1000

**Practical Considerations:**

```
Replicates    Runtime Impact    Confidence Precision    Recommended Use
100          +20% time         ±5-10%                  Quick screening
1000         +100% time        ±1-3%                   Publication standard
10000        +900% time        ±0.5-1%                 Ultra-high precision
```

**🚀 IQ-TREE's Ultra-Fast Bootstrap (UFBoot):**
Traditional bootstrap is very slow because it builds 1000 complete trees. IQ-TREE's UFBoot is much faster:

```
Traditional Bootstrap:
- Build 1000 complete trees from scratch
- Runtime: ~10-50x longer than tree building
- Memory: Standard

Ultra-Fast Bootstrap (UFBoot):
- Uses clever approximations and optimizations
- Runtime: ~2-5x longer than tree building
- Memory: ~25% more
- Accuracy: Nearly identical to traditional bootstrap
```

**⚙️ Practical Bootstrap Configuration:**

**For Different Research Scenarios:**

```yaml
# Quick exploratory analysis
iqtree_args: "-ninit 5 -n 5 -bb 100"
# → Fast, rough bootstrap estimates

# Standard publication
iqtree_args: "-ninit 10 -n 10 -bb 1000"
# → Reliable, widely accepted

# High-confidence study
iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
# → Bootstrap + additional statistical test

# Ultra-precise analysis
iqtree_args: "-ninit 50 -n 50 -bb 10000"
# → Maximum confidence, very slow
```

**🎯 When Bootstrap Values Matter Most:**

- **Manuscript submission**: Reviewers expect ≥80% support for key conclusions
- **Outbreak investigations**: High confidence needed for transmission chains
- **Evolutionary studies**: Branch support critical for phylogeographic conclusions
- **Vaccine design**: Tree topology must be highly reliable

**⚠️ Common Bootstrap Misconceptions:**

```
❌ "Bootstrap = probability the tree is correct"
✅ "Bootstrap = frequency this branch appears in resampled data"

❌ "Higher bootstrap always means better tree"
✅ "Bootstrap measures support for individual branches"

❌ "Low bootstrap means wrong tree"
✅ "Low bootstrap means insufficient data to resolve this relationship"
```

#### **📊 Configuration Examples by Use Case:**

**Pipeline/Development (Current Default):**

```yaml
iqtree_args: "-ninit 2 -n 2"
# Runtime: ~5-15 minutes
# Accuracy: Adequate for testing
# Memory: ~2-4GB
```

**Standard Research Analysis:**

```yaml
iqtree_args: "-ninit 10 -n 10 -bb 1000"
# Runtime: ~20-45 minutes
# Accuracy: Publication quality
# Memory: ~4-8GB
```

**High-Confidence Publication:**

```yaml
iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
# Runtime: ~45-90 minutes
# Accuracy: Maximum confidence
# Memory: ~6-12GB
```

**Large Dataset (>500 sequences):**

```yaml
iqtree_args: "-ninit 5 -n 5 -bb 100 -nt AUTO"
# Runtime: ~30-60 minutes
# Accuracy: Good for large datasets
# Memory: ~8-16GB
```

**Ultra-Thorough Analysis:**

```yaml
iqtree_args: "-ninit 100 -n 50 -bb 10000 -alrt 1000 -bnni"
# Runtime: ~2-6 hours
# Accuracy: Maximum possible
# Memory: ~8-20GB
```

#### **🚀 Performance Impact Matrix:**

| Setting              | Sequences | Runtime | RAM Usage | Accuracy Score |
| -------------------- | --------- | ------- | --------- | -------------- |
| `-ninit 2 -n 2`      | 100       | 5 min   | 2GB       | 85/100         |
| `-ninit 2 -n 2`      | 500       | 15 min  | 4GB       | 85/100         |
| `-ninit 10 -n 10`    | 100       | 15 min  | 3GB       | 92/100         |
| `-ninit 10 -n 10`    | 500       | 45 min  | 6GB       | 92/100         |
| `-ninit 20 -bb 1000` | 100       | 25 min  | 4GB       | 96/100         |
| `-ninit 20 -bb 1000` | 500       | 75 min  | 8GB       | 96/100         |

### **substitution_model** (default: "GTR") - DNA Evolution Models

```yaml
tree:
  substitution_model: "GTR" # General Time Reversible
  # Alternatives: "AUTO", "HKY", "K80", "JC69"
```

#### **🧬 Comprehensive Model Explanation:**

**Model Evolution Hierarchy:**

```
JC69 → F81 → K80 → HKY → TIM → TrN → TPM → GTR
Simple ←----------------------------------→ Complex
Fast   ←----------------------------------→ Accurate
Few Parameters ←--------------------------→ Many Parameters
```

#### **📚 Detailed Model Descriptions:**

**JC69 (Jukes-Cantor 1969):**

```yaml
substitution_model: "JC69"
```

- **Assumptions:** All substitutions equally likely, equal base frequencies
- **Parameters:** 1 (substitution rate)
- **Best for:** Very closely related sequences, quick analysis
- **Example:** Recently diverged viral strains
- **Runtime impact:** Fastest

**K80 (Kimura 1980):**

```yaml
substitution_model: "K80"
```

- **Assumptions:** Transitions ≠ transversions, equal base frequencies
- **Parameters:** 2 (transition rate, transversion rate)
- **Best for:** Moderately diverged sequences
- **Example:** Seasonal influenza analysis
- **Runtime impact:** Very fast

**HKY (Hasegawa-Kishino-Yano):**

```yaml
substitution_model: "HKY"
```

- **Assumptions:** Transitions ≠ transversions, unequal base frequencies
- **Parameters:** 5 (transition/transversion ratio + 3 base frequencies)
- **Best for:** Most viral datasets, good balance
- **Example:** HIV, hepatitis virus analysis
- **Runtime impact:** Fast

**GTR (General Time Reversible):**

```yaml
substitution_model: "GTR" # Current default
```

- **Assumptions:** All substitution types can have different rates
- **Parameters:** 9 (6 substitution rates + 3 base frequencies)
- **Best for:** Diverse datasets, publication-quality analysis
- **Example:** Cross-species viral analysis
- **Runtime impact:** Moderate

**AUTO (Automatic Model Selection):**

```yaml
substitution_model: "AUTO" # Best choice for important analysis
```

- **What it does:** IQ-TREE tests multiple models and picks the best
- **Models tested:** JC, F81, K80, HKY, TIM, TrN, TPM, GTR (+ variants)
- **Selection criteria:** AIC, BIC, or corrected AIC
- **Best for:** When you want optimal model choice
- **Runtime impact:** Slower (model testing overhead)

#### **🎯 Model Selection Guidelines:**

**By Dataset Characteristics:**

```
Closely related (>95% identity):     JC69 or K80
Moderately related (90-95%):         HKY
Distantly related (<90%):            GTR or AUTO
Mixed divergence levels:             AUTO
Publication analysis:                AUTO
Quick screening:                     GTR
Large datasets (>1000 seqs):        HKY or GTR
```

**By Research Goals:**

```
Exploratory analysis:                GTR (good default)
Hypothesis testing:                  AUTO (best model)
Comparative studies:                 AUTO (consistent methodology)
Real-time surveillance:              HKY (fast, adequate)
Manuscript preparation:              AUTO (reviewer expectations)
```

#### **⚡ Performance vs. Accuracy Trade-offs:**

| Model | Parameters | Relative Speed | Accuracy | Recommended Use         |
| ----- | ---------- | -------------- | -------- | ----------------------- |
| JC69  | 1          | 100%           | ★★☆☆☆    | Very similar sequences  |
| K80   | 2          | 95%            | ★★★☆☆    | Moderate diversity      |
| HKY   | 5          | 90%            | ★★★★☆    | Most viral datasets     |
| GTR   | 9          | 85%            | ★★★★★    | High diversity, default |
| AUTO  | Variable   | 60%            | ★★★★★    | Best possible choice    |

## 💻 **Resource Requirements - Detailed Specifications**

### **🔧 Configuration Structure:**

```python
threads: config["resources"].get("tree", {}).get("threads", 4)
resources:
    mem_mb = config["resources"].get("tree", {}).get("mem_mb", 8000)  # 8GB RAM
```

### **📊 Resource Scaling by Dataset Size:**

| Sequences | Recommended Threads | Recommended RAM | Typical Runtime |
| --------- | ------------------- | --------------- | --------------- |
| 50-100    | 2                   | 2-4 GB          | 5-15 min        |
| 100-300   | 4                   | 4-8 GB          | 10-30 min       |
| 300-500   | 6-8                 | 8-12 GB         | 20-60 min       |
| 500-1000  | 8-12                | 12-16 GB        | 30-120 min      |
| 1000+     | 16+                 | 16-32 GB        | 1-6 hours       |

### **⚙️ Hardware Optimization Guidelines:**

**CPU Considerations:**

```yaml
# Conservative (development/testing)
tree:
  threads: 2  # Good for laptops, shared systems

# Standard (research analysis)
tree:
  threads: 4  # Current default, good balance

# High-performance (production)
tree:
  threads: 8  # For dedicated servers

# Maximum (large datasets)
tree:
  threads: 16  # For high-end workstations
```

**Memory Scaling:**

```yaml
# Memory configurations by analysis type
resources:
  tree:
    # Minimal (small datasets, <100 sequences)
    mem_mb: 2000  # 2GB

    # Standard (current default)
    mem_mb: 8000  # 8GB - good for most analyses

    # Large datasets (>500 sequences)
    mem_mb: 16000  # 16GB

    # Very large (>1000 sequences)
    mem_mb: 32000  # 32GB
```

**� Resource Optimization Tips:**

- **IQ-TREE scales well** with multiple threads up to ~8-16 cores
- **Memory usage** increases with sequence count and model complexity
- **Bootstrap analysis** requires additional memory (~25% more)
- **Automatic model selection** needs extra RAM for model testing

## �🚀 **Execution Command - Complete Breakdown**

### **📋 Full Command Structure:**

```bash
python scripts/run_tree.py \
    --alignment {input.alignment} \
    --output {output.tree} \
    --log {log} \
    --method {params.method} \
    --threads {threads} \
    --substitution-model {params.substitution_model} \
    --iqtree-args "{params.iqtree_args}"
```

### **🔍 Parameter Expansion Examples:**

**Example 1 - Default Configuration:**

```bash
# With default parameters:
python scripts/run_tree.py \
    --alignment results/segments/L/aligned/rvf_L_aligned.fasta \
    --output results/segments/L/tree/rvf_L_tree.nwk \
    --log logs/tree_rvf_L.log \
    --method iqtree \
    --threads 4 \
    --substitution-model GTR \
    --iqtree-args "-ninit 2 -n 2"
```

**Example 2 - Publication Quality:**

```bash
# With enhanced parameters:
python scripts/run_tree.py \
    --alignment results/segments/L/aligned/rvf_L_aligned.fasta \
    --output results/segments/L/tree/rvf_L_tree.nwk \
    --log logs/tree_rvf_L.log \
    --method iqtree \
    --threads 8 \
    --substitution-model AUTO \
    --iqtree-args "-ninit 10 -n 10 -bb 1000 -alrt 1000"
```

**Example 3 - Fast Screening:**

```bash
# With FastTree for speed:
python scripts/run_tree.py \
    --alignment results/segments/L/aligned/rvf_L_aligned.fasta \
    --output results/segments/L/tree/rvf_L_tree.nwk \
    --log logs/tree_rvf_L.log \
    --method fasttree \
    --threads 4 \
    --substitution-model GTR \
    --iqtree-args ""  # Ignored for FastTree
```

### **🔧 What `scripts/run_tree.py` Does:**

**Internal Processing Steps:**

1. **Input Validation:**

   - Checks if alignment file exists and is non-empty
   - Counts number of sequences in alignment
   - Validates output directory paths

2. **Method-Specific Execution:**

   ```python
   if method == "iqtree":
       # Calls IQ-TREE with specified parameters
       cmd = f"iqtree -s {alignment} -nt {threads} -m {model} {iqtree_args}"
   elif method == "fasttree":
       # Calls FastTree (iqtree_args ignored)
       cmd = f"fasttree -gtr -nt < {alignment}"
   ```

3. **Error Handling:**

   - Logs execution details
   - Handles missing input files gracefully
   - Creates empty output files if no sequences available

4. **Output Generation:**
   - Converts output to standard Newick format
   - Ensures consistent file naming across methods
   - Logs runtime and resource usage

### **🗂️ File Path Resolution:**

**Single Segment Mode:**

```bash
# Input alignment:
{input.alignment} → results/aligned.fasta

# Output tree:
{output.tree} → results/tree/rvf_tree.nwk

# Log file:
{log} → logs/tree_rvf.log
```

**Multi-Segment Mode (RVF L, M, S segments):**

```bash
# Input alignments:
{input.alignment} → results/segments/L/aligned/rvf_L_aligned.fasta
                  → results/segments/M/aligned/rvf_M_aligned.fasta
                  → results/segments/S/aligned/rvf_S_aligned.fasta

# Output trees:
{output.tree} → results/segments/L/tree/rvf_L_tree.nwk
              → results/segments/M/tree/rvf_M_tree.nwk
              → results/segments/S/tree/rvf_S_tree.nwk

# Log files:
{log} → logs/tree_rvf_L.log
      → logs/tree_rvf_M.log
      → logs/tree_rvf_S.log
```

### **📊 Monitoring Execution:**

**Real-time Monitoring:**

```bash
# Watch log file during execution:
tail -f logs/tree_rvf_L.log

# Monitor system resources:
htop  # or Task Manager on Windows

# Check intermediate files:
ls -la results/segments/L/tree/
```

**Expected Log Output:**

```
Building tree from 245 aligned sequences
Using IQ-TREE method with GTR model
Starting tree search with 4 threads...
Checkpoint: Initial parsimony trees completed
Checkpoint: Likelihood optimization in progress...
Checkpoint: Bootstrap analysis starting...
Tree building completed successfully
Total runtime: 23.4 minutes
Peak memory usage: 6.2 GB
```

---

# ⏰ RULE 2: `refine` - Adding Time and Dates (Time Calibration)

## 🎯 **Simple Explanation**

This rule transforms your basic family tree into a **time-aware evolutionary timeline**. It's like taking a family tree that just shows "who's related to whom" and adding **birth dates, lifespans, and generational timing** to create a complete historical narrative.

## 🔄 **What Actually Happens (Step-by-Step)**

### **🌳 Before Refinement (Raw Tree):**

```
Basic tree from Rule 1:
- Shows evolutionary relationships only
- Branch lengths = number of mutations
- No timing information
- No dates on internal nodes

Example: "Virus A and B are closely related, but when did they diverge?"
```

### **⏰ After Refinement (Time-Calibrated Tree):**

```
Refined tree:
- Shows evolutionary relationships + timing
- Branch lengths = time periods (years/months)
- Internal nodes have estimated dates
- Molecular clock calibrated to real time

Example: "Virus A and B diverged around January 2023"
```

## 📊 **Visual Transformation Example**

**Raw Tree (mutations only):**

```
                    ┌─ Virus_A (collected: 2023-03-15)
    ┌───────────────┤
    │               └─ Virus_B (collected: 2023-03-20)
────┤
    │               ┌─ Virus_C (collected: 2023-02-10)
    └───────────────┤
                    └─ Virus_D (collected: 2023-01-05)

Branch lengths = mutations (0.001, 0.003, etc.)
Internal nodes = no dates
```

**Refined Tree (time-calibrated):**

```
                         ┌─ Virus_A (2023-03-15)
    ┌────────────────────┤
    │    ~2023-02-01     └─ Virus_B (2023-03-20)
────┤
    │    ~2022-12-15     ┌─ Virus_C (2023-02-10)
    └────────────────────┤
                         └─ Virus_D (2023-01-05)

Branch lengths = time periods (days/months)
Internal nodes = estimated divergence dates
```

## 📥 **Inputs - Detailed Breakdown**

```python
input:
    tree = "results/tree/{prefix}_tree.nwk",          # Raw tree from Rule 1
    alignment = get_alignment_input,                   # Same alignment file
    metadata = get_input_metadata                      # Sample collection dates/locations
```

### **🌳 Tree Input (from Rule 1):**

- **What it contains:** Basic phylogenetic relationships
- **Branch lengths:** Represent evolutionary distance (mutations per site)
- **What's missing:** No timing, no dates, no molecular clock calibration

### **🧬 Alignment Input:**

- **Purpose:** Calculate mutation rates along branches
- **How it's used:** Compare sequences to estimate evolutionary speed
- **Critical for:** Calibrating the molecular clock

### **📋 Metadata Input (The Key to Time Calibration):**

**Essential Metadata Columns:**

```tsv
strain          date        country     division    host
virus_001      2023-01-15   Kenya       Rift_Valley human
virus_002      2023-02-20   Tanzania    Arusha      bovine
virus_003      2023-03-10   Sudan       Kassala     human
virus_004      2022-11-05   Egypt       Cairo       human
virus_005      2023-04-02   Ethiopia    Addis_Ababa bovine
```

**🗓️ Date Column Requirements:**

- **Format:** YYYY-MM-DD (preferred) or YYYY-MM or YYYY
- **Completeness:** More dates = better time calibration
- **Date range:** Wider time span = more accurate molecular clock
- **Quality:** Accurate collection dates are crucial

**📍 Geographic Columns:**

- **Purpose:** Track virus movement through space and time
- **Examples:** country, division, region
- **Usage:** Combined with dates to infer transmission patterns

## 📤 **Outputs - What You Get**

```python
output:
    tree = "results/tree/{prefix}_refined.nwk",           # Time-calibrated tree
    node_data = "results/node_data/{prefix}_branch_lengths.json"  # Branch length data
```

### **🌳 Refined Tree (.nwk file):**

**What's Different:**

```newick
# Before (mutations):
((virus_A:0.001,virus_B:0.002):0.0005,virus_C:0.003);

# After (time-calibrated):
((virus_A:45.2,virus_B:23.7):67.8,virus_C:123.4);
# Numbers now represent days/years from root
```

**New Features:**

- **Branch lengths = time units** (days, months, years)
- **Proportional to calendar time** not just mutation count
- **Root positioned** at estimated common ancestor date

### **📊 Node Data (.json file):**

**Contains:**

```json
{
  "nodes": {
    "NODE_001": {
      "numdate": 2023.123, // Decimal year estimate
      "date": "2023-02-15", // Human-readable date
      "confidence": [2023.1, 2023.2], // Date confidence interval
      "clock_length": 0.0023 // Mutations on this branch
    }
  }
}
```

## ⚙️ **Configuration Parameters - Comprehensive Guide**

### **coalescent** (default: "opt") - Population Dynamics Model

```yaml
refine:
  coalescent: "opt" # Options: "opt", "const", "skyline"
```

#### **🧬 What is Coalescent Theory?**

**Simple Explanation:** Coalescent theory describes how lineages in a population "coalesce" (merge) backward in time as you trace ancestry. It's like asking: "If I go back in time, when do these viral lineages find their common ancestor?"

#### **📊 Detailed Options:**

**1. "opt" (Optimize) - RECOMMENDED:**

```yaml
coalescent: "opt"
```

- **What it does:** Automatically finds the best population model for your data
- **How it works:** Tests different population size assumptions and picks optimal
- **Best for:** Most analyses, when you don't know population dynamics
- **Pros:** Flexible, data-driven choice
- **Cons:** Slightly slower computation

**2. "const" (Constant Population):**

```yaml
coalescent: "const"
```

- **Assumption:** Population size stayed the same throughout time
- **When to use:** Stable endemic circulation, long-term surveillance data
- **Example:** Seasonal flu in established populations
- **Pros:** Fast computation, simple model
- **Cons:** May be unrealistic for outbreaks

**3. "skyline" (Variable Population):**

```yaml
coalescent: "skyline"
```

- **Assumption:** Population size can change over time
- **When to use:** Outbreaks, epidemics, population bottlenecks
- **Example:** COVID-19 pandemic, Ebola outbreaks
- **Pros:** Most realistic for epidemics
- **Cons:** Slower, needs larger datasets

#### **🎯 Choosing the Right Model:**

| Scenario                 | Recommended Setting | Why                        |
| ------------------------ | ------------------- | -------------------------- |
| **Unknown dynamics**     | `opt`               | Let algorithm decide       |
| **Stable endemic virus** | `const`             | Population likely stable   |
| **Outbreak/epidemic**    | `skyline`           | Population growing rapidly |
| **Long surveillance**    | `opt`               | Mixed dynamics likely      |
| **Small dataset**        | `const`             | Simple model more reliable |

### **date_inference** (default: "marginal") - Date Estimation Method

```yaml
refine:
  date_inference: "marginal" # Options: "marginal", "joint"
```

#### **🎯 What This Controls:**

How the algorithm estimates dates for internal nodes (ancestors) in your tree.

#### **📊 Detailed Comparison:**

**"marginal" (Independent Estimation) - DEFAULT:**

```yaml
date_inference: "marginal"
```

- **How it works:** Estimates each internal node's date independently
- **Process:** For each ancestor, finds most likely date given its descendants
- **Speed:** ⚡ Fast computation
- **Memory:** 💾 Low memory usage
- **Accuracy:** ⭐⭐⭐ Good for most purposes
- **Best for:** Large datasets, exploratory analysis, pipeline runs

**"joint" (Simultaneous Estimation):**

```yaml
date_inference: "joint"
```

- **How it works:** Estimates all internal node dates simultaneously
- **Process:** Considers interactions between all ancestral dates
- **Speed:** 🐌 Slower computation (2-5x longer)
- **Memory:** 💾💾 Higher memory usage
- **Accuracy:** ⭐⭐⭐⭐⭐ More accurate, especially for complex trees
- **Best for:** Publication-quality analysis, final results

#### **🔍 When Does the Difference Matter?**

**Marginal is Fine When:**

- Large, well-sampled datasets (>100 sequences)
- Good temporal coverage (samples across time range)
- Simple evolutionary patterns
- Exploratory analysis

**Joint is Better When:**

- Small datasets (<50 sequences)
- Poor temporal sampling (clumped dates)
- Complex evolutionary patterns
- Publication-quality analysis needed
- High precision required

#### **📈 Performance Comparison:**

| Dataset Size    | Marginal Runtime | Joint Runtime | Accuracy Difference |
| --------------- | ---------------- | ------------- | ------------------- |
| 50 sequences    | 2 min            | 5 min         | Significant         |
| 200 sequences   | 5 min            | 15 min        | Moderate            |
| 500 sequences   | 10 min           | 35 min        | Minimal             |
| 1000+ sequences | 20 min           | 90 min        | Negligible          |

### **clock_filter_iqd** (default: 4) - Outlier Detection and Removal

```yaml
refine:
  clock_filter_iqd: 4 # Interquartile distance multiplier
```

#### **🎯 Purpose: Molecular Clock Quality Control**

Some sequences evolve unusually fast or slow compared to the average rate. These "temporal outliers" can distort your molecular clock calibration.

#### **📊 How It Works:**

**Step 1: Calculate Divergence Rates**

```
For each sequence, calculate: mutations ÷ time_since_ancestor
```

**Step 2: Find the Normal Range**

```
Calculate quartiles of all rates:
Q1 = 25th percentile
Q3 = 75th percentile
IQR = Q3 - Q1 (interquartile range)
```

**Step 3: Define Outlier Boundaries**

```
Lower bound = Q1 - (clock_filter_iqd × IQR)
Upper bound = Q3 + (clock_filter_iqd × IQR)
```

**Step 4: Remove Outliers**

```
Keep sequences within bounds, remove others
```

#### **🔧 Practical Examples:**

**clock_filter_iqd: 2 (Strict Filtering):**

```yaml
clock_filter_iqd: 2
```

- **Effect:** Removes more sequences (aggressive filtering)
- **Keeps:** Only sequences with very consistent evolution rates
- **Risk:** May remove legitimate sequence variation
- **Use when:** Suspected poor-quality sequences, small outbreaks

**clock_filter_iqd: 4 (Standard - DEFAULT):**

```yaml
clock_filter_iqd: 4
```

- **Effect:** Balanced filtering (recommended for most cases)
- **Keeps:** Most sequences while removing clear outliers
- **Balance:** Good trade-off between quality and retention
- **Use when:** General analysis, unsure about data quality

**clock_filter_iqd: 6 (Permissive Filtering):**

```yaml
clock_filter_iqd: 6
```

- **Effect:** Removes fewer sequences (gentle filtering)
- **Keeps:** Almost all sequences, only extreme outliers removed
- **Risk:** May retain problematic sequences
- **Use when:** High-quality data, want to keep maximum sequences

#### **🚨 Sequences That Get Filtered (Examples):**

**Too Fast Evolution:**

- Laboratory-passaged viruses (adapted to cell culture)
- Sequencing errors creating false mutations
- Mislabeled collection dates (actually older than claimed)

**Too Slow Evolution:**

- Archival samples with degraded RNA
- Mixed infections (consensus of multiple strains)
- Mislabeled collection dates (actually newer than claimed)

#### **📊 Filtering Impact Example:**

```
Original dataset: 100 sequences
clock_filter_iqd: 4

Result: 94 sequences retained, 6 removed
Removed sequences:
- virus_047: evolving 8x faster than average (lab strain?)
- virus_083: evolving 12x slower than average (degraded sample?)
- virus_091: extreme date inconsistency
... (3 more outliers)
```

### **📈 Additional Optional Parameters**

#### **clock_rate** (optional):

```yaml
refine:
  clock_rate: 0.001 # Mutations per site per year
```

- **Purpose:** Fix evolutionary rate if known from literature
- **When to use:** Well-studied viruses with established rates
- **Example:** Influenza A: ~2.3×10⁻³, HIV: ~1.4×10⁻³

#### **clock_std_dev** (optional):

```yaml
refine:
  clock_std_dev: 0.0001 # Uncertainty in clock rate
```

- **Purpose:** Specify uncertainty in the fixed clock rate
- **Use with:** clock_rate parameter
- **Effect:** Allows some flexibility around the fixed rate

#### **root** (optional):

```yaml
refine:
  root: "best" # or specific strain name like "virus_001"
```

- **Purpose:** Control which sequence/node becomes the tree root
- **Options:**
  - `"best"`: Algorithm chooses optimal root
  - `"strain_name"`: Use specific sequence as root
- **When to specify:** Known ancestral sequence, specific rooting needed

### **Optional Clock Parameters:**

```yaml
refine:
  clock_rate: 0.001 # Mutations per site per year
  clock_std_dev: 0.0001 # Uncertainty in clock rate
  root: "best" # or specific strain name
```

## 💻 **Resource Requirements**

```python
threads: config["resources"].get("refine", {}).get("threads", 2)
resources:
    mem_mb = config["resources"].get("refine", {}).get("mem_mb", 4000)  # 4GB RAM
```

## 🚀 **Execution Command**

```bash
augur refine \
    --tree {input.tree} \
    --alignment {input.alignment} \
    --metadata {input.metadata} \
    --output-tree {output.tree} \
    --output-node-data {output.node_data} \
    --coalescent {params.coalescent} \
    --date-inference {params.date_inference} \
    --clock-filter-iqd {params.clock_filter_iqd} \
    --date-confidence \
    --keep-root
```

**Key Flags:**

- `--date-confidence`: Calculate uncertainty in date estimates
- `--keep-root`: Don't change the root position

---

# 🧬 RULE 3: `ancestral` - Reconstructing Ancient Sequences (Ancestral State Reconstruction)

## 🎯 **Simple Explanation**

This rule is like **genetic archaeology** - it reconstructs what ancient viral sequences looked like and maps exactly how viruses evolved over time. Think of it as building a complete evolutionary storybook that shows not just family relationships, but **what changed when** and **how**.

## 🔬 **What Actually Happens (Detailed Process)**

### **🧬 The Reconstruction Challenge:**

```
Your tree shows relationships:
- Virus A is related to Virus B
- They shared a common ancestor

But you want to know:
- What did that ancestor look like?
- Which mutations happened on which branches?
- When did specific changes occur?
```

### **🎯 The Solution - Ancestral Sequence Reconstruction:**

```
Given:
- Time-calibrated tree (with branch lengths as time)
- Modern sequences at tree tips
- Evolutionary model

Infer:
- Sequences at all internal nodes (ancestors)
- Mutations along each branch
- Timing of evolutionary changes
```

## 📊 **Visual Example: From Tips to Ancestors**

**What You Have (Tree Tips):**

```
                           ┌─ Virus_A: ATGCGATCGTAA...
    ┌─ Ancestor_1: ???     │
────┤                      └─ Virus_B: ATGCGTTCGTAA...
    │
    └─ Ancestor_2: ???     ┌─ Virus_C: ATGCGATCGTAG...
                           │
                           └─ Virus_D: ATGCGATCGTAG...
```

**What You Get (Complete Reconstruction):**

```
                           ┌─ Virus_A: ATGCGATCTTAA... (G→T at pos 6)
    ┌─ Ancestor_1: ATGCGATCGTAA     │
────┤                      └─ Virus_B: ATGCGTTCGTAA... (no change)
    │
    └─ Ancestor_2: ATGCGATCGTAG     ┌─ Virus_C: ATGCGATCGTAG... (no change)
           (A→G at pos 11)          │
                           └─ Virus_D: ATGCGATCGTAG... (no change)
```

## 📥 **Inputs - Detailed Breakdown**

```python
input:
    tree = "results/tree/{prefix}_refined.nwk",    # Time-calibrated tree
    alignment = get_alignment_input                # Aligned sequences
```

### **🌳 Time-Calibrated Tree (from Rule 2):**

**What it provides:**

- **Tree topology**: Who's related to whom
- **Branch lengths**: Time periods between nodes
- **Internal node dates**: When ancestors lived
- **Root position**: The most ancient common ancestor

**Critical for reconstruction:**

- Determines the **order of mutations** (which happened first)
- Provides **time constraints** for evolutionary inference
- Guides **likelihood calculations** for ancestral states

### **🧬 Aligned Sequences (Original Input):**

**What it provides:**

- **Modern sequences** at tree tips (what we observe today)
- **Sequence positions** where mutations can occur
- **Variation patterns** that inform evolutionary models

**How it's used:**

- **Reference points** for working backward to ancestors
- **Mutation rate estimation** for different sites
- **Evolutionary constraint** information (some positions change rarely)

### **� Data Requirements for Quality Reconstruction:**

| Data Aspect              | Minimum   | Good       | Excellent | Impact                    |
| ------------------------ | --------- | ---------- | --------- | ------------------------- |
| **Sequence Length**      | >500 bp   | >1000 bp   | >3000 bp  | Mutation signal strength  |
| **Number of Sequences**  | >10       | >50        | >100      | Phylogenetic resolution   |
| **Temporal Span**        | >6 months | >2 years   | >5 years  | Clock calibration         |
| **Geographic Diversity** | 1 region  | 3+ regions | Global    | Evolutionary completeness |

## �📤 **Outputs - What You Get**

```python
output:
    node_data = "results/node_data/{prefix}_nt_muts.json"  # Mutations on each branch
```

### **🔍 Complete Output Structure:**

```json
{
  "version": "v2",
  "meta": {
    "updated": "2023-06-19",
    "title": "Ancestral sequence reconstruction"
  },
  "nodes": {
    "NODE_0000001": {
      "sequence": "ATGCGATCGTAAGCCTTAG...",
      "muts": ["A123G", "T456C", "G789A"],
      "aa_muts": {
        "ORF1": ["K41E", "L152P"],
        "ORF2": ["V23I"]
      }
    },
    "virus_sample_001": {
      "muts": ["G234A"],
      "aa_muts": {
        "ORF1": ["A78T"]
      }
    }
  }
}
```

### **📋 Key Output Components:**

#### **1. Ancestral Sequences:**

```json
"sequence": "ATGCGATCGTAAGCCTTAG..."
```

- **Complete genome sequence** for each internal node
- **Reconstructed based on descendants** and evolutionary model
- **Most likely sequence** given the evidence

#### **2. Nucleotide Mutations (muts):**

```json
"muts": ["A123G", "T456C", "G789A"]
```

- **Format**: `[original][position][new]`
- **A123G**: Adenine at position 123 changed to Guanine
- **Lists all changes** that occurred on the branch leading to this node

#### **3. Amino Acid Mutations (aa_muts):**

```json
"aa_muts": {
  "ORF1": ["K41E", "L152P"],
  "ORF2": ["V23I"]
}
```

- **Protein-level changes** organized by gene/ORF
- **K41E**: Lysine at position 41 changed to Glutamic acid
- **Critical for understanding** functional impacts

#### **4. Branch-Specific Information:**

- **Each branch** gets its own mutation list
- **Mutations are mapped** to specific time periods
- **Allows tracking** of evolutionary pressures

## 🔬 **Reconstruction Algorithms Explained**

### **🎯 The Core Challenge:**

Given sequences at tree tips and a tree structure, what's the most likely sequence at each internal node?

### **📊 Maximum Likelihood Approach:**

```
For each position in the sequence:
1. Consider all possible nucleotides (A, T, G, C) at ancestor
2. Calculate likelihood of observing descendant sequences
3. Choose nucleotide with highest likelihood
4. Account for branch lengths and substitution model
```

### **🧮 Mathematical Foundation:**

```
L(ancestor) = P(descendants | ancestor, tree, model)

Where:
- L = likelihood of ancestral state
- P = probability function
- descendants = observed sequences at tips
- tree = topology and branch lengths
- model = substitution model (GTR, HKY, etc.)
```

## ⚙️ **Configuration Parameters - Comprehensive Guide**

### **inference** (default: "joint") - Reconstruction Strategy

```yaml
ancestral:
  inference: "joint" # Options: "joint", "marginal"
```

#### **🔍 What This Controls:**

How the algorithm coordinates ancestral sequence reconstruction across the entire tree.

#### **📊 Detailed Method Comparison:**

**"joint" (Simultaneous Reconstruction) - RECOMMENDED:**

```yaml
inference: "joint"
```

**How it works:**

1. **Considers entire tree** simultaneously
2. **Accounts for dependencies** between ancestral nodes
3. **Optimizes consistency** across all reconstructions
4. **Uses dynamic programming** for efficiency

**Process:**

```
Bottom-up phase:
- Start at tree tips (known sequences)
- Calculate likelihoods for each internal node
- Consider all possible ancestral states

Top-down phase:
- Start at root with most likely state
- Propagate optimal states down the tree
- Ensure consistency between ancestors and descendants
```

**Advantages:**

- ✅ **Most accurate** - considers global tree structure
- ✅ **Statistically rigorous** - proper likelihood framework
- ✅ **Consistent reconstructions** - ancestors fit with descendants
- ✅ **Publication quality** - standard method in literature

**Disadvantages:**

- ⚠️ **Slower computation** - especially for large trees
- ⚠️ **Higher memory usage** - stores intermediate calculations

**"marginal" (Independent Reconstruction):**

```yaml
inference: "marginal"
```

**How it works:**

1. **Each node reconstructed independently**
2. **Uses only immediate descendants**
3. **Faster but less coordinated**
4. **No global optimization**

**Process:**

```
For each internal node:
- Look only at direct descendants
- Find most likely ancestral state
- Ignore broader tree context
- No consistency checking
```

**Advantages:**

- ✅ **Faster computation** - parallelizable
- ✅ **Lower memory usage** - no global storage
- ✅ **Good for exploration** - quick results

**Disadvantages:**

- ⚠️ **Less accurate** - misses global patterns
- ⚠️ **Potential inconsistencies** - ancestors may not fit descendants
- ⚠️ **Not publication standard** - less rigorous

#### **🎯 When to Use Each Method:**

| Scenario                        | Recommended          | Reasoning                           |
| ------------------------------- | -------------------- | ----------------------------------- |
| **Final analysis**              | `joint`              | Maximum accuracy needed             |
| **Publication**                 | `joint`              | Scientific rigor required           |
| **Small datasets (<100 seqs)**  | `joint`              | Manageable computation              |
| **Large datasets (>1000 seqs)** | `marginal` → `joint` | Start fast, refine later            |
| **Exploratory analysis**        | `marginal`           | Quick feedback needed               |
| **Real-time surveillance**      | `marginal`           | Speed priority                      |
| **Functional analysis**         | `joint`              | Accuracy crucial for interpretation |

#### **📈 Performance Comparison:**

| Dataset Size    | Joint Runtime | Marginal Runtime | Accuracy Difference |
| --------------- | ------------- | ---------------- | ------------------- |
| 50 sequences    | 30 seconds    | 10 seconds       | High                |
| 200 sequences   | 2 minutes     | 30 seconds       | Moderate            |
| 500 sequences   | 8 minutes     | 1.5 minutes      | Low                 |
| 1000+ sequences | 20+ minutes   | 3 minutes        | Minimal             |

### **📊 Advanced Configuration Options**

#### **Handling Ambiguous Sites:**

```yaml
ancestral:
  inference: "joint"
  keep_ambiguous: false # How to handle uncertain reconstructions
```

#### **Output Customization:**

```yaml
ancestral:
  output_sequences: true # Include full ancestral sequences
  output_vcf: false # Generate VCF format mutations
  infer_gtr: true # Infer GTR model parameters
```

## 💻 **Resource Requirements - Scaling Guidelines**

### **🔧 Configuration Structure:**

```python
threads: config["resources"].get("ancestral", {}).get("threads", 1)
resources:
    mem_mb = config["resources"].get("ancestral", {}).get("mem_mb", 2000)  # 2GB RAM
```

### **📊 Resource Scaling by Dataset:**

| Sequences | Genome Length | Joint Memory | Joint Runtime | Marginal Memory | Marginal Runtime |
| --------- | ------------- | ------------ | ------------- | --------------- | ---------------- |
| 50        | 3kb           | 500 MB       | 30 sec        | 200 MB          | 10 sec           |
| 100       | 3kb           | 800 MB       | 1 min         | 300 MB          | 15 sec           |
| 200       | 3kb           | 1.5 GB       | 2 min         | 500 MB          | 30 sec           |
| 500       | 3kb           | 3 GB         | 8 min         | 1 GB            | 1.5 min          |
| 100       | 10kb          | 2 GB         | 3 min         | 800 MB          | 45 sec           |
| 100       | 30kb          | 6 GB         | 10 min        | 2 GB            | 2 min            |

### **🎯 Optimization Tips:**

**For Large Datasets:**

```yaml
ancestral:
  inference: "marginal" # Start with faster method
# Then optionally run joint for final analysis
```

**For Long Genomes:**

```yaml
resources:
  ancestral:
    mem_mb: 8000 # Increase memory for long sequences
```

**For High-Resolution Analysis:**

```yaml
ancestral:
  inference: "joint" # Maximum accuracy
  threads: 1 # Single-threaded but thorough
```

## 🚀 **Execution Command - Complete Breakdown**

### **📋 Full Command Structure:**

```bash
augur ancestral \
    --tree {input.tree} \
    --alignment {input.alignment} \
    --output-node-data {output.node_data} \
    --inference {params.inference} \
    --infer-ambiguous \
    --keep-root
```

### **🔧 Parameter Explanations:**

**Core Parameters:**

- `--tree`: Time-calibrated tree with branch lengths
- `--alignment`: Original sequence alignment
- `--output-node-data`: JSON file with ancestral sequences and mutations
- `--inference`: Method (joint/marginal)

**Quality Control:**

- `--infer-ambiguous`: Handle uncertain positions appropriately
- `--keep-root`: Maintain root sequence from tree

### **📊 Expected Runtime by Configuration:**

**Default (joint inference, 100 sequences, 3kb genome):**

```bash
# Expected output:
Reconstructing ancestral sequences
Reading alignment...
Reading tree...
Inferring ancestral sequences using joint reconstruction
Analyzing 3247 variable sites across 100 sequences
Completed ancestral sequence reconstruction
Runtime: 1.2 minutes
Memory peak: 1.1 GB
```

## 🎯 **Interpreting Results - Practical Guide**

### **🔍 Understanding Mutation Patterns:**

**High-Confidence Mutations:**

```json
"muts": ["A123G", "T456C"]
```

- **Single nucleotide changes** with high likelihood
- **Reliable for downstream analysis**
- **Can be used for phylogeographic inference**

**Complex Changes:**

```json
"muts": ["A123G", "T124C", "G125A"]
```

- **Clustered mutations** may indicate recombination
- **Selection pressure** at specific genomic regions
- **Potential sequencing artifacts** to investigate

### **📈 Biological Interpretation:**

**Synonymous vs. Non-synonymous:**

```json
"aa_muts": {
  "ORF1": [],           // No amino acid changes (synonymous)
  "ORF2": ["K41E"]      // Functional change (non-synonymous)
}
```

**Evolutionary Pressure Indicators:**

- **Many synonymous changes**: Neutral evolution
- **Many non-synonymous changes**: Selection or relaxed constraint
- **Clustered changes**: Potential recombination or strong selection

## ⚠️ **Common Issues and Solutions**

### **🚨 Low-Quality Reconstructions:**

**Symptoms:**

- Many ambiguous positions (N's in sequences)
- Inconsistent mutation patterns
- High uncertainty scores

**Solutions:**

1. **Improve tree quality** (better bootstrapping in Rule 1)
2. **Add more sequences** for better phylogenetic signal
3. **Check alignment quality** for systematic errors
4. **Use joint inference** for better accuracy

### **🚨 Memory Issues:**

**Error:** `OutOfMemoryError during ancestral reconstruction`

**Solutions:**

```yaml
# Reduce memory usage:
ancestral:
  inference: "marginal"

# Or increase memory allocation:
resources:
  ancestral:
    mem_mb: 8000
```

### **🚨 Runtime Issues:**

**Problem:** Reconstruction taking too long

**Solutions:**

1. **Start with marginal inference** for quick results
2. **Subset sequences** for initial analysis
3. **Check for very long branches** that slow computation

## 🎯 **Quality Assessment**

### **📊 Reconstruction Confidence Metrics:**

**High-Quality Indicators:**

- **Low ambiguity**: Few 'N' characters in ancestral sequences
- **Consistent patterns**: Similar mutation rates across branches
- **Biological plausibility**: Reasonable amino acid changes

**Warning Signs:**

- **High ambiguity**: Many uncertain positions
- **Rate variation**: Extremely fast/slow branches
- **Nonsense mutations**: Stop codons in essential genes

### **🔍 Validation Strategies:**

1. **Cross-validation**: Remove sequences and test reconstruction accuracy
2. **Sensitivity analysis**: Compare joint vs. marginal results
3. **Biological validation**: Check if mutations are functionally plausible

This comprehensive ancestral sequence reconstruction transforms your time-calibrated tree into a complete evolutionary narrative, showing not just relationships but the specific genetic changes that occurred throughout viral evolution.

---

# 🌍 RULE 4: `traits` - Geographic and Metadata Reconstruction

## 🎯 **Simple Explanation**

Figures out where ancestral viruses came from - like mapping your family tree onto a world map.

## 📥 **Inputs**

```python
input:
    tree = "results/tree/{prefix}_refined.nwk",    # Time-calibrated tree
    metadata = get_input_metadata                  # Geographic/trait information
```

## 📤 **Outputs**

```python
output:
    node_data = "results/node_data/{prefix}_traits.json"  # Ancestral trait reconstructions
```

**Output Contains:**

- **Ancestral locations** for internal nodes
- **Migration events** between geographic regions
- **Confidence intervals** for trait assignments

**Example JSON Structure:**

```json
{
  "nodes": {
    "internal_node_1": {
      "country": {
        "value": "Kenya",
        "confidence": 0.85
      }
    }
  }
}
```

## ⚙️ **Configuration Parameters**

### **columns** (default: ["country"])

```yaml
traits:
  columns: ["country", "division", "host"] # Which metadata columns to reconstruct
```

**Common Trait Types:**

- **Geographic**: country, division, region
- **Host**: human, animal, vector
- **Temporal**: year, season
- **Clinical**: severity, symptom

**Example Metadata Columns:**

```tsv
strain      country    division      host     season
virus_001   Kenya      Rift_Valley   human    wet
virus_002   Tanzania   Arusha        bovine   dry
```

## 💻 **Resource Requirements**

```python
threads: config["resources"].get("traits", {}).get("threads", 1)
resources:
    mem_mb = config["resources"].get("traits", {}).get("mem_mb", 2000)  # 2GB RAM
```

## 🚀 **Execution Command**

```bash
augur traits \
    --tree {input.tree} \
    --metadata {input.metadata} \
    --metadata-id-columns Accession \
    --output {output.node_data} \
    --columns {params.columns}
```

**Key Flags:**

- `--metadata-id-columns Accession`: Which column contains sequence IDs
- `--columns`: Which traits to reconstruct

---

# 🔄 Complete Pipeline Flow

## **Sequential Dependencies:**

```
1. tree      → 2. refine    → 3. ancestral
   ↓             ↓             ↓
   Basic tree    + Time       + Mutations
                 ↓
              4. traits
                 ↓
              + Geography
```

## **Resource Summary:**

| Rule      | CPU Threads | RAM (MB) | Typical Runtime |
| --------- | ----------- | -------- | --------------- |
| tree      | 4           | 8000     | 10-60 min       |
| refine    | 2           | 4000     | 5-20 min        |
| ancestral | 1           | 2000     | 2-10 min        |
| traits    | 1           | 2000     | 1-5 min         |

## **Total Resource Requirements:**

- **Peak RAM**: ~8GB (during tree building)
- **Total Time**: 20-95 minutes (depending on dataset size)
- **Disk Space**: ~500MB-2GB (depending on number of sequences)

---

# 🛠️ Configuration Examples

## **Fast Pipeline (Testing/Development):**

```yaml
tree:
  method: "fasttree"
  substitution_model: "GTR"

refine:
  coalescent: "const"
  date_inference: "marginal"

ancestral:
  inference: "marginal"

traits:
  columns: ["country"]
```

## **Publication Quality:**

```yaml
tree:
  method: "iqtree"
  iqtree_args: "-ninit 10 -n 10 -bb 1000 -alrt 1000"
  substitution_model: "AUTO"

refine:
  coalescent: "opt"
  date_inference: "joint"
  clock_filter_iqd: 3

ancestral:
  inference: "joint"

traits:
  columns: ["country", "division", "host"]
```

## **Memory-Constrained Environment:**

```yaml
resources:
  tree:
    threads: 2
    mem_mb: 4000
  refine:
    threads: 1
    mem_mb: 2000
  ancestral:
    threads: 1
    mem_mb: 1000
  traits:
    threads: 1
    mem_mb: 1000
```

---

# 🎯 Key Takeaways

## **What This File Accomplishes:**

1. **Builds evolutionary relationships** from sequence data
2. **Adds temporal dimension** using collection dates
3. **Reconstructs ancestral states** at internal nodes
4. **Maps geographic/trait evolution** through time

## **When to Modify These Rules:**

- **Change tree method** for different speed/accuracy tradeoffs
- **Adjust resource limits** based on available hardware
- **Modify trait columns** based on your specific research questions
- **Tune filtering parameters** for different data quality levels

## **Common Issues & Solutions:**

- **Out of memory**: Reduce `mem_mb` values or use fewer threads
- **Too slow**: Switch to `fasttree` or reduce `iqtree_args`
- **Poor date estimates**: Check metadata date formats and completeness
- **Missing traits**: Verify metadata column names match configuration

This pipeline forms the **analytical backbone** of Nextstrain's phylogenetic inference, transforming raw sequence data into rich evolutionary narratives with temporal and spatial context.

---

# 📚 **What is "Publication Quality Analysis"? - Comprehensive Guide**

## 🎯 **Definition and Context**

**Publication Quality Analysis** refers to phylogenetic analyses that meet the **rigorous standards expected by peer-reviewed scientific journals**. This means your analysis must be sufficiently robust, accurate, and well-documented to withstand scientific scrutiny and support your research conclusions.

## 🔬 **Key Requirements for Publication Quality**

### **1. Statistical Rigor**

```yaml
# Publication-quality configuration:
tree:
  method: "iqtree"
  iqtree_args: "-ninit 10 -n 10 -bb 1000 -alrt 1000"
  substitution_model: "AUTO"
```

**What this achieves:**

- **Multiple tree searches** (`-ninit 10 -n 10`) reduce chance of getting stuck in local optima
- **Bootstrap support** (`-bb 1000`) provides confidence measures for every branch
- **Additional statistical tests** (`-alrt 1000`) give independent validation
- **Optimal model selection** (`AUTO`) ensures best evolutionary model for your data

### **2. Reproducibility Standards**

```yaml
# All parameters explicitly documented:
refine:
  coalescent: "opt" # Document WHY you chose this
  date_inference: "joint" # Explain your methodology choice
  clock_filter_iqd: 3 # Justify your filtering threshold
  clock_rate: 0.001 # If known, cite literature source
```

### **3. Comprehensive Analysis Depth**

```yaml
ancestral:
  inference: "joint" # More accurate than "marginal"

traits:
  columns: ["country", "division", "host", "date"] # Multiple relevant traits
```

## 📊 **Publication vs. Exploratory Analysis Comparison**

| Aspect                   | Exploratory/Pipeline | Publication Quality            | Impact                  |
| ------------------------ | -------------------- | ------------------------------ | ----------------------- |
| **Bootstrap Replicates** | 0-100                | 1000+                          | Confidence in results   |
| **Tree Searches**        | 2-5                  | 10-50                          | Finding optimal tree    |
| **Model Selection**      | Fixed (GTR)          | AUTO                           | Best evolutionary model |
| **Statistical Tests**    | None                 | Multiple (bootstrap + SH-aLRT) | Robust support          |
| **Runtime**              | 5-30 min             | 1-6 hours                      | Investment in accuracy  |
| **Documentation**        | Minimal              | Comprehensive                  | Reproducibility         |
| **Validation**           | Internal use         | Peer review ready              | Scientific credibility  |

## 🎖️ **Standards by Research Field**

### **🦠 Virology/Epidemiology Publications:**

```yaml
# Minimum standards for virology journals:
tree:
  iqtree_args: "-ninit 10 -n 10 -bb 1000"
  # Branch support ≥80% required for epidemiological conclusions
  # Bootstrap ≥1000 replicates for transmission chains
```

### **🧬 Evolutionary Biology:**

```yaml
# Higher standards for evolution journals:
tree:
  iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
  substitution_model: "AUTO"
  # Multiple statistical tests required
  # Model justification essential
```

### **🏥 Clinical/Medical Research:**

```yaml
# Strictest standards for medical journals:
tree:
  iqtree_args: "-ninit 50 -n 50 -bb 10000 -alrt 1000"
  # Ultra-high confidence needed for clinical decisions
  # Extensive sensitivity analysis required
```

## 🔍 **Reviewer Expectations by Journal Tier**

### **🌟 Top-Tier Journals (Nature, Science, Cell):**

**Requirements:**

- **Bootstrap ≥1000** with support values clearly displayed
- **Multiple statistical tests** (bootstrap + SH-aLRT + others)
- **Model selection justification** (why AUTO or specific model chosen)
- **Sensitivity analysis** (how robust are results to parameter changes?)
- **Methodological comparison** (why IQ-TREE vs. alternatives?)
- **Code and data availability** (full reproducibility)

**Example Configuration:**

```yaml
tree:
  method: "iqtree"
  iqtree_args: "-ninit 50 -n 50 -bb 1000 -alrt 1000 -bnni"
  substitution_model: "AUTO"
refine:
  coalescent: "opt"
  date_inference: "joint"
ancestral:
  inference: "joint"
```

### **🎯 Field-Specific Journals (Virology, Mol Biol Evol):**

**Requirements:**

- **Bootstrap ≥1000** for key conclusions
- **Appropriate evolutionary model** (justified choice)
- **Clear methodology description**
- **Statistical support thresholds** defined
- **Biological interpretation** of results

**Example Configuration:**

```yaml
tree:
  method: "iqtree"
  iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
  substitution_model: "AUTO"
```

### **📋 Specialized Journals (Epidemiology, Bioinformatics):**

**Requirements:**

- **Bootstrap ≥100-1000** depending on application
- **Method comparison** when novel approaches used
- **Parameter sensitivity** analysis
- **Computational efficiency** considerations

## ⚠️ **Common Reasons for Rejection**

### **🚫 Insufficient Statistical Support:**

```yaml
# REJECTED configuration:
tree:
  iqtree_args: "-ninit 2 -n 2" # Too few searches
  # No bootstrap analysis
  # Fixed model without justification
```

**Reviewer comment:** _"The phylogenetic analysis lacks sufficient statistical support. Branch support values should be provided using bootstrap analysis with ≥1000 replicates."_

### **🚫 Poor Methodology Documentation:**

```yaml
# REJECTED approach:
# Using defaults without explanation
# No parameter justification
# Missing model selection rationale
```

**Reviewer comment:** _"The methods section does not adequately describe the phylogenetic reconstruction methodology. Please provide rationale for parameter choices and model selection."_

### **🚫 Inadequate Validation:**

```yaml
# REJECTED validation:
# Single statistical test
# No sensitivity analysis
# No comparison with alternative methods
```

## ✅ **Publication-Ready Checklist**

### **📝 Methodology Requirements:**

- [ ] **Bootstrap analysis** with ≥1000 replicates
- [ ] **Multiple tree searches** (≥10 independent runs)
- [ ] **Model selection** documented and justified
- [ ] **Statistical thresholds** clearly defined (e.g., ≥80% bootstrap)
- [ ] **Software versions** specified
- [ ] **Parameter choices** explained and justified

### **📊 Results Requirements:**

- [ ] **Branch support values** displayed on trees
- [ ] **Confidence intervals** for date estimates (if temporal analysis)
- [ ] **Statistical significance** of key relationships
- [ ] **Sensitivity analysis** demonstrating robustness
- [ ] **Alternative hypotheses** tested and excluded

### **📚 Documentation Requirements:**

- [ ] **Complete methods description** in manuscript
- [ ] **Supplementary files** with detailed parameters
- [ ] **Data availability** statement
- [ ] **Code availability** (scripts, configuration files)
- [ ] **Reproducibility instructions**

## 🚀 **Transitioning from Pipeline to Publication**

### **Step 1: Enhance Statistical Rigor**

```yaml
# From pipeline settings:
tree:
  iqtree_args: "-ninit 2 -n 2"

# To publication settings:
tree:
  iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
```

### **Step 2: Add Model Selection**

```yaml
# From fixed model:
tree:
  substitution_model: "GTR"

# To optimal model:
tree:
  substitution_model: "AUTO"
```

### **Step 3: Increase Analysis Depth**

```yaml
# From quick analysis:
refine:
  date_inference: "marginal"
ancestral:
  inference: "marginal"

# To thorough analysis:
refine:
  date_inference: "joint"
ancestral:
  inference: "joint"
```

### **Step 4: Document Everything**

```yaml
# Add comprehensive documentation:
tree:
  method: "iqtree" # Cite: Nguyen et al. 2015
  iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
  substitution_model: "AUTO" # ModelFinder: Kalyaanamoorthy et al. 2017

refine:
  coalescent: "opt" # TreeTime: Sagulenko et al. 2018
  clock_filter_iqd: 3 # Removes temporal outliers >3 IQD
```

## 💰 **Cost-Benefit Analysis**

| Investment              | Benefit                    | Journal Acceptance        |
| ----------------------- | -------------------------- | ------------------------- |
| **+2-10x runtime**      | Robust statistical support | ✅ Reviewer confidence    |
| **+Higher memory**      | Comprehensive analysis     | ✅ Methodological rigor   |
| **+Documentation time** | Full reproducibility       | ✅ Scientific credibility |
| **+Validation effort**  | Defendable conclusions     | ✅ Publication success    |

## 🎯 **Bottom Line for Your RVF Analysis**

**Current pipeline configuration** is perfect for:

- ✅ Development and testing
- ✅ Exploratory analysis
- ✅ Real-time surveillance
- ✅ Internal reporting

**For publication**, upgrade to:\*\*

```yaml
# Publication-ready RVF analysis:
tree:
  method: "iqtree"
  iqtree_args: "-ninit 20 -n 20 -bb 1000 -alrt 1000"
  substitution_model: "AUTO"

refine:
  coalescent: "opt"
  date_inference: "joint"
  clock_filter_iqd: 3

ancestral:
  inference: "joint"

traits:
  columns: ["country", "division", "host"]
```

This ensures your Rift Valley Fever phylogenetic analysis meets the standards expected by virology and epidemiology journals for peer-reviewed publication.
