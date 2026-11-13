# Data Science for Biology — Workshop Environment

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=main)

> **Multiomics, Microbiome and Host Gene Expression**  
> Nairobi, 9–13 November 2026

A cloud-based learning environment for integrated analysis of microbiome and host transcriptomics data. This workshop teaches participants to apply the same visualization and statistical pipelines to both microbial community and human gene expression data using R, Quarto, and the microViz/DESeq2 ecosystem.

---

## Table of Contents

- [Overview](#overview)
- [Getting Started](#getting-started)
- [Project Structure](#project-structure)
- [Data Architecture](#data-architecture)
- [Analysis Workflow](#analysis-workflow)
- [Technical Methods](#technical-methods)
  - [Data Transformations](#data-transformations)
  - [Distance Metrics](#distance-metrics)
  - [Ordination Methods](#ordination-methods)
  - [Differential Testing with DESeq2](#differential-testing-with-deseq2)
  - [Enrichment Analysis](#enrichment-analysis)
- [The Unified Pipeline Concept](#the-unified-pipeline-concept)
- [Working with the AI Assistant](#working-with-the-ai-assistant)
- [Group Project](#group-project)
- [Important Considerations](#important-considerations)
- [References](#references)

---

## Overview

This workshop demonstrates that **microbiome and transcriptomics data share the same mathematical structure**: both are sparse count matrices where features (species or genes) are rows and samples are columns. By recognizing this, we can apply a single analytical framework to both data types.

### Learning Objectives

1. Understand the anatomy of count tables and their statistical properties
2. Perform ordination and distance-based analyses on multiomics data
3. Conduct differential abundance/expression testing with DESeq2
4. Run pathway and gene set enrichment analyses
5. Integrate findings across microbiome and host response data

---

## Getting Started

### Prerequisites

- A web browser
- A free GitHub account
- An [OpenRouter](https://openrouter.ai/) API key for the AI assistant

### Setup Instructions

1. **Create your copy**: Click **"Use this template"** → **"Create a new repository"**

2. **Add your API key**: Go to GitHub **Settings → Codespaces → Secrets**
   - Create a secret named `OPENROUTER_API_KEY`
   - Grant it access to your new repository

3. **Launch the Codespace**: On your repository, click **"Code" → "Codespaces" → "Create codespace on main"**

4. **Verify your setup**: Open `notebooks/00-hello-r.qmd` and run the code cell
   > ⏱️ The first build takes several minutes. Do this **before the workshop starts**.

5. **Start the AI assistant**: Open a terminal and type `pi`

📖 [Full step-by-step walkthrough on the workshop website](https://elsherbini.github.io/data-science-for-biology/materials/2-workshop2/getting-started.html)

---

## Project Structure

```
.
├── notebooks/                    # Learning modules (Quarto documents)
│   ├── 00-hello-r.qmd           # Environment check
│   ├── 02-count-table-anatomy.qmd   # Module 2: Understanding count data
│   ├── 03-ordination.qmd        # Module 3: PCA, PCoA, distance metrics
│   ├── 04-differential.qmd      # Module 4: DESeq2 differential testing
│   └── 05-enrichment.qmd        # Module 5: Pathway analysis
│
├── group_project/               # Group project workspace
│   ├── 01-qc.qmd               # Quality control report
│   ├── 02-*.qmd                # Analysis documents
│   └── README.md               # Project prompts and workflow
│
├── data/2-workshop2/            # Multiomics datasets
│   ├── instructional_dataset/   # Pregnancy cohort (157 participants, 3 visits)
│   ├── group_a_probiotic/       # Probiotic trial (before/after)
│   └── group_b_contraception/   # Contraceptive methods comparison
│
├── R/
│   └── count_tables.R          # Data loading and phyloseq construction
│
├── .devcontainer/              # Codespace configuration
├── CLAUDE.md                   # AI assistant context
└── README.md
```

---

## Data Architecture

### Dataset Structure

Each of the three datasets contains six files with identical schemas:

| File | Description | Dimensions | Format |
|------|-------------|------------|--------|
| `01_metadata.csv` | Sample information (IDs, outcomes, covariates) | n_samples × variables | CSV |
| `02_microbiome_counts.csv` | Microbial species abundance (16S/metagenomics) | 220 species × n_samples | CSV, integer counts |
| `03_microbiome_taxonomy.csv` | Taxonomic classification (Domain→Species) | 220 × 7 ranks | CSV |
| `04_rnaseq_counts.csv.gz` | Human gene expression (RNA-seq) | 14,226 genes × n_samples | gzipped CSV, integer counts |
| `05_gene_annotation.csv` | Gene metadata (Ensembl ID, symbol, locus group) | 14,226 × annotations | CSV |
| `06_cytokines.csv` | Cytokine/chemokine measurements | n_samples × cytokines | CSV |

### Dataset Characteristics

| Dataset | Samples | Design | Outcome Variable |
|---------|---------|--------|------------------|
| `instructional` | 471 (157 participants × 3 visits) | Longitudinal pregnancy cohort | Term vs. preterm birth |
| `group_a` | varies | Probiotic intervention trial | Pre vs. post treatment |
| `group_b` | varies | Cross-sectional comparison | Contraceptive method |

### Data Orientation

Both count tables follow the convention: **features in rows, samples in columns** (transposed from the "tidy" format).

```
                Sample_001  Sample_002  Sample_003  ...  Sample_471
Species_001         142          0          89            203
Species_002           0         34           0              0
Species_003        1203        892        1547           2891
...
Species_220          45          0          12             67
```

### Sparsity Characteristics

| Table | Typical Zero Fraction | Explanation |
|-------|----------------------|-------------|
| Microbiome | ~89% | Many species absent from most samples; dominated by few taxa |
| RNA-seq | ~15–25% | Most genes detected in most samples; zeros often from low expression |

### Loading Data

The `R/count_tables.R` helper provides two interfaces:

#### 1. Raw Matrices and Tibbles

```r
source(here::here("R/count_tables.R"))

d <- read_counts("instructional")   # "instructional", "group_a", or "group_b"

# Returns a list:
d$metadata    # tibble: sample_id, participant_id, time_point, outcome, ...
d$species     # integer matrix: 220 rows × 471 columns, rownames = taxon IDs
d$genes       # integer matrix: 14,226 rows × 471 columns, rownames = Ensembl IDs
d$taxonomy    # tibble: taxon, Domain, Phylum, Class, Order, Family, Genus, species
d$gene_info   # tibble: ensembl_gene_id, symbol, locus_group, ...
d$cytokines   # tibble: sample_id + cytokine columns
```

#### 2. phyloseq Objects for microViz

```r
# Microbiome phyloseq object
ps_microbiome <- as_counts(
  dataset = "instructional",
  table = "microbiome",
  min_count = 0,      # minimum count threshold
  min_samples = 0,    # minimum samples a feature must appear in
  fix = TRUE          # run tax_fix() to fill blank labels
)

# RNA-seq phyloseq object (with recommended filtering)
ps_rnaseq <- as_counts(
  dataset = "instructional",
  table = "rnaseq",
  min_count = 10,     # require ≥10 counts
  min_samples = 100,  # in ≥100 samples
  fix = TRUE
)
```

#### phyloseq Object Structure

```r
ps_rnaseq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8934 taxa and 471 samples ]
# sample_data() Sample Data:       [ 471 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 8934 taxa by 2 taxonomic ranks ]
```

| Slot | Microbiome | RNA-seq |
|------|------------|---------|
| `otu_table()` | Species counts | Gene counts |
| `sample_data()` | Metadata | Metadata |
| `tax_table()` | Domain → Phylum → ... → species | locus_group → symbol |

---

## Analysis Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           ANALYSIS PIPELINE                                  │
└─────────────────────────────────────────────────────────────────────────────┘

     ┌──────────────┐         ┌──────────────┐         ┌──────────────┐
     │  Microbiome  │         │   RNA-seq    │         │  Cytokines   │
     │    Counts    │         │    Counts    │         │              │
     │  (220 × n)   │         │ (14,226 × n) │         │   (n × k)    │
     └──────┬───────┘         └──────┬───────┘         └──────┬───────┘
            │                        │                        │
            └───────────┬────────────┴────────────────────────┘
                        │
                        ▼
             ┌──────────────────┐
             │   read_counts()  │  ← Unified data loading
             └────────┬─────────┘
                      │
          ┌───────────┴───────────┐
          ▼                       ▼
 ┌─────────────────┐    ┌─────────────────┐
 │   as_counts()   │    │   as_counts()   │
 │table="microbiome│    │ table="rnaseq"  │
 └────────┬────────┘    └────────┬────────┘
          │                      │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────┐
          │  phyloseq object │  ← Unified data structure
          │                  │
          │  otu_table()     │  counts matrix
          │  sample_data()   │  metadata
          │  tax_table()     │  feature hierarchy
          └────────┬─────────┘
                   │
     ┌─────────────┼─────────────┬─────────────────┐
     ▼             ▼             ▼                 ▼
┌─────────┐  ┌──────────┐  ┌───────────┐    ┌────────────┐
│microViz │  │  DESeq2  │  │   fgsea   │    │ PERMANOVA  │
│         │  │          │  │           │    │            │
│• PCA    │  │• NB GLM  │  │• Ranked   │    │• Distance  │
│• PCoA   │  │• Wald    │  │  gene set │    │  partition │
│• NMDS   │  │  test    │  │  analysis │    │  by groups │
└─────────┘  └──────────┘  └───────────┘    └────────────┘
```

### Module Progression

| Module | Topic | Methods | Key Functions |
|--------|-------|---------|---------------|
| **02** | Count Table Anatomy | Sparsity analysis, library size, distributions | `summary()`, `hist()`, `rowSums()` |
| **03** | Ordination | Distance metrics, dimensionality reduction | `tax_transform()`, `dist_calc()`, `ord_calc()` |
| **04** | Differential Testing | Negative binomial GLM, shrinkage | `DESeqDataSetFromMatrix()`, `DESeq()`, `results()` |
| **05** | Enrichment | Ranked gene set enrichment | `fgsea()`, `msigdbr()` |

---

## Technical Methods

### Data Transformations

Count data requires transformation before certain analyses. The choice depends on the downstream method:

| Transform | Formula | Use Case | Properties |
|-----------|---------|----------|------------|
| **Raw counts** | None | DESeq2, edgeR | Required for NB models |
| **Relative abundance** | x_ij / Σ_i x_ij | Descriptive, stacked bars | Sums to 1 per sample |
| **CLR** (Centered log-ratio) | log(x_ij / g(x_j)) | PCA, Euclidean distance | Handles compositionality |
| **VST** (Variance stabilizing) | DESeq2::vst() | Heatmaps, clustering | Removes mean-variance relationship |
| **rlog** (Regularized log) | DESeq2::rlog() | Small sample sizes | Shrinks low-count variance |

#### CLR Transform Details

The centered log-ratio addresses compositionality by expressing each value relative to the geometric mean:

```
CLR(x_i) = log(x_i / g(x))

where g(x) = (∏ x_i)^(1/n) = exp(mean(log(x)))
```

**Zero handling**: CLR is undefined for zeros. Common approaches:
- Add pseudocount (e.g., +1 or +0.5)
- Bayesian replacement (e.g., `zCompositions::cmultRepl()`)
- microViz default: zeros replaced with small value before transform

```r
# microViz CLR transform
ps |> tax_transform("clr", rank = "species", zero_replace = "halfmin")
```

### Distance Metrics

Distance/dissimilarity metrics quantify sample differences for ordination and PERMANOVA:

| Metric | Formula | Range | Best For | Sensitive To |
|--------|---------|-------|----------|--------------|
| **Bray-Curtis** | Σ\|x_i - y_i\| / Σ(x_i + y_i) | [0, 1] | Abundance data, ecology | Abundant taxa |
| **Jaccard** | 1 - \|A ∩ B\| / \|A ∪ B\| | [0, 1] | Presence/absence | Rare taxa equally |
| **Euclidean** | √Σ(x_i - y_i)² | [0, ∞) | CLR-transformed data | Large differences |
| **Aitchison** | Euclidean on CLR | [0, ∞) | Compositional data | Log-ratios |
| **UniFrac** | Phylogenetic branch length | [0, 1] | 16S microbiome | Phylogenetic structure |

#### Choosing a Distance Metric

```
Decision tree:
├── Microbiome data?
│   ├── Care about phylogeny? → UniFrac (weighted or unweighted)
│   ├── Focus on abundant taxa? → Bray-Curtis
│   └── Focus on presence/absence? → Jaccard
├── RNA-seq data?
│   └── After CLR transform → Euclidean (= Aitchison distance)
└── General compositional data?
    └── CLR + Euclidean
```

```r
# microViz distance calculation
ps |>
  tax_transform("clr", rank = "species") |>
  dist_calc("euclidean")  # or "bray", "jaccard", "gunifrac", "wunifrac"
```

### Ordination Methods

Ordination reduces high-dimensional data to 2–3 axes for visualization:

| Method | Input | Preserves | Axes Meaning |
|--------|-------|-----------|--------------|
| **PCA** | Features × samples matrix | Global variance | Principal components (linear combinations) |
| **PCoA** (MDS) | Distance matrix | Pairwise distances | Coordinates in reduced space |
| **NMDS** | Distance matrix | Rank order of distances | Non-metric (relative positions) |
| **RDA/CCA** | Features + constraints | Variance explained by constraints | Constrained axes |

#### PCoA (Principal Coordinates Analysis)

1. Compute pairwise distance matrix D (n × n)
2. Double-center: B = -½ J D² J, where J = I - 11ᵀ/n
3. Eigendecompose B = UΛUᵀ
4. Coordinates: X = U√Λ

```r
# microViz ordination pipeline
ps |>
  tax_transform("clr", rank = "species") |>
  dist_calc("euclidean") |>
  ord_calc(method = "PCoA") |>  # or "PCA", "NMDS", "RDA"
  ord_plot(colour = "outcome", shape = "time_point")
```

#### PERMANOVA (Permutational MANOVA)

Tests whether group centroids differ in multivariate space:

```
H₀: Group centroids are equal
H₁: At least one group centroid differs

Test statistic: pseudo-F = (SS_between / df_between) / (SS_within / df_within)
P-value: proportion of permuted F ≥ observed F
```

```r
# PERMANOVA with vegan
library(vegan)
adonis2(distance_matrix ~ outcome + time_point, 
        data = metadata, 
        permutations = 999,
        by = "margin")  # Type III SS
```

**Assumptions**:
- Groups have similar multivariate dispersion (check with `betadisper()`)
- Observations are independent (violated by repeated measures)

### Differential Testing with DESeq2

DESeq2 models count data with a negative binomial distribution and performs shrinkage estimation of dispersions.

#### Statistical Model

```
K_ij ~ NegativeBinomial(μ_ij, α_i)

log₂(μ_ij) = β_i0 + β_i1 × condition_j + ...

where:
  K_ij = count for gene i in sample j
  μ_ij = expected count
  α_i  = gene-wise dispersion
  β    = log₂ fold changes
```

#### Normalization: Size Factors

DESeq2 uses the median-of-ratios method:

```
s_j = median_i (K_ij / (∏_j K_ij)^(1/n))
```

This accounts for library size and composition differences.

#### Workflow

```r
library(DESeq2)

# 1. Create DESeqDataSet (RAW COUNTS ONLY)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,      # genes × samples, integers
  colData = metadata,            # must match column order!
  design = ~ covariate + condition
)

# 2. Pre-filter low counts (optional, improves speed)
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]

# 3. Run the pipeline
dds <- DESeq(dds)  # estimate size factors → estimate dispersions → fit NB GLM

# 4. Extract results
res <- results(dds, 
               contrast = c("condition", "treatment", "control"),
               alpha = 0.05,           # FDR threshold
               lfcThreshold = 0)       # test against this LFC

# 5. Shrink LFC estimates (for visualization/ranking)
res_shrunk <- lfcShrink(dds, 
                        coef = "condition_treatment_vs_control",
                        type = "apeglm")  # recommended shrinkage estimator
```

#### Key Output Columns

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | Effect size (log₂ scale) |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald statistic (LFC / lfcSE) |
| `pvalue` | Raw p-value from Wald test |
| `padj` | Benjamini-Hochberg adjusted p-value |

#### Multiple Testing Correction

DESeq2 uses Benjamini-Hochberg FDR correction:

```
padj_i = p_i × n / rank(p_i)

Interpretation: Of genes with padj < 0.05, expect ≤5% false discoveries
```

### Enrichment Analysis

#### Gene Set Enrichment Analysis (GSEA) with fgsea

Tests whether genes in predefined sets are enriched at the top or bottom of a ranked list.

```r
library(fgsea)
library(msigdbr)

# 1. Create ranked gene list (use DESeq2 stat column)
ranks <- setNames(res$stat, rownames(res))
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

# 2. Get gene sets (e.g., Hallmark collection)
gene_sets <- msigdbr(species = "Homo sapiens", category = "H") |>
  split(x = .$ensembl_gene, f = .$gs_name)

# 3. Run fgsea
fgsea_res <- fgsea(
  pathways = gene_sets,
  stats = ranks,
  minSize = 15,    # minimum genes in set
  maxSize = 500,   # maximum genes in set
  nPermSimple = 10000
)
```

#### Enrichment Score Calculation

1. Walk down ranked gene list
2. When gene is IN set: ES increases by √((N - k) / k)
3. When gene is NOT in set: ES decreases by √(k / (N - k))
4. ES = maximum deviation from zero
5. NES = ES normalized by set size

```
NES > 0: genes enriched at TOP of list (upregulated)
NES < 0: genes enriched at BOTTOM of list (downregulated)
```

#### Gene Set Collections (MSigDB)

| Category | Code | Description | Sets |
|----------|------|-------------|------|
| Hallmark | H | Curated, coherent biological states | 50 |
| Positional | C1 | Chromosome cytogenetic bands | 299 |
| Curated | C2 | Canonical pathways (KEGG, Reactome, etc.) | 6,456 |
| Regulatory | C3 | microRNA and TF targets | 3,726 |
| Computational | C4 | Cancer modules, neighborhoods | 858 |
| Ontology | C5 | GO terms (BP, CC, MF) | 15,719 |
| Oncogenic | C6 | Cancer gene signatures | 189 |
| Immunologic | C7 | Immune cell signatures | 5,219 |
| Cell type | C8 | Cell type markers | 830 |

---

## The Unified Pipeline Concept

The key insight of this workshop: **microViz was built for microbiome data, but it works identically on gene expression data** when we reframe the annotation:

| Microbiome | RNA-seq | Role |
|------------|---------|------|
| Species ID | Ensembl Gene ID | Feature identifier |
| Taxonomy ranks | `locus_group` → `symbol` | Feature hierarchy |
| Domain → Species | Protein-coding → TP53 | Coarse → fine labels |

The `as_counts()` function constructs a phyloseq object from either table, mapping gene annotations to the taxonomy slot:

```r
# From count_tables.R - how gene annotations become "taxonomy"
ranks <- rank_matrix(info, "ensembl_gene_id", c("locus_group", "symbol"))

#                      locus_group       symbol
# ENSG00000000003      protein_coding    TSPAN6
# ENSG00000000005      protein_coding    TNMD
# ENSG00000000419      protein_coding    DPM1
```

This allows the same microViz pipeline to work unchanged:

```r
# Identical code structure, different data
as_counts(table = "microbiome") |>
  tax_transform("clr", rank = "species") |>
  dist_calc("euclidean") |>
  ord_calc("PCoA") |>
  ord_plot(colour = "outcome")

as_counts(table = "rnaseq") |>
  tax_transform("clr", rank = "symbol") |>
  dist_calc("euclidean") |>
  ord_calc("PCoA") |>
  ord_plot(colour = "outcome")
```

---

## Working with the AI Assistant

Pi is your AI coding assistant, billed through your OpenRouter account.

### Model Selection

| Model | Best For | Speed | Cost |
|-------|----------|-------|------|
| `z-ai/glm-5.3-flash` | Routine coding tasks | Fast | Low |
| `z-ai/glm-5.3` | Complex reasoning | Slower | Higher |

Switch models with `/model` or `Ctrl+L` in Pi.

### Usage Tips

- Pi reads `CLAUDE.md` automatically to understand your context
- Ask it to explain code, not biology
- Request verification steps: "How can I check this is correct?"
- If stuck, switch to the stronger model

📊 Track spending at [OpenRouter Activity](https://openrouter.ai/activity)

---

## Group Project

On Wednesday and Thursday, participants work in groups on one of:
- **Study A**: Probiotic trial (before/after intervention)
- **Study B**: Contraceptive methods comparison

### Workflow

1. **QC Report** (`01-qc.qmd`): Data quality assessment
2. **Analysis Documents** (`02-*.qmd`): One per research question
3. **Documentation**: Each document ends with a "Choices Made" section

### Key Principle: Avoiding Circular Analysis

> **Do not define a group from one count table and then test that same table.**

**Correct workflow**:
```
Microbiome data → Define "Lactobacillus-dominant" (>50% relative abundance)
                         ↓
Gene expression data → Test: Which genes differ between community types?
```

**Why this matters**: If you define groups based on microbiome composition, those same features will trivially differ between groups. Cross-omic testing asks a meaningful biological question: "Does the host transcriptome differ by microbial community state?"

📖 [Group Project prompts and guidance](group_project/README.md)

---

## Important Considerations

### Statistical Best Practices

| ✅ Do | ❌ Don't |
|-------|---------|
| Feed **raw integer counts** to DESeq2 | Use CPM, TPM, RPKM, or log-transformed data |
| Verify `colnames(counts) == rownames(metadata)` | Assume alignment is automatic |
| Pre-filter low-count features before testing | Test all 14,226 genes (wastes power on noise) |
| Report effect sizes (LFC) alongside p-values | Rely on p-values alone |
| Use `lfcShrink()` for visualization/ranking | Use raw LFC for noisy genes |
| Check dispersion estimates plot | Ignore QC diagnostics |

### Compositional Data Awareness

Both count tables are **compositional**: values represent relative abundances constrained by sequencing depth (library size).

#### The Compositionality Problem

If Species A truly doubles while everything else stays constant:
- True abundance: A goes from 1000 → 2000 copies
- Observed: A goes from 10% → 18% (not 20%, because total increased)
- Other species: appear to decrease even though they didn't change

#### Implications

| Statement | Validity |
|-----------|----------|
| "A has higher relative abundance in group 1" | ✓ Valid (describes data) |
| "There is more A in group 1" | ✗ Invalid (claims absolute abundance) |
| "B decreased in group 1" | ✗ Invalid (may be compositional artifact) |

#### Mitigation Strategies

1. **CLR transform**: Expresses changes as log-ratios, reducing compositional bias
2. **Spike-ins**: Add known quantities before sequencing (rarely available)
3. **Careful language**: "enriched" or "elevated relative abundance" rather than "increased"

### Sample Filtering Recommendations

| Data Type | Suggested Filter | Rationale |
|-----------|------------------|-----------|
| Microbiome | Total reads > 1000 | Below this, composition estimates are unstable |
| RNA-seq | Total reads > 1M | Standard QC threshold for bulk RNA-seq |
| Per-feature | ≥10 counts in ≥n samples | Remove noise, improve multiple testing power |

### Common Pitfalls

| Pitfall | Problem | Solution |
|---------|---------|----------|
| Treating zeros as "not present" | May be present below detection | Use zero-aware methods or presence threshold |
| Ignoring repeated measures | Inflated significance | Use mixed models or aggregate |
| Cherry-picking features | Confirmation bias | Pre-register analysis plan |
| Batch effects | Confound biological signal | Include batch in model or use ComBat |

### Data Provenance

> ⚠️ **Important**: These datasets are **constructed for teaching**. The counts are real sequencing data from real samples, but the study designs and outcomes were created for this workshop. Analyze them properly, but do not cite the results in publications.

---

## References

### Methods Papers

- **DESeq2**: Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. [DOI:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

- **fgsea**: Korotkevich G, Sukhov V, Sergushichev A (2019). Fast gene set enrichment analysis. *bioRxiv*. [DOI:10.1101/060012](https://doi.org/10.1101/060012)

- **microViz**: Barnett DJ, Arts IC, Penders J (2021). microViz: an R package for microbiome data visualization and statistics. *JOSS* 6(63):3201. [DOI:10.21105/joss.03201](https://doi.org/10.21105/joss.03201)

- **phyloseq**: McMurdie PJ, Holmes S (2013). phyloseq: An R package for reproducible interactive analysis of microbiome census data. *PLoS ONE* 8(4):e61217. [DOI:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)

- **Compositional data**: Gloor GB et al. (2017). Microbiome datasets are compositional: And this is not optional. *Frontiers in Microbiology* 8:2224. [DOI:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)

### Documentation

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [microViz documentation](https://david-barnett.github.io/microViz/)
- [fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [MSigDB gene sets](https://www.gsea-msigdb.org/gsea/msigdb/)

---

## License & Attribution

This workshop environment was developed for the Durban Data Science for Biology workshop series.

---

<div align="center">

**Questions?** Ask Pi, your AI coding assistant.

**Issues?** Open a GitHub issue on this repository.

</div>
