---
editor_options: 
  markdown: 
    wrap: 72
---

# Intro to Ecological Genomics 2026

### Author: Stephanie Iwinski

### 9/8/26

### Introduction to Notetaking

#### notess

### 9/17/26

### Data intro

### accesing our data in a z file

#### commands I used on the VACC

`cd /gpfs1/cl/biol3990`

`pwd`

`ll`

`ls`

`cd Transcriptomics/`

`ll`

`cd Cleandata/`

`ll`

`zcat AA_F0_Rep3_2_clean.fq.gz  | head -n 8`

`zcat AA_F0_Rep3_2_clean.fq.gz  | head -n 4`

### word count `zcat AA_F0_Rep3_2_clean.fq.gz  | wc -l`

### wc: 93638040

## 9/15/2026 - Setting up lab notebook and learning markdown

- setting up transcriptomics notebook

- learn how to take notes in markdown

- push notes to github

**Working Directory:**

`/gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/transcriptomics`

**Input files:**

`none`

**Output files:**

`/gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and Dependencies:**

- `R version 4.5.1`

- `R-studio`

**Scripts:**

`none`

**Code:**

``` r
print ("Hello World") 
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

**Notes/Observations:**

- oh cool graph! it makes sense! and some interpretation of your data

------------------------------------------------------------------------

## 9/15/2026 - Setting up lab notebook and learning markdown

- setting up transcriptomics notebook

- learn how to take notes in markdown

- push notes to github

**Working Directory:**

`/gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/transcriptomics`

**Input files:**

`none`

**Output files:**

`/gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and Dependencies:**

- `R version 4.5.1`

- `R-studio`

**Scripts:**

`none`

**Code:**

``` r
print ("Hello World") 
```

**Table:**

## 11/17/2026 - Group Project Notes

Copy 16S data and scripts into final_project folder (using terminal):

```{bash}
cp -r /gpfs1/cl/biol3990/PopulationGenomics/16S /gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/final_project
```

Created a new script titled `PCA_zebrafish.R`

```{r}
# DATA DIMENSIONS
cat("OTU table:", nrow(otu), "taxa x", ncol(otu), "samples\n")
# OTU table: 36756 taxa x 473 samples
cat("Taxonomy:", nrow(tax), "taxa\n")
# Taxonomy: 73171 taxa
cat("Metadata:", nrow(metadata), "samples\n")
# Metadata: 477 samples
cat("Tree:", length(phy_tree$tip.label), "tips\n")
# Tree: 36756 tips
# OTU taxa found in taxonomy: 36756 of 36756 
# OTU samples found in metadata: 473 of 473 
# Number of taxa: 36756 
# Number of samples: 473 


# "1] Min. number of reads = 33"

# "2] Max. number of reads = 76697"

# "3] Total number of reads = 17882766"

# "4] Average number of reads = 37807.1162790698"

# "5] Median number of reads = 40897"

# "7] Sparsity = 0.993582558151039"

# "6] Any OTU sum to 1 or less? NO"

# "8] Number of singletons = 0"

# "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0"

# "10] Number of sample variables are: 11"

# > ntaxa(ps)
# [1] 36756
# > nsamples(ps)
# [1] 473
```

```{r}
# Number of samples: 435 
# Number of features: 36728 
# Minimum reads/sample: 33 
# Median reads/sample: 39395 
# Maximum reads/sample: 76697 
```

Rarification, removing the samples with low amounts of reads.

```         
> print(sample_retention)
  depth samples_retained samples_total proportion_retained percent_retained samples_lost
1  1000              432           435          0.99310345        99.310345            3
2  2500              431           435          0.99080460        99.080460            4
3  5000              429           435          0.98620690        98.620690            6
4 10000              409           435          0.94022989        94.022989           26
5 15000              385           435          0.88505747        88.505747           50
6 20000              353           435          0.81149425        81.149425           82
7 30000              284           435          0.65287356        65.287356          151
8 60000               25           435          0.05747126         5.747126          410
```
