---
title: "The Global Threat of ST147 Klebsiella pneumoniae
"
layout: post
date: 2024-07-02
tag:
- Klebsiella pneumoniae
- ST147
- MDR (Multidrug-resistant)
- Nosocomial infections
- Antimicrobial resistance (AMR)
- Time-scaled Phylogenetics 
- R programming



projects: true
description: "The Global Threat of ST147 Klebsiella pneumoniae"
---

**Bangalore, India 2023**

---
## Introduction
**ST147 *Klebsiella pneumoniae* has emerged as a significant global health concern due to its multidrug resistance and rapid spread. While extensively studied in European and North American settings, its prevalence and characteristics in India and Southeast Asia remain under-explored. This study aimed to comprehensively characterize ST147 isolates from India and compare them to a global dataset to understand their unique features and contribution to the global burden of MDR *K. pneumoniae*.**

## Background
- ST147's Global Rise: ST147 *K. pneumoniae* has been identified as a major contributor to the global burden of antimicrobial resistance, particularly in healthcare settings. <br>
- India's Unique Challenges: India faces a unique challenge due to the high prevalence of ST147 strains carrying various carbapenemase genes.<br>
- Knowledge Gaps: Despite previous studies, a comprehensive understanding of ST147's epidemiology, transmission pathways, and molecular mechanisms in India remains elusive.<br>

## My Role in the Study
As a member of the Bioinformatics team, I played a crucial role in:

- Literature Review: Conducting a thorough literature review to identify knowledge gaps and inform our research design.<br>
- Methodology Development: Contributing to the development of the research methodology, including data collection and analysis techniques.<br>
- Data Analysis: Actively participating in the analysis of genomic data to characterize ST147 isolates and identify their unique features.
- Manuscript Preparation: Contributing to the writing and editing of the research manuscript<br>
## Methods
- **Genome Sequencing:** Whole-genome sequencing of ST147 isolates using the Illumina MiSeq platform.<br>
- **Genomic Characterization:** De novo assembly and genomic characterization using Kleborate.<br>
- **Antimicrobial Resistance Analysis:** In silico prediction of antimicrobial resistance using Kleborate and AMRfinderPlus.<br>
- **Plasmid Analysis:** Identification of plasmid replicon types and their association with resistance genes.<br>
- **Global Comparison:** Comparison of Indian ST147 isolates with a global dataset to identify unique features.<br>

## Results
**Global Phylogeny and Time-Scaled Analysis**<br>
I played a pivotal role in constructing the time-scaled phylogeny using the GHRU-SNP phylogeny pipeline and removing recombinations using the Gubbins software. By carefully analyzing the genomic data, we were able to identify multiple introduction events of ST147 into India, suggesting ongoing clonal expansion. Additionally, we employed the BactDating R package to estimate the time of the most recent common ancestor (MRCA), finding it to be around 1987.<br>
**Plasmid Replicon Types and Resistance Gene Associations**
I also contributed to the analysis of plasmid replicon types and their association with resistance genes. Using tools like Plasmidfinder and Mobsuite, We were able to identify the nature of contigs containing AMR genes and determine whether they were carried on plasmids or chromosomes. This analysis provided valuable insights into the molecular mechanisms underlying the observed multidrug resistance phenotype facilitated by hybrid plasmids carrying carbapenemase genes.

## Here are some R scripts that I have written to contribute to this study:
### Phylogenetic Tree Plot with Metadata Bars

In this section, we load the phylogenetic tree data and metadata, create a plot using the `ggtree` package, and overlay bars for each gene from the metadata.

### Required Libraries
```r
# Load required libraries
library(ggtree)
library(ape)
library(tidyr)
```
Libraries: We load ggtree for phylogenetic tree visualization, ape for handling tree objects, and tidyr for data manipulation.

### Reading Metadata and Tree File
```r
# Read the metadata from the CSV file
metadata <- read.csv("/data/internship_data/srikanth_kpn/global_st147_fastqs/metadata/lj.csv", stringsAsFactors = FALSE)

# Read the tree file (assuming it's in Newick format)
tree <- read.tree("/data/internship_data/srikanth_kpn/global_st147_fastqs/new_fastqs/snp_phylogeny_output/2023/gubbins_new/mid_point_rooted_ape.nwk")
```
- Metadata: We read the metadata file containing information about the phylogenetic tree tips.<br>
- Tree File: The phylogenetic tree is read in Newick format using read.tree.<br>
### Preparing the Tree Data Frame
```r
# Initialize vectors to store tip labels and branch lengths
tip_labels <- character(length = 2 * tree$Nnode + 1)
branch_lengths <- numeric(length = length(tip_labels))

# Assign tip labels and branch lengths
tip_labels[tree$edge[, 1]] <- tree$tip.label
branch_lengths[tree$edge[, 1]] <- tree$edge.length

# Create the tree_df data frame
tree_df <- data.frame(tip.label = tip_labels[tip_labels != ""], branch = branch_lengths[tip_labels != ""])

# Save the data frame to a CSV file
write.csv(tree_df, file = "/data/internship_data/srikanth_kpn/global_st147_fastqs/metadata/tree_df.csv")
```
- Tree Data: We extract the tip labels and branch lengths from the tree and store them in a data frame tree_df.<br>
- Save Data: The data frame is saved as a CSV file for future use.<br>
### Merging Metadata with Tree Data
```r
# Calculate y-coordinates for the bars
y_coords <- seq(0, -1, length.out = nrow(metadata))

# Merge filtered metadata with the tree data based on sample IDs
tree_with_metadata <- left_join(tree_df, metadata, by="tip.label")
```
- Y-Coordinates: We generate y_coords for plotting bars based on the number of metadata rows.<br>
- Merge Data: The tree_df and metadata are merged on the tip.label column.<br>
### Plotting the Tree with Metadata Bars
```r
# Make the original tree plot
p <- ggtree(tree_with_metadata)

# Create bars for each gene
for (gene in colnames(metadata)[-1]) {  # Exclude the 'ids' column
  tree_with_metadata <- mutate(tree_with_metadata, !!paste0(gene, "_y") := y_coords)
  
  p <- p + geom_segment(data = tree_with_metadata,
                        aes(x = branch, xend = branch, y = !!sym(paste0(gene, "_y")),
                            yend = !!sym(paste0(gene, "_y")), color = !!sym(gene)),
                        size = 3) +
    scale_color_manual(name = gene, values = c("yes" = "red", "no" = "white"), guide = "none")
}

# Show the plot
print(p)
```
- Plot: A base tree plot is created using ggtree.<br>
- Bars: For each gene in the metadata, colored bars (red for 'yes', white for 'no') are added to the plot, representing the presence or absence of the gene.<br

## Script 2: BactDating Phylogenetic Analysis
## BactDating Phylogenetic Analysis

This section demonstrates how to run **BactDating** to estimate the divergence times in a phylogenetic tree, including loading the tree, merging metadata, and visualizing the results.

### Required Libraries and Reading Data
```r
# Load libraries
library(ape)
library(BactDating)

# Read the tree file and metadata
t <- read.tree(file="~/pw_fastqs/snp_phylogeny_output/gubbin_out/aligned_pseudogenome.node_labelled.final_tree.tre")
metadata <- read.csv(file="~/Kpn_data/bactdating/dates.csv")
```
- Libraries: ape is used for handling phylogenetic trees, and BactDating is used for dating trees.<br>
- Tree and Metadata: The phylogenetic tree and metadata with date information are loaded.<br>
### Merging Metadata with Tree Labels
```r
# Extract tree tip labels and merge with metadata
tree_labels <- as.data.frame(t$tip.label)
merged <- merge(tree_labels, metadata, by.x='t$tip.label', by.y='name', all.x = TRUE, all.y = FALSE)

# Reorder to match the tree tip order
merged2 <- merged[match(t$tip.label, merged$`t$tip.label`),]
```
- Merge: The tip labels from the tree are merged with the metadata based on sample names.
- Reordering: The merged data is reordered to match the order of the tree tip labels.
### Initializing and Running BactDating
```r
# Initialize the rooted tree
rooted <- initRoot(t, merged2$date)

# Root-to-tip regression for time estimates
r <- roottotip(rooted, merged2$date)

# Run BactDating with a relaxed gamma model
res <- bactdate(unroot(t), merged2$date, nbIts=10000, initSigma=0.000005, 
                updateSigma=TRUE, updateRoot=TRUE, updateAlpha=TRUE, 
                updateMu=TRUE, model="relaxedgamma")
```
- Rooting the Tree: We initialize the rooted tree with the date metadata.<br>
- Root-to-Tip Regression: roottotip() performs a regression to estimate divergence times.<br>
- BactDating: bactdate() is used to estimate dates using the relaxedgamma model, running for 10,000 iterations.<br>
### Visualizing the BactDating Results
```r
# Plot the results with confidence intervals
plot(res, 'treeCI', show.tip.label = TRUE)

# Additional BactDating run on previously rooted tree
res0 <- bactdate(rooted, merged2$date)

# Run BactDating again on unrooted tree from root-to-tip regression
res2 <- bactdate(unroot(r), merged2$date, nbIts=10000, initSigma=0.000005, 
                 updateSigma=TRUE, updateRoot=TRUE, updateAlpha=TRUE, 
                 updateMu=TRUE, model="relaxedgamma")

# Plot second set of results
plot(res0, 'treeCI', show.tip.label = TRUE)
```
- Plot: The BactDating results are visualized with confidence intervals around the divergence time estimates.
- Additional Runs: Additional BactDating analyses are performed and plotted, including both rooted and unrooted tree versions.