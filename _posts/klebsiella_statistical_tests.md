---
title: "K and O loci
"
layout: post
date: 2024-07-02
tag:
- machine learning
- R
- ESKAPE pathogens
- public health
- healthcare
- epidemiology
- global health
- infectious diseases
- bacteria
- antibiotic resistance



projects: true
description: "Predicting Resistance Trends of ESKAPE Pathogens using R"
---

**Bangalore, India 2023**

---
## Introduction
**This project contains R scripts associated with a study ocused on understanding the geographical distribution and diversity of *Klebsiella pneumoniae* serotypes in India. Using genomic data from *Klebsiella pneumoniae* isolates, we employed Fisher’s exact test to examine the relationship between specific serotypes and geographical regions. I employed Simpson’s diversity index to assess serotype diversity within different regions. By combining these statistical analyses with the geographical data, I was able to identify patterns in serotype distribution and uncover potential regional variations in the pathogen’s virulence and antibiotic resistance.**

This work is published in the journal Microbial Genomics under the DOI: [https://doi.org/10.1099/mgen.0.001271]

## Step-by-step Process to use Fisher's exact test on my dataset

### Step 1: Load Required Libraries

We load the necessary R libraries to handle data manipulation and perform Fisher's exact tests.

```r
# Load necessary libraries
library(dplyr)
library(tidyverse)
library(magrittr)
```
### Step 2: Read the Data
The CSV file contains information on the Klebsiella pneumoniae isolates, including virulence factors and K/O loci data.
```r
# Load dataset from CSV
mr_data_2 <- read.csv("~/Downloads/MR_data_K & O Loci _1072.csv")
```
### Step 3: Perform Fisher's Exact Test for Virulence Factors
Here, we perform Fisher's exact test to determine if there is a statistically significant association between the presence of specific virulence factors (Yersiniabactin, Colibactin, Aerobactin, Salmochelin, RmpADC, rmpA2) and invasive disease.
```r
# Fisher's exact test for Yersiniabactin
contingency_table_invas_yer <- table(mr_data_2$Invasive, mr_data_2$Yersiniabactin)
result_invas_yer <- fisher.test(contingency_table_invas_yer)
print(result_invas_yer)

# Fisher's exact test for Colibactin
contingency_table_invas_coli <- table(mr_data_2$Invasive, mr_data_2$Colibactin)
result_invas_coli <- fisher.test(contingency_table_invas_coli)
print(result_invas_coli)

# Fisher's exact test for Aerobactin
contingency_table_invas_aero <- table(mr_data_2$Invasive, mr_data_2$Aerobactin)
result_invas_aero <- fisher.test(contingency_table_invas_aero)
print(result_invas_aero)

# Fisher's exact test for Salmochelin
contingency_table_invas_salmo <- table(mr_data_2$Invasive, mr_data_2$Salmochelin)
result_invas_salmo <- fisher.test(contingency_table_invas_salmo)
print(result_invas_salmo)

# Fisher's exact test for RmpADC
contingency_table_invas_rmpADC <- table(mr_data_2$Invasive, mr_data_2$RmpADC)
result_invas_rmpADC <- fisher.test(contingency_table_invas_rmpADC)
print(result_invas_rmpADC)

# Fisher's exact test for rmpA2
contingency_table_invas_rmpa2 <- table(mr_data_2$Invasive, mr_data_2$rmpA2)
result_invas_rmpa2 <- fisher.test(contingency_table_invas_rmpa2)
print(result_invas_rmpa2)
```
### Step 4: Fisher's Test for K/O Loci vs Invasive Disease
Here, we test whether there is an association between specific K and O loci (e.g., KL64, KL51, O1/O2v1) and invasive disease.
```r
# Define K and O loci of interest
k_types_of_interest <- c("KL64", "KL51","KL2","KL10")
o_types_of_interest <- c("O1/O2v1", "O1/O2v2", "OL101")

# Fisher's exact test for K loci
for (k_locus in k_types_of_interest) {
  contingency_table_k <- table(mr_data_2$Invasive, mr_data_2$K_locus == k_locus)
  fisher_result_k <- fisher.test(contingency_table_k)
  cat("Fisher's exact test for K-type", k_locus, ":\n")
  cat("P-value:", fisher_result_k$p.value, "\n\n")
}

# Fisher's exact test for O loci
for (o_locus in o_types_of_interest) {
  contingency_table_o <- table(mr_data_2$Invasive, mr_data_2$O_locus == o_locus)
  fisher_result_o <- fisher.test(contingency_table_o)
  cat("Fisher's exact test for O-type", o_locus, ":\n")
  cat("P-value:", fisher_result_o$p.value, "\n\n")
}
```
### Step 5: Fisher's Test for K Loci vs Sequence Type (ST)
We perform Fisher's test to determine the relationship between Klebsiella pneumoniae sequence types (ST) and specific K loci.
```r
# Fisher's exact test for K Loci vs ST
fisher_test_results_st <- data.frame()

for (st in unique(mr_data_2$ST)) {
  st_subset <- mr_data_2 %>% filter(ST == st)
  unique_combinations <- st_subset %>% distinct(ST, K_locus)
  
  for (i in 1:nrow(unique_combinations)) {
    k_locus <- unique_combinations[i, "K_locus"]
    
    subset_df <- st_subset %>% filter(ST == st, K_locus == k_locus)
    
    if (nrow(subset_df) >= 2) {
      contingency_table <- matrix(c(nrow(subset_df), 0, 0, 0), nrow = 2)
      contingency_table[1, "Present"] <- sum(subset_df$K_locus == k_locus)
      contingency_table[2, "Absent"] <- sum(subset_df$K_locus != k_locus)
      
      result <- fisher.test(contingency_table)
      fisher_test_results_st <- bind_rows(fisher_test_results_st, data.frame(ST = st, K_locus = k_locus, p_value = result$p.value))
    }
  }
}

print(fisher_test_results_st)
```
### Step 6: Fisher's Test for K/O Loci vs Carbapenem Genes
Here, we perform Fisher's exact test to determine the association between K and O loci and the presence of carbapenemase genes.
```r
# Fisher's exact test for K Loci vs Carbapenem Genes
results_KL <- data.frame(K_locus = character(), P_value = numeric())
unique_k_locus <- unique(mr_data_2$K_locus)

for (k_locus in unique_k_locus) {
  contingency_table_k <- table(mr_data_2$Bla_Carb_acquired, mr_data_2$K_locus == k_locus)
  fisher_result_k <- fisher.test(contingency_table_k, simulate.p.value = TRUE)
  results_KL <- rbind(results_KL, data.frame(K_locus = k_locus, P_value = fisher_result_k$p.value))
}

print(results_KL)
# Fisher's exact test for O Loci vs Carbapenem Genes
results_O <- data.frame(O_locus = character(), P_value = numeric())
unique_o_locus <- unique(mr_data_2$O_locus)

for (o_locus in unique_o_locus) {
  contingency_table_o <- table(mr_data_2$Bla_Carb_acquired, mr_data_2$O_locus == o_locus)
  fisher_result_o <- fisher.test(contingency_table_o, simulate.p.value = TRUE)
  results_O <- rbind(results_O, data.frame(O_locus = o_locus, P_value = fisher_result_o$p.value))
}

print(results_O)
```

### Conclusion
This analysis provides insights into the relationship between virulence factors, K/O loci, and key characteristics in Klebsiella pneumoniae, leveraging Fisher's exact test to assess statistical significance.

## Step-by-step process employed to use Simpson's diversity index on my data
This analysis calculates **Simpson's Diversity Index** for the K-locus and O-locus across different regions and age groups. The script performs calculations for both the full dataset and deduplicated data, and it confirms the results using the `vegan` package.

### Libraries and Data
```r
# Load necessary libraries
library(dplyr)
library(vegan)
library(magrittr)
library(ggplot2)
```
### Full Dataset: Simpson's Diversity Index for K-locus by Region
- Data: Full dataset of 1072 entries <br>
- Calculation: Group the data by region and calculate Simpson's Diversity Index using the K-locus data.<br>
```r
# Load the full dataset
mr_all_data <- read.csv("~/Downloads/MR_data_K & O Loci _1072.csv")

# Group by region and calculate Simpson's Diversity Index for K-locus
diversity_by_region_K_all_data <- mr_all_data %>%
  group_by(Regions) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(K_locus) / sum(table(K_locus)))^2))

# Print the results
print(diversity_by_region_K_all_data)
```
### Full Dataset: Simpson's Diversity Index for K-locus by Age Group
- Data: Full dataset of 1072 entries <br>
- Calculation: Group the data by age group and calculate Simpson's Diversity Index for the K-locus data. <br>
```r
# Group by age and calculate Simpson's Diversity Index for K-locus
diversity_by_age_K_all_data <- mr_all_data %>%
  group_by(Age_group) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(K_locus) / sum(table(K_locus)))^2))

# Print the results
print(diversity_by_age_K_all_data)
```
### Full Dataset: Simpson's Diversity Index for O-locus by Region
- Data: Full dataset of 1072 entries <br>
- Calculation: Group the data by region and calculate Simpson's Diversity Index using O-locus data. <br>
```r
# Group by region and calculate Simpson's Diversity Index for O-locus
diversity_by_region_O_all_data <- mr_all_data %>%
  group_by(Regions) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(O_locus) / sum(table(O_locus)))^2))

# Print the results
print(diversity_by_region_O_all_data)
```
### Deduplicated Dataset: Simpson's Diversity Index for K-locus by Region
- Data: Deduplicated dataset <br>
- Calculation: Group the deduplicated data by region and calculate Simpson's Diversity Index for K-locus. <br>
```r
# Load deduplicated data
derep_data <- read.csv("~/Downloads/Derplicated_MR.csv")

# Group by region and calculate Simpson's Diversity Index for K-locus
diversity_by_region_K_derep_data <- derep_data %>%
  group_by(Regions) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(K_locus) / sum(table(K_locus)))^2))

# Print the results
print(diversity_by_region_K_derep_data)
```
### Deduplicated Dataset: Simpson's Diversity Index for K-locus by Age Group
- Data: Deduplicated dataset <br>
- Calculation: Group the deduplicated data by age group and calculate Simpson's Diversity Index for K-locus. <br>
```r
# Group by age and calculate Simpson's Diversity Index for K-locus
diversity_by_age_K_derep_data <- derep_data %>%
  group_by(Age_group) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(K_locus) / sum(table(K_locus)))^2))

# Print the results
print(diversity_by_age_K_derep_data)
```
### Deduplicated Dataset: Simpson's Diversity Index for O-locus by Region
- Data: Deduplicated dataset <br>
- Calculation: Group the deduplicated data by region and calculate Simpson's Diversity Index for O-locus. <br>
```r
# Group by region and calculate Simpson's Diversity Index for O-locus
diversity_by_region_O_derep <- derep_data %>%
  group_by(Regions) %>%
  summarize(Simpson_Diversity_Index = 1 - sum((table(O_locus) / sum(table(O_locus)))^2))

# Print the results
print(diversity_by_region_O_derep)
```
### Reconfirmation Using the Vegan Package
- Purpose: Reconfirm the Simpson’s Diversity Index using the vegan package.
- Steps:
  - Create a table of counts of K-locus per region.
  - Use the diversity() function from the vegan package.
```r
# Create a table of K-locus counts by region
k_locus_table_derep <- table(derep_data$Regions, derep_data$K_locus)

# Calculate Simpson's Diversity Index for each region using the vegan package
simpson_diversity_derep <- apply(k_locus_table_derep, 1, function(x) diversity(x, index = "simpson"))

# Print the results
print(simpson_diversity_derep)
```
### Visualizing Simpson's Diversity Index (Optional)
You can use the ggplot2 library to visualize the diversity indices by region. <br>
```r
# Visualize Simpson's Diversity Index by region for the deduplicated K-locus data
ggplot(diversity_by_region_K_derep_data, aes(x = Regions, y = Simpson_Diversity_Index)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Regions", y = "Simpson's Diversity Index") +
  ggtitle("Simpson's Diversity Index by Region for K-locus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
### Conclusion
In this analysis, we calculated the Simpson’s Diversity Index across different groupings (region and age group) for both the full dataset and the deduplicated data. We confirmed the results using the vegan package, and visualized them using ggplot2.