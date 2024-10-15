---
title: "Predicting Resistance Trends of ESKAPE Pathogens using R
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

**This project is aimed at tackling the escalating threat of antimicrobial resistance (AMR) posed by ESKAPE pathogens. We analyze a the Pfizer ATLAS dataset (8 lakh+ isolates) from the Vivli AMR Register competition (2023) to understand their prevalence, geographical distribution, and resistance trends. Additionally, we leverage machine learning to predict future resistance patterns and the analysis is conducted on a High-Performance Computing (HPC) environment using R Studio to handle the massive dataset size effectively.**


## Sample figures generated from the project

## 1. Data Cleaning and Preparation
This code utilizes the dplyr and tidyr libraries to filter and manipulate the raw data based on specific criteria. You can adjust the filter conditions within the filter() function to select desired species, countries, etc. The filtered data is then saved as a new CSV file.
```R
# Sorting the raw data based on requirement
library(dplyr)
library(magrittr)
library(tidyr)

# Read the data from the CSV file
raw_data <- read.csv("/path/to/your/csvfile.csv")

# Filter the data for a specific species and country
# Replace "Species" and "Country" with the appropriate column names
# Replace "Species_Name" and "Country_Name" with the desired species and country names
filtered_data <- raw_data %>%
  filter(Species == "Species_Name" & Country == "Country_Name")

# Save the filtered data to a new CSV file
# Replace "/path/to/your/output/file.csv" with the desired output path and file name
write.csv(filtered_data, '/path/to/your/output/file.csv')
```
- Removing Empty Columns:   
This snippet employs the dplyr library to identify and eliminate empty columns from the dataset. The which() function locates empty columns, and these are subsequently removed using [, -empty_cols]. Finally, the modified data is saved as a new CSV file.
```R
# Remove any empty columns
# Load required libraries
library(dplyr)

# Read the data from the CSV file
# Replace "/path/to/your/input/file.csv" with the actual path to your input CSV file
data <- read.csv("/path/to/your/input/file.csv")

# Identify and remove empty columns
empty_cols <- which(colSums(is.na(data)) == nrow(data))
data <- data[, -empty_cols]

# Save the modified data to a new CSV file
# Replace "/path/to/your/output/modified_file.csv" with the desired output path and file name
write.csv(data, "/path/to/your/output/modified_file.csv", row.names = FALSE)
```
## 2. Data Visualisation
- Distribution Plot of Samples:
This code snippet leverages the ggplot2 and maps libraries to create a world map illustrating the geographical distribution of samples across the globe. The data is merged with a world map using inner_join(), and the plot is customized using various ggplot2 functions. Finally, the plot is saved as a PNG image.
```R
# Creating a distribution plot of samples across the world from the provided dataset 
# Load required libraries
library(ggplot2)
library(maps)

# Read the CSV file into a data frame
# Replace "/path/to/your/input/file.csv" with the actual path to your input CSV file
data <- read.csv("/path/to/your/input/file.csv")

# Create a world map
world <- map_data("world")

# Check for missing or non-unique country names in the data
if (any(is.na(data$Country))) {
  stop("There are missing values in the 'Country' column.")
}

if (any(duplicated(data$Country))) {
  stop("There are non-unique values in the 'Country' column.")
}

# Merge the data frame with the world map
data_map <- inner_join(world, data, by = c("region" = "Country"))

# Create the plot using ggplot2
plot <- ggplot() +
  geom_polygon(data = data_map, aes(x = long, y = lat, group = group, fill = SampleID)) +
  coord_map("miller") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_viridis_c() +
  labs(title = "Distribution of Samples on World Map")

# Save the plot as a PNG file
# Replace "/path/to/your/output/sample_map.png" with the desired output path and file name
ggsave("/path/to/your/output/sample_map.png", plot, width = 10, height = 6, dpi = 300)
```
**Analyzing Prevalence by Species and Year: Key Steps**
- Load required libraries: ggplot2 and dplyr.
- Read CSV data: Ensure correct file path and column names.
- Check data quality: Verify for missing and duplicate values in 'Year' and 'n' columns.
- Create line plot: Use ggplot() to define aesthetics (x, y, group, color), add lines and points, set labels, and apply a theme.
- Display and save: Print the plot and save it as a PNG file.
```R
# To generate a plot for prevalence of species by year 
# Load required libraries
library(ggplot2)
library(dplyr)

# Read the CSV file into a data frame
# Replace "/path/to/your/input/file.csv" with the actual path to your input CSV file
data_counts <- read.csv("/path/to/your/input/file.csv")

# Check for missing or non-unique values in the data
if (any(is.na(data_counts$Year)) || any(is.na(data_counts$n))) {
  stop("There are missing values in the 'Year' or 'n' columns.")
}

if (any(duplicated(data_counts$Year))) {
  stop("There are non-unique values in the 'Year' column.")
}

# Create the line plot using ggplot2
line_plot <- ggplot(data_counts, aes(x = Year, y = n, group = Species, color = Species)) +
  geom_line() +
  geom_point() +
  labs(title = "Prevalence of ESKAPE pathogens by species and year",
       x = "Year",
       y = "Number of Samples",
       color = "Species") +
  theme_minimal()

# Display the line plot
print(line_plot)

# Save the plot as a PNG file
# Replace "/path/to/your/output/line_plot.png" with the desired output path and file name
ggsave("/path/to/your/output/line_plot.png", line_plot, width = 10, height = 6, dpi = 300)
```

## 3. Analyzing MDR Isolate Count
1. Data Preprocessing:

- Load Libraries: Import the dplyr library for data manipulation.
- Read Raw Data: Load the CSV file containing your data using read.csv().
- Convert Resistance Labels: Define a function to replace "Resistant", "Susceptible", and "Intermediate" with "R", "I", and "S" respectively.
- Apply Conversion: Use mutate and across to automatically apply the conversion function to relevant columns in the data frame.
- Save Modified Data (Optional): If desired, save the updated data to a new CSV file. <br>
2. MDR Analysis with AMR Package:

- Call mdro Function: Use the mdro function from the AMR package to calculate MDR isolate counts.
Specify Parameters:
- x: Provide the dataframe containing antibiotic columns like AMX or amox.
- guideline: Set the guideline for MDR determination (e.g., "CMI2012").
- col_mo: Specify the column name for the names or codes of microorganisms.
- info: Set interactive() to print progress to the console.
- pct_required_classes: Adjust the minimum required percentage of antimicrobial classes for an isolate to be considered MDR.
- combine_SI: If TRUE, combines "S" and "I" into one category, considering only "R" as resistant.
- verbose: Set to FALSE for less verbose output.
- only_sir_columns: If FALSE, detects all antibiotic columns that have been transformed to "SIR" format. <br>
3. Generating MDR Plot:

- Load Libraries: Import ggplot2 and tidyverse for data visualization.
- Prepare Data: Create a data frame with pathogen names, total counts, and MDR counts. Calculate the MDR proportion.
- Create Grouped Bar Plot:
- Use ggplot() to create the plot.
- Set x to the reordered pathogen names for better visualization.
- Use geom_bar() to create bars for total counts and MDR counts, using different colors for each.
- Add text labels for the proportion of MDR isolates.
- Customize the plot with titles, labels, and theme settings as desired.
```R
# preparing the data for MDR analysis using the AMR package 
# First convert the entries in the dataset columns into R,I,S format from Resistant,Intermediate,Susceptible format 
# Load required libraries
library(dplyr)

# Read the raw data from CSV file
# Replace "/path/to/your/input/raw_data.csv" with the actual path to your input CSV file
raw_data <- read.csv("/path/to/your/input/raw_data.csv")

# Function to replace "Resistant", "Susceptible", and "Intermediate" with "R", "S", and "I" respectively
replace_resistant_susceptible_intermediate <- function(x) {
  x <- gsub("Resistant", "R", x)
  x <- gsub("Susceptible", "S", x)
  gsub("Intermediate", "I", x)
}

# Automatically detect and replace columns containing "Resistant", "Susceptible", and "Intermediate"
raw_data <- raw_data %>%
  mutate(across(where(is.character), replace_resistant_susceptible_intermediate))

# Save the updated data to a new CSV file if needed
# Replace "/path/to/your/output/modified_data.csv" with the desired output path and file name
write.csv(raw_data, "/path/to/your/output/modified_data.csv", row.names = FALSE)



# Determine the MDR isolates count from the prepared dataset usinf AMR package https://msberends.github.io/AMR/reference/mdro.html

mdro(
  # A dataframe antibiotic columns like AMX or amox.
  x = NULL, 
  # A specefic supported National/International/Custom guideline to be followed 
  guideline = "CMI2012",
  # column name of the names or codes of the microorganisms
  col_mo = NULL,
  # a logical to indicate whether progress should be printed to the console
  info = interactive(),
  # minimal required percentage of antimicrobial classes that must be available per isolate
  pct_required_classes = 0.5,
  # a logical to indicate whether all values of S and I must be merged into one, so resistance is only considered when isolates are R, not I. 
  combine_SI = TRUE,
  # a logical to turn Verbose mode on and off (default is off).
  verbose = FALSE,
  # a logical to indicate whether only antibiotic columns must be detected that were transformed to class sir
  only_sir_columns = FALSE,
  ...
)



# Generating a plot for the MDR analysis 
# Load required libraries
library(ggplot2)
library(tidyverse)

# Example data: Replace these with your actual data
pathogen_names <- c("E. faecium", "S. aureus", "K. pneumoniae", "A. baumannii", "P. aeruginosa", "Enterobacter spp")
total_counts <- c(16661, 147921, 88712, 38036, 94460, 48041)
mdr_counts <- c(5710, 1596, 27327, 23077, 17165, 7549)

# Combine vectors into a data frame
mdr_data <- data.frame(pathogen = pathogen_names, total_count = total_counts, mdr_count = mdr_counts)

# Calculate the proportion of MDR isolates relative to total isolates
mdr_data$mdr_proportion <- mdr_data$mdr_count / mdr_data$total_count

# Create the grouped bar plot with text labels
ggplot(mdr_data, aes(x = reorder(pathogen, total_count), y = total_count, fill = "Total")) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_bar(aes(y = mdr_count, fill = "MDR"), stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(mdr_proportion * 100), "%")),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
  labs(title = "MDR Isolate Count Among Total Isolates for ESKAPE Pathogens",
       x = "Pathogen",
       y = "Count",
       fill = "Isolate Type") +
  scale_fill_manual(values = c("Total" = "gray", "MDR" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"))
```
## 4. Predicting Future Antibiotic Resistance Trends <br>
This code snippet explores the possibility of predicting future resistance patterns using the AMR package. Here's a breakdown of the steps involved:

**Required Libraries** <br>

dplyr for data manipulation tasks.
ggplot2 for creating visualizations.
AMR for resistance prediction functionalities.

**Data Loading:** <br>

- The code reads a CSV file containing your data into a data frame (data). Remember to replace the path ("/path/to/your/input/your_file_path.csv") with the actual location of your data file.
- Selecting Antibiotic for Prediction:
- The code defines the specific antibiotic (antibiotic_column) for which resistance prediction will be performed. 
- Replace "Your_Antibiotic_Column_Name" with the actual column name in your data that holds antibiotic resistance information.

**Resistance Prediction:** <br>

- The resistance_predict function from the AMR package is used to predict resistance patterns for the chosen antibiotic.
- It takes the data frame (data), specifies the column containing date information (col_date), the antibiotic column (col_ab), and sets the model type (model = "binomial") for this example.
- The output (predict_antibiotic) is a data frame containing predictions.

**Visualizing Predictions:**<br>

- The code offers two approaches to visualize the predictions:
1. plot(predict_antibiotic) generates a basic plot based on the AMR package's default settings.
2. ggplot_sir_predict(predict_antibiotic) leverages ggplot2 for a more customized plot displaying the predicted proportions of Susceptible, Intermediate, and Resistant categories for the chosen antibiotic over time.
```R
# Future antimicrobial resistance prediction using regression models 
# AMR package (https://msberends.github.io/AMR/articles/resistance_predict.html)
# Load the necessary libraries
library(dplyr)
library(ggplot2)
library(AMR)

# Read the CSV file into a data frame
# Replace "/path/to/your/input/your_file_path.csv" with the actual path to your input CSV file
data <- read.csv("/path/to/your/input/your_file_path.csv")

# Specify the name of the antibiotic for resistance prediction
# Replace "Your_Antibiotic_Column_Name" with the actual column name containing antibiotic resistance information
antibiotic_column <- "Your_Antibiotic_Column_Name"

# Perform resistance prediction for the specified antibiotic
# Replace "Year" and "Your_Antibiotic_Column_Name" with the actual column names for date and antibiotic resistance information
predict_antibiotic <- data %>%
  resistance_predict(
    col_date = "Year",
    col_ab = antibiotic_column,
    model = "binomial"
  )

# Plot the resistance prediction results
plot(predict_antibiotic)

# Plot the predicted SIR (Susceptible, Intermediate, Resistant) proportions
ggplot_sir_predict(predict_antibiotic)
```