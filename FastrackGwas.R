# Introduction
# This GWAS analysis uses the rMVP package to identify SNPs associated with a given phenotypic trait. 

# 1. Installation and loading of libraries
# 1.1. Install the libraries if necessary
if (!requireNamespace("rMVP", quietly = TRUE)) install.packages("rMVP")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
if (!requireNamespace("mgsub", quietly = TRUE)) install.packages("mgsub")
if (!requireNamespace("bigmemory", quietly = TRUE)) install.packages("bigmemory")
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")

# 1.2. Load the libraries
library(rMVP)
library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)
library(mgsub)
library(bigmemory)
library(vcfR)
library(outliers)


# 2. Set the working directory to the directory where the script is located
# This automatically sets the working directory to the folder where the script is located, no need for modification.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Print the working directory to confirm it's set correctly
cat("Working directory has been set to:\n")
cat(getwd(), "\n")


# 3. Loading data
# Please make sure that the following files are uploaded and placed in the "data" folder:
# - Phenotype_African.txt (mandatory)
# - African_SNPs.hmp.txt (mandatory)
# - Kinship matrix file (optional, if you have it)
# - PCA matrix file (optional, if you have it)
# If you do not have the kinship or PCA files, the script will automatically generate them.

pheno_file <- "data/Phenotype_African.txt"  # File with phenotype data (mandatory)
hmp_file <- "data/African_SNPs.hmp.txt"  # File with genotype data (mandatory)

# 3.1. Loading genotype data
geno <- read.csv(hmp_file, sep = "\t", header = TRUE, check.names = FALSE)  # Loading genotype data
str(geno)
head(geno)

# 3.2. Loading phenotypic data
pheno <- read.table(pheno_file, sep = "\t", header = TRUE)  # Loading phenotype data
str(pheno)
head(pheno)

## 3.3 Optional input:
# If you do not have the kinship or PCA files, the script will automatically generate them!
# read from file and convert my result kin
MVP.Data.Kin("data/mvp.kin.txt", out="data/mvp", maxLine=1e4) # Kinship file (optional)
# read from file and convert my result PCA 
#rMVP.Data.PC("data/mvp.pc.txt", out='data/mvp', sep='\t') # PCA file (optional)


# 4. Filtering SNPs with too low MAF (default <5%), as well as rows containing too many NAs (default >10%)
# This filtering process is performed separately for different input formats:
# - HapMap format
# - VCF format

# 4.1. For hapmap genotype file format
# 4.1.1. Function to filter hapmap genotype data
filter_hapmap <- function(hapmap_data, freq_threshold = 5, na_threshold = 10) {
  # Function to calculate the frequency of each letter
  calculate_frequencies <- function(genotypes) {
    # Remove NA values
    genotypes <- genotypes[!is.na(genotypes)]
    # Calculate the frequency of each genotype
    freq_table <- table(genotypes)
    # Sort frequencies in descending order
    sorted_freq <- sort(freq_table, decreasing = TRUE)
    return(sorted_freq)
  }
  
  # Apply the frequency calculation function to each row
  hapmap_data <- hapmap_data %>%
    rowwise() %>%
    mutate(
      freq = list(calculate_frequencies(c_across(12:ncol(hapmap_data)))),  # Apply frequency calculation across relevant columns
      na_count = sum(is.na(c_across(12:ncol(hapmap_data)))),  # Count NAs in the relevant columns
      total_count = ncol(hapmap_data) - 11,  # Calculate the total count of genotypes (after column 11)
      na_percentage = (na_count / total_count) * 100  # Calculate NA percentage
    ) %>%
    ungroup()  # Ungroup the data to avoid issues with rowwise()
  
  # Filter rows based on the NA percentage and the frequency of the second most frequent genotype
  filtered_hmp <- hapmap_data %>%
    filter(na_percentage < na_threshold) %>%  # Filter rows where NA percentage is less than threshold
    filter(map_lgl(freq, ~ {
      if (length(.x) > 1) {  # Ensure there is more than one genotype to check the second most frequent
        second_most_freq <- .x[2]  # Get the second most frequent genotype
        second_most_freq_value <- as.numeric(second_most_freq)  # Convert it to numeric
        total_genotypes <- sum(.x)  # Get the total number of genotypes
        freq_percentage <- (second_most_freq_value / total_genotypes) * 100  # Calculate the percentage
        return(freq_percentage > freq_threshold)  # Return TRUE if the percentage is above the threshold
      } else {
        return(FALSE)  # If only one genotype is present, return FALSE
      }
    })) %>%
    select(-freq, -na_count, -total_count, -na_percentage)  # Drop temporary columns used for filtering
  
  return(filtered_hmp)
}


# 4.1.2. Apply the function to your hapmap genotype file
filtered_hmp <- filter_hapmap(geno)

write.table(filtered_hmp, "data/filtered_hmp.txt", sep = "\t", row.names = FALSE, quote = FALSE)


############OR#################
# 4.2. For VCF genotype file format
# 4.2.1. Function to filter VCF genotype data
filter_vcfmap <- function(geno_vcf, freq_threshold = 5, na_threshold = 10) {
  # Calculate the MAF for each variant
  maf_values <- maf(geno_vcf)
  
  # Convert the MAF result to a data frame
  maf_df <- as.data.frame(maf_values)
  
  # Add variant IDs (CHROM and POS) to the MAF data frame
  maf_df <- maf_df %>%
    mutate(CHROM = geno_vcf@fix[, "CHROM"],
           POS = as.numeric(geno_vcf@fix[, "POS"]))
  
  # Calculate the percentage of NA values for each row
  maf_df <- maf_df %>%
    rowwise() %>%
    mutate(na_percentage = (`NA` / nAllele) * 100) %>%
    ungroup()
  
  # Filter rows based on the MAF and NA percentage
  filtered_data <- maf_df %>%
    filter(na_percentage < na_threshold) %>%
    filter(Frequency > (freq_threshold / 100))
  
  # Get the CHROM and POS of the filtered rows
  filtered_variants <- filtered_data %>%
    select(CHROM, POS)
  
  # Subset the original VCF object to keep only the filtered rows
  filtered_vcf <- geno_vcf[geno_vcf@fix[, "CHROM"] %in% filtered_variants$CHROM & 
                           as.numeric(geno_vcf@fix[, "POS"]) %in% filtered_variants$POS, ]
  
  return(filtered_vcf)
}

# 4.2.2. Apply the function to your VCF genotype file
filtered_vcf <- filter_vcfmap(geno)

write.table(filtered_vcf, "data/filtered_vcf.txt", sep = "\t", row.names = FALSE, quote = FALSE)
##################


# 5. Phenotypic Data Filtering
# 5.1. Function to remove outliers using IQR
process_phenotypic_data <- function(phenotypic_data, na_threshold = 10) {
  # Identify numeric columns (excluding the first column which is assumed to be the sample identifier)
  numeric_cols <- names(phenotypic_data)[sapply(phenotypic_data, is.numeric)]
  
  # If there are no numeric columns, return the original data
  if (length(numeric_cols) == 0) {
    stop("No numeric columns found in the dataset.")
  }
  
  # Process each numeric column
  for (column_name in numeric_cols) {
    # Remove rows with NA values in the specific column
    phenotypic_data <- phenotypic_data[!is.na(phenotypic_data[[column_name]]), ]
    
    # Calculate NA percentage for numeric columns
    phenotypic_data <- phenotypic_data %>%
      rowwise() %>%
      mutate(
        na_percentage = (sum(is.na(c_across(all_of(numeric_cols)))) / length(numeric_cols)) * 100
      ) %>%
      ungroup()
    
    # Filter rows with NA percentage less than the threshold
    phenotypic_data <- phenotypic_data %>% 
      filter(na_percentage < na_threshold)
    
    # Remove outliers using IQR for the current column
    Q1 <- quantile(phenotypic_data[[column_name]], 0.25, na.rm = TRUE)
    Q3 <- quantile(phenotypic_data[[column_name]], 0.75, na.rm = TRUE)
    IQR_value <- Q3 - Q1
    
    # Define lower and upper bounds
    lower_bound <- Q1 - 1.5 * IQR_value
    upper_bound <- Q3 + 1.5 * IQR_value
    
    # Filter data within bounds for the current column
    phenotypic_data <- phenotypic_data %>%
      filter(!!sym(column_name) >= lower_bound & !!sym(column_name) <= upper_bound)
  }
  
  # Drop temporary columns
  phenotypic_data <- phenotypic_data %>% select(-na_percentage)
  
  return(phenotypic_data)
}

# 5.2 Process the phenotypic data for all numeric columns (except the first column)
filtered_pheno <- process_phenotypic_data(pheno)  # Automatically filters all numeric columns

# Save the filtered data to a file
write.table(filtered_pheno, "data/filtered_pheno.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Check the results
head(filtered_pheno)  # Display the first few rows of the filtered data


# 6. Conversion HapMap -> format MVP
MVP.Data(
  fileHMP = "data/filtered_hmp.txt",  # Filtered HapMap file
  filePhe = "data/filtered_pheno.txt",  # Filtered phenotypic data file
  fileKin = TRUE,   # Set to TRUE if kinship file is not available, otherwise set to FALSE
  filePC = FALSE,    # Set to TRUE if PCA file is not available, otherwise set to FALSE
  out = "data/mvp_hmp"  # Output in MVP format
)

# Note: The fileKin and filePC parameters are related to step 3.3. 
# If you already have the kinship and PCA files (e.g., "data/mvp.kin.txt" and "data/mvp.pc.txt"),
# set these parameters to FALSE to avoid regenerating them.
# If the files are not available, leave these as TRUE, and the script will create them automatically.


# 7. Input file if want Gwas
genotype <- attach.big.matrix("data/mvp_hmp.geno.desc")
phenotype <- read.table("data/mvp_hmp.phe",head=TRUE)
map <- read.table("data/mvp_hmp.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("data/mvp_hmp.kin.desc")
#Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("data/mvp.hmp.pc.desc"))


# 8. GWAS Analysis for Single and Multiple Traits
# In this section, we perform two GWAS models simultaneously: one for a single trait and another for multiple traits.
# The key difference between these two is that the single-trait model analyzes one trait at a time, 
# whereas the multiple-trait model analyzes several traits simultaneously.
# The outputs for each will be saved in separate folders.


# Create output folders if they don't exist
if (!dir.exists("output")) {
  dir.create("output")
}
if (!dir.exists("output/single_trait_results")) {
  dir.create("output/single_trait_results")
}
if (!dir.exists("output/multi_trait_results")) {
  dir.create("output/multi_trait_results")
}

# Verify directories were created
cat("Directories created successfully:\n")
print(list.dirs("output"))

# 8.1 GWAS for a Single Trait
cat("Performing GWAS for a single trait:\n")
setwd("output/single_trait_results")  

imMVP_single_trait <- tryCatch({
  MVP(
    phe = phenotype,          # NA is acceptable in phenotype
    geno = genotype,
    map = map,
    nPC.GLM = 5,              # Number of PCs for GLM
    nPC.MLM = 3,              # Number of PCs for MLM
    nPC.FarmCPU = 3,          # Number of PCs for FarmCPU
    maxLine = 10000,          # Smaller value reduces memory cost
    vc.method = "BRENT",      # Method for MLM
    method.bin = "static",    # "FaST-LMM", "static" (only works for FarmCPU)
    threshold = 0.05,
    method = c("GLM", "MLM", "FarmCPU"),
    file.output = c("pmap", "pmap.signal", "plot", "log")  
  )
}, error = function(e) {
  cat("Error during single trait GWAS: ", e$message, "\n")
  NULL
})

setwd("../../")
cat("Checking if single trait results were saved:\n")
if (!is.null(imMVP_single_trait)) {
  cat("Single trait results files:\n")
  print(list.files("output/single_trait_results"))
} else {
  cat("No results generated for single trait GWAS.\n")
}

# 8.2 GWAS for Multiple Traits
cat("Performing GWAS for multiple traits:\n")
setwd("output/multi_trait_results") 

for (i in 2:ncol(phenotype)) {
  imMVP_multi_trait <- tryCatch({
    MVP(
      phe = phenotype[, c(1, i)],  # Using the first column as trait identifier and i-th trait
      geno = genotype,
      map = map,
      nPC.GLM = 5,
      nPC.MLM = 3,
      nPC.FarmCPU = 3,
      maxLine = 10000,
      vc.method = "BRENT",
      method.bin = "static",
      threshold = 0.05,
      method = c("GLM", "MLM", "FarmCPU"),
      file.output = c("pmap", "pmap.signal", "plot", "log")  
    )
  }, error = function(e) {
    cat("Error during multiple trait GWAS (trait ", i, "): ", e$message, "\n")
    NULL
  })
  
  # Check if output files are generated for multiple traits
  cat("Checking if results for trait ", i, " were saved:\n")
  if (!is.null(imMVP_multi_trait)) {
    cat("Multiple trait results files for trait ", i, ":\n")
    print(list.files("."))
  } else {
    cat("No results generated for trait ", i, ".\n")
  }
  
  # Clean up memory after each iteration
  gc()
}

# Return to main directory after processing
setwd("../../")
cat("GWAS for multiple traits completed.\n")

