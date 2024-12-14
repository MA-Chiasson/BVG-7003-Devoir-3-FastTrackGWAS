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

# 2. Set the working directory to the directory where the script is located
# This automatically sets the working directory to the folder where the script is located, no need for modification.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 2.1. Print the working directory to confirm it's set correctly
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

# Optional files:
#kin_file <- "data/kinship.txt"  # Kinship file (optional)
#pca_file <- "data/pca.txt"  # PCA file (optional)

# 3.1. Loading phenotypic data
pheno <- read.table(pheno_file, sep = "\t", header = TRUE)  # Loading phenotype data
str(pheno)
head(pheno)

# 3.2. Loading genotype data
geno <- read.csv(hmp_file, sep = "\t", header = TRUE, check.names = FALSE)  # Loading genotype data
str(geno)
head(geno)

# 4. Filtering SNPs with too low MAF (default <5%), as well as rows containing too many NAs (default >10%)
# 4.1. For hapmap genotype file format
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

# Apply the function to your hapmap genotype file
filtered_hmp <- filter_hapmap(geno)

write.table(filtered_hmp, "data/filtered_hmp.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#OR
# 4.2. For VCF genotype file format

filter_vcf <- function(geno_vcf, freq_threshold = 5, na_threshold = 10) {
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


# Conversion HapMap -> format MVP
MVP.Data(
  fileHMP = "data/filtered_hmp.txt",
  filePhe = "data/Phenotype_African.txt",
  fileKin=TRUE,
  filePC=TRUE,
  out = "data/mvp_hmp"
)



# 5. Running the GWAS analysis
genotype <- attach.big.matrix("data/mvp_hmp.geno.desc")
phenotype <- read.table("data/mvp_hmp.phe", header = TRUE)
map <- read.table("data/mvp_hmp.geno.map", header = TRUE)

# Running a simple GWAS with GLM
results <- MVP(
  phe = phenotype,
  geno = genotype,
  map = map,
  method = c("GLM"),
  nPC.GLM = 3,  # Number of principal components for GLM
  threshold = 0.05
)

# Extract results and save them
map <- results$map
glm_results <- as.data.frame(results$glm.results)
colnames(glm_results) <- c("Effect", "SE", "P.value")

combined_results <- cbind(map, glm_results)
write.table(combined_results, "Test_GWAS_results.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Prepare the data for visualization
vis_data <- combined_results[, c("SNP", "CHROM", "POS", "P.value")]
colnames(vis_data) <- c("SNP", "Chromosome", "Position", "P.value")

MVP.Report(
  vis_data,
  plot.type = c("qq", "manhattan"),
  threshold = 0.05 / nrow(vis_data) # Bonferroni threshold
)
