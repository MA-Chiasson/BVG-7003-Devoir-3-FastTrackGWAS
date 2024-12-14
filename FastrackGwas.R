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

# Define the working directory
setwd("C:/Users/tony7*OneDrive/Git/BVG-7003-Devoir-3-FastTrackGWAS") # To modify
# setwd("C:/path/to/your/data") # Replace with your actual path # I don't know if you'll want to change it only as a finishing touch up...

# cat("Please set your working directory to the folder containing your data files.\n") # Alternative interactive way to set working directory
# setwd(readline(prompt = "Enter the path to your working directory: "))


# 3. Loading data
pheno_file <- "data/Phenotype_African.txt"
hmp_file <- "data/African_SNPs.hmp.txt"

# Load and modify the phenotype

#3.1 Loading phenotypic data
pheno <- read.table(pheno_file, , sep = "\t", header = TRUE) 
str(pheno)
head(pheno)

geno <- read.csv(hmp_file, sep = "\t", header = TRUE)
str(geno)
head(geno)


pheno$Sample <- gsub("TGx ", "", pheno$Sample)

# Saving modified phenotype table
write.table(pheno, "data/pheno_modified.txt", sep = "\t", row.names = FALSE, quote = FALSE)


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

# Example usage
# Assuming geno_vcf is your VCF object
# filtered_vcf <- filter_vcf(geno_vcf, freq_threshold = 5, na_threshold = 10)

# 4.3. Should we create a MAF filtering function for other formats as well ?
# For PLINK and binary formats however this would mean building entirely new functions which may not be necessary for the assignment
# However if we do not provide any function for other formats, we simply need to specify that this pipeline only supports hmp and maybe vcf formats, (which would seem fine to me, as they are the most common anyway).
# I tried combining the functions in 4.1 and 4.2, but I get persistent errors which me and copilot aren't able to solve, so I'll leave it like in 2 separate function for now, but if you guys feel it would be necessary to merge them, feel free to try !


# Conversion HapMap -> format MVP
MVP.Data(
  fileHMP = filtered_hmp,
  filePhe = "data/Phenotype_African.txt",
  fileKin=TRUE,
  filePC=TRUE,
  out = "mvp_hmp"
)

# 4. Filtering SNPs with MAF < 0.05  
# Load the genotype data
genotype_data <- attach.big.matrix("mvp_hmp.geno.desc")

# Calculate MAF using the bigmemory object
maf_data <- apply(as.matrix(genotype_data), 2, function(x) {
  table_x <- table(x)  # Count the alleles
  maf <- min(table_x) / sum(table_x)  # Calculate the minor allele frequency
  return(maf)
})

# Print the first few values of maf_data to inspect
print(head(maf_data))

# Filter out the SNPs with MAF < 0.05
filtered_genotype_data <- genotype_data[, maf_data >= 0.05]

# Check the size of the data after filtering
print(dim(filtered_genotype_data))


# 5. Running the GWAS analysis
genotype <- attach.big.matrix("mvp_hmp.geno.desc")
phenotype <- read.table("mvp_hmp.phe", header = TRUE)
map <- read.table("mvp_hmp.geno.map", header = TRUE)

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
