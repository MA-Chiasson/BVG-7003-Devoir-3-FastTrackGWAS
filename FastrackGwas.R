# Introduction
# This GWAS analysis uses the rMVP package to identify SNPs associated with a given phenotypic trait. 

# 1. Installation and loading of libraries
# 1.1. Install the libraries if necessary
if (!requireNamespace("rMVP", quietly = TRUE)) install.packages("rMVP")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("mgsub", quietly = TRUE)) install.packages("mgsub")
if (!requireNamespace("bigmemory", quietly = TRUE)) install.packages("bigmemory")

# 1.2. Load the libraries
library(rMVP)
library(ggplot2)
library(data.table)
library(dplyr)
library(mgsub)
library(bigmemory)

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

# Conversion HapMap -> format MVP
MVP.Data(
  fileHMP = hmp_file,
  filePhe = "data/Phenotype_African.txt",
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
