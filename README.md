# BVG-7003-Assignment-3-FastTrackGWAS

## Table of Contents
1. [Introduction](#introduction)
2. [Installation of R and RStudio](#installation-of-r-and-rstudio)
   1. [Installing R](#1-installing-r)
   2. [Installing RStudio](#2-installing-rstudio)
   3. [Launching RStudio](#3-launching-rstudio)
3. [Script execution](#script-execution)
   1. [Input and Output File Formats](#1-input-and-output-file-formats)
   2. [Dependencies](#2-dependencies)
      1. [Installation of the packages](#2-1-installation-of-the-packages)
      2. [Loading the packages](#2-2-loading-the-packages)
   3. [Define working directory](#3-define-working-directory)
4. [Loading Data in R](#4-loading-data-in-r)
   1. [Loading phenotypic data](#4-1-loading-phenotypic-data)
   2. [Loading Genotypic Data](#4-2-loading-genotypic-data)
5. [Data Preprocessing](#5-data-preprocessing)
   1. [Formatting Genotype and Phenotype Files](#5-1-formatting-genotype-and-phenotype-files)
   2. [Quality Control](#5-2-quality-control)
      1. [SNP Filtering with MAF < 0.05](#5-2-1-snp-filtering-with-maf-005)
   3. [Generating Covariates](#5-3-generating-covariates)

***Complete table of contents once readme has its final structure***

---
## Introduction

Genome-wide association Studies (GWAS) are a powerful tool for uncovering the genetic basis of traits, and their application in agricultural science is pivotal for crop improvement.

We used GWAS to investigate the genetic foundation of traits associated with symbiotic nitrogen fixation (SNF) in African soybeans, identifying potential candidate genes. This research offers crucial insights into SNF-related traits and will expedite progress in SNF breeding initiatives.
The steps provided in this repository guide users through setting up the required tools and efficiently performing a GWAS analysis.

This project employs **R** and **RStudio** as the core computational environment to enable reproducibility. R is renowned for its statistical computing capabilities and robust library ecosystem, making it suitable for GWAS workflows. This repository demonstrates how to install and configure R and RStudio, manage dependencies, and execute the GWAS analysis pipeline.

The guide is designed to help users of all experience levels replicate the analysis and adapt it to other datasets or traits of interest.

--- 

## 1. Installation of R and RStudio

To perform a GWAS analysis using R, follow the steps below to install R and RStudio on your computer:

### 1.1. **Installing R**
R is the programming language required for the analysis:
1. Visit the official R website: [https://cran.r-project.org/](https://cran.r-project.org/).
2. Click on **Download R for [your operating system]** (Windows, macOS, or Linux).
3. Download the appropriate version for your system:
   - For Windows: Click on *base* and then download the installer.
   - For macOS: Download the version compatible with your system.
   - For Linux: Follow the instructions specific to your distribution.
4. Run the downloaded file and follow the installation instructions.

### 1.2. **Installing RStudio**
RStudio is a user-friendly interface for R:
1. Go to the RStudio website: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/).
2. Click on **RStudio Desktop** (free version).
3. Download the installer for your operating system.
4. Install RStudio by running the downloaded file.

### 1.3. **Launching RStudio**
1. Once installed, open RStudio.
2. Ensure that R is correctly configured with RStudio (RStudio automatically detects R if it is installed).
---
## 2. **Data Download and Preparation**
This analysis is based on the dataset provided in the paper:
**"Genome-Wide Association Analyses Reveal the Genetic Basis of Symbiotic Nitrogen Fixation in African Soybean"**
- **DOI:** [10.1007/s00122-019-03499-7](https://link.springer.com/article/10.1007/s00122-019-03499-7)
- **Data link:** [Figshare Dataset](https://figshare.com/projects/GenomeWide_Association_Analyses_Reveal_the_Genetic_Basis_of_Symbiotic_Nitrogen_Fixation_in_African_Soybean/65885)
### Test Data Used

For the test data in this analysis, the following were selected:

- **Genotypic data:** Information for all chromosomes.
- **Phenotypic data:** The trait SDWi/SDWf (SDW ratio) represents shoot dry weight in inoculated (SDWi) and fertilized (SDWf).
--- 

## 3. Script execution
### 3.1. Input and Output File Formats
#### 3.1.1. Input
#### 3.1.2. Output
**REPLACE BY THE RESULT THE USER SHOULD HAVE**

### 3.2. **Dependencies**
#### 3.2.1. Installation of the packages
Open the script and start to install the following R packages which are required to execute this analysis:
- **rMVP**: A memory-efficient, visualization-enhanced, and parallel-accelerated tool for GWAS.
- **ggplot2**: For creating high-quality graphics.
- **data.table**: For fast and flexible data manipulation.
- **dplyr**: For data manipulation with a consistent grammar.
- **purrr**: For working with functions and vectors.
- **mgsub**: For multiple, simultaneous string substitutions.
- **bigmemory**: For managing massive matrices.
- **vcfR**: For pre-treatment of VCF format genotype data.
- **outliers**: For outlier identification.
  **Adjust the remainder of this list once code is final**

Dependencies for these packages will be automatically installed when executing the script, ensuring a smooth setup process.
For this purpose, the installation of packages is conditional using the following function:

```r
if (!requireNamespace("package", quietly = TRUE)) install.packages("package")
```

This checks whether each package is already installed before attempting to install it. If the package is not found, it is installed. This approach ensures that packages are only installed when necessary, preventing redundant installations during each script iteration, which can slow down execution.

On the other hand, if a single line of code was used to install all the packages at once, like this:

```r
install.packages(c("rMVP", "ggplot2", "data.table", "dplyr", "mgsub", "bigmemory"))
```
This would lead to the installation of packages every time the script is run, even if they are already installed, which is inefficient. The if condition ensures that packages are only installed when required, avoiding unnecessary installation.

#### 3.2.2. Loading the packages
Once they've been installed, it's time to load them. To do this, use the library() command, inserting the package name between the brackets, as in the example below. 
```r
library(rMVP)
```
---
### 4. Define working directory
**Defining a working directory in R offers several key benefits:**

- **File Management**: It helps centralize project files (data, scripts, results) in one location, avoiding confusion.

- **Simplified File Access**: With a set working directory, you can use relative file paths (e.g., read.csv("data.csv")) instead of full paths.

- **Easier Data Import/Export**: Files can be easily imported or exported without needing to specify full file paths.

- **Reduced Path Errors**: Setting a working directory minimizes errors related to incorrect file paths, especially in team projects.

- **Project Organization**: It encourages better project structure, making it easier to manage and maintain files.

#### How to Set a Working Directory
Usually, the WD is the folder where our script and data are stored, and where we want to add the results of our analyses.
A very basic way of setting a working directory is using setwd():

```r
setwd("C:/Users/YourPath/ToYour/Folder")
```
Check the current working directory with getwd():
```r
getwd()
```

To simplify the use of the demo script, we have automated this step to set the working directory to the folder where you have placed the script with the following commands:

```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```
And also an interactive way to verify it has been done correctly :

```r
cat("Working directory has been set to:\n")
cat(getwd(), "\n")
```

***Adjust as the script gets final***

---

### 5. Loading data in R
There are several good ways to load data into R, depending on the data format. Below, we explain how to load both phenotypic and genotypic data in the R environment using specific commands.

#### 5.1. Loading genotype and phenotype data
To load data in R, we typically use `read.csv()` or `read.table()` functions to import the data. The key difference is that `read.csv()` is used when the data is separated by commas, while `read.table()` is more general and allows specifying different delimiters (e.g., tab-separated data). 

In our case, we use the following commands to load the phenotypic data:
```r
geno <- read.csv(hmp_file, sep = "\t", header = TRUE, check.names = FALSE)
```
```r
pheno <- read.table(pheno_file, sep = "\t", header = TRUE)
```
Explanation of the parameters:
- `hmp_file` or `pheno_file`: The file paths to the genotype and phenotype data.
- `sep = "\t"`: This specifies that the data is tab-separated.
- `header = TRUE`: This indicates that the first row contains column names.
- `check.names = FALSE`: This prevents R from altering column names to make them syntactically valid.

---

### 6. Data Preprocessing
In this section, we guide users through the following steps for data preprocessing. These steps are crucial for ensuring the quality and integrity of the data used in the analysis. Preprocessing helps to clean, format, and structure the data properly, minimizing errors and biases in the final results.

#### 6.1. Formatting Genotype and Phenotype Files
**IMPORTANT NOTE:** Before starting, it is important to note that a genotype data file cannot be easily modified directly in R, as it is a special file format. Therefore, any necessary modifications to harmonize the datasets should be made to the phenotype data file, as we will see below.

To check if the format of our data is correct, we can use two commonly used commands in R: `str()` and `head()`. These functions are useful for quickly exploring the structure and content of our data before proceeding with more complex analyses.

##### 6.1.1. `str()`
The str() (structure) function is used to display the internal structure of an object in R. It provides a quick overview of the data structure, including the type of each column (e.g., numeric, factor, etc.), the number of missing values, and a sample of the first few values in each column.

**Usefulness of `str()`:**

Chcking the data types in each column (e.g., whether it's a factor, character, numeric value).
Helping detect any formatting errors, such as a numeric column being incorrectly read as a factor.
Providing a summary of the data object, which is useful for verifying everything is in order before applying transformations or analyses.
Example:

```r
str(pheno)
```
This allows you to check that the phenotype data is in the correct format before using it in the analysis.

##### 6.1.2. `head()`
The `head()` function displays the first few rows of a dataset. By default, it shows the first six rows, but this number can be adjusted. This function gives a quick look at the actual values in the first rows, which is useful for checking the quality of the data and ensuring that they have been imported correctly.

**Usefulness of `head()`**:

It allows you to quickly view the first rows of a dataset to verify that the values are consistent and correctly formatted.
It provides a quick check for missing values or incorrect data in the first few rows.
It is particularly useful for large files, as it allows you to examine a subset of the data without loading the entire file.

**Example:**
```r
head(pheno)
```
This allows you to view the first few rows of the phenotype data and check everything is in order before continuing with the analysis.

**Conclusion**
These two commands are essential tools for quickly inspecting the structure and content of your data. They help identify potential problems before the analysis process, saving a lot of time and effort upfront. Using str() and head() before starting any detailed analysis ensures that the data is properly formatted and ready to be processed.

##### 6.1.3. Example of formatting
As observed in the phenotypic data (Figure 1), sample names often begin with "TGx," whereas this is not the case in the genotypic data (Figure 2). If we were to run the MVP function as is, only 18 matching samples would be identified between the two datasets, as shown in Figure 3 below.
**Name of the sample in the phenotypic data:**  
![image](https://github.com/user-attachments/assets/3e01bf96-879e-4a7c-90ad-50e8acc23228)

**Name of the sample in the genotypic data:**  
![image](https://github.com/user-attachments/assets/eb3511ab-2b10-4182-9645-abac278ffaa2)

**MVP without modifications:**  
![image](https://github.com/user-attachments/assets/9220e96c-394f-4423-acf8-e065da1cb288)

To address this, we will remove the "TGx" prefix from the phenotypic data using the following command:
```r
pheno$Sample <- gsub("TGx ", "", pheno$Sample)
```
Although the genotypic data uses a different naming convention, where it includes an "X" and replaces "-" with ".", it may seem intuitive to make similar changes in the phenotypic data for consistency. However, this step is unnecessary because RStudio automatically modifies column names when loading the data. Specifically, it prepends an "X" to column names starting with a number and replaces "-" with ".", ensuring compatibility with R’s column naming rules. As a result, if we run the MVP function now, we will have 281 matching samples between the two datasets.

**MVP after modifications**  
![image](https://github.com/user-attachments/assets/6b18b186-145f-4939-8d67-4b8352ea95d1)


#### 6.2. Quality Control

The quality control step involves filtering loci, for instance SNPs (Single Nucleotide Polymorphisms), based on their Minor Allele Frequency (MAF). We suggest removing loci with a MAF lower than 0.05 to ensure that only common variants are included in the analysis, but other tresholds may be used. 

Loci with low MAF (for example, lower than 0.05) carry an elevated risk of having genotyping errors (i.e. the minor allele may be a genotyping error and not even exist in reality). Further, loci with alleles present at very low frequencies are more prone to statistical errors, whether false positives or false negatives. Therefore, including them doesn't contribute significantly to the discovery of meaningful associations.

A second quality control step is also to eliminate loci with an important proportion of missing data. We suggest removing loci with 10% or more missing data, but other thresholds may be used.

Loci with too much missing data can potentially induce statistical errors.

The final quality control step implemented in this pipeline is that of missing data imputation. We recommend the users to impute missing data before importing their data into the pipeline, for example using LD-kNNi, FILLIN, FSHAP or BEAGLE software.  ***MAYBE THIS SHOULD BE MOVED TO THE DATA PREPROCESSING PART***

##### 6.2.1. Filtering out SNPs with MAF < 0.05

We calculate the MAF for each SNP and filter out those with a MAF below 0.05. At the same time, this is filtering out monomorphic loci that would not be useful for a GWAS analysis. Below is the R code used for this process:

***THIS PART OF THE CODE DOESN'T WORK FOR ME, I have created a new function to do this MAF filtering (as well as data completeness filtering) for hapmap format, I think it is better to do this filtering step before moving on to the MVP format which is quite complex.***
***NOTE TO NOT FORGET TO UPDATE THIS PART OF CODE ONCE THE FILTERING FUNCTION IS DONE***
```r
# Load genotype data
genotype_data <- attach.big.matrix("mvp_hmp.geno.desc")

# Calculate the allele frequency for each SNP
# Each column represents a SNP, and each row represents an individual
maf_data <- apply(genotype_data, 2, function(x) {
  # Count the alleles
  table_x <- table(x)  
  maf <- min(table_x) / sum(table_x)  # Calculate the minor allele frequency
  return(maf)
})

# Filter out SNPs with MAF < 0.05
filtered_genotype_data <- genotype_data[, maf_data >= 0.05]
```

This code filters out SNPs with a MAF lower than 0.05, ensuring only those with a sufficient minor allele frequency are retained for further analysis.

#### 6.3. Generating Covariates
***This step can be done automatically when doing the gwas analysis (see 6.3. below). It can also be done separately with other functions from the rMVP package but it seems to me less efficient to do so. Maybe we can still document this as a separate step, but otherwise this section should go in the 6.2. or 6.3. sections below for better consistency I think.***
To account for population structure or relatedness in the data, you may generate covariates such as Principal Component Analysis (PCA) scores or a relatedness (kinship) matrix. These covariates are used to adjust for the potential confounding factors in the GWAS analysis that are kinship or population structure. Indeed, population structure and relatedness may induce non-random distribution of alleles in the sampling pool. GWAS analysis over phenotypes that happen to covary with the population structure would then result in non-relevant associations with these non-randomly distributed alleles. Including a kinship matrix or a a PCA can help reduce this risk.

Generating a kinship matrix on its own

``` r
MVP.K.VanRaden(
M,
maxLine = 5000,
ind_idx = NULL,
cpu = 1,
verbose = TRUE,
checkNA = TRUE
)
```

Generating a pca for population structure on its own

``` r
MVP.PCA(
M = NULL,
K = NULL,
maxLine = 10000,
ind_idx = NULL,
pcs.keep = 5,
cpu = 1,
verbose = TRUE
)
```

Or generating both as part of the MVP.DATA() function at the same time as when generating MVP compatible genotype and phenotype files, done through setting the arguments `fileKin` and `filePC` as true.

``` r
MVP.DATA(fileKin=TRUE, filePC=TRUE)
```
---

### 7. Execution of the GWAS analysis

#### 7.1. Loading genotype and phenotype data in the rMVP environment 
The first step of running the GWAS analysis is to load properly the genotypic and phenotypic data

Once you have correctly imported your data in the R environment (see step 4) and preprocessed it (see step 5), you are ready to load it in the rMVP environment.
The sample code below shows a simple way to load hmp format genotypic data, using the attach.big.matrix() function.
If you have a genotype map, you can also add it here.

In order to generate rMVP usable files, we first need to use the `MVP.Data()` function from the rMVP package. Two of the most common types of data used in genomics are HapMap and VCF formats. We load each type of data with different parameters.

##### 7.1.1. Generating rMVP readable files, VCF data version: 
If the genotypic data is in VCF (Variant Call Format), we use the following command:
```r
MVP.Data(fileVCF = vcf_file, filePhe = pheno_file, out = "mvp_vcf")
```
- `fileVCF`: The file path to the VCF file.
- `filePhe`: The file path to the phenotypic data.
- `out` = "mvp_vcf": The output file name (the data will be saved in this format).

##### 7.1.2. Generating rMVP readable files, HapMap data version: 
If the genotypic data is in HapMap format, we use the following command:
``` r
MVP.Data(fileHMP = hmp_file, filePhe = pheno_file, out = "mvp_hmp")
```
- `fileHMP`: The file path to the HapMap file.
- `filePhe`: The file path to the phenotypic data.
- `out` = "mvp_hmp": The output file name.

##### 7.1.3. Loading the data in rMVP environment:
Once the rMVP readable files have been generated, we use the following commands to import them in the rMVP environment:
``` r
genotype <- attach.big.matrix("mvp_hmp.geno.desc")
phenotype <- read.table("mvp_hmp.phe", header = TRUE)
map <- read.table("mvp_hmp.geno.map", header = TRUE)

```
#### 7.2. Loading covariates in MVP function (optional)
***see 7.3.?***

Although not covered in the example code, it is possible to also include other environmental variables or covariates as fixed effects to account for various experimental designs (for example, breed, sex, weight, age, diet, socioeconomic status, batch effects, genotyping platform, etc.). The inclusion of such factors is done through the arguments CV.GLM, CV.MLM, CV.FarmCPU. See [https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#advanced](url) for more details.

#### 7.3. Running MVP function

Executing the GWAS analysis is done with the MVP() function. In addition to handling the genotype, phenotype, and mapping data, the MVP() function allows to :

1. Include covariates.
       If you do not provide your own kinship matrix, population structure data, or other covariates, the MVP() function can compute them for you while running its analysis. See parameters K, CV.GLM, CV.MLM, CV.FarmCPU, nPC.GLM, nPC.MLM=3, nPC.FarmCPU.
2. Choose the analytical method and the significance threshold.
       Three analytical methods can be used for the GWAS analysis, each with their pros ans cons. They are the following : "GLM", "MLM", "FarmCPU". Note that you can execute the function with all three methods and examine the quality of each.
       The default singificance threshold (before p-value adjustment for multiple testing) is set at 0.05 
3. And others

***The chunk of code and its comments below comes from rMVP github. Not sure how much we should change this as it is already pretty clear and straigthforward... Or how we should notify that it is directly from rMVP github***
```
imMVP <- MVP(
  phe=phenotype,          #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  #K=Kinship,             #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  #CV.GLM=Covariates,     #if you have environmental covariates, please keep all 'CV.*' open
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  nPC.GLM=5,              #if you have added PCs into covariates, please keep there closed
  nPC.MLM=3,              #if you don't want to add PCs as covariates, please comment out the parameter instead of setting it to 0.
  nPC.FarmCPU=3,
  maxLine=10000,          #smaller value would reduce the memory cost
  #ncpus=10,
  vc.method="BRENT",      #only works for MLM
  method.bin="static",    # "FaST-LMM", "static" (#only works for FarmCPU)
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

```

---

### 8. Example of results obtained from MVP

#### 8.1. QQ plots 
**Add some sample plots when script is final**
#### 8.2. Manhattan plots (association strength between individual loci and the studied phenotype)
**Add some sample plots when script is final**
#### 8.3. Filtering significant SNPs from the analysis
**Does the presence of this requirement in the assignment instructions mean we need to filter the genotype file as a final step of the pipeline ?**
---

### 9. Interpretation of results

#### 9.1. Interpreting QQ plots

QQ plots are generally produced to help assess the quality of the GWAS analysis. They can help identify if there are deviations from the expected distribution, which can indicate potential issues such as population stratification, genotyping errors, or true genetic associations.

If the points on the QQ plot follow the diagonal line, it suggests that the observed p-values match the expected distribution, indicating no systematic bias.
Deviations above the line, especially at the tail, suggest an excess of small p-values, which could indicate true associations or potential confounding factors.

#### 9.2. Identifying significant SNPs from Manhattan plots

Manhattan plots display the p-values of SNPs across the genome:

Peaks: Tall peaks indicate SNPs with low p-values, suggesting significant associations.

Threshold Line: The horizontal line represents the significance threshold set by the user (e.g., ( p < 0.05 )) that is usually adjusted for multiple testing. Due to the fact multiple testing increases the probability of false positive discoveries, the significance threshold is generally adjusted using multiple comparison p-value adjustment methods, such as Bonferonni's. SNPs above this significant threshold line are considered significant.

#### 9.3. Linking significant loci to candidate genes (bonus)

To link significant loci to candidate genes:

Identify Significant SNPs: Extract the SNPs that surpass the significance threshold. **See 8.3. as a starting step if we decide to complete the bonus requirement**

Gene Annotation: Use gene annotation tools or databases (e.g., Ensembl, NCBI) to find genes located near these SNPs.

Biological Relevance: Investigate the biological functions of these genes to understand their potential role in the trait of interest.







