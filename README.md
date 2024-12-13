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

![Explainer-close-up](https://github.com/user-attachments/assets/cd9c6e19-db32-427a-a9ea-43a81169bb9a)

We used GWAS to investigate the genetic foundation of traits associated with symbiotic nitrogen fixation (SNF) in African soybeans, identifying potential candidate genes. This research offers crucial insights into SNF-related traits and will expedite progress in SNF breeding initiatives.
The steps provided in this repository guide users through setting up the required tools and efficiently performing a GWAS analysis.

This project employs **R** and **RStudio** as the core computational environment to enable reproducibility. R is renowned for its statistical computing capabilities and robust library ecosystem, making it suitable for GWAS workflows. This repository demonstrates how to install and configure R and RStudio, manage dependencies, and execute the GWAS analysis pipeline.

The guide is designed to help users of all experience levels successfully replicate the analysis and adapt it to other datasets or traits of interest.

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

--- 

## 2. Script execution
### 2.1. Input and Output File Formats
#### 2.1.1. Input
#### 2.1.2. Output
**REPLACE BY THE RESULT THE USER SHOULD HAVE**

### 2.2. **Dependencies**
#### 2.2.1. Installation of the packages
Open the script and start to install the following R packages which are required to execute this analysis:
- **rMVP**: A memory-efficient, visualization-enhanced, and parallel-accelerated tool for GWAS.
- **ggplot2**: For creating high-quality graphics.
- **data.table**: For fast and flexible data manipulation.
- **dplyr**: For data manipulation with a consistent grammar.
- **mgsub**: For multiple, simultaneous string substitutions.

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

#### 2.2.2. Loading the packages
Once they've been installed, it's time to load them. To do this, use the library() command, inserting the package name between the brackets, as in the example below. 
```r
library(rMVP)
```
---
### 3. Define working directory
**Defining a working directory in R offers several key benefits:**

- File Management: It helps centralize project files (data, scripts, results) in one location, avoiding confusion.

- **Simplified File Access**: With a set working directory, you can use relative file paths (e.g., read.csv("data.csv")) instead of full paths.

- **Easier Data Import/Export**: Files can be easily imported or exported without needing to specify full file paths.

- **Reduced Path Errors**: Setting a working directory minimizes errors related to incorrect file paths, especially in team projects.

- **Project Organization**: It encourages better project structure, making it easier to manage and maintain files.

#### How to Set a Working Directory
Usually, the WD is the folder where our script and data are stored, and where we want to add the results of our analyses.
You can define a working directory using setwd():

```r
setwd("C:/Users/YourName/Documents/Project")
```
Check the current working directory with getwd():
```r
getwd()
```
**Do we want to do something more like Adrian exemple?** 
My take is we don't need to, and it is more steps where the user has to change info in the script. I'd leave it as you did up here ! -Edouard

#### Conclusion
Defining a working directory streamlines file management reduces errors, and helps organize your project efficiently.

---

### 4. Loading data in R
There are several good ways to load data into R, depending on the data format. Below, we explain how to load both phenotypic and genotypic data in the R environment using specific commands.

#### 4.1. Loading phenotypic data
For phenotypic data, we typically use `read.csv()` or `read.table()` to import the data. The key difference is that `read.csv()` is used when the data is separated by commas, while `read.table()` is more general and allows specifying different delimiters (e.g., tab-separated data).

In our case, we use the following code to load the phenotypic data:
```r
pheno <- read.table(pheno_file, sep = "\t", header = TRUE)
```
Explanation of the parameters:
- `pheno_file`: The file path to the phenotypic data.
- `sep = "\t"`: This specifies that the data is tab-separated.
- `header = TRUE`: This indicates that the first row contains column names.

#### 4.2. Loading Genotypic Data
For genotypic data, we use the `MVP.Data()` function from the rMVP package. There are two types of genotypic data commonly used: HapMap and VCF. We load each type of data with different parameters.

1. **Loading VCF Data**: If the genotypic data is in VCF (Variant Call Format), we use the following command:
```r
MVP.Data(fileVCF = vcf_file, filePhe = pheno_file, out = "mvp_vcf")
```
- `fileVCF`: The file path to the VCF file.
- `filePhe`: The file path to the phenotypic data.
- `out` = "mvp_vcf": The output file name (the data will be saved in this format).

2. **Loading HapMap Data**: If the genotypic data is in HapMap format, we use the following command:
``` r
MVP.Data(fileHMP = hmp_file, filePhe = pheno_file, out = "mvp_hmp")
```
- `fileHMP`: The file path to the HapMap file.
- `filePhe`: The file path to the phenotypic data.
- `out` = "mvp_hmp": The output file name.
  
**However, in many cases the MVP function won't work perfectly the first time, as you need to ensure that the data is formatted in the right way (see section 5).**

#### Conclusion
- **Phenotypic data** is loaded using `read.csv()` (or `read.table()`) with parameters like `sep` (separator) and `header` (column names).
- **Genotypic data** is loaded using the `MVP.Data()` function, which handles different formats (VCF or HapMap) by specifying the appropriate file parameter (`fileVCF` or `fileHMP`). The output can be saved in the specified format using the `out` parameter.

---

### 5. Data Preprocessing
In this section, we guide users through the following steps for data preprocessing. These steps are crucial for ensuring the quality and integrity of the data used in the analysis. Preprocessing helps to clean, format, and structure the data properly, minimizing errors and biases in the final results.

#### 5.1. Formatting Genotype and Phenotype Files
**IMPORTANT NOTE:** Before starting, it is important to note that the genotype data file cannot be modified directly in R, as it is a special file format. Therefore, any necessary modifications should be made to the phenotype data file, as we will see below.

To check if the format of our data is correct, we can use two commonly used commands in R: `str()` and `head()`. These functions are useful for quickly exploring the structure and content of our data before proceeding with more complex analyses.

##### 5.1.1. `str()`
The str() (structure) function is used to display the internal structure of an object in R. It provides a quick overview of the data structure, including the type of each column (e.g., numeric, factor, etc.), the number of missing values, and a sample of the first few values in each column.

**Usefulness of `str()`:**

It checks the data types in each column (e.g., whether it's a factor, character, numeric value).
It helps detect any formatting errors, such as a numeric column being incorrectly read as a factor.
It provides a summary of the data object, which is useful for verifying everything is in order before applying transformations or analyses.
Example:

```r
str(pheno)
```
This allows you to check that the phenotype data is in the correct format before using it in the analysis.

##### 5.1.2. `head()`
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
These two commands are essential tools for quickly inspecting the structure and content of your data. They help identify potential errors early in the analysis process, saving a lot of time and effort upfront. Using str() and head() before starting any detailed analysis ensures that the data is properly formatted and ready to be processed.

##### 5.1.3. Example of formatting
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


#### 5.2. Quality Control

The quality control step involves filtering SNPs (Single Nucleotide Polymorphisms) based on their Minor Allele Frequency (MAF). We suggest removing SNPs with a MAF lower than 0.05 to ensure that only common variations are included in the analysis, but other tresholds may be used. 

Loci with low MAF (for example, lower than 0.05) carry an elevated risk of having genotyping errors (i.e. the minor allele may be a genotyping error and not even exist in reality). Further, loci with alleles present at very low frequencies are more prone to statistical errors, whether false positives or false negatives. Therefore, including them doesn't contribute significantly to the discovery of meaningful associations.

A second quality control step is also to eliminate loci with an important proportion of missing data. We suggest removing loci with 10% or more missing data, but other thresholds may be used.

Loci with too much missing data can potentially induce statistical errors.

The final quality control step implemented in this pipeline is that of missing data imputation. We recommend the users to impute missing data before importing their data into the pipeline, for example using LD-kNNi, FILLIN, FSHAP or BEAGLE software.  ***MAYBE THIS SHOULD BE MOVED TO THE DATA PREPROCESSING PART***

##### 5.2.1. Filtering out SNPs with MAF < 0.05

We calculate the MAF for each SNP and filter out those with a MAF below 0.05. At the same time, this is filtering out monomorphic loci that would not be useful for a GWAS analysis. Below is the R code used for this process:

***THIS PART OF THE CODE DOESN'T WORK FOR ME, I have created a new function to do this MAF filtering (as well as data completeness filtering) for hapmap format, I think it is better to do this filtering step before moving on to the MVP format which is quite complex.***
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

#### 5.3. Generating Covariates

To account for population structure or relatedness in the data, you may generate covariates such as Principal Component Analysis (PCA) scores or a relatedness matrix. These covariates are used to adjust for potential confounding factors in the GWAS analysis such as kinship or population structure.

```

This section explains how we preprocess the data, including formatting the files, performing quality control, and generating covariates for the analysis.

```
---

### 6. Execution of the GWAS analysis

#### 6.1. Loading genotype and phenotype data in the rMVP environment 
The first step of running the GWAS analysis is to load properly the genotypic and phenotypic data

Once you have correctly imported your data in the R environment (see step 4) and preprocessed it (see step 5), you are ready to load it in the rMVP environment.
The sample code below shows a simple way to load hmp format genotypic data, using the attach.big.matrix() function.
If you have a genotype map, you can also add it here.

```
genotype <- attach.big.matrix("mvp_hmp.geno.desc")
phenotype <- read.table("mvp_hmp.phe", header = TRUE)
map <- read.table("mvp_hmp.geno.map", header = TRUE)

```
#### 6.2. Loading covariates in rMVP function (optional)

#### 6.3. Running rMVP function


---

### 7. Example of results obtained from rMVP

#### 7.1. QQ plots 

#### 7.2. Manhattan plots (association strength between individual loci and the studied phenotype)

#### 7.3. Filtering significant SNPs from the analysis

---

### 8. Interpretation of results

#### 8.1. Interpreting QQ plots

#### 8.2. Identifying significant SNPs from Manhattan plots

Manhattan plots display the p-values of SNPs across the genome:

Peaks: Tall peaks indicate SNPs with low p-values, suggesting significant associations.

Threshold Line: The horizontal line represents the significance threshold set by the user (e.g., ( p < 0.000000005 )). Due to multiple testing increasing the probability of false positive discoveries, the significance threshold is generally adjusted using multiple comparison p-value adjustment methods, such as Bonferonni's. SNPs above this significant threshold line are considered significant.

#### 8.3. Linking significant loci to candidate genes (bonus)

To link significant loci to candidate genes:

Identify Significant SNPs: Extract the SNPs that surpass the significance threshold.

Gene Annotation: Use gene annotation tools or databases (e.g., Ensembl, NCBI) to find genes located near these SNPs.

Biological Relevance: Investigate the biological functions of these genes to understand their potential role in the trait of interest.







