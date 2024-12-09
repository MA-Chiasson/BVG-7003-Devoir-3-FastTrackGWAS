# BVG-7003-Assignment-3-FastTrackGWAS

# Table of Contents
1. [Introduction](#introduction)
2. [Installation of R and RStudio](#installation-of-r-and-rstudio)
   1. [Installing R](#installing-r)
   2. [Installing RStudio](#installing-rstudio)
   3. [Launching RStudio](#launching-rstudio)
3. [Script Execution](#script-execution)
   1. [Input and Output File Formats](#input-and-output-file-formats)
   2. [Dependencies](#dependencies)
      1. [Installation](#installation)
      2. [Loading](#loading)
   3. [Define Working Directory](#define-working-directory)
   4. [Loading Data in R](#loading-data-in-r)
      1. [Loading Phenotypic Data](#loading-phenotypic-data)
      2. [Loading Genotypic Data](#loading-genotypic-data)

## Introduction

Genome-Wide Association Studies (GWAS) are a powerful tool for uncovering the genetic basis of traits, and their application in agricultural science is pivotal for crop improvement. This project focuses on analyzing the genetic basis of symbiotic nitrogen fixation in African soybean, using a publicly available dataset. The steps provided in this repository guide users through setting up the required tools and performing a GWAS analysis efficiently.

To enable reproducibility, this project employs **R** and **RStudio** as the core computational environment. R is renowned for its statistical computing capabilities and robust library ecosystem, making it a suitable choice for GWAS workflows. This repository demonstrates how to install and configure R and RStudio, manage dependencies, and execute the GWAS analysis pipeline.

The guide is designed to help users of all experience levels successfully replicate the analysis and adapt it to other datasets or traits of interest.

--- 

## Installation of R and RStudio

To perform a GWAS analysis using R, follow the steps below to install R and RStudio on your computer:

### 1. **Installing R**
R is the programming language required for the analysis:
1. Visit the official R website: [https://cran.r-project.org/](https://cran.r-project.org/).
2. Click on **Download R for [your operating system]** (Windows, macOS, or Linux).
3. Download the appropriate version for your system:
   - For Windows: Click on *base* and then download the installer.
   - For macOS: Download the version compatible with your system.
   - For Linux: Follow the instructions specific to your distribution.
4. Run the downloaded file and follow the installation instructions.

### 2. **Installing RStudio**
RStudio is a user-friendly interface for R:
1. Go to the RStudio website: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/).
2. Click on **RStudio Desktop** (free version).
3. Download the installer for your operating system.
4. Install RStudio by running the downloaded file.

### 3. **Launching RStudio**
1. Once installed, open RStudio.
2. Ensure that R is correctly configured with RStudio (RStudio automatically detects R if it is installed).

--- 

## Script execution
### 1. Input and Output File Formats
#### Input
#### Output

### 2. **Dependencies**
#### 2.1. Installation
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
Copier le code
install.packages(c("rMVP", "ggplot2", "data.table", "dplyr", "mgsub", "bigmemory"))
```
This would lead to the installation of packages every time the script is run, even if they are already installed, which is inefficient. The if condition ensures that packages are only installed when required, avoiding unnecessary installation.

#### 2.2. Loading
Once they've been installed, it's time to load them. To do this, use the library() command, inserting the package name between the brackets, as in the example below. 
```r
library(rMVP)
```
### 3. Define working directory
**Defining a working directory in R offers several key benefits:**

- File Management: It helps centralize project files (data, scripts, results) in one location, avoiding confusion.

- Simplified File Access: With a set working directory, you can use relative file paths (e.g., read.csv("data.csv")) instead of full paths.

- Easier Data Import/Export: Files can be easily imported or exported without needing to specify full file paths.

- Reduced Path Errors: Setting a working directory minimizes errors related to incorrect file paths, especially in team projects.

- Project Organization: It encourages better project structure, making it easier to manage and maintain files.

#### How to Set a Working Directory
Usually, the WD is the folder where our script and data are stored, and where we want to add the results of our analyses.
You can define a working directory using setwd():

```r
setwd("C:/Users/YourName/Documents/Project")
```
```r
Check the current working directory with getwd():
getwd()
```
#### Conclusion
Defining a working directory streamlines file management, reduces errors, and helps organize your project efficiently.

### 4. Loading Data in R
There are several good ways to load data into R, depending on the data format. Below, we explain how to load both phenotypic and genotypic data using specific commands.

#### 4.1. Loading Phenotypic Data
For phenotypic data, we typically use `read.csv()` or `read.table()` to import the data. The key difference is that `read.csv()` is used when the data is separated by commas, while `read.table()` is more general and allows specifying different delimiters (e.g., tab-separated data).

In our case, we use the following code to load the phenotypic data:

```r
pheno <- read.csv(pheno_file, sep = "\t", header = TRUE)
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

#### Conclusion
In summary:

- **Phenotypic data** is loaded using `read.csv()` (or `read.table()`) with parameters like `sep` (separator) and `header` (column names).
- **Genotypic data** is loaded using the `MVP.Data()` function, which handles different formats (VCF or HapMap) by specifying the appropriate file parameter (`fileVCF` or `fileHMP`). The output can be saved in the specified format using the `out` parameter.
