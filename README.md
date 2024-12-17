# BVG-7003-Assignment-3-FastTrackGWAS

## Table of Contents
- [1. Installation of R and RStudio](#1-installation-of-r-and-rstudio)
  * [1.1. Installing R](#11---installing-r--)
  * [1.2. Installing RStudio](#12---installing-rstudio--)
  * [1.3. Launching RStudio](#13---launching-rstudio--)
- [2. Data Download and Preparation](#2---data-download-and-preparation--)
  * [2.1. Test Data Used](#21-test-data-used)
  * [2.2. Other data recommandations](#22-other-data-recommandations)
- [3. Dependencies](#3---dependencies--)
  * [3.1. Installation of the packages](#31-installation-of-the-packages)
  * [3.2. Loading the packages](#32-loading-the-packages)
  * [4. Define working directory](#4-define-working-directory)
    + [How to Set a Working Directory](#how-to-set-a-working-directory)
  * [5. Loading data in R](#5-loading-data-in-r)
    + [5.1. Loading genotype and phenotype data](#51-loading-genotype-and-phenotype-data)
  * [6. Data Preprocessing](#6-data-preprocessing)
    + [6.1. Useful tools for inspecting data structure](#61-useful-tools-for-inspecting-data-structure)
      - [6.1.1. `str()`](#611--str---)
      - [6.1.2. `head()`](#612--head---)
      - [6.1.3. Example of formatting](#613-example-of-formatting)
    + [6.2. Quality Control](#62-quality-control)
      - [6.2.1. Genotype data quality control](#621-genotype-data-quality-control)
      - [6.2.2. Phenotype data quality control](#622-phenotype-data-quality-control)
    + [6.3. Generating Covariates](#63-generating-covariates)
  * [7. Execution of the GWAS analysis](#7-execution-of-the-gwas-analysis)
    + [7.1. Loading genotype and phenotype data in the rMVP environment](#71-loading-genotype-and-phenotype-data-in-the-rmvp-environment)
      - [7.1.1. Generating rMVP readable files, VCF data version:](#711-generating-rmvp-readable-files--vcf-data-version-)
      - [7.1.2. Generating rMVP readable files, HapMap data version:](#712-generating-rmvp-readable-files--hapmap-data-version-)
      - [7.1.3. Loading the data in rMVP environment:](#713-loading-the-data-in-rmvp-environment-)
    + [7.2. Loading covariates in MVP function (optional)](#72-loading-covariates-in-mvp-function--optional-)
    + [7.3. Running MVP function](#73-running-mvp-function)
  * [8. Example of results obtainable from MVP and their interpretation](#8-example-of-results-obtainable-from-mvp-and-their-interpretation)
    + [8.1. Q-Q plots](#81-q-q-plots)
    + [8.2. Manhattan plots (association strength between individual loci and the studied phenotype)](#82-manhattan-plots--association-strength-between-individual-loci-and-the-studied-phenotype-)
    + [8.3. Filtering significant SNPs from the analysis](#83-filtering-significant-snps-from-the-analysis)
  * [9. Expected results from this pipeline with the example data](#9-expected-results-from-this-pipeline-with-the-example-data)
    + [9.1. Interpreting QQ plots](#91-interpreting-qq-plots)
    + [9.2. Identifying significant SNPs from Manhattan plots](#92-identifying-significant-snps-from-manhattan-plots)
    + [9.3. Linking significant loci to candidate genes (bonus)](#93-linking-significant-loci-to-candidate-genes--bonus-)

<small><i><a href='http://github.com/3kh0/readme-toc/'>Table of contents generated with readme-toc</a></i></small>


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

### 1.1. Installing R
R is the programming language required for the analysis:
1. Visit the official R website: [https://cran.r-project.org/](https://cran.r-project.org/).
2. Click on **Download R for [your operating system]** (Windows, macOS, or Linux).
3. Download the appropriate version for your system:
   - For Windows: Click on *base* and then download the installer.
   - For macOS: Download the version compatible with your system.
   - For Linux: Follow the instructions specific to your distribution.
4. Run the downloaded file and follow the installation instructions.

### 1.2. Installing RStudio
RStudio is a user-friendly interface for R:
1. Go to the RStudio website: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/).
2. Click on **RStudio Desktop** (free version).
3. Download the installer for your operating system.
4. Install RStudio by running the downloaded file.

### 1.3. Launching RStudio
1. Once installed, open RStudio.
2. Ensure that R is correctly configured with RStudio (RStudio automatically detects R if it is installed).
---
## 2. Data Download and Preparation
This analysis is based on the dataset provided in the paper:
**"Genome-Wide Association Analyses Reveal the Genetic Basis of Symbiotic Nitrogen Fixation in African Soybean"**
- **DOI:** [10.1007/s00122-019-03499-7](https://link.springer.com/article/10.1007/s00122-019-03499-7)
- **Data link:** [Figshare Dataset](https://figshare.com/projects/GenomeWide_Association_Analyses_Reveal_the_Genetic_Basis_of_Symbiotic_Nitrogen_Fixation_in_African_Soybean/65885)
### 2.1. Test Data Used

For the test data in this analysis, the following were selected:

- **Genotypic data:** Information for all chromosomes.


  File `African_SNPs.hmp.txt` , of HapMap format, can be found in the data directory of this repository.

  Note that genotype data of VCF format could also be used with the example script if this was your case. 
  
- **Phenotypic data:** The trait SDWi/SDWf (SDW ratio) represents shoot dry weight in inoculated (SDWi) and fertilized (SDWf).

  File `Phenotype_African.txt`, in .txt format, can be found in the data directory of this repository.

  Note that phenotype data could be provided in other table formats, as long as the first column represents the genetic ID, and the other columns represent trait values. See [https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#phenotype]([url](https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#phenotype)) for more detailed information.

### 2.2. Other data recommendations

**Missing genotype data imputation** :
We recommend the users impute missing data **before** importing their data into the pipeline. Data imputation methods relying on inference from haplotypes are very accurate and should preferentially be used. Unfortunately, this can not be done through the rMVP package, but it can be done using external software such as LD-kNNi, FILLIN, FSHAP, or BEAGLE. 

--- 

## 3. Dependencies
### 3.1. Installation of the packages
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
  

Dependencies for these packages will be automatically installed when executing the script, ensuring a smooth setup process.
For this purpose, the installation of packages is conditional using the following function:

```r
if (!requireNamespace("package", quietly = TRUE)) install.packages("package")
```

This checks whether each package is already installed before attempting to install it. If the package is not found, it is installed. This approach ensures that packages are only installed when necessary, preventing redundant installations during each script iteration, which can slow down execution.

On the other hand, if a single line of code was used to install all the packages at once, like this:

```r
install.packages(c("rMVP", "ggplot2", "data.table", "dplyr", "purrr", "mgsub", "bigmemory", "vcfR", "outliers"))
```
This would lead to installing packages every time the script is run, even if they are already installed, which is inefficient. The if condition ensures that packages are only installed when required, avoiding unnecessary installation.

### 3.2. Loading the packages
Once they've been installed, it's time to load them. To do this, use the library() command, inserting the package name between the brackets, as in the example below. 
```r
library(rMVP, ggplot2, data.table, dplyr, purrr, mgsub, bigmemory, vcfR, outliers)
```
---
### 4. Define working directory
**Defining a working directory in R offers several key benefits:**

- **File Management**: It helps centralize project files (data, scripts, results) in one location, avoiding confusion.

- **Simplified File Access**: With a set working directory, you can use relative file paths (e.g., read.csv("data.csv")) instead of full paths.

- **Easier Data Import/Export**: Files can be easily imported or exported without specifying full file paths.

- **Reduced Path Errors**: Setting a working directory minimizes errors related to incorrect file paths, especially in team projects.

- **Project Organization**: It encourages better project structure, making managing and maintaining files easier.

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

#### 6.1. Useful tools for inspecting data structure
**IMPORTANT NOTE:** Before starting, it is important to note that genotype data files are sometimes complex to modify directly in R, as they have special file formats. Therefore, it is sometimes preferable that any necessary modifications to harmonize the datasets be made to the phenotype data file, as we will see below.

To check if the data is properly formatted, we can use two commonly used commands in R: `str()` and `head()`. These functions are useful for quickly exploring the structure and content of our data before proceeding with more complex analyses.

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

The final quality control step implemented in this pipeline is that of missing data imputation. **Show an example or remove this sentence**

##### 6.2.1. Genotype data quality control

1. **Filtering loci based on their Minor Allele Frequency (MAF).** We suggest removing loci with a MAF lower than 0.05 to ensure that only common variants are included in the analysis, but other tresholds may be used. 

Loci with low MAF (for example, lower than 0.05) carry an elevated risk of having genotyping errors (i.e. the minor allele may be a genotyping error and not even exist in reality). Further, loci with alleles present at very low frequencies are more prone to statistical errors, whether false positives or false negatives. Therefore, including them doesn't contribute significantly to the discovery of meaningful associations.

2. **Filtering loci with an important proportion of missing data.**  We suggest removing loci with 10% or more missing data, but other thresholds may be used.

A high level of missing data suggests lower genotyping quality for that loci and, even if missing data is imputed, lowers accuracy to some degree. Loci with too much missing data can additionally induce statistical errors.



Using `filter_hapmap()` or `filter_vcfmap()` functions, we can filter SNPs with too low MAF, with a MAF threshold modifiable, but set by default at 5%. At the same time, these functions filter out monomorphic loci that would not be useful for a GWAS analysis. The same functions allow for filtering SNPs with too much missing data, with the threshold modifiable, but set at 10% as default option. Below is the R code used for this process:

For Hapmap data
```r
filtered_hmp <- filter_hapmap(geno, freq_threshold = 5, na_threshold = 10)
```
Or for VCF data
```r
filtered_vcf <- filter_vcfmap(geno,  freq_threshold = 5, na_threshold = 10)
```
In both cases, thresholds can be either left out to use the default options, or modified with the threshold expressed in percentage.

##### 6.2.2. Phenotype data quality control 

1. **Verifying or removing outliers.** There is some debate concerning whether or not outlier observations should be removed from the analysis. Outliers can skew the results and lead to false associations in GWAS. By removing them, we ensure that the data accurately represents the population and improves the reliability of the findings. However, if you are dealing with your own experimental dataset, we suggest verifying if these outliers are true measures or result from measuring errors for example, to help you take a decision about removing them or not. For the sake of this example, and in the absence of real experimental information, we filter outlier data based on quantile distribution. 

2. **Verifying normality of distribution.** Many statistical tests used in GWAS assume that the data follows a normal distribution. Verifying and transforming the data to meet this assumption helps in obtaining valid and interpretable results. This step is not implemented in the pipeline currently as verifying normality shoudl involve proper scientific judgement, especially with large datasets where many formal normality tests can be very sensitive. Further, when original data is not normally distributed, there are many possibilities regarding the transformations that can be applied and which need to consider experimental design. For example, one could wish to transform the phenotype data using logarithmic, square root or inverse transformation, while the data may also be suited to add specific covariates and/or the transformation of these covariates. This step is therefore intentionally left out to the scientific judgement of the user.

3. **Handling missing data.** Missing data can reduce the power of the study and introduce bias. Properly addressing missing data through imputation or exclusion ensures that the analysis is robust and the results are not compromised. Here we will not discuss imputation of missing phenotypic data, as there are many different methods, which are outside the scope of this pipeline. We must however make sure the analysis doesn't get blocked by missing data, and we can even apply more severe quality tresholds by eliminating entries with too much missing phenotype data. Indeed, an individual or a variety showing many missing values may have been assessed under inappropriate experimental conditions, which leads to question its entire phenotype data. Under circumstances where detailed experimental information is not known, it is cautious to remove such individuals or varieties.

The `process_phenotypic_data()` is designed to remove outliers and to remove entries with too much missing data. Outliers are removed based on quantile distribution, and a treshold of 10% missing phenotype data for an entry gets it removed from the dataset. Note that individual missing data points can be handled directly by the `MVP()` function, so there is no need to preprocess those.

Example usage of the `process_phenotypic_data()` function is shown below.

```r
filtered_pheno <- process_phenotypic_data(pheno, na_threshold = 10)
```
***
#### 6.3. Generating Covariates

To account for population structure or relatedness in the data, you may generate covariates such as Principal Component Analysis (PCA) scores or a relatedness (kinship) matrix. These covariates are used to adjust for the potential confounding factors in the GWAS analysis that are kinship or population structure. Indeed, population structure and relatedness may induce non-random distribution of alleles in the sampling pool. GWAS analysis over phenotypes that happen to covary with the population structure would then result in non-relevant associations with these non-randomly distributed alleles. Including a kinship matrix or a a PCA can help reduce this risk.

The recommended option if you don't already have your own kinship matrix or PCA of population structure, is to directly do it in the GWAS analysis part which is done when using the function MVP(), through the arguments `K`and `nPC` (see section 7 below for complete analysis). Closing the argument `K`with a # will make the function generate it automatically. Setting the `nPC`arguments (depending on which method(s) you are using) with the number of principal components you would like to use for each will compute them automatically.

``` r
MVP(
    #K=kinship,        
    nPC.GLM=5,   
    nPC.MLM=3,             
    nPC.FarmCPU=3
    )
```
Alternatively, kinship matrices and population structure PCAs can be generated at the data conversion step, with the function `MVP.Data()`, by setting the arguments `fileKin` or `filePC` as `TRUE`. This is not covered in the sample script.

``` r
MVP.Data(
         fileKin=TRUE,
         filePC=TRUE,
         )
```

There is also the possibility of generating a kinship matrix from an already rMVP imported genotype file using the `MVP.Data.Kin()` function. This is not covered in the sample script.

``` r
MVP.Data.Kin(
             fileKin = TRUE,
             mvp_prefix = "mvp",
             out = NULL
             )
```

And even the possibility of generating a PCA for population structure, also from from an already rMVP imported genotype file, using the `MVP.Data.PC()` function. This is not covered in the sample script.

``` r
MVP.Data.PC(
            filePC = TRUE,
            mvp_prefix = "mvp",
            out = NULL,
            pcs.keep = 5, # this number representing the number of PCs to keep
            )
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

Although not covered in the example code, it is possible to also include other environmental variables or covariates as fixed effects to account for various experimental designs (for example, breed, sex, weight, age, diet, socioeconomic status, batch effects, genotyping platform, etc.). The inclusion of such factors is done through the arguments CV.GLM, CV.MLM, CV.FarmCPU. See [https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#advanced]([url](https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#advanced)), as the original creators of this package, for more details.

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

### 8. Example of results obtainable from MVP and their interpretation

`MVP()` function automatically outputs visuals from the GWAS analysis if you add the option `"plot"` to the vector of `file.output` argument (see above). Plots can also be viewed separately and customized with many options, see rMVP original github for more details [https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#output]([url](https://github.com/xiaolei-lab/rMVP?tab=readme-ov-file#output)).

#### 8.1. Q-Q plots 

Q-Q plots are generally produced to help assess the quality of the GWAS analysis. They can help identify if there are deviations from the expected distribution, which can indicate potential issues such as population stratification, genotyping errors, or true genetic associations.

![Sample_QQplots_combined](https://github.com/user-attachments/assets/150fcbe7-7740-471e-a6b7-e6c258a64e24)

If the points on the Q-Q plot follow the diagonal line, it suggests that the observed p-values match the expected distribution, indicating no systematic bias.
Deviations above the line, especially at the tail, suggest an excess of small p-values, which could indicate true associations or potential confounding factors.
So for instance, in the above sample Q-Q plots, the middle Q-Q plot would indicate that the method used may be the best suited out of the three shown. 

#### 8.2. Manhattan plots (association strength between individual loci and the studied phenotype)

Manhattan plots display the p-values of SNPs across the genome:

![Sample_ManhattanPlot](https://github.com/user-attachments/assets/5f911087-1cab-4ac9-9057-a57990ddf4c5)

Peaks: Tall peaks indicate SNPs with low p-values, suggesting significant associations.

Threshold Line: The horizontal line represents the significance threshold set by the user (e.g., ( p < 0.05 )) that is usually adjusted for multiple testing. Due to the fact multiple testing increases the probability of false positive discoveries, the significance threshold is generally adjusted using multiple comparison p-value adjustment methods, such as Bonferonni's. SNPs above this significant threshold line are considered significant.

#### 8.3. Filtering significant SNPs from the analysis
**Does the presence of this requirement in the assignment instructions mean we need to filter the genotype file as a final step of the pipeline ?**

---

### 9. Expected results from this pipeline with the example data

#### 9.1. Interpreting Q-Q plots

**Add actual QQplots from our pipeline + our example datasets**

#### 9.2. Identifying significant SNPs from Manhattan plots

**Add actual Manhattan plots from our pipeline + our example datasets**

#### 9.3. Linking significant loci to candidate genes (bonus)

To link significant loci to candidate genes:

Identify Significant SNPs: Extract the SNPs that surpass the significance threshold. **See 8.3. as a starting step if we decide to complete the bonus requirement**

Gene Annotation: Use gene annotation tools or databases (e.g., Ensembl, NCBI) to find genes located near these SNPs.

Biological Relevance: Investigate the biological functions of these genes to understand their potential role in the trait of interest.







