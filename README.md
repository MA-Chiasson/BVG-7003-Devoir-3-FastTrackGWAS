# BVG-7003-Assignment-3-FastTrackGWAS

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
Check the current working directory with getwd():
```
```r
getwd()
```
#### Conclusion
Defining a working directory streamlines file management, reduces errors, and helps organize your project efficiently.

### 4. Loading the dada
