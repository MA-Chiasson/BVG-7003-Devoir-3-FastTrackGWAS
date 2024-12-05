##############################################################
# Analyse_MAC
# Auteurs : Marc-Antoine Chiasson, Edouard Reed-Métayer, Mehdi Babaei, Ladan Ajdanian
# Date    : Sys.Date()
##############################################################

# 1. Introduction
# Cette analyse GWAS utilise le package rMVP pour identifier des SNPs
# associés à un trait phénotypique donné. Les données utilisées sont issues
# du fichier data.

# 2. Chargement des Bibliothèques
# Installer les bibliothèques si nécessaire
if (!requireNamespace("rMVP", quietly = TRUE)) install.packages("rMVP")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("mgsub", quietly = TRUE)) install.packages("mgsub")

# Charger les bibliothèques
library(rMVP)
library(ggplot2)
library(data.table)
library(dplyr)
library(mgsub)

# 3. Chargement des données
pheno_file <- "data/Phenotype_African.txt"
hmp_file <- "data/African_SNPs.hmp.txt"

# Charger et modifier le phénotype
pheno <- read.csv(pheno_file, sep = "\t", header = TRUE) %>% 
  mutate(Sample = mgsub(
    pattern = c("TGx", " "), 
    replacement = c("", ""), 
    string = Sample
  ))

# Sauvegarder le tableau phénotype modifié
write.table(pheno, "data/pheno_modified.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Conversion HapMap -> format MVP
MVP.Data(
  fileHMP = hmp_file,
  filePhe = "data/pheno_modified.txt",
  out = "mvp_hmp"
)

# 4. Exécution de l'Analyse GWAS
genotype <- attach.big.matrix("mvp_hmp.geno.desc")
phenotype <- read.table("mvp_hmp.phe", header = TRUE)
map <- read.table("mvp_hmp.geno.map", header = TRUE)

# Lancement du GWAS simple avec GLM
results <- MVP(
  phe = phenotype,
  geno = genotype,
  map = map,
  method = c("GLM"),
  nPC.GLM = 3,  # Nombre de composantes principales pour GLM
  threshold = 0.05
)

# Extraction des résultats et sauvegarde
map <- results$map
glm_results <- as.data.frame(results$glm.results)
colnames(glm_results) <- c("Effect", "SE", "P.value")

combined_results <- cbind(map, glm_results)
write.table(combined_results, "Test_GWAS_results.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Préparation des données pour la visualisation
vis_data <- combined_results[, c("SNP", "CHROM", "POS", "P.value")]
colnames(vis_data) <- c("SNP", "Chromosome", "Position", "P.value")

MVP.Report(
  vis_data,
  plot.type = c("qq", "manhattan"),
  threshold = 0.05 / nrow(vis_data) # Seuil de Bonferroni
)
