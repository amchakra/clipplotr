# Script to install necessary libraries for CLIPplotR
# A. M. Chakrabarti

packages <- c("optparse", "BiocManager", "ggplot2", "scales", "data.table", "patchwork")
for(package in packages) {

    if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    message("Installing ", package)
    install.packages(package, character.only = TRUE, repos = "https://cloud.r-project.org")
    
  }

}

biocpackages <- c("GenomeInfoDb", "rtracklayer", "GenomicFeatures", "txdbmaker")
for(package in biocpackages) {
  
  if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    
    message("Installing ", package)
    BiocManager::install(package)
    
  }  
  
}