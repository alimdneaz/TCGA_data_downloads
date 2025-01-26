# Library
# 0. Setting Up The Environment
if(!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("EDASeq")
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("maftools")
install.packages("pheatmap")
BiocManager::install("SummarizedExperiment")
BiocManager::install("sesameData")
BiocManager::install("sesame")


# Importing Libraries
library(TCGAbiolinks) 
library(dplyr)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesameData)
library(sesame)

# get a list of projects
gdcprojects <- getGDCprojects()
View(gdcprojects)

# Project Summary
getProjectSummary("TCGA-BRCA")
# Project Query
query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling")

output_brca <- getResults(query_TCGA)
class(output_brca)
View(output_brca)
# dplyr: grammer of dataframe: dataframe = tibble
sum(table(output_brca$sample_type))

cancer.table <- output_brca %>% 
  select(sample_type) %>% 
  table()
cancer.table


# Main Query: 
query_braca2 <- GDCquery(project = "TCGA-BRCA",
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = 'STAR - Counts',
                         access = "open")
output_braca2 <- getResults(query_braca2)
View(output_braca2)
output_braca2 %>% 
  select(sample_type) %>% 
  table()
# Select 25 samples of each category
Normal.type <- output_braca2 %>% 
  filter(sample_type == "Solid Tissue Normal") %>% 
  select(cases) %>% 
  head(25)
Normal.type
Primary.tumor <- output_braca2 %>% 
  filter(sample_type == "Primary Tumor") %>% 
  select(cases) %>% 
  head(25)

# Now Do the query again
query_braca2 <- GDCquery(project = "TCGA-BRCA",
                         data.category  = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = 'STAR - Counts',
                         access = "open",
                         barcode = c("TCGA-AN-A04A-01A-21R-A035-13"))

# Download Only One Normal Data
GDCdownload(query_braca2)