# Skeletal-Muscle-Organoid-Transcriptomics
"Transcriptomic analysis of human iPSC-derived skeletal muscle organoids (GSE147513) identifying key regulatory hub genes through RNA-seq and PPI network analysis."
# 🧬 Skeletal Muscle Organoid Transcriptomics

## Project Overview
This repository contains the full bioinformatics pipeline for analyzing the maturation of human iPSC-derived skeletal muscle organoids (Dataset: **GSE147513**). By comparing 4-week vs 16-week samples, this study identifies the genetic drivers behind muscle fiber development.

## 🛠 Analysis Workflow
- **Data Source:** NCBI GEO (GSE147513)
- **Normalization:** Log2 transformation of FPKM values.
- **DGE Analysis:** Identified **5,508 DEGs** using the `limma` package in R.
- **Pathways:** Significant enrichment in **Muscle Contraction** and **Sarcomere Organization**.
- **Systems Biology:** Protein-Protein Interaction (PPI) network construction via STRING.

## 🚀 Key Findings
- **Upregulated Genes:** 3,004
- **Downregulated Genes:** 2,504
- **Identified Hub Genes:** TNNI2, CAV3, SRL, ACTA1, MYL1.

## 📂 Repository Contents
- `analysis.R`: The complete R script for the analysis.
- `/plots`: High-resolution Volcano, Dotplots, and Boxplots.
