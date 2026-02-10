#!/usr/bin/env Rscript

# argument information
# 1 - Sample/genotype quality control markdown file

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# markdown template
markdown_rmd <- args[1]

# sample directory from nextflow module
sampleDir <- args[2]

rmarkdown::render(markdown_rmd, params = list(sampleDir = sampleDir))