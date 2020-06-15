#! /usr/bin/Rscript

suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default=NULL,
                help="input file [%default]"),
    make_option(c("-o", "--output"), type = "character", default=NULL,
                help="output file [%default]")
    )

opt <- parse_args(OptionParser(option_list=option_list))

a <- read.delim(opt$input, header = T, sep="\t")
write.xlsx2(a, file=opt$output, sheetName="Sheet1", 
            col.names=TRUE, row.names=FALSE, append=FALSE)
