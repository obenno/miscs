#! /usr/bin/Rscript

suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default=NULL,
                help="input file [%default]"),
    make_option(c("-o", "--output"), type = "character", default=NULL,
                help="output file [%default]"),
    make_option(c("--noheader"), action="store_false", default=TRUE,
                help="input file has no header")
    )

opt <- parse_args(OptionParser(option_list=option_list))

a <- read.delim(opt$input, header = opt$noheader, sep="\t")
write.xlsx2(a, file=opt$output, sheetName="Sheet1", 
            col.names= opt$noheader, row.names=FALSE, append=FALSE)
