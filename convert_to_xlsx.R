#! /usr/bin/Rscript

suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default=NULL,
                help="input file(s), could be comma seperated list [%default]"),
    make_option(c("-n", "--name"), type = "character", default="Sheet1",
                help="sheet names, could be comma seperated list [%default]"),
    make_option(c("-o", "--output"), type = "character", default=NULL,
                help="output file [%default]"),
    make_option(c("--header"), action="store_false", default=TRUE,
                help="read_tsv col_names option")
    )

opt <- parse_args(OptionParser(option_list=option_list))


inputFile <- unlist(strsplit(opt$input, ","))
sheetNames <- unlist(strsplit(opt$name, ","))
if(length(inputFile) != length(sheetNames)){
    stop("Please ensure sheet names length match number of inputs.")
}

d <- list()
for(i in 1:length(inputFile)){
    d[[sheetNames[i]]] <- read_tsv(inputFile[i],
                                   col_names = opt$header)
}

write_xlsx(d, path=opt$output,
           col_names= opt$header)
