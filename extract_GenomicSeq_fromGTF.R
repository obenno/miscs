#! /usr/bin/Rscript

library(optparse)

option_list <- list(
    make_option(c("-i", "--input"),
                help="input GTF file"),
    make_option(c("-g", "--genome"),
                help = "genome fasta file"),
    make_option(c("--outfmt"), type="character", default = "fasta",
                help = "output result format, tabular or fasta"),
    make_option(c("-t", "--thread"), type="integer", default = 4,
                help = "threads number to be used")
)

parser <- OptionParser(option_list = option_list,
                       description = "Extract exon and both intron sequences from transcript,\nlabel intron seq as lower case. Output to stdout.")

opt <- parse_args(parser)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doFuture))

registerDoFuture()
plan(multisession, workers = opt$thread)

options(future.globals.maxSize= 2*1024*1024^2)

txdb <- makeTxDbFromGFF(opt$input)
genome <- readDNAStringSet(opt$genome)

## Note intron order is genomic coordinates based,
## not transcript based
## ?intronsByTranscript for detail
introns <- intronsByTranscript(txdb, use.names = T)
exons <- exonsBy(txdb, use.names = T, by="tx")

exonSeqs <-  getSeq(genome, exons)

intronSeqs <- getSeq(genome, introns)

pasteSeqs <- function(exonSeqs, intronSeqs, strand){
    if(strand == "-"){
        intronSeqs <- rev(intronSeqs)
    }
    intronSeqs <- stringr::str_to_lower(intronSeqs)
    intronSeqs <- c(intronSeqs, "")
    outSeq <- paste0(exonSeqs, intronSeqs, collapse = "")
    return(outSeq)
}


outSeqList <- foreach(i=1:length(exonSeqs)) %dopar% {
    strand <- GenomicRanges::strand(exons[[i]]) %>%
        as.character() %>% unique()
    pasteSeqs(exonSeqs[[i]],
              intronSeqs[[i]],
              strand)
}

names(outSeqList) <- names(exons)

format_seq <- function(seq, width=80){
    seq <- as.character(seq) %>%
        strsplit(split = "") %>%
        unlist
    outSeq <- foreach(i=seq(from=1, to=length(seq), by=width), .combine=c) %do% {
        seq[i:(i+width-1)] %>%
            na.omit %>%
            paste(collapse = "")
    }
    paste(outSeq, collapse = "\n")
}

if(opt$outfmt == "tabular"){
    for(i in names(outSeqList)){
        cat(paste(i, outSeqList[[i]], sep="\t"), "\n")
    }
}else if(opt$outfmt == "fasta"){
    for(i in names(outSeqList)){
        cat(paste0(">", i, "\n", format_seq(outSeqList[[i]])), "\n")
    }
}else{
    stop("Output format only support tabular or fasta.")
}
