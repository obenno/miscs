#! /usr/bin/Rscript

library(optparse)

option_list <- list(
    make_option(c("-i", "--input"),
                help="input GTF file"),
    make_option(c("-g", "--genome"),
                help = "genome fasta file")
)

parser <- OptionParser(option_list = option_list,
                       description = "Extract exon and both intron sequences from transcript,\nlabel intron seq as lower case. Output to stdout.")

opt <- parse_args(parser)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome))


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

outSeqList <- list()

for(i in 1:length(exonSeqs)){
    strand <- strand(exons[[i]]) %>%
        as.character() %>% unique()
    outSeqList[i] <- pasteSeqs(exonSeqs[[i]], intronSeqs[[i]],
                           strand)
}

names(outSeqList) <- names(exons)

for(i in names(outSeqList)){
    cat(paste(i, outSeqList[[i]], sep="\t"), "\n")
}

