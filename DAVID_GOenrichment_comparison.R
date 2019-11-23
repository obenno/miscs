#! /usr/bin/Rscript

## This script is to extract GO enrichment
## information from multiple results of DAVID
## and draw clusterProfiler like dot plot
suppressPackageStartupMessages(library(tidyverse))
library(optparse)

## make options
##ars <- commandArgs(trailingOnly = T)
##option.list <- list(
##    make_option(c("-i", "--input"), type="character", default=NULL,
##                help="input Gene list files, separate by comma"),
##    make_option(c("-n", "--names"), type="character", default=NULL,
##                help="category names, separatre by comma"),
##    make_option(c("-o", "--output"), type="character", default=NULL,
##                help="output file name")
##)
##
##usage <- "usage: %prog -i A.lst,B.lst -n TypeA,TypeB -o output.pdf\n"
##
##parser <- OptionParser(option_list=option.list,
##                       description="Comparative GO enrichment",
##                       usage=usage)
##opt <- parse_args(parser, args=args, positional_arguments=T)

args <- commandArgs(trailingOnly = T)
option.list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input Gene list files, separate by comma"),
    make_option(c("-n", "--names"), type="character", default=NULL,
                help="category names, separatre by comma"),
    make_option(c("-f","--flip"), action="store_true", default=FALSE,
                help="Whether flip the plot [%default]"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output file name"),
    make_option(c("-k","--KEGG"), action="store_true", default=FALSE,
                help="Whether input is KEGG enrichment"),
    make_option(c("-w", "--width"), type="double", default=9,
                help="output file width"),
    make_option(c("-e", "--height"), type="double", default=6,
                help="output file height")
)
usage <- "usage: %prog -i A.lst,B.lst -n TypeA,TypeB -o output.pdf\n"
parser <- OptionParser(option_list=option.list,
                       description="Comparative GO enrichment",
                       usage=usage)
opt <- parse_args(parser, args=args, positional_arguments=T)
inputFiles <- unlist(strsplit(opt$options$input, ",",
                              fixed = TRUE))

CatNames <- unlist(strsplit(opt$options$names, ",",
                            fixed = TRUE))

message("reading input arguments...")
if(length(inputFiles)!=length(CatNames)){
    stop("The category names doesn't match input list length...")
}

res <- list()

for(i in 1:length(inputFiles)){
    res[[i]] <- read_tsv(inputFiles[i]) %>%
        mutate(Description=case_when(opt$options$KEGG==FALSE~str_sub(Term, 12),
                                     opt$options$KEGG==TRUE~str_sub(Term, 10)),
               GO_id=case_when(opt$options$KEGG==FALSE~str_sub(Term, start=1, 10),
                               opt$options$KEGG==TRUE~str_sub(Term, start=1, 8)),
               sample=CatNames[i]) %>%
        ##select(sample, GO_id, Description, Count, FDR)
        select(sample, GO_id, Description, Count, PValue) # Use P-value instead of FDR
}

##head(res)
combinedList <- res %>% bind_rows() %>%
    group_by(sample) %>%
    arrange(desc(Count), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Description=factor(Description,
                              levels = rev(unique(Description)),
                              ordered = TRUE)) %>%
    mutate(sample=factor(sample,
                         levels = CatNames))
##combinedList %>% select(sample, Description, Count) %>%
##    pivot_wider(names_from = sample, values_from = Count) %>%
##    arrange()
## Filter List with count and FDR
combinedList <- combinedList %>% filter(PValue<=0.01, Count>=2)
##head(combinedList)
p <- ggplot(data=combinedList, aes(y=Description, x=sample))
p <- p+geom_point(shape=21, color="grey20", stroke = 0.8,
                  aes(size=Count, fill=PValue))
#p <- p+facet_grid(sample~., scales="free", space="free")
p <- p+scale_fill_gradient(name="Pvalue",
                            high = "#700CBC", low = "#EA202C")
p <- p+theme_bw()+xlab("")+ylab("")

ggsave(filename=opt$options$output, width = opt$options$width, height = opt$options$height)

