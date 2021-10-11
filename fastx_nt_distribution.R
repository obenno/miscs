#! /usr/bin/Rscript


library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "input file, fastx_quality_stats output stats file."),
    make_option(c("-s", "--start"), type = "integer", default  = 1,
                help = "set start position of reads profile, [%default]"),
    make_option(c("-e", "--end"), type = "integer", default  = NULL,
                help = "set start position of reads profile, [%default]"),
    make_option("--plotN", action = "store_true", default = FALSE,
                help = "plot N nucleotides [%default]"),
    make_option(c("-o", "--output"), type = "character", default = "output.pdf",
                help = "output nucleotides distribution plot"),
    make_option("--width", type = "double", default = 10,
                help = "output plot width"),
    make_option("--height", type = "double", default = 6,
                help = "output plot height")
)

parser <- OptionParser(option_list=option_list,
                       description="Programm to plot nucleotide distribution from fastx_quality_stats result file")

opt <- parse_args(parser)

if(length(commandArgs(trailingOnly = TRUE))==0){
    print_help(parser)
    quit(save="no")
}
suppressPackageStartupMessages(library(tidyverse))

d <- read_tsv(opt$input, show_col_types = FALSE) %>%
    mutate(A_percent = A_Count/Max_count,
           T_percent = T_Count/Max_count,
           C_percent = C_Count/Max_count,
           G_percent = G_Count/Max_count,
           N_percent = N_Count/Max_count) %>%
    select(column,
           A_percent,T_percent,C_percent,G_percent,
           N_percent) %>%
    pivot_longer(-column,
                 names_to = "nucleotide",
                 values_to = "percent") %>%
    mutate(nucleotide = str_replace(nucleotide, "_percent", ""))

if(!is.null(opt$start)){
    d <- d %>%
        filter(column >= opt$start)
}

if(!is.null(opt$end)){
    d <- d %>%
        filter(column <= opt$end)
}

if(!opt$plotN){
    d <- d %>%
        filter(nucleotide != "N")
}

p <- ggplot(d, aes(x=column, y=percent, group = nucleotide,
                   color = nucleotide)) +
    geom_line(alpha = 0.8) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color="black"))

ggsave(opt$output,
       width = opt$width, height = opt$height)
