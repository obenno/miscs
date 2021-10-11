#! /usr/bin/env Rscript


suppressPackageStartupMessages(library(optparse))

args <- commandArgs(trailingOnly = T)
option.list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input agri GO files, separate by comma"),
    make_option(c("-n", "--names"), type="character", default=NULL,
                help="category names, separatre by comma"),
    make_option(c("-c", "--category"), type="character", default=NULL,
                help="GO category, BP, MF, CC or All"),
    ## deprecate simplify
    ##make_option(c("-s", "--simplify"), action="store_true", default=FALSE,
    ##            help="If perform simplify [%default]"),
    make_option(c("-f","--flip"), action="store_true", default=FALSE,
                help="Whether flip the plot [%default]"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output file name"),
    make_option(c("-w", "--width"), type="double", default=9,
                help="output file width"),
    make_option(c("-e", "--height"), type="double", default=6,
                help="output file height"),
    make_option(c("-y", "--legendY"), type="double", default=-4,
                help="Legend y asix position"),
    make_option(c("-x", "--legendX"), type="double", default=0.8,
                help="Legend X asix position")
)
usage <- "usage: %prog -i A_GO.txt,B_GO.txt -n TypeA,TypeB -o output.pdf\n"
parser <- OptionParser(option_list=option.list,
                       description="Comparative GO enrichment",
                       usage=usage)
opt <- parse_args(parser, args=args, positional_arguments=T)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))

set.seed(123)

inputFiles <- unlist(strsplit(opt$options$input, ",",
                       fixed = TRUE))
CatNames <- unlist(strsplit(opt$options$names, ",",
                     fixed = TRUE))

if(length(inputFiles)!=length(CatNames)){
    stop("The category names doesn't match input list length...")
}

if(is.null(opt$options$category)){
    message("GO category not indicated, use all three")
    opt$options$category <- "All"
}else if(!(opt$options$category %in% c("BP", "MF", "CC", "All"))){
    stop("Please indicate GO category")
}

GO_input <- list()
for(i in 1:length(inputFiles)){
    GO_input[[i]] <- read_tsv(inputFiles[i]) %>%
        mutate(sample=CatNames[i])
}

GO_result <- bind_rows(GO_input) %>%
    filter(FDR<=0.05) %>%
    select(sample, GO_acc, Term, term_type, FDR, queryitem) %>%
    mutate(FDR = -log10(FDR))

if(opt$options$category == "BP"){
    GO_result <- GO_result %>%
        filter(term_type == "P")
}else if(opt$options$category == "MF"){
    GO_result <- GO_result %>%
        filter(term_type == "F")
}else if(opt$options$category == "CC"){
    GO_result <- GO_result %>%
        filter(term_type == "C")
}

GO_result <- GO_result %>%
    arrange(sample, queryitem) %>%
    mutate(Term = factor(Term, levels = unique(Term)))

Term_List <- list()
for(i in 1:length(CatNames)){
    Term_List[[i]] <- GO_result %>%
        filter(sample == CatNames[i]) %>%
        pull(Term)
}
names(Term_List) <- CatNames

##Term_List
m1=make_comb_mat(Term_List)
term_order <- vector()
for(i in comb_name(m1)){
    term_order <- c(term_order,extract_comb(m1, i))
}

GO_result <- GO_result %>%
    mutate(Term = factor(Term, levels = term_order))

##p <- ggplot(GO_result, aes(x=Term, y=FDR, fill=sample))
##p <- p + geom_bar(stat = "identity",
##                  position = "dodge")
##p <- p + geom_hline(yintercept = -log10(0.05))
p <- ggplot(GO_result, aes(y = sample,
                           x = Term,
                           size = queryitem,
                           fill = FDR))
p <- p + geom_point(shape = 21)
p <- p + scale_fill_distiller(palette = 'Blues',
                              name = '-log10(FDR)')
p <- p + scale_size_binned(name = 'NumGenes')
##continuous(breaks=c(5, 100, 1000, 2000, 3000))

if(opt$options$category == "All"){
    p <- p + facet_grid(~term_type,
                        scales = "free",
                        space = "free")
}
p <- p + xlab("GO_terms [P]") + ylab("Category")
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank(),
               axis.text.x = element_text(angle=90,
                                          vjust = 0.5,
                                          hjust =1),
               axis.text = element_text(color = 'black'))

p <- p + theme(legend.position = c(opt$options$legendX, opt$options$legendY),
               legend.box.just = "bottom",
               legend.direction = "horizontal")
ggsave(file = opt$options$output,
       width = opt$options$width,
       height = opt$options$height)
