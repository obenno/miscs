#! /usr/bin/Rscript

library(ggtree)
args <- commandArgs(trailingOnly = T)

# input1 is tree file,
# input2 is taxa category of genes,
# input3 is output file name
tree <- read.tree(args[1])
group <- read.table(args[2], header = F)
names(group) <- c("id","taxa")
cls <- list()
for(i in 1:length(levels(group$taxa))){
    cls[[i]] <- as.vector(group[group$taxa==levels(group$taxa)[i],1])
}
tree <- groupOTU(tree, cls)
p <- ggtree(tree, aes(color=group),
            layout = "circular",
            branch.length="none")
p <- p+geom_tiplab2(size=2.5, aes(angle=angle), color="dimgrey")

ggsave(p, filename = args[3], width = 12, height = 12)
