#! /usr/bin/Rscript

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(doParallel))

args <- commandArgs(trailingOnly = T)
## args[1] GFF, args[2] is sample1 wig, args[3] is sample2 wig
## args[4] is mobile RNA list, args[5] is output pdf
txdb <- makeTxDbFromGFF(args[1])
sig <- import(args[2])
cds <- cdsBy(txdb)
fiveUTR <- fiveUTRsByTranscript(txdb)
threeUTR <- threeUTRsByTranscript(txdb)
mobileList <- read.table(args[4], header=F)



numCores <- detectCores()
registerDoParallel(numCores)
message("Used parallel with ", numCores, " thread(s).")

compute_binAverage <- function (inputGranges, signal_obj, binName_prefix, bin_number){
    ## inputGranges could be fiveUTR, threeUTR and cds Granges List
    ## signal file is imported wiggle/bigWiggle
    ## binName_prefix is "fiveUTR_", combined with bin_number to generate
    ## "fiveUTR_1", "fiveUTR_2"... bin name

    ## subset inputGranges according to singal_obj
    gr <- unlist(inputGranges)
    commonTrans <- intersect(seqnames(gr), seqnames(signal_obj))
    gr <- gr[seqnames(gr) %in% commonTrans]
    signal_obj <- signal_obj[seqnames(signal_obj) %in% commonTrans]
    outGranges <- foreach (i=1:length(gr), .combine=c) %dopar% {
        g0 <- gr[i]
        if(length(g0)==0 | sum(width(g0)<bin_number)){
            g2 <- g0[0, NULL]
        }else{
            g1 <- tile(g0, bin_number)[[1]]
            seqlevels(g1) <- unique(as.character(seqnames(g1)))
            sig_rle <- coverage(signal_obj[seqnames(signal_obj)==seqlevels(g1)], weight="score")[seqlevels(g1)]
            g2 <- binnedAverage(g1, sig_rle, "binScore")
            g2$bin <- paste0(binName_prefix, c(1:bin_number))
            g2
        }
    }
    return(outGranges)
}

message("Starting sample1", " [", Sys.time(), "]")
result_5UTR <- compute_binAverage(fiveUTR, sig, "fiveUTR_", 10)
message("Finished sample1 fiveUTR", " [", Sys.time(), "]")
result_cds <- compute_binAverage(cds, sig, "cds_", 15)
message("Finished sample1 cds", " [", Sys.time(), "]")
result_3UTR <- compute_binAverage(threeUTR, sig, "threeUTR_", 10)
message("Finished sample1 threeUTR", " [", Sys.time(), "]")


out_sample1 <- c(result_5UTR, result_cds, result_3UTR)
out_sample1 <- as.data.frame(out_sample1, row.names=NULL)
head(out_sample1)
out_sample1 <- out_sample1 %>%
    mutate(type=case_when(str_sub(seqnames,1,15) %in% mobileList$V1 ~ "mobile",
                          TRUE ~ "normal")) %>%
    mutate(bin=factor(bin, levels=c(paste0("fiveUTR_", c(1:10)),
                                    paste0("cds_", c(1:15)),
                                    paste0("threeUTR_", c(1:10))))) %>%
    group_by(type, bin) %>%
    summarise(mean = mean(binScore)) %>%
    ungroup()
out_sample1 <- out_sample1 %>%
    mutate(sample=rep("C08", nrow(out_sample1)))
head(out_sample1)

sig <- import(args[3])
message("Starting sample2", " [", Sys.time(), "]")
result_5UTR <- compute_binAverage(fiveUTR, sig, "fiveUTR_", 10)
message("Finished sample2 fiveUTR", " [", Sys.time(), "]")
result_cds <- compute_binAverage(cds, sig, "cds_", 15)
message("Finished sample2 cds", " [", Sys.time(), "]")
result_3UTR <- compute_binAverage(threeUTR, sig, "threeUTR_", 10)
message("Finished sample2 threeUTR", " [", Sys.time(), "]")


out_sample2 <- c(result_5UTR, result_cds, result_3UTR)
out_sample2 <- as.data.frame(out_sample2, row.names=NULL)
head(out_sample2)
out_sample2 <- out_sample2 %>%
    mutate(type=case_when(str_sub(seqnames,1,15) %in% mobileList$V1 ~ "mobile",
                          TRUE ~ "normal")) %>%
    mutate(bin=factor(bin, levels=c(paste0("fiveUTR_", c(1:10)),
                                    paste0("cds_", c(1:15)),
                                    paste0("threeUTR_", c(1:10))))) %>%
    group_by(type, bin) %>%
    summarise(mean = mean(binScore)) %>%
    ungroup()
out_sample2 <- out_sample2 %>%
    mutate(sample=rep("W05", nrow(out_sample1)))
head(out_sample2)

out <- rbind(out_sample1, out_sample2)

p <- ggplot(out, aes(x=bin, y=mean, group=interaction(sample, type), color=sample, linetype=type))
p <- p+geom_point()+geom_line()
p <- p + geom_vline(xintercept="cds_1") + geom_vline(xintercept="cds_15")
#p <- p + scale_x_discrete("", labels=NULL)
labels <- rep("", 35)
labels[1] <- "5' UTR"
labels[11] <- "ATG"
labels[25] <- "Stop_Codon"
labels[34] <- "3' UTR"
p <- p + scale_x_discrete("", labels=labels)
p <- p+theme_bw()+theme(panel.grid=element_blank(),
                        legend.position = "bottom",
                        axis.text.x=element_text(size=12, face="bold"))
p <- p+xlab("")+ylab("Fraction of 5mC")

ggsave(args[5], width=10, height=5)

write.table(out, file="tmp.txt", quote=F, sep="\t", row.names=F)

stopImplicitCluster()
