#!/usr/bin/Rscript
data <- read.table(file("stdin"))
ca <- commandArgs(trailingOnly=TRUE)
pdf(file=ca[1], width=12, height=6)
barplot(data[[3]])
dev.off()

