#!/usr/bin/Rscript
data <- read.table(file("stdin"))
ca <- commandArgs(trailingOnly=TRUE)
png(file=ca[1])
barplot(data[[3]])
dev.off()

