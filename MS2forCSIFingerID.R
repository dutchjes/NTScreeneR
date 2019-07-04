
## Script to prepare text files for batch CSI:FingerID analysis
## J.E.Schollee
## June 2018

library(MSnbase)
library(mzR)

dir <- "C:\\DATA\\NeugutOzonationNTS\\posttreat_clustering\\cluster6_NT"
setwd(dir)

text.dir <- "20150409-49\\MSMS"
out.dir <- "20150409-49\\MSfiles"
dir.create(paste(dir, out.dir, sep = "\\"))

files <- list.files(text.dir, pattern = ".txt", full.names = TRUE)
MSMS <- sapply(files, read.delim2)
sapply(files, writeMSData)


mat <- res[[55]]
spec <- new("Spectrum2", mz=mat[,1], intensity=mat[,2])
writeMgfData(spec, con = "testdata.mgf")

