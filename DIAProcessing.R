
### Evaluation of DIA data
## O3 Batch Experiments
## May 2018
## J.E.Schollee

library(MSnbase)
library(lattice)
library(reshape)
library(devtools)
library(RMassScreening)
library(XLConnect)
library(cluster)
library(pheatmap)
library(vegan)
library(ape)

pkg <- "C:/DATA/R/MSMSsim"
load_all(pkg)
document(pkg)

base_dir <- "C:\\DATA\\BAFU\\O3Batch\\"
#raw_dir <- paste(base_dir, "mzXML\\pos-RP-DIA", sep="");
raw_dir <- paste(base_dir, "mzXML\\neg-RP-DIA-rerun", sep="");

setwd(base_dir)
source("DIAFunctions.R")
source("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\O3 Batch Experiments\\6_NontargetScreening\\ProcessingFunctions_JES.R")

files <- list.files(raw_dir, pattern="mzXML",full.names=TRUE, recursive=FALSE)

#sampleList <- readWorksheetFromFile("input/pos-RP/sampleListNew_pos_JES.xlsx", 1)
#sampleGroup <- readWorksheetFromFile("input/pos-RP/sampleGroups_JES.xlsx", 1)

sampleList <- readWorksheetFromFile("input/neg-RP/sampleListNew_neg_JES.xlsx", 1)
sampleGroup <- readWorksheetFromFile("input/neg-RP/sampleGroups_JES.xlsx", 1)

sampleList <- assignSamples(files, sampleList)  

## DIA fragment peak picking from Michele
#pickedDir <- paste0(base_dir, "picked/pos-DIA")
pickedDir <- paste0(base_dir, "picked/neg-DIA")
#outDir <- paste0(base_dir, "results/pos-DIA")
outDir <- paste0(base_dir, "results/neg-DIA")
dir.create(pickedDir)
dir.create(outDir)

#loadRmsSettings("settings-iso.ini")
loadRmsSettings("settings-iso-negESI.ini")

# pick each DIA slice individually
runs <- batchPickDIA(files, pickedDir, writeData=TRUE, writeList = TRUE, multicore = FALSE)

# "downstream" stuff if you want to try it
# profile each DIA slice individually
profiles <- lapply(runs, function(run) fillProfiles(pickedDir, files, pattern = run))
profiles <- lapply(profiles, function(profile) computeProfiles(profile, dret=40, dmass=4))
save(profiles, file = paste0(outDir,"\\profiles.RData"))

profiles <- assignProfiles(profiles, assignSamples(files, sampleList))
save(profiles, file = paste0(outDir, "\\profilesAssigned.RData"))

# remove profiles with a single occurrence across all files (optional step to be done with many files to reduce data amount)
# profiles <- cutSingleProfiles(profiles)

# extract sample x profile table from profiles, combining the DIA slices together in one table
# you are probably most interested in "fullInfo", which has RT, mz and DIA slice number for each peak
# you can of course do the profiling with a single file if you don't want the issues from mixing info from different samples
profMatrix <- makeProfileMatrix(profiles, type = "graphviz")
fullMatrix <- do.call(cbind, profMatrix)
profInfo <- makeProfileInfo(profiles)
fullInfo <- do.call(rbind, profInfo)


# build time trends
sampleAssignment <- readWorksheetFromFile("input/pos-RP/cultureAssignment_JES.xlsx", 1)
sampleGroups <- readWorksheetFromFile("input/pos-RP/sampleGroups_JES.xlsx", 1)

sampleAssignment <- readWorksheetFromFile("input/neg-RP/cultureAssignment_JES.xlsx", 1)
sampleGroups <- readWorksheetFromFile("input/neg-RP/sampleGroups_JES.xlsx", 1)

# Check whether the right samples will be grouped (note mock=TRUE)

colnames(profiles$peaks)[1] <- "m/z"

summaryTest <- lapply(profiles, function(prof){
  groupSummaries(
    prof, sampleList, sampleAssignment, sampleGroups, groupBy="time", groups=c("t0", "t1", "t3", "t7"), mock=TRUE
  )
})

summaryTables <- lapply(profiles, function(prof){
  groupSummaries(
    prof, sampleList, sampleAssignment, sampleGroups, groupBy="time", groups=c("t0", "t1", "t3", "t7")
  )
})

save(summaryTables, file=paste0(outDir,"/summaryTables.RData"))


# compute additional values for each sample group, at will.
# Here the mean is computed (note: currently only one-column addons are supported with the viewer)
# summaryTables <- lapply(summaryTables, function(df)
# {
#   within(df, {
#     mean = (t1 + t3) / 2
#     
#   })
# })

# merge to the "tt.total" table
totalTable <- list()
for(i in 1:length(profiles)){
  totalTable[[i]] <- mergeGroups(profiles[[i]], sampleGroups, summaryTables[[i]])
}

save(totalTable, file=paste0(outDir,"/totalTable.RData"))

timepoints <- data.frame(
  name = c("t0", "t1", "t3", "t7"),
  t=c(0,15,75,210)
)

############### simple characteristic comparison
sink(paste0(outDir, "/nfeat.txt"))
nfeat <- lapply(profiles, function(x){
  
  res <- table(x$peaks$sample, x$peaks$time)
  res <- data.frame(cbind(rownames(res), res))
  
})
nfeat

sink()

#nfeat <- do.call(rbind, nfeat)
#nfeat <- data.frame(cbind(rownames(nfeat), nfeat))
#nfeat <- merge(nfeat, sampleAssignment, by = "sample")

for(i in 1:length(nfeat)){
  
  res <- nfeat[[i]]
  
  res <- merge(res, sampleAssignment, by.x = "V1", by.y = "sample")
  res[,2:5] <- apply(res[,2:5], 2, function(x) as.numeric(as.character(x)))
  
  feat.means <- t(apply(res[,2:5], 2, function(x) tapply(x, res$sampleIdentifier, mean)))
  feat.sds <- t(apply(res[,2:5], 2, function(x) tapply(x, res$sampleIdentifier, sd)))
  
  png(paste0(outDir,"/Nfeat_BarChart_Scan", i, ".png"), res = 300, width = 3300, height = 1800)
  barCenters <- barplot(height = feat.means,
                        beside = TRUE, las = 1,
                        ylim = c(0, ifelse(i==1, 20000, 6000)),
                        #xlab = "Buechi Enrichment",
                        ylab = "Number of Non-target Features",
                        legend.text = TRUE,
                        #names.arg = c("Dry", "Wet"),
                        #args.legend = list(cex = 1, legend = c("Before CAS", expression("After O"[3]), "After GAK5b")),
                        col = rev(c(blues9[3:6])))
  segments(barCenters, feat.means - feat.sds * 2, barCenters,
           feat.means + feat.sds * 2, lwd = 1.5)
  arrows(barCenters, feat.means - feat.sds * 2, barCenters,
         feat.means + feat.sds * 2, lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  dev.off()
  
}
  


## take only fragments not MS1
fragMatrix <- fullMatrix[,which(fullInfo$scan != 1)]
fragInfo <- subset(fullInfo, fullInfo$scan != 1)
fragInfo$profile_ID <- c(paste(fragInfo$profile_ID, fragInfo$scan, sep = "_"))

for(i in 1:6){
  totalTable[[i]]$profileIDs <- c(paste(totalTable[[i]]$profileIDs, i, sep = "_"))
}

fragTable <- do.call(rbind, totalTable)
fragTable <- subset(fragTable, fullInfo$scan != 1)


dim(fragMatrix) ## 51  98529
dim(fragInfo) ## 98529  5
dim(fragTable) ## 98529  41

### analyze only fragments that can be annotated by the MassBank list
known.frag <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Abgeschlossene\\EDA-Emerge\\Data\\MassBank\\17092015\\fragmentAnalysis\\CountofFragmentsAll.csv", header = TRUE)
#known.frag <- known.frag[known.frag$eng == 45,]

NCEs <- c(15, 30, 45, 60, 75, 90)
known.frag <- known.frag[known.frag$eng %in% NCEs == 1,]
known.frag <- known.frag[!duplicated(known.frag$form),]

suspects <- data.frame(mass = known.frag$mz,
                       name = known.frag$form)

hits <- screenMass(peaklist = fragInfo, suspects = suspects)
hitsInfo <- merge(hits, fragInfo, by.x = "profileID", by.y = "profile_ID")

form.hits <- fragInfo$profile_ID %in% hits$profileID
form.hits <- which(form.hits == 1)

length(form.hits) ## 40406

fragInfo <- fragInfo[form.hits,]
fragMatrix <- fragMatrix[,form.hits]
fragTable <- fragTable[form.hits,]

dim(fragMatrix) ## 51  40406
dim(fragInfo) ## 404006  5
dim(fragTable) ## 40406  41

merge <- sapply(fragInfo$profile_ID,  function(x) which(hits$profileID == x))

merge <- lapply(merge, function(x){
  #hits[x,]
  names(which.max(table(hits[x,]$name)))
  }
) 

fragInfo$form <- unlist(merge)


###analyze all fragments together and then each scan individually)
pdf(paste0(outDir, "//FragmentCharperScan.pdf"))
for(i in 2:6){
  plot(y=fragInfo$mean_mz[which(fragInfo$scan ==i)], x=fragInfo$mean_RT[which(fragInfo$scan==i)], col=log10(fragInfo$mean_int[which(fragInfo$scan==i)]),
       main = paste("Scan Number", i))
}
dev.off()

png(paste0(outDir, "//ColLegend.png"), res = 300, width = 500, height = 5000)
test <- seq(3, 10, 1)
plot.new()
plot.window(xlim = c(0,1), ylim = c(2.5,10))
#points(y = test, x = rep(0.5, length(test)), col = test, pch = 17, cex = 2)
#abline(h = test, col = "lightgrey", lwd = 3)
abline(h = test, col = test, lwd = 60)
axis(2, at = seq(2,10,1), labels = seq(2,10,1))
title(ylab = "log10(intensity)")
box()
dev.off()

#xyplot(mean_mz ~ mean_RT | scan, data = fullInfo, col = log10(fullInfo$mean_int))
#densityplot( ~ mean_mz | scan,  data = fullInfo)

detected <- fragMatrix > 0
detected <- detected[-1,]
frag.count <- colSums(detected)
plot(x=fragInfo$mean_mz, y=fragInfo$mean_RT, col=frag.count)

pdf(paste0(outDir, "/FragsvsCountperScan.pdf"))
for(i in 2:6){
  plot(x=fragInfo$mean_mz[which(fragInfo$scan==i)], y=frag.count[which(fragInfo$scan==i)], main = paste("Scan Number", i))
  plot(x=fragInfo$mean_RT[which(fragInfo$scan==i)], y=frag.count[which(fragInfo$scan==i)], main = paste("Scan Number", i))
  plot(x=log10(fragInfo$mean_int)[which(fragInfo$scan==i)], y=frag.count[which(fragInfo$scan==i)], main = paste("Scan Number", i))
}
dev.off()

pdf(paste0(outDir,"/PCAplotsperScan.pdf"), paper = "a4r")
for(i in 2:6){
  pca <- prcomp(scale(fragMatrix[-1,which(fragInfo$scan==i)])) 
  plot(pca)
  plot(pca$rotation[,1:2])
  plot(pca$x, col = as.factor(sampleList$time) , pch = as.numeric(as.factor(sampleList$sample)) #, cex = as.factor(sampleList$spike.level)
  )
  legend("topleft", legend = c(unique(sampleList$time), unique(sampleList$sample)), col = c(unique(as.factor(sampleList$time)), rep("black", length(unique(sampleList$sample)))), 
         pch = c(rep(16, length(unique(sampleList$time))), as.numeric(unique(as.factor(sampleList$sample)))), ncol = 2)
  pairs(pca$x[,1:5], col = as.factor(sampleList$time) , pch = as.numeric(as.factor(sampleList$sample)), cex = 1.5)
}
dev.off()


### blank subtraction

blank <- which(sampleList$sample == "blank")

blank.cor <- apply(fragMatrix[-1,], 2, function(feature){
  
  blank <- mean(feature[blank])
  feature > (blank * 3)
})

blank.cor <- blank.cor * fragMatrix[-1,]

NDs <- which(colSums(blank.cor) == 0) ## 1367
blank.cor <- blank.cor[,-NDs] ## 50  39039
blank.cor.info <- fragInfo[-NDs,] ## 39039  6
blank.cor.table <- fragTable[-NDs,]  ##39039  41
  
hist(colSums(blank.cor>0))
table(colSums(blank.cor>0))
summary(colSums(blank.cor>0))
hist(blank.cor.info$scan)

### frequency/replicate filter
time <- c("t0", "t1", "t3", "t7")
rep.cor <- blank.cor

for(i in 1:length(unique(sampleGroup$sampleIdentifier))){
  
  for(j in 1:length(time)){
    
    take <- intersect(grep(sampleList$sample, pattern = unique(sampleGroup$sampleIdentifier[i])), grep(sampleList$time, pattern = time[j]))
    reps <- which(colSums(rep.cor[take,] > 0) < 2)
    rep.cor[take,reps] <- 0.0
  }
}

NDs <- which(colSums(rep.cor) == 0) ## n=3 31810; n=2 23593
rep.cor <- rep.cor[,-NDs] ## 50 7229/15446
rep.cor.info <- blank.cor.info[-NDs,] ## 7229/15446  6
rep.cor.table <- blank.cor.table[-NDs,] ## 7229/15446  41

#rep.cor <- blank.cor
#rep.cor.info <- blank.cor.info
#rep.cor.table <- blank.cor.table

hist(colSums(rep.cor>0))
table(colSums(rep.cor>0))
summary(colSums(rep.cor>0))
hist(rep.cor.info$scan)

save(rep.cor, file = paste0(outDir, "\\ReplicateCorrected.Rdata"))
save(rep.cor.info, file = paste0(outDir, "\\ReplicateCorrectedInfo.Rdata"))

## general characteristics
sink(paste0(outDir, "\\FragsperScan_RepCor.txt"))
table(rep.cor.info$scan)
for(i in 2:6){
  print(summary(rep.cor.info[rep.cor.info$scan==i,]))
}
sink()


### general plotting of results
pdf(paste0(outDir,"/PCAplots_MS2perScan_RepCor.pdf"), paper = "a4r")
for(i in 2:6){
  pca <- prcomp(scale(rep.cor[,which(rep.cor.info$scan==i)])) 
  plot(pca, sub = paste("Scan Number", i))
  plot(pca$rotation[,1:2], sub = paste("Scan Number", i))
  plot(pca$x, col = as.factor(sampleList$time) , pch = as.numeric(as.factor(substr(sampleList$sample, start = 1, stop = 4))),
       sub = paste("Scan Number", i) #, cex = as.factor(sampleList$spike.level)
  )
  text(pca$x, labels = seq(nrow(pca$x)), pos = 3, cex = 0.7)
  legend("topright", legend = c(unique(sampleList$time), unique(substr(sampleList$sample,1,4))), col = c(unique(as.factor(sampleList$time)), rep("black", length(unique(sampleList$sample)))), 
         pch = c(rep(16, length(unique(sampleList$time))), as.numeric(as.factor(unique(substr(sampleList$sample, start = 1, stop = 4))))), ncol = 2)
  pairs(pca$x[,1:5], col = as.factor(sampleList$time) , pch = as.numeric(as.factor(sampleList$sample)), cex = 1.5)
}
dev.off()


## if using only known fragments, plotting of formula distributions
n.form <- length(unique(rep.cor.info$form))
w.form <- unique(rep.cor.info$form)

pdf(paste0(outDir,"/FormulaDistributions_RepCor.pdf"), paper = "a4r")
par(mfrow = c(2,2))
for(i in 1:n.form){
  data <- subset(rep.cor.info, rep.cor.info$form == w.form[i])
  plot(data$mean_RT, log10(data$mean_int), col = data$scan, main = paste(w.form[i], "\n", round(mean(data$mean_mz), 4), "+/-", round(sd(data$mean_mz),8), sep = " "), pch = 16)
  legend(x = "topleft", legend = unique(data$scan), col = data$scan, pch = 16, title = "Scan")
}
dev.off()

for(i in 2:6){
  scan.take <- subset(rep.cor.info, rep.cor.info$scan==i)
  print(length(unique(scan.take$form)))
}

### calculate dbe
rep.cor.info$dbe <- sapply(rep.cor.info$form, function(x) dbe(x))
rep.cor.info$scan2 <- as.character(rep.cor.info$scan)

library(lattice)
pdf(paste0(outDir, "/DBEDist_RepCor.pdf"), paper = "a4r")
densityplot(~dbe, groups=scan2, data = rep.cor.info, auto.key = TRUE)
histogram(~dbe|scan2, data = rep.cor.info)
for(i in 2:6){
  hist(rep.cor.info$dbe[rep.cor.info$scan==i], main = paste("Scan Number", i))
  
}
dev.off()


### Van Krevelan plots

miss.form <- which(apply(rep.cor.info,1,function(row) row[["form"]] == ""))
rep.cor.info$form[miss.form] <- "NA"

elems <-  lapply(rep.cor.info$form, function(x) formulastring.to.list(x))
elems_fixed = do.call(rbind, lapply(elems, function(e) c(unlist(e)[c('C','H','O')])))
elems_fixed[is.na(elems_fixed)] <- 0
colnames(elems_fixed) <- c("C", "H", "O")
elems_fixed <- as.data.frame(elems_fixed)

elems_fixed$ratioHC <- c(elems_fixed[,"H"] / elems_fixed[,"C"])
elems_fixed$ratioOC <- elems_fixed[,"O"] / elems_fixed[,"C"]

#miss.form <- which(apply(rep.cor.info,1,function(row) row[["form"]] == ""))
#elems_fixed <- cbind(elems_fixed, rep.cor.info[-miss.form, ])

elems_fixed <- cbind(elems_fixed, rep.cor.info)
elems_fixed$MD <- elems_fixed$mean_mz - round(elems_fixed$mean_mz, 0)


pdf(paste0(outDir, "//FormDist_RepCor.pdf"), paper = "a4r")
xyplot(ratioOC ~ ratioHC | scan2, data = elems_fixed, 
       xlab = list(label = list("H/C")), ylab = list(label = list("O/C"))
       #,par.settings = list(fill = as.factor(log10(elems_fixed$mean_int)))
       )
xyplot(dbe ~ C | scan2, data = elems_fixed)
histogram(~MD | scan2, data = elems_fixed)
xyplot(O ~ MD | scan2, data = elems_fixed)
xyplot(MD ~ mean_mz | scan2, groups = O, data = elems_fixed
       , key = list(space = "top", columns = 5 
                    ,text = list(as.character(c(1:9)))
                    ,points = list(col = as.factor(c(1:9)))
                    )
)
dev.off()

## need to replace zeros for log10 calc
rep.cor[rep.cor==0] <- 0.01

rownames(rep.cor) <- paste(sampleList$sample, sampleList$time, sep = "-")

pdf(paste0(outDir, "\\PheatMap3_perScan_RepCor.pdf"))
for(i in 2:6){
  pheatmap(t(log10(rep.cor[,which(rep.cor.info$scan==i)])), sub = paste("Scan Number", i))
}
dev.off()

# rep.cor <- rep.cor[,which(colSums(rep.cor>0.1)>3)]
# 
# pdf(paste0(outDir, "\\PheatMap_repCor_n6.pdf"))
# pheatmap(t(log10(rep.cor)))
# dev.off()
# 

### correlated fragments
most.cor <- mosthighlycorrelated(rep.cor, 1000)


# # build time trends
# sampleAssignment <- readWorksheetFromFile("input/pos-RP/cultureAssignment_JES.xlsx", 1)
# sampleGroups <- readWorksheetFromFile("input/pos-RP/sampleGroups_JES.xlsx", 1)
# 
# sampleAssignment <- readWorksheetFromFile("input/neg-RP/cultureAssignment_JES.xlsx", 1)
# sampleGroups <- readWorksheetFromFile("input/neg-RP/sampleGroups_JES.xlsx", 1)
# 
# # Check whether the right samples will be grouped (note mock=TRUE)
# 
# colnames(profiles$peaks)[1] <- "m/z"
# 
# summaryTest <- lapply(profiles, function(prof){
#   groupSummaries(
#     prof, sampleList, sampleAssignment, sampleGroups, groupBy="time", groups=c("t0", "t1", "t3", "t7"), mock=TRUE
#   )
# })
# 
# summaryTables <- lapply(profiles, function(prof){
#   groupSummaries(
#     prof, sampleList, sampleAssignment, sampleGroups, groupBy="time", groups=c("t0", "t1", "t3", "t7")
#   )
# })
# 
# save(summaryTables, file=paste0(outDir,"/summaryTables.RData"))
# 
# 
# # compute additional values for each sample group, at will.
# # Here the mean is computed (note: currently only one-column addons are supported with the viewer)
# # summaryTables <- lapply(summaryTables, function(df)
# # {
# #   within(df, {
# #     mean = (t1 + t3) / 2
# #     
# #   })
# # })
# 
# # merge to the "tt.total" table
# totalTable <- list()
# for(i in 1:length(profiles)){
#   totalTable[[i]] <- mergeGroups(profiles[[i]], sampleGroups, summaryTables[[i]])
# }
# 
# save(totalTable, file=paste0(outDir,"/totalTable.RData"))
# 
# timepoints <- data.frame(
#   name = c("t0", "t1", "t3", "t7"),
#   t=c(0,15,75,210)
# )


### hits table for most common fragments MassBank
## fragments with only CH
# suspects <- data.frame(
#   mass = c(91.05418, 77.03834, 58.06512, 65.0385, 57.06988, 105.0699, 79.05419),
#   name = c("C7H7", "C6H5", "C3H8N", "C5H5", "C4H9", "C8H9", "C6H7")
# )
# 
# ## fragments with CHO
# suspects <- data.frame(
#   mass = c(95.04914, 81.03349, 53.0022, 107.0491, 55.0178, 145.0647),
#   name = c("C6H7O", "C5H5O", "C3HO", "C7H7O", "C3H3O", "C10H9O")
# )
# 
# hits <- lapply(profiles, function(prof) screenProfiles(prof, suspects, "+", 5))
# do.call(rbind, hits)

#runViewer(totalTable, hits, timepoints ,sampleGroups, launch.browser=TRUE)


## find fragments with specific mz
spec.frag <- lapply(profiles, function(prof){
  plotTable(prof, mzrange = rev(RMassBank::ppm(mass = 176.0320, dppm=5, l=TRUE)), mode = "+")
})


runViewer(totalTable[[2]], spec.frag[[2]], timepoints ,sampleGroups, launch.browser=TRUE)


spec.frag <- do.call(rbind, spec.frag)
plot(spec.frag$RT, log10(spec.frag$int))


### clustering for fragment characterization
normmax_pos <- apply(rep.cor, 2, function(feature){
    feature/max(feature)
})

m <- as.matrix(normmax_pos)
my.dist <- function(c) dist(c, method="euclidean")
my.clust <- function(d) hclust(d, method="ward.D2")

# png(file = paste0(outDir, "//HeatMapFrags_pos.png"))
# heatmap(x=t(m), Colv=NA, distfun=my.dist, hclustfun=my.clust, scale='none',
#         col=topo.colors(100))
# dev.off()

dis_normmax_pos <- daisy(t(normmax_pos), metric = "euclidean", stand = FALSE)
cluster_ward_normmax_pos <- hclust(dis_normmax_pos, method="ward.D2")

### Cutting the tree either at specific hight or number of clusters 
kp <- 20                                        # No. of clusters
method_pos <- cluster_ward_normmax_pos
method_name_pos <- "cluster_ward_normmax_pos"
y_value_pos <- normmax_pos

r.cutree_pos <- cutree(method_pos, k=kp)
tab_pos <- table(r.cutree_pos)


# ##### Plotting the Clusters positive
# 
ids <- unlist(strsplit(colnames(y_value_pos), "_"))
ids <- ids[seq(1,length(ids), 2)]
 
 
 pdf(file=paste0(outDir,"\\", method_name_pos,"_",kp, ".pdf"),  paper="a4", pointsize=10, width=8, height=11)
 par(mfrow=c(4,2))
 
 plot(method_pos,  main = "Clusteranlyse", hang=0.1, labels=cutree(method_pos, k=kp), las=2)
 
 plot(method_pos,  main = "Clusteranlyse", hang=0.1, labels=cutree(method_pos, k=kp), las=2)
rect.hclust(method_pos, k=kp, border="red")
 
 for (i in 1:kp){
   plot(x=c(1,50), y=c(0,1), type="n", main=paste("normalized cluster",i),
        xlab="sample sequence", ylab="Intensity", axes = FALSE)
   axis(1, at = seq_len(50), cex = 0.2, las = 2)
   axis(2)
   
   list_in_i <- which(r.cutree_pos==i)
   
   for (j in 1:tab_pos[i]){
     lines(x=1:50, y=y_value_pos[1:50,list_in_i[j]])
   }
}
 dev.off()

 
pdf(file=paste0(outDir,"\\FragChar_", method_name_pos,"_",kp, ".pdf"),  paper="a4", pointsize=10, width=8, height=11)
par(mfrow = c(4,2))
for (i in 1:kp){
  plot(x=c(200,1500), y=c(50,1000), type="n", main=paste("normalized cluster",i),
       xlab="RT (sec)", ylab="m/z")
  #axis(1, at = seq_len(50), cex = 0.2, las = 2)
  #axis(2)
  list_in_i <- which(r.cutree_pos==i)
  points(x=rep.cor.info$mean_RT[list_in_i], y=rep.cor.info$mean_mz[list_in_i])
}

for(i in 1:kp){
  list_in_i <- which(r.cutree_pos==i)
  
  hist(rep.cor.info$mean_RT[list_in_i], breaks = 50, main = paste("RT - Cluster", i))
  hist(rep.cor.info$mean_mz[list_in_i], breaks = 100, main = paste("m/z - Cluster", i))
}

dev.off()

### remove blank samples from table
no.blank <- which(sampleList$sample != "blank")
NB.rep.cor <- rep.cor[no.blank,]
NB.rep.info <- rep.cor.info[,no.blank]

sampleList2 <- sampleList[no.blank,]

pdf(paste0(outDir,"/PCAplots_MS2perScan_NBRepCor.pdf"), paper = "a4r")
for(i in 2:6){
  pca <- prcomp(scale(NB.rep.cor[,which(rep.cor.info$scan==i)])) 
  plot(pca, sub = paste("Scan Number", i))
  plot(pca$rotation[,1:2], sub = paste("Scan Number", i))
  plot(pca$x, col = as.factor(sampleList2$time) , pch = as.numeric(as.factor(substr(sampleList2$sample, start = 1, stop = 4))),
       sub = paste("Scan Number", i) #, cex = as.factor(sampleList$spike.level)
  )
  text(pca$x, labels = seq(nrow(pca$x)), pos = 3, cex = 0.7)
  legend("topright", legend = c(unique(sampleList2$time), unique(substr(sampleList2$sample,1,4))), col = c(unique(as.factor(sampleList2$time)), rep("black", length(unique(sampleList2$sample)))), 
         pch = c(rep(16, length(unique(sampleList2$time))), as.numeric(as.factor(unique(substr(sampleList2$sample, start = 1, stop = 4))))), ncol = 2)
  pairs(pca$x[,1:5], col = as.factor(sampleList2$time) , pch = as.numeric(as.factor(sampleList2$sample)), cex = 1.5)
}
dev.off()

pdf(paste0(outDir, "/HeatMaps_NBRepCor.pdf"), paper = "a4r")
pheatmap(NB.rep.cor, scale = "column", main = "NBRepCor - Col Scaled")
pheatmap(log10(NB.rep.cor), main = "NBRepCor - log10 No Scale")
pheatmap(log10(NB.rep.cor), scale = "column", main = "NBRepCor - log10 Col Scaled")
dev.off()


##########

## adonis analysis

#no.blank <- which(sampleList$sample != "blank")
rep.dist <- vegdist(NB.rep.cor)

adonis.var <- sampleList2
#adonis.var$time[c(1,12)] <- "blank"
adonis.var$group <- substr(adonis.var$sample, 1,4)
#adonis.var$group[c(1,12)] <- "blank"
adonis.var$day <- substr(adonis.var$sample, 1, 2)
#adonis.var$day[c(1,12)] <- "blank"
adonis.var$spike <- substr(adonis.var$sample, 4,4)
#adonis.var$spike[c(1,12)] <- "blank"
adonis.var$label <- paste0(adonis.var$group, "-", adonis.var$time)

sink(paste0(outDir, "\\AdonisResults.txt"))
ado.res <- adonis2(rep.dist~time/group, strata = group, data = adonis.var)
ado.res
ado.res <- adonis2(rep.dist~day*time*spike, data = adonis.var)
ado.res
ado.res <- adonis2(rep.dist~label, data = adonis.var)
ado.res
mod <- with(adonis.var, betadisper(rep.dist, label, bias.adjust = FALSE))
permutest(mod)
anova(mod)
TukeyHSD(mod)
sink()

pdf(paste0(outDir, "\\AdonisResults.pdf"))
plot(mod)
boxplot(mod, col = as.factor(unique(adonis.var$time)))
dev.off()



#############
### other clustering options

pdf(paste0(outDir, "\\AsPhylo_NBRepCor.pdf"))
hc <- hclust(rep.dist, method = "ward.D2")
hc$labels <- adonis.var$label
plot(as.phylo(hc), type = "unrooted", cex = 0.6, no.margin = TRUE)
dev.off()


#### CCA
cca <- cca(NB.rep.cor~day+time+spike, data = adonis.var)

sink(paste0(outDir, "\\CCAResults.txt"))
vif.cca(cca)
permutest(cca)
cca
step(cca, scope = NB.rep.cor ~., test = "perm", perm.max  = 100)
anova(cca, by = "terms")
anova(cca, by = "mar")
anova(cca, by = "axis")
summary(cca)
sink()

pdf(paste0(outDir, "\\CCAResults.pdf"))
plot(cca)
dev.off()

### RDA when using the MassBank constrained data
rda.data <- NB.rep.cor
rda <- rda(rda.data~time+day+spike, adonis.var)
plot(rda)



rda.data <- t(NB.rep.cor)
rda <- rda(rda.data~scan2+dbe+C+O, elems_fixed)
plot(rda)




### plotting changes over ozonation

library(wesanderson)
library(RColorBrewer)

group <- unique(adonis.var$group)
colo <- wes_palette("Darjeeling1", n= 4, type = "discrete")
colo <- blues9[3:6]
colo <- brewer.pal(4,"BrBG")


for(k in 2:6){
  
  want.scan <- which(rep.cor.info$scan == k)
  pdf(paste0(outDir, "/vanKrevelen_Scan", k, ".pdf"), paper = "a4r", width = 1600)
  
  
  for(j in 1:length(group)){
    
    take.group <- group[j]
    
    plot.new()
    plot.window(xlim = c(0,2), ylim = c(0,4))
    grid(col = "lightgrey")
    for(i in 1:length(time)){
      
      take.time <- time[i]
      take <- intersect(which(adonis.var$group == take.group), which(adonis.var$time == take.time))
      #take <- NB.rep.cor[take,-miss.form]
      take <- NB.rep.cor[take, want.scan]
      
      detect <- which(colSums(take>0.01) > 1)
      points(y = jitter(elems_fixed$ratioHC[detect]), x = jitter(elems_fixed$ratioOC[detect]), col = colo[i],
             bg = alpha(colo[i], 0.4), pch = 21, cex = 1.2 #rev(seq(.8,1.2,.1))[i]
             )
      
    }
    axis(1, cex.axis = 1.5)
    axis(2, cex.axis = 1.5)
    box()
    title(main = paste("Van Krevelen \n", take.group, sep = ""), ylab = "H/C", xlab = "O/C", cex.lab = 1.5, sub = paste("Scan", k))
    legend("topleft", legend = time, col = colo, pch = 16)
    
  }
  
  dev.off()
}



elems_fixed$KendMass.O <- elems_fixed$mean_mz * (16 / 15.994915)
elems_fixed$KendMassDef.O <- round(elems_fixed$KendMass.O,0) - elems_fixed$KendMass.O 


for(k in 2:6){
  
  want.scan <- which(rep.cor.info$scan == k)
  pdf(paste0(outDir, "/KendrickMD_O_Scan", k, ".pdf"), paper = "a4r", width = 1600)
  
for(j in 1:length(group)){
  
  take.group <- group[j]
  
  plot.new()
  plot.window(xlim = c(-0.3,0.2), ylim = c(0,200))
  grid(col = "lightgrey")
  for(i in 1:length(time)){
    
    take.time <- time[i]
    take <- intersect(which(adonis.var$group == take.group), which(adonis.var$time == take.time))
    #take <- NB.rep.cor[take,-miss.form]
    take <- NB.rep.cor[take, want.scan]
    detect <- which(colSums(take>0.01) > 1)
    points(y = elems_fixed$KendMass.O[detect], x = elems_fixed$KendMassDef.O[detect], col = colo[i],
           bg = alpha(colo[i], 0.4), pch = 21, cex = 1.2)
    
  }
  axis(1, cex.axis = 1.5)
  axis(2, cex.axis = 1.5)
  box()
  title(main = paste("Kendrick O \n", take.group, sep = ""), ylab = "Kendrick Nominal Mass", 
        xlab = "Kendrick Mass Defect", cex.lab = 1.5)
  legend("topleft", legend = time, col = colo, pch = 16)
  
}
dev.off()

}


for(k in 2:6){
  
  want.scan <- which(rep.cor.info$scan == k)
  pdf(paste0(outDir, "/CvsDBE_Scan", k, ".pdf"), paper = "a4r", width = 1600)
  
  for(j in 1:length(group)){
    
    take.group <- group[j]
    
    plot.new()
    plot.window(xlim = c(0,15), ylim = c(-0.5,15))
    grid(col = "lightgrey")
    for(i in 1:length(time)){
      
      take.time <- time[i]
      take <- intersect(which(adonis.var$group == take.group), which(adonis.var$time == take.time))
      #take <- NB.rep.cor[take,-miss.form]
      take <- NB.rep.cor[take, want.scan]
      detect <- which(colSums(take>0.01) > 1)
      points(x = elems_fixed$C[detect]+ seq(-.225,.225,.15)[i], y = jitter(elems_fixed$dbe[detect]), col = colo[i],
             bg = alpha(colo[i], 0.4), pch = 21, cex = 1.2)
      
    }
    axis(1, cex.axis = 1.5, at = c(0:15))
    axis(2, cex.axis = 1.5)
    box()
    title(main = paste("Carbon Number vs. DBE \n", take.group, sep = ""), xlab = "Carbon Number", 
          ylab = "Double Bond Equilavent (DBE)", cex.lab = 1.5)
    legend("topleft", legend = time, col = colo, pch = 16)
    
  }
  dev.off()
  
}


for(k in 2:6){
  
  want.scan <- which(rep.cor.info$scan == k)
  pdf(paste0(outDir, "/HistO_Scan", k, ".pdf"), paper = "a4r", width = 1600)
  
  for(j in 1:length(group)){
    
    take.group <- group[j]
    
    #plot.new()
    #plot.window(xlim = c(-0.3,0.2), ylim = c(0,200))
    #grid(col = "lightgrey")
    for(i in 1:length(time)){
      
      take.time <- time[i]
      take <- intersect(which(adonis.var$group == take.group), which(adonis.var$time == take.time))
      #take <- NB.rep.cor[take,-miss.form]
      take <- NB.rep.cor[take, want.scan]
      detect <- which(colSums(take>0.01) > 1)
      #points(y = elems_fixed$KendMass.O[detect], x = elems_fixed$KendMassDef.O[detect], col = colo[i],
      #       bg = alpha(colo[i], 0.4), pch = 21, cex = 1.2)
      hist(elems_fixed$ratioOC[detect], breaks = c(seq(0,1,.02),5), xlim = c(0,1), ylim = c(0,30), main = paste(take.group, take.time))
    }
    #axis(1, cex.axis = 1.5)
    #axis(2, cex.axis = 1.5)
    #box()
    #title(main = paste("Kendrick O \n", take.group, sep = ""), ylab = "Kendrick Nominal Mass", 
    #      xlab = "Kendrick Mass Defect", cex.lab = 1.5)
    #legend("topleft", legend = time, col = colo, pch = 16)
    
  }
  dev.off()
  
}


library(gplots)

pdf(paste0(outDir, "\\Venn_NFrag_RepCor.pdf"), paper = "a4r")
WK.U.feats <- data.frame(t0 = rep.cor.table$`t0.WK-U` > 0,
                         t1 = rep.cor.table$`t1.WK-U` > 0,
                         t3 = rep.cor.table$`t3.WK-U` > 0,
                         t7 = rep.cor.table$`t7.WK-U` > 0
)


venn(WK.U.feats)
title(main = "WK-U")


WK.S.feats <- data.frame(t0 = rep.cor.table$`t0.WK-S` > 0,
                         t1 = rep.cor.table$`t1.WK-S` > 0,
                         t3 = rep.cor.table$`t3.WK-S` > 0,
                         t7 = rep.cor.table$`t7.WK-S` > 0
)


venn(WK.S.feats)
title(main = "WK-S")

WD.U.feats <- data.frame(t0 = rep.cor.table$`t0.WD-U` > 0,
                         t1 = rep.cor.table$`t1.WD-U` > 0,
                         t3 = rep.cor.table$`t3.WD-U` > 0,
                         t7 = rep.cor.table$`t7.WD-U` > 0
)

venn(WD.U.feats)
title(main = "WD-U")

WD.S.feats <- data.frame(t0 = rep.cor.table$`t0.WD-S` > 0,
                         t1 = rep.cor.table$`t1.WD-S` > 0,
                         t3 = rep.cor.table$`t3.WD-S` > 0,
                         t7 = rep.cor.table$`t7.WD-S` > 0
)


venn(WD.S.feats)
title(main = "WD-S")

dev.off()


sink(paste0(outDir,"/sessionInfo.txt"))
session_info()
sink()

save.image(paste0(outDir, "//", Sys.Date(), "-Rdata.RData"))



#all.data <- cbind(elems_fixed, t(NB.rep.cor))

clusters <- list()
clusters[[1]] <- NULL

for(i in 2:6){
 res <- pheatmap(t(log10(rep.cor[,which(rep.cor.info$scan==i)])))
 clusters[[i]] <- cbind(rep.cor.info[which(rep.cor.info$scan==i),], cluster = cutree(res$tree_row, k = 20))
 pheatmap(t(log10(rep.cor[,which(rep.cor.info$scan==i)])), annotation_row =  as.data.frame(cutree(res$tree_row, k=20)))
}

all.clusters <- do.call(rbind, clusters)

pdf(paste0(outDir, "\\FormperClust.pdf"), paper = "a4r")
par(mfrow = c(2,2))
for(k in 1:20){
  
  dat <- table(clusters[[5]][which(clusters[[5]]$cluster==k),"form"])
  n.feat <- nrow(clusters[[5]][which(clusters[[5]]$cluster==k),])
  #hist(dat, main = paste("Scan 5, Cluster", k))
  ifelse(length(which(dat>1))>0,
         pie.data <- subset(dat, dat > 1), pie.data <- dat)
  pie(pie.data, main = paste("Scan 5, Cluster", k), sub = paste("N Cluster = ", n.feat, "\n N Pie Data = ", sum(pie.data), "\n N form =", length(pie.data)))
}
dev.off()
# xyplot(mean_RT~mean_mz | cluster + scan2, groups = log10(mean_int), data = clusters)



#### Calculate neutral losses
# apply(data[[1]]$ms2, FUN = function(x) getFragMz(x))
# 
# specs <- sapply(data[[1]]$ms2, function(d) as.data.frame(d))
# 
# 
# all.frag <- do.call(rbind, all.frag)
# histogram( ~mz | file, data = all.frag)
# 
# mymean <- function(x) {ifelse(is.na(mean(x)), 0.01, mean(x))}
# pca.data <- cast(data = all.frag, formula = mz ~ file, value = "i", function(x) mymean(x))
# rownames(pca.data) <- pca.data$mz
# pca.data <- pca.data[,2:5]
# head(pca.data)
# 
# 
# pca <- prcomp(pca.data)
# plot(pca)
# biplot(pca)
# 
# 
# 
# # 
# all.frag <- list()
# frags <- list()
# 
# for(i in 1:length(data)){
# 
#   for(j in 1:10 #length(data[[1]]$ms2)
#       ){
# 
#     spec <- as.data.frame.Spectrum(data[[i]]$ms2[[j]])
#     spec$scan <- scanIndex(data[[i]]$ms2[j])
#     frags[[j]] <- spec
#   }
# 
#   all.frag[[i]] <- do.call(rbind, frags)
#   all.frag[[i]]$file <- i
# }
# 
# #optimize fragment extraction
# total.frags <- sum(data[[1]]$hd.ms2$originalPeaksCount)
# frags <- seq_len(total.frags)
# nfrags <- data[[1]]$hd.ms2$originalPeaksCount
# start <- seq_len(length(nfrags))
# end <- seq_len(length(nfrags))
# 
# #seq_along(data[[1]]$ms2)
# res <- lapply(1:50, function(i)
# {
#   as.data.frame(data[[1]]$ms2[[i]])
# })
# for(i in 1:length(nfrags)){
# 
#   ifelse(i == 1, start[i] <- 1, start[i] <- end[i-1] + 1)
#   ifelse(i == 1, end[i] <- nfrags[i], end[i] <- start[i] + nfrags[i] - 1)
# 
# }
# 
# MS2.spec <- list()
# for(i in 1:length(data[[1]]$ms2)){
#   MS2.spec[[i]] <- data[[1]]$ms2[[i]]
# }
# 
