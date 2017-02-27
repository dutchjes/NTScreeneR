
####

sample.info <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\EDA-Emerge\\Data\\Neugut\\post-treat\\pos\\sampleinfo.csv",
                        header = TRUE)
types <- unique(sample.info$SampleType)
types <- subset(types, types != "standard") ## removes the standards
sample.info <- subset(sample.info, sample.info$SampleType %in% types)

maindir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\EDA-Emerge\\Data\\Neugut\\post-treat\\pos\\03022016"
MSlistcsv <- "Q:/Abteilungsprojekte/uchem/Projekte/Aktuelle Projekte/EDA-Emerge/Data/Neugut/post-treat/pos/EnviPick/csv"

library(nontarget)


for(z in 2:length(types)){
  
  
posttreattype <- types[z]

if(posttreattype == "blind"){
  next
}

sample.info <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\EDA-Emerge\\Data\\Neugut\\post-treat\\pos\\sampleinfo.csv",
                        header = TRUE)
sample.info <- subset(sample.info, sample.info$SampleType == "blind" | sample.info$SampleType == posttreattype)

names <- as.data.frame(list.files(paste(MSlistcsv),full.names = FALSE))
#sample.info$measured <- paste(sample.info$Filename, ".mzXML.csv", sep = "") %in% names[,1]
names <- subset(names, rownames(names) %in% sample.info$rn)

#dir.create(paste(maindir, "\\", posttreattype, sep = ""))
setwd(paste(maindir, "\\", posttreattype, sep = ""))

peaklist <- read.csv("allpeaks.csv", header = TRUE)

####################################################################################################
####################  Start Nontarget script  ######################################################
####################################################################################################


#dir.create(paste(getwd(), "\\nontarget", sep = ""))
#setwd(paste(getwd(), "\\nontarget", sep = ""))

data.mono <- list()
## remove satellite peaks ##
#peaklist<-rm.sat(peaklist,dmz=0.3,drt=0.1,intrat=0.015,spar=0.8, corcut=-1000,plotit=TRUE);
#peaklist<-peaklist[peaklist[,4],1:3]

#for(i in 1:length(data)){
  
#  peaklist <- data[[i]]
nc <- ncol(peaklist)
peaklist$max_int <- c()

for(k in 1:nrow(peaklist)){
  peaklist$max_int[k] <- max(peaklist[k,5:nc])
}

  peaklist <- peaklist[,c("mz", "max_int", "rt")]
  colnames(peaklist) <- c("mass", "intensity", "rt")
  
  summary(peaklist)
  str(peaklist)
  
  ######## so run nontarget package first and try to find homologues with S
  data(adducts) # load adduct list
  data(isotopes) # load isotope list
  ######## isotope grouping - define isotopes and charge
  iso <- make.isos(isotopes,
                   use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C",
                                  "15N","34S","37Cl","81Br","41K"),
                   use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
  ######## isotope grouping - run on peaklist
  pattern <- pattern.search(peaklist, iso, cutint=10000, rttol=c(-0.5,0.5), mztol=10,
                            mzfrac=0.1, ppm=TRUE, inttol=0.2, rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                            deter=FALSE, entry=100);
  
  save(pattern, file = paste("pattern.rds", sep = ""))
  
  # ###################### isotopes: plotting ############################################
  # 
  pdf(paste("plotisopatt.pdf", sep = ""),paper="a4r")
  plotisotopes(pattern);
  dev.off()
  pdf(paste("plotdefectS.pdf", sep = ""),paper="a4r")
  plotdefect(pattern,elements=c("S"));
  dev.off()
  pdf(paste("plotdefectCl.pdf", sep = ""),paper="a4r")
  plotdefect(pattern,elements=c("Cl"));
  dev.off()
  pdf(paste("plotdefectN.pdf", sep = ""),paper="a4r")
  plotdefect(pattern,elements=c("N"));
  dev.off()
  
  ############################## adduct grouping and plotting ###################################
  ### if peaklist is from positive
  adduct <- adduct.search(peaklist, adducts, rttol=0.1, mztol=15,# massfrac=0.1, 
                          ppm=TRUE, use_adducts=c("M+K","M+H","M+Na","M+NH4"), 
                          ion_mode="positive" #, entry=50
  )
  
   save(adduct, file = paste("adduct.rds", sep = ""))
  # 
  # ############ maybe try adding more entries here to find the most relevant...
  # ################ more adducts in adduct_list
  pdf(paste("plotadduct.pdf", sep = ""),paper="a4r")
  plotadduct(adduct)
  dev.off()
  # ################ show a pattern group and relationship to adduct groups
  # plotall(pattern, adduct)
  
  
  ############ combine results into components
  ############ for those only with isotope pattern
  ############ without homologous series
  comp <- combine(pattern, adduct, homol=FALSE, dont=FALSE,  rules=c(TRUE,FALSE,FALSE))
  #plotcomp(comp,compoID=1,peakID=FALSE)
  
  #save(pattern, file = "pattern.rds")
  #save(adduct, file = "adduct.rds")
  save(comp, file = paste("comp.rds",sep = ""))
  
  #### export of results to have a look..
  fileNT <- paste("data_NT.txt", sep = "")
  filecomp1 <- paste("data_comp1.txt", sep = "")
  
  write.table(pattern[[1]],file=paste(fileNT),row.names=FALSE,quote=TRUE,sep=";")
  write.table(comp[[1]],file=paste(filecomp1),row.names=FALSE,quote=TRUE,sep=";")
  
  
  NTdata <- read.table(paste(fileNT), header = TRUE, sep = ";")
  #comp1data <- read.table(paste(filecomp1), header = TRUE, sep = ";")
  
 # mono <- which(NTdata$interaction.level == "0" | NTdata$interaction.level == "1")
  
  data.mono <- subset(NTdata, NTdata$interaction.level == "0" | NTdata$interaction.level == "1")
  
  save(data.mono, file = "datamono.rds")
  
  #data.mono[[i]] <- data[[i]][mono,]
  
  #NTdata <- NTdata[with(NTdata, order(mz)), ]
  #peaklist <- peaklist[with(peaklist, order(mz)), ]
  
  
}



