

#### setting up a profile list for enviMass profiling

library(devtools)
install_github("blosloos/enviMass", force = TRUE)
library(enviMass)

### load in your data. Structure of data is a list, where each element of the list is one peak-picked file, ie. enviPick output

load("data.rds")
picked <- data
makeProfileList <- function(picked){
  
  profiles<-list(0)
  profiles[[1]]<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
  colnames(profiles[[1]])<-c("peaks?","agglom?","profiled?","trends?")
  profiles[[2]]<-0  # peaks
  profiles[[3]]<-0  # datetime
  profiles[[4]]<-0  # time
  profiles[[5]]<-0  # place
  profiles[[6]]<-0  # index_agglom
  profiles[[7]]<-0  # index_prof
  profiles[[8]]<-0  # parameters
  profiles[[9]]<-0  # sample type
  names(profiles)<-c("state","peaks","datetime","sampleID","place",
                     "index_agglom","index_prof","parameters","type")
  
  data <- list()
  leng <- length(picked)
  at<-0
  
  for(i in 1:leng){ 
    at<-c(at+length(picked[[i]][,1]))  	
  }
  
  peaks<-matrix(nrow=at,ncol=8,0)
  colnames(peaks)<-c("m/z","intensity","RT","peakIDs","componentIDs","sampleIDs","partitionIDs","profileIDs")
  da1<-1;
  for(i in 1:leng){ 
    peaklist <- as.matrix(picked[[i]]);
    peaklist<-cbind(peaklist,peaklist[,c(1,4,5)])
    colnames(peaklist)[12]<-"m/z_corr";
    colnames(peaklist)[13]<-"sum_int_corr";
    colnames(peaklist)[14]<-"RT_corr";                
    da2<-(da1+length(peaklist[,1])-1)
    that<-c(length(peaklist[,1]))
    peaks[da1:da2,]<-as.matrix(cbind(peaklist[,c(12,13,14,10)],
                                     rep(0,that),rep(i,that),
                                     rep(0,that),rep(0,that)
    ));
    da1<-c(da2+1);
  }
  
  peaks<-peaks[order(peaks[,1],decreasing=FALSE),]
  profiles[[2]]<-peaks;
  datetime<-as.POSIXct.numeric(seq(1,leng,1),origin = "1960-01-01")
  profiles[[3]]<-datetime;
  rm(peaks)
  return(profiles)
}

profiles <- makeProfileList(picked)

profiles<-agglomer(
  profiles,
  dmass=15,
  ppm=TRUE,
  dret=50
)

profiles<-partcluster(
  profiles,
  dmass=15,
  ppm=TRUE,
  dret=50,
  from=FALSE,
  to=FALSE,
  progbar=FALSE,
  plot_it=FALSE,
  replicates=FALSE,
  IDs=FALSE,
  with_test=FALSE
)

matrix <- data.frame(profiletopeak(profiles, progbar = TRUE))


#### building data frame of peak intensities in samples

sampleIDs <- unique(profiles$peaks[,"sampleIDs"])
sampleIDs <- sampleIDs[order(sampleIDs)]
num <- max(matrix$profileID)

peaks <- as.data.frame(matrix(data = 0, ncol = length(sampleIDs)+3, nrow = num))
colnames(peaks) <- c("profileID", "mz", "rt", sampleIDs)
peaks$profileID <- matrix$profileID
peaks$mz <- matrix$mean_m.z
peaks$rt <- matrix$mean_RT

findId <- function(x) which(colnames(peaks) == x)

for(i in 1:num){
  
  if(i==1){
    print("Filling Peaks Table...")
    pb <- txtProgressBar(min = 0, max = num)
  }
  
  setTxtProgressBar(pb, i)
  prof <- peaks$profileID[i]
  detects <- subset(profiles$peaks, profiles$peaks[,"profileIDs"] == prof)
  
  samps <- sapply(detects[,"sampleIDs"], FUN = findId)
  peaks[i,samps] <- detects[,"intensity"]
  
  
  if(i==length(peaks$profileID)){close(pb)}
}




##### making a feature list for further analysis

num <- max(profiles$peaks[,"profileIDs"])
count <- table(profiles$peaks[,"profileIDs"])
feature <- list()

for(i in 1:num){
  
  if(i==1){
    print("Filling Features List...")
    pb <- txtProgressBar(min = 0, max = num)
  }
  
  setTxtProgressBar(pb, i)
  
  count_ID <- i
  prof <- as.data.frame(profiles[[2]][
    profiles[[7]][count_ID,1]:profiles[[7]][count_ID,2]
    ,])
  
  
  if(count[count_ID] > 1){
    prof <- prof[order(prof[,"sampleIDs"],decreasing = FALSE),]
    feature[[i]] <- prof[,c("peakIDs", "m/z", "intensity", "RT", "sampleIDs")]
    feature[[i]][,"RT"] <- feature[[i]][,"RT"]/60
    feature[[i]][,"sampleIDs"] <- feature[[i]][,"sampleIDs"]
    colnames(feature[[i]]) <- c("ID", "m/z", "intensity", "RT", "sample")
  }else{
    prof <- as.data.frame(t(prof))
    prof <- prof[order(prof[,"sampleIDs"],decreasing = FALSE),]
    feature[[i]] <- prof[,c("peakIDs", "m/z", "intensity", "RT", "sampleIDs")]
    feature[[i]][,"RT"] <- feature[[i]][,"RT"]/60
    feature[[i]][,"sampleIDs"] <- feature[[i]][,"sampleIDs"]
    colnames(feature[[i]]) <- c("ID", "m/z", "intensity", "RT", "sample")    
  }
  
  if(i==length(peaks$profileID)){close(pb)}
  
}

