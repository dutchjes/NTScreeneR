

#### making profiles list into feature table. Need profiled enviMass data to start


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

