
library(RChemMass)
## functions for processing HRMS data
## written by J.E.Schollee

## loading multiple packages simultaneously
Install_And_Load <- function(packages, repos = repos, type = type) {
  k <- packages[!(packages %in% installed.packages()[,"Package"])];
  if(length(k))
  {install.packages(k, repos=repos, type = type)};
  
  for(package_name in packages)
  {library(package_name,character.only=TRUE, quietly = TRUE);}
};


## wrapper for getting SMILES from chemical formula
## function from RChemMass
FormulafromSMILES <- function(SMILES){
  
  x <- getSuspectFormulaMass(as.character(SMILES))
  return(x$MolForm)
}

## wrapper for formatting for SMILES
## functions from rcdk
read.smiles <- function(smiles){
  
  x <- parse.smiles(as.character(smiles), kekulise = FALSE)[[1]]
  do.aromaticity(x)
  do.typing(x)
  do.isotopes(x)
  return(x)
  
}

## wrapper for calculating charge of compound
## functions from rcdk
calcCharge <- function(x){
  
  y <- get.atoms(x)
  l <- lapply(y, function(l) get.formal.charge(l))
  ch <- sum(unlist(l))
  return(ch)
  
}

## for calculating exact mass from SMILES
exactMass <- function(smiles){
  
  #x <- parse.smiles(as.character(smiles), kekulise = FALSE)[[1]]
  #do.aromaticity(x)
  #do.typing(x)
  #do.isotopes(x)
  x <- read.smiles(smiles)
  ch <- calcCharge(x)
  if(ch<0){y <- get.exact.mass(x) + p.mass * abs(ch)}
  if(ch>0){y <- get.exact.mass(x) - p.mass * abs(ch)}
  if(ch==0){y <- get.exact.mass(x)}
  
  #y <- get.exact.mass(x)
  
  return(y)
}

## calculate optimum collision energy based on compound mass
## developed in Uchem based on target compounds
## published in Gulde 2014
optCE <- function(mass){
  
  topCE <- c()
  if(mass < 350){
    topCE <- mass*(-0.41)+160
    topCE <- round(topCE*20, digits = -2)/20
  }else{
    topCE <- 15
  }
  
  return(topCE)
}

## wrapper function for predicting logP from smiles
## functions from rcdk
predlogP <- function(smiles){
  
  #x <- parse.smiles(as.character(x), kekulise = FALSE)[[1]]
  #do.aromaticity(x)
  #do.typing(x)
  #do.isotopes(x)
  x <- read.smiles(smiles)
  p <- get.xlogp(x)
  return(p)
}

## calculate RT match
## returns TRUE/FALSE if suspect match is within allow RT window
RTfilter <- function(expected.rt, measured.rt, rttol){
  
  match <- c()
  if(measured.rt > expected.rt - rttol && measured.rt < expected.rt + rttol)
    match <- TRUE
  else
    match <- FALSE
  
  return(match)
}

## exact mass suspect screening
## from Michele, probably built into RMassScreening also...
screenMass <- function(peaklist, suspects, ppmLimit = 5, neutral.mass = FALSE, polarity = "+"){
  
  if (is.null(ppmLimit)) 
    stop("No ppm limit specified. Specify either as a parameter or in your settings file.")
  potential <- suspects

  if(neutral.mass == TRUE){
    if (polarity == "+") 
      potential$mass <- potential$mass + 1.0072
    else if (polarity == "-") 
      potential$mass <- potential$mass - 1.0072
  }

  potential$X <- NULL
  ppmMax <- ppmLimit
  hits <- lapply(1:nrow(potential), function(n) {
    mass <- potential$mass[[n]]
    dppm <- ((peaklist$mean_mz/mass) - 1) * 1e+06
    hit <- peaklist[abs(dppm) < ppmMax, c("profile_ID", "mean_mz", 
                                          "mean_RT", "mean_int"), drop = FALSE]
    colnames(hit) <- c("profileID", "mz", "RT", "int")
    df <- cbind(potential[rep(n, nrow(hit)), , drop = FALSE], 
                hit)
  })
  hits.total <- do.call(rbind, hits)
  hits.total$dppm <- (hits.total$mz/hits.total$mass - 1) * 
    1e+06
  if (!is.null(hits.total$ret)) {
    hits.total$dRT <- (hits.total$RT - hits.total$ret)
    if (!is.null(rtLimit)) {
      rtStrict <- (sign(rtLimit) == 1)
      hits.total <- hits.total[(!is.na(hits.total$dRT) | 
                                  !rtStrict), , drop = FALSE]
      hits.total <- hits.total[(is.na(hits.total$dRT)) | 
                                 (abs(hits.total$dRT) < abs(rtLimit)), , drop = FALSE]
    }
  }
  return(hits.total)
  
}

## recalibration of picked files using internal standards
## functions from enviMass, adapted for mz
recalibrateMZ <- function(file, scan, pickedDir, recalDir, recalTable, tolmz = 10, tolret = 120)
{
  # determine in and out filenames
  suffix <- ".MSlist.RData"
  if(scan == "+")
    suffix <- ".MSlist.pos.RData"
  else if(scan == "-")
    suffix <- ".MSlist.neg.RData"
  
  
  filePicked <- paste0(pickedDir, "/", basename(file), suffix )
  fileOut <- paste0(recalDir, "/", basename(file), suffix )
  pathPlot.ret <- paste0(recalDir, "/", basename(file),  ".MZ.png")
  pathModel.ret <- paste0(recalDir, "/", basename(file), ".MZ.RData")
  
  load(filePicked)  
  
  res <- recalibCE(MSlist$Peaklist[,c("m/z","sum_int", "RT")]
                   , recalTable[,1], tolmz = tolmz, ppm=TRUE, ret = 60*recalTable[,2], tolret = tolret,
                   what="mass", plotit = TRUE, one=TRUE, path_1 = pathPlot.ret, path_2 = pathModel.ret)
  
  MSlist$Peaklist[,c("m/z","sum_int", "RT")] <- res
  save(MSlist, file=fileOut)
}


## function 1 for grouping profiles which should belong together
## applied after normal profiling to combat "peak splitting" and poor intensity repeatability

profileGroup <- function(profiles, dppm = 1, drt = 10){
  
  l <- nrow(profiles$index_prof)
  peaklist <- as.data.frame(profiles$index_prof[,c("mean_mz", "mean_RT", "max_int", "number_peaks_total")])
  peaklist$grouped <- NA
  rownames(peaklist) <- profiles$index_prof[,"profile_ID"]
  colnames(peaklist) <- c("mean_mz", "mean_RT", "max_int", "number_peaks_total", "grouped")
  peaklist <- peaklist[order(peaklist[,3], decreasing = TRUE),]
  n.peaks <- peaklist[,4]
  peaklist[,1] <- as.numeric(as.character(peaklist[,1]))
  peaklist[,2] <- as.numeric(as.character(peaklist[,2]))
  peaklist[,3] <- as.numeric(as.character(peaklist[,3]))
  peaklist[,4] <- as.numeric(as.character(peaklist[,4]))
  
  groupID <- list()
  
  pb <- txtProgressBar(min = 0, max = l)
  
  for(i in 1:l){
    
    #print(i)
    setTxtProgressBar(pb,i)
    
    if(!is.na(peaklist[i,5])){next} ## checks if profile already in a group
    
    if(n.peaks[i] == 1){  ## no grouping if profile only has 1 peak in it
      
      groupID[[i]] <- paste(rownames(peaklist[i,]))
      next
      } 

    mz.min <- peaklist[i,1] - RMassScreening::ppm(peaklist[i,1], dppm)
    mz.max <- peaklist[i,1] + RMassScreening::ppm(peaklist[i,1], dppm)
    
    take1 <- which(peaklist[,1] < mz.max & peaklist[,1] > mz.min)
    
    rt.min <- peaklist[i,2] - drt
    rt.max <- peaklist[i,2] + drt
    
    take2 <- which(peaklist[,2] < rt.max & peaklist[,2] > rt.min)
    
    take <- intersect(take1, take2)
    
    groupID[[i]] <- paste(rownames(peaklist[take,]))
    
    if(length(take) == 1){next}else{

      #peaklist[i,5] <- paste(rownames(peaklist[take,]), collapse = ';')
       ## indicates profileIDs to be grouped
      
      peaklist[take,5] <- "x" ## marks profiles which are already in a group
      
      }
    
  }
  
  close(pb)
  return(groupID)
  
}



## calculate the variability of profile characteristics in the composite samples
## select either m/z, rt, int
## select either relative standard deviation or standard deviation
compVariability <- function(profiles, sampleCat, value, type){
  
  #if(sum(c("mz", "rt", "int" %in% value != 1)))
  #  stop("Choose either mz, rt, or int as value to calculate variability")
  
  l <- length(unique(profiles$peaks$profileIDs))
  u <- unique(profiles$peaks$profileIDs)
  repIDs <- unique(sampleCat$sampleCat[which(sampleCat$rep)])
  n.rep <- length(unique(repIDs))
  
  varib <- data.frame(matrix(NA, nrow = l, ncol = n.rep))
  rownames(varib) <- u
  colnames(varib) <- repIDs
  
  pb <- txtProgressBar(min = 0, max = l)
  
  for(i in 1:l){
    
    setTxtProgressBar(pb,i)
    prof <- profiles$peaks[which(profiles$peaks$profileIDs == u[i]),]

    var <- rep(NA, n.rep)
    
    for(j in 1:n.rep){
      
      data <- prof[which(prof$sampleCat == repIDs[j]),]
      if(nrow(data)<3){next}
      else{
        
        if(value == "mz"){
          data <- data$mz
        }
        if(value == "rt"){
          data <- data$RT
        }
        if(value == "int"){
          data <- data$intensity
        }
        
      }
      
      if(type == "SD"){
        var[j] <- sqrt(var(data))
      }
      if(type == "RSD"){
        var[j] <- 100 * (sqrt(var(data)) / mean(data))
      }
      if(type == "mean"){
        var[j] <- mean(data)
      }
      if(type == "median"){
        var[j] <- median(data)
      }
      if(type == "max"){
        var[j] <- max(data)
      }
      if(type == "min"){
        var[j] <- min(data)
      }
    }
    varib[i,] <- var
  }
  
  close(pb)
  return(varib)
}



## filtering for features that are increased by a specific factor
## developed to look for OTPs in a batch experiment
IntIncrease <- function(totalTable, nm.OTP.cols, nm.parent.col, factor.increase = 2){
  
  OTP.sums <- rowSums(totalTable[,colnames(totalTable) %in% nm.OTP.cols])
  par.sum <- totalTable[,colnames(totalTable) %in% nm.parent.col]
  posOTP <- OTP.sums > par.sum * factor.increase 
  return(posOTP)
  
}

## returns the p-value of various significance tests
ks.p <- function(x=x, y=y, ...){
  ks <- wilcox.test(x=x, y=y, ...)
  # ks <- ks.test(x=x, y=y, ...)
  # ks <- t.text(x=x, y=y, ...)
  return(ks$p.value)
}

## calculates the p-value of intensity differences between groups of samples
## groups defined by the sampleList
getSignif <- function(profilesRc, sampleList, sample.x = "", time.x = "", sample.y = "", time.y = ""){
  
  take.x <- intersect(grep(sample.x, sampleList$sample), grep(time.x, sampleList$time))
  take.y <- intersect(grep(sample.y, sampleList$sample), grep(time.y, sampleList$time))
  
  samps.x <- profilesRc[take.x,]
  samps.y <- profilesRc[take.y,]
  
  samps <- rbind(samps.x, samps.y)
  
  #ks.pvalues <- apply(samps[,], 2, function(data) {ks.p(x=data[1:nrow(samps.x)], y=data[c(nrow(samps.x)+1):c(nrow(samps.y)+nrow(samps.x))])})
  
  y <- c(rep("1", length(take.x)), rep("2", length(take.y)))
  ks.pvalues <- apply(samps, 2, function(data) {
    t<-t.test(data~y) 
    return(t$p.value)}
  )
  
  return(ks.pvalues)
}

## calculates the fold change of intensity differences between groups of samples
## groups defined by the sampleList
getFC <- function(profilesRc, sampleList, sample.x = "", time.x = "", sample.y = "", time.y = ""){
  
  take.x <- intersect(grep(sample.x, sampleList$sample), grep(time.x, sampleList$time))
  take.y <- intersect(grep(sample.y, sampleList$sample), grep(time.y, sampleList$time))
  
  samps.x <- profilesRc[take.x,]
  samps.y <- profilesRc[take.y,]
  
  fc <- colMeans(samps.x) / colMeans(samps.y)
  
  return(fc)
}

##returns features in a particular quadrant of the volcano plot
VolcanoQuad <- function(pvalues, foldchange, p.thres = "", fc.thres = "", fc.quad = c("higher", "lower")){
  
  p.sig <- which(-log10(pvalues) > c(-log10(p.thres)))
  fc.sig <- if(fc.quad == "higher") {
    which(log2(foldchange) > log2(fc.thres))
  }else{
    which(log2(foldchange) < c(-log2(fc.thres)))
    }
  
  not.inf <- which(!is.infinite(log2(foldchange)) > 0)
  fc.sig <- intersect(not.inf, fc.sig)
  
  sig.feat <- intersect(p.sig, fc.sig)
  return(sig.feat)
}



## function for plotting the  4 Venn graph
## includes calculation of overlapping feature counts then plot output
## requires "VennDiagram" package
plot4Venn <- function(venn.data, ...){
  
  n12 <- length(intersect(venn.data[[1]], venn.data[[2]]))
  n13 <- length(intersect(venn.data[[1]], venn.data[[3]]))
  n14 <- length(intersect(venn.data[[1]], venn.data[[4]]))
  n23 <- length(intersect(venn.data[[2]], venn.data[[3]]))
  n24 <- length(intersect(venn.data[[2]], venn.data[[4]]))
  n34 <- length(intersect(venn.data[[3]], venn.data[[4]]))
  n123 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[3]])))
  n124 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[4]])))
  n134 <- length(intersect(venn.data[[1]], intersect(venn.data[[3]], venn.data[[4]])))
  n234 <- length(intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[4]])))
  n1234 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[4]]))))
  
  
  venn.plot <- draw.quad.venn(
    area1 = length(venn.data[[1]]),
    area2 = length(venn.data[[2]]),
    area3 = length(venn.data[[3]]), 
    area4 = length(venn.data[[4]]),
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23  = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    ...
  )
  
  return(venn.plot)
  
}

## function for plotting the  3 Venn graph
## includes calculation of overlapping feature counts then plot output
## requires "VennDiagram" package
plot3Venn <- function(venn.data, ...){
  
  n12 <- length(intersect(venn.data[[1]], venn.data[[2]]))
  n13 <- length(intersect(venn.data[[1]], venn.data[[3]]))
  n23 <- length(intersect(venn.data[[2]], venn.data[[3]]))
  n123 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[3]])))

  
  venn.plot <- draw.triple.venn(
    area1 = length(venn.data[[1]]),
    area2 = length(venn.data[[2]]),
    area3 = length(venn.data[[3]]), 
    n12 = n12,
    n13 = n13,
    n23  = n23,
    n123 = n123,
    ...
  )
  
  return(venn.plot)
  
}


## function for plotting the 5 venn graph
## includes calculation of overlapping feature counts then plot output
## requires "VennDiagram" package
plot5venn <- function(venn.data, ...){
  
  n12 <- length(intersect(venn.data[[1]], venn.data[[2]]))
  n13 <- length(intersect(venn.data[[1]], venn.data[[3]]))
  n14 <- length(intersect(venn.data[[1]], venn.data[[4]]))
  n15 <- length(intersect(venn.data[[1]], venn.data[[5]]))
  n23 <- length(intersect(venn.data[[2]], venn.data[[3]]))
  n24 <- length(intersect(venn.data[[2]], venn.data[[4]]))
  n25 <- length(intersect(venn.data[[2]], venn.data[[5]]))
  n34 <- length(intersect(venn.data[[3]], venn.data[[4]]))
  n35 <- length(intersect(venn.data[[3]], venn.data[[5]]))
  n45 <- length(intersect(venn.data[[4]], venn.data[[5]]))
  n123 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[3]])))
  n124 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[4]])))
  n125 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], venn.data[[5]])))
  n134 <- length(intersect(venn.data[[1]], intersect(venn.data[[3]], venn.data[[4]])))
  n135 <- length(intersect(venn.data[[1]], intersect(venn.data[[3]], venn.data[[5]])))
  n145 <- length(intersect(venn.data[[1]], intersect(venn.data[[4]], venn.data[[5]])))
  n234 <- length(intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[4]])))
  n235 <- length(intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[5]])))
  n245 <- length(intersect(venn.data[[2]], intersect(venn.data[[4]], venn.data[[5]])))
  n345 <- length(intersect(venn.data[[3]], intersect(venn.data[[4]], venn.data[[5]])))  
  n1234 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[4]]))))
  n1235	<- length(intersect(venn.data[[1]], intersect(venn.data[[2]], intersect(venn.data[[3]], venn.data[[5]]))))
  n1245	<- length(intersect(venn.data[[1]], intersect(venn.data[[2]], intersect(venn.data[[4]], venn.data[[5]]))))
  n1345	<- length(intersect(venn.data[[1]], intersect(venn.data[[3]], intersect(venn.data[[4]], venn.data[[5]]))))
  n2345	<- length(intersect(venn.data[[2]], intersect(venn.data[[3]], intersect(venn.data[[4]], venn.data[[5]]))))
  n12345 <- length(intersect(venn.data[[1]], intersect(venn.data[[2]], intersect(venn.data[[3]], intersect(venn.data[[4]], venn.data[[5]])))))	
  
  
  venn.plot <- draw.quintuple.venn(
    area1 = length(venn.data[[1]]),
    area2 = length(venn.data[[2]]),
    area3 = length(venn.data[[3]]), 
    area4 = length(venn.data[[4]]),
    area5 = length(venn.data[[5]]),
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n15 = n15,
    n23  = n23,
    n24 = n24,
    n25 = n25,
    n34 = n34,
    n35 = n35,
    n45 = n45,
    n123 = n123,
    n124 = n124,
    n125 = n125,
    n134 = n134,
    n135 = n135,
    n145 = n145,
    n234 = n234,
    n235 = n235,
    n245 = n245,
    n345 = n345,
    n1234 = n1234,
    n1235 = n1235,
    n1245 = n1245,
    n1345 = n1345,
    n2345 = n2345,
    n12345 = n12345,
    ...
  )
  
  return(venn.plot)
}


## plot any desired feature with RMassScreening
plotTable <- function(profiles, mzrange = c(10,1000), rtrange = c(0, 3000), mode = c("+", "-")){
  
  min.mz <- mzrange[1]
  max.mz <- mzrange[2]
  sel.mz <- which(profiles$index_prof[,"mean_mz"] > min.mz & profiles$index_prof[,"mean_mz"] < max.mz)
  
  min.rt <- rtrange[1]
  max.rt <- rtrange[2]
  sel.rt <- which(profiles$index_prof[,"min_RT"] > min.rt & profiles$index_prof[,"max_RT"] < max.rt)
  
  sel.feat <- intersect(sel.mz, sel.rt)
  
  hits <- as.data.frame(profiles$index_prof[sel.feat,])
  hits <- data.frame(
    profileID = hits$profile_ID,
    name = hits$profile_ID,
    mass = ifelse(mode == "+", hits$mean_mz + 1.0072, hits$mean_mz - 1.0072),
    dppm = 0,
    mz = hits$mean_mz,
    RT = hits$mean_RT,
    int = hits$mean_int
  )
  
  return(hits)
  
}

## function to fit the proper distribution fit
## developed for LOD estimation but ultimately not implemented.
fitData <- function(data, fit="gamma", sample=0.5){
  distrib = list()
  numfit <- length(fit)
  results = matrix(0, ncol=5, nrow=numfit)
  
  for(i in 1:numfit){
    if((fit[i] == "gamma") | 
       (fit[i] == "poisson") | 
       (fit[i] == "weibull") | 
       (fit[i] == "exponential") |
       (fit[i] == "logistic") |
       (fit[i] == "normal") | 
       (fit[i] == "geometric")
    ) 
      distrib[[i]] = fit[i]
    else stop("Provide a valid distribution to fit data" )
  }
  
  # take a sample of dataset
  n = round(length(data)*sample)
  data = sample(data, size=n, replace=F)
  
  for(i in 1:numfit) {
    if(distrib[[i]] == "gamma") {
      gf_shape = "gamma"
      fd_g <- fitdistr(data, "gamma")
      est_shape = fd_g$estimate[[1]]
      est_rate = fd_g$estimate[[2]]
      
      ks = ks.test(data, "pgamma", shape=est_shape, rate=est_rate)
      
      # add to results
      results[i,] = c(gf_shape, est_shape, est_rate, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "poisson"){
      gf_shape = "poisson"
      fd_p <- fitdistr(data, "poisson")
      est_lambda = fd_p$estimate[[1]]
      
      ks = ks.test(data, "ppois", lambda=est_lambda)
      # add to results
      results[i,] = c(gf_shape, est_lambda, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "weibull"){
      gf_shape = "weibull"
      fd_w <- fitdistr(data,densfun=dweibull,start=list(scale=1,shape=2))
      est_shape = fd_w$estimate[[1]]
      est_scale = fd_w$estimate[[2]]
      
      ks = ks.test(data, "pweibull", shape=est_shape, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_shape, est_scale, ks$statistic, ks$p.value) 
    }
    
    else if(distrib[[i]] == "normal"){
      gf_shape = "normal"
      fd_n <- fitdistr(data, "normal")
      est_mean = fd_n$estimate[[1]]
      est_sd = fd_n$estimate[[2]]
      
      ks = ks.test(data, "pnorm", mean=est_mean, sd=est_sd)
      # add to results
      results[i,] = c(gf_shape, est_mean, est_sd, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "exponential"){
      gf_shape = "exponential"
      fd_e <- fitdistr(data, "exponential")
      est_rate = fd_e$estimate[[1]]
      ks = ks.test(data, "pexp", rate=est_rate)
      # add to results
      results[i,] = c(gf_shape, est_rate, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "logistic"){
      gf_shape = "logistic"
      fd_l <- fitdistr(data, "logistic")
      est_location = fd_l$estimate[[1]]
      est_scale = fd_l$estimate[[2]]
      ks = ks.test(data, "plogis", location=est_location, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_location, est_scale, ks$statistic,    ks$p.value) 
    }
  }
  results = rbind(c("distribution", "param1", "param2", "ks stat", "ks    pvalue"),   results)
  #print(results)
  return(results)
}
