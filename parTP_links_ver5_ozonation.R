
#############################################################################################################
###################Element 1################################################################################
###################building of reaction table################################################################
#############################################################################################################
reactions.phase1 <- data.frame(
  index = c("oh", "deme", "deet", "deh2", "h2", "deh2o", 
            "clXh", "deco2", "deno2", "oxicooh", "oxicho", "h2o", "meoxi", 
            "deamin", "nitrored", "disf", "deet2", "disclox", "disfox", "disnhox", 
            "deipr", "acetyl", "deacteyl", "gluco", "degluco", "sulfo", "desulfo", "dealk", "descypr", "dioxi", "ozone", "oxideam", "deamket"), 
  reaction = c("Hydroxylation, N/S-Oxidation, Epoxidation (HOE)", 
               "Demethylation", "Di-demethy or Deethylation", "Dehydrogenation", 
               "Hydrogenation", "Dehydration", "Reductive displacement of chlorine ", 
               "Decarboxylation", "Loss of nitro group", "Alcohol to carboxylic acid or primary amine to nitro", 
               "Ketone formation", "Hydration", "Methyl to carboxylic acid", 
               "Deamination", "Nitro reduction", "Reductive displacement of fluorine", 
               "Desethylation", "Oxydative displacement of chlorine", "Oxydative displacement of fluorine", 
               "Oxydative displacement of amine", "Deisopropyl", "Acetylation", "De-acetylation", "Glucuronidation",
		           "De-glucuronidation", "Sulfation", "De-sulfation", "Dealkylation", "Descyclopropyl", "Oxygen addition", 
		           "Addition of ozone", "Oxidative deamination", "Deamination to ketone"),
  lossmass = as.numeric(c(0, 14.0157, 28.0313, 2.0157,0, 18.0106, 34.9689, 43.9898, 45.9929,2.0157,2.0157,0,2.0157, 15.0109,
                          31.9898, 18.9984,28.0313,34.9689,18.9984,15.0109,42.0470,0,42.01057,0,176.0320,0,79.9568, 30.0470,
                          40.0313, 0, 0, 17.0266, 17.0266)),
  gainmass = as.numeric(c(15.9949,0,0,0,2.0157,0,1.0078,0,1.0078,15.9949,0,18.0106,31.9898,0,2.0157,1.0078,0,17.0027,
                          17.0027, 17.0027,0,42.010565,0,176.0320,0,79.9568,0, 0, 0, 31.989830, 47.984745, 31.9898, 15.9949)),  
  loss = c("", 
           "CH2", "C2H4", "H2", "", "H2O", "Cl", "CO2", "NO2", "H2", "H2", 
           "", "H2", "NH", "O2", "F", "C2H4", "Cl", "F", "NH2", "C3H6", 
           "" , "C2H20", "", "C6H8O6", "", "SO3", "C2H6", "C3H4", "", "", "NH3", "NH3"), 
  gain = c("O", "", "", "", "H2", "", "H", "", "H", "O", 
           "", "H2O", "O2", "", "H2", "H", "", "OH", "OH", "OH", "", "C2H2O", "", "C6H8O6", "", "SO3", "", "", "", "O2", "O3", "O2", "O"), 
  #dmass = c(15.99491462, -14.01565006, -28.03130013, -2.015650064, 
  #  	2.015650064, -18.01056468, -33.96102765, -43.98982924, -44.98507821, 
  #		13.97926456, -2.015650064, 18.01056468, 29.97417917, -15.01089904, 
  #		-29.97417917, -17.99057819, -28.03130013, -17.96611303, -1.995663575, 
  #		0.984015583, -42.04695019, -150.0680796), 
  active = c("+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+",
             "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+") 
)


reactions.phase1$massdiff <- NULL
for(i in 1:nrow(reactions.phase1)){
  reactions.phase1$massdiff[i] <- reactions.phase1$gainmass[i]-reactions.phase1$lossmass[i]
}

act <- subset(reactions.phase1,reactions.phase1[,"active"]!="")
x <- act[,"index"]

############################################################################################################
######Linkage function using the determined parent and TP peaks#############################################
######par.peaks and TP.peaks come from tests_merge scripts##################################################
############################################################################################################
#load("data.rds")

####NOTE:: sample.info table must be in the working directory. CHANGE IF NECESSARY!!#######

linkages <- function(act, par.peaks, TP.peaks, tolmz = 10, ppm = TRUE, RTfilter = TRUE, filename = "links"){
  
  print(paste("Using working directory ", getwd(), sep=""))
#  dir <- list.files(path = getwd(), all.files = TRUE)
#  if("shortfeatlist.rds" %in% dir == FALSE){stop("Add shortfeatlist.rds to directory!")}
#  if("sampleinfo.csv" %in% dir == FALSE){stop("Add sampleinfo.csv to directory!")}
#  if("data.rds" %in% dir == FALSE){stop("Add data.rds to directory!")}
  
  x <- act[,"index"]
  options <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      options[j,k] <- paste(rownames(par.peaks[j,]), x[k])}
  }
  
  option.masses <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  rownames(option.masses) <- rownames(par.peaks)
  colnames(option.masses) <- act[,"index"]
  
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      option.masses[j,k] <- par.peaks[j,"mz"] + act[k,"massdiff"]}
  }
  
  save(option.masses, file="optionmasses2.rds")
  save(options, file="optionnames2.rds")
  print("Saved option masses")
  print("Saved option names")

  potential <- as.vector(as.matrix(option.masses))
  p <- as.vector(as.matrix(options))
  names(potential) <- p
  
  potential <- sort(potential, decreasing = FALSE)
  
  back <- strsplit(names(potential), split = " ")
  
  links <- data.frame(matrix(NA, nrow = length(potential), ncol = 11))
  rownames(links) = names(potential)
  colnames(links) = c("Link.data ID", "Parent Peak ID", "TP peak ID", "transformation", 
                      "parent m/z", "proposed TP m/z", "TP peaks m/z", "Parent RT", "TP RT", "Parent Int", "TP Int")
  target_up <- rep(NA, length(potential)) 
  target_down <- rep(NA, length(potential))
  
  prog<-winProgressBar("Calculate Potential TP Masses",min=0,max=length(potential));
  
  for(j in 1:length(potential)){   
    setWinProgressBar(prog, j, title = NULL, label = NULL)   
    # calculate bounds per target mass #################################
    if(ppm==TRUE){
      target_down[j]<-as.numeric(potential[j])-((as.numeric(tolmz)*(as.numeric(potential[j])/1e6))*2);
      target_up[j]<-as.numeric(potential[j])+(as.numeric(tolmz)*((as.numeric(potential[j])/1e6)*2));
    }else{
      target_down[j]<-as.numeric(potential[j])-((as.numeric(tolmz)/1000)*2);
      target_up[j]<-as.numeric(potential[j])+((as.numeric(tolmz)/1000)*2);
    }
  } 
  
  close(prog)
  
  prog<-winProgressBar("Linkage Analysis",min=0,max=length(potential));
  
  for(j in 1:length(potential)){  
    setWinProgressBar(prog, j, title = NULL, label = NULL) 
    for(k in 1:length(TP.peaks[,1])){
      if(as.numeric(TP.peaks[k,1]) > as.numeric(target_up[j])){
        next
      }else{ if(
        ( as.numeric(TP.peaks[k,1]) >= as.numeric(target_down[j])) &&
          ( as.numeric(TP.peaks[k,1]) <= as.numeric(target_up[j]))
      ){
        p <- back[[j]][1]
        links[j,2] <- back[[j]][1]
        links[j,3] <- row.names(TP.peaks[k,])
        links[j,4] <- back[[j]][2]
        links[j,5] <- par.peaks[paste(p), "mz"]
        links[j,6] <- potential[j]
        links[j,7] <- TP.peaks[k,"mz"]
        links[j,8] <- par.peaks[paste(p), "rt"]
        links[j,9] <- TP.peaks[k, "rt"]
        links[j,10] <- par.peaks[paste(p), "into"]
        links[j,11] <- TP.peaks[k, "into"]
      } else {next}
      }
    }  
  }
  
  close(prog)
  links <- subset(links, links[,2] != "NA")
  print(c("Number of links made", nrow(links)))
  print(table(links$transformation)) 

  if(RTfilter == TRUE){
    links <- subset(links, links[,8] > links[,9])
    print(c("Number of links after RT filter", nrow(links)))
    print(table(links$transformation))
  }

  if(nrow(links)>1){
    links[,1] <- c(1:nrow(links))
  }else{
    stop("No links found")
  }
  
  save(links, file=paste(filename, ".rds", sep = ""))
  write.csv(links, file = paste(filename, ".csv", sep = ""))
  print("Linkage completed")
  
 # load("shortfeatlist.rds")
 # sample.info <- read.csv("sampleinfo.csv", header = TRUE)
#  sample.info <- sample.info[with(sample.info, order(rn)), ]
#  rownames(sample.info) <- c(1:nrow(sample.info))
 # load("data.rds")
  
}
  ############################################Link.data table############################################
linkdata <- function(links, data, filename = "links"){
  
  print(paste("Using working directory ", getwd(), sep=""))
  dir <- list.files(path = getwd(), all.files = TRUE)
  if("shortfeatlist.rds" %in% dir == FALSE){stop("Add shortfeatlist.rds to directory!")}
  if("sampleinfo.csv" %in% dir == FALSE){stop("Add sampleinfo.csv to directory!")}
#  if("data.rds" %in% dir == FALSE){stop("Add data.rds to directory!")} 


 link.data <- list()
  
  for(i in 1:nrow(links)){
  #  cat("Percent done with Linkage data table ", paste(round(i/nrow(links)*100)), "%", sep = "")
    a <- links[i,"Parent Peak ID"]  #short.feat peak ID of parent
    b <- links[i,"TP peak ID"]  #short.feat peak ID of TP
    
    pt <- data.frame(short.feat[[as.numeric(a)]])
    pt$PredictedClass <- rep("par", nrow(pt))
    pt$SampleClass <- rep(NA, nrow(pt))
    pt$SampleType <- rep(NA, nrow(pt))
   # pt$IsoPresence <- rep(NA, nrow(pt))
    
    sample.info <- read.csv("sampleinfo.csv", header = TRUE)
  # sample.info$SampleNumber <- gsub(c("20140924_online_pos_", "20141010_online_pos_"), "", sample.info$Filename)
    sample.info$SampleNumber <- sample.info$rn
    
    for(j in 1:nrow(pt)){
      c <- as.numeric(paste(pt[j,"sample"])) #sample with detection in pt table, starting from top
      d <- subset(sample.info, sample.info$rn == c) #sample info for current pt detection
      pt[j, "SampleClass"] <- paste(d$SampleClass)
      pt[j, "SampleType"] <- paste(d$SampleType)
      
      z <- pt[j,"ID"] #peak ID for current pt detection
      w <- as.numeric(rownames(d)) #need this because the rowname corresponds to the number of the data in the "data" list
      y <- subset(data[[as.numeric(w)]], data[[as.numeric(w)]]$peak.ID == z)
      
    #  pt[j, "IsoPresence"] <- paste(y$isotope.s.)   
    }
    
    pt$AssignCorr <- rep(NA, nrow(pt))
    
    for(j in 1:nrow(pt)){
      if(pt[j,"PredictedClass"]=="par" && pt[j,"SampleClass"]=="AV"){
        pt[j,"AssignCorr"] <- paste("TRUE")
      } else {
        pt[j,"AssignCorr"] <- paste("FALSE")
      }
    }
    
    tp <- data.frame(short.feat[[as.numeric(b)]])
    tp$PredictedClass <- rep("TP", nrow(tp))
    tp$SampleClass <- rep(NA, nrow(tp))
    tp$SampleType <- rep(NA, nrow(tp))
  #  tp$IsoPresence <- rep(NA, nrow(tp))
    
    for(k in 1:nrow(tp)){
      e <- as.numeric(paste(tp[k,"sample"]))
      f <- subset(sample.info, sample.info$rn == e)
      tp[k, "SampleClass"] <- paste(f$SampleClass)
      tp[k, "SampleType"] <- paste(f$SampleType)
      
      x <- short.feat[[as.numeric(b)]][k,"ID"]
      u <- as.numeric(rownames(f)) #need this because the rowname corresponds to the number of the data in the "data" list
      w <- subset(data[[as.numeric(u)]], data[[as.numeric(u)]]$peak.ID == x)
      
  #    tp[k, "IsoPresence"] <- paste(w$isotope.s.)               
    }
    
    tp$AssignCorr <- rep(NA, nrow(tp))  
    for(k in 1:nrow(tp)){
      if(tp[k,"SampleClass"] == "AB" | tp[k,"SampleClass"] =="noISTD"){
        tp[k,"AssignCorr"] <- paste("TRUE")
      }else{
        tp[k,"AssignCorr"] <- paste("FALSE")
      }      
    }
    link.data[[i]] <- rbind(pt, tp) 
  }
  
  Rfile <- paste(filename, ".rds", sep = "")
  save(link.data, file=Rfile)
  load(Rfile)
  
  ####################################
#  iso.pt <- list()
#  for(i in 1:length(link.data)){
#    iso.pt[[i]] <- subset(link.data[[i]], link.data[[i]]$PredictedClass == "par")
#    iso.pt[[i]] <- table(iso.pt[[i]]$IsoPresence)
#  }
  
#  save(iso.pt, file = "isopt.rds")
  
#  iso.tp <- list()
#  for(i in 1:length(link.data)){
#    iso.tp[[i]] <- subset(link.data[[i]], link.data[[i]]$PredictedClass == "TP")
#    iso.tp[[i]] <- table(iso.tp[[i]]$IsoPresence)
#  }
  
#  save(iso.tp, file = "isotp.rds")
  
  ###################################
  links$ParentCount
  links$ParentStdev
  links$TPCount
  links$TPStdev
  links$sumCount
#  links$ParentIsoPercent
#  links$TPIsoPercent
#  links$sumIsoPercent
  links$CorrectClass
  
  
  for(i in 1:nrow(links)){
    cat("\t Percent done with prioritization", paste(round(i/nrow(links)*100)))
    
    c <- links[i,2]  #short.feat peak ID of parent
    d <- links[i,3]  #short.feat peak ID of TP
    
    pt <- short.feat[[as.numeric(c)]]
    links[i,"ParentStdev"] <- sd(pt[,"m/z"])
    tp <- short.feat[[as.numeric(d)]]
    links[i,"TPStdev"] <- sd(tp[,"m/z"])
    
    
    a <- table(link.data[[i]]$PredictedClass)
    links[i,"ParentCount"] <- as.integer(a[names(a)=="par"])
    links[i,"TPCount"] <- as.integer(a[names(a)=="TP"])
    links[i,"sumCount"] <- links[i,"ParentCount"] + links[i,"TPCount"]
    
    
 #   q <- sum(iso.pt[[i]][names(iso.pt[[i]]) != "none"])
 #   r <- sum(iso.pt[[i]])
 #   links[i,"ParentIsoPercent"] <- as.numeric(paste(q/r * 100))
        
 #   s <- sum(iso.tp[[i]][names(iso.tp[[i]]) != "none"])
 #   t <- sum(iso.tp[[i]])
 #   links[i,"TPIsoPercen"] <- as.numeric(paste(s/t * 100))
    
 #   links[i,"sumIsoPercent"] <- paste(as.numeric(paste(q/r * 100)) + as.numeric(paste(s/t * 100)))
    
    u <- sum(link.data[[i]]$AssignCorr == TRUE)
    v <- nrow(link.data[[i]])
    links[i,"CorrectClass"] <- paste(u/v * 100)
    
  }
  
#  print(head(rownames(links)))
  names(link.data) <- rownames(links)
  
 # links <- links[with(links, order(sumIsoPercent, decreasing = TRUE)), ]  
  
  print(head(links))
  print(table(links$transformation))
  write.csv(links, file=paste(filename, ".csv", sep = ""))
  save(links, file=paste(filename, ".rds", sep = ""))
  load(paste(filename, ".rds", sep = ""))
  print(paste("Saved links to ", getwd(), sep=""))
  
}

#####   function for calculating p-values of each feature     ########################
linksPvalues <- function(links, PCA.peaks, sample.info, filename = "linksPvalues"){

  print(paste("Using working directory ", getwd(), sep=""))
 # dir <- list.files(path = getwd(), all.files = TRUE)
 # if("PCApeaks.rds" %in% dir == FALSE){stop("Add PCApeaks.rds to directory!")}
#  if("sampleinfo.csv" %in% dir == FALSE){stop("Add sampleinfo.csv to directory!")}
  #if("links.csv" %in% dir == FALSE){stop("Add links.csv to directory!")}

linkPvaluesPT <- c()
linkPvaluesTP <- c()

AB <- subset(sample.info, SampleClass == "effluent")
AB <- which(rownames(PCA.peaks) %in% AB$rn)
AV <- subset(sample.info, SampleClass == "influent")
AV <- which(rownames(PCA.peaks) %in% AV$rn)

if((length(AB)+length(AV) == nrow(PCA.peaks)) == FALSE) {stop("Check sample.info table samples vs. samples in PCA.peaks table")}

  for(i in 1:nrow(links)){
 #   cat("Percent done with Linkage data table ", paste(round(i/nrow(links)*100)), "%", sep = "")
    a <- links[i,"Parent Peak ID"]  #short.feat peak ID of parent

    a <- which(colnames(PCA.peaks) == paste(a))
    
    x <- PCA.peaks[,as.numeric(paste(a))]   

  #  y <- as.vector(subset(x, names(x) %in% AB$rn))
    y <- x[c(AB)]
  #  z <- as.vector(subset(x, names(x) %in% AV$rn))
    z <- x[c(AV)]
  
    g <- factor(c(rep("AV", length(z)), rep("AB", length(y))))
    v <- c(z,y)
    p <- t.test(v~g)
#print(c("p-value parent link", i, " ", p$p.value))

linkPvaluesPT[i] <- p$p.value


 b <- links[i,"TP peak ID"]  #short.feat peak ID of TP
    
 b <- which(colnames(PCA.peaks) == paste(b))

         k <- PCA.peaks[,as.numeric(paste(b))]   

  #  l <- as.vector(subset(k, names(k) %in% AB$rn))
  l <- k[c(AB)]
  #  m <- as.vector(subset(k, names(k) %in% AV$rn))
  m <- k[c(AV)]

    n <- factor(c(rep("AV", length(m)), rep("AB", length(l))))
    q <- c(m,l)
    t <- t.test(q~n)
#print(c("p-value TP link", i, " ", t$p.value))

linkPvaluesTP[i] <- t$p.value
 
    }

tab <- cbind(c(1:nrow(links)), linkPvaluesPT, linkPvaluesTP)
colnames(tab) <- c("link", "parent p-value", "TP p-value")
write.csv(tab, file = paste(filename, ".csv", sep = ""))
save(tab, file = paste(filename, "_ptable.rds", sep = ""))
#return(tab)

}

################################################################################################################################
################################################################################################################################

# match.par <- c()
# 
# for(i in 1:length(links[,1])){
#   
#   props <- links[i,c("parent m/z", "Parent RT")]
#   
#   for(j in 1:length(NTdata[,1])){
#     
#     NT <- NTdata[j,c("mz", "rt")]
#     
#     
#     if(identical(as.numeric(round(NTdata[j,1], 4)), as.numeric(round(links[i,1], 4))) == TRUE
#        && identical(as.numeric(round(NTdata[j,2], 5)), as.numeric(round(links[i,2], 5))) == TRUE){
#       
#       match.par[i] <- as.numeric(rownames(NTdata[j,]))
#       
#     }
#     
#     
#   }
# }


################################################################################################################################
####NOTE:: In this function the short.feat and the sample.info table MUST be in the working directory######

#linkages(act, par.peaks, TP.peaks, tolmz = 5, ppm = TRUE, RTfilter = FALSE, filename = "linksAll3mgL_deamin_23022016")
#linkdata(links, data, filename = "linksAlldata2mgL_29012016")
#linksPvalues(links, PCA.peaks, sample.info, filename = "linksAll5mgL_deamin_Pvalues_23022016")

# load("links.rds")
# load("linkdata.rds")
# load("ptable.rds")
# load("isopt.rds")
# load("isotp.rds")

