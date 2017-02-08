

#' Starting a new profile list to fill
#'
#' @param picked load in your data. Structure of data is a list, where each element of the list is one peak-picked file, ie. enviPick output
#'
#' @return empty profile list
#' @export
#'
#' @examples


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