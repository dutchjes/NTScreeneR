

#### Simple suspect screening in R
#### written by M. Stravs

ppm <- function(mz, ppmlimit) (1e-6 * mz * ppmlimit)

suspScreen <- function(measured, suspects, ppmMax = 5, rtTol = 1){
  
  if("mz" %in% colnames(suspects) == FALSE) {stop("Need a 'mz' column in suspects table!")} 
  if("mz" %in% colnames(measured) == FALSE) {stop("Need a 'mz' column in measured table!")} 
  if("rt" %in% colnames(suspects) == FALSE) {stop("Need a 'rt' column in suspects table!")} 
  if("rt" %in% colnames(measured) == FALSE) {stop("Need a 'rt' column in measured table!")} 
  
  peaklist <- measured
  hits <- lapply(1:nrow(peaklist), function(n)
  {
    mass <- peaklist[n,"mz"]
    peaklist.hit <- suspects[abs(suspects$mz - mass) < ppm(mass, ppmMax),,drop=FALSE]
 #   peaklist.hit <- peaklist.hit[,c("mz", "rt")]
    
    peaklist.hit$dppm <- ((peaklist.hit$mz / mass)-1)*1e6
    peaklist.hit$mz_measured <- rep(peaklist[n, "mz"], nrow(peaklist.hit))
    peaklist.hit$rt_measured <- rep(peaklist[n, "rt"], nrow(peaklist.hit))
   # peaklist.hit$ID <- rep(peaklist[n, "name"], nrow(peaklist.hit))

    return(peaklist.hit)
  })
  
  hits.total <- do.call(rbind, hits)
  
#  if(rtTol == FALSE){return(hits.total)}
  
  hits.total <- hits.total[which(hits.total$rt < hits.total$rt_measured + rtTol),] # && hits.total$rt > hits.total$RT - 1
  hits.total <- hits.total[which(hits.total$rt > hits.total$rt_measured - rtTol),]
  
  return(unique(hits.total))
}

