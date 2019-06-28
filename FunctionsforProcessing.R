

#### Functions for HRMS(MS)
### Extensions to either MSnbase or enviMass

### MSnbase

mzWindow <- function(mz, mztol){
  max <- mz + ppm(mz, mztol, p = TRUE)
  min <- mz - ppm(mz, mztol, p = TRUE)
  wind <- max-min
  return(wind)
}

rtbounds <- function(exp.rt, rtwind){
  
  rtmin <- exp.rt - rtwind
  rtmax <- exp.rt + rtwind
  return(c(rtmin, rtmax))
}

findMS2 <- function(hd.ms1, hd.ms2, precursor.mz, mztol){
  
  mzlimits <- ppm(precursor.mz, mztol, l = TRUE)
  hd.ms1 <- as.data.frame(hd.ms1)
  hd.ms2 <- as.data.frame(hd.ms2)
  
  ms2.scans <- which(hd.ms2$precursorMZ < mzlimits[1] & hd.ms2$precursorMZ > mzlimits[2])
  if(length(ms2.scans)==0){return(NULL)}
  
  eic.int <- as.numeric(as.character(hd.ms2$precursorIntensity[ms2.scans]))
  ms2.rt <- as.numeric(as.character(hd.ms2$retentionTime[ms2.scans]))
  scan.nm <- hd.ms2$acquisitionNum[ms2.scans]
  return(as.data.frame(cbind(scan.nm  = scan.nm, ms2.rt = ms2.rt, eic.int = eic.int)))
  
}

selectIon <- function(all.file, polarity, msLevel){
  
  if(polarity == "positive"){
    pos.data <- which(names(all.file) == "1")
    if(length(pos.data)==0) {return(NULL)}
    else{
      ms.file <- which(names(all.file[[pos.data]]) == msLevel)
      plot.file <- all.file[[pos.data]][[ms.file]]
      return(plot.file)
      
    }
  }
  
  if(polarity == "negative"){
    neg.data <- which(names(all.file) == "0")
    if(length(neg.data)==0) {return(NULL)}
    else{
      ms.file <- which(names(all.file[[neg.data]]) == msLevel)
      plot.file <- all.file[[neg.data]][[ms.file]]
      return(plot.file)
      
    } 
  }
  
}

readmzXML <- function(file.path, mode){
  
  # ms1 <- readMSData(file.path, centroided. = TRUE, msLevel. = 1, mode = mode)
  # h.ms1 <- header(ms1)
  # 
  # ms2 <- readMSData(file.path, centroided. = TRUE, msLevel. = 2, mode = mode)
  # h.ms2 <- header(ms2)
  # 
  # return(list(list(ms1, h.ms1), list(ms2, h.ms2)))
  
  all.file <- readMSData(file.path, centroided. = TRUE, mode = mode)
  split.file <- split(all.file, polarity(all.file))
  split.file <- lapply(split.file, function(x) split(x, msLevel(x)))
  
  return(split.file)
  
}

t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
  
}

myXIC <- function(eic, ms2.data, exp.rt){
  
  plot(eic)
  abline(v = exp.rt, col = "blue")
  
  if(nrow(ms2.data) < 1 || is.null(ms2.data))
    legend("topleft", legend = paste("No MS2 scans"))
  else
    points(ms2.data[,2:3], col = t_col("red", 50), pch = 16, type = "h", lty = 2)
}

myHeadtoTail <- function(sample.spectrum, ref.spectrum, y, z, t, b, id.file, id.ref){
  
  my.spectra <- OrgMSSim(as.data.frame(sample.spectrum), as.data.frame(ref.spectrum), y = y, z = z, t = t, b = b)
  match <- subset(my.spectra[[2]], my.spectra[[2]]$match > 0)
  plot.new()
  plot.window(xlim = c(0,max(my.spectra[[2]]$mz)+50), ylim = c(-1,1))
  points(x = my.spectra[[2]]$mz, y = (my.spectra[[2]][,4])/max(my.spectra[[2]][,4]), type = "h", col = "grey")
  points(x = my.spectra[[2]]$mz, y = -(my.spectra[[2]][,5])/max(my.spectra[[2]][,5]), type = "h", col = "grey")
  points(x = match$mz, y = (match[,4])/max(my.spectra[[2]][,4]), type = "h", col = "red")
  points(x = match$mz, y = -(match[,5])/max(my.spectra[[2]][,5]), type = "h", col = "blue")
  abline(h = 0)
  abline(h = input$MSMS_b, lty = 2, col = "grey")
  abline(h = -input$MSMS_b, lty = 2, col = "grey")
  axis(side = 1)
  axis(side = 2)
  title(main = c(id.file, "\n", id.ref), sub = paste("similarity score = ", round(my.spectra$score, 3), sep = ""),
        xlab = "Fragment m/z", ylab = "Relative intensity")
  box()
}



#### enviMass

#### blank/blind detection

blankDetect <- function(profiles, blank.code){
  
  bk.det <- unique(unlist(sapply(blank.code, function(x){which(profiles$peaks$sample == x)})))
  bk.profs <- unique(profiles$peaks$profileIDs[bk.det])
  profiles$index_prof[bk.profs,"in_blind?"] <- TRUE
  
}


#### RAMClustR and RMassScreening
### get component associated with profile and plot

### for checking RT of suspects. just adds TRUE/FALSE if measured.rt within window of expected.rt
RTfilter <- function(expected.rt, measured.rt, rttol){
  
  match <- c()
  if(measured.rt > expected.rt - rttol && measured.rt < expected.rt + rttol)
    match <- TRUE
  else
    match <- FALSE
  
  return(match)
}

### extracts the spectra with RMassBank of the desired profile
getComponent <- function(profileID, spectra){
  
  loc <- which(spectra@aggregated$profile == profileID)
  loc2 <- as.numeric(spectra@aggregated[loc,]$cpdID)
  data <- getData(spectra@spectra[[loc2]]@children[[1]])
  return(data)
}

### plots MS1 or MS2
plotComponent <- function(component.data, ...){
  plot.new()
  plot.window(ylim = c(0, max(component.data$intensity)*1.2), xlim = c(0,max(component.data$mz)*1.5), xlab = "m/z", ylab = "Absolute intensity")
  points(component.data$mz, component.data$intensity, type = "h", lwd = 2)
  text(component.data$mz, component.data$intensity, labels = round(component.data$mz,4), pos = 3, offset = 1.5)
  text(component.data$mz, component.data$intensity, labels = component.data$profile, pos = 3)
  box()
  axis(1)
  axis(2)
  title("Plot RAMClust Component", xlab="m/z",
        ylab="Absolute intensity")
}


