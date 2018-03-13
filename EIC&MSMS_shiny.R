


#load packages for processing
library(RMassBank)
library(mzR)
library(Hmisc)
library(enviPat)
library(shiny)
library(MSnbase)
library(ggplot2)
library(devtools)
#install_github("dutchjes/MSMSsim")
library(MSMSsim)

ui <- fluidPage(
  
  titlePanel("Extracting EICs and MS/MS scans"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "P_data", label = "Choose file with precursor masses and rts:"), ## user upload csv with mz and rt
      selectInput(inputId = "masscol", label = "Column with precursor masses", choices = c()),
      selectInput(inputId = "rtcol", label = "Column with RTs", choices = c()),
      selectInput(inputId = "rtunit", label = "Retention time units", choices = c("minutes", "seconds")),
      checkboxInput(inputId = "useRT", "Use retention time for feature extraction"),
      selectInput(inputId = "IDcol", label = "Column with ID", choices = c()),
      selectInput(inputId = "ioncol", label = "Column with ionization mode", choices = c()),
      selectInput(inputId = "ionization", label = "Which Ionization mode?", choices = c()),
      numericInput(inputId = "masstol", label = "ppm tolerence", value = c(10)),
      numericInput(inputId = "rtwind", label = "rt window (sec)", value = c(60)), ## user selected
      textInput(inputId = "filedir", label = "Directory with mzXML files:"), ## path to RAW files
      textOutput(outputId = "status")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", tableOutput(outputId = "userP_data")),
        
        tabPanel("TIC", 
                 selectInput(inputId = "RAWFileTIC", label = "Current mzXML file:", choices = c(), multiple = FALSE), ##  filled with names in filedir
                 actionButton(inputId = "createTIC", label = "Create TIC"),
                 plotOutput(outputId = "TIC"),
                 actionButton(inputId = "create2d", label = "Create 2d plot"),
                 plotOutput(outputId = "plot2d"),
                 selectInput(inputId = "coloring", label = "Color points according to:", choices = c("ionCount", "file", "peaks.count", "charge"))
                 ),
        
        tabPanel("EIC and MS2", 
                 selectInput(inputId = "feature", label = "feature to extract", choices = c()), ## filled with masses from P_data
                 selectInput(inputId = "RAWFileEIC", label = "Current mzXML file:", choices = c(), multiple = FALSE), ##  filled with names in filedir
                 actionButton(inputId = "createEIC", label = "Create EIC"),
                 textOutput(outputId = "Compound"),
                 plotOutput(outputId = "EIC"),
                 tableOutput(outputId = "availMSMSscans"),
                 selectInput(inputId = "selectedMSMS", label = "Which MSMS to plot?", choices = c()),
                 plotOutput(outputId = "MSMS")
                 ),
        
        tabPanel("EIC Output", 
                 checkboxGroupInput(inputId = "RAWFileOutput", label = "Files to process:", choices = c(), inline = TRUE), ##  filled with names in filedir
                 actionButton("selectall", "Select All Files"),
                 #checkboxInput(inputId = "pdf", label = "Save images to pdf?"), 
                 textInput(inputId = "outputdir", label = "Directory to save pdf:"),
                 actionButton(inputId = "process", label = "Save to pdf")
                 ),
        
        tabPanel("MS2 Output",
                 selectInput(inputId = "MS2filecol", label = "Column with mzXMLs to extract MS2", choices = c()),
                 selectInput(inputId = "reffilecol", label = "Column with mzXMLs to extract reference MS2", choices = c()),
                 checkboxInput(inputId = "htot", label = "Plot head to tail between File and Reference"),
                 radioButtons(inputId = "whichMSMS", label = "Which MS2 scans to plot", choices = c("only most intense precursor in selected RT range",
                                                                                                    "all MS2 scans in selected RT range")),
                 numericInput(inputId = "MSMS_y", label = "m/z weighting factor (y)", value = 0),
                 numericInput(inputId = "MSMS_z", label = "intensity weighting factor (z)", value = 1),
                 numericInput(inputId = "MSMS_t", label = "m/z tolerence for merging fragments (Da)", value = 0.005),
                 numericInput(inputId = "MSMS_b", label = "relative intensity cutoff", value = 0.01),
                 textInput(inputId = "MS2outputdir", label = "Directory to save MS2 pdfs:"),
                 actionButton(inputId = "processMS2", label = "Save MS2s to pdf")
                 )
      )
    )
  )
)

server <- function(input, output, session){
  
  P_data <- eventReactive(input$P_data, {
    inFile <- input$P_data
    newFile <- read.csv(inFile$datapath, header = TRUE)
    
  })
  
  # mzWindow <- function(mz, mztol){
  #   max <- mz + ppm(mz, mztol, p = TRUE)
  #   min <- mz - ppm(mz, mztol, p = TRUE)
  #   wind <- max-min
  #   return(wind)
  # }
  
  rtbounds <- function(exp.rt, rtwind){
    
    rtmin <- exp.rt - rtwind
    rtmax <- exp.rt + rtwind
    return(c(rtmin, rtmax))
  }
  
  findMS2 <- function(hd.ms1, hd.ms2, precursor.mz, mztol){
    
    mzlimits <- ppm(precursor.mz, mztol, l = TRUE)
    ms2.scans <- which(hd.ms2$precursorMZ < mzlimits[1] & hd.ms2$precursorMZ > mzlimits[2])
    if(length(ms2.scans)==0){return(NULL)}

    eic.int <- as.numeric(as.character(hd.ms2$precursorIntensity[ms2.scans]))
    ms2.rt <- as.numeric(as.character(hd.ms2$retentionTime[ms2.scans]))
    scan.nm <- hd.ms2$acquisitionNum[ms2.scans]
    return(as.data.frame(cbind(scan.nm  = scan.nm, ms2.rt = ms2.rt, eic.int = eic.int)))
    
  }
  
  readmzXML <- function(file.path, mode){
    
    ms1 <- readMSData(file.path, centroided. = TRUE, msLevel. = 1, mode = mode)
    h.ms1 <- header(ms1)
    
    ms2 <- readMSData(file.path, centroided. = TRUE, msLevel. = 2, mode = mode)
    h.ms2 <- header(ms2)
    
    return(list(list(ms1, h.ms1), list(ms2, h.ms2)))
    
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
  
  output$userP_data <- renderTable({
    P_data()
  })
  
  rv <- reactiveValues(masses = c(), ## masses to populate mass dropdown
                       rts_min = c(), ## rts to populate rt dropdown
                       rts_sec = c(), ## rts for chromatograph function
                       names = c(), ## names of NTs
                       ESImodes = c(),
                       col = c(),
                       features = c(),
                       filepath = c(),
                       filename = c(),
                       f = c(),
                       r = c(),
                       allcalc = c(),
                       msms = list(),
                       ms2.data = c(),
                       ms2.files = c(),
                       ref.files = c(),
                       spec.file = list(),
                       spec.ref = list())
  
  
  observeEvent(input$P_data,{
    rv$col <- colnames(P_data())
    
    guess.mass <- c(which(grepl("Mass", colnames(P_data()))), which(grepl("mz", colnames(P_data()))))
    updateSelectInput(session, "masscol", choices = rv$col, selected = rv$col[guess.mass[1]])
    
    guess.rt <- c(which(grepl("RT", colnames(P_data()))), which(grepl("retention", colnames(P_data()))))
    updateSelectInput(session, "rtcol", choices = rv$col, selected = rv$col[guess.rt[1]])
    
    guess.name <- c(which(grepl("ID", colnames(P_data()))), which(grepl("name", colnames(P_data()))))
    updateSelectInput(session, "IDcol", choices = rv$col, selected = rv$col[guess.name[1]])
    
    guess.mode <- c(which(grepl("Polarity", colnames(P_data()))), which(grepl("ESI", colnames(P_data()))))
    updateSelectInput(session, "ioncol", choices = rv$col, selected = rv$col[guess.mode[1]])
    
    guess.file <- c(which(grepl("Sample", colnames(P_data()))), which(grepl("File", colnames(P_data()))))
    updateSelectInput(session, "MS2filecol", choices = rv$col, selected = rv$col[guess.file[1]])
    
    guess.ref <- c(which(grepl("Reference", colnames(P_data()))), which(grepl("ref", colnames(P_data()))))
    updateSelectInput(session, "reffilecol", choices = rv$col, selected = rv$col[guess.ref[1]])
    
  })
 
  
#  observeEvent(c(input$ionization, input$ioncol),{
#       pol <- which(colnames(P_data()) == rv$ioncol)
#       P_data() <- subset(P_data(), P_data()[,pol] == input$ionization)
#  }) 
  
  
  observeEvent(c(input$masscol, input$rtcol, input$IDcol, input$ioncol, input$rtunit), {

    tuke <- which(colnames(P_data()) == input$ioncol)
    rv$ESImodes <- unique(P_data()[,tuke])
    updateSelectInput(session, "ionization", choices = rv$ESImodes)

    take <- which(colnames(P_data()) == input$masscol)
    rv$masses <- as.numeric(as.character(P_data()[,take]))

    tike <- which(colnames(P_data()) == input$rtcol)
    if(input$rtunit == "seconds"){
      rv$rts_min <- as.numeric(as.character(P_data()[,tike])) / 60
      rv$rts_sec <- as.numeric(as.character(P_data()[,tike]))
    }
    if(input$rtunit == "minutes"){
      rv$rts_min <- as.numeric(as.character(P_data()[,tike]))
      rv$rts_sec <- as.numeric(as.character(P_data()[,tike])) * 60
    }
    
    rv$features <- paste(rv$masses, rv$rts_min, sep = "_")
    updateSelectInput(session, "feature", choices = rv$features)

    toke <- which(colnames(P_data()) == input$IDcol)
    rv$names <- P_data()[,toke]
    
    # tuke <- which(colnames(P_data()) == input$ioncol)
    # rv$ESImodes <- unique(P_data()[,tuke])

  })
  

  observeEvent(input$filedir,{
    
    rv$filepath <- list.files(input$filedir, pattern = "mzXML", full.names = TRUE, recursive = FALSE)
    rv$filename <- list.files(input$filedir, pattern = "mzXML", full.names = FALSE, recursive = FALSE)
    
    updateSelectInput(session, "RAWFileTIC", choices = rv$filename)
    updateSelectInput(session, "RAWFileEIC", choices = rv$filename)
    updateCheckboxGroupInput(session, "RAWFileOutput", choices = rv$filename)
    
  })
  
  ##########################
  ### Plotting the TIC in Shiny
  
  observeEvent(input$createTIC,{

      process.file1 <- which(rv$filename == input$RAWFileTIC)
      
      if(length(process.file1) == 0){output$status <- renderText({"Waiting for mzXML files..."})}
      
      else{
        
        output$status <- renderText({"Plotting TIC"})
        rv$f <- try(readmzXML(rv$filepath[process.file1], mode = "onDisk"))

        output$status <- renderText({"Read in file"})
        output$TIC <- renderPlot({
          plot(chromatogram(rv$f[[1]][[1]]))
        })
        
        output$status <- renderText({"Finished"})

        observeEvent(input$create2d,{
          spec <- tryCatch(readMSData(rv$filepath[process.file1], msLevel. = 1, centroided. = TRUE))
          
          output$status <- renderText({"Plotting 2d"})
          
          output$plot2d <- renderPlot({
            plot2d(spec, z  = input$coloring)
          })
          
          output$status <- renderText({"Finished 2d"})

          #  plot(MSmap(rv$ftic, lowMz = min(rv$masses), highMz = max(rv$masses), resMz = ppm(max(rv$masses), input$masstol, p = TRUE)), aspect = 1, allTicks = FALSE)
        })
      }
  });
  
  ##########################
  ### Plotting the EIC in Shiny
  
  observeEvent(input$createEIC,{

    process.file2 <- which(rv$filename == input$RAWFileEIC)
    process.feature <- which(rv$features == input$feature) ### will need to fix the problem of multiple features which are the same
    
    if(length(process.file2) == 0){output$status <- renderText({"Waiting for mzXML files..."})}
    if(length(process.feature) == 0){output$status <- renderText({"Waiting for suspect list..."})}
    

    if(length(process.file2) == 1){
      
      output$status <- renderText({"Ready to plot EIC..."})
      rv$f <- readmzXML(rv$filepath[process.file2], mode = "onDisk")

    if(input$useRT == TRUE){
      rt.lims <- rtbounds(rv$rts_sec[process.feature[1]], input$rtwind)
    }else{
      rt.lims <- c(min(rtime(rv$f[[1]][[1]])),max(rtime(rv$f[[1]][[1]])))
    }
      
      eic <- chromatogram(rv$f[[1]][[1]], rt = rt.lims, mz = ppm(rv$masses[process.feature[1]], input$masstol, l = TRUE))
      rv$ms2.data <- findMS2(rv$f[[1]][[2]], rv$f[[2]][[2]], rv$masses[process.feature[1]], input$masstol) 
      rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rt.lims[1] & rv$ms2.data[,2] < rt.lims[2])
    output$EIC <- renderPlot({myXIC(eic = eic, ms2.data = rv$ms2.data, exp.rt = rv$rts_sec[process.feature[1]])})

        if(!is.null(rv$ms2.data) & nrow(rv$ms2.data) > 0 ){
      updateSelectInput(session, "selectedMSMS", choices = rv$ms2.data[,1])
      output$availMSMSscans <- renderTable({rv$ms2.data})
      observeEvent(input$selectedMSMS, 
                   output$MSMS <- renderPlot({
                     plot(rv$f[[2]][[1]][[which(rv$f[[2]][[2]]$acquisitionNum == input$selectedMSMS)]], reporters = NULL, full = TRUE) + theme_bw()
                   })
      )
    }
    output$Compound <- renderText({paste(rv$names[process.feature[1]])})
    output$status <- renderText({"Finished"})
    }
  })


  observeEvent(c(input$RAWFileOutput,input$selectall), {

    if(input$selectall == 0){
      rv$allcalc <- input$RAWFileOutput

    }else if(input$selectall%%2 == 0){
      updateCheckboxGroupInput(session, "RAWFileOutput", choices = rv$filename)
      rv$allcalc <- input$RAWFileOutput

    }else{
      updateCheckboxGroupInput(session, "RAWFileOutput", choices = rv$filename, selected = rv$filename)
      rv$allcalc <- input$RAWFileOutput
    }
  })


  ##########################
  ### Output the MS/MS in pdfs
  
  observeEvent(c(input$MS2filecol, input$reffilecol),{
    
    get <- which(colnames(P_data()) == input$MS2filecol)
    rv$ms2.files <- P_data()[,get]
    
    git <- which(colnames(P_data()) == input$reffilecol)
    rv$ref.files <- P_data()[,git]
    
  })
  
  observeEvent(input$processMS2, {
    if(input$htot == TRUE){
      pdf(file = paste(input$MS2outputdir, "\\", "HeadtoTailMS2.pdf", sep = ""))
      
      for(i in 1:length(rv$ms2.files)){
        
        print(i)
        ### sample file
        if(grepl(".mzXML", rv$ms2.files[i]))
          process.file <- which(grepl(paste(rv$ms2.files[i]), rv$filename) > 0)
        else
          process.file <- which(paste(rv$ms2.files[i], ".mzXML", sep = "") == rv$filename)
        
        output$status <- renderText({paste("Working with file", rv$ms2.files[i], sep = " ")})
        rv$f <- readmzXML(file.path = rv$filepath[process.file], mode = "onDisk")
        
        if(input$useRT == TRUE){
          rt.lims <- rtbounds(rv$rts_sec[i], input$rtwind)
        }else{
          rt.lims <- c(min(rtime(rv$f[[1]][[1]])),max(rtime(rv$f[[1]][[1]])))
        }
        
        rv$ms2.data <- findMS2(rv$f[[1]][[2]], rv$f[[2]][[2]], rv$masses[i], input$masstol)
        rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rt.lims[1] & rv$ms2.data[,2] < rt.lims[2])
        
        if(is.null(rv$ms2.data)==TRUE || nrow(rv$ms2.data) < 1)
          rv$spec.file <- NA
        else{
          
          ## input for intense or all MS2
          if(input$whichMSMS == "only most intense precursor in selected RT range"){
            get.scan <- which.max(rv$ms2.data[,3])
            
            if(length(get.scan)<1)
              rv$spec.file <- NA
            else
              rv$spec.file <- rv$f[[2]][[1]][[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][get.scan])]]
            
            ce <- rv$f[[2]][[2]]$collisionEnergy[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][get.scan])]
            id.file <- paste(rv$features[i], "_NCE", ce, "_", formatC(max(as.data.frame(rv$spec.file)$i), format = "e", 1), "_", rv$ms2.files[i], sep = "")
          }
 
          if(input$whichMSMS=="all MS2 scans in selected RT range"){

              ce <- c()
              id.file <- c()
              rv$spec.file <- list(c())
              
              for(k in 1:nrow(rv$ms2.data)){
                
                rv$spec.file[[k]] <- rv$f[[2]][[1]][[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][k])]]
                ce[k] <- rv$f[[2]][[2]]$collisionEnergy[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][k])]
                id.file[k] <- paste(rv$features[i], "_NCE", ce[k], "_", formatC(max(as.data.frame(rv$spec.file[[k]])$i), format = "e", 1), "_", rv$ms2.files[i], sep = "")
              }
          }
        }
        
        ### reference file
        if(grepl(".mzXML", rv$ref.files[i]))
          ref.file <- which(grepl(paste(rv$ref.files[i]), rv$filename) > 0)
        else
          ref.file <- which(paste(rv$ref.files[i], ".mzXML", sep = "") == rv$filename)
        
        output$status <- renderText({paste("Working with file", rv$ref.files[i], sep = " ")})
        rv$r <- readmzXML(file.path = rv$filepath[ref.file], mode = "onDisk")
        
        if(input$useRT == TRUE){
          rt.lims <- rtbounds(rv$rts_sec[i], input$rtwind)
        }else{
          rt.lims <- c(min(rtime(rv$r[[1]][[1]])),max(rtime(rv$r[[1]][[1]])))
        }
        
        rv$ms2.data <- findMS2(rv$r[[1]][[2]], rv$r[[2]][[2]], rv$masses[i], input$masstol)
        rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rt.lims[1] & rv$ms2.data[,2] < rt.lims[2])
        
        if(is.null(rv$ms2.data)==TRUE || nrow(rv$ms2.data) < 1)
          rv$spec.ref <- NA
        else{
          
          ## input for intense or all MS2
          if(input$whichMSMS == "only most intense precursor in selected RT range"){
            get.scan <- which.max(rv$ms2.data[,3])
          
            if(length(get.scan)<1)
             rv$spec.ref <- NA
            else
             rv$spec.ref <- rv$r[[2]][[1]][[which(rv$r[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][get.scan])]]
            
          ce <- rv$r[[2]][[2]]$collisionEnergy[which(rv$r[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][get.scan])]
          id.ref <- paste(rv$features[i], "_NCE", ce, "_", formatC(max(as.data.frame(rv$spec.ref)$i), format = "e", 1), "_", rv$ref.files[i], sep = "")
          }
          
          if(input$whichMSMS=="all MS2 scans in selected RT range"){
            ce <- c()
            id.ref <- c()
            rv$spec.ref <- list(c())
            
            for(k in 1:nrow(rv$ms2.data)){
              
              rv$spec.ref[[k]] <- rv$r[[2]][[1]][[which(rv$r[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][k])]]
              ce[k] <- rv$r[[2]][[2]]$collisionEnergy[which(rv$r[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][k])]
              id.ref[k] <- paste(rv$features[i], "_NCE", ce[k], "_", formatC(max(as.data.frame(rv$spec.ref[[k]])$i), format = "e", 1), "_", rv$ref.files[i], sep = "")
            }
          }
        }
        
        if(!sum(!is.na(rv$spec.file) & !is.na(rv$spec.ref))==0){
          if(input$whichMSMS == "only most intense precursor in selected RT range"){
            
            myHeadtoTail(sample.spectrum = as.data.frame(rv$spec.file), ref.spectrum = as.data.frame(rv$spec.ref),
                         y = input$MSMS_y, z = input$MSMS_z, t = input$MSMS_t, b = input$MSMS_b, id.file = id.file, id.ref = id.ref)
            
            sink(paste(input$MS2outputdir,"\\", id.file, ".txt", sep = ""))
            print(rv$spec.file)
            print(as.data.frame(rv$spec.file))
            print(rv$spec.ref)
            print(as.data.frame(rv$spec.ref))
            sink()
          }
          
          if(input$whichMSMS == "all MS2 scans in selected RT range"){
            
            sink(paste(input$MS2outputdir,"\\AllMS2scans.txt", sep = ""))
            
            for(k in 1:length(rv$spec.file)){
              
              for(l in 1:length(rv$spec.ref)){
                
                myHeadtoTail(sample.spectrum = as.data.frame(rv$spec.file[[k]]), ref.spectrum = as.data.frame(rv$spec.ref[[l]]),
                             y = input$MSMS_y, z = input$MSMS_z, t = input$MSMS_t, b = input$MSMS_b, 
                             id.file = id.file[k], id.ref = id.ref[l])
                
                print(rv$spec.file[[k]])
                print(as.data.frame(rv$spec.file[[k]]))
                print(rv$spec.ref[[l]])
                print(as.data.frame(rv$spec.ref[[l]]))
              }
            }
            sink()
          }
        }
        else
          next
      }
      
      output$status <- renderText({"Finished building MS2 pdfs"})
      dev.off()
      
    }else{
      
      pdf(file = paste(input$MS2outputdir, "\\", "ExtractedMS2s.pdf", sep = ""))
      
      for(i in 1:length(rv$ms2.files)){
        print(i)
        if(grepl(".mzXML", rv$ref.files[i]))
          process.file <- which(grepl(paste(rv$ms2.files[i]), rv$filename) > 0)
        else
          process.file <- which(paste(rv$ms2.files[i], ".mzXML", sep = "") == rv$filename)
        
        output$status <- renderText({paste("Working with file", rv$ms2.files[i], sep = " ")})
        rv$f <- readmzXML(rv$filepath[process.file], mode = "onDisk")
        
        if(input$useRT == TRUE){
          rt.lims <- rtbounds(rv$rts_sec[process.feature], input$rtwind)
        }else{
          rt.lims <- c(min(rtime(rv$f[[1]][[1]])),max(rtime(rv$f[[1]][[1]])))
        }
        
        eic <- chromatogram(rv$f[[1]][[1]], rt = rt.lims, mz = ppm(rv$masses[i], input$masstol, l = TRUE))
        rv$ms2.data <- findMS2(rv$f[[1]][[2]], rv$f[[2]][[2]], rv$masses[i], input$masstol)
        if(is.null(rv$ms2.data)[1]==TRUE)
          next
        else{
          rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rt.lims[1] & rv$ms2.data[,2] < rt.lims[2])
          
          if(input$whichMSMS == "only most intense precursor in selected RT range"){
            
          get.scan <- which.max(rv$ms2.data[,3])
          if(length(get.scan)<1)
            next
          else
            plot(rv$f[[2]][[1]][[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][get.scan])]], reporters = NULL, full = TRUE) + theme_bw()
          }
          
          if(input$whichMSMS == "all MS2 scans in selected RT range"){
              ce <- c()
              id.ref <- c()
              
              for(k in 1:nrow(rv$ms2.data)){
                plot(rv$f[[2]][[1]][[which(rv$f[[2]][[2]]$acquisitionNum == rv$ms2.data[,1][k])]], reporters = NULL, full = TRUE) + theme_bw()
              }
            }
          }
        }
      };
      
      dev.off()
      output$status <- renderText({"Finished building MS2 pdfs"})
  })
  
  
  ##################################
  ### EIC output tab
  ### add compound name to plot
  
  observeEvent(input$process,{

    for(i in 1:length(rv$allcalc)){

      calcfile <- which(rv$filename == rv$allcalc[i])

      output$status <- renderText({paste("Working on file ", rv$allcalc[i], sep = "")})

      pdf(file = paste(input$outputdir, "\\", rv$filename[calcfile], "_PlottedEICs.pdf", sep = ""))
      plotting.file <- readmzXML(rv$filepath[calcfile], mode = "onDisk")
      
      sapply(rv$masses, function(x){

        index <- which(rv$masses == x)
        
        if(input$useRT == TRUE)
          rt.lims <- rtbounds(rv$rts_sec[index], input$rtwind)
        else
          rt.lims <- c(min(rtime(plotting.file[[1]][[1]])),max(rtime(plotting.file[[1]][[1]])))
          
         eic <- chromatogram(plotting.file[[1]][[1]], rt = rt.lims, mz = ppm(rv$masses[index], input$masstol, l = TRUE))
         ms2.data <- findMS2(plotting.file[[1]][[2]], plotting.file[[2]][[2]], rv$masses[index], input$masstol)
         myXIC(eic = eic, ms2.data = ms2.data, exp.rt = rv$rts_sec[index])
      
      }
    )
      dev.off()

      if(i == length(rv$allcalc)){
        output$status <- renderText({paste("Finished building pdfs")})
      }
    }
  })
}

shinyApp(ui = ui, server = server)


