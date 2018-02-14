


#load packages for processing
library(RMassBank)
library(mzR)
library(Hmisc)
library(enviPat)
library(shiny)
library(MSnbase)

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
      
   
 #     checkboxInput(inputId = "MSMS2", label = "Display MS/MS scans?") ## user selected
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
        
        tabPanel("MS2 Output",
                # fileInput(inputId = "MS2topull", label = "Upload file with precursor masses and data files"),
                 selectInput(inputId = "MS2filecol", label = "Column with mzXMLs to extract MS2", choices = c()),
                 selectInput(inputId = "reffilecol", label = "Column with mzXMLs to extract reference MS2", choices = c()),
                 checkboxInput(inputId = "htot", label = "Plot head to tail between File and Reference"),
                 textInput(inputId = "MS2outputdir", label = "Directory to save MS2 pdfs:"),
                 #tableOutput(outputId = "userMS2data"),
                 actionButton(inputId = "processMS2", label = "Save MS2s to pdf")
                 #selectInput(inputId = "MS2feature", label = "feature to extract", choices = c()), ## filled with masses from P_data
                 #plotOutput(outputId = "MS2plot"),
                 #selectInput(inputId = "RAWFileMS2", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
                 ),
        
        tabPanel("Data Output", 
                 checkboxGroupInput(inputId = "RAWFileOutput", label = "Files to process:", choices = c(), inline = TRUE), ##  filled with names in filedir
                 actionButton("selectall", "Select All Files"),
                 #checkboxInput(inputId = "pdf", label = "Save images to pdf?"), 
                 textInput(inputId = "outputdir", label = "Directory to save pdf:"),
                 actionButton(inputId = "process", label = "Save to pdf")
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
  
  findMS2 <- function(hd.ms1, hd.ms2, eic, precursor.mz, mztol){
    
    mzlimits <- ppm(precursor.mz, mztol, l = TRUE)
    ms2.scans <- which(hd.ms2$precursorMZ < mzlimits[1] & hd.ms2$precursorMZ > mzlimits[2])
    if(length(ms2.scans)==0){return(NA)}
    ms1.scans <- sapply(ms2.scans, function(x) which.min(abs(hd.ms1$acquisitionNum - hd.ms2$acquisitionNum[x])))
    
   # eic.int <- as.numeric(as.character(eic@.Data[[1]]@intensity[ms1.scans]))
    eic.int <- as.numeric(as.character(hd.ms2$precursorIntensity[ms2.scans]))
    ms2.rt <- as.numeric(as.character(hd.ms2$retentionTime[ms2.scans]))
    scan.nm <- hd.ms2$acquisitionNum[ms2.scans]
    return(as.data.frame(cbind(scan.nm  = scan.nm, ms2.rt = ms2.rt, eic.int = eic.int)))
    
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
                       ftic = c(),
                       hdtic = c(),
                       feic.ms1 = c(),
                       hdeic.ms1 = c(),
                       feic.ms2 = c(),
                       hdeic.ms2 = c(),
                       rfeic.ms1 = c(),
                       rhdeic.ms1 = c(),
                       rfeic.ms2 = c(),
                       rhdeic.ms2 = c(),
                       allcalc = c(),
                      # mzmin = c(),
                      # mzmax = c(),
                       msms = list(),
                       ms2.data = c(),
                       ms2.files = c(),
                       ref.files = c(),
                       spec.file = c(),
                       spec.ref = c())
  
  
  observeEvent(input$P_data,{
    rv$col <- colnames(P_data())
    
    updateSelectInput(session, "masscol", choices = rv$col)
    updateSelectInput(session, "rtcol", choices = rv$col)
    updateSelectInput(session, "IDcol", choices = rv$col)
    updateSelectInput(session, "ioncol", choices = rv$col)
    updateSelectInput(session, "MS2filecol", choices = rv$col)
    updateSelectInput(session, "reffilecol", choices = rv$col)
    
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
      
      #browser()
      
      if(length(process.file1) == 0){output$status <- renderText({"Waiting for mzXML files..."})}
      
      else{
        
        output$status <- renderText({"Plotting TIC"})
        rv$ftic <- try(readMSData(rv$filepath[process.file1], mode = "onDisk"))
        rv$htic <- tryCatch(header(rv$ftic))
        
        output$status <- renderText({"Read in file"})
        output$TIC <- renderPlot({
          plot(chromatogram(rv$ftic))
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
    process.feature <- which(rv$features == input$feature)
    
    if(length(process.file2) == 0){output$status <- renderText({"Waiting for mzXML files..."})}
    if(length(process.feature) == 0){output$status <- renderText({"Waiting for suspect list..."})}
    

    if(length(process.file2) == 1){
      
      output$status <- renderText({"Ready to plot EIC..."})
    
    rv$feic.ms1 <- readMSData(rv$filepath[process.file2], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
    rv$hdeic.ms1 <- header(rv$feic.ms1)

    rv$feic.ms2 <- readMSData(rv$filepath[process.file2], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
    rv$hdeic.ms2 <- header(rv$feic.ms2)

    output$EIC <- renderPlot({

    if(input$useRT == TRUE){

     ## first the whole eic to get the correct ms2 acquisitionScans 
      eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      rv$ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[process.feature], input$masstol) 

     ## now the shortened eic to plot 
      eic <- chromatogram(rv$feic.ms1, rt = rtbounds(rv$rts_sec[process.feature], input$rtwind), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      plot(eic)
      abline(v = rv$rts_sec[process.feature], col = "blue")
      if(is.na(rv$ms2.data)[1] == TRUE)
        legend("topleft", legend = paste("No MS2 scans"))
      else
        points(rv$ms2.data[,2:3], col = "red", pch = 16)
      
      rtlims <- rtbounds(rv$rts_sec[process.feature], input$rtwind)
      rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rtlims[1] & rv$ms2.data[,2] < rtlims[2])


    }else{
      
      eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      rv$ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[process.feature], input$masstol) ### problem with plotting ms2 scans comes from shortened eic
      
      plot(eic)
      abline(v = rv$rts_sec[process.feature], col = "blue")
      if(is.na(rv$ms2.data)[1] == TRUE)
        legend("topleft", legend = paste("No MS2 scans"))
      else
        points(rv$ms2.data[,2:3], col = "red", pch = 16)

    }

    })
    
    
    updateSelectInput(session, "selectedMSMS", choices = rv$ms2.data[,1])
    output$availMSMSscans <- renderTable({
      rv$ms2.data
    })
    
    output$Compound <- renderText({paste(rv$names[process.feature])})
    output$status <- renderText({"Finished"})
        
    }
    
    observeEvent(input$selectedMSMS, 
                 output$MSMS <- renderPlot({
                   plot(rv$feic.ms2[[which(rv$hdeic.ms2$acquisitionNum == input$selectedMSMS)]], reporters = NULL, full = TRUE) + theme_bw()

    })
    )
    
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
        
        ### sample file
        process.file <- which(grepl(paste(rv$ms2.files[i]), rv$filename) > 0)
        output$status <- renderText({paste("Working with file", rv$ms2.files[i], sep = " ")})
        
        rv$feic.ms1 <- readMSData(rv$filepath[process.file], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
        rv$hdeic.ms1 <- header(rv$feic.ms1)
        
        rv$feic.ms2 <- readMSData(rv$filepath[process.file], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
        rv$hdeic.ms2 <- header(rv$feic.ms2)
        
        ## first the whole eic to get the correct ms2 acquisitionScans
        eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[i], input$masstol, l = TRUE))
        #plot(eic)
        rv$ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[i], input$masstol)
        if(is.na(rv$ms2.data)[1]==TRUE)
          rv$spec.file <- NA
        else{
          
          rtlims <- rtbounds(rv$rts_sec[i], input$rtwind)
          rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rtlims[1] & rv$ms2.data[,2] < rtlims[2])
          get.scan <- which.max(rv$ms2.data[,3])
          if(length(get.scan)<1)
            rv$spec.file <- NA
          else
            rv$spec.file <- rv$feic.ms2[[which(rv$hdeic.ms2$acquisitionNum == rv$ms2.data[,1][get.scan])]]
           ce <- rv$hdeic.ms2$collisionEnergy[which(rv$hdeic.ms2$acquisitionNum == rv$ms2.data[,1][get.scan])]
           id.file <- paste(rv$features[i], ce, rv$ms2.files[i], sep = "_")
          
          
        }
        
        
        ### reference file
        ref.file <- which(grepl(paste(rv$ref.files[i]), rv$filename) > 0)
        output$status <- renderText({paste("Working with file", rv$ref.files[i], sep = " ")})
        
        rv$rfeic.ms1 <- readMSData(rv$filepath[ref.file], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
        rv$rhdeic.ms1 <- header(rv$rfeic.ms1)
        
        rv$rfeic.ms2 <- readMSData(rv$filepath[ref.file], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
        rv$rhdeic.ms2 <- header(rv$rfeic.ms2)
        
        ## first the whole eic to get the correct ms2 acquisitionScans
        eic <- chromatogram(rv$rfeic.ms1, rt = c(min(rtime(rv$rfeic.ms1)),max(rtime(rv$rfeic.ms1))), mz = ppm(rv$masses[i], input$masstol, l = TRUE))
        #plot(eic)
        rv$ms2.data <- findMS2(rv$rhdeic.ms1, rv$rhdeic.ms2, eic, rv$masses[i], input$masstol)
        if(is.na(rv$ms2.data)[1]==TRUE)
          rv$spec.ref <- NA
        else{
          
          rtlims <- rtbounds(rv$rts_sec[i], input$rtwind)
          rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rtlims[1] & rv$ms2.data[,2] < rtlims[2])
          get.scan <- which.max(rv$ms2.data[,3])
          if(length(get.scan)<1)
            rv$spec.ref <- NA
          else
            rv$spec.ref <- rv$rfeic.ms2[[which(rv$rhdeic.ms2$acquisitionNum == rv$ms2.data[,1][get.scan])]]
          ce <- rv$rhdeic.ms2$collisionEnergy[which(rv$rhdeic.ms2$acquisitionNum == rv$ms2.data[,1][get.scan])]
          id.ref <- paste(rv$features[i], ce, rv$ref.files[i], sep = "_")
          
        }
        
        if(!is.na(rv$spec.file) & !is.na(rv$spec.ref)){
          plot(rv$spec.file, rv$spec.ref, main = c(id.file, "/n", id.ref))
          
          browser()
          
          sink(paste(input$MS2outputdir,"\\", id.file, ".txt", sep = ""))
          print(rv$spec.file)
          print(as.data.frame(rv$spec.file))
          print(rv$spec.ref)
          print(as.data.frame(rv$spec.ref))
          sink()
        }
        else
          next
        
      }
      
      output$status <- renderText({"Finished building MS2 pdfs"})
      dev.off()
      
    }else{
      
      pdf(file = paste(input$MS2outputdir, "\\", "ExtractedMS2s.pdf", sep = ""))
      
      for(i in 1:length(rv$ms2.files)){
        
        process.file <- which(grepl(paste(rv$ms2.files[i]), rv$filename) > 0)
        output$status <- renderText({paste("Working with file", rv$ms2.files[i], sep = " ")})
        
        rv$feic.ms1 <- readMSData(rv$filepath[process.file], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
        rv$hdeic.ms1 <- header(rv$feic.ms1)
        
        rv$feic.ms2 <- readMSData(rv$filepath[process.file], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
        rv$hdeic.ms2 <- header(rv$feic.ms2)
        
        ## first the whole eic to get the correct ms2 acquisitionScans
        eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[i], input$masstol, l = TRUE))
        plot(eic)
        rv$ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[i], input$masstol)
        if(is.na(rv$ms2.data)[1]==TRUE)
          next
        else{
          
          rtlims <- rtbounds(rv$rts_sec[i], input$rtwind)
          rv$ms2.data <- subset(rv$ms2.data, rv$ms2.data[,2] > rtlims[1] & rv$ms2.data[,2] < rtlims[2])
          get.scan <- which.max(rv$ms2.data[,3])
          if(length(get.scan)<1)
            next
          else
            plot(rv$feic.ms2[[which(rv$hdeic.ms2$acquisitionNum == rv$ms2.data[,1][get.scan])]], reporters = NULL, full = TRUE) + theme_bw()
          
        }
        
      };
      
      dev.off()
      output$status <- renderText({"Finished building MS2 pdfs"})
      
    }
    
      
  })

  
  
  ##################################
  ### EIC output tab
  ### add compound name to plot
  
  observeEvent(input$process,{


    for(i in 1:length(rv$allcalc)){

      calcfile <- which(rv$filename == rv$allcalc[i])

      output$status <- renderText({paste("Working on file ", rv$allcalc[i], sep = "")})

      pdf(file = paste(input$outputdir, "\\", rv$filename[calcfile], "_PlottedEICs.pdf", sep = ""))

      ms1 <- readMSData(rv$filepath[calcfile], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
      ms1.hd <- header(ms1)
      
      ms2 <- readMSData(rv$filepath[calcfile], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
      ms2.hd <- header(ms2)
      
      sapply(rv$masses, function(x){

        if(input$useRT == TRUE){

          index <- which(rv$masses == x)
          
          ## plot eic over total runtime to get correct ms2 aquisition numbers
          eic <- chromatogram(ms1, rt = c(min(rtime(ms1)),max(rtime(ms1))), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          ms2.data <- findMS2(ms1.hd, ms2.hd, eic, rv$masses[index], input$masstol)
          
          ## get shortened eic for plotting
          eic <- chromatogram(ms1, rt = rtbounds(rv$rts_sec[index], input$rtwind), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          plot(eic)
          abline(v = rv$rts_sec[index], col = "blue")
        #  ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[index], input$masstol)
          if(is.na(ms2.data))
            legend("topleft", legend = paste("No MS2 scans"))
          else
            points(ms2.data, col = "red", pch = 16)

        }else{

         index <- which(rv$masses == x)
        
          eic <- chromatogram(ms1, rt = c(min(rtime(ms1)),max(rtime(ms1))), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          ms2.data <- findMS2(ms1.hd, ms2.hd, eic, rv$masses[index], input$masstol)

          plot(eic)
          abline(v = rv$rts_sec[index], col = "blue")
          if(is.na(ms2.data))
            legend("topleft", legend = paste("No MS2 scans"))
          else
            points(ms2.data, col = "red", pch = 16)
          
        }
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


