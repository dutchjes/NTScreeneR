


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
                 actionButton(inputId = "createTIC", label = "Create TIC"),
                 plotOutput(outputId = "TIC"),
                 actionButton(inputId = "create2d", label = "Create 2d plot"),
                 plotOutput(outputId = "plot2d"),
                 selectInput(inputId = "coloring", label = "Color points according to:", choices = c("ionCount", "file", "peaks.count", "charge")),
                 selectInput(inputId = "RAWFileTIC", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
                 ),
        
        tabPanel("EIC", 
                 actionButton(inputId = "createEIC", label = "Create EIC"),
                 selectInput(inputId = "feature", label = "feature to extract", choices = c()), ## filled with masses from P_data
                 textOutput(outputId = "Compound"),
                 plotOutput(outputId = "EIC"),
                 selectInput(inputId = "RAWFileEIC", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
                 ),
        
        tabPanel("MS2 Spectra",
                 selectInput(inputId = "MS2feature", label = "feature to extract", choices = c()), ## filled with masses from P_data
                 plotOutput(outputId = "MS2plot"),
                 selectInput(inputId = "RAWFileMS2", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
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
  
  findMS2 <- function(hd.ms1, hd.ms2, eic, precursor.mz, mztol){
    
    mzlimits <- ppm(precursor.mz, mztol, l = TRUE)
    ms2.scans <- which(hd.ms2$precursorMZ < mzlimits[1] & hd.ms2$precursorMZ > mzlimits[2])
    if(length(ms2.scans)==0){return("NA")}
    ms1.scans <- sapply(ms2.scans, function(x) which.min(abs(hd.ms1$acquisitionNum - hd.ms2$acquisitionNum[x])))
    
    eic.int <- eic@.Data[[1]]@intensity[ms1.scans]
    ms2.rt <- hd.ms2$retentionTime[ms2.scans]
    return(cbind(ms2.rt, eic.int))
    
  }
  
  
  output$userP_data <- renderTable({
    P_data()
  })
  
  rv <- reactiveValues(masses = c(), ## masses to populate mass dropdown
                       rts_min = c(), ## rts to populate rt dropdown
                       rts_sec = c(),
                       names = c(),
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
                       allcalc = c(),
                      # mzmin = c(),
                      # mzmax = c(),
                       msms = list())
  
  
  observeEvent(input$P_data,{
    rv$col <- colnames(P_data())
    
    updateSelectInput(session, "masscol", choices = rv$col)
    updateSelectInput(session, "rtcol", choices = rv$col)
    updateSelectInput(session, "IDcol", choices = rv$col)
    updateSelectInput(session, "ioncol", choices = rv$col)
    
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
    
  #  rv$mzmax <- sapply(rv$masses, function(x) (x + ppm(x, input$masstol, p = TRUE)))
  #  rv$mzmin <- sapply(rv$masses, function(x) (x - ppm(x, input$masstol, p = TRUE)))
  })
  
  # observeEvent(input$ionization,{
  #   
  #   
  # })
  # 
  
  observeEvent(input$filedir,{
    
    rv$filepath <- list.files(input$filedir, pattern = "mzXML", full.names = TRUE, recursive = FALSE)
    rv$filename <- list.files(input$filedir, pattern = "mzXML", full.names = FALSE, recursive = FALSE)
    
    updateSelectInput(session, "RAWFileTIC", choices = rv$filename)
    updateSelectInput(session, "RAWFileEIC", choices = rv$filename)
    updateCheckboxGroupInput(session, "RAWFileOutput", choices = rv$filename)
    
  })
  

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
  
  
  observeEvent(input$createEIC,{

    process.file2 <- which(rv$filename == input$RAWFileEIC)
    process.feature <- which(rv$features == input$feature)
    
    if(length(process.file2) == 0){output$status <- renderText({"Waiting for mzXML files..."})}
    if(length(process.feature) == 0){output$status <- renderText({"Waiting for suspect list..."})}
    
    
  #  rtlim.small <- rv$rts_sec[process.feature] - input$rtwind
  #  rtlim.large <- rv$rts_sec[process.feature] + input$rtwind
  #  ppm <- ppm(rv$masses[process.feature], input$masstol, p = T)
  #  mzwidth <- as.numeric(rv$mzmax[process.feature]) - as.numeric(rv$mzmin[process.feature])

    if(length(process.file2) == 1){
      
      output$status <- renderText({"Ready to plot EIC..."})
    
    rv$feic.ms1 <- readMSData(rv$filepath[process.file2], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
    rv$hdeic.ms1 <- header(rv$feic.ms1)

    rv$feic.ms2 <- readMSData(rv$filepath[process.file2], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
    rv$hdeic.ms2 <- header(rv$feic.ms2)

    output$EIC <- renderPlot({

    if(input$useRT == TRUE){

      
      eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[process.feature], input$masstol) ### problem with plotting ms2 scans comes from shortened eic
      ### with the smaller rtbounds. need to adjust the function
      ### findMS2 to look for the right acquisitionScan rather than
      ### with a which function
      
      eic <- chromatogram(rv$feic.ms1, rt = rtbounds(rv$rts_sec[process.feature], input$rtwind), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      plot(eic)
      abline(v = rv$rts_sec[process.feature], col = "blue")
      #ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[process.feature], input$masstol)
      if(is.na(ms2.data)){legend("topleft", legend = paste("No MS2 scans"))}
      points(ms2.data, col = "red", pch = 16)

    }else{
      
      eic <- chromatogram(rv$feic.ms1, rt = c(min(rtime(rv$feic.ms1)),max(rtime(rv$feic.ms1))), mz = ppm(rv$masses[process.feature], input$masstol, l = TRUE))
      ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[process.feature], input$masstol) ### problem with plotting ms2 scans comes from shortened eic
      
      plot(eic)
      abline(v = rv$rts_sec[process.feature], col = "blue")
      if(is.na(ms2.data)){legend("topleft", legend = paste("No MS2 scans"))}
      points(ms2.data, col = "red", pch = 16)

    }

    })

    # output$EIC <- renderPlot({
    #
    #   xic(rv$feic, mz = rv$masses[process.feature],
    #      # rtlim = c(rtlim.small, rtlim.large),
    #       width = mzWindow(rv$masses[process.feature], input$masstol),
    #      points = TRUE, legend = TRUE, main = rv$features[process.feature])
    #   abline(v = rv$rts_sec[process.feature], col = "blue", lwd = 2)
    #
    # }, height = 1000, width = 1200
    # )

        output$Compound <- renderText({paste(rv$names[process.feature])})
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


  observeEvent(input$process,{

    #output.dir <- input$outputdir

    for(i in 1:length(rv$allcalc)){

      calcfile <- which(rv$filename == rv$allcalc[i])

      output$status <- renderText({paste("Working on file ", rv$allcalc[i], sep = "")})

      pdf(file = paste(input$outputdir, "\\", rv$filename[calcfile], "_PlottedEICs.pdf", sep = ""))


      ms1 <- readMSData(rv$filepath[calcfile], centroided. = TRUE, msLevel. = 1, mode = "onDisk")
      ms1.hd <- header(ms1)
      
      ms2 <- readMSData(rv$filepath[calcfile], centroided. = TRUE, msLevel. = 2, mode = "onDisk")
      ms2.hd <- header(ms2)
      
      sapply(rv$masses, function(x){
        # {
        # index <- which(rv$masses == x);
        # try(xic(ms, x, width = mzWindow(x, input$masstol), points = TRUE, main = rv$features[index]));
        # try(abline(v = as.numeric(as.character(rv$rts_sec[index])), col = "blue", lwd = 2))
        # }
        
        if(input$useRT == TRUE){

          index <- which(rv$masses == x)
          
          eic <- chromatogram(ms1, rt = c(min(rtime(ms1)),max(rtime(ms1))), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          ms2.data <- findMS2(ms1.hd, ms2.hd, eic, rv$masses[index], input$masstol)
          
          
          eic <- chromatogram(ms1, rt = rtbounds(rv$rts_sec[index], input$rtwind), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          plot(eic)
          abline(v = rv$rts_sec[index], col = "blue")
        #  ms2.data <- findMS2(rv$hdeic.ms1, rv$hdeic.ms2, eic, rv$masses[index], input$masstol)
          if(is.na(ms2.data)){legend("topleft", legend = paste("No MS2 scans"))}
          points(ms2.data, col = "red", pch = 16)

        }else{

         index <- which(rv$masses == x)
        
          eic <- chromatogram(ms1, rt = c(min(rtime(ms1)),max(rtime(ms1))), mz = ppm(rv$masses[index], input$masstol, l = TRUE))
          ms2.data <- findMS2(ms1.hd, ms2.hd, eic, rv$masses[index], input$masstol)

          plot(eic)
          abline(v = rv$rts_sec[index], col = "blue")
          if(is.na(ms2.data)){legend("topleft", legend = paste("No MS2 scans"))}
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


