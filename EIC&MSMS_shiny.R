


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
      selectInput(inputId = "IDcol", label = "Column with ID", choices = c()),
      selectInput(inputId = "ioncol", label = "Column with ionization mode", choices = c()),
      selectInput(inputId = "ionization", label = "Which Ionization mode?", choices = c()),
      
      selectInput(inputId = "feature", label = "feature to extract", choices = c()), ## filled with masses from P_data
      numericInput(inputId = "masstol", label = "ppm tolerence", value = c(5)),
      numericInput(inputId = "rtwind", label = "rt window (sec)", value = c(60)), ## user selected
      
      textInput(inputId = "filedir", label = "Directory with mzXML files:"), ## path to RAW files
      
      textOutput(outputId = "status")
      
   
 #     checkboxInput(inputId = "MSMS2", label = "Display MS/MS scans?") ## user selected
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", tableOutput(outputId = "userP_data")),
        
        tabPanel("TIC", 
          plotOutput(outputId = "TIC"),
          actionButton(inputId = "create2d", label = "Create 2d plot"),
          plotOutput(outputId = "plot2d"),
          selectInput(inputId = "coloring", label = "Color points according to:", choices = c("ionCount", "file", "peaks.count", "charge")),
          selectInput(inputId = "RAWFileTIC", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
          ),
        
        tabPanel("EIC", 
                 textOutput(outputId = "Compound"),
                 plotOutput(outputId = "EIC", inline = TRUE),
                 selectInput(inputId = "RAWFileEIC", label = "Current mzXML file:", choices = c(), multiple = FALSE) ##  filled with names in filedir
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
  
  output$userP_data <- renderTable({
    P_data()
  })
  
  rv <- reactiveValues(masses = c(), ## masses to populate mass dropdown
                       rts = c(), ## rts to populate rt dropdown
                       names = c(),
                       ESImodes = c(),
                       col = c(),
                       features = c(),
                       filepath = c(),
                       filename = c(),
                       ftic = c(),
                       hdtic = c(),
                       fteic = c(),
                       hdeic = c(),
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
    rv$rts <- as.numeric(as.character(P_data()[,tike]))
    if(input$rtunit == "seconds"){
      rv$rts <- as.numeric(as.character(P_data()[,tike])) / 60
    }

    rv$features <- paste(rv$masses, rv$rts, sep = "_")
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
  

  observeEvent(input$RAWFileTIC,{

      process.file1 <- which(rv$filename == input$RAWFileTIC)

      rv$ftic <- try(openMSfile(rv$filepath[process.file1]))
      rv$htic <- try(header(rv$ftic))

      output$TIC <- renderPlot({
        chromatogram(rv$ftic)
      })
      
      observeEvent(input$create2d,{
        spec <- try(readMSData(rv$filepath[process.file1]))
        output$plot2d <- renderPlot({
          plot2d(spec, z  = input$coloring)
      })

      #  plot(MSmap(rv$ftic, lowMz = min(rv$masses), highMz = max(rv$masses), resMz = ppm(max(rv$masses), input$masstol, p = TRUE)), aspect = 1, allTicks = FALSE)
      })
  })
  
  
  observeEvent(c(input$RAWFileEIC, input$feature),{
    
    process.file2 <- which(rv$filename == input$RAWFileEIC)
    process.feature <- which(rv$features == input$feature)
    rtlim.small <- rv$rts[process.feature] - input$rtwind
    rtlim.large <- rv$rts[process.feature] + input$rtwind
  #  ppm <- ppm(rv$masses[process.feature], input$masstol, p = T)
  #  mzwidth <- as.numeric(rv$mzmax[process.feature]) - as.numeric(rv$mzmin[process.feature])
    
    
    rv$feic <- try(openMSfile(rv$filepath[process.file2]))
    rv$heic <- try(header(rv$feic))

        output$EIC <- renderPlot({
      
      xic(rv$feic, mz = rv$masses[process.feature], 
         # rtlim = c(rtlim.small, rtlim.large), 
          width = mzWindow(rv$masses[process.feature], input$masstol),
         points = TRUE, legend = TRUE, main = rv$features[process.feature])
      abline(v = rv$rts[process.feature]*60, col = "blue", lwd = 2)
          
    }, height = 1000, width = 1200)
        
        output$Compound <- renderText({paste(rv$names[process.feature])})
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
      
      calcfile <- which(rv$filename == input$RAWFileOutput[i])
      
      output$status <- renderText({paste("Working on file ", input$RAWFileOutput[i], sep = "")})
      
      pdf(file = paste(input$outputdir, "\\", rv$filename[calcfile], "_PlottedEICs.pdf", sep = ""))
      
      
      ms <- try(openMSfile(rv$filepath[calcfile]))
      hd <- try(header(ms))

      sapply(rv$masses, try(function(x) 
        {
        index <- which(rv$masses == x);
        try(xic(ms, x, width = mzWindow(x, input$masstol), points = TRUE, main = rv$features[index]));
        try(abline(v = as.numeric(as.character(rv$rts[index])) * 60, col = "blue", lwd = 2))
        })
        )
      
      dev.off()
      
      if(i == length(rv$allcalc)){
        
        output$status <- renderText({paste("Finished building pdfs")})
      }
    }
      
    })

}

shinyApp(ui = ui, server = server)


