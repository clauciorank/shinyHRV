library(shiny)
library(dplyr)
library(pracma)
library(ggplot2)
library(RHRV)
library(tidyverse)
library(gridExtra)
library(shinythemes)
library(DT)



ui <- fluidPage(theme = shinytheme(theme = "cosmo"),
titlePanel("Time domain HRV Analysis"),
  tabsetPanel(
    tabPanel("RR Peak selection",
             sidebarPanel(
               em("Make sure that in your .csv file the time column is in"),strong("seconds"),
                  fileInput("file1", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               uiOutput("colnamestime"),
               uiOutput("colnamesvolt"),
               checkboxInput("invert", "Invert Voltage"),
               checkboxInput("normalize", "Normalize" ),
               actionButton("render", "Render ECG"),
               
                  sliderInput(inputId = "voltageslider", 
                              label= "Voltage Threshold (mV)", 
                              min = 0, 
                              max = 1.5, value = 0.7)
             ),
             mainPanel(
                  tableOutput('head'),
                  plotOutput("EKG"))
             ),
    tabPanel("HRV Analysis",
             sidebarPanel(
               em("The points out of the red lines will be removed"),
             sliderInput("sliderHRV", label = ("Adjust HR range"), min = 0, 
                         max = 250, value = c(60, 120))),
             mainPanel(plotOutput("HRVplot"),
                       tableOutput("HRVtable"))
             ),
    tabPanel("Download File",
             sidebarPanel(
             downloadButton("downloadcsv", "Download .csv"),
             downloadButton("downloadpdf", "Download .pdf")
             ))
  )
)
  


server <- function(input, output, session) {
  
  ############ INPUTS
  
  output$colnamestime <- renderUI({
    req(input$file1)
    colnames <- colnames(read.csv(input$file1$datapath))
    selectInput("colunametime", "Choose time column", colnames)
  })
  
  output$colnamesvolt <- renderUI({
    req(input$file1)
    colnames <- colnames(read.csv(input$file1$datapath))
    selectInput("colunamevoltage", "Choose voltage column", colnames)
  })
  
  ########### REACTIVES
  
  #ECG PROCESSING
  
  ecg <- eventReactive(input$render, {
    req(input$colunamevoltage,input$colunametime)
    elet <- read.csv(input$file1$datapath)
    time.sec <- elet[,input$colunametime]
    voltage.mv <- elet[,input$colunamevoltage]
    if(input$invert) {
      voltage.mv <- elet[,input$colunamevoltage]*-1
    } 
    if(input$normalize) {
      normalize <- function(x) {
        return ((x - min(x)) / (max(x) - min(x)))
      }
      voltage.mv <- normalize(elet[,input$colunamevoltage])
      
      if(input$invert){
        voltage.mv <- normalize(elet[,input$colunamevoltage]*-1)
      }
    }
    eletro <- data.frame(time.sec, voltage.mv)
    eletro$X <- c(1:length(eletro$voltage.mv))
    eletro
  })
  
  # Find Peaks
  peaks <- reactive({
    peaks <- data.frame(findpeaks(ecg()$voltage.mv, 
                                  minpeakheight = input$voltageslider, minpeakdistance = 30, nups = 0))
    sliceposi <- ecg() %>% filter(ecg()$X %in% peaks$X2)
    sliceposi
  })
  
  # HRV
  hrv.data <- reactive({
    rrint<- diff(peaks()$time.sec) 
    
    hrv.data <- NULL
    hrv.data <- CreateHRVData(Verbose = F)
    hrv.data <- LoadBeatVector(hrv.data, peaks()$time.sec)
    hrv.data <- BuildNIHR(hrv.data)
    hrv.data$Beat <- hrv.data$Beat[-1,]
    hrv.data <- CreateTimeAnalysis(hrv.data, size = 40)
    hrv.data
  })
  
  hrv.save <- reactive({
    to_save <- hrv.data()
    to_save$Beat <- to_save$Beat %>% filter(niHR>= min(input$sliderHRV) & niHR <=max(input$sliderHRV))
    to_save <- CreateTimeAnalysis(to_save, size = 40)
    to_save
  })
  
  ########### OUTPUT RENDERS
  
  ## Table EKG
  output$head <- renderTable({
    req(input$file1)
    table <- read.csv(input$file1$datapath)
    head(table,3)
  })
  
  ## EKG PLOT
  
  output$EKG <- renderPlot({
      ggplot(ecg(), aes(x = time.sec, y = voltage.mv))+
      geom_line()+
      geom_hline(yintercept = input$voltageslider)+
      geom_point(data = peaks(), aes(x=time.sec, y =voltage.mv), color = 'red')+
      labs(x= "Time(sec)", y = "Voltage (mV)")+
      theme_light()
  })
  
  
  ## Output HRV PLOT
  output$HRVplot <- renderPlot({
    ggplot(hrv.data()$Beat, aes(x = hrv.data()$Beat$Time, y = hrv.data()$Beat$niHR))+
      geom_line()+
      geom_point()+
      geom_hline(yintercept = min(input$sliderHRV), color = "red")+
      geom_hline(yintercept = max(input$sliderHRV), color = "red")+
      labs(x="Time(sec)", y = "Heart Rate")+
      theme_light()
    
  })
  
  ## Output HRV TABLE
  output$HRVtable <- renderTable({
    TimeAnaysis <- hrv.save()$TimeAnalysis[[1]]
    HeartRate <- 60/(mean(hrv.save()$Beat$RR)/1000)
    data.frame(HeartRate, TimeAnaysis)
  })
  
  
  ## Download CSV
  
  output$downloadcsv <- downloadHandler(
  filename = function(){
    paste("HRV", Sys.Date(),Sys.time(), ".csv", sep = "")
  }, 
  content = function(file) {
    TimeAnaysis <- data.frame(hrv.save()$TimeAnalysis[[1]])
    HeartRate <- 60/(mean(hrv.save()$Beat$RR)/1000)
    save <- data.frame(HeartRate, TimeAnaysis)
    
    write.csv(save, file)
  })
  
  ## Download PDF
  output$downloadpdf <- downloadHandler(filename = function(){
                          paste("HRV",Sys.Date(), Sys.time(),".pdf", sep = "" )},
                          content = function(file){
                            
                            TimeAnaysis <- data.frame(hrv.save()$TimeAnalysis[[1]])
                            HeartRate <- 60/(mean(hrv.save()$Beat$RR)/1000)
                            save <- round(data.frame(HeartRate, TimeAnaysis),2)
                            save1 <-save[,c(1,3,6,7,8,9,10,11,12)]
                  
                            tbl1 <- tableGrob(save1)
                         
                            p2 <- ggplot(hrv.save()$Beat, aes(x = hrv.save()$Beat$Time, y = hrv.save()$Beat$niHR))+
                              geom_line()+
                              geom_hline(yintercept = min(input$sliderHRV), color = "red")+
                              geom_hline(yintercept = max(input$sliderHRV), color = "red")+
                              labs(x="Time(sec)", y = "Heart Rate")+
                              theme_light()
                            p1 <- ggplot(ecg(), aes(x = time.sec, y = voltage.mv))+
                              geom_line()+
                              geom_hline(yintercept = input$voltageslider)+
                              geom_point(data = peaks(), aes(x=time.sec, y =voltage.mv), color = 'red')+
                              labs(x= "Time(sec)", y = "Voltage (mV)")+
                              labs(title = "HRV ANALYSIS")+
                              theme_light()
                            
                        
                            pdf(file, paper = "a4", height = 10)
                            grid.arrange(p1,p2,tbl1,nrow = 3)
                            dev.off()
                          })
}

shinyApp(ui, server)