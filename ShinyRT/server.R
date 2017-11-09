#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)

source('load.data.R')

# cap uploads at 1 GB
options(shiny.maxRequestSize=1000*1024^2) 

# load STAN function
fit_RT <- readRDS('./../fit_RT.rds')

shinyServer(function(input, output, session) {
  
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  load.data(inFile$datapath)
  
  #output$ev <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    #inFile <- input$file1
    
    # catch null file
    #if (is.null(inFile))
    #  return(NULL)
    
    #read_tsv(inFile$datapath)
    
  #})
})
