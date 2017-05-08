
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)
library(rbokeh)
#library(parallel)

# Calculate the number of cores
#no_cores <- detectCores() - 1

# Initiate cluster
#cl <- makeCluster(no_cores)

shinyUI(fluidPage(

  # Application title
  titlePanel("Single cell data viewer"),

  # Sidebar with input files and data metrics
  sidebarLayout(
    sidebarPanel(
           
        #fileInput('file1', 'RaceID2_class_ext.R to source', 
             #         accept = '.R'),
            #fileInput('dataFile', 'Select .Rdata File', 
           #           accept = '.Rdata'),
       
      textInput("gene",
                 "Gene Expression Tsne - Enter Gene:", value = 'Cx3cr1'),
      checkboxInput("Log",
                "Expression on Logarithmic Scale", value = FALSE, width = NULL),
      numericInput("cluster",
                "Differential Gene Expression - Enter Cluster number:", value = 1, min = 1, max = 30, step=1),
      tableOutput('clusterTableUp'),
      tableOutput('clusterTableDown')
      ),

    # Show plots from the sc data
    mainPanel(
            plotOutput("tsneExpPlot"),
            rbokehOutput("bokehPlot"),
      plotOutput("tsnelabelPlot", width = "84%"),
      plotOutput("tsnePlot", width = "84%"),
      
      numericInput("cluster1",
                   "DIFFERENTIAL GENE EXPRESSION ANALYSIS First Cluster:", value = 1, min = 1, max = 30, step=1),
      numericInput("cluster2",
                   "Second Cluster:", value = 2, min = 1, max = 30, step=1),
      plotOutput('maPlot', width = "84%")
      
      
    )
  )
))

