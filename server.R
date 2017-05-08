
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(RColorBrewer)
library(tidyverse)
library(rbokeh)

#Please upload your Rdata file here
load('Your file.Rdata')

#Paste the path to your RaceID class_ext file here
source('RaceID2_class_ext.R')

#extract cluster genes
clust <- clustdiffgenes(sc, pvalue=0.05)

##Extract tsne variables
cells <- as.data.frame(sc@tsne)
cell_ids <- cbind(as.data.frame(names(sc@cpart)), as.data.frame(sc@cpart))
head(cells)
cell_df <- cbind(cells, cell_ids)
head(cell_df)
colnames(cell_df) <- c('x', 'y', 'Cell_ID', 'Cluster')

#Please update your Cell ID pattern and distinguish your groups based on the underlying condition - here 'Sick' or 'Healthy'
cell_df$Condition <- ifelse(grepl("Your_Cell_ID_Pattern", cell_df$Cell_ID), 'Sick', 'Healthy')

#get mean coordinates for clusters
cell_df_mean <- cell_df %>%
        group_by(Cluster) %>%
        summarise(mean_x = median(x),
                  mean_y = median(y))



shinyServer(function(input, output) {
        
        #RaceID class file to source:
    
        #output$contents <- renderTable({
        #        inFile <- input$file1
        #        if (is.null(inFile))
        #                return(NULL)
        #        else{
        #                tdata <- as.matrix(   read.table(inFile$datapath))
        #                head(tdata, n = 2)
         #               source("./CODE/run_myApp.r", local = TRUE)
        #        }
       # })
        
        #output$contents <- renderTable({
         #       inFile <- input$dataFile
         #       if (is.null(inFile))
        #                return(NULL)
        #        else{
        #                tdata <- as.matrix(   read.table(inFile$datapath))
        #                head(tdata, n = 2)
         #               load("./CODE/run_myApp.r", local = TRUE)
         #       }
        #})
        
        output$bokehPlot <- renderRbokeh({
                
                figure(width = 1050, height = 600) %>%
                        ly_points(x, y, data = cell_df,
                                  color = factor(Cluster), glyph = Condition,
                                  hover = c(Cluster, Condition, Cell_ID)) %>% 
                        ly_text(mean_x, mean_y, data = cell_df_mean, text = Cluster)
        })
        #plot Tsne with cluster number
        output$tsnePlot <- renderPlot({
        
                plottsne(sc)
        })
        
        #Plot Tsne with group labels
        output$tsnelabelPlot <- renderPlot({
                
                # generate bins based on input$bins from ui.R
                #x    <- faithful[, 2]
                #bins <- seq(min(x), max(x), length.out = input$bins + 1)
                
                # draw the histogram with the specified number of bins
                #hist(x, breaks = bins, col = 'darkgray', border = 'white')
                plotsymbolstsne(sc,types=sub("(\\_|\\.).+","", names(sc@ndata)))
        })
       
        #Plot tsne plot Gene expression
        output$tsneExpPlot <- renderPlot({
                
                plotexptsne(sc, name2id (input$gene, rownames(sc@ndata)), logsc=input$Log) # to plot expression
        })
        
        #MA plot 
        output$maPlot <- renderPlot({
                cl1 <- names(sc@cpart[sc@cpart %in% input$cluster1])
                cl2 <- names(sc@cpart[sc@cpart %in% input$cluster2])
                diffexp <- diffexpnb(sc@ndata,cl1,cl2,norm=F,vfit=sc@background$vfit) # 
                plotdiffgenesnb(diffexp,show_names = T)
                
        })
        
        #plot upregulated genes
        output$clusterTableUp <- renderTable({
                
                        a <- cbind(as.data.frame(rownames(as.data.frame(clust[input$cluster]))), as.data.frame(clust[input$cluster]))
                        colnames(a) <- c('Gene', 'Mean overall', 'Mean in Clust', 'fc', 'p-Value')
                        a[a$fc>1,]
                        

        }, caption=paste("Upregulated Genes"),
        caption.placement = getOption("xtable.caption.placement", "top"),
        caption.width = getOption("xtable.caption.width")
        )
        
        #Plot downregulated genes
        output$clusterTableDown <- renderTable({
                
                        a <- cbind(as.data.frame(rownames(as.data.frame(clust[input$cluster]))), as.data.frame(clust[input$cluster]))
                        colnames(a) <- c('Gene', 'Mean overall', 'Mean in Clust', 'fc', 'p-Value')
                        a[a$fc<1,]
             
                
        }, caption=paste("Downregulated Genes"),
        caption.placement = getOption("xtable.caption.placement", "top"),
        caption.width = getOption("xtable.caption.width", NULL)
        )
        
})
