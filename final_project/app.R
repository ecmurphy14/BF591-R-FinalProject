## Author: Liz Murphy
## ecmurphy@bu.edu
## BU BF591
## Final Project

library(shiny)
library(ggplot2)
library(tidyverse)
library(colourpicker) # you might need to install this.


ui <- fluidPage( 
    titlePanel(h1('BF591 Final Project', h5("This application explores ."))),
    tabsetPanel(
        tabPanel('Samples',
                 sidebarLayout(
                     sidebarPanel(fileInput('sample_file', label='Load sample information csv.', buttonLabel = 'Browse...', accept='.csv', placeholder = 'sample_info.csv')),
                     mainPanel(
                        tabPanel("Summary",
                                tableOutput('summary_table')
                        ),
                        tabPanel("Table",
                                tableOutput('sample_table')        
                        ), 
                        tabPanel('Plots',
                                plotOutput('sample_plot'))
                        ) #close main panel
                     ) #close sidebarLayout
                 ),  #close samples tabPanel
        tabPanel("Counts",
                 sidebarLayout(
                     sidebarPanel(sliderInput('slider1', min = 0, max = 100, 'Select variance filter', value = 20, step = 1),
                                  sliderInput('slider2', min = 0, max = 500, 'Select non-zero filter', value = 50, step= 5)
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Filter Summary",
                                      tableOutput('filter_summary')
                             ),
                             tabPanel("Volcano Plots",
                                      plotOutput('volcano1')
                             ),
                             tabPanel('Heatmap',
                                      plotOutput('heatmap')
                             ),
                             tabPanel('PCA',
                                      sidebarLayout(
                                          sidebarPanel(radioButtons('xbutton', 'Choose a PC for the x-axis', choices = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'), selected = 'PC1'),
                                                       radioButtons('ybutton', 'Choose a PC for the y-axis', choices =c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'), selected = 'PC2')
                                          ),
                                          mainPanel(plotOutput('pca'))
                                      )
                             )
                         ) #close main tabset panel
                     ) #close main panel
                 ) #close sidebar layout
        ), #close counts tabpanel 
        tabPanel('DE',
                 tabsetPanel(
                     tabPanel('Table', 
                              tableOutput('de_table')),
                     tabPanel('Plots',
                              sidebarLayout(
                                  sidebarPanel( fileInput('file', label='Load differential expression results.', buttonLabel = 'Browse...', accept='.csv', placeholder = 'deseq_res.csv'),
                                                p("A volcano plot can be generated with 'log2 fold-change' on the x-axis and 'p-adjusted' on the y-axis."),
                                                radioButtons('xbutton', 'Choose a column for the x-axis', choices = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'log2FoldChange'),
                                                radioButtons('ybutton', 'Choose a column for the y-axis', choices =c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'padj'),
                                                colourInput('base', 'Base point color', value = '#22577A'),
                                                colourInput('highlight', 'Highlight point color', value = '#FFCF56'),
                                                sliderInput('slider', min = -300, max = 0, 'Select the magnitude of the p adjusted coloring:', value = -150, step = 1),
                                                submitButton('Plot', icon=icon(name='chart-line', lib='font-awesome'), width = '100%')
                                  ),
                                  mainPanel(
                                      tabsetPanel(
                                          tabPanel("Plot", plotOutput("volcano2")),
                                          tabPanel("Table", tableOutput("table"))
                                      )
                                  )#close main panel
                              ) #close sidebar layout
                     ) #close plots tabpanel
                 ) #close tabset panel
        ), #close DE tab panel
        tabPanel('GSEA',
                 sidebarLayout(
                     sidebarPanel(fileInput('file', label='Load fgsea results.', buttonLabel = 'Browse...', accept='.csv', placeholder = 'fgsea_res.csv')
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel('Top Results',
                                      sidebarPanel(sliderInput('slider3', min = 0.00001, max = 0.1, 'Select adjusted p-value threshold', value=0.05, step = 0.00001)),
                                      mainPanel(plotOutput('gsea_bar'))),
                             tabPanel('Table',
                                      sidebarPanel(sliderInput('slider4', min = 0.00001, max = 0.1, 'Select adjusted p-value threshold', value=0.05, step = 0.00001),
                                                   radioButtons('radio_gsea', 'Choose an option', choices = c('All genes', "Negative NES", "Positive NES", selected='All genes'))),
                                      mainPanel(tableOutput('gsea'))
                             ),
                             tabPanel('Plots',
                                      sidebarPanel(sliderInput('slider5', min = 0.00001, max = 0.1, 'Select adjusted p-value threshold', value=0.05, step = 0.00001)),
                                      mainPanel(plotOutput('gsea_scatter'))
                             )
                         )#close tabset panel
                     )#close main panel
                 ) #close sidebar layout
        )#CLose GSEA tab panel
    )#Close outer tabset panel
) #Close fluidPage 

# Define server logic required to load different files, and return the various tables and plots
server <- function(input, output, session) {
    
    #' This first set of functions work within the Samples tab
    #' 
    #' load_sample_data
    #'
    #' @details A reactive function to load the first file input, the sample_info.csv
    
    load_sample_data <- reactive({
        data <- read.csv(input$sample_file$datapath, header = TRUE)
        return(data)
    })
    
    #' format_summary_table
    #' @param df Data frame loaded by load_sample_data()
    #'
    #' @return Data frame summarizing sample info for Samples subtab 'Summary'
    #' @details A function to generate a summary table of sample information from sample_info.csv
    
    format_summary_table <- function(df){
        
    }
    
    #' format_sample_table
    #' @param df Data frame loaded by load_sample_data()
    #'
    #' @return Data frame of sample info for Samples subtab 'Table'
    #' @details A function to generate a table of sample information from sample_info.csv
    
    format_sample_table <- function(df){
        
    }
    
    #' plot_continuous
    #' @param df Data frame loaded by load_sample_data()
    #'
    #' @return Histograms of continuous variables included in the sample information
    #' @details A function to generate plots of continuous variables from sample_info.csv
    
    plot_continuous <- function(df){
        
    }
    
    
    #' Outputs for Samples tab
    #' 
    #' 
    output$summary_table <- renderTable(format_summary_table(load_sample_data()),
                                        striped = TRUE) 
    
    output$sample_table <- renderTable(format_sample_table(load_sample_data()),
                                        striped = TRUE) 
    
    output$sample_plot <- renderPlot(plot_continuous(load_sample_data()),
                                     width = 900,
                                     height = 700)
    
}

shinyApp(ui = ui, server = server)