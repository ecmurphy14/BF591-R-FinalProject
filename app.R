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
        sidebarPanel(fileInput('file', label='Load sample information csv.', buttonLabel = 'Browse...', accept='.csv', placeholder = 'sample_info.csv'))
            ),
      mainPanel(
        tabPanel("Summary",
              tableOutput('summary_table')
                ),
        tabPanel("Table",
                  tableOutput('sample_table')        
                  ), 
        tabPanel('Plots',
                 plotOutput())                 
                ) #close main panel
      ), #close samples tab panel
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

