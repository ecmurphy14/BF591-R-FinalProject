## Author: Liz Murphy
## ecmurphy@bu.edu
## BU BF591
## Final Project

library(shiny)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(colourpicker) 
library(RColorBrewer)

# increase max file upload size to 30 MB
options(shiny.maxRequestSize=30*1024^2) 


ui <- fluidPage( theme = shinythemes::shinytheme('flatly'),
    titlePanel(h1('BF591 Final Project', h5(HTML("<p>This tool explores RNA-Seq data of post-mortem human brain tissue from both healthy individuals and those with with Huntington's Disease. <br>
                                                 The normalized counts matrix and DESeq2 results are accessible using the GEO Acession number <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810'>GSE64810</a>. <br>
                                                 See the R script 'processing.R' from this app's <a href='https://github.com/ecmurphy14/BF591-R-FinalProject'>repository</a> to see how the sample metadata was extracted from the GSE SOFT files, as well as how GSEA results were generated.</p>")))),
    tabsetPanel(
        tabPanel('Samples',
                 sidebarLayout(
                     sidebarPanel(fileInput('sample_file', label='Load sample information file', buttonLabel = 'Browse...', accept=c('.csv', '.tsv'), placeholder = 'sample_info.csv')),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Summary",
                                      tableOutput('summary_table')
                                      ),
                             tabPanel("Table",
                                      DT::dataTableOutput('sample_table')
                                      ),
                             tabPanel('Plots',
                                      fluidRow(
                                          splitLayout(cellWidths = c("50%", "50%"), plotOutput("age_death_beeswarm"), plotOutput("age_onset_hist")
                                                      )
                                          )#close fluidRow
                         )#close tabsetPanel
                        ) #close main panel
                     ) #close sidebarLayout
                   )#close samples tabPanel
        ),
        tabPanel("Counts",
                 sidebarLayout(
                     sidebarPanel(p('After changing the filters, please click submit to update the table or plot output.'),
                                  fileInput('counts_file', label='Load normalized counts csv.', buttonLabel = 'Browse...', accept=c('.csv', '.tsv'), placeholder = 'counts.csv'),
                                  sliderInput('slider1', min = 0, max = 100, 'Select variance filter', value = 20, step = 1),
                                  sliderInput('slider2', min = 0, max = 40, 'Select non-zero filter', value = 5, step= 1),
                                  submitButton('Submit', width = '100%')
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Filter Summary",
                                      tableOutput('filter_summary')
                             ),
                             tabPanel("Plots",
                                      fluidRow(
                                          splitLayout(cellWidths = c("50%", "50%"), plotOutput("median_vs_variance"), plotOutput("median_vs_zeros"))
                                      )
                             ),
                             tabPanel('Heatmap',
                                      plotOutput('heatmap')
                             ),
                             tabPanel('PCA',
                                      sidebarLayout(
                                          sidebarPanel(p('After changing the radio buttons, please click submit to update the PCA plot.'),
                                                       radioButtons('xPCA', 'Choose a PC for the x-axis', choices = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'), selected = 'PC1'),
                                                       radioButtons('yPCA', 'Choose a PC for the y-axis', choices =c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'), selected = 'PC2'),
                                                       submitButton('Submit', width = '100%')
                                          ),
                                          mainPanel(plotOutput('pca'))
                                          )
                                      )
                         ) #close main tabset panel
                     ) #close main panel
                 ) #close sidebar layout
        ), #close counts tabpanel 
        tabPanel('DE',
                 sidebarLayout(
                     sidebarPanel(p('After changing the slider value, radio buttons or colors, please click submit to update the table or plot output.'),
                                  fileInput('deseq_file', label='Load differential expression results.', buttonLabel = 'Browse...', accept=c('.csv','.tsv'), placeholder = 'deseq_res.csv'),
                                  sliderInput('slider3', min = -40, max = 0, 'Set a threshold for adjusted p value. This slider determines table filtering and volcano plot coloring.', value = -15, step = 1),
                                  radioButtons('xbutton', 'Choose a column for the x-axis', choices = c('baseMean', 'HD.mean', 'Control.mean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'log2FoldChange'),
                                  radioButtons('ybutton', 'Choose a column for the y-axis', choices =c('baseMean', 'HD.mean', 'Control.mean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'padj'),
                                  colourInput('base', 'Base point color', value = '#22577A'),
                                  colourInput('highlight', 'Highlight point color', value = '#FFCF56'),
                                  submitButton('Submit', icon=icon(name='chart-line', lib='font-awesome'), width = '100%')   
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel('Table', 
                                      DT::dataTableOutput('de_table')),
                             tabPanel('Plots',
                                      plotOutput("volcano_deseq"))
                     )
                     )#close main panel
                 ) #close sidebar layout
                ), #close tabpanel
        tabPanel('GSEA',
                 sidebarLayout(
                     sidebarPanel(p('After changing the adjusted p-value threshold, please click submit to update the table or plot output'),
                         fileInput('fgsea_file', label='Load fgsea results.', buttonLabel = 'Browse...', accept=c('.csv', '.tsv'), placeholder = 'fgsea_res.csv'),
                                  sliderInput('slider4', min = 0.00000001, max = 0.05, 'Select adjusted p-value threshold', value=0.01),
                                  submitButton('Submit', width = '100%')
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel('Top Results',
                                      plotOutput('gsea_bar')),
                             tabPanel('Table',
                                      sidebarLayout(
                                        sidebarPanel(radioButtons('radio_gsea', 'Choose an option', choices = c('All genes', "Negative NES", "Positive NES"), selected='All genes'),
                                                     submitButton('Submit', width = '100%'),
                                                     downloadButton('download_gsea', 'Download Filtered GSEA Results')),
                                        mainPanel(DT::dataTableOutput('gsea'))
                                        )
                             ),
                             tabPanel('Plots',
                                      plotOutput('gsea_scatter'))
                             )#close tabset panel
                         ) #close main panel
                     ) #close sidebar layout
                 ) #Close GSEA tab panel 
        ) #Close outer tabset panel
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
        df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], 
                                               as.factor)
        df$vonsattel_grade <- factor(df$vonsattel_grade)
        summary <- data.frame(`Column Name` = colnames(df)[-1:-2], 
                              Type = sapply(df[-1:-2], class), 
                              Mean_or_Distinct_Value = c(levels(df$type), levels(df$source_name), 
                                                            levels(df$organism), formatC(mean(df$age_at_death), digits=2, format='f'), 
                                                            formatC(mean(df$age_at_onset,na.rm=T), digits = 2, format = 'f'), 
                                                            formatC(mean(df$cag, na.rm=T), digits = 2, format = 'f'), 
                                                            toString(c(levels(df$diagnosis))), formatC(mean(df$duration, na.rm=T), digits = 2, format = 'f'), 
                                                            formatC(mean(df$hv_cortical_score, na.rm=T), digits = 2, format = 'f'), 
                                                            formatC(mean(df$hv_striatal_score, na.rm=T), digits = 2, format = 'f'), 
                                                            formatC(mean(df$mrna_seq_reads), digits=2, format = 'f'), 
                                                            formatC(mean(df$pmi, na.rm=T), digits = 2, format='f'), 
                                                            formatC(mean(df$rin), digits = 2, format = 'f'), levels(df$tissue), 
                                                            toString(c(levels(df$vonsattel_grade)))) )
        summary <- tibble(summary)
        return(summary)
        
    }
    
    
    #' plot_continuous1
    #' @param df Data frame loaded by load_sample_data()
    #'
    #' @return Beeswarm plot of age_at_death variable included in the sample information

    plot_continuous1 <- function(df){
        plot <- ggplot(df, aes(x=diagnosis, y=age_at_death)) +
            geom_beeswarm(cex=3) +
            labs(title='Beeswarm Plot: Age at Death vs. Diagnosis', x='Diagnosis', y='Age at Death')
        return(plot)
    }
    
    #' plot_continuous2
    #' @param df Data frame loaded by load_sample_data()
    #'
    #' @return Histogram plot of age_at_onset variable included in the sample information
        
    plot_continuous2 <- function(df){
        plot <- ggplot(df, aes(x=age_at_onset)) + 
            geom_histogram(binwidth = 4) +
            labs(title="Histogram of Age at Onset of Huntington's Disease Samples", x='Age at Onset', y='Counts')
        return(plot)
    }
    
    #' Outputs for Samples tab
    #' 
    #' 
    output$summary_table <- renderTable({
        req(input$sample_file)
        data <- load_sample_data()
        return(format_summary_table(data))},
        striped = TRUE) 
    
    output$sample_table <- DT::renderDataTable({
        req(input$sample_file)
        data <- load_sample_data()
        return(DT::datatable(data, options = list(pageLength=25)))},
        striped = TRUE) 
    
    output$age_death_beeswarm <- renderPlot({
        req(input$sample_file)
        data <- load_sample_data()
        return(plot_continuous1(data))})
    
    output$age_onset_hist <- renderPlot({
        req(input$sample_file)
        data <- load_sample_data()
        return(plot_continuous2(data))})
    
    
    
    
    #' This following set of functions work within the Counts tab
    #' 
    #' load_counts_data
    #'
    #' @details A reactive function to load the file input for the Counts tab.
    
    load_counts_data <- reactive({
        data <- read_delim(input$counts_file$datapath, delim='\t', col_names = TRUE)
        return(data)
    })
    
    #' filtered_counts_table
    #' @param data Data frame loaded by load_counts_data()
    #' @param slider1 slider input from slider1 input, variance percentile threshold
    #' @param slider2 slider input from slider2 input, number of 0 counts threshold 
    #'
    #' @return A summary table of genes after filtering
    #' @details A function to generate a summary table of genes from the counts matrix after filtering
    
    filtered_counts_table <- function(data, slider1, slider2){
        data_filtered <- dplyr::mutate(data, variance=apply(data[-1], 1, var)) %>%
            dplyr::arrange(desc(variance))%>%
            dplyr::top_frac(n=(1-slider1),wt=variance)
        data_filtered <-dplyr::select(data_filtered, -variance) %>%
            dplyr::filter(rowSums(data_filtered == 0) < slider2)
        tib <- tibble(Samples = ncol(data)-1,
                      Total_Genes = nrow(data),
                      Genes_Pass = nrow(data_filtered),
                      genes_pass_pct = paste(formatC(nrow(data_filtered)/nrow(data)*100, digits=2, format='f'),'%', sep=''),
                      genes_filtered = nrow(data)-nrow(data_filtered),
                      genes_filtered_pct = paste(formatC(((nrow(data)-nrow(data_filtered))/nrow(data))*100, digits=2, format='f'),'%', sep='') 
        )
        transposed <- t(tib) %>%
            as_tibble() %>%
            dplyr::rename(Value = V1) %>%
            dplyr::mutate(Statistic= c("Sample Number", 'Total Genes', 'Passing Genes', '% Passing Genes',"Filtered Genes", "% Filtered Genes" ), 
                          .before = Value)
        return(transposed) 
    }
    #' plot_variance_vs_median
    #' @param data Data frame loaded by load_counts_data()
    #' @param scale_y_axis whether or not to use a log scale on the y axis
    #' @param slider1 slider input from slider1 input, variance percentile threshold
    #' 
    #'
    #' @return A plot of variance vs median counts for each gene, colored by slider input

    
    plot_variance_vs_median <- function(data, scale_y_axis=FALSE, slider1) {
        new_data <- dplyr::mutate(data, median_exp = apply(data[-1],1, median, na.rm = TRUE), variance = apply(data[-1], 1, var))%>%
            dplyr::arrange(median_exp) %>%
            dplyr::mutate(rank_median = rank(median_exp))
        quants <- quantile(new_data$variance, slider1/100)
        new_data$quant <- with(new_data, factor(ifelse(variance < quants[1], 0, 1)))
        if (scale_y_axis==TRUE) {
            plot <- ggplot(new_data)+
                geom_point(aes(x=rank_median, y=variance, color=quant), alpha =0.3)+
                scale_y_log10()+
                scale_color_manual(labels=c(paste('Below ', slider1, '%', ''), paste('Above ', slider1,'%', '')), values = c('gray', 'black'))+
                labs(title='Variance vs. Median Expression', color='Variance', x='Rank of Median Expression', y='Variance')+
                geom_smooth(aes(x=rank_median, y=variance))
        }
        return(plot)
    }
    
    #' plot_zeros_vs_median
    #' @param data Data frame loaded by load_counts_data()
    #' @param scale_y_axis whether or not to use a log scale on the y axis
    #' @param slider2 slider input from slider2, number of zeros threshold
    #' 
    #'
    #' @return A plot of zeros vs median counts for each gene, colored by slider input
    
    plot_zeros_vs_median <- function(data, scale_y_axis=FALSE, slider2){
        new_data <- dplyr::mutate(data, zeros = rowSums(data[-1] == 0),  median_exp = apply(data[-1],1, median, na.rm = TRUE)) %>%
            dplyr::arrange(median_exp) %>%
            dplyr::mutate(rank_median = rank(median_exp))
        if (scale_y_axis==TRUE) {
            plot <- ggplot(new_data)+
                geom_point(aes(x=rank_median, y=zeros, color=zeros<slider2), alpha =0.3)+
                scale_y_log10()+
                scale_color_manual(labels=c(paste('>', slider2, ' '), paste('<', slider2, ' ')), values = c('gray', 'black'))+
                labs(title='Number of Zeros vs. Median Expression', color='Number of Zeros', x='Rank of Median Expression')+
                geom_smooth(aes(x=rank_median, y=zeros))
        }
        return(plot)
    }
    
    #' counts_heatmap
    #' @param data Data frame loaded by load_counts_data()
    #' @param slider1 slider input from slider1, variance percentile threshold
    #' @param slider2 slider input from slider2, number of zeros threshold
    #'
    #' @return A heatmap of filtered counts
    
    counts_heatmap <- function(data, slider1, slider2){
        data_filtered <- dplyr::mutate(data, variance=apply(data[-1], 1, var)) %>%
            dplyr::arrange(desc(variance))%>%
            dplyr::top_frac(n=slider1/100,wt=variance)
        data_filtered <-dplyr::select(data_filtered, -variance) %>%
            dplyr::filter(rowSums(data_filtered[-1] == 0) < slider2)
        data_filtered <- data.matrix(data_filtered[-1])
        pal <- brewer.pal(n = 11, name = 'RdBu')
        heatmap(data_filtered, col = pal, main="Clustered Heatmap of RNA-Seq Counts")
        legend(x = "right", legend = c("low", "medium", "high"),cex = 0.9, fill = pal)
        return()
    }
    
    #' counts_heatmap
    #' @param data tibble loaded by load_counts_data()
    #' @param metadata dataframe loaded by load_sample_data(), to color PCA plot by diagnosis
    #' @param slider1 slider input from slider1, variance percentile threshold
    #' @param slider2 slider input from slider2, number of zeros threshold
    #' @param xPCA radio button input from xPCA, user choice of PC to plot on x axis
    #' @param yPCA radio button input from yPCA, user choice of PC to plot on y axis
    #'
    #' @return A scatter plot of chosen PCs colored by sample diagnosis
    
    counts_pca <- function(data, metadata, slider1, slider2, xPCA, yPCA) {
        data_filtered <- dplyr::mutate(data, variance=apply(data[-1], 1, var)) %>%
            dplyr::arrange(desc(variance))%>%
            dplyr::top_frac(n=slider1/100,wt=variance)
        data_filtered <- dplyr::select(data_filtered, -variance) %>%
            dplyr::filter(rowSums(data_filtered[-1] == 0) < slider2)
        pca_results <- prcomp(t(data_filtered[-1]))
        pca_x <- pca_results$x %>%
            as_tibble(rownames = 'sample_name')%>%
            dplyr::left_join( . , metadata, by = c('sample_name'='sample_name') ) %>%
            dplyr::mutate(diagnosis=as_factor(diagnosis))
        variance <- data.frame(summary(pca_results)$importance)
        variance1 <- toString(variance[2, xPCA]*100)
        variance2 <- toString(variance[2, yPCA]*100)
        plot <- ggplot(pca_x, aes_string(x=xPCA, y=yPCA))+
            geom_point(aes(color=diagnosis))+
            labs(title='PCA on Counts Matrix', x = paste(xPCA, ' (', variance1, '%)'), y =paste(yPCA, ' (', variance2, '%)'), color="Diagnosis") 
        return(plot)
    }
    
    #' Outputs for Counts tab
    #' 
    #'   
    
    output$filter_summary <- renderTable({
        req(input$counts_file)
        data <- load_counts_data()
        colnames(data)[1] <- "gene"
        return(filtered_counts_table(data, input$slider1, input$slider2))}
        , striped = TRUE) 
    
    
    output$median_vs_variance <- renderPlot({
        req(input$counts_file)
        data <- load_counts_data()
        colnames(data)[1] <- "gene"
        return(plot_variance_vs_median(data, scale_y_axis = TRUE, input$slider1))
    })
    
    
    output$median_vs_zeros <- renderPlot({
        req(input$counts_file)
        data <- load_counts_data()
        colnames(data)[1] <- "gene"
        return(plot_zeros_vs_median(data, scale_y_axis = TRUE, input$slider2))
    })
    
    
    output$heatmap <- renderPlot({
        req(input$counts_file)
        data <- load_counts_data()
        colnames(data)[1] <- "gene"
        return(counts_heatmap(data, input$slider1, input$slider2))},
        width = 900,
        height = 700)
    
    
    output$pca <- renderPlot({
        req(input$counts_file)
        data <- load_counts_data()
        metadata <- load_sample_data()
        colnames(data)[1] <- "gene"
        return(counts_pca(data, metadata, input$slider1, input$slider2, input$xPCA, input$yPCA))
    })
    
    #' This following set of functions work within the DE tab
    #' 
    #' load_deseq_data
    #'
    #' @details A reactive function to load the file input for the DE tab.
    
    load_deseq_data <- reactive({
        data <- read.table(input$deseq_file$datapath, sep='\t', header = TRUE)
        return(data)
    })
    
    #' volcano_plot
    #' @param dataf Data frame loaded by load_deseq_data()
    #' @param x_name radio button input to select variable to plot on x axis
    #' @param y_name radio button input to select variable to plot on y axis
    #' @param slider slider input from slider3, adjusted p-value threshold
    #' @param color1 base color input for volcano plot
    #' @param color2 highlight color input for volcano plot
    #'
    #' @return A volcano plot of DESeq2 reults, colored by adjusted p-value
    
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        dataf['neglog10'] <- -log(dataf[y_name], base=10)
        dataf['x_name'] <- dataf[x_name]
        
        plot <-  dataf %>%
            ggplot(aes(x= x_name, y=neglog10)) +
            geom_point(aes(color=padj < 1 * 10 ^ (as.numeric(slider)))) +
            scale_color_manual(name = paste('padj < 1 * 10 ^', slider, sep=' '), values=c(color1, color2)) +
            labs(x= x_name, y= paste('-log10(', y_name, ')', sep = '' )) + 
            theme_bw() +
            theme(legend.position="bottom")
        
        return(plot)
    }
    
    #' draw_table_deseq
    #' @param dataf Data frame loaded by load_deseq_data()
    #' @param slider slider input from slider3, adjusted p-value threshold
    #'
    #' @return A data frame of deseq2 results, filtered by slider input
    
    draw_table_deseq <- function(dataf, slider) {
        dataf <- dataf %>%
            dplyr::filter(padj< (1 * 10^slider)) %>%
            dplyr::mutate(padj = formatC(padj), pvalue = formatC(pvalue))
        return(dataf)
    }
    
    #' Outputs for DE
    #' 
    #' 
    
    output$volcano_deseq <- renderPlot({
        req(input$deseq_file)
        data <- load_deseq_data()
        colnames(data)[1] <- "gene"
        return(volcano_plot(data, input$xbutton, input$ybutton, input$slider3, input$base, input$highlight))
    })
    
    output$de_table <- DT::renderDataTable({
        req(input$deseq_file)
        data <- load_deseq_data()
        colnames(data)[1] <- "gene"
        return(DT::datatable(draw_table_deseq(data, input$slider3), options = list(pageLength=25)))}
        , striped = TRUE) 
    
    #' This last set of functions work within the GSEA tab
    #' 
    #' load_fgsea
    #'
    #' @details A reactive function to load the file input for fgsea_results.tsv
    
    load_fgsea <- reactive({
        data <- read_csv(input$fgsea_file$datapath, col_names = TRUE)
        return(data)
    })
    
    #' plot_gsea_bar
    #' @param fgseaData tibble loaded by load_fgsea()
    #' @param slider4 slider input from slider4, p- adj threshold for pathways to plot
    #' 
    #'
    #' @return A bar chart of NES values from the fgsea reults, filtered by p-adj slider
     
    plot_gsea_bar <- function(fgseaData, slider4){
        fgseaData <- fgseaData %>%
            dplyr::filter(padj<slider4) %>%
            dplyr::mutate(pathway=str_replace_all(pathway, '_', ' ')) %>%
            dplyr::mutate(pathway = str_wrap(pathway, width = 100)) %>%
            dplyr::mutate(pathway=factor(pathway, levels=pathway))
        plot <- ggplot(fgseaData, aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill=NES<0)) +
            coord_flip() +
            labs(x="Pathway", y="Normalized Enrichment Score",
                 title="C2 Canonical Pathways NES from GSEA") + 
            scale_fill_discrete('NES', labels=c('Positive', 'Negative'))
            theme_minimal()
        return(plot)
    }
    
    #' plot_gsea_bar
    #' @param fgseaData tibble loaded by load_fgsea()
    #' @param slider4 slider input from slider4, p- adj threshold for pathways to plot
    #' @param radio_gsea input from the radio button, selecting pathways by NES sign
    #' 
    #'
    #' @return a tibble of fgsea results filtered by padj as well as positive or negative NES
    
    fgsea_res_table <- function(fgseaData, slider4, radio_gsea){
        fgseaData <- dplyr::filter(fgseaData, padj < input$slider4) %>%
            dplyr::arrange(padj) 
        if (radio_gsea == 'Positive NES') {
            fgseaData <- dplyr::filter(fgseaData, NES>0)
        } else if (radio_gsea == 'Negative NES') {
            fgseaData <- dplyr::filter(fgseaData, NES<0)
        } 
        return(fgseaData)
    }
    
    #' plot_gsea_scatter
    #' @param fgseaData tibble loaded by load_fgsea()
    #' @param slider4 slider input from slider4, p- adj threshold for pathways to plot
    #' 
    #'
    #' @return a scatter plt of fgsea results -log10 adjusted p-value vs NES, colored by slider4 padj threshold
    #' 
    
    plot_gsea_scatter <- function(fgseaData, slider4){
        fgseaData <- fgseaData %>%
            dplyr::mutate(neglog10padj = -log(padj, base=10))
        plot <- ggplot(fgseaData, aes(x=NES, y=neglog10padj, color=padj<slider4)) +
            geom_point() +
            labs(title = 'Scatter Plot of GSEA Results', 
                 x = 'Normalized Enrichment Score', 
                 y='Negative log10(padj)') +
            scale_color_manual('Adjusted P-value', 
                               labels = c(paste('>', slider4, sep = ' '), paste('<', slider4, sep = ' ')),
                               values = c('#22577A', '#FFCF56'))
        return(plot)
    }
    
    
    #' Outputs for GSEA tab
    #' 
    #' 
    
    output$gsea_bar <- renderPlot({
        req(input$fgsea_file)
        data <- load_fgsea()
        return(plot_gsea_bar(data, input$slider4))
    }, height = 2500, width = 900)
    
    output$gsea <- DT::renderDataTable({
        req(input$fgsea_file)
        data <- load_fgsea()
        return(DT::datatable(
            fgsea_res_table(data, input$slider4, input$radio_gsea)))},
        striped = TRUE)
    
    output$gsea_scatter <- renderPlot({
        req(input$fgsea_file)
        data <- load_fgsea()
        return(plot_gsea_scatter(data, input$slider4))
    })
    
    output$download_gsea <- downloadHandler(
        filename = 'filtered_gsea.csv',
        content = function(file) {
            data <- 
            write_csv(fgsea_res_table(load_fgsea(), input$slider4, input$radio_gsea), file)
        }
    )
    
}


runApp(shinyApp(ui = ui, server = server), launch.browser = TRUE)