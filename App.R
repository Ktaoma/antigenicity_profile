library(ggplot2)
library(data.table)
library(shiny)
library(plotly)
library(DT)
library(r3dmol)
library(bio3d)
library(Biostrings)
data(BLOSUM100)
source("R/config.R")
source("R/utils.R")
options(shiny.maxRequestSize = 1024 * 1024 * 2048) # 2 GB

# Define UI ----
ui <- fluidPage(
  titlePanel("FluWatch + FluInsight"),
  
  # Main panel content first
  tabsetPanel(type = "tabs",
              tabPanel("Demo",
                       fluidRow(
                         column(
                           width = 7,
                           hr(),
                           div(style = "font-size:17px; text-align:center; border:2px solid #A9A9A9; background-color:#F5F5F5; padding:15px; border-radius:10px; width:100%; margin:0 auto;", textOutput("demo_report")),
                           div(
                             style = "width: 100%; margin: 0 auto;",
                             plotlyOutput("plot_antigencity", height = "500px"),
                           ),
                         ),
                         column(
                           width = 5,
                           #div(style = "font-size: 17px;",style = "text-align: center;",textOutput("num_epitope")),
                           r3dmolOutput("mol", height = "300px"),
                           hr(),
                           plotlyOutput("plot_epitope", height = "300px", width = "100%"),
                         )
                       ),
                       fluidRow(
                         column(
                           width = 12,
                           dataTableOutput("segment_summary")
                         )
                       )),
              tabPanel("User's analysis",
                       br(),
                       fluidRow(
                         column(
                           width = 3,
                           radioButtons( 
                             inputId = "radio", 
                             label = "Analysis Type", 
                             choices = list( 
                               "Genome Assembly with Antigenicity Prediction" = 1, 
                               "Antigenicity Prediction Only" = 2
                             ),
                           )
                         ),
                         column(
                           width = 3,
                           conditionalPanel(
                             condition = "input.radio == 1",
                             fileInput(
                               "fastq_files", 
                               "Upload FASTQ File", 
                               multiple = FALSE,
                               accept = c(".fastq",".fastq.gz", ".fq",".fq.gz","application/gzip",
                                          "application/x-gzip")
                             )
                           ),
                           conditionalPanel(
                             condition = "input.radio == 2",
                             fileInput(
                               "fasta_files", 
                               "Upload FASTA File", 
                               multiple = FALSE,
                               accept = c(".fasta",".fasta.gz", ".fa",".fa.gz")
                             )
                           )
                         ),
                         column(
                           width = 2,
                           dateInput("date1", "Collected Date:", value = "2012-02-29")
                         ),
                         column(
                           width = 2,
                           br(), # adds spacing before button
                           actionButton("run_pipeline", "Start", class = "btn-primary")
                         )
                       ),
                       hr(),
                       fluidRow(dataTableOutput("segment_summary_usr"))
              )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  
  
  etext <- eventReactive(input$run_pipeline,{
    type  <- "both"
    leng <- 300
    system2("bash", args = c("script/01_run.sh",input$fastq_files$datapath, type,leng), stdout = TRUE)
  })
  
  
  
  
  output$segment_summary_usr <- renderDataTable({
    
    
      #sample <- "data/03_tmp/b49.fq.gz"
    
      etext()
    
      mydf <- read.table("result/sample.tsv")
      
      datatable(mydf,
                extensions = 'Buttons',
                options = list(
                  dom = 'Bfrtip',   # show buttons
                  buttons = c('copy', 'csv', 'excel'),
                  pageLength = 4
                ))
      
    })
  
  output$segment_summary <- renderDataTable({
    
    mydf <- data.frame(
      `genomic segment` = paste0("Segment_", 1:8),
      `Consensus Length` = sample(1000:2000, 8, replace = TRUE),
      `Read count` = sample(100:500, 8, replace = TRUE),
      `Averaged genomic depth` = sample(1000:5000, 8, replace = TRUE),
      `N number` = sample(10:50, 8, replace = TRUE)
    )
    
    
    datatable(mydf,
              extensions = 'Buttons',
              options = list(
                dom = 'Bfrtip',   # show buttons
                buttons = c('copy', 'csv', 'excel'),
                pageLength = 4
              ))
    
  })
  
  
  clicked_strain <- reactive({
    
    d <- event_data("plotly_click", source = "A")  # optional: add `source` for isolation
    if (is.null(d)) return("Demo_sample") # default
    strain <- d$customdata
    strain
    
  })
  

  output$num_epitope <- renderText({
    paste0(clicked_epitope()," epitope changes \n in ",clicked_strain())
  })
  
  
  output$demo_report <- renderText({
    
    demo_df <- fread("data/03_tmp/demo_antigenicty.tsv")
    demo_df$group <- ifelse(demo_df$type == "Demo_sample","Demo samples","Circulating strains")
    
    E <- round(unique(0.53-demo_df[demo_df$type == clicked_strain(),2]),3)*100
    
    paste0("Vaccine efficacy of ",clicked_strain(),"\n",
           "compared with A/Victoria/2570/2019 strain is reduced by ",E," percent")
  })
  
  
  clicked_epitope <- reactive({
    d <- event_data("plotly_click", source = "B")  # optional: add `source` for isolation
    if (is.null(d)) return("A") # default
    epitope <- d$customdata[1]
    epitope
    
  })
  
  
  output$plot_antigencity <- renderPlotly({
    
      demo_df <- fread("data/03_tmp/demo_antigenicty.tsv")
      demo_df$group <- ifelse(demo_df$type == "Demo_sample","Demo samples","Circulating strains from GISAID database")
      p <- demo_df |> ggplot(aes(date,E,
                                 color=group,shape=group,
                                 text = paste(
                                   "Date:", date,
                                   "<br>Name:", type,
                                   "<br>Vaccine efficacy:", round(E, 2)
                                 ),customdata = type))+
        geom_point(position = position_dodge(width = .9),alpha=0.7,size=3)+
        geom_hline(yintercept=0.53, linetype="dashed", 
                   color = "red")+
        scale_colour_manual(values = c("grey","red"))+
        scale_x_date(
          date_breaks = "3 month",
          date_labels = "%Y-%m"
        )+
        labs(color="",shape="")+
        ylab("Vaccine efficacy ((u - v)/v)")+
        xlab("Date (Year-month)")+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      
      ggplotly(p, tooltip = "text",source = "A") |>
        layout(
          legend = list(
            orientation = "v",   # vertical legend
            x = 0.05,            # distance from left (0 = far left, 1 = far right)
            y = 0.05,            # distance from bottom (0 = bottom, 1 = top)
            xanchor = "left",
            yanchor = "bottom",
            bgcolor = "rgba(255,255,255,0.6)", # semi-transparent background
            bordercolor = "black",
            borderwidth = 1
          )
        )

  })
  
  
  output$plot_epitope <- renderPlotly({
    
    strain <- clicked_strain()
    epitopes <- c("A", "B", "C", "D", "E")
    mutations <- round(runif(5,min = 0,max = 5))
    
    # Plot
    fig <- plot_ly(
      x = epitopes,
      y = mutations,
      customdata = epitopes,
      type = 'bar',
      marker = list(color = '#4C78A8'),
      source = "B"
    ) %>%
      layout(
        xaxis = list(title = "Epitope"),
        yaxis = list(title = "Number of Mutations"),
        plot_bgcolor = '#ffffff',   # plot area white
        paper_bgcolor = '#ffffff',  # outside plot white
        font = list(size = 12),
        bargap = 0.3
      )
    
    fig
  })
  
  output$mol <- renderR3dmol({
    
    epitope <- clicked_epitope()
  
    # Read your legacy-format PDB file
    pdb <- read.pdb("data/03_tmp/3LZG.pdb")
    pdb_trim <- trim.pdb(pdb, chain = c("A","B","C","D","E","F"))
    pdb_text <- paste(capture.output(write.pdb(pdb_trim)), collapse = "\n")  
  
    pdb_lines <- pdb_trim$atom
    pdb_text <- paste0(
      apply(pdb_lines, 1, function(row) {
        sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %2s",
                as.integer(row["eleno"]), row["elety"], row["resid"],
                row["chain"], as.integer(row["resno"]),
                as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
                row["elety"])
      }),
      collapse = "\n"
    )
    
    
    # Visualize it
    idx <- list("A"=c(132,133,134,135),
                "B"=c(54,155,156,157,160),
                "C"=c(38,40,41,43,44,45),
                "D"=c(89,94,95,96,113,117,163),
                "E"=c(258,259,260,261,263,267))
    
    r3dmol() %>%
      m_add_model(data = pdb_text, format = "pdb") %>% 
      m_set_style(sel = m_sel(chain = c("A","B","C","D","E","F")),
                  style = m_style_cartoon(color = "lightgray")) %>%
      m_zoom_to() %>%
      m_set_style(sel = m_sel(chain="A", resi = idx[[epitope]]),
                  style = m_style_sphere(color = "red",radius = 1)) %>%
      m_add_res_labels(
        m_sel(
          resi = idx[[epitope]],
          chain = "A"
        ),
        style = m_style_label(inFront = T,font = list(size = 6),
                              backgroundOpacity = 0.7,
                              fontColor = "white",showBackground = T,
                              alignment = "bottomRight")
      ) 
      
  }) 
}

# Run the app ----
shinyApp(ui = ui, server = server)

