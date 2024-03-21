### Global Options ####
options(shiny.maxRequestSize = 75*1024^2, # this sets the max file size for uploads to 75MB
        spinner.color = "#8e44ad", # this sets the options for all of the loading spinners
        spinner.type = 4,
        spinner.color.background = "#fdfefe")

### Libraries ####
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(forcats)
library(scales)
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinythemes)
library(shinyWidgets)
library(data.table)
library(purrr)
library(stringr)
library(zip)
library(Rspc)
library(naniar)
library(scales)
library(shinyalert)
library(DT)
library(parallel)
library(future)
library(parallelly)
library(furrr)
library(microbenchmark)
library(lubridate)
library(hablar)

### UI Start ####
ui <- fluidPage(

### STYLING ####
  theme = shinytheme("cerulean"), 
  
  useShinyjs(), #shinyjs is needed to use hide/show on the download output widget
  useShinyalert(), #shinyalert needed to generate alert windows

  titlePanel(tags$b("shinyMBA"), # CSS tag that bolds the title text
             windowTitle = "shinyMBA"), # set the text that appears in browser windows/tabs

  
  tags$style(type = "text/css", ".navbar{ font-size: 17px; }"), # CSS tag used to change the text size in the navigation bar

  tags$head(
    tags$style(HTML(" .shiny-output-error-merge{ color: crimson; font-weight: bold; }", # customizing font for file merging error messages
                    " .shiny-output-error-ct{font-size: large; font-weight: bold; }" ))), # customizing font for control tracking select antigen/controls warning

  tagList(
   tags$head(tags$script(type="text/javascript", src = "cdcalign.js")), # javascript code that adds in CDC's logo with the correct alignment
   navbarPage(title = "",
              inverse = TRUE, # inverse color scheme

### MODULE: Upload ####             
             tabPanel("Upload",
                      sidebarLayout( # this dictates what will go in the side panel of the app. I'm choosing to put most of the widgets here.
                        sidebarPanel(awesomeRadio(inputId = "master_source",
                                                  label = "Data source",
                                                  choices = c("Only xPONENT/BPM", 
                                                              "Only shinyMBA datasets", 
                                                              "Append xPONENT/BPM to shinyMBA datasets"),
                                                  checkbox = TRUE),
                                     conditionalPanel(condition = "input.master_source == 'Only xPONENT/BPM' || input.master_source == 'Append xPONENT/BPM to shinyMBA datasets'",
                                       fileInput(inputId = "file", # This widget allows users to upload their own xPONENT output .csv files into the app
                                                 label = "Select xPONENT (.csv) and/or Bio-Plex Manager (.xlsx) files",
                                                 accept = c(".csv", ".xlsx"), #the file must be an xPONENT .csv or multi-tab BioPlex Manager .xlsx
                                                 multiple = TRUE)), #allows for uploading multiple files
                                     conditionalPanel(condition = "input.master_source == 'Only shinyMBA datasets' || input.master_source == 'Append xPONENT/BPM to shinyMBA datasets'",
                                       fileInput(inputId = "prev_datasets",
                                                 label = "Select previous shinyMBA datasets (.csv)",
                                                 accept = ".csv",
                                                 multiple = TRUE))),
                        mainPanel())),

### MODULE: Bead Count / MFI ####             
             tabPanel("Bead Count / MFI",
                      sidebarLayout(
                        sidebarPanel(radioButtons(inputId = "qc_choice", # widget that lets users switch between bead count and MFI QC 
                                                  label = "Choose QC metric",
                                                  choices = c("Bead Count", "MFI"),
                                                  selected = "Bead Count"), # bead count selected by default
                                     
                                     uiOutput("batches"), #reactive widget that populates as a list of plate ID's for the user to choose from once the xPONENT files are loaded
                                     uiOutput("ab_var"), # reactive widget that populates as a list of antigens/antibodies for the user to choose from once they upload their xPONENT files
                                     conditionalPanel(condition = "input.qc_choice == 'Bead Count'", # section displayed only when bead count QC is selected
                                                      numericInputIcon(inputId = "uthresh", # This allows users to manipulate the upper threshold value
                                                                   label = "Upper threshold",
                                                                   value = 35, # default value is 35 beads/well
                                                                   min = 1, # minimum value is 1 bead/well
                                                                   max = 1000, # maximum value is 1,000 beads/well,
                                                                   icon = list(NULL, "beads/well"),
                                                                   width = '200px'),
                                                      numericInputIcon(inputId = "lthresh", # This allows users to manipulate the lower threshold value
                                                                   label = "Lower threshold",
                                                                   value = 20, # default value is 20 beads/well
                                                                   min = 1, # minimum value is 1 bead/well
                                                                   max = 999, # maximum value is 999 beads/well
                                                                   icon = list(NULL, "beads/well"),
                                                                   width = '200px'),
                                                      radioButtons(inputId = "utlt_platefail", # widget that lets users define well failure as bead count either < upper threshold or < lower threshold
                                                                   label = "Set failed well criteria",
                                                                   choices = c("< Upper threshold", "< Lower threshold"), 
                                                                   selected = "< Lower threshold"), # < upper threshold selected by default
                                                      numericInputIcon(inputId = "bc_platefail", # widget that lets users define antigen-plate failure based on a % of failed wells
                                                                       label = "Set failed plate criteria",
                                                                       value = 30, # 30% set as the default
                                                                       min = 0,
                                                                       max = 100,
                                                                       icon = list(NULL, "% failed wells"),
                                                                       width = '200px')),
                                     conditionalPanel(condition = "input.qc_choice == 'MFI'", # section displayed only when MFI QC is selected
                                                      radioGroupButtons(inputId = "mfi_plot_switch", # button widget that switches between MFI plate plot and MFI heat map
                                                                        label   = "Plot type",
                                                                        choices = c("Background/Warning", "Heat Map"),
                                                                        individual = TRUE,
                                                                        checkIcon = list(
                                                                            yes = tags$i(class = "fa fa-circle", 
                                                                                         style = "color: steelblue"),
                                                                            no = tags$i(class = "fa fa-circle-o", 
                                                                                        style = "color: steelblue"))),
                                                      numericInputIcon(inputId = "mfi_warning", # widget that lets users flag wells based on an arbitrary MFI value
                                                                        label = "Set well warning threshold",
                                                                        value = 100, # < 100 MFI selected as default
                                                                        min = 1,
                                                                        max = 1000000,
                                                                        icon = list("<", "MFI"),
                                                                        width = '225px')),
                                     width = 3),
                                                       

                        mainPanel(
                          tabsetPanel(
                            tabPanel("Plots", # QC plots displayed here
                                     conditionalPanel(condition = "input.qc_choice == 'Bead Count'", # plots only displayed when bead count QC is selected
                                                      
                                                      tabPanel("Bead Count Plate Plot", plotOutput("plateplot") %>% withSpinner()), # 96-well bead count plot displayed here
                                                      tabPanel("QC Plot", plotlyOutput("qcplot") %>% withSpinner())), # bead count fluctuation plot displayed here
                                     conditionalPanel(condition = "input.qc_choice == 'MFI'", # plots only displayed when MFI QC is selected
                                                      tabPanel("MFI Plate Plot", plotOutput("plateplot_mfi") %>% withSpinner()))), # 96-well MFI plot displayed here
                            
                            tabPanel("Plate-Antigen Repeats Summary ",
                                     conditionalPanel(condition = "input.qc_choice == 'Bead Count'", # summary table only displayed when bead count QC is selected
                                                      tabPanel("Bead Count Repeats by Plate", dataTableOutput("bc_repeats") %>% withSpinner()))), # plate-antigen repeats table displayed here
                            
                            tabPanel("All Plates Repeats Summary ",
                                     conditionalPanel(condition = "input.qc_choice == 'Bead Count'", # summary table only displayed when bead count QC is selected
                                                      tabPanel("Bead Count Repeats (Whole Plates)", dataTableOutput("bc_repeats2") %>% withSpinner())))), # entire plate antigen failture table displayed here
                          width = 9)
                          )), 

### MODULE: Control Tracking ####             
             tabPanel("Control Tracking",
                      sidebarLayout(
                        sidebarPanel(radioButtons(inputId = "lj_choice", # widget that lets users choose their data source for generating LJ plots
                                                  label = "Data Source",
                                                  choices = c("xPONENT/Bio-Plex Manager", "Formatted Custom Data"),
                                                  selected = "xPONENT/Bio-Plex Manager" ), # xPONENT/Bio-Plex Manager selected by default
                                     conditionalPanel(condition = "input.lj_choice == 'Formatted Custom Data'", # widget section only appears when formatted custom data is selected
                                                      fileInput(inputId = "ljfile", # widget that allows users to upload their custom formatted .csv files
                                                                label = "Upload .csv formatted data",
                                                                accept = ".csv")),
                                    
                                     pickerInput(inputId = "chart_type", # drop down list widget that switches between control chart types
                                                 label = "Select chart type",
                                                 choices = c("Levey-Jennings", "Shewhart", "Moving Range"),
                                                 options = list(
                                                   style = "btn-primary")),
                                    
                                     conditionalPanel(condition = "input.chart_type == 'Shewhart'", # widget section that appears when the Shewhart chart option is selected 
                                                      numericInputIcon(inputId = "shew_ref", # widget that lets users set how many of the initial observations are used to determine confidence limits
                                                                   label = "Number of reference points",
                                                                   value = 2,
                                                                   min = 2, 
                                                                   max = NULL,
                                                                   icon = list(NULL, "data points"),
                                                                   width = '200px')),
                                     conditionalPanel(condition = "input.lj_choice == 'xPONENT/Bio-Plex Manager'", # widget section only appears when xPONENT/BPM is selected
                                                      pickerInput(inputId = "samples2", label = "Specify controls", choices = NULL, options = list(`actions-box` = TRUE), multiple = TRUE), # widget that allows users to tell the app what their controls are
                                                      pickerInput(inputId = "ag_spec", label = "Specify antigens/antibodies", choices = NULL, options = list(`actions-box` = TRUE), multiple = TRUE), # widget that allows users to tell the app which antigens should be paired with the controls
                                                          prettyCheckbox(inputId = "lj_avg", # checkbox widget that tells the app whether to average plate replicates before plotting
                                                                         label = "Average replicates?",
                                                                         value = FALSE,
                                                                         icon = icon("check"),
                                                                         status = "success",
                                                                         animation = "smooth")),
                                     
                                     radioButtons(inputId = "lj_type", # widget that lets users specify one of two LJ violation types 
                                                  label = "Violation Type",
                                                  choices = c("Nelson Rules", "2CL and 3CL"), # violation types
                                                  selected = "2CL and 3CL"), # 2CL and 3CL rules selected by default
                                     width = 3),
                                     
                        mainPanel(
                          tabsetPanel(id = "ljpanels",
                                      tabPanel("Plots", # interactive LJ plot displayed here
                                               conditionalPanel(condition = "input.lj_choice == 'xPONENT/Bio-Plex Manager'", # widget section only appears when xPONENT/BPM is selected
                                                                dropdownButton(inputId = "lj_xp_dd", # drop down list widget that lets users change which control + antigen combination is displayed in the plot area
                                                                               selectInput(inputId = "ab2", label = "Choose an antigen/antibody", choices = NULL), # widget that populates with a list of the antigens/antibodies from the user's xPONENT/BPM files.
                                                                               selectInput(inputId = "samples3", label = "Choose a control", choices = NULL), # widget that populates with all the controls/samples specified by the user
                                                                               status = "info",
                                                                               circle = TRUE,
                                                                               icon = icon("cog"), width = "300px",
                                                                               tooltip = tooltipOptions(title = "Click to change inputs!"))),
                                               conditionalPanel(condition = "input.lj_choice == 'Formatted Custom Data'", # widget section only appears when formatted custom data is selected
                                                                dropdownButton(inputId = "lj_custom_dd",  # drop down list widget that lets users change which control + antigen combination is displayed in the plot area
                                                                               selectInput(inputId = "ab2_custom", label = "Choose an antigen/antibody", choices = NULL), # widget that populates with a list of the antigens/antibodies from the custom dataset
                                                                               selectInput(inputId = "samples2_custom", label = "Choose a control", choices = NULL), # widget that populates with a list of the controls from the custom dataset
                                                                               status = "info",
                                                                               circle = TRUE,
                                                                               icon = icon("cog"), width = "300px",
                                                                               tooltip = tooltipOptions(title = "Click to change inputs!"))),
                                               plotlyOutput("lj") %>% withSpinner(), # displayed plots must be plotly objects
                                               conditionalPanel(condition = "input.chart_type == 'Moving Range'",
                                                                plotlyOutput("mr") %>% withSpinner()), # moving range chart only displayed if moving range chart type is selected
                                               conditionalPanel(condition = "input.lj_type == 'Nelson Rules'", # Nelson rules table only displayed when Nelson Rules is selected
                                                                tableOutput("nelson"))), # Nelson rules table explaining each of the violation types
                                      tabPanel("Flagged plates", # summary table that displays the control tracking analysis results
                                               dataTableOutput("lj_summary") %>% withSpinner()),
                                      tabPanel("Summary", # summary table that displays summary statistics for each antigen/antibody-control combination
                                               dataTableOutput("lj_big_summary") %>% withSpinner()),
                                      tabPanel("Custom Data format", # data table that shows users how to format their data for the custom data format option
                                               tableOutput("lj_format"))),
                          width = 9)
                        )),
                          
### MODULE: Downloads ####             
             navbarMenu("Downloads", # drop down menu for the various download types
                        
                        tabPanel("Raw and Summary Data Output", # This lets users download a merged raw data .csv and an .xlsx containing the bead count repeat summaries from the QC Plots module
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action3", "Generate Raw and Summary Data"), 
                                       downloadButton("saver", "Download Raw and Summary Data"))),
                                   mainPanel())), 
                        
                        
                        tabPanel("Bead QC Plots", # this lets users download high resolution bead count fluctuation plots for each plate faceted by antigen/antibody
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action", "Generate bead count QC plots"),
                                       downloadButton("bqcplots", "Download bead QC plots"))),
                                   mainPanel())), 
                        
                        tabPanel("Bead Count Plate Plots", # this lets users download high resolution bead count 96-well plate plots for each plate faceted by antigen/antibody
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action2", "Generate bead count plate plots"),
                                       downloadButton("plateplotsaver", "Download bead count plate plots"))),
                                   mainPanel())), 
                        
                        tabPanel("MFI Plate Plots", # this lets users download high resolution MFI 96-well plate plots for each plate faceted by antigen/antibody
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action5", "Generate MFI plate plots"),
                                       downloadButton("plateplotsaver_mfi", "Download MFI plate plots"))),
                                   mainPanel())), 
                        
                        tabPanel("MFI Plate Heat Maps", # this lets users download high resolution MFI 96-well plate heat maps for each plate faceted by antigen/antibody
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action7", "Generate MFI plate heat maps"),
                                       downloadButton("heatmapsaver_mfi", "Download MFI plate heat maps"))),
                                   mainPanel())), 
                        
                        tabPanel("Control Tracking Plots", # this lets users download high resolution LJ plots for each antigen/antibody faceted by control
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action4", "Generate Control Tracking Plots"),
                                       downloadButton("ljsaver", "Download Control Tracking Plots"))),
                                   mainPanel())),
                        
                        tabPanel("Control Tracking Summary Results", # this lets users download the LJ summary data tables from the Control Tracking module
                                 sidebarLayout(
                                   sidebarPanel(
                                     verticalLayout(
                                       actionButton("action6", "Generate Control Tracking Summary Results"),
                                       downloadButton("lj_summary_saver", "Download Control Tracking Summary Results"))),
                                   mainPanel())))
             
### UI End ####             
  )))
