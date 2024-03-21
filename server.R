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
library(furrr)
library(microbenchmark)
library(lubridate)
library(profvis)
library(hablar)

#plan(multicore) # Using multicore will greatly improve download speeds if using a Unix system with multiple CPU's
plan(sequential) # Use sequential if running app on Windows or if multiple CPU's are not available 

#### SERVER ####
server <- function(input, output, session) {
  
  #### DISPLAY MANIPULATION ####
  
    # disable these widgets at startup #
  disable("saver")
  disable("bqcplots")
  disable("plateplotsaver")
  disable("plateplotsaver_mfi")
  disable("ljsaver")
  disable("lj_summary_saver")
 
  
  observeEvent(input$file, {
    
    execute_safely(file_merger(), # starts merging the user's xPONENT and/or BPM files as soon as they're uploaded
                   title = "File merging error",
                   message = "An improperly formatted output file was submitted.")
    
    disable("saver")
    disable("bqcplots")
    disable("plateplotsaver")
    disable("plateplotsaver_mfi")
    disable("ljsaver")
  })
  
  observe({ # this hides all of the widgets until the user uploads xPONENT/BPM file(s) or previous shinyMBA datasets
    
    if (is.null(input$file) & is.null(input$prev_datasets)) {
      disable("plate_choice")
    } else {
      enable("plate_choice")
    }

    if (is.null(input$ljfile) & is.null(input$file) & is.null(input$prev_datasets)) {
      disable("action4")
    } else {
      enable("action4")
    }

    if (is.null(input$file) & is.null(input$prev_datasets)) {
      disable("action3")
    } else {
      enable("action3")
    }

    if (is.null(input$file) & is.null(input$prev_datasets)) {
      disable("action")
    } else {
      enable("action")
    }

    if (is.null(input$file) & is.null(input$prev_datasets)) {
      disable("action2")
    } else {
      enable("action2")
    }

    if (is.null(input$file) & is.null(input$prev_datasets)) {
      disable("action5")
    } else {
      enable("action5")
    }
    
    if (is.null(input$samples2)) {
      disable("ab_var2")
    } else {
      enable("ab_var2")
    }
    
    if (is.null(input$samples2) | is.null(input$ag_spec)) {
      disable("samples3")
      disable("ab2")
    } else {
      enable("samples3")
      enable("ab2")
    }
    
    if (input$mfi_plot_switch == "Heat Map") {
      disable("mfi_warning")
    } else {
      enable("mfi_warning")
    }
    
    toggle(id = "lj_xp_dd", condition = !is.null(input$ag_spec) & !is.null(input$samples2))
    toggle(id = "lj_custom_dd", condition = !is.null(input$ljfile))
    
  })
  
  #### DATA MERGING #### 
  
    # file_merger() combines cleans and merges individual xPONENT and BPM files into a single dataset
  file_merger <- reactive({ #this reactive function will take all of the user's uploaded xPONENT + BPM files and merge them into one dataset
    
    req(input$file) # prevents file_merger() from running unless the user has uploaded their files
    
    is_csv <- str_detect(input$file$datapath, ".csv") # detects which files are xPONENT csv's 
    
    nfiles = nrow(input$file) # detects the # of files that have been uploaded
    
    plate_list = list() # creates a blank list that will be used to store each plate's data
    
    for (n in 1:nfiles){ #this loops reads in each of the uploaded files, grabs the data, and stores each plate's data as an element in "plate_list"
      
      if (is_csv[n] == TRUE) { #xPONENT data cleaning starts here
        
        plate <- read.csv(input$file$datapath[n], # reading in the xPONENT file
                          header = FALSE, # prevents R from reading in the first row of data as the header
                          fill = TRUE, # pads shorter rows with blank cells 
                          stringsAsFactors = FALSE, # prevents R from reading in strings as factors
                          blank.lines.skip = FALSE, # VERY IMPORTANT. This keeps blank rows from being cut out when read in.
                          col.names = paste0("V", # naming the columns
                                             seq_len(max(count.fields(input$file$datapath[n], sep = ","), na.rm = TRUE)))) %>% # detects the maxmium number of data fields in the xPONENT file, and then creates a vector to name the columns in sequence
          slice(which(.$V1 == "Results")+2:n()) %>% #detects where the actual results (not metadata) are and slices the dataset there
          select_if(~!(all(. == ""))) # removes any lingering blank columns
        
        plate_metadata <- read.csv(input$file$datapath[n], # reads in the metadata as a separate object
                                   skip = 2, 
                                   nrows = 4, 
                                   header = FALSE, 
                                   stringsAsFactors = FALSE, 
                                   blank.lines.skip = FALSE) %>% select(V1:V3) 
        
        # extracting median FI data #
        med <- plate %>% 
          slice(((which(.$V2 == "Median"))+1):n()) %>% # detects where the median FI data starts and slices the dataset there
          filter(V1 != "") %>%  # removes any blank rows 
          slice(-((which(.$V1 =="DataType:")[1]):n())) # slices out the next DataType segment
        
        colnames(med) <- med[1,] #adjusting the column names
        med <- med[-1,] # ^see above
        
        agab_xp <- med %>% #pulling all antigens/antibodies as a single vector
          select(-c(Location, Sample, `Total Events`)) %>% # removes columns that aren't antigen/antibody names
          colnames() %>% sort() #creates the vector and sorts it alphabetically 
        
        med2 <- med %>% 
          pivot_longer(all_of(agab_xp), #pivots the dataset to a long format
                       names_to = "antigen", #creates and assigns the variable "antigen"
                       values_to = "median_fi") #creates and assigns the variable "median_fi"
        
        # extracting Net MFI data #
        ## see above code for extracting median FI
        netmfi <- plate %>% 
          slice(((which(.$V2 == "Net MFI"))+1):n()) %>% 
          filter(V1 != "") %>% 
          slice(-((which(.$V1 =="DataType:")[1]):n()))
        
        colnames(netmfi) <- netmfi[1,]
        netmfi <- netmfi[-1,]
        
        netmfi2 <- netmfi %>% 
          pivot_longer(all_of(agab_xp),
                       names_to = "antigen",
                       values_to = "net_mfi")
        
        # extracting bead count data #
        ## see above code for extracting median FI
        bcount <- plate %>% 
          slice(((which(.$V2 == "Count"))+1):n()) %>% 
          filter(V1 != "") %>% 
          slice(-((which(.$V1 =="DataType:")[1]):n()))

        colnames(bcount) <- bcount[1,]
        bcount <- bcount[-1,]
        
        
        bcount2 <- bcount %>% 
          pivot_longer(all_of(agab_xp),
                       names_to = "antigen",
                       values_to = "bead_count")
        
        # merging various data types #
        
        merged <- list(med2, netmfi2, bcount2) %>%
          reduce(left_join) 
        
        # assigning metadata #
        merged$date <- plate_metadata$V2[1]
        merged$time <- plate_metadata$V3[1]
        merged$datetime <- mdy_hm(paste(merged$date, merged$time)) # creates a datetime variable that is easily recognized as a date by R
        merged$SN <- plate_metadata$V2[3]
        merged$batch <- plate_metadata$V2[4]
        merged$file_name <- input$file$name[n]
        
        colnames(merged) <- tolower(colnames(merged)) #change all column names to lowercase for easy concatenation
        
        merged$location <- str_replace(merged$location, "\\d+\\(\\d+,", "") %>% # cleaning up well location strings
          str_replace("\\)", "")
        
        plate_list[[n]] <- merged # the cleaned xPONENT file is stored as an element in "plate_list"
        
      } else { # BPM data cleaning starts here
        
        validate(need(getSheetNames(input$file$datapath[n]) %in% "FI", 
                      message = paste("BPM file", input$file$name[n], "is missing FI data")),
                 need(getSheetNames(input$file$datapath[n]) %in% "FI - Bkgd", 
                      message = paste("BPM file", input$file$name[n], "is missing FI - bkgd data")),
                 need(getSheetNames(input$file$datapath[n]) %in% "Bead Count", 
                      message = paste("BPM file", input$file$name[n], "is missing Bead Count data")),
                 errorClass = "merge")
        
        ### reading in net mfi data ### 
        bplex1 <- read.xlsx(input$file$datapath[n], 
                            sheet = "FI - Bkgd",
                            colNames = FALSE,
                            skipEmptyRows = FALSE)
        
        sets_l <- length(which(bplex1$X1 == "Type")) #vectorizing the number of appearances of the string "Type" 
        
        bplex2 <- bplex1 %>% 
          slice(ifelse(sets_l == 2, 
                       which(.$X1 == "Type")[2], #pulls data from the second set of data if "Type" appears twice
                       which(.$X1 == "Type")) # pulls data from the first set of data if "Type" doesn't appear twice
                +(-1):n())
        
        ### cleaning steps ###
        run_date1 <- bplex1[2,1] %>% 
          str_replace("Acquisition Date: ", "") %>% #cleaning the date string to match the xPONENT format
          str_replace( "\\,.+", "")
        
        run_time1 <- bplex1[2,1] %>%
          str_replace(".+\\,", "") %>% # same as above ^
          str_trim()
        
        sn1 <- bplex1[3,1] %>% 
          str_replace("Reader Serial Number: ", "") # cleaning the serial number string to match xPONENT format
        
        bplex2[1, 1:3] <- bplex2[2, 1:3] 
        
        colnames(bplex2) <- bplex2[1,]
        
        bplex2 <- bplex2[-1:-2,] %>% 
          filter(!is.na(Well)) # remove well data that is read in as NA
        
        validate(need(colnames(bplex2) %in% "Description", paste("BPM file", input$file$name[n], "is missing sample IDs")),
                 errorClass = "merge") # generates an error message if the bpm file does not have a "Description" column (i.e. sample IDs)

        #       #       #       #       #
        
        bplex_bwell <- bplex2 %>% # cleaning well location strings to match xPONENT format
          filter(str_detect(.$Well, ",") == TRUE) %>% 
          mutate(Well = str_trim(str_replace(Well, "\\,\\w+", "")))
        
        bplex3 <- bind_rows(bplex_bwell, bplex2) %>% # same as above ^
          mutate(Well = str_trim(str_replace(Well, "\\w+\\,", "")))
        
        
        agab_bpm <- bplex3 %>% # pulling all antigens/antibodies as a single vector
          select(-c(Type, Well, Description)) %>% # removing any columns that aren't antigen/antibody names
          colnames() %>% sort() # creates the vector and sorts is alphabetically
        
        bplex4 <- bplex3 %>% 
          pivot_longer(all_of(agab_bpm),
                       names_to = "antigen",
                       values_to = "net_mfi") %>% 
          mutate(antigen = str_trim(str_replace(antigen, "\\(\\d+\\)", "" )), # cleaning antigen strings to match xPONENT format
                 Description = str_trim(str_replace(Description, "\\w+\\s\\w+\\:", "" )), # cleaning Description strings to match xPONENT format
                 date = run_date1,
                 time = run_time1,
                 sn = sn1) %>% 
          select(location = Well,
                 sample = Description,
                 everything(),
                 -(Type))
        
        ### reading in median FI ###
        
        
        bplex_mfi1 <- read.xlsx(input$file$datapath[n], 
                               sheet = "FI",
                               colNames = FALSE,
                               skipEmptyRows = FALSE)
        
        sets_mfil <- length(which(bplex_mfi1$X1 == "Type")) # see net mfi cleaning section
        
        bplex_mfi2 <- bplex_mfi1 %>% 
          slice(ifelse(sets_mfil == 2, 
                       which(.$X1 == "Type")[2],
                       which(.$X1 == "Type"))
                +(-1):n())
        
        ### cleaning steps ###
        
        bplex_mfi2[1, 1:3] <- bplex_mfi2[2, 1:3]
        
        colnames(bplex_mfi2) <- bplex_mfi2[1,]
        
        bplex_mfi2 <- bplex_mfi2[-1:-2,] %>% 
          filter(!is.na(Well))
        
        #       #       #       #       #
        
        bplex_mfiwell <- bplex_mfi2 %>% 
          filter(str_detect(.$Well, ",") == TRUE) %>% 
          mutate(Well = str_trim(str_replace(Well, "\\,\\w+", "")))
        
        bplex_mfi3 <- bind_rows(bplex_mfiwell, bplex_mfi2) %>% 
          mutate(Well = str_trim(str_replace(Well, "\\w+\\,", "")))
        
        
        bplex_mfi4 <- bplex_mfi3 %>% 
          pivot_longer(all_of(agab_bpm), #agab created during net_mfi extraction
                       names_to = "antigen",
                       values_to = "median_fi") %>% 
          mutate(antigen = str_trim(str_replace(antigen, "\\(\\d+\\)", "" )),
                 Description = str_trim(str_replace(Description, "\\w+\\s\\w+\\:", "" )))  %>% 
          select(location = Well,
                 sample = Description,
                 everything(),
                 -(Type))
        
        #  #    #      #      #
        ### reading in bead count data ###
        
        bplex_bc1 <- read.xlsx(input$file$datapath[n], 
                               sheet = "Bead Count",
                               colNames = FALSE,
                               skipEmptyRows = FALSE)
        
        sets_bcl <- length(which(bplex_bc1$X1 == "Type")) # see net mfi cleaning section
        
        bplex_bc2 <- bplex_bc1 %>% 
          slice(ifelse(sets_bcl == 2, 
                       which(.$X1 == "Type")[2],
                       which(.$X1 == "Type"))
                +(-1):n())
        
        ### cleaning steps ###
        
        bplex_bc2[1, 1:3] <- bplex_bc2[2, 1:3]
        
        colnames(bplex_bc2) <- bplex_bc2[1,]
        
        bplex_bc2 <- bplex_bc2[-1:-2,] %>% 
          filter(!is.na(Well))
        
        #       #       #       #       #
        
        bplex_bcwell <- bplex_bc2 %>% 
          filter(str_detect(.$Well, ",") == TRUE) %>% 
          mutate(Well = str_trim(str_replace(Well, "\\,\\w+", "")))
        
        bplex_bc3 <- bind_rows(bplex_bcwell, bplex_bc2) %>% 
          mutate(Well = str_trim(str_replace(Well, "\\w+\\,", "")))
        
        
        bplex_bc4 <- bplex_bc3 %>% 
          pivot_longer(all_of(agab_bpm), #agab created during net_mfi extraction
                       names_to = "antigen",
                       values_to = "bead_count") %>% 
          mutate(antigen = str_trim(str_replace(antigen, "\\(\\d+\\)", "" )),
                 Description = str_trim(str_replace(Description, "\\w+\\s\\w+\\:", "" )))  %>% 
          select(location = Well,
                 sample = Description,
                 everything(),
                 -(Type))
        
        ### merging BPM data ###
        
         bpm_merge <- list(bplex4, bplex_mfi4, bplex_bc4) %>% 
          reduce(left_join) %>% 
          mutate(datetime = dmy_hm(paste(date, time)), 
                 file_name = input$file$name[n],
                 batch = str_replace(file_name, "\\_\\d{8}\\_.+", ""), # removing trailing digits/characters from batch name
                 `total events` = NA) %>% 
          select(location:net_mfi,
                 median_fi,
                 bead_count,
                 `total events`,
                 datetime,
                 date:sn, 
                 batch,
                 file_name)
        
        plate_list[[n]] <- bpm_merge
        
      }
      
    } #End of the loop. All xPONENT and BPM outputs have been cleaned and inserted as individual elements in "plate_list"
    
    xp_bpm_data <- do.call(rbind, plate_list) %>% # this binds each element (the cleaned xPONENT and/or BPM data) of "plate_list" into one dataframe
                   mutate(date = as.factor(date), # converting these variables to factors makes certain analyses easier to perform later
                          sn = as.factor(sn),
                          batch = as.factor(batch),
                          file_name = as.factor(file_name),
                          sample = as.factor(str_trim(sample, "both")),
                          antigen = as.factor(tolower(.$antigen)),
                          median_fi = as.numeric(median_fi), # this makes sure that R treats these as numerical variables
                          bead_count = as.numeric(bead_count),
                          `total events` = as.numeric(`total events`),
                          net_mfi = as.numeric(net_mfi))%>%  
                   mutate_at(c("median_fi", "bead_count", "net_mfi"), ~ifelse(is.nan(.), NA, .)) %>% # changing NaN values to NA
                   mutate(bead_count = ifelse(is.na(bead_count), 0, bead_count))
    
    xp_bpm_data
  }) # end of the file_merger() function
  
    # prev_file_merger() reads in previous shinyMBA datasets and will merge them together if >1 are uploaded 
  prev_file_merger <- reactive({
    
    req(input$prev_datasets) # prevents prev_file_merger() from running unless the user has uploaded a previous shinyMBA dataset(s)
    
    data1 <- map_dfr(input$prev_datasets$datapath, read.csv, # reads in previous shinyMBA datasets and merges them into a single dataset
                     check.names = FALSE,
                     colClasses = c("factor", "factor", "numeric", "factor", "numeric", "numeric",
                                    "numeric", "numeric","numeric", "numeric","numeric", "numeric",
                                    "character","character","character", "factor", "factor", "factor")) 
    
    data2 <- data1 %>% 
      select(-(ut:over_ut)) %>%  #removing bead count dummy variables to make merging with individual xPONENT/BPM files easier
      mutate(datetime = as_datetime(datetime))
    
    data2
  })
    # master_data() dictates the structure of the central dataset that will be used for analysis in the app
  master_data <- reactive({
    
    if (input$master_source == "Only xPONENT/BPM") {
      master <- file_merger()
    } else if (input$master_source == "Only shinyMBA datasets") {
      master <- prev_file_merger()
    } else {
      master <- bind_rows(prev_file_merger(), file_merger())
    }
    
    master
  })
  
  beadcount <- reactive({ #creates a separate dataset with just bead count data
    
    master_data() %>% 
      select(-c(median_fi, net_mfi)) %>%
      pivot_wider(names_from = antigen,
                  values_from = bead_count)
  })
  
  #### CONTROL TRACKING DATA PREP ####
  lj_xp_data <- reactive({ # this reactive function performs some preliminary data cleaning on xPONENT/BPM data before proceeding with the control tracking analyses

    data1 <- master_data()
    
    lj_preclean <- data1 %>%
      select(location:median_fi, date:batch, datetime, -c(`total events`, date, time)) %>% # selecting only the needed variables
      group_by(sample, antigen) %>% # groups observations by sample ID's and antigens
      arrange(datetime) # orders data by datetime (earliest to latest)
    
    sum1 <- lj_preclean %>% # acquires the number of observations for each sample ID + antigen group
      summarize(n = n()) 
    
    rejoin1 <- left_join(lj_preclean, sum1) # adding observation number metadata back to the initial dataset
    
    lj_clean <- rejoin1 %>%
      ungroup() %>% #removing sample ID + antigen grouping
      filter(n >= 2) %>% # removes any data with a single observation (control tracking analyses cannot be performed unless there are >1 observations)
      select(-(n)) %>% # keep all variables except observation number
      mutate(antigen = as.factor(as.character(antigen)), # refactoring the antigen and sample variables
             sample = as.factor(as.character(sample))) 
    
    lj_clean
  })  
  
  lj_data <- reactive({ #this function cleans the user's data into a format needed to perform Levey-Jennings (LJ) analyses
      
      if (input$lj_choice == "Formatted Custom Data"){ # custom formatted data cleaning starts here
        
        req(input$ljfile)
        
        data1 <- read.csv(input$ljfile$datapath[1], # reading in custom data
                          header = TRUE,
                          stringsAsFactors = FALSE)
        
        colnames(data1) <- c("observation", "antigen", "sample", "median_fi") #renaming the column names in case the user doesn't name them correctly
        
        data2 <- data1 %>% 
          mutate(antigen = as.factor(antigen), #factoring these variables
                 sample = as.factor(sample),
                 observation = as.character(observation))
        
        data2$observation <- str_replace(data2$observation, "\\(\\d+\\)", "") %>% # cleaning observation strings in case they have excessive characters
          str_trim(.) 
        
        lj_clean <- data2 %>% # data cleaning for use in the control tracking analyses
          mutate(observation = as.numeric(observation)) %>% 
          group_by(sample, antigen) %>% 
          arrange(observation, .by_group = TRUE) %>% # orders data by observation within antigen-sample groups
          mutate(observation = as.character(observation)) %>% 
          nest() # nesting the dataframe so iterative map functions can be used later
        # custom formatted data cleaning function ends here
        
      } else { # xPONENT/Bio-Plex Manager cleaning starts here

       if (input$lj_avg == FALSE) { # NOT averaging sample duplicates for each plate
          lj_clean <- lj_xp_data() %>% 
            filter(sample %in% input$samples2, # subset data to include specified controls
                   antigen %in% input$ag_spec) %>% # subset data to include specified antigens/antibodies
            mutate(antigen = as.factor(as.character(antigen)), # refactors antigen variable
                   sample = as.factor(as.character(sample))) %>% # refactors sample variable
          group_by(sample, antigen) %>% # grouping sample ID's and antigens before chronologically ordering the data
          arrange(datetime) %>% # arranging the data in chronological order (this is required to properly conduct control tracking analyses)
          nest() # nesting the dataframe so iterative map functions can be used later
       } else { # averaging sample duplicates for each plate
          lj_clean <- lj_xp_data() %>% 
            filter(sample %in% input$samples2, # subset data to include specified controls
                   antigen %in% input$ag_spec) %>% # subset data to include specified antigens/antibodies
            mutate(antigen = as.factor(as.character(antigen)), # refactors antigen variable
                   sample = as.factor(as.character(sample))) %>% # refactors sample variable
            group_by(sample, antigen, batch, datetime) %>% # groups sample duplicates for each plate together
            summarize(reps = n(),
                      avg_mfi = mean(median_fi, na.rm = TRUE)) %>% # averages sample duplicates for each plate
            select(everything(), median_fi = avg_mfi) %>% 
            group_by(sample, antigen) %>% # grouping sample ID's and antigens before chronologically ordering the data
            arrange(datetime) %>% # arranging the data in chronological order (this is required to properly conduct control tracking analyses)
            nest() # nesting the dataframe so iterative map functions can be used later
       }
    } # xPONENT/Bio-Plex Manager cleaning ends here
    
      lj_cleanup <- function(df) { # this function will be used to conduct process control chart analyses on either the custom or xPONENT/BPM data that was cleaned earlier 
        
        if (input$chart_type == "Shewhart") {
          d1 <- slice(df, 1:input$shew_ref) # Shewhart analyses calculates CL's from a subset of data points selected by the user
        } else {
          d1 <- df # LJ and I-MR analyses calculates CL's from all data points
        }
        
        d1 <- d1 %>% 
          mutate(mr = abs(median_fi - lag(median_fi))) # creates a variable for the point-to-point moving range
        
        mfi_avg <- mean(d1$median_fi, na.rm = TRUE) #mfi average
        mr_avg <- mean(d1$mr, na.rm = TRUE) # moving range average
        mr_sd <- mr_avg / 1.128 # I-MR standard deviation calculated using n = 2 observations constant of 1.128
        
        if (input$chart_type == "Moving Range"){ # confidence limit determination for I-MR charts
          
          lim3 <- list(ucl = mfi_avg + (3 * mr_sd), # I-MR 3-sigma limits
                       cl  = mfi_avg,
                       lcl = mfi_avg - (3 * mr_sd))
          
          lim2 <- list(ucl = mfi_avg + (2 * mr_sd), # I-MR 2-sigma limits
                       cl  = mfi_avg,
                       lcl = mfi_avg - (2 * mr_sd))
          
          lim1 <- list(ucl = mfi_avg + mr_sd, # I-MR 1-sigma limits
                       cl  = mfi_avg,
                       lcl = mfi_avg - mr_sd)
        } else { # LJ and Shewhart confidence limit calculations
          
        lim3 <- CalculateLimits(d1$median_fi, 
                                type = "i") # calculates the 3-sigma limits
        lim2 <- CalculateLimits(d1$median_fi, 
                                controlLimitDistance = 2, 
                                type = "i") # calculates the 2-sigma limits
        lim1 <- CalculateLimits(d1$median_fi, 
                                controlLimitDistance = 1,
                                type = "i") # calculates the 1-sigma limits
        }
        
        d2 <- EvaluateRules(df$median_fi, # Control tracking analysis using Rspc package. See https://cran.r-project.org/web/packages/Rspc/index.html
                            type = "i",
                            controlLimitDistance = 3,
                            ucl = round(lim3[[1]], 0),
                            cl = round(lim3[[2]], 0),
                            lcl = round(lim3[[3]], 0)) # LJ analyses 
        
        d3 <- d2 %>% # cleaning up resulting control tracking analyses for use in upcoming plots + tables
          mutate(obs = row_number(),
                 mr = abs(x - lag(x)),
                 mr_avg = mr_avg,
                 mr_ucl = (3.267 * mr_avg),
                 cl = round(lim3[[2]], 0), #calculating confidence limits
                 ucl3 = round(lim3[[1]], 0),
                 lcl3 = round(lim3[[3]], 0),
                 ucl2 = round(lim2[[1]], 0),
                 lcl2 = round(lim2[[3]], 0),
                 ucl1 = round(lim1[[1]], 0),
                 lcl1 = round(lim1[[3]], 0),
                 # r1-r8 will assess whether any of the Nelson rules were violated
                 r1 = ifelse(Rule1 == 1, "1,", ""), 
                 r2 = ifelse(Rule2 == 1, "2,", ""),
                 r3 = ifelse(Rule3 == 1, "3,", ""),
                 r4 = ifelse(Rule4 == 1, "4,", ""),
                 r5 = ifelse(Rule5 == 1, "5,", ""),
                 r6 = ifelse(Rule6 == 1, "6,", ""),
                 r7 = ifelse(Rule7 == 1, "7,", ""),
                 r8 = ifelse(Rule8 == 1, "8", ""),
                 status = case_when(
                   rowSums(select(., contains("Rule"))) > 0 & is.na(x) == FALSE  ~ lj_status()[2], # criteria for out of control points
                   rowSums(select(., contains("Rule"))) == 0 & is.na(x) == FALSE ~ lj_status()[1]),
                 vio = ifelse(status == lj_status()[1], NA, paste(r1,r2,r3,r4,r5,r6,r7,r8, sep = ""))) # concatenating all of the violated rules as a single string
        
        d3
      }
      
      lj_final <- lj_clean %>% 
        mutate(lj = map(data, lj_cleanup)) %>% # map lj_cleanup() to iteratively run LJ analyses for each sample+antigen combination. DO NOT USE future_map
        unnest(cols = c(data, lj)) # unnesting to restore the dataframe
      
      lj_final
    })
  
  lj_data_23sd <- reactive({ # creates a separate dataset for CL2-CL3 violation types
    
    data1 <- lj_data() %>%  
      mutate(status = case_when( # changing the status variable to fit the 2CL-3CL analysis
        Rule1 == 1                     ~ lj_status()[2], # criteria for "Out of Control" datapoints
        Rule1 == 0 & median_fi >= ucl2 ~ lj_status()[3], # criteria for "2CL Warning" datapoints
        Rule1 == 0 & median_fi <= lcl2 ~ lj_status()[3],
        Rule1 == 0 & ucl2 >= median_fi & median_fi >= lcl2 ~ lj_status()[1])) # criteria for "In Control" datapoints
    
    data1
  })
  
  mr_data <- reactive({ # separate dataset that be used for the moving range plot (NOT individuals)
    
    data1 <- lj_data() %>% 
      mutate(status = ifelse(mr >= mr_ucl, # points are defined as "Out of Control" if moving range is higher than the upper confidence limit
                             lj_status()[2],
                             lj_status()[1]))
    data1
  })
    
  #### DYNAMIC VECTORS (QC PLOTS) ####
  
  ag_names <- reactive({ # this reactive function will generate a vector containing all of the antigens/antibodies tested in the plate chosen by the user
    
    req(input$plates)
    
    subset1 <- master_data() %>% 
      filter(batch == input$plates) %>%
      droplevels()
    
    levels(subset1$antigen) %>% sort()
  }) 
  
  ag_names_all <- reactive({ # this reactive function will generate a vector containing all of the antigens/antibodies tested in the user's xPONENT and/or BPM  file(s)
    levels(master_data()$antigen) %>% sort()
  })
  
  output$ab_var <- renderUI({ # BEAD COUNT QC: this generates a dynamic drop-down list of potential antigens/antibodies for the user to choose from
    selectInput("ab", "Choose an antigen/antibody", choices = ag_names())
  })
  
  output$ab_var_mfi <- renderUI({ # MFI QC: this generates a dynamic drop-down list of potential antigens/antibodies for the user to choose from
    selectInput("ab", "Choose an antigen/antibody", choices = ag_names())
  })
  
  batch_names <- reactive ({ # this reactive function will generate a vector containing all of the plate names from the user's xPONENT/BPM file(s)
    
    b_names1 <- master_data()
    
    levels(b_names1$batch) %>% str_sort(., numeric = TRUE) # creates a vector of the batch names (AKA plate names)
  })
  
  output$batches <- renderUI({ #BEAD QC: this generates a dynamic drop-down list of potential plate names for the user to choose from
    selectInput("plates", "Choose a plate", choices = batch_names())
  })
  
  output$batches_mfi <- renderUI({ #MFI QC: this generates a dynamic drop-down list of potential plate names for the user to choose from
    selectInput("plates", "Choose a plate", choices = batch_names())
  })
  
  mfi_status <- reactive({ # vector containing status factors for the MFI module analyses
    
    mfi_levels <- c("< Background or No Signal", paste("\u2265", "Background"), "Background Well", paste0("Warning ", "(< ", input$mfi_warning, " MFI)"))#, "No Signal")
    
    mfi_levels
  })
  
  #### DYNAMIC VECTORS (CONTROL TRACKING) ####
  
  lj_antigens <- reactive({ # creates a dynamic vector of the antigen names for use in LJ analyses
    
    data1 <- lj_data()
    
    levels(data1$antigen) %>% sort()
  })
  
  lj_controls <- reactive({ # creates a dynamic vector of the sample/control names for use in LJ analyses
    
    data1 <- lj_data()
    
    levels(data1$sample) %>% sort()
  })
  
   observe({
    if (!is.null(input$ljfile) & input$lj_choice == "Formatted Custom Data") { #updates drop-down lists when custom datasets are used
      updateSelectInput(session = session, inputId = "ab2_custom", label = "Choose an antigen/antibody", choices = lj_antigens())
      updateSelectInput(session = session, inputId = "samples2_custom", label = "Choose a control", choices = lj_controls())
     }
   })
  
  
  observe({ # tells the app to update the "samples2" and "ab2" UI widgets when the user chooses xPONENT/BPM files as their data source
    
    controls <- levels(lj_xp_data()$sample) %>% sort()
    
    ag <- levels(lj_xp_data()$antigen) %>% sort
    
    if (!is.null(input$file) | !is.null(input$prev_datasets) & input$lj_choice == "xPONENT/Bio-Plex Manager") {
       updatePickerInput(session = session, inputId = "samples2", label = "Specify controls", choices = controls)
       updatePickerInput(session = session, inputId = "ag_spec", label = "Specify antigens/antibodies", choices = ag)
    }
  })
   
  observe({ # tells the app to update the "samples3" UI widget once the user specifies their controls

   
    if (!is.null(input$samples2) & !is.null(input$ag_spec)) {
      updateSelectInput(session = session, inputId = "samples3", label = "Choose a control", choices = input$samples2)
      updateSelectInput(session = session, inputId = "ab2", label = "Choose an antigen/antibody", choices = input$ag_spec)
      
    }
  })
  
  lj_status <- reactive({ # character vector used for control tracking flagging 
    
    statuses <- c("In control", "Out of control", "CL2 Warning")
    statuses
  })
  
  #### MERGED AND SUMMARY DATA OUTPUT  ####
  
  raw_data <- reactive({ #cleaning up the raw data that will be downloaded by the user
    
    raw1 <- left_join(master_data(), all_plots_data()) %>% 
      select(-c(observation, mean, per_under_ut, per_under_lt), status2 = utlt) %>%
      select(location:net_mfi, median_fi, bead_count, ut:over_ut, datetime, date:file_name) 
    
    raw1
  })
  
  raw_data_wide <- reactive({ # reformatting raw data in a wide format for users to download
    
    data1 <- raw_data()
    
    bc_wide1 <- data1 %>% # this creates a wide dataframe for the bead count data
      select(-c(median_fi, net_mfi, ut:over_ut)) %>% 
      pivot_wider(names_from = "antigen",
                  values_from = "bead_count",
                  names_prefix = "bead_count_")
    
    mfi_wide1 <- data1 %>% # this creates a wide dataframe for the mfi data
      select(-c(median_fi, bead_count, ut:over_ut)) %>% 
      mutate(net_mfi = round(net_mfi, 0)) %>% 
      pivot_wider(names_from = "antigen",
                  values_from = "net_mfi",
                  names_prefix = "mfi_")
    
    wide_merge <- left_join(bc_wide1, mfi_wide1) #merging the wide bead count and mfi data into a single dataframe
    
    wide_merge
  })
  
  psum_ut <- reactive({ #this reactive function creates a % < the upper threshold summary dataframe for wells 
    
    psum_ut1 <- beadcount()
    
    psum_ut2 <- psum_ut1 %>% 
      mutate_at(vars(ag_names_all()),
                list(~ ifelse( . > input$uthresh, 0, 1))) %>% 
      select(ag_names_all(), batch, everything()) %>% 
      mutate(sumrow = rowSums(.[1:length(ag_names_all())], na.rm = TRUE)) %>% 
      select(batch, ag_names_all(), -(location:sumrow)) %>% 
      group_by(batch) %>% 
      summarize_at(vars(ag_names_all()), mean) %>% 
      mutate_if(is.numeric, percent)
    
    psum_ut2
  })
  
  psum_lt <- reactive({ #this reactive function creates a % < the lower threshold summary dataframe for wells 
    
    psum_lt1 <- beadcount()
    
    psum_lt2 <- psum_lt1 %>% 
      mutate_at(vars(ag_names_all()),
                list(~ ifelse( . > input$lthresh, 0, 1))) %>% 
      select(ag_names_all(), batch, everything()) %>% 
      mutate(sumrow = rowSums(.[1:length(ag_names_all())], na.rm = TRUE)) %>% 
      select(batch, ag_names_all(), -(location:sumrow)) %>% 
      group_by(batch) %>% 
      summarize_at(vars(ag_names_all()), mean) %>% 
      mutate_if(is.numeric, ~percent(., accuracy = 1.0)) 
    
    psum_lt2
  })
  
  bc_summary <- reactive({ # generates a dataframe summarizing batch + antigen bead count failures 
    
    data1 <- all_plots_data() %>% 
      select(antigen, batch, location, bead_count, under_lt, under_ut) %>% 
      filter(!is.na(bead_count)) %>% 
      group_by(batch, antigen)
    
    summary1 <- data1 %>% 
      summarize(wells = n(),
                fails = ifelse(input$utlt_platefail == "< Upper threshold", sum(under_ut), sum(under_lt)), # user sets thresholds
                per_fail = fails/wells,
                status = as.character(ifelse(per_fail >= (input$bc_platefail / 100), "Repeat", "Safe")), # user sets threshold
                status_flag = ifelse(status == "Repeat", 1, 0))
    
    summary1
  })
  
  bc_summary2 <- reactive({ # generates a dataframe summarizing batch bead count failures
    
    data1 <- bc_summary() %>% 
      group_by(batch) %>% 
      nest() # nesting the data by batch so the following functions can be iteratively applied
    
    batch_fail <- function (df) { # function that will subset antigen bead count fails within a plate and then concatenate the failed antigen names together
      
      fdata1 <- df
      
      fdata2  <- fdata1 %>%  
        filter(status_flag == 1) %>% 
        mutate(antigen2 = as.factor(as.character(antigen)))
      
      paste(fdata2$antigen, collapse = ", ")
    }
    
    n_badag <- function (df) { # function that sums the number of antigens within a plate whose bead count failed 
      
      ffdata1 <- df
      
      sum(ffdata1$status_flag)
    }
    
    data2 <- data1 %>% 
      mutate(badags = map(data, batch_fail), 
             n_failed_ag = map(data, n_badag)) %>% 
      select(batch, n_failed_ag, badags) %>% 
      rename(Batch = batch, `n Failed Antigens` = n_failed_ag, `Failed Antigens` = badags) %>% # renaming variable so they look nice in a table
      unnest() # remove nesting to return a normal dataframe
    
    data2
  })
  
  #### CONTROL TRACKING SUMMARY DATA OUTPUT ####
  
  nelson_summary <- reactive({ #creates a summary statistics table for Nelson rule violations
   
   if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
    controls <- input$samples2
    ag <- input$ag_spec
   } else {
     controls <- lj_controls()
     ag <- lj_antigens()
   }
    
    data1 <- lj_data()
    
    data2 <- data1 %>% 
      filter(sample %in% controls,
             antigen %in% ag) %>% # filters out samples that are not controls
      mutate(ooc = ifelse(status == lj_status()[2], 1, 0)) %>% 
      group_by(sample, antigen) # summary statistics applied to each control-antigen/antibody combination
    
    summary1 <- summarize(data2, 
                          n = n(),
                          n_missing = sum(is.na(median_fi)),
                          mean_mfi = round(mean(median_fi, na.rm = TRUE), 1),
                          mean_mr = round(mean(mr, na.rm = TRUE), 1),
                          mr_ucl = round(mean(mr_ucl, na.rm = TRUE), 1),
                          median_mfi = round(median(median_fi, na.rm = TRUE), 1),
                          max_mfi = max(median_fi, na.rm = TRUE),
                          min_mfi = min(median_fi, na.rm = TRUE),
                          sd_mfi = round(sd(median_fi, na.rm = TRUE), 2),
                          upper_cl3 = mean(ucl3),
                          lower_cl3 = mean(lcl3),
                          upper_cl2 = mean(ucl2),
                          lower_cl2 = mean(lcl2),
                          upper_cl1 = mean(ucl1),
                          lower_cl1 = mean(lcl1),
                          out_of_control = percent((sum(ooc)/n), accuracy = 0.01),
                          total_violations = sum(Rule1,
                                                 Rule2,
                                                 Rule3,
                                                 Rule4,
                                                 Rule5,
                                                 Rule6,
                                                 Rule7,
                                                 Rule8),
                          rule1_violations = sum(Rule1),
                          rule2_violations = sum(Rule2),
                          rule3_violations = sum(Rule3),
                          rule4_violations = sum(Rule4),
                          rule5_violations = sum(Rule5),
                          rule6_violations = sum(Rule6),
                          rule7_violations = sum(Rule7),
                          rule8_violations = sum(Rule8))
    
    summary1
  })
  
  nelson_summary_small <- reactive({ # creates a data table for control tracking flagging results using the Nelson Rules
    
    if (input$lj_choice == "Formatted Custom Data"){ # custom formatted data
      table1 <- lj_data() %>% 
        filter(sample %in% lj_controls()) %>% 
        select(obs, 
               Control = sample, 
               Antigen = antigen, 
               MFI = median_fi, 
               MR = mr,
               Status = status, 
               Violations = vio)
    } else { # xPONENT/BPM data
        table1 <- lj_data() %>% 
          select(Obs = obs, 
                 Date = datetime, 
                 Control = sample, 
                 Antigen = antigen, 
                 Batch = batch, 
                 MFI = median_fi,
                 MR = mr,
                 Status = status, 
                 Violations = vio)
    }
    
    table1
  })
  
  cl23_summary <- reactive({ #creates a summary statistics table for CL2-CL3 violations
    
    if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
      controls <- input$samples2
      ag <- input$ag_spec
    } else {
      controls <- lj_controls()
      ag <- lj_antigens()
    }
    
    data1 <- lj_data_23sd() %>%
      filter(sample %in% controls,
             antigen %in% ag) # filters out samples that are not controls
    
    data2 <- data1 %>% 
      mutate(ooc3sd = ifelse(status == lj_status()[2], 1, 0),
             ooc2sd = ifelse(status == lj_status()[3], 1, 0)) %>% 
      group_by(sample, antigen) # summary statistics applied to each control-antigen/antibody combination
    
    summary1 <- summarize(data2, 
                                 n = n(),
                                 n_missing = sum(is.na(median_fi)),
                                 mean_mfi = round(mean(median_fi, na.rm = TRUE), 1),
                                 mean_mr = round(mean(mr, na.rm = TRUE), 1),
                                 mr_ucl = round(mean(mr_ucl, na.rm = TRUE), 1),
                                 median_mfi = round(median(median_fi, na.rm = TRUE), 1),
                                 max_mfi = max(median_fi, na.rm = TRUE),
                                 min_mfi = min(median_fi, na.rm = TRUE),
                                 sd_mfi = round(sd(median_fi, na.rm = TRUE), 2),
                                 upper_cl3 = mean(ucl3),
                                 lower_cl3 = mean(lcl3),
                                 upper_cl2 = mean(ucl2),
                                 lower_cl2 = mean(lcl2),
                                 n_outside_CL3 = sum(ooc3sd, na.rm = TRUE),
                                 per_outside_CL3 = percent((n_outside_CL3/n), accuracy = 0.01),
                                 n_between_CL2_CL3 = sum(ooc2sd, na.rm = TRUE),
                                 per_between_CL2_CL3 = percent((n_between_CL2_CL3/n), accuracy = 0.01))
    
    summary1
  })
  
  cl23_summary_small <- reactive({ # creates a data table for control tracking flagging results using the 2CL-3CL flagging criteria
    
    if (input$lj_choice == "Formatted Custom Data"){ #custom dataset
      table1 <- lj_data_23sd() %>% 
        filter(sample %in% lj_controls()) %>% 
        select(obs, 
               Control = sample, 
               Antigen = antigen, 
               MFI = median_fi,
               MR = mr,
               Status = status)
      
    } else { #xPONENT/BPM data
      table1 <- lj_data_23sd() %>% 
        select(Obs = obs, 
               Date = datetime, 
               Control = sample, 
               Antigen = antigen, 
               Batch = batch, 
               MFI = median_fi,
               MR = mr,
               Status = status)
    }
    
    table1
  })
  
  mr_summary_small <- reactive({ # control tracking data table specific for the point-to-point moving range chart (NOT individuals)
    
    if (input$lj_choice == "Formatted Custom Data"){
      
      table1 <- mr_data() %>% 
        filter(sample %in% lj_controls()) %>% 
        select(obs, 
               Control = sample, 
               Antigen = antigen, 
               MFI = median_fi,
               MR = mr,
               Status = status)
      
    } else {
      
      table1 <- mr_data() %>% 
        select(Obs = obs, 
               Date = datetime, 
               Control = sample, 
               Antigen = antigen, 
               Batch = batch, 
               MFI = median_fi,
               MR = mr,
               Status = status)
    }
    
    table1
  })

  #### DATA MANIPULATION FOR QC PLOTS ####
  all_plots_data <- reactive({ # this reactive object creates a dataset that will be used for generating plots in the app
    
    plots1 <- beadcount()
    
    plots2 <- plots1 %>% # modifying the bead count data so that it can be easily plotted later
      select(any_of(ag_names_all()), everything()) %>% 
      pivot_longer(ag_names_all(), # returning data to a long format
                   names_to = "antigen",
                   values_to = "bead_count") %>% 
      filter(!is.na(bead_count)) %>% # removing any missing data
      group_by(antigen, batch) %>% # grouping so that the following variable creation is performed for each batch + antigen combination
      mutate(observation = 1:n(),
             mean = round(mean(bead_count, na.rm = TRUE), 0), # mean bead count
             ut = input$uthresh, # upper bead count threshold (set by user)
             lt = input$lthresh, # lower bead count threshold (set by user)
             under_ut     = ifelse(bead_count <  ut,  1, 0), # flags wells that are < upper threshold
             under_lt     = ifelse(bead_count <  lt,  1, 0), # flags wells that are < lower threshold
             over_ut      = ifelse(bead_count >= ut, 1, 0), # flags wells that are > upper threshold
             per_under_ut = percent(sum(under_ut) / n(), accuracy = 0.1), # calculates the percentage of wells that are < upper threshold for that batch + antigen combination
             per_under_lt = percent(sum(under_lt) / n(), accuracy = 0.1),# calculates the percentage of wells that are < lower threshold for that batch + antigen combination
             utlt         = ifelse(bead_count < lt, paste("<", input$lthresh, "beads", sep = " "), # creates a categorical variable that describes the bead count for a well
                                   ifelse((ut > bead_count) & (bead_count >= lt), paste("<", input$uthresh, "beads", sep = " "), 
                                          paste("\u2265", input$uthresh, "beads", sep = " ")))) %>%
      arrange(antigen, batch) %>% 
      ungroup() %>% 
      mutate(antigen = as.factor(antigen)) # refactoring antigen variable
    
    # MFI plots data cleaning #
    mfi_data1 <- master_data() %>% 
      group_by(batch, antigen)
    
    bkgd_max <- mfi_data1 %>% 
      filter(sample == "Background0") %>% 
      summarize(max_bkgd = max(s(net_mfi))) # creates a variable of the highest background MFI value for each plate + antigen combination (NA's are not included unless all background wells are NA)
    
    merged <- list(mfi_data1, bkgd_max, plots2) %>% 
      reduce(left_join) %>% # joining the bead count, mfi, and bkgd value datasets
      ungroup()
    
    merged2 <- merged %>% # this object defines the MFI status for every well
      mutate(mfi_well_status = as.factor(case_when( 

           sample == "Background0"                                                    ~ mfi_status()[3], # status = Background Well
           sample != "Background0" & is.na(max_bkgd) & !is.na(net_mfi)                ~ mfi_status()[2], # status = > Background Well
           sample != "Background0" &  net_mfi >= input$mfi_warning                    ~ mfi_status()[2], # status = > Background Well
           sample != "Background0" & net_mfi < input$mfi_warning & net_mfi > max_bkgd ~ mfi_status()[4], # status = Warning Well
           sample != "Background0" & net_mfi <= max_bkgd                              ~ mfi_status()[1], # status = < Background Well or No Signal
           sample != "Background0" & is.na(net_mfi)                                   ~ mfi_status()[1]  # status = < Background Well or No Signal
         )))  

    merged2 
  })
  
  plate_plot_data <- reactive({ #this reactive object will create a dataset that generates a 96-well plate template and then merge the user's uploaded xPONENT/BPM files based on which wells are present
    
    ############################ 96-well plate template  ####################################
    location <- c( "A1",  "A2",  "A3",  "A4",  "A5",  "A6",  "A7",  "A8",  "A9",  "A10", "A11", "A12", "B1",  "B2",  "B3",  "B4",  "B5",  "B6",  "B7", 
                   "B8",  "B9",  "B10", "B11", "B12", "C1",  "C2",  "C3",  "C4",  "C5",  "C6",  "C7",  "C8",  "C9",  "C10", "C11", "C12", "D1",  "D2", 
                   "D3",  "D4",  "D5",  "D6",  "D7",  "D8",  "D9",  "D10", "D11", "D12", "E1",  "E2",  "E3",  "E4",  "E5",  "E6",  "E7",  "E8",  "E9", 
                   "E10", "E11", "E12", "F1",  "F2",  "F3",  "F4",  "F5",  "F6",  "F7",  "F8",  "F9",  "F10", "F11", "F12", "G1",  "G2",  "G3",  "G4", 
                   "G5",  "G6",  "G7",  "G8",  "G9",  "G10", "G11", "G12", "H1",  "H2",  "H3",  "H4",  "H5",  "H6",  "H7",  "H8",  "H9",  "H10", "H11",
                   "H12")
    wellnum <- 1:96
    
    welltonum <- data.frame(location, wellnum, stringsAsFactors = FALSE)
    
    #                                 #                                 #
    
    cols <- rep(c("1", "2", "3", "4", "5", "6", "7", "8", "9","10", "11", "12"), times = 8)
    cols_n <- rep(1:12, times = 8)
    rows <- rep(c("A","B", "C", "D", "E", "F", "G", "H"), each = 12)
    rows_n <- rep(1:8, each = 12)
    
    plate_format <- data.frame(rows, rows_n, cols, cols_n) %>% 
      mutate(cols = fct_relevel(cols, c("1", "2", "3", "4", "5", "6", "7", "8", "9","10", "11", "12"))) 
    #############################################################################
    
    pdata1 <- all_plots_data() %>% 
      left_join(welltonum, "location") %>% # joining the user's data to the 96-well template
      select(location, wellnum, everything())

    pdata2 <- pdata1 %>% 
      mutate(rows = as.factor(str_sub(location, 1,1)), # cleaning the row strings in the user's data
             cols = fct_relevel(as.factor(str_sub(location, 2, 3)), # cleaning the column strings in the user's data
                                c("1", "2", "3", "4", "5", "6", "7", "8", "9","10", "11", "12"))) %>% 
      left_join(plate_format) %>% # joining the column and row formats to the user's cleaned up data
      select(location, rows, rows_n, cols, cols_n, wellnum, everything()) 
    
    pdata2
  })
  
  display_plate_plot <- reactive({ # BEAD COUNT QC: this reactive function displays a  single plate+antigen/antibody combination 96-well plate plot in the app
    
    req(input$ab,
        input$plates)
    
    pdata1 <- plate_plot_data()
    
    pdata2 <- pdata1 %>% # prepping the data to be plotted
      filter(antigen == input$ab,
             batch == input$plates) %>% 
      mutate(utlt2 = ifelse(utlt == paste("<", input$lthresh, "beads", sep = " "), # utlt2 will be used to label the legend
                            paste0(utlt, " ", "(", "n = ", sum(under_lt), ")"),
                            ifelse(
                              utlt == paste("<", input$uthresh, "beads", sep = " "),
                              paste0(utlt, " ", "(", "n = ", sum(under_ut), ")"),
                              paste0(utlt, " ", "(", "n = ", sum(over_ut), ")"))))
    
    
   plot1 <- ggplot(pdata2, aes(x = cols_n, y = rows_n))+ # generating the display plot
      geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color = "grey50",
                 fill = "white", shape = 21, size = 12)+ 
      geom_point(aes(color = utlt2), size = 15)+
      geom_text(aes(label = bead_count), size = 5)+
      scale_y_reverse(breaks = seq(1, 8),
                      limits = c(8.25, 0.75),
                      labels = LETTERS[1:8])+
      scale_x_continuous(breaks = seq(1, 12),
                         limits = c(0.75, 12.25),
                         position = "top")+
      labs(title = paste(input$plates, input$ab, sep = " "),
           x     = "",
           y     = "",
           color = "")+
      scale_color_manual(limits = c(paste0("<", " ", input$lthresh, " ", "beads", " ", "(", "n = ", sum(pdata2$under_lt), ")"), 
                                    paste0("<", " ", input$uthresh, " ", "beads", " ", "(", "n = ", sum(pdata2$under_ut), ")"),
                                    paste0("\u2265", " ", input$uthresh, " ", "beads", " ", "(", "n = ", sum(pdata2$over_ut), ")")),
                         values = c("red", "yellow", "springgreen3"))+
      theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
            legend.position = "right",
            legend.text = element_text(size = 16),
            panel.border = element_rect(color = "black", fill = NA, size = 2),
            strip.text.x = element_text(size = 16),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
   
   plot1
  })
  
  display_qc_plot <- reactive({ #this reactive function generates the single bead fluctuation QC plot in the app
    
    req(input$plates, 
        input$ab)
    
    if (input$utlt_platefail == "< Lower threshold"){
      bc_flag <- input$lthresh
    } else {
      bc_flag <- input$uthresh
    }
    
    data1 <- all_plots_data()
    
    data2 <- data1 %>% 
      mutate(batch2 = as.character(batch)) %>% #,
      group_by(batch) %>% 
      nest()
    
    pdata1 <- plate_plot_data() # reads in the merged dataset
    
    pdata2 <- pdata1 %>%
      mutate(status = ifelse(bead_count >= bc_flag, "Sufficient bead count", "Low bead count")) %>%
      filter(batch == input$plates, 
             antigen == input$ab) %>% 
      rename(`bead count` = bead_count)
    
    plot1 <- ggplot(pdata2, aes(label = location))+ # this block of code makes the bead QC plot that is outputed when plot_output() is called
      geom_hline(yintercept = input$lthresh, linetype = "dashed", color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = input$uthresh, linetype = "dashed", color = "forestgreen", size = 1, alpha = 0.5)+
      geom_hline(yintercept = mean(pdata2$`bead count`), color = "blue", size = 1, alpha = 0.3)+
      geom_line(aes(x = observation, y = `bead count`), alpha = 0.3)+
      geom_point(aes(x = observation, y = `bead count`, color = fct_rev(status)), size = 2)+
      scale_y_continuous(limits = c(0, NA))+
      scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60, 72, 84, 96))+
      labs(x = "Well order",
           y = "Bead count",
           color = "",
           title = "")+
      scale_color_manual(limits = c("Sufficient bead count", "Low bead count"),
                         values = c("royalblue4", "red"))+
      theme(legend.position = "bottom",
            legend.text = element_text(size = 16),
            panel.border = element_rect(color = "black", fill = NA))
    
    ggplotly(plot1, tooltip = c("label", "y")) # ggplot converted to plotly format
  })
  
  display_plate_plot_mfi <- reactive({ # MFI QC: this reactive function displays the single plate+antigen/antibody 96-well plate plot in the app
    
    req(input$ab,
        input$plates)
    
    pdata1 <- plate_plot_data()
    
    pdata2 <- pdata1 %>% 
      filter(antigen == input$ab,
             batch == input$plates) %>% 
      mutate(net_mfi = prettyNum(round(net_mfi, 0), big.mark = ",")) # formatting the net_mfi values to display with comma separators 
    
    ggplot(pdata2, aes(x = cols_n, y = rows_n))+ # generating the display plot
      geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color = "grey50",
                 fill = "white", shape = 21, size = 12)+
      geom_point(aes(color = mfi_well_status), size = 15)+
      geom_text(aes(label = net_mfi), size = 4)+
      scale_y_reverse(breaks = seq(1, 8),
                      limits = c(8.25, 0.75),
                      labels = LETTERS[1:8])+
      scale_x_continuous(breaks = seq(1, 12),
                         limits = c(0.75, 12.25),
                         position = "top")+
      labs(title = paste(input$plates, input$ab, sep = " "),
           x     = "",
           y     = "",
           color = "Well Status")+
      scale_color_manual(limits = mfi_status(),
                         values = c("red", "springgreen3", "mediumpurple", "yellow"))+
      theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            panel.border = element_rect(color = "black", fill = NA, size = 2),
            strip.text.x = element_text(size = 16),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  })
  
  display_mfi_heatmap <- reactive({ # MFI QC: this reactive function displays the single plate+antigen/antibody heat map plate plot in the app
    
    req(input$ab,
        input$plates)
    
    pdata1 <- plate_plot_data()
    
    pdata2 <- pdata1 %>% 
      filter(antigen == input$ab,
             batch == input$plates) %>% 
      mutate(median_fi2 = prettyNum(round(median_fi, 0), big.mark = ","), # formatting the net_mfi values to display with comma separators
             z_score = scale(net_mfi),
             log_mfi = log10(median_fi))  
   
   plot1 <- ggplot(pdata2, aes(x = cols_n, y = rows_n, fill = log_mfi))+ # generating the display plot
     geom_tile(color = "black")+
     geom_text(aes(label = median_fi2), size = 4)+
     scale_y_reverse(breaks = seq(1, 8),
                     labels = LETTERS[1:8])+
     scale_x_continuous(breaks = seq(1, 12),
                        position = "top")+
     labs(title = paste(input$plates, input$ab, sep = " "),
          x     = "",
          y     = "",
          fill = "log(MFI)")+
     scale_fill_gradientn(limits = c(0,6),#c(-4,4),
                          breaks = c(0,1,2,3,4,5,6),
                          labels = format(c(0,1,2,3,4,5,6)),
                          colors = viridis_pal()(7))+
     theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
           legend.position = "right",
           legend.title = element_text(size = 18),
           legend.text = element_text(size = 16),
           panel.background = element_rect(fill = "white"),
           panel.grid = element_blank(),#element_line(color = "gray75"),
           panel.border = element_blank(),#element_rect(color = "black", fill = NA, size = 1),
           strip.text.x = element_text(size = 16),
           axis.title = element_text(size = 16),
           axis.text = element_text(size = 14))+
     guides(fill = guide_legend(override.aes = list(size = 10)))
   
   plot1
    
  })
  
  all_plate_plots <- reactive({ # BEAD COUNT QC: this reactive function creates a list containing all of the 96-well plate layout plots for each plate faceted by antigen
    
    data1 <- plate_plot_data()
    
    data2 <- data1 %>% 
      mutate(batch2 = as.character(batch)) %>% 
      group_by(batch) %>% 
      nest() # nesting by batch so plot generation can be performed iteratively
    
    plateplot_maker <- function(df) { #function that will be applied iteratively to generate the plots
      
      fdata1 <- df %>%
        mutate(batch2 = as.factor(batch2))
      
      batch_name <- levels(fdata1$batch2)
      
      b1_plot <- ggplot(fdata1, aes(x = cols_n, y = rows_n))+
        facet_wrap(.~antigen, drop = TRUE)+ # faceting the batch plots by antigen and dropping any antigen levels that are not present in the batch
        geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color = "grey50",
                   fill = "white", shape = 21, size = 4)+ #5 4 3
        geom_point(aes(color = utlt), size = 7.5)+ #8 7 6
        geom_text(aes(label = bead_count), size = 2.5)+ #3
        scale_y_reverse(breaks = seq(1, 8),
                        limits = c(8.25, 0.75),
                        labels = LETTERS[1:8])+
        scale_x_continuous(breaks = seq(1, 12),
                           limits = c(0.75, 12.25),
                           position = "top")+
        labs(title = paste(batch_name, "Bead Count"),
             x     = "",
             y     = "",
             color = "")+
        scale_color_manual(limits = c(paste("<", input$lthresh, "beads", sep = " "), paste("<", input$uthresh, "beads", sep = " "), paste("\u2265", input$uthresh, "beads", sep = " ")),
                           values = c("red", "yellow", "springgreen3"))+
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = "bold"), 
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 20), 
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              strip.text.x = element_text(size = rel(1.5)), #28 #12
              axis.text = element_blank(),
              axis.ticks = element_blank())
      
      b1_plot
    }
    
    data3 <- data2 %>% 
      mutate(plots = map(data, plateplot_maker)) %>%  # runs plateplot_maker() iteratively over each batch and adds the output as a new list-column ("plots") to the nested dataframe
      ungroup() %>% 
      select(plots)
    
    data3
  })
  
  all_qc_plots <- reactive({# BEAD COUNT QC: this reactive function creates a list containing all of the bead QC plots for each plate faceted by antigen
    
    if (input$utlt_platefail == "< Lower threshold"){
      bc_flag <- input$lthresh
    } else {
      bc_flag <- input$uthresh
    }
    
    data1 <- all_plots_data()
    
    data2 <- data1 %>% 
      mutate(batch2 = as.character(batch),
             status = ifelse(bead_count >= bc_flag, "Sufficient bead count", "Low bead count")) %>% 
      group_by(batch) %>% 
      nest()
    
    qcplot_maker <- function(df) {
      
      fdata1 <- df %>% 
        mutate(batch2 = as.factor(batch2))
      
      batch_name <- levels(fdata1$batch2)
      
      b1_plot <- ggplot(df)+
        facet_wrap(.~antigen)+
        #coord_fixed(xlim = c(1,96))+
        geom_hline(aes(yintercept = lt), linetype = "dashed", color = "red", size = 0.5, alpha =0.5)+
        geom_hline(aes(yintercept = ut), linetype = "dashed", color = "forestgreen", size = 0.5, alpha = 0.5)+
        geom_line(aes(x = observation, y = bead_count), alpha = 0.25)+
        geom_hline(aes(yintercept = mean), color = "blue", size = 0.5, alpha = 0.3)+
        geom_point(aes(x = observation, y = bead_count, color = fct_rev(status)), size = 1.5)+
        scale_y_continuous(limits = c(0,NA))+
        labs(x = "Well order",
             y = "Bead count",
             color = "",
             title = paste(batch_name, "Bead QC Plots"))+
        scale_color_manual(limits = c("Sufficient bead count", "Low bead count"),
                           values = c("royalblue4", "red"))+
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
              legend.position = "bottom",
              legend.text = element_text(size = 18),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              strip.text.x = element_text(size = rel(1.5)),
              axis.title = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()
              )
      
      b1_plot
    }
    
    data3 <- data2 %>% 
      mutate(plots = map(data, qcplot_maker),
             filename = paste0(batch, ".png"),
             panels = ) %>% 
      ungroup() %>% 
      select(plots, filename)
    
    data3
  })
  
  all_plate_plots_mfi <- reactive({ # MFI QC: this reactive function creates a list containing all of the 96-well plate layout plots for each plate faceted by antigen
    
    # IDENTICAL TO all_plate_plots PROCEDURE EXCEPT MFI PLATE PLOTS ARE CREATED INSTEAD OF BEAD COUNT #
    
    data1 <- plate_plot_data()
    
    data2 <- data1 %>% 
      mutate(batch2 = as.character(batch),
             net_mfi = prettyNum(round(net_mfi, 0), big.mark = ",")) %>% 
      group_by(batch) %>% 
      nest()
    
    plateplot_maker_mfi <- function(df) {
      
      fdata1 <- df %>%
        mutate(batch2 = as.factor(batch2))
      
      batch_name <- levels(fdata1$batch2)
      
      plot1 <- ggplot(fdata1, aes(x = cols_n, y = rows_n))+
        facet_wrap(.~antigen, drop = TRUE)+
        geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color = "grey50",
                   fill = "white", shape = 21, size = 4) +
        geom_point(aes(color = mfi_well_status), size = 7.5)+
        geom_text(aes(label = net_mfi), size = 1.5)+
        scale_y_reverse(breaks = seq(1, 8),
                        limits = c(8.25, 0.75),
                        labels = LETTERS[1:8])+
        scale_x_continuous(breaks = seq(1, 12),
                           limits = c(0.75, 12.25),
                           position = "top")+
        labs(title = paste(batch_name, "MFI"),
             x     = "",
             y     = "",
             color = "")+
        scale_color_manual(limits = mfi_status(),
                           values = c("red", "springgreen3", "mediumpurple", "yellow"))+
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 18), 
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              strip.text.x = element_text(size = rel(1.5)), 
              axis.text = element_blank(),
              axis.ticks = element_blank())
      
      plot1
    }
    
    data3 <- data2 %>% 
      mutate(plots = map(data, plateplot_maker_mfi)) %>% 
      ungroup() %>% 
      select(plots)
    
    data3
  })
  
  all_mfi_heatmaps <- reactive({
    
    data1 <- plate_plot_data()
    
    data2 <- data1 %>% 
      mutate(batch2 = as.character(batch),
             median_fi2 = prettyNum(round(median_fi, 0), big.mark = ","), # formatting the net_mfi values to display with comma separators
             log_mfi = log10(median_fi)) %>% 
      group_by(batch) %>% 
      nest()
    
    heatmap_maker_mfi <- function(df) {
      
      fdata1 <- df %>%
        mutate(batch2 = as.factor(batch2))
      
      batch_name <- levels(fdata1$batch2)
      
      plot1 <- ggplot(fdata1, aes(x = cols_n, y = rows_n, fill = log_mfi))+ # generating the display plot
        facet_wrap(.~antigen, drop = TRUE)+
        geom_tile(color = "black")+
        geom_text(aes(label = median_fi2), size = 2)+
        scale_y_reverse(breaks = seq(1, 8),
                        labels = LETTERS[1:8])+
        scale_x_continuous(breaks = seq(1, 12),
                           position = "top")+
        labs(title = paste(batch_name, "MFI Heat Maps"),
             x     = "",
             y     = "",
             fill = "log(MFI)")+
        scale_fill_gradientn(limits = c(0,6),
                             breaks = c(0,1,2,3,4,5,6),
                             labels = format(c(0,1,2,3,4,5,6)),
                             colors = viridis_pal()(7))+
        theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 16),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              strip.text.x = element_text(size = rel(1.5)),
              strip.background = element_rect(size = 0.1),
              axis.text = element_blank(),
              axis.ticks = element_blank())+
        guides(fill = guide_legend(override.aes = list(size = 10)))
      
      plot1
    }
    
    data3 <- data2 %>% 
      mutate(plots = map(data, heatmap_maker_mfi)) %>% 
      ungroup() %>% 
      select(plots)
    
    data3
  })

  #### DATA MANIPULATION FOR CONTROL TRACKING PLOTS ####
  
  display_lj_plot <- reactive({ # creates an interactive LJ plot including Nelson rule violations
      
    data1 <- lj_data()
    
    vio_status <- c(lj_status()[1], lj_status()[2])
    status_colors <- c("forestgreen", "red")
    
     if (input$lj_choice == "Formatted Custom Data") {
      set_sample <- input$samples2_custom
      set_ag <- input$ab2_custom
     } else {
      set_sample <- input$samples3
      set_ag <- input$ab2
      }
      
    data2 <- data1 %>% 
      filter(sample == set_sample,
             antigen == set_ag) %>% 
      rename(`Median FI` = median_fi, Violations = vio)
      
    if (input$lj_choice == "Formatted Custom Data") {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = `Median FI`, label = Violations))
    } else {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = `Median FI`, label = Violations, text = batch))
    }
    
    lj_plot2 <- lj_plot1 +
      geom_line(aes(group = 1), alpha = 0.5)+ # must specify group = 1 or else the lines won't show up
      geom_hline(aes(yintercept = cl, linetype = "Mean"), color = "blue", size = 1, alpha = 0.4)+
      geom_hline(aes(yintercept = ucl3, linetype = "cl3"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl3, color = "red", size = 1, alpha = 0.5)+
      geom_hline(aes(yintercept = ucl2, linetype = "cl2"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl2, linetype = "longdash", color = "red", size = 1, alpha = 0.5)+
      geom_hline(aes(yintercept = ucl1, linetype = "cl1"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl1, linetype = "dotted", color = "red", size = 1, alpha = 0.5)+
      geom_point(x = data2$obs , y = data2$`Median FI`, aes(color = status), size = 3)+
      labs(x = "Observation",
           y = "MFI",
           color = "",
           title = paste(set_sample, " ", input$chart_type, " Individuals"," Chart ", "(", set_ag, ")", sep = ""))+
      scale_y_continuous(labels = comma)+
      scale_color_manual(limits = vio_status,
                         values = status_colors)+
      scale_linetype_manual(name = "", values = c(3,2,1,1),
                            guide = guide_legend(override.aes = list(color = c("red", "red", "red", "blue"))))+
      theme(plot.title = element_text(size = 16, vjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA))
    
    if (input$chart_type == "Shewhart") {
      
      lj_plot2 <- lj_plot2 +
        geom_vline(xintercept = input$shew_ref,
                   color = "purple4",
                   size = 1,
                   alpha = 0.75)
    }

  
    lj_plot3 <- ggplotly(lj_plot2, # converts non-interactive ggplot2 object into an interactive plotly object
                  tooltip = c("text", "x", "y", "label"))

      # fixing legend inconsistencies caused by ggplot2 ---> plotly conversion #
    
    for (i in 1:length(lj_plot3$x$data)) {
      
      if (!is.null(lj_plot3$x$data[[i]]$name)) {
        
        lj_plot3$x$data[[i]]$name <- str_replace(lj_plot3$x$data[[i]]$name, "\\(", "") %>% 
          str_replace(",1\\)", "")
      }
    }
    
    lj_plot3 
               
  })
  
  display_lj_plot_23sd <- reactive({ # creates an interactive LJ plot including CL2-CL3 violations
    
      data1 <- lj_data_23sd() 
      
      vio_status <- lj_status()
      status_colors <- c("forestgreen", "red",  "yellow1")
    
    if (input$lj_choice == "Formatted Custom Data") {
      set_sample <- input$samples2_custom
      set_ag <- input$ab2_custom
    } else {
      set_sample <- input$samples3
      set_ag <- input$ab2
    }
    
    data2 <- data1 %>% 
      filter(sample == set_sample,
             antigen == set_ag) %>% 
      rename(`Median FI` = median_fi, Violations = vio)
    
    if (input$lj_choice == "Formatted Custom Data") {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = `Median FI`))
    } else {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = `Median FI`, text = batch))
    }
    
    lj_plot2 <- lj_plot1 +
      geom_line(aes(group = 1), alpha = 0.5)+ # must specify group = 1 or else the lines won't show up
      geom_hline(aes(yintercept = cl, linetype = "Mean"), color = "blue", size = 1, alpha = 0.4)+
      geom_hline(aes(yintercept = ucl3, linetype = "cl3"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl3, color = "red", size = 1, alpha = 0.5)+
      geom_hline(aes(yintercept = ucl2, linetype = "cl2"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl2, linetype = "longdash", color = "red", size = 1, alpha = 0.5)+
      geom_hline(aes(yintercept = ucl1, linetype = "cl1"), color = "red", size = 1, alpha = 0.5)+
      geom_hline(yintercept = data2$lcl1, linetype = "dotted", color = "red", size = 1, alpha = 0.5)+
      geom_point(x = data2$obs , y = data2$`Median FI`, aes(color = status), size = 3)+
      labs(x = "Observation",
           y = "MFI",
           color = "",
           title = paste(set_sample, " ", input$chart_type, " Individuals", " Chart ", "(", set_ag, ")", sep = ""))+
      scale_y_continuous(labels = comma)+
      scale_color_manual(limits = vio_status,
                         values = status_colors)+
      scale_linetype_manual(name = "", values = c(3,2,1,1),
                            guide = guide_legend(override.aes = list(color = c("red", "red", "red", "blue"))))+
      theme(plot.title = element_text(size = 16, vjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA))
    
    if (input$chart_type == "Shewhart") {
      
      lj_plot2 <- lj_plot2 +
        geom_vline(xintercept = input$shew_ref,
                   color = "purple4",
                   size = 1,
                   alpha = 0.75)
    }
    
    lj_plot3 <- ggplotly(lj_plot2, # converts non-interactive ggplot2 object into an interactive plotly object
                         tooltip = c("text", "x", "y"))
    
    # fixing legend inconsistencies caused by ggplot2 ---> plotly conversion #
    
    for (i in 1:length(lj_plot3$x$data)) {
      
      if (!is.null(lj_plot3$x$data[[i]]$name)) {
        
        lj_plot3$x$data[[i]]$name <- str_replace(lj_plot3$x$data[[i]]$name, "\\(", "") %>% 
          str_replace(",1\\)", "")
      }
    }
    
    lj_plot3
    
  })  
  
  display_mr_plot <- reactive({ # creates the point-to-point moving range plot
    
    data1 <- mr_data()
    
    vio_status <- c(lj_status()[1], lj_status()[2])
    status_colors <- c("forestgreen", "red")
    
    if (input$lj_choice == "Formatted Custom Data") {
      set_sample <- input$samples2_custom
      set_ag <- input$ab2_custom
    } else {
      set_sample <- input$samples3
      set_ag <- input$ab2
    }
    
    data2 <- data1 %>% 
      filter(sample == set_sample,
             antigen == set_ag) %>% 
      rename(MR = mr)
    
    if (input$lj_choice == "Formatted Custom Data") {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = MR))
    } else {
      lj_plot1 <- ggplot(data2, aes(x = obs, y = MR, text = batch))
    }
    
    lj_plot2 <- lj_plot1 +
      geom_line(aes(group = 1), alpha = 0.5)+ # must specify group = 1 or else the lines won't show up
      geom_hline(aes(yintercept = mr_avg, linetype = "Mean MR"), color = "blue", size = 1, alpha = 0.4)+
      geom_hline(aes(yintercept = mr_ucl, linetype = "UCL"), color = "red", size = 1, alpha = 0.5)+
      geom_point(x = data2$obs , y = data2$MR, aes(color = status), size = 3)+
      labs(x = "Observation",
           y = "Moving range (MFI)",
           color = "",
           title = paste(set_sample, " Moving Range Chart ", "(", set_ag, ")", sep = ""))+
      scale_y_continuous(labels = comma)+
      scale_color_manual(limits = vio_status,
                         values = status_colors)+
      theme(plot.title = element_text(size = 16, vjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 16),
            legend.title = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA))
    
    lj_plot3 <- ggplotly(lj_plot2, # converts non-interactive ggplot2 object into an interactive plotly object
                         tooltip = c("text", "x", "y"))
    
    # fixing legend inconsistencies caused by ggplot2 ---> plotly conversion #
    
    for (i in 1:length(lj_plot3$x$data)) {
      
      if (!is.null(lj_plot3$x$data[[i]]$name)) {
        
        lj_plot3$x$data[[i]]$name <- str_replace(lj_plot3$x$data[[i]]$name, "\\(", "") %>% 
          str_replace(",1\\)", "")
      }
    }
    
    lj_plot3
  })
  
  all_lj_plots <- reactive({ # creates plots for each antigen/antibody faceted by control
    
    if (input$lj_type == "Nelson Rules") {
      
      data1 <- lj_data()
      vio_status <- c(lj_status()[1], lj_status()[2])
      status_colors <- c("forestgreen", "red")
    } else {
      
      data1 <- lj_data_23sd()
      vio_status <- lj_status()
      status_colors <- c("forestgreen", "red",  "yellow1")
    }
    
    data2 <- data1
    
    data3 <- data2 %>% 
      mutate(antigen2 = as.character(antigen)) %>% 
      group_by(antigen) %>% 
      nest()
    
    ljplot_maker <- function (df) { # function that will be applied iteratively to generate the plots
      
      fdata1 <- df %>% 
        mutate(antigen2 = as.factor(antigen2))
      
      ag_name <- levels(fdata1$antigen2)
      
      b1_plot <- ggplot(df, aes(x = obs , y = median_fi))+
        facet_wrap(.~sample,
                   scales = "free")+
        geom_line(alpha = 0.5)+
        geom_hline(aes(yintercept = cl, linetype = "Mean"), color = "blue", size = 1, alpha = 0.4)+
        geom_hline(aes(yintercept = ucl3, linetype = "cl3"), color = "red", size = 1, alpha = 0.5)+
        geom_hline(aes(yintercept = lcl3, linetype = "cl3"), color = "red", size = 1, alpha = 0.5)+
        geom_hline(aes(yintercept = ucl2, linetype = "cl2"), color = "red", size = 1, alpha = 0.5)+
        geom_hline(aes(yintercept = lcl2, linetype = "cl2"), color = "red", size = 1, alpha = 0.5)+
        geom_hline(aes(yintercept = ucl1, linetype = "cl1"), color = "red", size = 1, alpha = 0.5)+
        geom_hline(aes(yintercept = lcl1, linetype = "cl1"), color = "red", size = 1, alpha = 0.5)+
        geom_point(aes(x = obs , y = median_fi, color = status), size = 4)+
        labs(x = "Observation",
             y = "MFI",
             color = "",
             title = paste0(input$chart_type, " Individuals ", "(", ag_name, ")"))+
        scale_y_continuous(labels = comma)+
        scale_color_manual(limits = vio_status, 
                           values = status_colors)+
        scale_linetype_manual(name = "", values = c(3,2,1,1),
                              guide = guide_legend(override.aes = list(color = c("red", "red", "red", "blue"))))+
        theme(plot.title = element_text(size = 40, vjust = 0.5, face = "bold"),
              legend.position = "bottom",
              legend.text = element_text(size = 24),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 24),
              strip.text.x = element_text(size = 26),
              panel.border = element_rect(color = "black", fill = NA))
      
      if (input$chart_type == "Shewhart") {
        
        b1_plot <- b1_plot +
          geom_vline(xintercept = input$shew_ref,
                     color = "purple4",
                     size = 1,
                     alpha = 0.75)
        }
      
      b1_plot
    }
    
    data4 <- data3 %>% 
      mutate(plots = map(data, ljplot_maker)) %>% # runs ljplot_maker() iteratively over each antigen/antibody and adds the output as a list-column ("plots") to the nested dataframe
      ungroup() %>% 
      select(plots)
    
    data4
  })
  
  all_mr_plots <- reactive({
    
    data1 <- mr_data()
    
    vio_status <- c(lj_status()[1], lj_status()[2])
    status_colors <- c("forestgreen", "red")
    
    data2 <- data1
    
    data3 <- data2 %>% 
      mutate(antigen2 = as.character(antigen)) %>% 
      group_by(antigen) %>% 
      nest()
    
    mrplot_maker <- function (df) { # function that will be applied iteratively to generate the plots
      
      fdata1 <- df %>% 
        mutate(antigen2 = as.factor(antigen2))
      
      ag_name <- levels(fdata1$antigen2)
      
      b1_plot <- ggplot(df, aes(x = obs , y = mr))+
        facet_wrap(.~sample,
                   scales = "free")+
        geom_line(alpha = 0.5)+
        geom_hline(aes(yintercept = mr_avg, linetype = "Mean MR"), color = "blue", size = 1, alpha = 0.4)+
        geom_hline(aes(yintercept = mr_ucl, linetype = "UCL"), color = "red", size = 1, alpha = 0.5)+
        geom_point(aes(x = obs , y = mr, color = status), size = 4)+
        labs(x = "Observation",
             y = "Moving Range (MFI)",
             color = "",
             title = paste0("Moving Range Chart ", "(", ag_name, ")"))+
        scale_y_continuous(labels = comma)+
        scale_color_manual(limits = vio_status, 
                           values = status_colors)+
        theme(plot.title = element_text(size = 40, vjust = 0.5, face = "bold"),
              legend.position = "bottom",
              legend.text = element_text(size = 24),
              legend.title = element_blank(),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 24),
              strip.text.x = element_text(size = 26),
              panel.border = element_rect(color = "black", fill = NA))
      
      b1_plot
    }
    
    data4 <- data3 %>% 
      mutate(plots = map(data, mrplot_maker)) %>% # runs ljplot_maker() iteratively over each antigen/antibody and adds the output as a list-column ("plots") to the nested dataframe
      ungroup() %>% 
      select(plots)
    
    data4
  })
  
  #### DISPLAY OUTPUTS ####
  
  ## MODULE: QC Plots ##
  output$plateplot <- renderPlot({display_plate_plot()}) # this tells shiny what to use for the displaying the bead count 96-well plate plot
  
  output$plateplot_mfi <- renderPlot({# this tells shiny what to use for the displaying the MFI 96-well plate plot
    
    if (input$mfi_plot_switch == "Background/Warning"){
      display_plate_plot_mfi()
    } else {
      display_mfi_heatmap()
    }
    
    }) 
  
  output$qcplot <- renderPlotly({display_qc_plot()}) # this tells shiny what to use for the displaying the bead count fluctuation plot
  
  output$bc_repeats <- renderDataTable({ # this tells shiny what to use for displaying the bead count well failure summary table
    
    table1 <- bc_summary() %>% 
      select(batch, antigen, per_fail, status) %>% 
      filter(status == "Repeat") %>% 
      mutate(per_fail = percent(per_fail, accuracy = 1)) %>% 
      rename(`Failed Wells %` = per_fail,
             Batch = batch,
             Antigen = antigen,
             Status = status)
    
    table1
  }, rownames = FALSE)
  
  output$bc_repeats2 <- renderDataTable({bc_summary2()}, rownames = FALSE) # this tells shiny what to use for displaying the bead count plate failures summary table
  
  ## MODULE: Control Tracking ##
  output$lj <- renderPlotly({ # this tells shiny what to use for the interactive LJ plot
    
    if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
      validate(
        need(!is.null(input$samples2) & !is.null(input$ag_spec), "Please specify controls and antigens/antibodies"),
        errorClass = "ct")
    }
    
    if (input$lj_type == "Nelson Rules"){
      display_lj_plot()
    } else {
      display_lj_plot_23sd()
    }
  })
  
  output$mr <- renderPlotly({
    if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
      validate(
        need(!is.null(input$samples2) & !is.null(input$ag_spec), ""))
    }
    
    display_mr_plot()
    })
   
  output$nelson <- renderTable({ # this tells shiny what to use for displaying the Nelson rules explanation table
    
    data.table(`Nelson Rules` = paste("Rule", seq_len(8)),
               Violation = c("Point is outside CL3.",
                             "9 points in a row are on one side of the central line.",
                             "6 points in a row steadily increasing or decreasing.",
                             "14 or more points in a row alternate in direction, increasing then decreasing.",
                             "2/3 consecutive points beyond CL2 on the same side of the center line.",
                             "4/5 or 5/5 points in a row are beyond 1SD on the same side of the mean.",
                             "15 points in a row are all within 1SD on either side of the mean",
                             "8 points in a row outside 1SD of the mean in both directions."))
  })
  
  output$lj_summary <- renderDataTable({ # this tells shiny what to use for displaying the out of control points summary table
    if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
      validate(
        need(!is.null(input$samples2) & !is.null(input$ag_spec), "Please specify controls and antigens/antibodies"),
        errorClass = "ct")
    }
    
    if (input$lj_type == "Nelson Rules") {
      
      nelson_summary_small()
    } else {
      cl23_summary_small()
    }
      
  }, rownames = FALSE)
  
  output$lj_big_summary <- renderDataTable({ # this tells shiny what to use for displaying the summary statistics table
    
    if (input$lj_choice == "xPONENT/Bio-Plex Manager") {
      validate(
        need(!is.null(input$samples2) & !is.null(input$ag_spec), "Please specify controls and antigens/antibodies"),
        errorClass = "ct")
    }
    
    if (input$chart_type == "Levey-Jennings") { #Levey-Jennings analysis column names
      nelson_cols <- c("Control", "Antigen", "n", "n Missing", "Mean MFI", "Mean MR", "MR UCL", "Median MFI", "Max MFI", "Min MFI", "MFI SD",
                       "UCL3", "LCL3", "UCL2", "LCL2", "UCL1", "LCL1", "% OOC ", "Total Violations",
                       "n R1", "n R2", "n R3", "n R4", "n R5", "n R6", "n R7", "n R8")
    } else { #Shewart analysis column names
      nelson_cols <- c("Control", "Antigen", "n", "n Missing", "Mean MFI", "Mean MR", "MR UCL", "Median MFI", "Max MFI", "Min MFI", "MFI SD",
                       "UCL3", "LCL3", "UCL2", "LCL2", "UCL1", "LCL1", "% OOC ", "Total Violations",
                       "n R1", "n R2", "n R3", "n R4")
    }
    
    if (input$lj_type == "Nelson Rules") { # uses this if Nelson rule violations are selected
      datatable(nelson_summary(), 
                colnames = nelson_cols,
                rownames = FALSE,
                options = list( scrollCollapse = TRUE))
    } else { # uses this if CL2-CL3 violations are selected
      datatable(cl23_summary(), 
                colnames = c("Control", "Antigen", "n", "n Missing", "Mean MFI", "Mean MR", "MR UCL", "Median MFI", "Max MFI", "Min MFI", "MFI SD",
                             "UCL3", "LCL3", "UCL2", "LCL2", "n Outside CL3", "% Outside CL3", "n CL2-CL3", "% CL2-CL3"),
                rownames = FALSE,
                options = list(scrollCollapse = TRUE))
    }
  })
  
  output$lj_format <-renderTable({ # this tells shiny what to display for the custom formatted data example
    
    t1 <- data.table(observation = 1:10, antigen = "antigen 1", sample = "sample 1", median_fi = round(rnorm(10, 10000, 1000), 0))
    
    t2 <- data.table(observation = 1:10, antigen = "antigen 1", sample = "sample 2", median_fi = round(rnorm(10, 5000, 500), 0))
    
    t3 <- data.table(observation = 1:10, antigen = "antigen 2", sample = "sample 1", median_fi = round(rnorm(10, 15000, 1000), 0))
    
    t4 <- data.table(observation = 1:10, antigen = "antigen 2", sample = "sample 2", median_fi = round(rnorm(10, 1000, 150), 0))
    
    tf <- bind_rows(t1, t2, t3, t4)
    
    tf
    })
  
  #### DOWNLOADS ####
    ## Merged data and repeat summaries ##
  observeEvent(input$action3, { # this tells shiny what to do when the "action3" button is clicked
    
    shinyalert(title = "Now Generating Datasets",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/uprwwjptZW4Za/giphy.gif", #spongebob and patrick assembling the data
               #imageWidth = 335,
               #imageHeight = 250
               )
    
    
    withProgress(message = "Downloading merged data set...", value = 0, { #this adds in a download progress bar
      
      disable("uthresh") #disabling widgets to stop the user from confusing the app while it's assembling the data
      disable("lthresh")
      disable("ab")
      disable("plates")
      disable("saver")
      
      incProgress(amount = 1/5, detail = "Assembling raw data...") #updates download progress bar when it starts working on merging the uploaded xPONENT files
      
      # raw data (long format) compiling #
      data1 <- as.data.frame(raw_data())
      
      fwrite(data1,
             file = paste0(tempdir(), "/", "raw_data.csv"),
             row.names = FALSE)#,
      
      zipr(zipfile = paste0(tempdir(), "/", "raw_sum_data.zip"), # creates a zip folder to store the output
           files = paste0(tempdir(), "/", "raw_data.csv"))
      
      
      # raw data (wide format) compiling #
      
      data2 <- as.data.frame(raw_data_wide()) # wide format merged data
      
      fwrite(data2, 
             file = paste0(tempdir(), "/", "raw_data_wide.csv"),
             row.names = FALSE)
      
      zipr_append(zipfile = paste0(tempdir(), "/", "raw_sum_data.zip"), # adds the wide format data to the zip folder
                  files = paste0(tempdir(), "/", "raw_data_wide.csv"))
      
      # summary data compiling #
      
      output_workbook <- createWorkbook() #creates a blank excel workbook
      
      incProgress(amount = 1/5, detail = paste("Assembling % Wells < Upper Threshold ","(", input$uthresh, ")", sep = "")) #updates download progress bar when it starts working on the % wells < the upper threshold summary data
      sheet.1 <- addWorksheet(output_workbook, sheetName = paste("% Wells < Upper Threshold ","(", input$uthresh, ")", sep = "")) # creates a blank sheet with a name in the workbook
      writeData(output_workbook, as.data.frame(psum_ut()), sheet = sheet.1) # adds the % wells < upper threshold output from psum_ut() to the sheet
      
      incProgress(amount = 1/5, detail = paste("Assembling % Wells < Lower Threshold ","(", input$lthresh, ")", sep = ""))  #updates download progress bar when it starts working on the % wells < the lower threshold summary data
      sheet.2 <- addWorksheet(output_workbook, sheetName = paste("% Wells < Lower Threshold ","(", input$lthresh, ")", sep = "")) # creates a 2nd blank sheet with a name in the workbook
      writeData(output_workbook, as.data.frame(psum_lt()), sheet = sheet.2) # adds the % wells < lower threshold output from psum_lt() to the sheet
      
      # repeats #
      incProgress(amount = 1/5, detail = "Assembling plate repeats summary data...")
      
      repeats1 <- bc_summary() %>% 
        mutate(per_fail = percent(per_fail, 1)) %>% 
        rename(`Failed Wells %` = per_fail) %>% 
        select(-c(status_flag))
      
      sheet.3 <- addWorksheet(output_workbook, sheetName = paste0("Repeats plate-ag", " (", input$bc_platefail, "%=Fail)"))  # creates a 3rd blank sheet with a name in the workbook
      writeData(output_workbook, as.data.frame(repeats1), sheet = sheet.3)
      
      sheet.4 <- addWorksheet(output_workbook, sheetName = "Repeats all plates")
      writeData(output_workbook, as.data.frame(bc_summary2()), sheet = sheet.4)  # creates a 4th blank sheet with a name in the workbook
      
      
      saveWorkbook(output_workbook, 
                   file = paste0(tempdir(), "/", "summary_results.xlsx"),
                   overwrite = TRUE) # returns the saved workbook
      
      incProgress(amount = 1/5, detail = "Finalizing...") # updates the download progress bar when finished
      
      zipr_append(zipfile = paste0(tempdir(), "/", "raw_sum_data.zip"), # adds the summary results .xlsx file to the zip folder
                  files = paste0(tempdir(), "/", "summary_results.xlsx"))
      
      
      
      shinyalert(title = "Done!",
                 text = "Datasets are now ready to be downloaded.",
                 type = "success",
                 confirmButtonCol = "#5dade2"#,
                 #imageUrl = "https://media.giphy.com/media/26u4lOMA8JKSnL9Uk/giphy.gif", #spongebob wiping his hands
                 #imageWidth = 335,
                 #imageHeight = 250
                 )
      
      enable("uthresh") # user widgets enabled again once the download is complete
      enable("lthresh")
      enable("ab")
      enable("plates")
      enable("saver")
    })
  })
  
  output$saver <- downloadHandler( # this tells the app what to do when the "saver" download button is clicked
    filename = "Raw_sum_data.zip",# default file name
    content = function(file){ # this tells the app how to save the raw data output
      
      file.copy(paste0(tempdir(), "/", "raw_sum_data.zip"), file) # returns the zip folder containing the merged and summary datasets to the user
    })
  
    ## Bead count fluctuation plots ##
  observeEvent(input$action, { # this tells shiny what to do when the "action" button is clicked
    
    shinyalert(title = "Now Generating Plots",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/xT1R9ZORSvL6jtqOeQ/giphy.gif",
               #imageWidth = 320,
               #imageHeight = 180
               )
    
    
    disable("bqcplots")
    
    withProgress(message = "Assembling Bead QC plots", value = 0, {
      
      data1 <- all_qc_plots()
      
      dir.create(temp <- tempfile()) # creates a new temporary directory where the plots will be stored
      
      saver <- function(df){ # function that will be run in parallel to save the high resolution plots
        
        ggsave(paste0(df$labels$title, ".png"), path = temp, plot = df, 
               width = ifelse(length(unique(df$data$antigen)) > 6, 26, 10.6), #specifying plot sizes based on the number of panels
               height = ifelse(length(unique(df$data$antigen)) > 6, 18, 6),
               units = "in",
               device = "png")
      }
      
      incProgress(amount = 1/2, detail = "Generating Plots...")
      
      future_map(data1$plots, saver, .options = furrr_options(packages = "forcats"))#, #future_map runs saver() in parallel to speed up the downloads if futures plan is set to multicore
      
      incProgress(amount = 1, detail = "Zipping files...")
      
      plot_list <- list.files(path = temp, # returns all of the QC plot file names
                              pattern = "*.png",
                              full.names = TRUE) 
      
      zipr(zipfile = paste0(tempdir(), "/", "Bead QC plots.zip"), # zips all of the QC plots together
           files = plot_list)
      
      shinyalert(title = "Done!",
                 text = "Plots are now ready to be downloaded!",
                 type = "success",
                 confirmButtonCol = "#5dade2"#,
                 #imageUrl = "https://media.giphy.com/media/xUNd9MHCLPVyBIpNbW/giphy.gif",
                 #imageWidth = 320,
                 #imageHeight = 180
                 )
      
      enable("bqcplots")
    })
  })
  
  output$bqcplots <- downloadHandler(
    filename = "Bead-QC-plots.zip",
    content = function(file){
      
      file.copy(paste(tempdir(), "/", "Bead QC plots.zip", sep = ""), file) # returns the zip folder containing the QC plots to the user
    }
  )
  
    ## Bead count 96-well plate plots ##
          # See the "Bead count fluctuation plots" section for code explanations #
  observeEvent(input$action2, {
    
    shinyalert(title = "Now Generating Plots",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/HUplkVCPY7jTW/giphy.gif", # magic computer doing all the things
               #imageWidth = 250,
               #imageHeight = 350
               )
    
    disable("plateplotsaver")
    
    withProgress(message = "Assembling Bead Plate plots", value = 0, {
      
      data1 <- all_plate_plots()
      
      dir.create(temp <- tempfile())
      
      saver <- function(df){
        
        ggsave(paste0(df$labels$title, " plate", ".png"), path = temp, plot = df, 
               width = ifelse(length(unique(df$data$antigen)) > 6, 26, 10.6), 
               height = ifelse(length(unique(df$data$antigen)) > 6, 18, 6),
               device = "png")
      }
      
      incProgress(amount = 1/2, detail = "Generating Plots...")
      
      future_map(data1$plots, saver, .options = furrr_options(packages = "forcats"))
      
      incProgress(amount = 1, detail = "Zipping files...")
      
      plot_list <- list.files(path = temp,
                              pattern = "*.png",
                              full.names = TRUE)
      
      zipr(zipfile = paste0(tempdir(), "/", "Bead QC plate plots.zip"),
           files = plot_list)
      
      shinyalert(title = "Done!",
                 text = "Plots are now ready to be downloaded!",
                 type = "success",
                 confirmButtonCol = "#5dade2",
                 #imageUrl = "https://media.giphy.com/media/3ohryiYkE0DVwdLAys/giphy.gif",
                 #imageWidth = 250,
                 #imageHeight = 250
                 )
      
      enable("plateplotsaver")
    })
  })
  
  
  output$plateplotsaver <- downloadHandler( # this tells the app what to do when the "plateplotsaver" download button is clicked
    filename = "Bead-QC-plate-plots.zip", #default file name
    content = function(file){
      
      file.copy(paste(tempdir(), "/", "Bead QC plate plots.zip", sep = ""), file)
    }
  )
  
    ## MFI 96-well plate plots ##
          # See the "Bead count fluctuation plots" section for code explanations #
  observeEvent(input$action5, {
    
    shinyalert(title = "Now Generating Plots",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/DIbzujHh2PCbm/giphy.gif", # magic computer doing all the things
               #imageWidth = 281,
               #imageHeight = 225
               )
    
    disable("plateplotsaver_mfi")
    
    withProgress(message = "Assembling MFI Plate plots", value = 0, {
      
      data1 <- all_plate_plots_mfi()
      
      dir.create(temp <- tempfile())
      
      saver <- function(df){
        
        ggsave(paste0(df$labels$title, " plate", ".png"), path = temp, plot = df, 
               width = ifelse(length(unique(df$data$antigen)) > 6, 26, 10.6), 
               height = ifelse(length(unique(df$data$antigen)) > 6, 18, 6),  
               device = "png")
      }# crazy hiarious gif: https://media.giphy.com/media/DIbzujHh2PCbm/giphy.gif
      # good pairing gif: https://media.giphy.com/media/mJHSkWKziszrkcNJPo/giphy.gif 
      
      incProgress(amount = 1/2, detail = "Generating Plots...")
      
      future_map(data1$plots, saver)
      
      incProgress(amount = 1, detail = "Zipping files...")
      
      plot_list <- list.files(path = temp,
                              pattern = "*.png",
                              full.names = TRUE)
      
      zipr(zipfile = paste0(tempdir(), "/", "MFI QC plate plots.zip"),
           files = plot_list)
      
      shinyalert(title = "Done!",
                 text = "Plots are now ready to be downloaded!",
                 type = "success",
                 confirmButtonCol = "#5dade2",
                 imageUrl = "https://media.giphy.com/media/mJHSkWKziszrkcNJPo/giphy.gif",
                 imageWidth = 350,
                 imageHeight = 200)
      
      enable("plateplotsaver_mfi")
    })
  })
  
  output$plateplotsaver_mfi <- downloadHandler( # this tells the app what to do when the "plateplotsaver" download button is clicked
    filename = "MFI-QC-plate-plots.zip", #default file name
    content = function(file){# this tells the app how to assemble the plots
      
      file.copy(paste(tempdir(), "/", "MFI QC plate plots.zip", sep = ""), file)
    }
  )
  
    ## MFI heat maps ##
  observeEvent(input$action7, {
    
    shinyalert(title = "Now Generating Plots",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/DIbzujHh2PCbm/giphy.gif", # magic computer doing all the things
               #imageWidth = 281,
               #imageHeight = 225
               )
    
    disable("heatmapsaver_mfi")
    
    withProgress(message = "Assembling MFI Heat Maps", value = 0, {
      
      data1 <- all_mfi_heatmaps()
      
      dir.create(temp <- tempfile())
      
      saver <- function(df){
        
        ggsave(paste0(df$labels$title, " plate", ".png"), path = temp, plot = df, 
               width = ifelse(length(unique(df$data$antigen)) > 6, 32, 12.5), 
               height = ifelse(length(unique(df$data$antigen)) > 6, 18, 6), 
               device = "png")
      }# crazy hiarious gif: https://media.giphy.com/media/DIbzujHh2PCbm/giphy.gif
      # good pairing gif: https://media.giphy.com/media/mJHSkWKziszrkcNJPo/giphy.gif 
      
      incProgress(amount = 1/2, detail = "Generating Plots...")
      
      future_map(data1$plots, saver)
      
      incProgress(amount = 1, detail = "Zipping files...")
      
      plot_list <- list.files(path = temp,
                              pattern = "*.png",
                              full.names = TRUE)
      
      zipr(zipfile = paste0(tempdir(), "/", "MFI Heat Maps.zip"),
           files = plot_list)
      
      shinyalert(title = "Done!",
                 text = "Plots are now ready to be downloaded!",
                 type = "success",
                 confirmButtonCol = "#5dade2"#,
                 #imageUrl = "https://media.giphy.com/media/mJHSkWKziszrkcNJPo/giphy.gif",
                 #imageWidth = 350,
                 #imageHeight = 200
                 )
      
      enable("heatmapsaver_mfi")
    })
  })
  
  output$heatmapsaver_mfi <- downloadHandler( # this tells the app what to do when the "plateplotsaver" download button is clicked
    filename = "MFI-Heat-Maps.zip",#.xlsx",#function(){"Bead QC plate plots.xlsx"}, #default file name
    content = function(file){# this tells the app how to assemble the plots
      
      file.copy(paste(tempdir(), "/", "MFI Heat Maps.zip", sep = ""), file)
    }
  )
    ## LJ Plots ##
          # See the "Bead count fluctuation plots" section for code explanations #
  observeEvent(input$action4, {
    
    shinyalert(title = "Now Generating Plots",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/3oKIPEqDGUULpEU0aQ/giphy.gif",
               #imageWidth = 300,
               #imageHeight = 300
               )
    
    withProgress(message = "Assembling control tracking plots", value = 0, {
      
      disable("ljsaver")
      
      data1 <- all_lj_plots()
      
      dir.create(temp <- tempfile())
      
      saver <- function(df){
        
        ggsave(paste0(df$labels$title, ".png"), path = temp, plot = df, 
               width = 24, 
               height = 18, 
               device = "png")
      }
      
      incProgress(amount = 1/2, detail = "Generating Plots...")
      
      future_map(data1$plots, saver, .options = furrr_options(packages = "forcats"))
      
      if (input$chart_type == "Moving Range"){
        
        data2 <- all_mr_plots()
        
        future_map(data2$plots, saver, .options = furrr_options(packages = "forcats"))
      } 
      
      incProgress(amount = 1, detail = "Zipping files...")
      
      plot_list <- list.files(path = temp,
                              pattern = "*.png",
                              full.names = TRUE)
      
      zipr(zipfile = paste0(tempdir(), "/", "Control tracking plots.zip"),
           files = plot_list)
      
      shinyalert(title = "Done!",
                 text = "Control tracking plots are now ready to be downloaded.",
                 type = "success",
                 confirmButtonCol = "#5dade2"#,
                 #imageUrl = "https://media.giphy.com/media/3og0IUzdgwVczU67eg/giphy.gif",
                 #imageWidth = 300,
                 #imageHeight = 225
                 )
      
      enable("ljsaver")
    })
  })
  
  output$ljsaver <- downloadHandler(
    filename = "control_tracking_plots.zip",
    content = function(file){
      file.copy(paste0(tempdir(), "/", "Control tracking plots.zip"), file)
    }
  )
  
    ## LJ summary tables ##
        # See "Merged data and repeat summaries" section for code explanation #
  observeEvent(input$action6, {
    
    shinyalert(title = "Now Generating Control Tracking Summary Tables",
               text  = "Please wait until the app is finished before using other features.",
               type = "",
               confirmButtonCol = "#5dade2"#,
               #imageUrl = "https://media.giphy.com/media/uprwwjptZW4Za/giphy.gif", #spongebob and patrick assembling the data
               #imageWidth = 335,
               #imageHeight = 250
               )
    
    
    withProgress(message = "Assembling control tracking summary results...", value = 0, { #this adds in a download progress bar
      
      disable("uthresh") #disabling widgets to stop the user from confusing the app while it's assembling the data
      disable("lthresh")
      disable("ab")
      disable("plates")
      disable("saver")
      disable("lj_summary_saver")
      
      incProgress(amount = 1/3, detail = "Assembling Nelson rules summary table...")
      
      fwrite(as.data.table(nelson_summary()), 
             file = paste0(tempdir(), "/", "nelson_summary.csv"),
             row.names = FALSE)#,
      
      
      zipr(zipfile = paste0(tempdir(), "/", "control_tracking_summary.zip"),
           files = paste0(tempdir(), "/", "nelson_summary.csv"))
      
      incProgress(amount = 1/3, detail = "Assembling 2CL-3CL summary table...")
      ###
      
      fwrite(as.data.table(nelson_summary_small()), 
             file = paste0(tempdir(), "/", "nelson_repeat_summary.csv"),
             row.names = FALSE)#,
      
      
      zipr_append(zipfile = paste0(tempdir(), "/", "control_tracking_summary.zip"),
           files = paste0(tempdir(), "/", "nelson_repeat_summary.csv"))
      
      ###
    
      fwrite(as.data.table(cl23_summary()), 
             file = paste0(tempdir(), "/", "2cl_3cl_summary.csv"),
             row.names = FALSE)
      
      
      zipr_append(zipfile = paste0(tempdir(), "/", "control_tracking_summary.zip"),
                  files = paste0(tempdir(), "/", "2cl_3cl_summary.csv"))
      
      ###
      
      fwrite(as.data.table(cl23_summary_small()), 
             file = paste0(tempdir(), "/", "2cl_3cl_repeat_summary.csv"),
             row.names = FALSE)#,
      
      
      zipr_append(zipfile = paste0(tempdir(), "/", "control_tracking_summary.zip"),
           files = paste0(tempdir(), "/", "2cl_3cl_repeat_summary.csv"))
      
      ###
      
      if (input$chart_type == "Moving Range"){
        
        fwrite(as.data.table(mr_summary_small()), 
               file = paste0(tempdir(), "/", "mr_repeat_summary.csv"),
               row.names = FALSE)#,
        
        
        zipr_append(zipfile = paste0(tempdir(), "/", "control_tracking_summary.zip"),
                    files = paste0(tempdir(), "/", "mr_repeat_summary.csv"))
      }
      
      ###
      
      shinyalert(title = "Done!",
                 text = "Summary tables are now ready to be downloaded.",
                 type = "success",
                 confirmButtonCol = "#5dade2"#,
                 # imageUrl = "https://media.giphy.com/media/26u4lOMA8JKSnL9Uk/giphy.gif", 
                 # imageWidth = 335,
                 # imageHeight = 250
                 )
      
      enable("uthresh") # user widgets enabled again once the download is complete
      enable("lthresh")
      enable("ab")
      enable("plates")
      enable("saver")
      enable("lj_summary_saver")
      
    })
  })
  
  output$lj_summary_saver <- downloadHandler( # this tells the app what to do when the "lj_summary_saver" download button is clicked
    filename = "control tracking summary results.zip",# default file name
    content = function(file){ # this tells the app how to save the raw data output
      file.copy(paste0(tempdir(), "/", "control_tracking_summary.zip"), file)
    })
  
  #### End  ####
}

