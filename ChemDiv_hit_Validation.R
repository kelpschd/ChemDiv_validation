list.of.packages <- c("xml2", "openxlsx", "tidyverse", "ggpubr", "ggbeeswarm", "plater", "gridExtra", "cowplot", "scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

Library_decoder <- as_tibble(read.csv("CHEMdiv_rescreen_DJK.csv"))

##round_any() is a function in the package plyr, but plyr has some 
##compatibility issues with the tidyverse, thus just writing in the 
##function manually
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

Define_outliers <- function(x) {
  Validation.outliers <<- x %>% 
    filter(Drug_name == i) %>% 
    group_by(Treatment) %>% 
    mutate(RLU_Q1 = (quantile(RLU)[[2]]), 
           RLU_Q3 = (quantile(RLU)[[4]]), 
           avg = mean(RLU),
           RLU_IQR = IQR(RLU),
           Outlier = if_else(RLU > RLU_Q3+1.5*(RLU_IQR) | RLU < RLU_Q1-1.5*(RLU_IQR), "Outlier", "Normal"))}

Define_y_limit <- function(x) {
  y_lim <<- x %>%
    filter(Outlier != "Outlier") %>% 
    pull(RLU) %>% 
    max(.) %>%
    round_any(.,50000, ceiling)}

plot_validation_data <- function(x) {
  m <<- ggplot(data = x, 
               mapping = aes(x = factor(Treatment, levels = c(
                 "Vehicle",
                 "5 µM lomitapide",
                 "0.0625 µM",
                 "0.125 µM",
                 "0.25 µM",
                 "0.5 µM",
                 "1 µM",
                 "2 µM",
                 "4 µM",
                 "8 µM")), 
                 y = RLU)) +
    geom_boxplot() +
    geom_point() +
    scale_y_continuous(
      name = "Relative Luminescence Units (RLU)",
      labels = comma, limits = c(0,y_lim), #c(lower limit, upper limit)
      breaks = seq(0,y_lim,y_lim/10)) +  #seq(lower limit, upper limit, tick frequency)
    theme(
      axis.text.x = element_text(face = 'bold', angle = 45, hjust = 1), 
      axis.title.x = element_blank(),
      axis.text.y = element_text(face = 'bold'), 
      axis.title.y = element_text(face = 'bold'),
      legend.title = element_text(face = 'bold'),
      legend.text = element_text(face = 'bold'),
      strip.text.x = element_text(face = 'bold')) +
    labs(title = paste("Validation of ",i, sep=""), 
         subtitle = paste(Outlier.count, " total outliers removed", sep = ""),
         caption = paste("Experiment: ", a, sep = ""))
}

##following are the functions to calculate SSMD
SSMD_hs <- function(d) {
  N = length(d[!is.na(d)])
  (gamma((N-1)/2)/gamma((N-2)/2))*(sqrt(2/(N-1)))*(mean(d, na.rm = TRUE)/sd(d, na.rm = TRUE))
}

SSMD_qc <- function(a, b) {
  #a = sample
  #b = negative control
  (median(a, na.rm = TRUE)-median(b, na.rm = TRUE))/(1.4826*sqrt(((mad(a, na.rm = TRUE))^2)+((mad(b, na.rm = TRUE))^2)))
}

#Fix up these tibbles, they have a lot of uncessary info
All_plates <- tibble(Row = character(),
                     Column = numeric(),
                     RLU = numeric(),
                     Dose = character(),
                     Drug = character(),
                     Date_of_read = as.Date(NA),
                     Filename = character(),
                     Library = character(),
                     ChemDiv_rescreen_plate = numeric(),
                     JHDL_test_plate = numeric(),
                     Original_plate_ID = numeric(),
                     JHDL_stock_location = numeric(),
                     JHDL_well_location = character(),
                     Compound_name = character(), 
                     Fold_Change = numeric(),
                     log2_RLU = numeric(), 
                     Difference = numeric(),
                     log2_Fold_Change = numeric())

Summary_data <- tibble(Column = numeric(),
                       Sample_size = integer(),
                       hsSSMD = numeric(),
                       qcSSMD = numeric(),
                       Average_RLU = numeric(),
                       log2_Average_Fold_Change = numeric(),
                       Dose = character(),
                       Drug = character(),
                       Date_of_read = as.Date(NA),
                       Filename = character(), 
                       Library = character(),
                       ChemDiv_rescreen_plate = numeric(),
                       JHDL_test_plate = numeric(),
                       Original_plate_ID = numeric(),
                       JHDL_stock_location = character(),
                       JHDL_well_location = character(),
                       Compound_name = character())

target <- c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12",
            "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12",
            "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12",
            "D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12",
            "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
            "F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
            "G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12",
            "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")








files <- list.files("./input_data", full.names = TRUE) #generates a list of all the files within the working directory

#WORK ON MAKING THIS INTO A SINGLE OR LARGER FUNCTION? MAKE IT JUST A CALL
#STATUS BAR WOULD BE FUN TOO
for (document in files){
  name <- document #assigns a filename to the object 'name'
  XMLfile <- read_xml(name) #reads the file
  
  #Well label extraction
  wells <- xml_find_all(XMLfile, "///Well") 
  well_labels <- xml_attr(wells, "Pos")
  
  #RLU extraction
  values <- xml_find_all(XMLfile, "////Single")
  RLU <- xml_text(values) #identifies the luminescence reads and assigns them to the object 'RLU'
  RLU <- as.numeric(RLU)
  
  #Date extraction
  Date_node <- xml_child(XMLfile, search = 5)
  Read_start <- xml_attr(Date_node, "Time_Start") #assigns the start time of reading to the object Read_start
  Date_of_read <- as.Date(strptime(Read_start, "%Y-%m-%d"))
  #Read_start <- strptime(Read_start, "%Y-%m-%dT%H:%M:%S")
  
  #Generate and fill tibble
  Plate <- tibble(well_labels, RLU) #generates a dataframe with well labels and RLU paired
  Plate <- Plate %>% arrange(match(target, well_labels)) #reorganizes the dataframe according to target
  Plate <- Plate %>% separate(well_labels, into = c("Row", "Column"), sep = "(?<=[[:upper:]])")
  Plate <- Plate %>% mutate_at(vars(Column), as.numeric)
  Dose <- c(NA, "5 µM", "8 µM", "4 µM", "2 µM", "1 µM", "0.5 µM", "0.25 µM", "0.125 µM", "0.0625 µM", "Vehicle", NA)
  Drug <- c(NA, "Lomitapide", "A", "A", "A", "A", "A", "A", "A", "A", "Vehicle", NA)
  
  ##Remove drug/dose NAs 
  
  name <- name %>% str_remove(".xml") %>% str_remove('./input_data/') #removes filename directory information and '.xml'
  Plate <- Plate %>% add_column(Dose = rep(Dose, 8),
                                Drug = rep(Drug, 8),
                                Date_of_read = rep(Date_of_read),
                                Filename = rep(name)) %>% separate(col = Filename, 
                                                                   into = c("ChemDiv_rescreen_plate"), 
                                                                   sep = "_", 
                                                                   remove = FALSE, 
                                                                   convert = TRUE, 
                                                                   extra = "drop") %>% 
    left_join(Library_decoder, by = c("ChemDiv_rescreen_plate"), keep = TRUE) #this works but its being a pain in the ass right now, need to incorporate the columns into All_data and Summary_data
  
  ##Outlier handling here, could help wiht Vehicle values issues!
    
  #Calculate the mean and median of the 8 vehicle samples
  Vehicle_Values <- Plate %>% group_by(Column) %>% filter(Column == 11) %>% pull(RLU)
  
  #Plate analysis
  Plate <- Plate %>% mutate(Fold_Change = RLU/mean(Vehicle_Values, na.rm = TRUE), 
                            log2_RLU = log2(RLU), 
                            Difference = log2_RLU-log2(median(Vehicle_Values, na.rm = TRUE)), 
                            log2_Fold_Change = log2(Fold_Change))
  
  #Summary statistics
  Plate_summary <- Plate %>% group_by(Column) %>% summarise(Sample_size = length(log2_RLU[!is.na(log2_RLU)]),
                                                            hsSSMD = SSMD_hs(Difference),
                                                            qcSSMD = SSMD_qc(RLU, Vehicle_Values),
                                                            Average_RLU = mean(RLU, na.rm = TRUE), 
                                                            log2_Average_Fold_Change = log2(mean(Fold_Change, na.rm = TRUE)))
 
  Plate_summary <- Plate_summary %>% add_column(Dose = Dose,
                                                Drug = Drug,
                                                Date_of_read = rep(Date_of_read),
                                                Filename = rep(name)) %>% separate(col = Filename, 
                                                                                   into = c("ChemDiv_rescreen_plate"), 
                                                                                   sep = "_", 
                                                                                   remove = FALSE,
                                                                                   convert = TRUE,
                                                                                   extra = "drop")
  
  #Individual plate graphing
  #RLU boxplot
  #p1 <- ggplot(Plate %>% filter(!is.na(Drug)), aes(x = factor(Dose, levels = c("Vehicle","5 µM","8 µM","4 µM","2 µM","1 µM")),y = RLU)) + 
  #  scale_x_discrete(labels = c('label 1' = expression("Vehicle"))) +
  #  geom_boxplot() +
  #  labs(title = paste(name),
  #       subtitle = "Boxplot of RLU",
  #       x = "Drug treatment",
  #       y = "Relative Luminescence Units (RLU)") +
  #  facet_grid(cols = vars(factor(Drug, levels = c("Vehicle", "Lomitapide", "A", "B"))), scales = "free", space = "free")
  #ggsave(filename = paste(name, "_graph.png", sep = ""), width = 7, height = 5)
  #Fold change boxplot
  #p2 <- ggplot(Plate %>% filter(!is.na(Drug)), aes(x = factor(Dose, levels = c("Vehicle","5 µM","8 µM","4 µM","2 µM","1 µM")),y = Fold_Change)) + 
  #  scale_x_discrete(labels = c('label 1' = expression("Vehicle"))) +
  #  geom_boxplot() +
  #  labs(title = paste(name),
  #       subtitle = "Boxplot of fold change",
  #       x = "Drug treatment",
  #       y = "Relative fold change") +
  #  facet_grid(cols = vars(factor(Drug, levels = c("Vehicle", "Lomitapide", "A", "B"))), scales = "free", space = "free")
  #ggsave(filename = paste(name, "_graphs.png", sep = ""), width = 15, height = 7.5, arrangeGrob(p1, p2, ncol = 2, widths = c(1,1)))
  
  #Add to C...
  All_plates <- All_plates %>% add_row(Plate)
  Summary_data <- Summary_data %>% add_row(Plate_summary)
  
  #matrix generation - this is mostly wrong at the momment, but still helpful to do, worth fixing up
  RLUmatrix <- matrix(Plate$RLU, ncol = 12, byrow = T)
  rownames(RLUmatrix) <- c("1","2","3","4","5","6","7","8")
  colnames(RLUmatrix) <- c(paste(name,"empty"),
                           paste(name,"5 µM lomitapide"),
                           paste(name,"8 µM Drug A"),
                           paste(name,"4 µM Drug A"),
                           paste(name,"2 µM Drug A"),
                           paste(name,"1 µM Drug A"),
                           paste(name,"8 µM Drug B"),
                           paste(name,"4 µM Drug B"),
                           paste(name,"2 µM Drug B"),
                           paste(name,"1 µM Drug B"),
                           paste(name,"vehicle"),
                           paste(name,"empty"))
  
  ##Export to excel file below ##Necessary? i think this could be really helpful for whoever goes back to this? Maybe just print out as is?
  #Plate_sheet <- createWorkbook()
  #addWorksheet(Plate_sheet, "Raw data matrix")
  #addWorksheet(Plate_sheet, "Data frame")
  #addWorksheet(Plate_sheet, "Summary statistics")
  #writeData(Plate_sheet, sheet = "Raw data matrix", x = RLUmatrix, rowNames = TRUE, keepNA = TRUE)
  #writeData(Plate_sheet, sheet = "Data frame", x = Plate, rowNames = FALSE, keepNA = TRUE)
  #writeData(Plate_sheet, sheet = "Summary statistics", x = Plate_summary, rowNames = FALSE, keepNA = TRUE)
  #FN <- paste(name, ".xlsx", sep = "")
  #saveWorkbook(Plate_sheet, file = FN)
}

#All_plates <- All_plates %>% drop_na(Drug)
#Separate All plates and summary data? New function for graphing, just work from All_plates, will be perfectly fine to do and out of the same loop
#