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
  All_plates.outliers <<- x %>% 
    filter(Drug_name == i) %>% 
    group_by(Dose) %>% 
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

plot_All_plates_data <- function(x) {
  m <<- ggplot(data = x, 
               mapping = aes(x = factor(Dose, levels = c(
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
    labs(title = paste("Treatment of ",i, sep=""), 
         subtitle = paste(Outlier.count, " total outliers removed", sep = ""))
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
                     ChemDiv_rescreen_plate = numeric())

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
  Dose <- c(NA, "5 µM lomitapide", "8 µM", "4 µM", "2 µM", "1 µM", "0.5 µM", "0.25 µM", "0.125 µM", "0.0625 µM", "Vehicle", NA)
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
                                                                   extra = "drop")
  #Add to C...
  All_plates <- All_plates %>% add_row(Plate)
  }

All_plates <- All_plates %>% left_join(Library_decoder, by = c("ChemDiv_rescreen_plate"), keep = TRUE) %>%
  filter(Drug != "NA") %>%
  filter(RLU != "NA")

dir.create("./Output_data", showWarnings = FALSE)

Drug_list <- All_plates$Drug_name[!duplicated(All_plates$Drug_name)]

for (i in Drug_list) {
  if (nrow(All_plates %>% filter(Drug_name == i)) == 0) {next}
  
  dir.create(path = paste("./Output_data/",i, sep = ""), showWarnings = FALSE)
  All_plates.i <- All_plates %>% filter(Drug_name == i)
  Define_outliers(All_plates.i)
  Outlier.count <- All_plates.outliers %>% filter(Outlier == "Outlier") %>% nrow(.)
  Define_y_limit(All_plates.outliers)
    
  sink(file = paste("./Output_data/",i,"/",i,".txt", sep = ""))
  cat("ANOVA results\n")
  print(All_plates.outliers %>% 
          filter(Outlier == "Normal") %>% 
          compare_means(formula = RLU~Dose, 
                        method = "anova"))
  cat("\nPairwise t test results\n")
  print(All_plates.outliers %>% 
          filter(Outlier == "Normal") %>% 
          compare_means(formula = RLU~Dose,
                        ref.group = "Vehicle", 
                        p.adjust.method = "bonf", 
                        method = "t.test") %>% 
          arrange(factor(group2, levels = c(
            "5 µM lomitapide",
            "0.0625 µM",
            "0.125 µM",
            "0.25 µM",
            "0.5 µM",
            "1 µM",
            "2 µM",
            "4 µM",
            "8 µM"))))
  cat("\nns: p > 0.05\n*: p <= 0.05\n**: p <= 0.01\n***: p <= 0.001\n****: p <= 0.0001\n")
  sink(file = NULL)
  
  plot_All_plates_data(All_plates.outliers %>% filter(Outlier == "Normal"))
  ggsave(m, filename = paste(i, ".eps", sep = ""), width = 4, height = 5, path = paste("./Output_data/",i, sep = ""))
}









