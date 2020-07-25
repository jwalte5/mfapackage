# HEADER --------------------------------------------------------------
# SCRIPT NAME: Primary_Analysis_Template.R
#
# DATE CREATED: 28 June 2020
#
# AUTHOR: Joshua Walters
#
# DESCRIPTION: This script contains high level workflow for extracting
# data from XML files for the microplate feeder assay.
#
# NOTES
# -
# -
# SET WORKING DIRECTORY AND LOAD LIBRARIES INTO MEMORIES --------------
dir_path = choose.dir(caption = "Select Location of Data")
data_path = choose.files(caption = "Select Location of XML File")

dir_path = "E:\\CHG\\Work Directory\\Research\\Microplate Feeder Assay\\C - Plate Parsing"
data_path =  paste0(dir_path,"/Evap Test with double perforation.xml")

# Libraries
library("lubridate")
library("tidyverse")
library("xml2")
library("abind")
library("mfapackage")

# Use filepaths to set working directory, import functions, and
# import xml data
setwd(dir_path)
data = read_xml(data_path)

# Check all plate reads available
all_reads = import_plates(data)

# Calculate consumption time for each sample ID
line_reads = con_time(data)

# Import plate reads for calculations
plate_reads = read_plates(data)

# Calculate consumption
data_check = consumption_formula(vol = 10, plate1 = plate_reads[,,1], plate2 = plate_reads[,,2])

# Transforms samples into Tidy format instead of arrays
data_2check = collate_samples(data_check)

# Perform all calculations from xml document
layout_matrix = matrix(c("C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                         "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
                         nrow = 8,
                         ncol = 12,
                         byrow = TRUE)

data_df = consumption_from_xml(data, matrix_type = "custom", layout_matrix = layout_matrix)

#
stats = cbind(get_stats(df4),get_stats(df5, keep_row_names = FALSE))

write.csv(stats, file = "./evap_stats.csv")

colnames(stats) = c("Hardtop", "Hardtop_w_tape")


data.all = data.frame(Hardtop = df4, Hardtop_w_tape = df5)

data.all = pivot_longer(data.all, cols = c("Hardtop", "Hardtop_w_tape"), names_to = "Sealing", values_to = "Vol")

ggplot(data.all)+
  geom_boxplot(aes(x = Sealing, y = Vol))+
  geom_point(aes(x = Sealing, y = Vol))+
  labs(x = NULL, y = "Consumption (?L)", title = "Microplate CAFE Consumption" )

ggplot(data.all)+
  geom_boxplot(aes(x = Sealing, y = Vol))+
  geom_point(aes(x = Sealing, y = Vol))+
  labs(x = NULL, y = "Consumption (?L)", title = "Microplate CAFE Consumption" )+
  ylim(-0.5,1)

dev.copy(svg, "Hardtop_evap_range.svg")
dev.off()








plot_con = function(X) {
  ggplot(data.all)+
    geom_boxplot(aes(x = Sealing, y = Vol))+
    geom_point(aes(x = Sealing, y = Vol))+
    labs(x = NULL, y = "Consumption (?L)", title = "Microplate CAFE Consumption" )+
    ylim(-10,0)
}






