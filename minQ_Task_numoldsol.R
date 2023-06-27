# Load required libraries
library(dplyr)
library(readr)

#### 1. task number when the lowest Q-value is generated ####
# Access the input file name passed as an argument
args <- commandArgs(trailingOnly = TRUE)
base_file <- args[1]

# Read the content of the input base file
base_output = readLines("base_file")

# Determine the number of lines to read
end.line = length(base_output) - 2 
start.line = length(base_output) - 21 

# Extract the lines including Q values and task numbers
Q_lines <- base_output[start.line:end.line]

# Convert the selected lines (with Q values & task number) into a data frame
Q_task <- read.table(text = Q_lines, 
                     header = F, 
                     quote = "\"'")[, 1:4]
colnames(Q_task) = c("TaskNum", "Qmain", "Qaux", "Q.XX")

# Find the number of task when the value of Qm is the lowest
lowest_Qm_task <- 
  Q_task %>% 
  filter(Qmain == min(Qmain)) %>% 
  select(TaskNum) %>% 
  pull()

# Output the the number of selected task
cat("The number of Task when the value of Qm is the lowest:", lowest_Qm_task)

#### 2. Update the value of numoldsol for other iniparams.txt files ####
# Define the file paths
param_files <- c("iniparams_BS.txt", 
                 "iniparams_DISP.txt", 
                 "iniparams_before_dual.txt", 
                 "iniparams_BS_DISP.txt"
                 )

# Replace the value of numoldsol in each file
for (file_path in param_files) {
  # Read the file
  param_lines <- read_lines(file_path)
  
  # Replace the value of numoldsol in the 15th line, default is now 0
  param_lines[15] = gsub("0", 
                         lowest_Qm_task, 
                         param_lines[15])
  
  # Write the updated lines back to the file
  write.table(param_lines, 
              file = file_path, 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}