# Load required libraries
library(dplyr)
library(readr)
library(base)

# Access the input file name passed as an argument
args <- commandArgs(trailingOnly = TRUE)
input_ini <- args[1]

# Read the content of the input file
new_ini = readLines("input_ini")
# new_ini = readLines("/Users/TingZhang/Documents/HEI HAQ PMF/CSN_IMPROVE_comp/CSN_CMD_txt/Cluster_1/Factor_6/iniparams_base.txt")

# Edit the 30th line of the file content
new_ini[30] = gsub("3890",
                   "000",
                   new_ini[30])

write.table(new_ini, 
            file = "ini_new.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
