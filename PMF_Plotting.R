##clear environment
# rm(list=ls())

##set working directory
setwd("/Users/TingZhang/Documents/HEI HAQ PMF/PMF_Results")
getwd()
data.dir <- "/Users/TingZhang/Documents/HEI HAQ PMF/PMF_Results"


setwd("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/National_SA_PMF/")
getwd()
data.dir <- "/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/National_SA_PMF/"


##packages in need
library(tidyr) # separate{tidyr}, gather{tidyr}, spread{tidyr},  spread is VIP function, str_split_fixed{stringr} is better than separate
library(dplyr)
library(stats) # aggregate{stats}, VIP function
library(scales) # percent{}
library(ggrepel)
library(ggsci)
library(ggpattern)
library(ggplot2)

# read species-class dataset
species_class = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/National_SA_PMF/CSN_Species_class_sub.csv")
species_class$X = NULL

#SPECIATE-Source species
Top_5species_source = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/SPECIATE_2018/SPECIATE_Top_5_All_species_source.csv")
Top_5Element_source = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/SPECIATE_2018/SPECIATE_Top_5_Element_species_source.csv")
Top_5noCsub_source = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/SPECIATE_2018/SPECIATE_Top_5_noCsub_species_source.csv")
Top_5species_source$X = Top_5Element_source$X = Top_5noCsub_source$X = NULL

##
them_text_speciesName = theme(axis.title.x = element_text(color="grey25", size = 12,
                                                          family = "Arial Unicode MS", 
                                                          vjust=0, margin=margin(0,0,0,300)), 
                              axis.title.y = element_text(color="grey25", size = 12,
                                                          family = "Arial Unicode MS", 
                                                          vjust=1, margin=margin(0,2,0,0)),
                              axis.text.x = element_text(color="grey25", size = 11,
                                                         family = "Arial Unicode MS",
                                                         angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                              axis.text.y = element_text(color="grey25", size = 11,
                                                         family = "Arial Unicode MS",
                                                         angle = 0, hjust = 0.5))
                              

####### for bar start from a value less than 0 ####### 
# https://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
require(scales)
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

####### setting 30 colors for plotting clusters ####### 
npg.color = c("firebrick3", "indianred2", "lightskyblue", 
              "royalblue4", "lightseagreen", "skyblue3", 
              "goldenrod3", "burlywood4", 
              "burlywood3", "sandybrown", 
              "darkcyan")

########### information of the Site/Cluster #################
cluster_info_all = read.csv("PMF_NoGUI_cluster/CSN_PMF_CMD_StrongWeakBad_Cluster.csv")
cluster_info_all = read.csv("PMF_NoGUI_NoCsub_cluster/CSN_noCsub_PMF_CMD_StrongWeakBad_Cluster.csv")
cluster_info_all = read.csv("PMF_NoGUI_NoCsub_NoExtreme_cluster/CSN_noCsub_noExtreme_PMF_CMD_StrongWeakBad_Cluster.csv")
cluster_info_all$X = NULL

#### NO longer used 
cluster_info_row_extreme.vs.no = select(cluster_info_all, 
                                        Finaly.Decision, cluster.row, cluster.row.1, X.1)
colnames(cluster_info_row_extreme.vs.no) = c("Cluster", "Row_allData", "Row_noExtreme", "Decrease_Ratio")

plot(cluster_info_row_extreme.vs.no$Row_allData, cluster_info_row_extreme.vs.no$Row_noExtreme,
     xlab = "Row_allData", ylab = "Row_noExtreme")
cluster_info_row_extreme.vs.no$Decrease_Ratio = 
  (cluster_info_row_extreme.vs.no$Row_allData-
     cluster_info_row_extreme.vs.no$Row_noExtreme)/
  cluster_info_row_extreme.vs.no$Row_allData
summary(cluster_info_row_extreme.vs.no$Decrease_Ratio)
mean(cluster_info_row_extreme.vs.no$Decrease_Ratio)
sd(cluster_info_row_extreme.vs.no$Decrease_Ratio)
#### NO longer used 


cluster_info_all = plyr::rename(
  cluster_info_all, 
  c("K." = "KIon",
    "Na." = "NaIon", 
    "NH4." = "NH4Ion",
    "NO3" = "NO3Ion",
    "SO4" = "SO4Ion",
    "PM25" = "PM2.5"))

org_conc = read.csv("PMF_NoGUI_NoCsub_cluster/CSN_noCsub_C_6_PMF_CMD.csv")
org_conc = read.csv("PMF_GUI_NoCsub_NoExtreme_cluster/CSN_noCsub_noExtreme_C_6_unc.csv")
org_conc$X = NULL
org_conc$Date = as.Date(org_conc$Date)

base_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_C_6_F_9_base_result.txt")
base_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_C_6_F_9_2011-17_base_result.txt")

base_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_noCsub_C_6_F_9_base_result.txt")
disp_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_noCsub_C_6_F_9_DISP_DISPres1.txt")

base_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_noCsub_noExtreme_C_6_F_9_base_result.txt")
disp_output = readLines("CSN_PMF_noGUI_noCsub_AllData/CSN_noCsub_noExtreme_C_6_F_9_DISP_DISPres1.txt")


cluster.No = 6
factor.No = 9
Factor.serial = paste0("Factor", 1:factor.No)

cluster_info = subset(cluster_info_all, 
                      Finaly.Decision == cluster.No)
cluster.data.row = cluster_info$cluster.row

####### Select weak & strong variables by the value, All variable
cluster.weak.strong = 
  colnames(cluster_info)[which(
    !is.na(cluster_info[1, 1:48]))]
cluster.w.s.count = length(cluster.weak.strong)

cluster.strong = 
  colnames(cluster_info)[which(
    cluster_info[1, 1:48] == 1)]
cluster.str.count = length(cluster.strong)

cluster.weak = 
  colnames(cluster_info)[which(
    cluster_info[1, 1:48] == 0)]
cluster.weak.count = length(cluster.weak)

####### Select weak & strong variables by the value, NO C-sub
cluster.weak.strong = 
  colnames(cluster_info)[which(
    !is.na(cluster_info[1, 1:41]))]
cluster.w.s.count = length(cluster.weak.strong)

cluster.strong = 
  colnames(cluster_info)[which(
    cluster_info[1, 1:41] == 1)]
cluster.str.count = length(cluster.strong)

cluster.weak = 
  colnames(cluster_info)[which(
    cluster_info[1, 1:41] == 0)]
cluster.weak.count = length(cluster.weak)

#############################################################################
########### NO C-Subgroups Prepare data generated from PMF No-GUI CMD########### 
#############################################################################
###### 1. BASE run ######
# Determine the number of lines to read
end.bs.Q = length(base_output) - 2 
start.bs.Q = length(base_output) - 21 

# Extract the lines including Q values and task numbers
Q_lines <- base_output[start.bs.Q : end.bs.Q]

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

# Extract Factor matrix AA & BB from base result
# the regular expression pattern
pattern_start <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task, "\\s+---------------------------")
pattern_end <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task+1, "\\s+---------------------------")

# detect the line number where the pattern appears
base.start.line = grep(pattern_start, base_output)
base.end.line = grep(pattern_end, base_output)

#### 1.1 base - minumn Q ####
# the base run results for the selected task
base_selectQtask = 
  base_output[
    (base.start.line+3):
      (base.end.line-4)]

#### 1.2 base - G correlation ####
# fitted G vs. reference G, regression
G.correl.start = 10 + (factor.No+2)*2 + 3
base_G_cor_txt = 
  base_selectQtask[
    G.correl.start : 
      (G.correl.start + factor.No - 1)]

# split the lines into individual elements
base_G_cor_txt = unlist(strsplit(base_G_cor_txt, "\\s+"))
# create a matrix to store the individual elements and convert to data.frame
base_G_cor = data.frame(
  matrix(base_G_cor_txt, 
         nrow = factor.No, 
         ncol = factor.No+2, 
         byrow = T))
base_G_cor$X1 = NULL
base_G_cor = mutate_all(base_G_cor, as.numeric)

# rename the columns and replace the first column values
colnames(base_G_cor)[1] = "Factors"
colnames(base_G_cor)[2:(factor.No+1)] = Factor.serial
base_G_cor$Factors = Factor.serial

#### 1.3 base - time series ####
ts.start = G.correl.start + factor.No + 4
base_ts_txt = 
  base_selectQtask[
    ts.start : 
      (ts.start + cluster.data.row - 1)]

# split the lines into individual elements
base_ts_txt = unlist(strsplit(base_ts_txt, "\\s+"))

# create a matrix to store the individual elements and convert to data.frame
base_ts = data.frame(
  matrix(base_ts_txt, 
         nrow = cluster.data.row, 
         ncol = factor.No+2, 
         byrow = T))
base_ts$X1 = NULL
base_ts = mutate_all(base_ts, as.numeric)

# rename the columns and replace the first column values
colnames(base_ts)[1] = "Serial.No"
colnames(base_ts)[2:(factor.No+1)] = Factor.serial

# match with SiteCode & Date
base_ts_date = base_ts_PM = base_ts_all = base_ts
base_ts_date$Date = org_conc$Date
base_ts_date$SiteCode = org_conc$SiteCode
base_ts_PM$PM2.5 = org_conc$conc_PM25
base_ts_PM$Date = org_conc$Date
base_ts_PM$SiteCode = org_conc$SiteCode
base_ts_date$Serial.No = base_ts_PM$Serial.No = NULL

base_ts_all$Date = org_conc$Date
base_ts_all$State = org_conc$State
base_ts_all$SiteCode = org_conc$SiteCode
base_ts_all$PM2.5 = org_conc$conc_PM25
base_ts_all$species.sum = rowSums(base_ts_all[, 2:(factor.No+1)])
sapply(base_ts_all, class)
cor(base_ts_all$PM2.5, base_ts_all$species.sum)

##### 1.3.1 base - time series plot ##### 

# gather the data for plotting
base_ts_plot = gather(base_ts_date, 
                      "Factor", 
                      "Contribution", 
                      -Date, -SiteCode)
summary(base_ts_plot)

##### 1.3.2 base - linear, factor contribution ##### 
# factor contribution - linear
# ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM[, c("PM2.5", paste0("Factor", 1:factor.No))])
ts_PM_lm = lm(PM2.5 ~ ., 
              data = subset(base_ts_PM, 
                            SiteCode == 
                              unique(base_ts_PM$SiteCode)[3])[
                                , c("PM2.5", paste0("Factor", 1:factor.No))])
ts_PM_lm_beta = summary(ts_PM_lm)$coefficients[, 1]
ts_PM_lm_beta = data.frame(ts_PM_lm_beta[2:length(ts_PM_lm_beta)])
colnames(ts_PM_lm_beta) = "lm.beta.site-SiteBased"
ts_PM_lm_beta$Factor = rownames(ts_PM_lm_beta)
ts_PM_lm_beta$Factor.contribution = (ts_PM_lm_beta$lm.beta.site/
                                       sum(ts_PM_lm_beta$lm.beta.site))*100
# keep three significant digits
ts_PM_lm_beta$Factor.contr = paste0(
  signif(ts_PM_lm_beta$Factor.contribution, 
         3),
  "%")

# ts_PM_lm_beta_cluster = ts_PM_lm_beta

##### 1.4 base - conc & percent ##### 
conc.start = ts.start + cluster.data.row + 4
base_conc_txt = 
  base_selectQtask[
    conc.start : 
      (conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
base_conc_txt = unlist(strsplit(base_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
base_conc = data.frame(
  matrix(base_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+2, 
         byrow = T))
sapply(base_conc, class)
base_conc = mutate_all(base_conc, as.numeric)
base_conc$X1 = NULL

# rename the columns and replace the first column values
colnames(base_conc)[1] = "Species"
colnames(base_conc)[2:(factor.No+1)] = Factor.serial
base_conc$Species = cluster.weak.strong

# get the percent contribution - percent of species
base_percent_species = data.frame(Species = base_conc$Species)
base_percent_value = 
  signif(
    base_conc[, -1] *100 / 
      rowSums(base_conc[, -1]), 
    2)
base_percent = cbind(base_percent_species, 
                     base_percent_value)

######### JUMP how PMF calculate % of total variable does not work for SPECIATE ###########
### below is how PMF calculate % of total variable, but it's not always working
## not working reason 1, there are 0 in PM2.5 in base_conc
# base_percent_variable[Factor.serial] <- base_conc[Factor.serial] / subset(base_conc, Species == 'PM2.5')[Factor.serial]
base_percent_variable = base_percent_species
base_conc_PM = subset(base_conc, Species == 'PM2.5')[Factor.serial]
base_percent_variable[Factor.serial] = 0
for(i in 1:factor.No){
  factor.i = paste0("Factor", i)
  if(base_conc_PM[, i] != 0) {
    base_percent_variable[1:(cluster.w.s.count-1), i+1] = 
      base_conc[1:(cluster.w.s.count-1), i+1] / base_conc_PM[, i] * 100
  }
}

## not working reason 2, the % could be > 100%
# here is an example for csnC6F9_Alldata_allVariables
base_conc$Factor1[1]
base_conc$Factor1[1]/base_conc$Factor1[32]*100
base_percent_variable$Factor1[1]

base_conc$Factor8[14]/base_conc$Factor8[32]*100

######### JUMP how PMF calculate % of total variable does not work for SPECIATE ###########

## use adjusted method to estimate % of total variable
base_percent_variable = base_percent_species
base_conc_species = base_conc[-cluster.w.s.count, ][Factor.serial]
base_percent_var =
  round(
    base_conc_species * 100/ 
      rep(colSums(base_conc_species), 
          each = cluster.w.s.count - 1), 
    2)
base_percent_var[cluster.w.s.count, ] = 0
base_percent_variable = cbind(base_percent_variable, base_percent_var)

base_conc_plot = gather(base_conc,
                   "Factor", 
                   "Concentration", 
                   -Species)

base_percent_plot = gather(base_percent,
                      "Factor", 
                      "Percent", 
                      -Species)

#### 2. DISP run ####
##### 2.1 DISP - conc & percent - down lim ##### 
disp.down.conc.start = 9
disp_down_conc_txt = 
  disp_output[
    disp.down.conc.start : 
      (disp.down.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_down_conc_txt = unlist(strsplit(disp_down_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_down_conc = data.frame(
  matrix(disp_down_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_down_conc, class)
disp_down_conc = mutate_all(disp_down_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_down_conc)[1] = "Species"
colnames(disp_down_conc)[2:(factor.No+1)] = Factor.serial
disp_down_conc$Species = cluster.weak.strong

# get the percent contribution
disp_down_percent_species = data.frame(Species = disp_down_conc$Species)
disp_down_percent_value = 
  signif(
    disp_down_conc[, -1] *100 / 
      rowSums(disp_down_conc[, -1]), 
    2)
disp_down_percent = cbind(disp_down_percent_species, 
                     disp_down_percent_value)

disp_down_percent_plot = gather(disp_down_percent,
                                "Factor", 
                                "Percent.down", 
                                -Species)
summary(disp_down_percent_plot)


##### 2.2 DISP - conc & percent - up lim ##### 
disp.up.conc.start = disp.down.conc.start + cluster.w.s.count + 2
disp_up_conc_txt = 
  disp_output[
    disp.up.conc.start : 
      (disp.up.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_up_conc_txt = unlist(strsplit(disp_up_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_up_conc = data.frame(
  matrix(disp_up_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_up_conc, class)
disp_up_conc = mutate_all(disp_up_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_up_conc)[1] = "Species"
colnames(disp_up_conc)[2:(factor.No+1)] = Factor.serial
disp_up_conc$Species = cluster.weak.strong

# get the percent contribution
disp_up_percent_species = data.frame(Species = disp_up_conc$Species)
disp_up_percent_value = 
  signif(
    disp_up_conc[, -1] *100 / 
      rowSums(disp_up_conc[, -1]), 
    2)
disp_up_percent = cbind(disp_up_percent_species, 
                          disp_up_percent_value)

disp_up_percent_plot = gather(disp_up_percent,
                              "Factor", 
                              "Percent.up", 
                              -Species)
summary(disp_up_percent_plot)

#### plotting - Concentration & percent contribution ####
conc_percent_bsDisp = merge(base_conc_plot, base_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_down_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_up_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, species_class)

sapply(conc_percent_bsDisp, class)
head(conc_percent_bsDisp)
unique(conc_percent_bsDisp$Species)
conc_percent_bsDisp = subset(conc_percent_bsDisp,
                             Species != "PM25")
conc_percent_bsDisp_noCsub = conc_percent_bsDisp

source_first = read.csv("percent_source_determin_C_6_F_9_noCsub.csv")
source_first = select(source_first, Percent.data, Source1:Source9)
colnames(source_first)[2:10] = paste0("Factor", 1:factor.No)
source_first_use = data.frame(Factor = paste0("Factor", 1:factor.No),
                              Source = NA,
                              Source.No = NA)
source_first_use$Source[source_first_use$Factor == "Factor1"] = "F2-Secondary Nitrate"
source_first_use$Source[source_first_use$Factor == "Factor2"] = "F8-Biomass"
source_first_use$Source[source_first_use$Factor == "Factor3"] = "F9-Soil/Dust"
source_first_use$Source[source_first_use$Factor == "Factor4"] = "F3-Secondary Sulfate"
source_first_use$Source[source_first_use$Factor == "Factor5"] = "F4-Aged Sea Salt"
source_first_use$Source[source_first_use$Factor == "Factor6"] = "F7-"
source_first_use$Source[source_first_use$Factor == "Factor7"] = "F6-Fresh Sea Salt"
source_first_use$Source[source_first_use$Factor == "Factor8"] = "F5-Industry"
source_first_use$Source[source_first_use$Factor == "Factor9"] = "F1-Vehicle"

source_first_use$Source.No = sapply(strsplit(source_first_use$Source, "-"), "[", 1)

conc_percent_bsDisp_noCsub = merge(conc_percent_bsDisp_noCsub, source_first_use)


ggplot(conc_percent_bsDisp_noCsub, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species"
  ) +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_bsDisp_noCsub, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), 
              color = "black", shape = 15) +
  geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), 
                width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species",
    limits = c(-5, 110)) +
  theme_classic() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_bsDisp_noCsub, 
       aes(x = reorder(Species, sequence), 
           group = Source, fill = Source)) +
  facet_grid(Source ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6, alpha = 0.8)+
  scale_y_log10() + 
  scale_y_continuous(trans = mylog_trans(base=10, from=-5),
                     limits = c(1E-5, max(conc_percent_bsDisp$Concentration))) +
  scale_fill_npg() +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_bsDisp_noCsub, 
       aes(x = reorder(Species, sequence), 
           group = Factor, fill = Factor)) +
  facet_grid(Factor ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6, alpha = 0.8)+
  scale_y_log10() + 
  scale_y_continuous(trans = mylog_trans(base=10, from=-5),
                     limits = c(1E-5, max(conc_percent_bsDisp$Concentration))) +
  scale_fill_npg() +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

#### plotting - time series ####
base_ts_plot_noCsub = base_ts_plot
base_ts_plot_noCsub = merge(base_ts_plot_noCsub, source_first_use)
head(base_ts_plot_noCsub)
dim(base_ts_plot_noCsub)

ggplot(subset(base_ts_plot_noCsub,
              SiteCode = unique(base_ts_plot_noCsub$SiteCode)[3]), 
       aes(x = Date, y = Contribution, 
           group = Source, color = Source)) +
  facet_grid(Source ~., scales = "free_y") +
  geom_line(linewidth = 0.3, alpha = 0.4)+
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_npg() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.x = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

#### plotting - factor percent contribution to PM2.5 ####
# stripe, https://coolbutuseless.github.io/package/ggpattern/articles/pattern-stripe.html
ts_PM_lm_beta_noCsub = ts_PM_lm_beta
ts_PM_lm_beta_noCsub = merge(ts_PM_lm_beta_noCsub, source_first_use)

ggplot(ts_PM_lm_beta_noCsub, 
       aes(x = Source.No, y = Factor.contribution, fill = Source.No)) +
  geom_bar(stat="identity", width = 0.2) +    
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0),
        axis.title.y = element_text(color="grey25", size = 13, angle = 90, hjust = 0.5, vjust = 0))

ggplot(ts_PM_lm_beta_noCsub, 
       aes(x = Source.No, y = Factor.contribution)) +
  geom_bar_pattern(aes(fill = Source.No),
                   stat = "identity",
                   pattern_color = "white",
                   pattern_fill = "white",
                   # alpha = 0.8,
                   width = 0.3)  +
  geom_text(aes(label = paste(Source, Factor.contr)), 
            size = 6, angle = 90, hjust = 0, vjust = -3,
            position = position_stack(vjust = 0)) + # start from same bottom level
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color="grey25", size = 16, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 16, angle = 90, hjust = 5, vjust = -0.5),
        axis.title.x = element_text(color="grey25", size = 22, angle = 180, hjust = 0.5),
        axis.title.y = element_text(color="grey25", size = 0, angle = -90, vjust = -1))

#############################################################################
########### NO C-Sub & NO Extreme --- Prepare data generated from PMF No-GUI CMD########### 
#############################################################################
###### 1. BASE run ######
# Determine the number of lines to read
end.bs.Q = length(base_output) - 2 
start.bs.Q = length(base_output) - 21 

# Extract the lines including Q values and task numbers
Q_lines <- base_output[start.bs.Q : end.bs.Q]

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

# Extract Factor matrix AA & BB from base result
# the regular expression pattern
pattern_start <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task, "\\s+---------------------------")
pattern_end <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task+1, "\\s+---------------------------")

# detect the line number where the pattern appears
base.start.line = grep(pattern_start, base_output)
base.end.line = grep(pattern_end, base_output)

#### 1.1 base - minumn Q ####
# the base run results for the selected task
base_selectQtask = 
  base_output[
    (base.start.line+3):
      (base.end.line-4)]

#### 1.2 base - G correlation ####
# fitted G vs. reference G, regression
G.correl.start = 10 + (factor.No+2)*2 + 3
base_G_cor_txt = 
  base_selectQtask[
    G.correl.start : 
      (G.correl.start + factor.No - 1)]

# split the lines into individual elements
base_G_cor_txt = unlist(strsplit(base_G_cor_txt, "\\s+"))
# create a matrix to store the individual elements and convert to data.frame
base_G_cor = data.frame(
  matrix(base_G_cor_txt, 
         nrow = factor.No, 
         ncol = factor.No+2, 
         byrow = T))
base_G_cor$X1 = NULL
base_G_cor = mutate_all(base_G_cor, as.numeric)

# rename the columns and replace the first column values
colnames(base_G_cor)[1] = "Factors"
colnames(base_G_cor)[2:(factor.No+1)] = Factor.serial
base_G_cor$Factors = Factor.serial

#### 1.3 base - time series ####
ts.start = G.correl.start + factor.No + 4
base_ts_txt = 
  base_selectQtask[
    ts.start : 
      (ts.start + cluster.data.row - 1)]

# split the lines into individual elements
base_ts_txt = unlist(strsplit(base_ts_txt, "\\s+"))

# create a matrix to store the individual elements and convert to data.frame
base_ts = data.frame(
  matrix(base_ts_txt, 
         nrow = cluster.data.row, 
         ncol = factor.No+2, 
         byrow = T))
base_ts$X1 = NULL
base_ts = mutate_all(base_ts, as.numeric)

# rename the columns and replace the first column values
colnames(base_ts)[1] = "Serial.No"
colnames(base_ts)[2:(factor.No+1)] = Factor.serial

# match with SiteCode & Date
base_ts_date = base_ts_PM = base_ts_all = base_ts
base_ts_date$Date = org_conc$Date
base_ts_date$SiteCode = org_conc$SiteCode
base_ts_PM$PM2.5 = org_conc$PM25
base_ts_PM$Date = org_conc$Date
base_ts_PM$SiteCode = org_conc$SiteCode
base_ts_date$Serial.No = base_ts_PM$Serial.No = NULL

base_ts_all$Date = org_conc$Date
base_ts_all$State = org_conc$State
base_ts_all$SiteCode = org_conc$SiteCode
base_ts_all$PM2.5 = org_conc$PM25
base_ts_all$species.sum = rowSums(base_ts_all[, 2:(factor.No+1)])
sapply(base_ts_all, class)
cor(base_ts_all$PM2.5, base_ts_all$species.sum)

##### 1.3.1 base - time series plot ##### 
# gather the data for plotting
base_ts_plot = gather(base_ts_date, 
                      "Factor", 
                      "Contribution", 
                      -Date, -SiteCode)
summary(base_ts_plot)

##### 1.3.2 base - linear, factor contribution ##### 
# factor contribution - linear
# ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM)
ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM[, c("PM2.5", paste0("Factor", 1:factor.No))])
ts_PM_lm_beta = summary(ts_PM_lm)$coefficients[, 1]
ts_PM_lm_beta = data.frame(ts_PM_lm_beta[2:length(ts_PM_lm_beta)])
colnames(ts_PM_lm_beta) = "lm.beta.site-SiteBased"
ts_PM_lm_beta$Factor = rownames(ts_PM_lm_beta)
ts_PM_lm_beta$Factor.contribution = (ts_PM_lm_beta$lm.beta.site/
                                       sum(ts_PM_lm_beta$lm.beta.site))*100
# keep three significant digits
ts_PM_lm_beta$Factor.contr = paste0(
  signif(ts_PM_lm_beta$Factor.contribution, 
         3),
  "%")

##### 1.4 base - conc & percent ##### 
conc.start = ts.start + cluster.data.row + 4
base_conc_txt = 
  base_selectQtask[
    conc.start : 
      (conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
base_conc_txt = unlist(strsplit(base_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
base_conc = data.frame(
  matrix(base_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+2, 
         byrow = T))
sapply(base_conc, class)
base_conc = mutate_all(base_conc, as.numeric)
base_conc$X1 = NULL

# rename the columns and replace the first column values
colnames(base_conc)[1] = "Species"
colnames(base_conc)[2:(factor.No+1)] = Factor.serial
base_conc$Species = cluster.weak.strong

# get the percent contribution
base_percent_species = data.frame(Species = base_conc$Species)
base_percent_value = 
  signif(
    base_conc[, -1] *100 / 
      rowSums(base_conc[, -1]), 
    2)
base_percent = cbind(base_percent_species, 
                     base_percent_value)

## use adjusted method to estimate % of total variable
base_percent_variable = base_percent_species
base_conc_species = base_conc[-cluster.w.s.count, ][Factor.serial]
base_percent_var =
  round(
    base_conc_species * 100/ 
      rep(colSums(base_conc_species), 
          each = cluster.w.s.count - 1), 
    2)
base_percent_var[cluster.w.s.count, ] = 0
base_percent_variable = cbind(base_percent_variable, base_percent_var)

base_conc_plot = gather(base_conc,
                        "Factor", 
                        "Concentration", 
                        -Species)

base_percent_plot = gather(base_percent,
                           "Factor", 
                           "Percent", 
                           -Species)

#### 2. DISP run ####
##### 2.1 DISP - conc & percent - down lim ##### 
disp.down.conc.start = 9
disp_down_conc_txt = 
  disp_output[
    disp.down.conc.start : 
      (disp.down.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_down_conc_txt = unlist(strsplit(disp_down_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_down_conc = data.frame(
  matrix(disp_down_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_down_conc, class)
disp_down_conc = mutate_all(disp_down_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_down_conc)[1] = "Species"
colnames(disp_down_conc)[2:(factor.No+1)] = Factor.serial
disp_down_conc$Species = cluster.weak.strong

# get the percent contribution
disp_down_percent_species = data.frame(Species = disp_down_conc$Species)
disp_down_percent_value = 
  signif(
    disp_down_conc[, -1] *100 / 
      rowSums(disp_down_conc[, -1]), 
    2)
disp_down_percent = cbind(disp_down_percent_species, 
                          disp_down_percent_value)

disp_down_percent_plot = gather(disp_down_percent,
                                "Factor", 
                                "Percent.down", 
                                -Species)
summary(disp_down_percent_plot)


##### 2.2 DISP - conc & percent - up lim ##### 
disp.up.conc.start = disp.down.conc.start + cluster.w.s.count + 2
disp_up_conc_txt = 
  disp_output[
    disp.up.conc.start : 
      (disp.up.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_up_conc_txt = unlist(strsplit(disp_up_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_up_conc = data.frame(
  matrix(disp_up_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_up_conc, class)
disp_up_conc = mutate_all(disp_up_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_up_conc)[1] = "Species"
colnames(disp_up_conc)[2:(factor.No+1)] = Factor.serial
disp_up_conc$Species = cluster.weak.strong

# get the percent contribution
disp_up_percent_species = data.frame(Species = disp_up_conc$Species)
disp_up_percent_value = 
  signif(
    disp_up_conc[, -1] *100 / 
      rowSums(disp_up_conc[, -1]), 
    2)
disp_up_percent = cbind(disp_up_percent_species, 
                        disp_up_percent_value)

disp_up_percent_plot = gather(disp_up_percent,
                              "Factor", 
                              "Percent.up", 
                              -Species)
summary(disp_up_percent_plot)

#### plotting - Concentration & percent contribution ####
conc_percent_bsDisp = merge(base_conc_plot, base_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_down_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_up_percent_plot)
conc_percent_bsDisp = merge(conc_percent_bsDisp, species_class)

sapply(conc_percent_bsDisp, class)
head(conc_percent_bsDisp)
unique(conc_percent_bsDisp$Species)
conc_percent_bsDisp = subset(conc_percent_bsDisp,
                             Species != "PM25")


source_first = read.csv("percent_source_determin_C_6_F_9_noCsub_noExtreme.csv")
source_first = select(source_first, Percent.data, Source1:Source9)
colnames(source_first)[2:10] = paste0("Factor", 1:factor.No)
source_first_use = data.frame(Factor = paste0("Factor", 1:factor.No),
                              Source = NA,
                              Source.No = NA)

source_first_use$Source[source_first_use$Factor == "Factor1"] = "F9-Soil/Dust"
source_first_use$Source[source_first_use$Factor == "Factor2"] = "F4-Aged Sea Salt"
source_first_use$Source[source_first_use$Factor == "Factor3"] = "F7-"
source_first_use$Source[source_first_use$Factor == "Factor4"] = "F5-Industry"
source_first_use$Source[source_first_use$Factor == "Factor5"] = "F1-Vehicle"
source_first_use$Source[source_first_use$Factor == "Factor6"] = "F2-Secondary Nitrate"
source_first_use$Source[source_first_use$Factor == "Factor7"] = "F8-Biomass"
source_first_use$Source[source_first_use$Factor == "Factor8"] = "F3-Secondary Sulfate"
source_first_use$Source[source_first_use$Factor == "Factor9"] = "F6-Fresh Sea Salt"

source_first_use$Source.No = sapply(strsplit(source_first_use$Source, "-"), "[", 1)

conc_percent_bsDisp_noCsub_noExtreme = conc_percent_bsDisp
conc_percent_bsDisp_noCsub_noExtreme = merge(conc_percent_bsDisp_noCsub_noExtreme, source_first_use)


ggplot(conc_percent_bsDisp_noCsub_noExtreme, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species"
  ) +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_bsDisp_noCsub_noExtreme, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species",
    limits = c(-5, 110)) +
  theme_classic() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_bsDisp_noCsub_noExtreme, 
       aes(x = reorder(Species, sequence), 
           group = Source, fill = Source)) +
  facet_grid(Source ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6, alpha = 0.8)+
  scale_y_log10() + 
  scale_y_continuous(trans = mylog_trans(base=10, from=-5),
                     limits = c(1E-5, max(conc_percent_bsDisp$Concentration))) +
  scale_fill_npg() +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

#### plotting - time series ####
base_ts_plot_noCsub_noExtreme = base_ts_plot
base_ts_plot_noCsub_noExtreme = merge(base_ts_plot_noCsub_noExtreme, source_first_use)

ggplot(subset(base_ts_plot_noCsub_noExtreme,
              SiteCode = unique(base_ts_plot_noCsub_noExtreme$SiteCode)[3]), 
       aes(x = Date, y = Contribution, 
           group = Source, color = Source)) +
  facet_grid(Source ~., scales = "free_y") +
  geom_line(linewidth = 0.3, alpha = 0.8)+
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_npg() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.x = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

#### plotting - factor percent contribution to PM2.5 ####
# stripe, https://coolbutuseless.github.io/package/ggpattern/articles/pattern-stripe.html
ts_PM_lm_beta_noCsub_noExtreme = ts_PM_lm_beta
ts_PM_lm_beta_noCsub_noExtreme = merge(ts_PM_lm_beta_noCsub_noExtreme, source_first_use)

ggplot(ts_PM_lm_beta_noCsub_noExtreme, 
       aes(x = Source.No, y = Factor.contribution)) +
  geom_bar_pattern(aes(fill = Source.No),
                   stat = "identity",
                   pattern_color = "white",
                   pattern_fill = "white",
                   # alpha = 0.8,
                   width = 0.3)  +
  geom_text(aes(label = paste(Source, Factor.contr)), 
            size = 6, angle = 90, hjust = 0, vjust = -3,
            position = position_stack(vjust = 0)) + # start from same bottom level
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color="grey25", size = 16, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 16, angle = 90, hjust = 5, vjust = -0.5),
        axis.title.x = element_text(color="grey25", size = 22, angle = 180, hjust = 0.5),
        axis.title.y = element_text(color="grey25", size = 0, angle = -90, vjust = -1))

ggplot(ts_PM_lm_beta_noCsub_noExtreme, 
       aes(x = Source, y = Factor.contribution, fill = Source)) +
  geom_bar(stat="identity", width = 0.2) +    
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0),
        axis.title.y = element_text(color="grey25", size = 13, angle = 90, hjust = 0.5, vjust = 0))



#############################################################################
########### ALL Variables -- Prepare data generated from PMF No-GUI CMD########### 
#############################################################################
###### 1. BASE run ######
# Determine the number of lines to read
end.bs.Q = length(base_output) - 2 
start.bs.Q = length(base_output) - 21 

# Extract the lines including Q values and task numbers
Q_lines <- base_output[start.bs.Q : end.bs.Q]

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

# Extract Factor matrix AA & BB from base result
# the regular expression pattern
pattern_start <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task, "\\s+---------------------------")
pattern_end <- paste0("Results written by postprocessing for task #\\s+", lowest_Qm_task+1, "\\s+---------------------------")

# detect the line number where the pattern appears
base.start.line = grep(pattern_start, base_output)
base.end.line = grep(pattern_end, base_output)

#### 1.1 base - minumn Q ####
# the base run results for the selected task
base_selectQtask = 
  base_output[
    (base.start.line+3):
      (base.end.line-4)]

#### 1.2 base - G correlation --- NOT ALWAYS DETECTED ####
# fitted G vs. reference G, regression
G.correl.start = 10 + (factor.No+2)*2 + 3
base_G_cor_txt = 
  base_selectQtask[
    G.correl.start : 
      (G.correl.start + factor.No - 1)]

# split the lines into individual elements
base_G_cor_txt = unlist(strsplit(base_G_cor_txt, "\\s+"))
# create a matrix to store the individual elements and convert to data.frame
base_G_cor = data.frame(
  matrix(base_G_cor_txt, 
         nrow = factor.No, 
         ncol = factor.No+2, 
         byrow = T))
base_G_cor$X1 = NULL
base_G_cor = mutate_all(base_G_cor, as.numeric)

# rename the columns and replace the first column values
colnames(base_G_cor)[1] = "Factors"
colnames(base_G_cor)[2:(factor.No+1)] = Factor.serial
base_G_cor$Factors = Factor.serial

#### 1.3 base - time series ####
ts.start = 10 + (factor.No+2)*2 + 4
ts.start = 10 + (factor.No+2)*3 + 5 # 2011-17
cluster.data.row = 2934  # 2011-17

base_ts_txt = 
  base_selectQtask[
    ts.start : 
      (ts.start + cluster.data.row - 1)]

# split the lines into individual elements
base_ts_txt = unlist(strsplit(base_ts_txt, "\\s+"))

# create a matrix to store the individual elements and convert to data.frame
base_ts = data.frame(
  matrix(base_ts_txt, 
         nrow = cluster.data.row, 
         ncol = factor.No+2, 
         byrow = T))
base_ts$X1 = NULL
base_ts = mutate_all(base_ts, as.numeric)

# rename the columns and replace the first column values
colnames(base_ts)[1] = "Serial.No"
colnames(base_ts)[2:(factor.No+1)] = Factor.serial

# match with SiteCode & Date
base_ts_date = base_ts_PM = base_ts_all = base_ts

org_conc = subset(org_conc, Date < as.Date("2018-01-01")) # 2011-17

base_ts_date$Date = org_conc$Date
base_ts_date$SiteCode = org_conc$SiteCode
base_ts_PM$PM2.5 = org_conc$conc_PM25
base_ts_PM$Date = org_conc$Date
base_ts_PM$SiteCode = org_conc$SiteCode
base_ts_date$Serial.No = base_ts_PM$Serial.No = NULL

base_ts_all$Date = org_conc$Date
base_ts_all$State = org_conc$State
base_ts_all$SiteCode = org_conc$SiteCode
base_ts_all$PM2.5 = org_conc$conc_PM25
base_ts_all$species.sum = rowSums(base_ts_all[, 2:(factor.No+1)])
sapply(base_ts_all, class)
cor(base_ts_all$PM2.5, base_ts_all$species.sum)

##### 1.3.1 base - time series plot ##### 

# gather the data for plotting
base_ts_plot = gather(base_ts_date, 
                      "Factor", 
                      "Contribution", 
                      -Date, -SiteCode)
summary(base_ts_plot)

##### 1.3.2 base - linear, factor contribution ##### 
# factor contribution - linear
# ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM)
ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM[, c("PM2.5", paste0("Factor", 1:factor.No))])

base_ts_PM_17_in_20 = subset(base_ts_PM, Date < as.Date("2018-01-01"))
ts_PM_lm = lm(PM2.5 ~ ., data = base_ts_PM_17_in_20[, c("PM2.5", paste0("Factor", 1:factor.No))])

ts_PM_lm_beta = summary(ts_PM_lm)$coefficients[, 1]
ts_PM_lm_beta = data.frame(ts_PM_lm_beta[2:length(ts_PM_lm_beta)])
colnames(ts_PM_lm_beta) = "lm.beta.site-SiteBased"
ts_PM_lm_beta$Factor = rownames(ts_PM_lm_beta)
ts_PM_lm_beta$Factor.contribution = (ts_PM_lm_beta$lm.beta.site/
                                       sum(ts_PM_lm_beta$lm.beta.site))*100
# keep three significant digits
ts_PM_lm_beta$Factor.contr = paste0(
  signif(ts_PM_lm_beta$Factor.contribution, 
         3),
  "%")

# ts_PM_lm_beta_17_in_20 = ts_PM_lm_beta
##### 1.4 base - conc & percent ##### 
conc.start = ts.start + cluster.data.row + 4
base_conc_txt = 
  base_selectQtask[
    conc.start : 
      (conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
base_conc_txt = unlist(strsplit(base_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
base_conc = data.frame(
  matrix(base_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+2, 
         byrow = T))
sapply(base_conc, class)
base_conc = mutate_all(base_conc, as.numeric)
base_conc$X1 = NULL

# rename the columns and replace the first column values
colnames(base_conc)[1] = "Species"
colnames(base_conc)[2:(factor.No+1)] = Factor.serial
base_conc$Species = cluster.weak.strong

# get the percent contribution
base_percent_species = data.frame(Species = base_conc$Species)
base_percent_value = 
  signif(
    base_conc[, -1] *100 / 
      rowSums(base_conc[, -1]), 
    2)
base_percent = cbind(base_percent_species, 
                     base_percent_value)

## use adjusted method to estimate % of total variable
base_percent_variable = base_percent_species
base_conc_species = base_conc[-cluster.w.s.count, ][Factor.serial]
base_percent_var =
  round(
    base_conc_species * 100/ 
      rep(colSums(base_conc_species), 
          each = cluster.w.s.count - 1), 
    2)
base_percent_var[cluster.w.s.count, ] = 0
base_percent_variable = cbind(base_percent_variable, base_percent_var)

base_conc_plot = gather(base_conc,
                        "Factor", 
                        "Concentration", 
                        -Species)

base_percent_plot = gather(base_percent,
                           "Factor", 
                           "Percent", 
                           -Species)

#### 2. DISP run ####
##### 2.1 DISP - conc & percent - down lim ##### 
disp.down.conc.start = 9
disp_down_conc_txt = 
  disp_output[
    disp.down.conc.start : 
      (disp.down.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_down_conc_txt = unlist(strsplit(disp_down_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_down_conc = data.frame(
  matrix(disp_down_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_down_conc, class)
disp_down_conc = mutate_all(disp_down_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_down_conc)[1] = "Species"
colnames(disp_down_conc)[2:(factor.No+1)] = Factor.serial
disp_down_conc$Species = cluster.weak.strong

# get the percent contribution
disp_down_percent_species = data.frame(Species = disp_down_conc$Species)
disp_down_percent_value = 
  signif(
    disp_down_conc[, -1] *100 / 
      rowSums(disp_down_conc[, -1]), 
    2)
disp_down_percent = cbind(disp_down_percent_species, 
                          disp_down_percent_value)

disp_down_percent_plot = gather(disp_down_percent,
                                "Factor", 
                                "Percent.down", 
                                -Species)
summary(disp_down_percent_plot)


##### 2.2 DISP - conc & percent - up lim ##### 
disp.up.conc.start = disp.down.conc.start + cluster.w.s.count + 2
disp_up_conc_txt = 
  disp_output[
    disp.up.conc.start : 
      (disp.up.conc.start + cluster.w.s.count - 1)]

# split the lines into individual elemenconc
disp_up_conc_txt = unlist(strsplit(disp_up_conc_txt, "\\s+"))

# create a matrix to store the individual elemenconc and convert to data.frame
disp_up_conc = data.frame(
  matrix(disp_up_conc_txt, 
         nrow = cluster.w.s.count, 
         ncol = factor.No+1, 
         byrow = T))
sapply(disp_up_conc, class)
disp_up_conc = mutate_all(disp_up_conc, as.numeric)

# rename the columns and replace the first column values
colnames(disp_up_conc)[1] = "Species"
colnames(disp_up_conc)[2:(factor.No+1)] = Factor.serial
disp_up_conc$Species = cluster.weak.strong

# get the percent contribution
disp_up_percent_species = data.frame(Species = disp_up_conc$Species)
disp_up_percent_value = 
  signif(
    disp_up_conc[, -1] *100 / 
      rowSums(disp_up_conc[, -1]), 
    2)
disp_up_percent = cbind(disp_up_percent_species, 
                        disp_up_percent_value)

disp_up_percent_plot = gather(disp_up_percent,
                              "Factor", 
                              "Percent.up", 
                              -Species)
summary(disp_up_percent_plot)

#### plotting - Concentration & percent contribution ####
conc_percent_bsDisp = merge(base_conc_plot, base_percent_plot)

# conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_down_percent_plot)
# conc_percent_bsDisp = merge(conc_percent_bsDisp, disp_up_percent_plot)
 
conc_percent_bsDisp = merge(conc_percent_bsDisp, species_class)
sapply(conc_percent_bsDisp, class)
head(conc_percent_bsDisp)
unique(conc_percent_bsDisp$Species)
conc_percent_bsDisp = subset(conc_percent_bsDisp,
                             Species != "PM25")
unique(conc_percent_bsDisp$Species)

source_first = read.csv("percent_source_determin_C_6_F_9_Alldata_AllVaraible.csv")
source_first = read.csv("percent_source_determin_C_6_F_9_Alldata_AllVaraible_2011-17.csv")

source_first = select(source_first, Percent.data, Source1:Source9)
colnames(source_first)[2:10] = paste0("Factor", 1:factor.No)

source_first_use = data.frame(Factor = paste0("Factor", 1:factor.No),
                              Source = NA,
                              Source.No = NA)

# "F1-Vehicle" "F2-Secondary Nitrate" "F3-Secondary Sulfate" 
# "F4-Aged Sea Salt"  "F5-Industry" "F6-Fresh Sea Salt" "F7-"  "F8-Biomass"
source_first_use$Source[source_first_use$Factor == "Factor1"] = "F1-"
source_first_use$Source[source_first_use$Factor == "Factor2"] = "F4-Aged Sea Salt"
source_first_use$Source[source_first_use$Factor == "Factor3"] = "F2-Secondary Nitrate"
source_first_use$Source[source_first_use$Factor == "Factor4"] = "F7-"
source_first_use$Source[source_first_use$Factor == "Factor5"] = "F5-"
source_first_use$Source[source_first_use$Factor == "Factor6"] = "F3-Secondary Sulfate"
source_first_use$Source[source_first_use$Factor == "Factor7"] = "F8-Biomass"
source_first_use$Source[source_first_use$Factor == "Factor8"] = "F9-Soil/Dust"
source_first_use$Source[source_first_use$Factor == "Factor9"] = "F6-Fresh Sea Salt"

# 2011-17
source_first_use$Source[source_first_use$Factor == "Factor1"] = "F6-Fresh Sea Salt (??)"
source_first_use$Source[source_first_use$Factor == "Factor2"] = "F3-Secondary Sulfate"
source_first_use$Source[source_first_use$Factor == "Factor3"] = "F8-Biomass"
source_first_use$Source[source_first_use$Factor == "Factor4"] = "F4-Aged Sea Salt"
source_first_use$Source[source_first_use$Factor == "Factor5"] = "F1-"
source_first_use$Source[source_first_use$Factor == "Factor6"] = "F9-Soil/Dust "
source_first_use$Source[source_first_use$Factor == "Factor7"] = "F5-"
source_first_use$Source[source_first_use$Factor == "Factor8"] = "F2-Secondary Nitrate"
source_first_use$Source[source_first_use$Factor == "Factor9"] = "F7-"

source_first_use$Source.No = sapply(strsplit(source_first_use$Source, "-"), "[", 1)

conc_percent_bsDisp_AllVariable = conc_percent_bsDisp
conc_percent_bsDisp_AllVariable = merge(conc_percent_bsDisp_AllVariable, source_first_use)

conc_percent_2011.17_AllVariable = conc_percent_bsDisp
conc_percent_2011.17_AllVariable = merge(conc_percent_2011.17_AllVariable, source_first_use)

ggplot(conc_percent_2011.17_AllVariable, #conc_percent_bsDisp_AllVariable, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species"
  ) +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

ggplot(conc_percent_2011.17_AllVariable, #conc_percent_bsDisp_AllVariable, 
       aes(x = reorder(Species, sequence), 
           group = Source)) +
  facet_grid(Source ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  # geom_errorbar(aes(ymin = Percent.down, ymax = Percent.up), width = 0.4) +
  scale_y_continuous(
    position = "right",
    name = "% of Species",
    limits = c(-5, 110)) +
  theme_classic() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

ggplot(conc_percent_2011.17_AllVariable, #conc_percent_bsDisp_AllVariable, 
       aes(x = reorder(Species, sequence), 
           group = Source, fill = Source)) +
  facet_grid(Source ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6, alpha = 0.8)+
  scale_y_log10() + 
  scale_y_continuous(trans = mylog_trans(base=10, from=-5),
                     limits = c(1E-5, max(conc_percent_bsDisp$Concentration))) +
  scale_fill_npg() +
  theme_bw() +
  xlab(format_variable("PM25 Species")) +
  scale_x_discrete(labels = function(x) format_variable(x)) +
  them_text_speciesName +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")

#### plotting - time series ####
base_ts_plot_AllVariable = base_ts_plot
base_ts_plot_AllVariable = merge(base_ts_plot_AllVariable, source_first_use)

base_ts_2011.17_AllVariable = base_ts_plot
base_ts_2011.17_AllVariable = merge(base_ts_2011.17_AllVariable, source_first_use)

base_ts_plot_AllVariable_20 = base_ts_plot_AllVariable
base_ts_plot_AllVariable_17 = base_ts_2011.17_AllVariable

#base_ts_plot_AllVariable_20$period = "2011-20"
#base_ts_plot_AllVariable_17$period = "2011-17"

#colnames(base_ts_plot_AllVariable_20)[4] = "Contribution_20"
#colnames(base_ts_plot_AllVariable_17)[4] = "Contribution_17"
colnames(base_ts_plot_AllVariable_20) = c("Factor_20", "Date", "SiteCode", 
                                          "Contribution_20", "Source_20", "Source.No_20")
colnames(base_ts_plot_AllVariable_17) = c("Factor_17", "Date", "SiteCode", 
                                          "Contribution_17", "Source_17", "Source.No_17")

base_ts_plot_AllVariable_20 = subset(base_ts_plot_AllVariable_20,
                                     !(Source.No_20 %in% c("F1", "F5", "F7")))
base_ts_plot_AllVariable_17 = subset(base_ts_plot_AllVariable_17,
                                     !(Source.No_17 %in% c("F1", "F5", "F7")))
base_ts_plot_AllVariable_20 = subset(base_ts_plot_AllVariable_20,
                                     Date < as.Date("2018-01-01"))

base_ts_plot_AllVariable_20$Factor = base_ts_plot_AllVariable_17$Factor = NULL
base_ts_plot_AllVariable_20$Source_20[
  base_ts_plot_AllVariable_20$Source_20 == "F6-Fresh Sea Salt (??)"] = 
  "F6-Fresh Sea Salt"
  
base_ts_plot_AllVariable_20 = base_ts_plot_AllVariable_20[with(base_ts_plot_AllVariable_20, 
                                                               order(SiteCode, Date, Source_20)), ]
base_ts_plot_AllVariable_17 = base_ts_plot_AllVariable_17[with(base_ts_plot_AllVariable_17, 
                                                               order(SiteCode, Date, Source_17)), ]

summary(base_ts_plot_AllVariable_20$Date == base_ts_plot_AllVariable_17$Date)
summary(base_ts_plot_AllVariable_20$SiteCode == base_ts_plot_AllVariable_17$SiteCode)
summary(base_ts_plot_AllVariable_20$Source.No_20 == base_ts_plot_AllVariable_17$Source.No_17)

base_ts_AllVariable_17_20_rb = rbind(base_ts_plot_AllVariable_17, 
                                     base_ts_plot_AllVariable_20)
base_ts_AllVariable_17_20_cb = cbind(base_ts_plot_AllVariable_17, 
                                     base_ts_plot_AllVariable_20)
cor(base_ts_AllVariable_17_20_cb$Contribution_17, base_ts_AllVariable_17_20_cb$Contribution_20)


ggplot(base_ts_2011.17_AllVariable, #base_ts_plot_AllVariable, 
       aes(x = Date, y = Contribution, 
           group = Source, color = Source)) +
  facet_grid(Source ~., scales = "free_y") +
  geom_line(alpha = 0.5, linewidth = 0.3)+
  geom_point(size = 0.3) +
  scale_color_npg() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.x = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

#### plotting - factor percent contribution to PM2.5 ####
# stripe, https://coolbutuseless.github.io/package/ggpattern/articles/pattern-stripe.html
ts_PM_lm_beta_AllVariable = ts_PM_lm_beta
ts_PM_lm_beta_AllVariable = merge(ts_PM_lm_beta_AllVariable, source_first_use)

ts_PM_lm_2011.17_AllVariable = ts_PM_lm_beta
ts_PM_lm_2011.17_AllVariable = merge(ts_PM_lm_2011.17_AllVariable, source_first_use)
ts_PM_lm_2011.17_AllVariable$Factor.contr[ts_PM_lm_2011.17_AllVariable$Factor == "Factor9"] = "0.446%"


ts_PM_lm_beta_17_in_20 = ts_PM_lm_beta
ts_PM_lm_beta_17_in_20 = merge(ts_PM_lm_beta_17_in_20, source_first_use)

ggplot(ts_PM_lm_beta_17_in_20, #ts_PM_lm_beta_AllVariable, 
       aes(x = Source.No, y = Factor.contribution)) +
  geom_bar_pattern(aes(fill = Source.No),
                   stat = "identity",
                   pattern_color = "white",
                   pattern_fill = "white",
                   # alpha = 0.8,
                   width = 0.3)  +
  geom_text(aes(label = paste(Source, Factor.contr)), 
            size = 6, angle = 90, hjust = 0, vjust = -3,
            position = position_stack(vjust = 0)) + # start from same bottom level
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color="grey25", size = 16, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 16, angle = 90, hjust = 5, vjust = -0.5),
        axis.title.x = element_text(color="grey25", size = 22, angle = 180, hjust = 0.5),
        axis.title.y = element_text(color="grey25", size = 0, angle = -90, vjust = -1))

ggplot(ts_PM_lm_beta, 
       aes(x = Factor, y = Factor.contribution)) +
  # geom_bar(stat="identity", width = 0.2) +    
  geom_col_pattern(aes(pattern_fill = Factor),  # fill = Factor
                   width = 0.2, fill = "white", 
                   pattern = 'stripe', pattern_density = 0.5) +
  scale_fill_nejm() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5))


ggplot(ts_PM_lm_beta, 
       aes(x = Factor, y = Factor.contribution, fill = Factor)) +
  geom_bar(stat="identity", width = 0.2) +    
  scale_fill_npg() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0),
        axis.title.y = element_text(color="grey25", size = 13, angle = 90, hjust = 0.5, vjust = 0))


#####################################################################################
############# First determination of sources based on two percent files ############# 
#####################################################################################
## generate dataframe to store all first-step source interpretation info
percent_source = setNames(data.frame(matrix(NA, ncol = 27, nrow = 0)), 
                          c("Dataset", "Version", "Cluster", 
                            "Factor.No", "Percent.data",
                            paste0("Factor", 1:11), 
                            paste0("Source", 1:11)))

##### 1. all species ######
base_percent_ofSpecies = base_percent[-cluster.w.s.count, ]
base_percent_ofVariable = base_percent_variable[-cluster.w.s.count, ]

## Detected the main species for each factor based on % of Species
rownames(base_percent_ofSpecies) = base_percent_ofSpecies$Species
base_percent_ofSpecies$Species = NULL

# get the rownames of which the %ofSpecies ranks top 5 of the column
main5_Species = data.frame(
  apply(
    base_percent_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_ofSpecies)
    [order(x, decreasing = T)[1:5]]
  ))

main3_Species = data.frame(
  apply(
    base_percent_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_ofSpecies)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_Species = apply(main5_Species, 2, 
                      function(x) 
                        paste(x, collapse = " "))
main3_Species = apply(main3_Species, 2, 
                      function(x) 
                        paste(x, collapse = " "))

## Detected the main Variable for each factor based on % of all Variable
rownames(base_percent_ofVariable) = base_percent_ofVariable$Species
base_percent_ofVariable$Species = NULL

# get the rownames of which the %ofVariable ranks top 5 of the column
main5_Variable = data.frame(
  apply(
    base_percent_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_ofVariable)
    [order(x, decreasing = T)[1:5]]
  ))

main3_Variable = data.frame(
  apply(
    base_percent_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_ofVariable)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_Variable = apply(main5_Variable, 2, 
                       function(x) 
                         paste(x, collapse = " "))
main3_Variable = apply(main3_Variable, 2, 
                       function(x) 
                         paste(x, collapse = " "))

## combine results based on % of Species & % of all Variable
main5_main3 = 
  do.call(rbind,
          list(
            main5_Species, 
            main3_Species,
            main5_Variable,
            main3_Variable))

#  
base_percent_source = data.frame(Dataset = rep("CSN", 4), 
                                 Version = rep("allVariable_AllData", 4), 
                                 Cluster = rep(cluster.No, 4), 
                                 Factor.No = rep(factor.No, 4),
                                 Percent.data = rep(c("%_Species", "%_Variable"),
                                                    each = 2))
base_percent_source = cbind(base_percent_source, main5_main3)

percent_source = rbind.fill(percent_source, base_percent_source)

##### 2. No C-subgroup ######
C.sub.PM = c("OC1", "OC2", "OC3", "OC4", 
          "EC1", "EC2", "EC3", "PM2.5")

base_percent_noCsub_ofSpecies = 
  subset(base_percent, 
         !(Species %in% C.sub.PM))
base_percent_noCsub_ofVariable = 
  subset(base_percent_variable, 
         !(Species %in% C.sub.PM))

## Detected the main species for each factor based on % of Species
rownames(base_percent_noCsub_ofSpecies) = base_percent_noCsub_ofSpecies$Species
base_percent_noCsub_ofSpecies$Species = NULL

# get the rownames of which the %ofSpecies ranks top 5 of the column
main5_noCsub_Species = data.frame(
  apply(
    base_percent_noCsub_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_noCsub_ofSpecies)
    [order(x, decreasing = T)[1:5]]
  ))

main3_noCsub_Species = data.frame(
  apply(
    base_percent_noCsub_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_noCsub_ofSpecies)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_noCsub_Species = apply(main5_noCsub_Species, 2, 
                             function(x) 
                               paste(x, collapse = " "))
main3_noCsub_Species = apply(main3_noCsub_Species, 2, 
                             function(x) 
                               paste(x, collapse = " "))

## Detected the main Variable for each factor based on % of all Variable
rownames(base_percent_noCsub_ofVariable) = base_percent_noCsub_ofVariable$Species
base_percent_noCsub_ofVariable$Species = NULL

# get the rownames of which the %ofVariable ranks top 5 of the column
main5_noCsub_Variable = data.frame(
  apply(
    base_percent_noCsub_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_noCsub_ofVariable)
    [order(x, decreasing = T)[1:5]]
  ))

main3_noCsub_Variable = data.frame(
  apply(
    base_percent_noCsub_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_noCsub_ofVariable)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_noCsub_Variable = apply(main5_noCsub_Variable, 2, 
                              function(x) 
                                paste(x, collapse = " "))
main3_noCsub_Variable = apply(main3_noCsub_Variable, 2, 
                              function(x) 
                                paste(x, collapse = " "))

## combine results based on % of Species & % of all Variable
main5_main3_noCsub = 
  do.call(rbind,
          list(
            main5_noCsub_Species, 
            main3_noCsub_Species,
            main5_noCsub_Variable,
            main3_noCsub_Variable))

# 
base_percent_noCsub_source = data.frame(Dataset = rep("CSN", 4), 
                                        Version = rep("noCsub_AllData", 4), 
                                        Cluster = rep(cluster.No, 4), 
                                        Factor.No = rep(factor.No, 4),
                                        Percent.data = rep(c("%_Species", "%_Variable"),
                                                           each = 2))
base_percent_noCsub_source = cbind(base_percent_noCsub_source, main5_main3_noCsub)

percent_source = rbind.fill(percent_source, base_percent_noCsub_source)

##### 2. ONLY element ######
C.sub.PM = c("OC", "OC1", "OC2", "OC3", "OC4", 
             "EC", "EC1", "EC2", "EC3", "PM2.5")

base_percent_element_ofSpecies = 
  subset(base_percent, 
         !(Species %in% C.sub.PM |
             grepl("Ion", Species, fixed = T)))
base_percent_element_ofVariable = 
  subset(base_percent_variable, 
         !(Species %in% C.sub.PM |
             grepl("Ion", Species, fixed = T)))

## Detected the main species for each factor based on % of Species
rownames(base_percent_element_ofSpecies) = base_percent_element_ofSpecies$Species
base_percent_element_ofSpecies$Species = NULL

# get the rownames of which the %ofSpecies ranks top 5 of the column
main5_element_Species = data.frame(
  apply(
    base_percent_element_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_element_ofSpecies)
    [order(x, decreasing = T)[1:5]]
  ))

main3_element_Species = data.frame(
  apply(
    base_percent_element_ofSpecies, 
    2, 
    function(x) 
      rownames(base_percent_element_ofSpecies)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_element_Species = apply(main5_element_Species, 2, 
                             function(x) 
                               paste(x, collapse = " "))
main3_element_Species = apply(main3_element_Species, 2, 
                             function(x) 
                               paste(x, collapse = " "))

## Detected the main Variable for each factor based on % of all Variable
rownames(base_percent_element_ofVariable) = base_percent_element_ofVariable$Species
base_percent_element_ofVariable$Species = NULL

# get the rownames of which the %ofVariable ranks top 5 of the column
main5_element_Variable = data.frame(
  apply(
    base_percent_element_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_element_ofVariable)
    [order(x, decreasing = T)[1:5]]
  ))

main3_element_Variable = data.frame(
  apply(
    base_percent_element_ofVariable, 
    2, 
    function(x) 
      rownames(base_percent_element_ofVariable)
    [order(x, decreasing = T)[1:3]]
  ))

# combine the selected rownames into one cell
main5_element_Variable = apply(main5_element_Variable, 2, 
                              function(x) 
                                paste(x, collapse = " "))
main3_element_Variable = apply(main3_element_Variable, 2, 
                              function(x) 
                                paste(x, collapse = " "))

## combine results based on % of Species & % of all Variable
main5_main3_element = 
  do.call(rbind,
          list(
            main5_element_Species, 
            main3_element_Species,
            main5_element_Variable,
            main3_element_Variable))

# 
base_percent_element_source = data.frame(Dataset = rep("CSN", 4), 
                                        Version = rep("elementOnly_AllData", 4), 
                                        Cluster = rep(cluster.No, 4), 
                                        Factor.No = rep(factor.No, 4),
                                        Percent.data = rep(c("%_Species", "%_Variable"),
                                                           each = 2))
base_percent_element_source = cbind(base_percent_element_source, main5_main3_element)

percent_source = rbind.fill(percent_source, base_percent_element_source)

##### 4. source assignment - SPECIATE ######
percent_source_variable = subset(percent_source, 
                                 Percent.data == "%_Variable")
percent_source_all = subset(percent_source_variable, 
                            Version == "allVariable_AllData")
percent_source_noCsub = subset(percent_source_variable, 
                               Version == "noCsub_AllData")
percent_source_element = subset(percent_source_variable, 
                                Version == "elementOnly_AllData")

Top_5noCsub_source$AircraftVessel = Top_5noCsub_source$Fuel = 
  Top_5noCsub_source$Agriculture = Top_5noCsub_source$DustCons = 
  Top_5noCsub_source$NotClass = NULL
Top_5Element_source$AircraftVessel = Top_5Element_source$Cook = 
  Top_5Element_source$DustCons = Top_5Element_source$NitrateSulfate = 
  Top_5Element_source$NotClass = NULL
Top_5species_source$Fuel = Top_5species_source$Agriculture = 
  Top_5species_source$NotClass = NULL
# Top_5species_source; Top_5Element_source; Top_5noCsub_source

source.All = ncol(Top_5species_source) - 1
source.noCsub = ncol(Top_5noCsub_source) - 1
source.element = ncol(Top_5Element_source) - 1

for(i in 1:11) {
  # colnames(percent_source_all)[6]
  if(!is.na(percent_source_all[1, i+6-1])){
    main.species.sample = percent_source_all[1, i+6-1]
    main.species.sample = as.character(str_split(main.species.sample, "\\s+")[[1]])
    for(l in 1:source.All) {
      source.name = colnames(Top_5species_source)[l+1]
      main.species.source = Top_5species_source[, l+1]
      percent_source_all[1, i+16] = 
        ifelse(
          sum(main.species.sample %in% main.species.source) >= 3,
          ifelse(is.na(percent_source_all[1, i+16]), 
                 source.name, 
                 paste(source.name, as.character(percent_source_all[1, i+16]))),
          percent_source_all[1, i+16])
    }
  }
}

for(i in 1:11) {
  # colnames(percent_source_noCsub)[6]
  if(!is.na(percent_source_noCsub[1, i+6-1])){
    main.species.sample = percent_source_noCsub[1, i+6-1]
    main.species.sample = as.character(str_split(main.species.sample, "\\s+")[[1]])
    for(l in 1:source.noCsub) {
      source.name = colnames(Top_5noCsub_source)[l+1]
      main.species.source = Top_5noCsub_source[, l+1]
      percent_source_noCsub[1, i+16] = 
        ifelse(
          sum(main.species.sample %in% main.species.source) >= 3,
          ifelse(is.na(percent_source_noCsub[1, i+16]), 
                 source.name, 
                 paste(source.name, as.character(percent_source_noCsub[1, i+16]))),
          percent_source_noCsub[1, i+16])
    }
  }
}

for(i in 1:11) {
  # colnames(percent_source_element)[6]
  if(!is.na(percent_source_element[1, i+6-1])){
    main.species.sample = percent_source_element[1, i+6-1]
    main.species.sample = as.character(str_split(main.species.sample, "\\s+")[[1]])
    for(l in 1:source.element) {
      source.name = colnames(Top_5Element_source)[l+1]
      main.species.source = Top_5Element_source[, l+1]
      percent_source_element[1, i+16] = 
        ifelse(
          sum(main.species.sample %in% main.species.source) >= 3,
          ifelse(is.na(percent_source_element[1, i+16]), 
                 source.name, 
                 paste(source.name, as.character(percent_source_element[1, i+16]))),
          percent_source_element[1, i+16])
    }
  }
}

percent_source_var_result = 
  do.call(rbind,
          list(
            percent_source_all, 
            percent_source_noCsub,
            percent_source_element))

##### 5.1 source assignment - %_of_Species - noCsub ######
percent_source_species = subset(percent_source, 
                                 Percent.data == "%_Species")

percent_source_species$Source5 = "F4-Aged Sea Salt"
percent_source_species$Source4 = "F2-Secondary Nitrate"
percent_source_species$Source1 = "F3-Secondary Sulfate"
percent_source_species$Source2 = "F8-Biomass"
percent_source_species$Source3 = "F9-Soil/Dust"
percent_source_species$Source7 = "F6-Fresh Sea Salt"

percent_source_determin = rbind(percent_source_var_result, percent_source_species)
write.csv(percent_source_determin, "percent_source_determin_C_6_F_9_noCsub.csv")

##### 5.2 source assignment - %_of_Species - noCsub noExtreme######
percent_source_species = subset(percent_source, 
                                Percent.data == "%_Species")

percent_source_species$Source2 = "F4-Aged Sea Salt"
percent_source_species$Source6 = "F2-Secondary Nitrate"
percent_source_species$Source8 = "F3-Secondary Sulfate"
percent_source_species$Source7 = "F8-Biomass"
percent_source_species$Source1 = "F9-Soil/Dust"
percent_source_species$Source9 = "F6-Fresh Sea Salt"

percent_source_determin = rbind(percent_source_var_result, percent_source_species)
write.csv(percent_source_determin, "percent_source_determin_C_6_F_9_noCsub_noExtreme.csv")


##### 5.3 source assignment - %_of_Species - ALL data ALL variable ######
percent_source_species = subset(percent_source, 
                                Percent.data == "%_Species")

percent_source_species$Source2 = "F4-Aged Sea Salt"
percent_source_species$Source3 = "F2-Secondary Nitrate"
percent_source_species$Source6 = "F3-Secondary Sulfate"
percent_source_species$Source7 = "F8-Biomass"
percent_source_species$Source8 = "F9-Soil/Dust"
percent_source_species$Source9 = "F6-Fresh Sea Salt"

percent_source_determin = rbind(percent_source_var_result, percent_source_species)
write.csv(percent_source_determin, "percent_source_determin_C_6_F_9_Alldata_AllVaraible.csv")

##### 5.4 source assignment - %_of_Species - 2011-17 ALL data ALL variable ######
percent_source_species = subset(percent_source, 
                                Percent.data == "%_Species")

percent_source_species$Source4 = "F4-Aged Sea Salt"
percent_source_species$Source8 = "F2-Secondary Nitrate"
percent_source_species$Source2 = "F3-Secondary Sulfate"
percent_source_species$Source3 = "F8-Biomass"
percent_source_species$Source6 = "F9-Soil/Dust"
percent_source_species$Source1 = "F6-Fresh Sea Salt??"

percent_source_determin = rbind(percent_source_var_result, percent_source_species)
write.csv(percent_source_determin, "percent_source_determin_C_6_F_9_Alldata_AllVaraible_2011-17.csv")

#############################################################################
###########  Prepare data generated from PMF GUI ########### 
#############################################################################
library(readxl)

site_ts = read_excel("/Users/ztttttt/Downloads/Hopper_VM_PMF_compare_2023.06.09.xlsx", sheet = "MatrixA")
site_conc = read_excel("/Users/ztttttt/Downloads/Hopper_VM_PMF_compare_2023.06.09.xlsx", sheet = "MatrixB")

mutate(site_ts, across(all_of(paste0("Factor", 1:7), as.numeric)))

#### Factor - concentration contributions ####
# read dataset of concentration contributions for each factor
site_conc = read.csv("C6S3F8_base_Concentration.csv")
colnames(site_conc)[1] = "Species"




# replace species names
site_conc$Species[site_conc$Species == "K."] = "KIon"
site_conc$Species[site_conc$Species == "Na."] = "NaIon"
site_conc$Species[site_conc$Species == "NH4."] = "NH4Ion"
site_conc$Species[site_conc$Species == "NO3"] = "NO3Ion"
site_conc$Species[site_conc$Species == "SO4"] = "SO4Ion"

# gather the data for plotting
site_conc = gather(site_conc, 
                   "Factor", 
                   "Concentration", 
                   -Species)

# match species class
site_conc = merge(site_conc, 
                  species_class, 
                  all.x = T)

# site_conc$data = "Concentration"

#### Factor - Percent contributions ####
# read dataset of percent contributions for each factor
site_perc = read.csv("C6S3F8_base_Percent.csv")
colnames(site_perc)[1] = "Species"

# replace species names
site_perc$Species[site_perc$Species == "K."] = "KIon"
site_perc$Species[site_perc$Species == "Na."] = "NaIon"
site_perc$Species[site_perc$Species == "NH4."] = "NH4Ion"
site_perc$Species[site_perc$Species == "NO3"] = "NO3Ion"
site_perc$Species[site_perc$Species == "SO4"] = "SO4Ion"

# gather the data for plotting
site_perc = gather(site_perc, 
                   "Factor", 
                   "Percent", 
                   -Species)

# match species class
site_perc = merge(site_perc, 
                  species_class, 
                  all.x = T)

# site_perc$data = "Percent"

#### Factor - time-series ####
# read dataset of concentration contributions for each factor
site_ts = read.csv("C6S3F8_base_TimeSeries.csv")
colnames(site_ts)[1] = "Date"
site_ts$Date = as.Date(site_ts$Date, format = "%m/%d/%y")

# match with PM2.5 concentration
site_ts_PM = site_ts
site_PM = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/Top_PM_cluster&site/CSN_C6_S3_2011-14_conc.csv")
site_ts_PM$PM2.5 = site_PM$PM25
site_ts_PM$Date = NULL

# gather the data for plotting
site_ts = gather(site_ts, 
                 "Factor", 
                 "Contribution", 
                 -Date)
summary(site_ts)

# factor contribution - linear
site_ts_PM_lm = lm(PM2.5 ~ ., data = site_ts_PM)
site_ts_PM_lm_beta = summary(site_ts_PM_lm)$coefficients[, 1]
site_ts_PM_lm_beta = data.frame(site_ts_PM_lm_beta[2:length(site_ts_PM_lm_beta)])
colnames(site_ts_PM_lm_beta) = "lm.beta.site-SiteBased"
site_ts_PM_lm_beta$Factor = rownames(site_ts_PM_lm_beta)
site_ts_PM_lm_beta$Factor.contribution = (site_ts_PM_lm_beta$lm.beta.site/
                                                   sum(site_ts_PM_lm_beta$lm.beta.site))*100

#### plotting - Concentration & percent contribution ####

site_conc_perc = merge(site_conc, site_perc)
sapply(site_conc_perc, class)
head(site_conc_perc)
site_conc_perc = subset(site_conc_perc,
                        Species != "PM25")


ggplot(site_conc_perc, 
       aes(x = reorder(Species, sequence), 
           group = Factor)) +
  facet_grid(Factor ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  scale_y_continuous(
    position = "right",
    name = "% of Species"
    # Features of the first axis
    #name = "Concentration of Species",
    # Add a second axis and specify its features
    #sec.axis = sec_axis(~., name="% of Species") # second y-axis, dual y-axis
    ) +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

ggplot(site_conc_perc, 
       aes(x = reorder(Species, sequence), 
           group = Factor)) +
  facet_grid(Factor ~.) +
  geom_point( aes(y = Percent), color = "black", shape = 15) +
  scale_y_continuous(
    position = "right",
    name = "% of Species"
  ) +
  theme_classic() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.line.x.bottom=element_line(color=NA),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

# for bar start from a value less than 0
# https://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
require(scales)
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

ggplot(site_conc_perc, 
       aes(x = reorder(Species, sequence), 
           group = Factor, fill = Factor)) +
  facet_grid(Factor ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6)+
  scale_y_log10() + 
  scale_y_continuous(trans = mylog_trans(base=10, from=-5)) +
  scale_fill_nejm() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))


ggplot(site_conc_perc, 
       aes(x = reorder(Species, sequence), 
           group = Factor, fill = Factor)) +
  facet_grid(Factor ~.) +
  geom_bar( aes(y = Concentration), 
            stat="identity", width = 0.6)+
  geom_point( aes(y = value.fake), color = NA) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


#### plotting - time series ####

ggplot(site_ts, 
       aes(x = Date, y = Contribution, 
           group = Factor, color = Factor)) +
  facet_grid(Factor ~.) +
  geom_line(alpha = 0.8)+
  geom_point(size = 0.5) +
  scale_color_nejm() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.05, vjust = -25, size = 16),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.x = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

#### plotting - factor percent contribution to PM2.5 ####
# stripe, https://coolbutuseless.github.io/package/ggpattern/articles/pattern-stripe.html

library(ggpattern)

ggplot(site_ts_PM_lm_beta, 
       aes(x = Factor, y = Factor.contribution)) +
 # geom_bar(stat="identity", width = 0.2) +    
  geom_col_pattern(aes(pattern_fill = Factor),  # fill = Factor
                   width = 0.2, fill = "white",
                   pattern = 'stripe', pattern_density = 0.5) +
  scale_fill_nejm() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5))



ggplot(site_ts_PM_lm_beta, 
       aes(x = Factor, y = Factor.contribution, fill = Factor)) +
  geom_bar(stat="identity", width = 0.2, pattern_fill = "slashes") +    
  scale_fill_nejm() +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 11, angle = 90, hjust = 0.5))

######## Below MDL (BDL) values. ########
# PM2.5 species classification, Strong, Weak, or Bad
# HEI_PM_Source_2nd_Committee_Meeting_20230606_Ting.pptx
# randomly extrac Cr concentration from csn_C6 to explain how BDL criteria work
bdl_Cr = read.csv("/Users/ztttttt/Downloads/Cr_fakeValue_for_BDL_criteria_2023.06.06.csv")
head(bdl_Cr)

ggplot(bdl_Cr, 
       aes(x = Serial.No, y = Cr)) +
  geom_point(size = 3, alpha = 0.8) +    
  facet_grid(. ~ Group) +
  geom_hline(yintercept = 0.00107*2, color = "red")+
  ylab(format_variable("Cr Concentration µg/m3")) +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour="grey60", linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 16),
        axis.text.x = element_text(color="grey25", size = 14, angle = 0, hjust = 0.5, vjust = 0), 
        axis.text.y = element_text(color="grey25", size = 14, angle = 90, hjust = 0.5),
        axis.title.x = element_text(color="grey25", size = 16, angle = 0, hjust = 0.5, vjust = 0),
        axis.title.y = element_text(color="grey25", size = 16, angle = 90, hjust = 0.5))


