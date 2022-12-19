##clear environment
# rm(list=ls())

##set working directory
setwd("XXXX")
getwd()
data.dir <- "XXXX"

##packages in need
require(tidyr) 
require(stats) 
require(ggplot2)
require(scales) 
require(stringr) 
require(dplyr)
require(plyr)
require(lubridate)
require(gridExtra) 
require(grid) 
library(data.table)

#### import & prepare data to use ####
imp_data = read.csv("IMPROVE component only 10092022.csv")
imp_data$X = imp_data$X.1 = NULL
head(imp_data)
setDT(imp_data)

imp_data$Date = as.Date(imp_data$Date)
# imp_data$CompName[is.na(imp_data$CompName)] = "Na"

imp_data$Qualifier = imp_data$Status
imp_data_compare = select(imp_data, Dataset, State, SiteCode, Date, CompName, Val, Qualifier)

csn_data = read.csv("CSN data for analysis 12122022.csv") ## with extracted collection, analysis methods
csn_data$X = NULL

head(csn_data)
setDT(csn_data)

csn_data$Date = as.Date(csn_data$Date)
csn_data = subset(csn_data, Date > as.Date("2010-12-31"))

csn_data$Dataset = "EPACSN"
excluded.variables.csn = c("Soil", "CS2", 
                           "Levoglucosan", "Mannosan", "Galactosan")
csn_data = subset(csn_data, !(CompName %in% excluded.variables.csn))
csn_data = subset(csn_data, Unit == "ug/m^3")
csn_data.1 = csn_data

# csn_ocec = subset(csn_data, grepl("OC1", CompName, fixed = T) | grepl("EC1", CompName, fixed = T))

csn_data$class[csn_data$CompName %in% c("Accept.PM2.5", "PM2.5RC")] = 0
csn_data$class[csn_data$CompName %in% c("Cl-", "K+", "NH4+", "Na+", "NO3", "SO4")] = "Ion"
csn_data.2 = csn_data

csn_data$Status = "V0"
csn_data$Qualifier = paste(csn_data$Qualifier1, csn_data$Qualifier2, csn_data$Qualifier3)
csn_data_compare = select(csn_data, Dataset, State, SiteCode, Date, CompName, Val, Qualifier)

# exclude sites according to previous summary 
# these sites a)not in mainland US, b)lack many species data, or c)short sampling duration
csn.exclude.site = c(800020014, 61072003, 20900034, 20904101, 
                     560210100, 720210010, 150030010) # site from CSN = 156-7 = 150
csn_data_compare = subset(csn_data_compare, 
                          !(SiteCode %in% csn.exclude.site))
dim(csn_data_compare)

#### 1. IMPROVE data preparation ####
head(imp_data_compare)
#imp_data_try = select(imp_data_compare, Dataset, State, SiteCode, Date, Qualifier, CompName, Val)[1:192, ]
#nrow(imp_data_try)/length(unique(imp_data_try$CompName))
#imp_data_try_sp = imp_data_try %>% spread(CompName, Val)
#head(imp_data_try_sp)

imp_data_use = select(imp_data_compare, Dataset, State, SiteCode, Date, Qualifier, CompName, Val)
nrow(imp_data_use)/length(unique(imp_data_use$CompName))

imp_state_site = select(imp_data_use, State, SiteCode)
imp_state_site = imp_state_site[!duplicated(imp_state_site), ]
imp_state_site$dup.site = duplicated(imp_state_site$SiteCode)
summary(imp_state_site$dup.site)

# create a data.frame with full dates and Components
imp_date = unique(imp_data_use$Date)
imp_Comp = unique(imp_data_use$CompName)
imp_site = unique(imp_data_use$SiteCode)

imp_date_comp_site = data.frame(Date = rep(imp_date, each = length(imp_Comp) * length(imp_site)), 
                                CompName = rep(imp_Comp, each = length(imp_site)), 
                                SiteCode = rep(imp_site))

# combine the full date-site-component list with the concentration values
imp_date_comp_site_data = merge(imp_date_comp_site, imp_data_use, all.x = T)

##### detect the duplicated values from POC!!! Parameter Occurrence Code  ####
# detect where the extracted duplicated date-site-component groups are from
imp_date_comp_site_data$dup.date.site.comp = duplicated(select(imp_date_comp_site_data, SiteCode, Date, CompName))
imp_date_comp_site_data_dup = subset(imp_date_comp_site_data, dup.date.site.comp)
length(unique(imp_date_comp_site_data_dup$Date)); length(imp_date)
length(unique(imp_date_comp_site_data_dup$SiteCode))
length(unique(imp_date_comp_site_data_dup$CompName))
# imp_date_comp_site_data_dup = imp_date_comp_site_data_dup[with(imp_date_comp_site_data_dup, order(State, SiteCode, Date, CompName)), ]
colnames(imp_date_comp_site_data_dup)[6:8] = c("Qualifier.dup", "Val.dup", "dup")

imp_date_comp_site_data_NOdup = subset(imp_date_comp_site_data, dup.date.site.comp == FALSE)
colnames(imp_date_comp_site_data_NOdup)[6:8] = c("Qualifier.NONEdup", "Val.NONEdup", "NONEdup")

imp_dup_date_comp_site = merge(imp_date_comp_site_data_dup, imp_date_comp_site_data_NOdup, all.x = T)
imp_dup_date_comp_site$oneValMissing = ifelse(imp_dup_date_comp_site$Val.dup == -999 | 
                                                imp_dup_date_comp_site$Val.NONEdup == -999,
                                              1, 0)
sum(imp_dup_date_comp_site$oneValMissing); nrow(imp_dup_date_comp_site)
imp_dup_date_comp_site_noMissing = subset(imp_dup_date_comp_site, oneValMissing == 0)
length(unique(imp_dup_date_comp_site_noMissing$Date)); length(imp_date)
length(unique(imp_dup_date_comp_site_noMissing$SiteCode))
length(unique(imp_dup_date_comp_site_noMissing$CompName))
imp_dup_date_comp_site_noMissing$dup = imp_dup_date_comp_site_noMissing$NONEdup = imp_dup_date_comp_site_noMissing$oneValMissing = NULL
head(imp_dup_date_comp_site_noMissing)
nrow(imp_dup_date_comp_site_noMissing)

##### check if the strange duplicated component values from the same date-site-component groups exist in original dataset #####
imp_data_org = read.csv("/Users/ztttttt/Documents/HEI PMF/IMPROVE & CSN original/ailsa2be_20221008_205315_uNQ0v IMPROVE.txt",
                        sep = ",", dec = ".")
head(imp_data_org)
dim(imp_data_org) # 18003337, 25
imp_meta_sites_use = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE/IMPROVE metadata 196 sample sites info 2010-20.csv")
imp_data_org = subset(imp_data_org, SiteCode %in% as.character(imp_meta_sites_use$SiteCode)) # some sites have no gps/location information
nrow(imp_data_org) # 17730878 
imp_data_org = subset(imp_data_org, Unit == "ug/m^3")
dim(imp_data_org) # 12659632, 25
as.Date(imp_data_org$Date[1], format = "%m/%d/%Y")
imp_data_org$Date = as.Date(imp_data_org$Date, format = "%m/%d/%Y")

imp_data_org_check_1 = subset(imp_data_org, Date == imp_dup_date_comp_site_noMissing$Date[1] &
                                SiteCode == imp_dup_date_comp_site_noMissing$SiteCode[1] &
                                ParamCode == "ALf")
imp_data_org_check_1

imp_data_org_check_2 = subset(imp_data_org, Date == imp_dup_date_comp_site_noMissing$Date[10000] &
                                SiteCode == imp_dup_date_comp_site_noMissing$SiteCode[10000] &
                                ParamCode == "ALf")
imp_data_org_check_2

subset(imp_data, Date == imp_dup_date_comp_site_noMissing$Date[10000] &
         SiteCode == imp_dup_date_comp_site_noMissing$SiteCode[10000] &
         ParamCode == "ALf")

##### the duplicated data are due to POC, Parameter Occurrence Code !!!!

##### compare concentrations between POC from same date-site-component groups ####
imp_poc = imp_dup_date_comp_site
imp_poc_noMissing = imp_dup_date_comp_site_noMissing
colnames(imp_poc_noMissing)[7] = "Val.POC.1"
colnames(imp_poc_noMissing)[9] = "Val.POC.2"

theme.poc = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                  strip.text.x = element_text(size = 12, colour = "grey25", angle = 0),
                  axis.title.x = element_text(color="grey25", size = 13, vjust=0, margin=margin(0,0,0,300)), 
                  axis.title.y = element_text(color="grey25", size = 13, vjust=1, margin=margin(0,2,0,0)),
                  axis.text.x = element_text(color="grey25", size = 12, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                  axis.text.y = element_text(color="grey25", size = 12, angle = 0, hjust = 0.5))

ggplot(subset(imp_poc_noMissing, grepl("OC", CompName, fixed = T) | grepl("EC", CompName, fixed = T)), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(imp_poc_noMissing, grepl("OP", CompName, fixed = T)), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(imp_poc_noMissing, CompName %in% c("Al", "Ca", "Cl",  "Fe", "K",  "Mg", "Na", "S", "Si")), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(imp_poc_noMissing, CompName %in% c("As", "Br", "Cr", "Cu", "Mn", "Ni", "P",  "Pb", "Rb","Se", 
                                                 "Sr", "Ti", "V",  "Zn", "Zr")), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(imp_poc_noMissing, CompName %in% c("Cl-", "NO3",  "SO4")), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

imp_poc_noMissing_P_Zn = subset(imp_poc_noMissing, CompName == "P" | CompName == "Zn")
unique(imp_poc_noMissing_P_Zn$Qualifier.dup); unique(imp_poc_noMissing_P_Zn$Qualifier.NONEdup)
unique(imp_poc_noMissing$Qualifier.dup); unique(imp_poc_noMissing$Qualifier.NONEdup)
subset(imp_poc_noMissing_P_Zn, Val.POC.2 < 0.01 & Val.POC.1 > 0.01)

ggplot(subset(imp_poc_noMissing, CompName == "P" | CompName == "Zn"), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc
ggplot(subset(imp_poc_noMissing, CompName == "S" | CompName == "SO4"), 
       aes(Val.POC.1, Val.POC.2, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc
# 
##### calculate component concentrations for those with POC = 2 ####
imp_data_use$POC = imp_data$POC
head(imp_data_use)
head(imp_dup_date_comp_site)
sapply(imp_dup_date_comp_site, class)
imp_dup_date_comp_site$Val.opc = (imp_dup_date_comp_site$Val.dup + imp_dup_date_comp_site$Val.NONEdup)/2
imp_dup_date_comp_site$Val.opc = ifelse(imp_dup_date_comp_site$Val.dup == -999, imp_dup_date_comp_site$Val.NONEdup, imp_dup_date_comp_site$Val.opc)
imp_dup_date_comp_site$Val.opc = ifelse(imp_dup_date_comp_site$Val.NONEdup == -999, imp_dup_date_comp_site$Val.dup, imp_dup_date_comp_site$Val.opc)
imp_dup_date_comp_site_conc = select(imp_dup_date_comp_site, Date, SiteCode, CompName, Val.opc)
imp_data_poc1 = join(imp_data_use, imp_dup_date_comp_site_conc)
# imp_data_poc1.1 = imp_data_poc1
imp_data_poc1$Val = ifelse(is.na(imp_data_poc1$Val.opc), imp_data_poc1$Val, imp_data_poc1$Val.opc)
imp_data_poc1 = subset(imp_data_poc1, POC == 1)
imp_poc_noNA = subset(imp_data_poc1, Val != -999)

imp_date_comp_site_noDup_noNA = merge(imp_date_comp_site, imp_poc_noNA, all.x = T)
nrow(imp_date_comp_site_noDup_noNA)
imp_date_comp_site_noDup_noNA$Val.opc = imp_date_comp_site_noDup_noNA$POC = 
  imp_date_comp_site_noDup_noNA$Dataset = imp_date_comp_site_noDup_noNA$State = NULL

imp_date_comp_site_noDup_noNA$site.date.qualifier = paste(imp_date_comp_site_noDup_noNA$Date,
                                                          imp_date_comp_site_noDup_noNA$SiteCode,
                                                          imp_date_comp_site_noDup_noNA$Qualifier)
imp_date_comp_site_noDup_noNA$Date = imp_date_comp_site_noDup_noNA$SiteCode = imp_date_comp_site_noDup_noNA$Qualifier = NULL
head(imp_date_comp_site_noDup_noNA)
write.csv(imp_date_comp_site_noDup_noNA, "IMPROVE_Component_to_be_spread.csv")

imp_date_comp_site_noDup_noNA = read.csv("IMPROVE_Component_to_be_spread.csv")
imp_daily_comp = imp_date_comp_site_noDup_noNA %>% spread(CompName, Val)
head(imp_daily_comp)
nrow(imp_daily_comp)
imp_site_date_qualifier = data.frame(select(imp_daily_comp, site.date.qualifier))
imp_site_date_qualifier_sep = imp_site_date_qualifier %>% 
  separate(site.date.qualifier, c("Date", "SiteCode", "Qualifier"),  
           sep = "\\s+")
head(imp_site_date_qualifier_sep)
imp_daily_comp_use = cbind(imp_site_date_qualifier_sep, imp_daily_comp)
imp_daily_comp_use_Q.Na = subset(imp_daily_comp_use, Qualifier == "NA")
summary(imp_daily_comp_use_Q.Na)
imp_daily_comp_use = subset(imp_daily_comp_use, Qualifier != "NA")

imp_meta_sites = read.csv("IMPROVE metadata 196 sample sites info 2010-20.csv")
imp_state_site = select(imp_meta_sites, State, SiteCode)

imp_daily_comp_use = merge(imp_daily_comp_use, imp_state_site, all.x = T)
imp_daily_comp_use$site.date.qualifier = imp_daily_comp_use$State
imp_daily_comp_use$State = NULL
colnames(imp_daily_comp_use)[4] = "State"
write.csv(imp_daily_comp_use, "IMPROVE_Component_with_missing.csv")

#### 2. CSN data preparation ####

##### 2.1 prepare dataset with duplicates and for reference #####
head(csn_data_compare)

# dataset with measurement duplicates
csn_data_use = select(csn_data_compare, Dataset, State, SiteCode, Date, Qualifier, CompName, Val)
csn_data_use = subset(csn_data_use, Date > as.Date("2010-12-31"))
# nrow(csn_data_use)/length(unique(csn_data_use$CompName))
dim(csn_data_use)
head(csn_data_use)

# csn_state_site = select(csn_data_use, State, SiteCode)
# csn_state_site = csn_state_site[!duplicated(csn_state_site), ]
# csn_state_site$dup.site = duplicated(csn_state_site$SiteCode)
# summary(csn_state_site$dup.site)

# create a data.frame with all dates and Components for each site
csn_date = unique(csn_data_use$Date)
csn_Comp = unique(csn_data_use$CompName)
csn_site = unique(csn_data_use$SiteCode)

# check & clear the messy dates
csn_date_dt = setDT(data.frame(csn_date))
colnames(csn_date_dt) = "Date"

# check the date different between one row and the row below - every-3-day sample?
csn_date_dt = csn_date_dt %>%
  mutate(diff = Date - lag(Date, default = first(Date)))

# after observation, separate csn_date_dt into
csn_date_normalSamp = csn_date_dt[1:1143, ]
# write.csv(csn_date_normalSamp, "CSN_cira_normal_sampling_dates.csv")
csn_date_oddSamp = csn_date_dt[1144:nrow(csn_date_dt), ]
csn_oddDateSamp = subset(csn_data_use, 
                         Date %in% csn_date_oddSamp$Date)
csn_normalDateSamp = subset(csn_data_use, 
                            !(Date %in% csn_date_oddSamp$Date))

unique(csn_oddDateSamp$SiteCode)
summary(unique(csn_oddDateSamp$SiteCode) %in% 
          csn_normalDateSamp$SiteCode)
# all sites with odd sampling dates have the normal sampling dates
# then just remove the the odd date data

# ONLY include the normal dates
csn_date = csn_date[1:1143]

csn_date_comp_site_ref = data.frame(Date = rep(csn_date, each = length(csn_Comp) * length(csn_site)), 
                                    CompName = rep(csn_Comp, each = length(csn_site)), 
                                    SiteCode = rep(csn_site))
dim(csn_date_comp_site_ref)
head(csn_date_comp_site_ref)
unique(csn_date_comp_site_ref$Date)[1:15]

##### 2.2 check POC repeatability (duplicated measurements) #####
csn_data$Method = csn_data$Status = csn_data$Elevation = 
  csn_data$Qualifier1 = csn_data$Qualifier2 = csn_data$Qualifier3 = 
  csn_data$Unit = csn_data$Unc = csn_data$MDL = csn_data$SampleDuration = 
  csn_data$year = csn_data$month = csn_data$day = NULL
csn_data$date.site.comp = paste0(csn_data$Date, csn_data$SiteCode, csn_data$CompName)

# check the type and number of different POC in CSN
table(csn_data$POC) # POC = 5, 6, or 7

csn_data_poc5 = subset(csn_data, POC == 5)
csn_data_poc6 = subset(csn_data, POC == 6)
csn_data_poc7 = subset(csn_data, POC == 7)

# check if there could be >2 measurements for one species from the same date & site
summary(csn_data_poc7$date.site.comp %in% csn_data_poc6$date.site.comp)
summary(csn_data_poc6$date.site.comp %in% csn_data_poc5$date.site.comp)
summary(csn_data_poc7$date.site.comp %in% csn_data_poc5$date.site.comp)
# the results suggest there are up to 2 measurements for a given date/site/specie group

# prepare data for comparison between duplicated measurement 
csn_data_poc5_dup_use = select(csn_data_poc5, date.site.comp, Val,  
                               Collection, Analysis, Qualifier)
csn_data_poc5_dup_use$poc5 = "Y"
csn_data_poc7_dup_use = select(csn_data_poc7, date.site.comp, Val, 
                               Collection, Analysis, Qualifier)
csn_data_poc7_dup_use$poc7 = "Y"

# duplication type 1, POC 5 & 6
csn_poc56_dup = merge(csn_data_poc6, csn_data_poc5_dup_use, by = "date.site.comp", all.x = T)
csn_poc56_dup = subset(csn_poc56_dup, poc5 == "Y")
dim(csn_poc56_dup)
# check if collection methods and qualifiers for the duplicates are same
summary(csn_poc56_dup$Collection.x == csn_poc56_dup$Collection.y)
summary(csn_poc56_dup$Analysis.x == csn_poc56_dup$Analysis.y)
summary(csn_poc56_dup$Qualifier1.x == csn_poc56_dup$Qualifier1.y)
summary(csn_poc56_dup$Qualifier2.x == csn_poc56_dup$Qualifier2.y)
summary(csn_poc56_dup$Qualifier3.x == csn_poc56_dup$Qualifier3.y)

# duplication type 2, POC 6 & 7
csn_poc67_dup = merge(csn_data_poc6, csn_data_poc7_dup_use, by = "date.site.comp", all.x = T)
csn_poc67_dup = subset(csn_poc67_dup, poc7 == "Y")
dim(csn_poc67_dup)
# check if collection methods and qualifiers for the duplicates are same
summary(csn_poc67_dup$Collection.x == csn_poc67_dup$Collection.y)
summary(csn_poc67_dup$Analysis.x == csn_poc67_dup$Analysis.y)
summary(csn_poc67_dup$Qualifier1.x == csn_poc67_dup$Qualifier1.y)
summary(csn_poc67_dup$Qualifier2.x == csn_poc67_dup$Qualifier2.y)
summary(csn_poc67_dup$Qualifier3.x == csn_poc67_dup$Qualifier3.y)

# prepare the duplicate data for merging   
csn_poc56_dup$Val.dup = (csn_poc56_dup$Val.x + csn_poc56_dup$Val.y)/2 
# attention, above function, divided by 2 is different from using "mean"
dup_info56 = select(csn_poc56_dup, 
                    date.site.comp, POC, 
                    Date, SiteCode, 
                    CompName, ParamCode, Val.dup)
dup_info56$POC = 56

csn_poc67_dup$Val.dup = (csn_poc67_dup$Val.x + csn_poc67_dup$Val.y)/2
dup_info67 = select(csn_poc67_dup,  
                    date.site.comp, POC, 
                    Date, SiteCode, 
                    CompName, ParamCode, Val.dup)
dup_info67$POC = 67

# double check if there are duplicates 
dup_info_56_67 = rbind(dup_info56, dup_info67)
dup_info_56_67$dup = duplicated(dup_info_56_67$date.site.comp)
summary(dup_info_56_67$dup)
dup_info_56_67 = select(dup_info_56_67, 
                        date.site.comp, 
                        Val.dup,
                        POC)
head(dup_info_56_67)
write.csv(dup_info_56_67, "CSN_Duplicates_POC.csv")

# check how the duplicates match 
ggplot(subset(csn_poc56_dup, grepl("OC", CompName, fixed = T) & 
                !grepl("OC.", CompName, fixed = T) & CompName != "OC"), 
       aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, grepl("EC", CompName, fixed = T) & 
                !grepl("EC.", CompName, fixed = T) & CompName != "EC"), 
       aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, CompName %in% c("SO4", "NO3", "Cl-", "NH4+", "Na+", "K+")), 
       aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, CompName %in% unique(csn_poc56_dup$CompName[
  csn_poc56_dup$class == "element.ion.PM"])[1:9]), 
  aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, CompName %in% unique(csn_poc56_dup$CompName[
  csn_poc56_dup$class == "element.ion.PM"])[14:22]), 
  aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, CompName %in% unique(csn_poc56_dup$CompName[
  csn_poc56_dup$class == "element.ion.PM"])[23:32] & CompName != "NO3"), 
  aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

ggplot(subset(csn_poc56_dup, CompName %in% unique(csn_poc56_dup$CompName[
  csn_poc56_dup$class == "element.ion.PM"])[34:38]), 
  aes(Val.x, Val.y, color = CompName)) +
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.poc

##### 2.3 mark the duplicated measurement and prepare dataset for interpolation #####
csn_date_comp_site_measured = csn_data_use

csn_date_comp_site_measured$date.site.comp = paste0(csn_date_comp_site_measured$Date, 
                                                    csn_date_comp_site_measured$SiteCode, 
                                                    csn_date_comp_site_measured$CompName)
dim(csn_date_comp_site_measured)

# merge the information of whether it is a duplicated measurement into dataset
csn_date_comp_with_dup = merge(csn_date_comp_site_measured, 
                               dup_info_56_67)
dim(csn_date_comp_with_dup)
csn_date_comp_with_dup$Val = ifelse(is.na(csn_date_comp_with_dup$Val.dup), 
                                    csn_date_comp_with_dup$Val,
                                    csn_date_comp_with_dup$Val.dup)
# each line has one duplicate, remove the duplicates via using row number
csn_date_comp_with_dup$row.No = 1:nrow(csn_date_comp_with_dup)
csn_date_comp_with_dup = subset(csn_date_comp_with_dup, 
                                row.No %% 2 == 0)
csn_date_comp_with_dup$row.No = NULL



# extract the dataset without duplicates
csn_date_comp_no_dup = subset(merge(csn_date_comp_site_measured, 
                                    dup_info_56_67, all.x = T), 
                              is.na(POC))

# to check if we have miss any data group from original csn_data
length(unique(csn_date_comp_no_dup$date.site.comp)) + # 6425441
  length(unique(csn_date_comp_with_dup$date.site.comp)) == # 179081
  length(unique(csn_date_comp_site_measured$date.site.comp)) # 6604522
nrow(csn_date_comp_site_measured) == # 6783603
  nrow(csn_date_comp_no_dup) + # 6425441
  2*nrow(csn_date_comp_with_dup) # 179081

# combine those with & without duplicated measurements 
csn_date_comp_site_data = rbind(csn_date_comp_with_dup, 
                                csn_date_comp_no_dup)
summary(csn_date_comp_site_data)

# only keep data sampled in normal dates
csn_date_comp_site_data = subset(csn_date_comp_site_data,
                                 Date %in% csn_date)

##### 2.4 merge data for reference (include all date-site-component groups) & Val #####

# the detected non-negative Val is 5e-06 , set all non-positive Val to 0.000009
# the 0.000009 will be transferred to MDL/2 for each component later
csn_date_comp_site_data$Val[csn_date_comp_site_data$Val <= 5e-06] = 0.000009

# prepare information before spreading the dataset
csn_date_comp_site_data$site.date.qualifier = paste(csn_date_comp_site_data$Date,
                                                    csn_date_comp_site_data$State,
                                                    csn_date_comp_site_data$SiteCode,
                                                    csn_date_comp_site_data$Qualifier)

# delete columns not needed for merging or not to be used in the spread dataset
csn_date_comp_site_data$Qualifier = 
  csn_date_comp_site_data$Val.dup = csn_date_comp_site_data$POC = 
  csn_date_comp_site_data$Dataset = csn_date_comp_site_data$State = NULL

# double check if the date, site, and component match in the reference dataset and the extracted data
summary(unique(csn_date_comp_site_ref$Date) %in% unique(csn_date_comp_site_data$Date))
summary(unique(csn_date_comp_site_ref$SiteCode) %in% unique(csn_date_comp_site_data$SiteCode))
summary(unique(csn_date_comp_site_ref$CompName) %in% unique(csn_date_comp_site_data$CompName))

summary(unique(csn_date_comp_site_data$Date) %in% unique(csn_date_comp_site_ref$Date))
summary(unique(csn_date_comp_site_data$SiteCode) %in% unique(csn_date_comp_site_ref$SiteCode))
summary(unique(csn_date_comp_site_data$CompName) %in% unique(csn_date_comp_site_ref$CompName))

# combine the full date-site-component list with the concentration values
csn_date_with_missing = merge(csn_date_comp_site_ref, 
                              csn_date_comp_site_data, all.x = T)

length(unique(csn_date_comp_site_measured$date.site.comp)) ==
  length(unique(csn_date_with_missing$date.site.comp))
length(unique(csn_date_with_missing$site.date.qualifier))

csn_site_date_qualifier = data.frame(select(csn_date_with_missing, site.date.qualifier))
setDT(csn_site_date_qualifier)
csn_site_date_qualifier_sep = csn_site_date_qualifier %>% 
  separate(site.date.qualifier, c("Date", "State", "SiteCode", "Qualifier"),  
           sep = "\\s+")
csn_site_date_qualifier_sep_noNA = subset(csn_site_date_qualifier_sep,
                                          !is.na(SiteCode))
csn_site_date_qualifier_sep_noNA$Date = as.Date(csn_site_date_qualifier_sep_noNA$Date)
csn_site_date_qualifier_sep_noNA$SiteCode = as.integer(csn_site_date_qualifier_sep_noNA$SiteCode)
csn_site_date_qualifier_sep_noNA = 
  csn_site_date_qualifier_sep_noNA[!duplicated(csn_site_date_qualifier_sep_noNA), ]

# 
csn_date_miss  = merge(csn_date_with_missing, 
                       csn_site_date_qualifier_sep_noNA)
setDT(csn_date_miss)
summary(csn_date_miss)
summary(is.na(csn_date_miss$State))
summary(is.na(csn_date_miss$CompName))
summary(is.na(csn_date_miss$Qualifier))

csn_date_miss$site.date.qualifier = NULL
csn_date_miss$site.date.qualifier = paste(csn_date_miss$Date, csn_date_miss$State, 
                                          csn_date_miss$SiteCode, csn_date_miss$Qualifier)

# check if qualifiers for same day-site group are same
# csn_miss_sample = subset(csn_date_with_missing, 
#                          Date == csn_date_with_missing$Date[1] &
#                            SiteCode == csn_date_with_missing$SiteCode[1])
# length(unique(csn_miss_sample$site.date.qualifier))

# write.csv(csn_date_miss, "CSN_Component_to_be_spread_2022-12.csv")

##### 2.5 spreading - data for interpolation - before & after 2015 #####
library(psych) # corr.test{}
library(corrplot) # corrplot.mixed{}

# csn_miss_to_spread = read.csv("CSN_Component_to_be_spread_2022-12.csv")
# csn_miss_to_spread$X = NULL
csn_miss_to_spread = csn_date_miss

# spreading the data to daily-site-component style
csn_miss_to_spread$Date = csn_miss_to_spread$SiteCode = 
  csn_miss_to_spread$State = csn_miss_to_spread$Qualifier = NULL
head(csn_miss_to_spread)
dim(csn_miss_to_spread)

csn_spread_try = csn_miss_to_spread[1:888, ]
csn_spread_try %>% spread(CompName, Val)

csn_daily_comp = csn_miss_to_spread %>% spread(CompName, Val)
head(csn_daily_comp)
dim(csn_daily_comp)
csn_site_date_qualifier_sep_noNA$site.date.qualifier = paste(csn_site_date_qualifier_sep_noNA$Date, 
                                                             csn_site_date_qualifier_sep_noNA$State, 
                                                             csn_site_date_qualifier_sep_noNA$SiteCode, 
                                                             csn_site_date_qualifier_sep_noNA$Qualifier)
csn_daily_comp = merge(csn_daily_comp, csn_site_date_qualifier_sep_noNA, all.x = T)
summary(csn_daily_comp)
csn.d.col = ncol(csn_daily_comp)
csn_daily_comp_use = cbind(csn_daily_comp[, (csn.d.col-3):csn.d.col],
                           csn_daily_comp[, 2:(csn.d.col-4)])
write.csv(csn_daily_comp_use, "CSN_Component_with_missing.csv")

##### 2.6 check correlation & missing pattern #####
# checking the distributinon of OC.unadjusted.88 vs. OC.TOR.unadjusted.88
plot(csn_daily_comp_use$OC.unadjusted.88, 
     csn_daily_comp_use$OC.TOR.unadjusted.88)

cor(csn_daily_comp_use$OC.unadjusted.88, 
    csn_daily_comp_use$OC.TOR.unadjusted.88,
    use="complete.obs", 
    method = "pearson") # 0.9937977
cor(csn_daily_comp_use$OC.unadjusted.88, 
    csn_daily_comp_use$OC.TOR.unadjusted.88,
    use="complete.obs", 
    method = "spearman") # 0.994943

# checking the distributinon of OC vs. OC.88
plot(csn_daily_comp_use$OC, 
     csn_daily_comp_use$OC.88)

cor(csn_daily_comp_use$OC, 
    csn_daily_comp_use$OC.88,
    use="complete.obs", 
    method = "pearson") # 0.9942456
cor(csn_daily_comp_use$OC, 
    csn_daily_comp_use$OC.88,
    use="complete.obs", 
    method = "spearman") # 0.9949892

# checking the distributinon of OC.unadjusted.88 vs. OC.88
plot(csn_daily_comp_use$OC.unadjusted.88, 
     csn_daily_comp_use$OC.88)

cor(csn_daily_comp_use$OC.unadjusted.88, 
    csn_daily_comp_use$OC.88,
    use="complete.obs", 
    method = "pearson") # 0.9733577
cor(csn_daily_comp_use$OC.unadjusted.88, 
    csn_daily_comp_use$OC.88,
    use="complete.obs", 
    method = "spearman") # 0.9991176

# checking the overall correlations for OC/EC groups
csn_ocec = data.frame(csn_daily_comp_use)[, grepl("OC", names(data.frame(csn_daily_comp_use))) |
                                            grepl("EC", names(data.frame(csn_daily_comp_use)))]
csn_ocec_before_2015 = csn_ocec[1:62120, ]
csn_ocec_after_2015 = csn_ocec[62121:nrow(csn_ocec), ]

# convert the dataset to matrix for cor function - before 2015
before_ocec_m <- as.matrix(csn_ocec_before_2015[, 8:15])
after_ocec_m <- as.matrix(csn_ocec_after_2015[, 8:15])

csn_ocec_m = before_ocec_m
csn_ocec_m = after_ocec_m

# calculate pearson & spearman correlations - before 2015
corPsn_ocec = cor(csn_ocec_m, 
                  use="pairwise.complete.obs", 
                  method = "pearson")
corSpmn_ocec = cor(csn_ocec_m, 
                   use="pairwise.complete.obs", 
                   method = "spearman")

p_psn_ocec <- as.data.frame(corr.test(csn_ocec_m, 
                                      adjust = "none", 
                                      method = "pearson")[4])
p_spm_ocec <- as.data.frame(corr.test(csn_ocec_m, 
                                      adjust = "none", 
                                      method = "spearman")[4])

# only keep value with significant correlations (p < 0.05)
NA_psn <- p_psn_ocec > 0.05
NA_spm <- p_spm_ocec > 0.05
corPsn_ocec[NA_psn] <- 0
corSpmn_ocec[NA_spm] <- 0

# plot pearson & spearman correlations
corEF_Psn <- corrplot.mixed(as.matrix(corPsn_ocec)) # ,mar=c(0,0,0,0), title = "XX"
corEF_Spmn <- corrplot.mixed(as.matrix(corSpmn_ocec)) # ,mar=c(0,0,0,0), title = "XX"

# calculate the missing data rate
csn_daily_before_2015 = subset(csn_daily_comp_use, 
                               Date <= as.Date("2015-12-31"))
csn_daily_after_2015 = subset(csn_daily_comp_use, 
                              Date > as.Date("2015-12-31"))
dim(csn_daily_before_2015)
dim(csn_daily_after_2015)

# overall NA rate
p_miss_with_neg <- data.frame(unlist(lapply(csn_daily_comp[, 2:(csn.d.col-4)], 
                                            function(x) 
                                              sum(is.na(x))))/
                                nrow(csn_daily_comp))
p_miss_with_neg$CompCode = rownames(p_miss_with_neg) 
colnames(p_miss_with_neg)[1] = "Missing_Rate"
rownames(p_miss_with_neg) = 1:nrow(p_miss_with_neg)
ggplot(p_miss_with_neg, aes(CompCode, Missing_Rate)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

write.csv(p_miss_with_neg, "CSN_Overall_missing_rate_NA.csv")

# before 2015 - NA rate
p_miss_with_neg_before2015 <- data.frame(unlist(lapply(
  csn_daily_before_2015[, 5:ncol(csn_daily_before_2015)], 
  function(x) 
    sum(is.na(x))))/
    nrow(csn_daily_before_2015))

p_miss_with_neg_before2015$CompCode = rownames(p_miss_with_neg_before2015) 
colnames(p_miss_with_neg_before2015)[1] = "Missing_Rate"
rownames(p_miss_with_neg_before2015) = 1:nrow(p_miss_with_neg_before2015)
ggplot(p_miss_with_neg_before2015, aes(CompCode, Missing_Rate)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

write.csv(p_miss_with_neg_before2015, "CSN_Overall_missing_rate_NA_before_2015.csv")

# after 2015 - NA rate
p_miss_with_neg_after2015 <- data.frame(unlist(lapply(
  csn_daily_after_2015[, 5:ncol(csn_daily_after_2015)], 
  function(x) 
    sum(is.na(x))))/
    nrow(csn_daily_after_2015))

p_miss_with_neg_after2015$CompCode = rownames(p_miss_with_neg_after2015) 
colnames(p_miss_with_neg_after2015)[1] = "Missing_Rate"
rownames(p_miss_with_neg_after2015) = 1:nrow(p_miss_with_neg_after2015)
ggplot(p_miss_with_neg_after2015, aes(CompCode, Missing_Rate)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
        axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
        axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
        axis.text.x = element_text(color="grey25", size = 11, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
        axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

write.csv(p_miss_with_neg_after2015, "CSN_Overall_missing_rate_NA_after_2015.csv")

# using mice to check the missing data pattern
library(mice)

plot.window(xlim=c(-1, ncol(csn_daily_before_2015) + 1), 
            ylim=c(-1, nrow(csn_daily_before_2015) + length_of_longest_colname), asp=1)
mice::md.pattern(csn_daily_before_2015[, 5:ncol(csn_daily_before_2015)], rotate.names = T) 

plot.window(xlim=c(-1, ncol(csn_daily_after_2015) + 1), 
            ylim=c(-1, nrow(csn_daily_after_2015) + length_of_longest_colname), asp=1)
mice::md.pattern(csn_daily_after_2015[, 5:ncol(csn_daily_after_2015)], rotate.names = T) 

