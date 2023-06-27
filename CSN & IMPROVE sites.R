##clear environment
# rm(list=ls())

##set working directory
setwd("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE")
getwd()
data.dir <- "Users/ztttttt/Documents/HEI PMF/R - original IMPROVE"

##packages in need
require(tidyr) # separate{tidyr}, gather{tidyr}, spread{tidyr},  spread is VIP function, str_split_fixed{stringr} is better than separate
require(stats) # aggregate{stats}, VIP function
require(ggplot2)
require(scales) # percent{}
require(stringr) # str_split_fixed, separate one column into multiple
require(dplyr)
require(plyr)
require(lubridate)
require(gridExtra) #grid.arrange{}
require(grid) #textGrob{}
library(data.table)

library(maps) 
library(usmap)
library(ggrepel)

#### input IMPROVE data ####
imp_data = read.csv("IMPROVE component only 10092022.csv")
imp_data$X = imp_data$X.1 = NULL
head(imp_data)
imp_data$Date = as.Date(imp_data$Date)
setDT(imp_data)

# imp_data$CompName[is.na(imp_data$CompName)] = "Na"
imp_data$Method = imp_data$MethodSimple
imp_data$class = 0
imp_data$class[imp_data$Method == "A-XRF"] = "Element"
imp_data$class[imp_data$Method == "B-IC"] = "Ion"
imp_data$class[imp_data$Method == "C-TOR"] = "OC/EC subgroup"

imp_data$Qualifier = imp_data$Status
imp_data_compare = select(imp_data, Dataset, State, SiteCode, Date, Method, class, CompName, Val, Qualifier, 
                                 Unc, MDL, year, month)

#### input CSN data ####

# csn_data = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original CSN/CSN method collection analysis use 10032022.csv") ## with extracted collection, analysis methods
# csn_data$Date = as.Date(csn_data$Date, format="%m/%d/%Y")
# csn_data = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original CSN/CSN data for analysis 10232022.csv") ## with extracted collection, analysis methods
csn_data = read.csv("CSN data for analysis 12122022.csv") ## with extracted collection, analysis methods
csn_data$X = NULL

head(csn_data)
setDT(csn_data)

csn_data$Date = as.Date(csn_data$Date)
csn_data = subset(csn_data, Date > as.Date("2010-12-31"))

csn_data$Dataset = "EPACSN"
excluded.variables.csn = c("FlowRate", "FlowRate.tf", "FlowRate.ny", "FlowRate.qz", 
                           "Volumn.tf", "Volumn.qz", "avgT.URG", "avgP.URG", "Volumn.ny", 
                           "Soil", "Volume", "CS2", 
                           "Levoglucosan", "Mannosan", "Galactosan", 
                           "MinT", "MaxT", "avgT", "MinP", "MaxP", "avgP")
csn_data = subset(csn_data, !(CompName %in% excluded.variables.csn))
# csn_data.1 = csn_data

# csn_ocec = subset(csn_data, grepl("OC1", CompName, fixed = T) | grepl("EC1", CompName, fixed = T))

csn_data$class[csn_data$class == "element.ion.PM"] = "Element"
csn_data$class[csn_data$CompName %in% c("Accept.PM2.5", "PM2.5RC")] = 0
csn_data$class[csn_data$CompName %in% c("Cl-", "K+", "NH4+", "Na+", "NO3", "SO4")] = "Ion"
csn_data$class[csn_data$class == "OC.EC"] = "OC/EC subgroup"
# csn_data.2 = csn_data

csn_data$Method = csn_data$Analysis
csn_data$Status = "V0"
csn_data$Qualifier = paste(csn_data$Qualifier1, csn_data$Qualifier2, csn_data$Qualifier3)
csn_data_compare = select(csn_data, Dataset, State, SiteCode, Date, Method, class, CompName, Val, Qualifier, 
                                 Unc, MDL, year, month)

# exclude sites according to previous summary 
# these sites a)not in mainland US, b)lack many species data, or c)short sampling duration
csn.exclude.site = c(800020014, 61072003, 20900034, 20904101, 
                     560210100, 720210010, 150030010) # site from CSN = 156-7 = 150
csn_data_compare = subset(csn_data_compare, 
                          !(SiteCode %in% csn.exclude.site))
dim(csn_data_compare)

#### merge IMPROVE & CSN dataset ####

imp_csn_data = rbind(imp_data_compare, csn_data_compare)

imp_meta_sites = read.csv("IMPROVE metadata 196 sample sites info 2010-20.csv")
imp_meta_sites$StartDate = as.Date(imp_meta_sites$StartDate)
imp_meta_sites$EndDate = as.Date(imp_meta_sites$EndDate)
imp_sites_use = subset(imp_meta_sites, Longitude > -999)
head(imp_sites_use)
imp_sites_use = select(imp_sites_use, Dataset, State, SiteCode, Latitude, Longitude, StartDate, EndDate)

csn_meta_sites = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original CSN/CSN metadata sample sites 2010-20 use.csv")
csn_sites_use = subset(csn_meta_sites, 
                       !(SiteCode %in% 
                           c("800020014", "61072003", "20900034", 
                             "20904101", "560210100", "720210010",
                             "150030010")))
csn_sites_use$StartDate = as.Date(csn_sites_use$StartDate, format = "%m/%d/%Y")
csn_sites_use$EndDate = as.Date(csn_sites_use$EndDate, format = "%m/%d/%Y")
head(csn_sites_use)
csn_sites_use = select(csn_sites_use, Dataset, State, SiteCode, Latitude, Longitude, StartDate, EndDate)

imp_csn_site = rbind(imp_sites_use, csn_sites_use)

##########################################################################################
######## check the site distance between CSN & IMPROVE ########
##########################################################################################
library(FNN)

imp_gps = select(imp_sites_use, Latitude, Longitude)
csn_gps = select(csn_sites_use, Latitude, Longitude)

near_points = get.knnx(imp_gps, csn_gps, k=1)
near_csn_in_imp = imp_gps[near_points$nn.index[,1], ]
near_csn_in_imp$distance = near_points$nn.dist[,1]
near_csn_in_imp$csn.site = csn_sites_use$SiteCode
near_csn_in_imp$distance = round(near_csn_in_imp$distance, 3)
colnames(near_csn_in_imp)[1:2] = c("imp.Latitude", "imp.Longitude")
near_csn_in_imp$imp.site = imp_sites_use$SiteCode[near_points$nn.index[,1]]
near_csn_in_imp$csn.Latitude = csn_sites_use$Latitude
near_csn_in_imp$csn.Longitude = csn_sites_use$Longitude
near_csn_in_imp$State = csn_sites_use$State
near_csn_in_imp = near_csn_in_imp[with(near_csn_in_imp, order(distance)), ]
write.csv(near_csn_in_imp, "The closest sampling points of CSN in IMPROVE.csv")

MainStates <- map_data("state")
theme.3 = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                axis.title.x = element_text(color="grey25", size = 12, vjust=0, margin=margin(0,0,0,300)), 
                axis.title.y = element_text(color="grey25", size = 12, vjust=1, margin=margin(0,2,0,0)),
                axis.text.x = element_text(color="grey25", size = 11, angle = 0, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                axis.text.y = element_text(color="grey25", size = 11, angle = 0, hjust = 0.5))

# subset(imp_csn_site, Longitude < -140), sites from AK & HI
# subset(imp_csn_site, Latitude < 20), sites from VI (Idaho?) & HI
ggplot(subset(imp_csn_site, Longitude > -140 & Latitude > 20), 
       aes(Longitude, Latitude, color= Dataset)) + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="darkgrey", fill="lightgrey", alpha = 0.6) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("yellow", "purple")) +
  ggtitle("Distribution of IMPROVE & CSN sites") + 
  theme.3

ggplot(subset(imp_csn_site, Longitude > -140 & Latitude > 20 & Dataset == "IMPAER"), 
       aes(Longitude, Latitude, color= Dataset)) + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="darkgrey", fill="lightgrey", alpha = 0.6) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = "purple") +
  ggtitle("Distribution of IMPROVE sites") + 
  theme.3

ggplot(subset(imp_csn_site, Longitude > -140 & Latitude > 20 & Dataset == "EPACSN"), 
       aes(Longitude, Latitude, color= Dataset)) + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="darkgrey", fill="lightgrey", alpha = 0.6) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = "yellow") +
  ggtitle("Distribution of IMPROVE sites") + 
  theme.3

##########################################################################################
######## check component concentration relationship between CSN & IMPROVE for sites of identical GPS ########
##########################################################################################
#### extract data for the nearest sites ####
near_csn_in_imp = read.csv( "The closest sampling points of CSN in IMPROVE.csv"); near_csn_in_imp$X = NULL
nearest_csn_imp = near_csn_in_imp[1:7, ]
unique(nearest_csn_imp$csn.site)
unique(nearest_csn_imp$imp.site)

nearest_csn_imp_sites = select(nearest_csn_imp, csn.site, imp.site)
nearest_csn_imp_sites$SameSite = 1:nrow(nearest_csn_imp_sites)
# colnames(nearest_csn_imp_sites)[1] = "SiteCode"

# double check if it is the nearest site
i
unique(subset(csn_data, SiteCode == nearest_csn_imp_sites$csn.site[i])$Longitude)
unique(subset(imp_data, SiteCode == nearest_csn_imp_sites$imp.site[i])$Longitude)
unique(subset(csn_data, SiteCode == nearest_csn_imp_sites$csn.site[i])$Latitude)
unique(subset(imp_data, SiteCode == nearest_csn_imp_sites$imp.site[i])$Latitude)

nearest.csn.sites = as.character(unique(nearest_csn_imp$csn.site))
nearest.imp.sites  = unique(nearest_csn_imp$imp.site)
nearest.sites = c(nearest.csn.sites, nearest.imp.sites)

# imp_csn_nearest_site_data = subset(imp_csn_data, SiteCode %in% nearest.sites)
## output component info for the side-by-side sampling in CSN & IMPROVE
# write.csv(imp_csn_nearest_site_data, "IMPROVE_CSN_Nearest_site_comparison.csv") 

#### JUMP to the subtask - overall component distribution check between nearest sites - JUMP to the subtask of by date & measurement ####
imp_csn_nearest_site_data = read.csv("IMPROVE_CSN_Nearest_site_comparison.csv") 

imp_csn_nearest_site_data = subset(imp_csn_nearest_site_data, Val != -999)
unique(subset(imp_csn_nearest_site_data, Dataset == "IMPAER")$Date)[1:10]
unique(subset(imp_csn_nearest_site_data, Dataset == "EPACSN")$Date)[1:10]

## to check if the sampling date of nearest site is the same
imp_date = data.frame(table(subset(imp_csn_nearest_site_data, Dataset == "IMPAER")$Date))
csn_date = data.frame(table(subset(imp_csn_nearest_site_data, Dataset == "EPACSN")$Date))
colnames(imp_date) = c("Date", "Freq.IMP")
colnames(csn_date) = c("Date", "Freq.CSN")
imp_csn_date = join(csn_date, imp_date)
csn_extra_date_check = subset(imp_csn_date, is.na(Freq.IMP))
csn_extra_date = subset(select(imp_csn_nearest_site_data, Dataset, State, SiteCode, Date, CompName), 
                        as.character(Date) %in% as.character(csn_extra_date_check$Date))
# write.csv(csn_extra_date, "Extra_date_only_in_CSN_data_during_nearest_site_comparison.csv")

imp_csn_nearest_Val_data = select(imp_csn_nearest_site_data, Dataset, SiteCode, Date, CompName, Val, class, Method, Qualifier)
imp_nearest_data = subset(imp_csn_nearest_Val_data, Dataset == "IMPAER")
csn_nearest_data = subset(imp_csn_nearest_Val_data, Dataset == "EPACSN")

nearest_csn_imp_sites_csn = select(nearest_csn_imp_sites, csn.site, SameSite)
colnames(nearest_csn_imp_sites_csn)[1] = "SiteCode"
nearest_csn_imp_sites_imp = select(nearest_csn_imp_sites, imp.site, SameSite)
colnames(nearest_csn_imp_sites_imp)[1] = "SiteCode"

csn_nearest_data = join(csn_nearest_data, nearest_csn_imp_sites_csn)
csn_nearest_data$CompName[csn_nearest_data$CompName == "OC.88" ] = "OC"
csn_nearest_data$CompName[csn_nearest_data$CompName == "OC1.88" ] = "OC1"
csn_nearest_data$CompName[csn_nearest_data$CompName == "OC2.88" ] = "OC2"
csn_nearest_data$CompName[csn_nearest_data$CompName == "OC3.88" ] = "OC3"
csn_nearest_data$CompName[csn_nearest_data$CompName == "OC4.88" ] = "OC4"
imp_nearest_data = join(imp_nearest_data, nearest_csn_imp_sites_imp)
imp_nearest_data$Dataset = csn_nearest_data$Dataset = NULL
colnames(imp_nearest_data) = c("SiteCode.IMP", "Date", "CompName", "IMPROVE", "class", "Method.IMP", "Qualifier.IMP", "SameSite")
colnames(csn_nearest_data) = c("SiteCode.CSN", "Date", "CompName", "CSN", "class", "Method.CSN", "Qualifier.CSN", "SameSite")
imp_csn_nearest_data = join(imp_nearest_data, csn_nearest_data)
head(imp_csn_nearest_data)
summary(imp_csn_nearest_data)

# write.csv(imp_csn_nearest_data, "IMPROVE_CSN_Nearest_side-by-side_comparison.csv")

###### component - nearest sites - by date & by measurement method #######
imp_csn_nearest_data = read.csv("IMPROVE_CSN_Nearest_side-by-side_comparison.csv")
imp_csn_nearest_data$X = NULL
imp_csn_nearest_data$Date = as.Date(imp_csn_nearest_data$Date)

## elements from both datasets tested by XRF till 2015.11
## ions alwasy IC
## OC, EC from CSN using IMPROVE protocal from 2018.10
## subset data with identical species measurement method
imp_csn_nearest_data = 
  subset(imp_csn_nearest_data,
         (Date < as.Date("2015-11-01") & class == "Element") |
           (Date > as.Date("2018-10-01") & class == "OC/EC subgroup") |
           class == "Ion")
dim(imp_csn_nearest_data)

# remove NA
imp_csn_nearest_data = na.omit(imp_csn_nearest_data)
dim(imp_csn_nearest_data)

sapply(imp_csn_nearest_data, class)
# Calculate correlation by component
correl.mean_by_comp <- imp_csn_nearest_data %>%
  group_by(CompName) %>%
  # there is summarize in plyr{}, need to specify here
  dplyr::summarize(correlation = cor(IMPROVE, CSN, 
                                     use = "complete.obs",
                                     method = "pearson"),
                   mean.IMPROVE = mean(IMPROVE),
                   mean.CSN = mean(CSN),
                   median.IMPROVE = median(IMPROVE),
                   median.CSN = median(CSN))

correl.mean_by_comp$mean.diff = 
  label_percent()(
    abs((correl.mean_by_comp$mean.IMPROVE - 
           correl.mean_by_comp$mean.CSN)/
          correl.mean_by_comp$mean.IMPROVE))

correl.mean_by_comp$median.diff = 
  label_percent()(
    abs((correl.mean_by_comp$median.IMPROVE - 
           correl.mean_by_comp$median.CSN)/
          correl.mean_by_comp$median.IMPROVE))

correl.mean_by_comp$mean.diff1 = 
  abs((correl.mean_by_comp$mean.IMPROVE - 
           correl.mean_by_comp$mean.CSN)/
          correl.mean_by_comp$mean.IMPROVE)

correl.mean_by_comp$median.diff1 = 
    abs((correl.mean_by_comp$median.IMPROVE - 
           correl.mean_by_comp$median.CSN)/
          correl.mean_by_comp$median.IMPROVE)

write.csv(correl.mean_by_comp, "CSN_IMPROVE_Nearest7sites_Compnenet_summary.csv")
summary(correl.mean_by_comp)

imp_csn_nearest_low_S = subset(imp_csn_nearest_data, 
                               (CompName == "S" &
                                 CSN < 0.15 & 
                                 IMPROVE > 0.25) |
                                 (CompName == "SO4" &
                                    CSN < 0.15 & 
                                    IMPROVE > 1) )
write.csv(imp_csn_nearest_low_S, "S_SO4_low.CSN_high.IMPORVE.csv")

###### component - nearest sites - by date & by measurement method - Plotting #######

theme.comp = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                   strip.text.x = element_text(size = 14, colour = "grey25", angle = 0),
                   axis.title.x = element_text(color="grey25", size = 18, vjust=0, margin=margin(0,0,0,300)), 
                   axis.title.y = element_text(color="grey25", size = 18, vjust=1, margin=margin(0,2,0,0)),
                   axis.text.x = element_text(color="grey25", size = 16, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                   axis.text.y = element_text(color="grey25", size = 16, angle = 0, hjust = 0.5))

element.group = unique(imp_csn_nearest_data$CompName[imp_csn_nearest_data$class == "Element"])
length(element.group)
#ion.group = unique(imp_csn_nearest_data$CompName[imp_csn_nearest_data$class == "Ion"])
ion.group = c("NO3", "SO4")
oc.ec.group = unique(imp_csn_nearest_data$CompName[imp_csn_nearest_data$class == "OC/EC subgroup"])
# element.group = c("Al", "As", "Cr", "Cu", "Fe", "K", "Mn", "Na", "Ni", "V", "Zn")
# ion.group = c("Cl-", "NO3", "SO4")
# oc.ec.group = c("EC1", "EC2", "EC3", "OC1", "OC2", "OC3", "OC4")


element.1 = ggplot(subset(imp_csn_nearest_data, CompName %in% element.group[1:6]), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
element.2 = ggplot(subset(imp_csn_nearest_data, CompName %in% element.group[7:12]), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
element.3 = ggplot(subset(imp_csn_nearest_data, CompName %in% element.group[13:18]), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
element.4 = ggplot(subset(imp_csn_nearest_data, CompName %in% element.group[19:24]), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp

# multiplot(element.1, element.2, element.3, element.4, cols = 1)

ggplot(subset(imp_csn_nearest_data, CompName %in% c("Al", "Na", "S", "Si")),
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle("Elements test with EDXRF only in CSN") +
  theme.comp

ggplot(subset(imp_csn_nearest_data, CompName %in% c("Br", "Ca", "K")),
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle("Elements test with XRF or more than EDXRF in CSN") +
  theme.comp

ion.dis = ggplot(subset(imp_csn_nearest_data, CompName %in% c("Cl-", "NO3", "SO4")), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
# subset(imp_csn_nearest_data, CSN < 0.5 & IMPROVE > 1 & CompName == "SO4")

ec.dis = ggplot(subset(imp_csn_nearest_data, CompName %in% c("EC1", "EC2", "EC3")), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
oc.dis = ggplot(subset(imp_csn_nearest_data, CompName %in% c("OC1", "OC2", "OC3", "OC4")), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp

imp_csn_nearest_2015on = subset(imp_csn_nearest_data, Date >= as.Date("2015-01-01"))
imp_csn_nearest_2015on$Period[imp_csn_nearest_2015on$Date < as.Date("2018-01-01")] = "2015-18"
imp_csn_nearest_2015on$Period[imp_csn_nearest_2015on$Date >= as.Date("2018-01-01")] = "2018-20"

ggplot(subset(imp_csn_nearest_2015on, CompName %in% c("EC1", "EC2", "EC3")), 
       aes(IMPROVE, CSN, color = Period)) + 
  geom_point(size = 1.5, alpha = 0.5) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme_bw() +
  theme.comp

ggplot(subset(imp_csn_nearest_2015on, CompName %in% c("OC1", "OC2", "OC3", "OC4")), 
       aes(IMPROVE, CSN, color = Period)) + 
  geom_point(size = 1.5, alpha = 0.5) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme_bw() +
  theme.comp

#### EC/OC subgroups between nearest sites ####
imp_csn_nearest_C = subset(imp_csn_nearest_data, class == "OC/EC subgroup")
imp_csn_nearest_C_noNA = subset(imp_csn_nearest_C, !is.na(Method.CSN))

###### 1. in CSN, EC-subgroup using STN-TOT & OC with IMPROVE_TOR ######
## these OC/EC-subs in CSN was shown as OC1-4, EC1-3
imp_csn_nearest_C_csnTOT = subset(imp_csn_nearest_C_noNA, Method.CSN == "STN-TOT")
imp_csn_nearest_C_csnTOR = subset(imp_csn_nearest_C_noNA, Method.CSN == "IMPROVE_TOR")
range(imp_csn_nearest_C_csnTOT$Date) # OC, from "2010-01-02" to "2010-04-29"
range(imp_csn_nearest_C_csnTOR$Date) # EC, from "2015-11-20" to "2020-12-29"

ggplot(subset(imp_csn_nearest_C_csnTOT, CompName %in% c("OC1", "OC2", "OC3", "OC4")), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle(paste("OC-Subgroup_CSN-TOT_vs._IMPROVE-TOR: ", min(imp_csn_nearest_C_csnTOT$Date), "~", max(imp_csn_nearest_C_csnTOT$Date))) +
  theme.comp

ggplot(subset(imp_csn_nearest_C_csnTOR, CompName %in% c("EC1", "EC2", "EC3")), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ CompName, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle(paste("EC-Subgroup_CSN-TOR_vs._IMPROVE-TOR: ", min(imp_csn_nearest_C_csnTOR$Date), "~", max(imp_csn_nearest_C_csnTOR$Date))) +
  theme.comp

###### 2. in CSN, EC/OC-subgroup with other methods ######
csn_near_C = subset(csn_data_compare, SiteCode %in% nearest.csn.sites & class == "OC/EC subgroup")
csn_near_C_other = subset(csn_near_C, !(CompName %in% c("EC1", "EC2", "EC3", "OC1", "OC2", "OC3", "OC4")))
unique(csn_near_C_other$CompName)

csn_near_OC = subset(csn_near_C, grepl("OC1", CompName, fixed = T) | 
                       grepl("OC2", CompName, fixed = T) | 
                       grepl("OC3", CompName, fixed = T) |
                       grepl("OC4", CompName, fixed = T))
csn_near_EC = subset(csn_near_C, grepl("EC1", CompName, fixed = T) |
                       grepl("EC2", CompName, fixed = T) | 
                       grepl("EC3", CompName, fixed = T) |
                       grepl("OPC", CompName, fixed = T))
write.csv(csn_near_OC, "nearest_IMP_CSN_OC.csv")
write.csv(csn_near_EC, "nearest_IMP_CSN_EC.csv")
# to be continued 2022-11-05

#### unmatched Elements between nearest sites ####
imp_csn_near_ele = subset(imp_csn_nearest_data, class == "Element")
imp_csn_near_ele_noNA = subset(imp_csn_near_ele, !is.na(Method.CSN))
imp_csn_near_ele_NA = subset(imp_csn_near_ele, is.na(Method.CSN))

# Al
imp_csn_near_Al = subset(imp_csn_near_ele_noNA, CompName == "Al")
min(imp_csn_near_Al$Date); max(imp_csn_near_Al$Date)
unique(imp_csn_near_Al$Method.CSN)

# Br
imp_csn_near_Br = subset(imp_csn_near_ele_noNA, CompName == "Br")
min(imp_csn_near_Br$Date); max(imp_csn_near_Br$Date)
unique(imp_csn_near_Br$Method.CSN)

ggplot(imp_csn_near_Br, 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, color = "red", alpha = 0.4) + 
  facet_grid(. ~ Method.CSN, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle("Br") +
  theme.comp

# Ca
imp_csn_near_Ca = subset(imp_csn_near_ele_noNA, CompName == "Ca")
min(imp_csn_near_Ca$Date); max(imp_csn_near_Ca$Date)
unique(imp_csn_near_Ca$Method.CSN)

ggplot(imp_csn_near_Ca, 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, color = "forestgreen", alpha = 0.6) + 
  facet_grid(. ~ Method.CSN, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle("Ca") +
  theme.comp

# K
imp_csn_near_K = subset(imp_csn_near_ele_noNA, CompName == "K")
min(imp_csn_near_K$Date); max(imp_csn_near_K$Date)
unique(imp_csn_near_K$Method.CSN)

ggplot(imp_csn_near_K, 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, color = "dodgerblue", alpha = 0.6) + 
  facet_grid(. ~ Method.CSN, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  ggtitle("K") +
  theme.comp

# Na
imp_csn_near_Na = subset(imp_csn_near_ele_noNA, CompName == "Na")
min(imp_csn_near_Na$Date); max(imp_csn_near_Na$Date)
unique(imp_csn_near_Na$Method.CSN)

# S
imp_csn_near_S = subset(imp_csn_near_ele_noNA, CompName == "S")
min(imp_csn_near_S$Date); max(imp_csn_near_S$Date)
unique(imp_csn_near_S$Method.CSN)

# Si
imp_csn_near_Si = subset(imp_csn_near_ele_noNA, CompName == "Si")
min(imp_csn_near_Si$Date); max(imp_csn_near_Si$Date)
unique(imp_csn_near_Si$Method.CSN)

imp_csn_near_csnXRF = subset(imp_csn_near_ele_noNA, Method.CSN == "XRF")
unique(imp_csn_near_csnXRF$CompName)

imp_csn_near_csnXRF = subset(imp_csn_near_ele_noNA, Method.CSN == "EDXRF")
unique(imp_csn_near_csnXRF$CompName)

#### Ions between nearest sites ####
imp_csn_near_ion = subset(imp_csn_nearest_data, class == "Ion")
imp_csn_near_ion_noNA = subset(imp_csn_near_ion, !is.na(Method.CSN))
imp_csn_near_ion_NA = subset(imp_csn_near_ion, is.na(Method.CSN))

# Cl-
imp_csn_near_Cl_ion = subset(imp_csn_near_ion_noNA, CompName == "Cl-")
min(imp_csn_near_Cl_ion$Date); max(imp_csn_near_Cl_ion$Date)
unique(imp_csn_near_Cl_ion$Method.CSN)

# NO3
imp_csn_near_NO3 = subset(imp_csn_near_ion_noNA, CompName == "NO3")
min(imp_csn_near_NO3$Date); max(imp_csn_near_NO3$Date)
unique(imp_csn_near_NO3$Method.CSN)

# SO4
imp_csn_near_SO4 = subset(imp_csn_near_ion_noNA, CompName == "SO4")
min(imp_csn_near_SO4$Date); max(imp_csn_near_SO4$Date)
unique(imp_csn_near_SO4$Method.CSN)

ggplot(imp_csn_near_SO4, 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8, color = "dodgerblue", alpha = 0.8) + 
  facet_grid(. ~ Method.CSN, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp

##########################################################################################
########### Preparing data for PMF ###########
##########################################################################################
#### 1. IMPROVE data preparation ####
head(imp_data_compare)
#imp_data_try = select(imp_data_compare, Dataset, State, SiteCode, Date, Qualifier, CompName, Val)[1:192, ]
#nrow(imp_data_try)/length(unique(imp_data_try$CompName))
#imp_data_try_sp = imp_data_try %>% spread(CompName, Val)
#head(imp_data_try_sp)

imp_data_use = select(imp_data_compare, Dataset, State, SiteCode, Date, Qualifier, CompName, Val)
nrow(imp_data_use)/length(unique(imp_data_use$CompName))
setDT(imp_data_use)

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

# change data to data.table to be more efficient
setDT(imp_date_comp_site)
setDT(imp_data_use)

# combine the full date-site-component list with the concentration values
imp_date_comp_site_data = merge(imp_date_comp_site, imp_data_use, all.x = T)


####!!!JUMP to 1.4 - convert data to day-sitecode-all.compnents style for PMF!!


##### 1.1 detect the duplicated values from POC!!! Parameter Occurrence Code  ####
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

##### 1.2 check if the strange duplicated component values from the same date-site-component groups exist in original dataset #####
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

##### 1.3 compare concentrations between POC from same date-site-component groups ####
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
##### 1.4 calculate component concentrations for those with POC = 2 ####
imp_data_use$POC = imp_data$POC
head(imp_data_use)
head(imp_dup_date_comp_site)
dim(imp_dup_date_comp_site)
dim(imp_data_use)
sapply(imp_dup_date_comp_site, class)

# estimate the Val for those with duplicated measurements
imp_dup_date_comp_site$Val.opc = (imp_dup_date_comp_site$Val.dup + 
                                    imp_dup_date_comp_site$Val.NONEdup)/2

# with NA in one of the duplicated measurement, use the none NA value
imp_dup_date_comp_site$Val.opc = ifelse(imp_dup_date_comp_site$Val.dup == -999, 
                                        imp_dup_date_comp_site$Val.NONEdup, 
                                        imp_dup_date_comp_site$Val.opc)
imp_dup_date_comp_site$Val.opc = ifelse(imp_dup_date_comp_site$Val.NONEdup == -999, 
                                        imp_dup_date_comp_site$Val.dup, 
                                        imp_dup_date_comp_site$Val.opc)
imp_dup_date_comp_site_conc = select(imp_dup_date_comp_site, 
                                     Date, SiteCode, CompName, Val.opc)

# match dataset for duplicates with whole dataset, and determine the final Val
imp_data_poc1 = join(imp_data_use, imp_dup_date_comp_site_conc)
# imp_data_poc1.1 = imp_data_poc1
imp_data_poc1$Val = ifelse(is.na(imp_data_poc1$Val.opc), 
                           imp_data_poc1$Val, 
                           imp_data_poc1$Val.opc)

# remove the duplicates
imp_data_poc1 = subset(imp_data_poc1, POC == 1)

# extract those with real records & match with dataset of full dates & sites
imp_poc_noNA = subset(imp_data_poc1, Val != -999)
table(imp_poc_noNA$Qualifier)

##### 1.5 deal with qualifier, and convert data to the style for PMF ####
# remove marked with given flags as NA
# the flags were selected based on description & the influence on Val
imp_poc_noNA = subset(imp_poc_noNA, 
                      !(Qualifier %in% c("M2", "V2", "V4", "V5")))
imp_poc_noNA$Qualifier = imp_poc_noNA$Val.opc = imp_poc_noNA$POC = 
  imp_poc_noNA$Dataset = NULL

# match those of real records with dataset of full dates & sites
imp_date_comp_site_noDup_noNA = merge(imp_date_comp_site, imp_poc_noNA, all.x = T)
dim(imp_date_comp_site_noDup_noNA)

imp_date_comp_site_noDup_noNA$site.date = paste(imp_date_comp_site_noDup_noNA$Date,
                                                          imp_date_comp_site_noDup_noNA$State,
                                                          imp_date_comp_site_noDup_noNA$SiteCode)
imp_date_comp_site_noDup_noNA$Date = imp_date_comp_site_noDup_noNA$State = 
  imp_date_comp_site_noDup_noNA$SiteCode = NULL

head(imp_date_comp_site_noDup_noNA)
# write.csv(imp_date_comp_site_noDup_noNA, "IMPROVE_Component_to_be_spread.csv")

# imp_date_comp_site_noDup_noNA = read.csv("IMPROVE_Component_to_be_spread.csv")
imp_date_comp_site_noDup_noNA$X = NULL
setDT(imp_date_comp_site_noDup_noNA)
# spread the dataset to day-every.component style
imp_daily_comp = imp_date_comp_site_noDup_noNA %>% spread(CompName, Val)
head(imp_daily_comp)

# double check if there are duplicated date&site group
imp_daily_comp$dup.site.date = duplicated(imp_daily_comp$site.date)
summary(imp_daily_comp$dup.site.date)
imp_daily_comp$dup.site.date = NULL

# get the separate column of site and date info
imp_site_date = data.frame(select(imp_daily_comp, site.date))
imp_site_date_sep = imp_site_date %>% 
  separate(site.date, c("Date", "State", "SiteCode"),  
           sep = "\\s+")
head(imp_site_date_sep)
imp_daily_comp_use = cbind(imp_site_date_sep, imp_daily_comp)
imp_daily_comp_use$site.date = NULL

# No longer save qualifier records since 2023
#imp_daily_comp_use_Q.Na = subset(imp_daily_comp_use, Qualifier == "NA")
#summary(imp_daily_comp_use_Q.Na)
#imp_daily_comp_use = subset(imp_daily_comp_use, Qualifier != "NA")

# already exist since 2023 
#imp_meta_sites = read.csv("IMPROVE metadata 196 sample sites info 2010-20.csv")
#imp_state_site = select(imp_meta_sites, State, SiteCode)

#imp_daily_comp_use = merge(imp_daily_comp_use, imp_state_site, all.x = T)
#imp_daily_comp_use$site.date = imp_daily_comp_use$State
#imp_daily_comp_use$State = NULL
#colnames(imp_daily_comp_use)[4] = "State"

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
# from row 1144, no longer every-3-day measurement
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

# prepare the file with all sampling dates for each site and each PM species
csn_date_comp_site_ref = data.frame(Date = rep(csn_date, each = length(csn_Comp) * length(csn_site)), 
                                CompName = rep(csn_Comp, each = length(csn_site)), 
                                SiteCode = rep(csn_site))
dim(csn_date_comp_site_ref)
head(csn_date_comp_site_ref)
unique(csn_date_comp_site_ref$Date)[1:15]
setDT(csn_date_comp_site_ref)

#### !!! jump to 2.3 to prepare data for PMF !!!

##### 2.1.1 get the monthly MDL for PM species #####
# according to "CSN_contract_transition_data_advisory_final_508.pdf", 
## MDL in CSN were evaluated each month instead of each year

# replace negative values with NA 
csn_data$MDL.not.neg = csn_data$MDL
csn_data$MDL.not.neg = replace(csn_data$MDL.not.neg, 
                               which(csn_data$MDL.not.neg < 0), 
                               NA)

# get the montly average MDL for each site
csn_mdl = ddply(csn_data, .(SiteCode, CompName, year, month), 
                summarise,
                MDL = mean(MDL.not.neg, na.rm = T))
csn_data$MDL.not.neg = NULL

csn_mdl$site.date = paste(csn_mdl$SiteCode, csn_mdl$year, csn_mdl$month)
csn_mdl_gathered = select(csn_mdl, CompName, MDL, site.date)

# convert to site-day-allSpecies format
csn_mdl_spread = csn_mdl_gathered %>% spread(CompName, MDL)
head(csn_mdl_spread)
nrow(csn_mdl_spread)
csn_SiteDate = data.frame(select(csn_mdl_spread, site.date))
csn_SiteDate = csn_SiteDate %>% 
  separate(site.date, c("SiteCode", "year", "month"),  
           sep = "\\s+")
csn_mdl_use = cbind(csn_SiteDate, csn_mdl_spread)
csn_mdl_use$site.date = NULL

# get the MDL for PM2.5 before & after year 2015
csn_mdl_before_2015 = subset(csn_mdl_use, 
                             year <= 2015)
csn_mdl_after_2015 = subset(csn_mdl_use, 
                            year > 2015)

##### 2.1.2-Version A: monthly MDL for PM species - OC/EC #####

# variables not to be used for PMF
csn_bfr_remove = c("EC1.unadjusted.88", "EC2.unadjusted.88", "EC3.unadjusted.88",  
                   "OC1.unadjusted.88", "OC2.unadjusted.88", "OC3.unadjusted.88", "OC4.unadjusted.88",
                   "OC1.88", "OC2.88", "OC3.88", "OC4.88",
                   "EC1.88", "EC2.88", "EC3.88", 
                   "OC1", "OC2", "OC3", "OC4",
                   "EC1", "EC2", "EC3", 
                   "EC.unadjusted.88", "EC.88", # "EC.TOR.unadjust.88"
                   "OC.unadjusted.88", "OC", # "OC.TOR.unadjusted.88"
                   "OPC.88", "OPC.unadjusted.88", # "OP.TOR.unadjusted.88"
                   "EC.TOR.88", "OC.88", "OPC.TOR.88",
                   "EC", "OC.unadjusted", "OPC")
# "Ag", "Ba", "Cd", "Ce", "Co", "Cs", "In", 
# "K+", "Na+", "NH4+", "Sb", "Sn", 
# "Rb", "Zr", "Cl-")
# Above sepecies were kept since 2023.02

csn_aft_remove = c("EC1.unadjusted.88", "EC2.unadjusted.88", "EC3.unadjusted.88",  
                   "OC1.unadjusted.88", "OC2.unadjusted.88", "OC3.unadjusted.88", "OC4.unadjusted.88",
                   "OC1.88", "OC2.88", "OC3.88", "OC4.88",
                   "EC1.88", "EC2.88", "EC3.88", 
                   "OC1", "OC2", "OC3", "OC4",
                   "EC1", "EC2", "EC3", 
                   "EC.unadjusted.88", "EC.TOR.unadjust.88", "EC.88",
                   "OC.unadjusted.88", "OC.TOR.unadjusted.88", "OC", 
                   "OPC.88", "OP.TOR.unadjusted.88", "OPC.unadjusted.88",
                   "EC", "OC.unadjusted", "OPC")
# "Ag", "Ba", "Cd", "Ce", "Co", "Cs", "In", 
# "K+", "Na+", "NH4+", "Sb", "Sn", 
# "Rb", "Zr", "Cl-")
# Above sepecies were kept since 2023.02

csn_mdl_before_2015[ ,csn_bfr_remove] <- list(NULL)
csn_mdl_after_2015[ ,csn_aft_remove] <- list(NULL)

# Accept.PM2.5 as PM concentrations before 2015, and PM2.5RC for after
csn_mdl_before_2015$PM25 = csn_mdl_before_2015$Accept.PM2.5
csn_mdl_after_2015$PM25 = csn_mdl_after_2015$PM2.5RC
csn_mdl_before_2015$Accept.PM2.5 = csn_mdl_before_2015$PM2.5RC = 
  csn_mdl_after_2015$Accept.PM2.5 = csn_mdl_after_2015$PM2.5RC = NULL

# reset colnames, OC/EC were selected based on test protocol
colnames(csn_mdl_before_2015)[c(18, 30, 31)] = c("EC", "OC", "OP")
# "EC.TOR.unadjust.88" "OC.TOR.unadjusted.88" "OP.TOR.unadjusted.88" used till 2015
colnames(csn_mdl_after_2015)[c(18, 30, 31)] = c("EC", "OC", "OP")
# "EC.TOR.88" "OC.88" "OPC.TOR.88" are used after 2015 

# check if columns from data before & after 2015 match
summary(colnames(csn_mdl_after_2015) == colnames(csn_mdl_before_2015))

# combine and output final MDL file for CSN
csn_mdl_combine = rbind(csn_mdl_before_2015, csn_mdl_after_2015)
summary(csn_mdl_combine)
dim(csn_mdl_combine)
sapply(csn_mdl_combine, class)


##### 2.1.2-Version B: monthly MDL for PM species - C-Subgroup #####
# remove multiple columns, name of which would be used
C.before.col.remove = c("EC", "EC1", "EC2", "EC3", 
                        "OC", "OC1", "OC2", "OC3", "OC4", 
                        "OP")
C.after.col.remove = c("EC", "OP", 
                       "OC", "OC1", "OC2", "OC3", "OC4")
csn_mdl_before_2015[ , C.before.col.remove] <- list(NULL)
csn_mdl_after_2015[ , C.after.col.remove] <- list(NULL)

# Changed selected to OC, EC & OP for PMF
csn_mdl_before_2015 = plyr::rename(csn_mdl_before_2015, 
                                   c("EC1.unadjusted.88" = "EC1", 
                                     "EC2.unadjusted.88" = "EC2",
                                     "EC3.unadjusted.88" = "EC3",
                                     "OC1.unadjusted.88" = "OC1",
                                     "OC2.unadjusted.88" = "OC2", 
                                     "OC3.unadjusted.88" = "OC3",
                                     "OC4.unadjusted.88" = "OC4",
                                     "Accept.PM2.5" = "PM25",
                                     "EC.TOR.unadjust.88" = "EC", 
                                     "OC.TOR.unadjusted.88" = "OC",
                                     "OP.TOR.unadjusted.88" = "OP"))

csn_mdl_after_2015 = plyr::rename(csn_mdl_after_2015, 
                                  c("OC1.88" = "OC1",
                                    "OC2.88" = "OC2", 
                                    "OC3.88" = "OC3",
                                    "OC4.88" = "OC4",
                                    "PM2.5RC" = "PM25",
                                    "EC.TOR.88" = "EC", 
                                    "OC.88" = "OC",
                                    "OPC.TOR.88" = "OP"))

# get shared colnames
conc_common_cols <- intersect(names(csn_mdl_before_2015), 
                              names(csn_mdl_after_2015))

# Subset data frames using shared colnames
csn_mdl_before_2015 <- csn_mdl_before_2015[, conc_common_cols]
csn_mdl_after_2015 <- csn_mdl_after_2015[, conc_common_cols]

# combine the dataset
csn_mdl_OrigOrder = rbind(csn_mdl_before_2015, csn_mdl_after_2015)

# Remove variables not needed for PMF (C-subgroups)
csn_remove = c("OC.unadjusted.88", "OC.unadjusted",
               "EC.unadjusted.88", "EC.88",
               "OPC.unadjusted.88", "OPC", "OP", "OPC.88" )
csn_mdl_OrigOrder[ ,csn_remove] <- list(NULL)
colnames(csn_mdl_OrigOrder)

# Move grouped columns to the right of the dataframe
# match(), select columns that match a specific pattern in their names
# everything(), include all remaining columns not selected by matches()
OC.EC = c("EC", "EC1", "EC2", "EC3", 
          "OC", "OC1", "OC2", "OC3", "OC4")
ions = c("Na.", "K.", "NH4.", "NO3", "SO4")
csn_mdl_combine = csn_mdl_OrigOrder %>%
  select(!(matches(OC.EC)), 
         everything())
csn_mdl_combine = csn_mdl_combine %>%
  select(!(matches(ions)), 
         everything())

# Replace other columns
csn_mdl_combine = csn_mdl_combine %>% relocate(PM25, .after = SO4)
csn_mdl_combine = csn_mdl_combine %>% relocate(SiteCode, .before = year)

colnames(csn_mdl_combine)
dim(csn_mdl_combine)

##### 2.1.3 Replace NAs in monthly MDL with annual median or others for OC/EC #####
# get the median MDL of each year for each site
csn_mdl_median_annual = csn_mdl_combine  %>% 
  group_by(SiteCode, year, month) %>%
  summarise_if(is.numeric,
    median, na.rm = T)
dim(csn_mdl_median_annual)

# get the 10th percentile MDL of each year for each site
#csn_mdl_10th_annual = csn_mdl_combine  %>% 
#  dplyr::group_by(SiteCode, year) %>%
#  dplyr::summarise(
#    dplyr::across(
#    where(is.numeric), 
#    ~ quantile(., 0.1, na.rm = T)))
#dim(csn_mdl_10th_annual)

# For OC/EC & subgroups, replaced with MDL of 2011-17 & 2018-20
csn_mdl_OCEC = select(csn_mdl_combine, 
                      contains("EC") |
                        contains("OC") |
                        contains("year") |
                        contains("month") ) 
csn_mdl_OCEC$Duration = "2011-17"
csn_mdl_OCEC$Duration[csn_mdl_OCEC$year >= 2018] = "2018-20"

csn_mdl_C_median = csn_mdl_OCEC  %>% 
  group_by(SiteCode, Duration) %>%
  summarise_if(is.numeric,
               median, na.rm = T)
summary(csn_mdl_C_median)

# detect site-year groups of data with NA in MDL
summary(csn_mdl_median_annual)
csn_mdl_med_annual_NA = subset(csn_mdl_median_annual, is.na(Al))
unique(csn_mdl_med_annual_NA$SiteCode)
## detect NA exists for 22 site-year groups of data (ions & elements) from 7 sites
### "132950002" "171190024" "20904101"  "390350065" "390350076" "560210100" "720210010"
### when including C-subgroups, the number increased?

# detect sites of data with NA in MDL
csn_mdl_median_site = csn_mdl_combine  %>% 
  group_by(SiteCode) %>%
  summarise_if(is.numeric,
               median, na.rm = T)
summary(csn_mdl_median_site)
unique(subset(csn_mdl_median_site, is.na(Al))$SiteCode)
## detect NA exists for 5 sites across the whole sampling period
### "132950002" "20904101"  "390350065" "390350076" "560210100"

# those with high OC MDL
csn_mdl_highOC = subset(csn_mdl_combine, SiteCode %in% 
                          unique(subset(csn_mdl_combine, OC3 > 1.5)$SiteCode))
# "220330009" "20904101" 

# prepare datasets for replace
csn_mdl_combine_date = select(csn_mdl_OCEC, SiteCode, year, month)
csn_mdl_annual_combine = merge(csn_mdl_combine_date, 
                               csn_mdl_median_annual,
                               all.x = T) 
# extend it to yearly
csn_mdl_C_annual = select(csn_mdl_OCEC, SiteCode, Duration, year, month)
csn_mdl_C_median_annual = merge(csn_mdl_C_annual, 
                                csn_mdl_C_median, 
                                all.x = T)
csn_mdl_C_median_annual$Duration = NULL

csn_mdl_final = csn_mdl_combine
# reorder files for data matching
csn_mdl_annual_combine = csn_mdl_annual_combine[with(csn_mdl_annual_combine, 
                                                     order(SiteCode, year, month)), ]
csn_mdl_C_median_annual = csn_mdl_C_median_annual[with(csn_mdl_C_median_annual, 
                                                       order(SiteCode, year, month)), ]
csn_mdl_final = csn_mdl_final[with(csn_mdl_final, 
                                   order(SiteCode, year, month)), ]

# manually check if the rows matched for random selected sites
csn_mdl_C_cb_1 = subset(csn_mdl_C_median_annual, SiteCode == "132950002")
csn_mdl_cb_1 = subset(csn_mdl_annual_combine, SiteCode == "132950002")
csn_mdl_fl_1 = subset(csn_mdl_final, SiteCode == "132950002")
summary(rownames(csn_mdl_cb_1) == rownames(csn_mdl_fl_1) & 
          rownames(csn_mdl_C_cb_1) == rownames(csn_mdl_fl_1) )
View(csn_mdl_cb_1)
View(csn_mdl_fl_1)

csn_mdl_final = data.frame(csn_mdl_final)
csn_mdl_annual_combine = data.frame(csn_mdl_annual_combine)
csn_mdl_C_median_annual = data.frame(csn_mdl_C_median_annual)

# replace NAs for sites having MDL values with annual median
csn_mdl_final[is.na(csn_mdl_final)] = csn_mdl_annual_combine[is.na(csn_mdl_final)]
sapply(csn_mdl_final, class)
colnames(csn_mdl_final)

# replace all OC.EC columns by intersect columns
C_col_names <- intersect(colnames(csn_mdl_final), 
                         colnames(csn_mdl_C_median_annual)[4:12])
csn_mdl_final[OC.EC] <- csn_mdl_C_median_annual[OC.EC]

# convert species concentration to numeric
csn_mdl_final <- csn_mdl_final %>% 
  mutate_if(is.character, as.numeric)
summary(csn_mdl_final)

# replace NAs for sites having no MDL values with mean of all sites
indx <- which(is.na(csn_mdl_final), arr.ind = T)
csn_mdl_final[indx] <- colMeans(csn_mdl_final, na.rm = T)[indx[, 2]]
csn_mdl_final$SiteCode = as.character(csn_mdl_final$SiteCode)

# output MDL file
# write.csv(csn_mdl_final, "CSN_MDL_monthly_2023.02.27.csv")
write.csv(csn_mdl_final, "CSN_MDL_C-Sub_monthly_2023.05.csv")

##### 2.2 check POC repeatability (duplicated measurements) #####
csn_data$Method = csn_data$Status = csn_data$Elevation = 
  csn_data$Qualifier1 = csn_data$Qualifier2 = csn_data$Qualifier3 = 
  csn_data$Unit = csn_data$Unc = csn_data$MDL = csn_data$SampleDuration = 
  csn_data$year = csn_data$month = csn_data$day = NULL
csn_data$site.date = paste(csn_data$Date, 
                           csn_data$State,
                           csn_data$SiteCode, 
                           csn_data$CompName)

# check the type and number of different POC in CSN
table(csn_data$POC) # POC = 5, 6, or 7

csn_data_poc5 = subset(csn_data, POC == 5)
csn_data_poc6 = subset(csn_data, POC == 6)
csn_data_poc7 = subset(csn_data, POC == 7)

# check if there could be >2 measurements for one species from the same date & site
summary(csn_data_poc7$site.date %in% csn_data_poc6$site.date)
summary(csn_data_poc6$site.date %in% csn_data_poc5$site.date)
summary(csn_data_poc7$site.date %in% csn_data_poc5$site.date)
# the results suggest there are up to 2 measurements for a given date/site/specie group

# prepare data for comparison between duplicated measurement 
csn_data_poc5_dup_use = select(csn_data_poc5, site.date, Val,  
                               Collection, Analysis, Qualifier)
csn_data_poc5_dup_use$poc5 = "Y"
csn_data_poc7_dup_use = select(csn_data_poc7, site.date, Val, 
                               Collection, Analysis, Qualifier)
csn_data_poc7_dup_use$poc7 = "Y"

# duplication type 1, POC 5 & 6
csn_poc56_dup = merge(csn_data_poc6, csn_data_poc5_dup_use, by = "site.date", all.x = T)
csn_poc56_dup = subset(csn_poc56_dup, poc5 == "Y")
dim(csn_poc56_dup)
# check if collection methods and qualifiers for the duplicates are same
summary(csn_poc56_dup$Collection.x == csn_poc56_dup$Collection.y)
summary(csn_poc56_dup$Analysis.x == csn_poc56_dup$Analysis.y)
summary(csn_poc56_dup$Qualifier1.x == csn_poc56_dup$Qualifier1.y)
summary(csn_poc56_dup$Qualifier2.x == csn_poc56_dup$Qualifier2.y)
summary(csn_poc56_dup$Qualifier3.x == csn_poc56_dup$Qualifier3.y)

# duplication type 2, POC 6 & 7
csn_poc67_dup = merge(csn_data_poc6, csn_data_poc7_dup_use, by = "site.date", all.x = T)
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
                    site.date, POC, 
                    Date, SiteCode, 
                    CompName, ParamCode, Val.dup)
dup_info56$POC = 56

csn_poc67_dup$Val.dup = (csn_poc67_dup$Val.x + csn_poc67_dup$Val.y)/2
dup_info67 = select(csn_poc67_dup,  
                    site.date, POC, 
                    Date, SiteCode, 
                    CompName, ParamCode, Val.dup)
dup_info67$POC = 67

# double check if there are duplicates 
dup_info_56_67 = rbind(dup_info56, dup_info67)
dup_info_56_67$dup = duplicated(dup_info_56_67$site.date)
summary(dup_info_56_67$dup)
dup_info_56_67 = select(dup_info_56_67, 
                        site.date, 
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
csn_date_comp_site_measured = csn_normalDateSamp

# import the duplicated measurement info generated from step 2.2
dup_info_56_67 = read.csv("CSN_Duplicates_POC.csv")
dup_info_56_67$X = NULL
setDT(dup_info_56_67)
  
csn_date_comp_site_measured$site.date = paste(csn_date_comp_site_measured$Date, 
                                              csn_date_comp_site_measured$State,
                                              csn_date_comp_site_measured$SiteCode, 
                                              csn_date_comp_site_measured$CompName)
dim(csn_date_comp_site_measured)
# check the existance of duplicated date-site-component groups
summary(duplicated(csn_date_comp_site_measured$site.date))

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
length(unique(csn_date_comp_no_dup$site.date)) + # 6425441
  length(unique(csn_date_comp_with_dup$site.date)) == # 179081
  length(unique(csn_date_comp_site_measured$site.date)) # 6604522
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
# double check if there are still duplicated records
summary(duplicated(csn_date_comp_site_data$site.date))
# double check if there are NAs 
summary(is.na(csn_date_comp_site_data$site.date))

##### 2.4 deal with low & qualifier, merge with reference (include all date-site-component groups) #####

# the detected non-negative Val is 5e-06 , set all non-positive Val to 0.000009
# the 0.000009 will be transferred to MDL/2 for each component later
csn_date_comp_site_data$Val[csn_date_comp_site_data$Val <= 5e-06] = 0.000009

# set data with given qualifier as NA
csn_date_comp_site_data$Val[grepl("1", csn_date_comp_site_data$Qualifier) |
                              grepl("6", csn_date_comp_site_data$Qualifier) |
                              grepl("QX", csn_date_comp_site_data$Qualifier)] = NA

# set data with given qualifier as 0.000009 for later convertion to MDL/2
csn_date_comp_site_data$Val[grepl("MD", csn_date_comp_site_data$Qualifier) |
                              grepl("ND", csn_date_comp_site_data$Qualifier) |
                              grepl("SQ", csn_date_comp_site_data$Qualifier)] = 0.000009

# randomly check one site that was detected as have many qualifiers
data_site_manyNA_bf_qualifier = subset(csn_data_compare, SiteCode == 11130003)
data_site_manyNA_aft_qualifier = subset(csn_date_comp_site_data, SiteCode == 11130003)

# save the qualifier info
csn_qualifier = select(csn_date_comp_site_data, 
                       Date, SiteCode, CompName, Qualifier)
dim(csn_qualifier)
write.csv(csn_qualifier, "CSN_Qulifier_no_missing_Val.csv")

# generating information related qualifier
csn_qualifier = select(csn_qualifier, Date, SiteCode, Qualifier)
csn_qualifier_info = subset(csn_qualifier,
                            grepl("I", Qualifier, fixed = T))
csn_qualifier_info = csn_qualifier_info[!duplicated(csn_qualifier_info), ]

csn_qualifier_info = csn_qualifier_info[with(csn_qualifier_info, 
                                             order(SiteCode, Date)), ]
csn_qualifier_info$Date = as.Date(csn_qualifier_info$Date)
write.csv(csn_qualifier_info, "CSN_SiteCode_Information_Qualifier.csv")

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

(length(unique(csn_date_comp_site_data$site.date)) + 1) ==
  length(unique(csn_date_with_missing$site.date))
length(unique(csn_date_with_missing$site.date))

csn_site_date_qualifier_sep = 
  select(csn_date_with_missing, site.date) %>% 
  separate(site.date, c("Date", "State", "SiteCode", "CompName"),  
           sep = "\\s+")

# prepare the date-site-compname list with no NA
csn_site_date_qualifier_sep_noNA = subset(csn_site_date_qualifier_sep,
                                          !is.na(SiteCode))
csn_site_date_qualifier_sep_noNA$Date = 
  as.Date(csn_site_date_qualifier_sep_noNA$Date)
csn_site_date_qualifier_sep_noNA$SiteCode = 
  as.integer(csn_site_date_qualifier_sep_noNA$SiteCode)
csn_site_date_qualifier_sep_noNA = 
  csn_site_date_qualifier_sep_noNA[!duplicated(
    csn_site_date_qualifier_sep_noNA), ]
dim(csn_site_date_qualifier_sep_noNA)

# merge the date-site-compname list with no NA with component values
csn_date_miss  = merge(csn_date_with_missing, 
                       csn_site_date_qualifier_sep_noNA,
                       all.x = T)
# setDT(csn_date_miss)
summary(csn_date_miss)
summary(is.na(csn_date_miss$State))
summary(is.na(csn_date_miss$CompName))
# summary(is.na(csn_date_miss$Qualifier))

## Below is no longer needed, cause we have added state in the site.date info
# csn_date_miss$site.date = NULL
# csn_date_miss$site.date = paste(csn_date_miss$Date, csn_date_miss$State, 
#                                csn_date_miss$SiteCode, csn_date_miss$Qualifier)

# check if qualifiers for same day-site group are same
# csn_miss_sample = subset(csn_date_with_missing, 
#                          Date == csn_date_with_missing$Date[1] &
#                            SiteCode == csn_date_with_missing$SiteCode[1])
# length(unique(csn_miss_sample$site.date))

# write.csv(csn_date_miss, "CSN_Component_to_be_spread_2022-12.csv")
# write.csv(csn_date_miss, "CSN_Component_to_be_spread_2022-12_2023.02.csv")

##### 2.5 spreading - data for interpolation - before & after 2015 #####
library(psych) # corr.test{}
library(corrplot) # corrplot.mixed{}

# csn_date_miss = read.csv("CSN_Component_to_be_spread_2022-12.csv")
# csn_date_miss = read.csv("CSN_Component_to_be_spread_2022-12_2023.02.csv")
# csn_date_miss$X = NULL

# this old site.date info contains lots of NA, need a new one with complete list
csn_date_miss$site.date = NULL 
csn_date_miss$site.date = paste(csn_date_miss$Date, 
                                csn_date_miss$SiteCode)

csn_miss_to_spread = csn_date_miss
# check if there is still NA in the site.date info
summary(is.na(csn_miss_to_spread$site.date))

# remove extral columns
csn_miss_to_spread$Date = csn_miss_to_spread$SiteCode = 
  csn_miss_to_spread$State = NULL
head(csn_miss_to_spread)
dim(csn_miss_to_spread)

# spreading the data to daily-site-component style
csn_spread_try = csn_miss_to_spread[1:900, ]
csn_spread_try_sp = csn_spread_try %>% spread(CompName, Val)

csn_daily_comp = csn_miss_to_spread %>% spread(CompName, Val)
head(csn_daily_comp)
dim(csn_daily_comp)

### when running all steps in section 2.4
csn_daily_comp = merge(csn_daily_comp, csn_site_date_qualifier_sep_noNA, all.x = T)
summary(csn_daily_comp)
csn.d.col = ncol(csn_daily_comp)
csn_daily_comp_use = cbind(csn_daily_comp[, (csn.d.col-3):csn.d.col],
                           csn_daily_comp[, 2:(csn.d.col-4)])

### when jumping previous steps in 2.4, then we need new info of date, site, state
csn_site_date_daily = 
  select(csn_daily_comp, site.date) %>% 
  separate(site.date, c("Date", "SiteCode"),  
           sep = "\\s+")
# merge date, sitecode with species concentrations
csn_daily_comp_use = cbind(csn_site_date_daily, 
                           csn_daily_comp[, 2:ncol(csn_daily_comp)])
# match state info
csn_sites_use = select(csn_sites_use, State, SiteCode)
csn_daily_comp_use = merge(csn_daily_comp_use, 
                           csn_sites_use, 
                           all.x = T)

# write.csv(csn_daily_comp_use, "CSN_Component_with_missing.csv")
write.csv(csn_daily_comp_use, "CSN_Component_with_missing_2023.02.csv")

##### 2.6 check correlation & missing pattern #####

csn_daily_comp_use = read.csv("CSN_Component_with_missing_2023.02.csv")
csn_daily_comp_use$X = NULL
csn_daily_comp_use$Date = as.Date(csn_daily_comp_use$Date)



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

# check the relationship between species groups suggested by Phil on 2023-02-21

plot(csn_daily_comp_use$`K+`, csn_daily_comp_use$K, xlim = c(0,40), ylim = c(0,40))
plot(csn_daily_comp_use$`Na+`, csn_daily_comp_use$Na, xlim = c(0,8), ylim = c(0,8))
plot(csn_daily_comp_use$`Cl-`, csn_daily_comp_use$Cl, xlim = c(0,10), ylim = c(0,10))
plot(csn_daily_comp_use$SO4, csn_daily_comp_use$S, xlim = c(0,45), ylim = c(0,15))


# overall NA rate
p_miss_with_neg <- data.frame(unlist(lapply(
  csn_daily_comp[,2:ncol(csn_daily_comp)], 
  function(x) 
    sum(is.na(x))))/
    nrow(csn_daily_comp))
p_miss_with_neg$CompCode = rownames(p_miss_with_neg) 
colnames(p_miss_with_neg)[1] = "Missing_Rate"
rownames(p_miss_with_neg) = 1:nrow(p_miss_with_neg)
missing_rate = ggplot(p_miss_with_neg, aes(CompCode, Missing_Rate)) +
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
  csn_daily_before_2015[, 4:(ncol(csn_daily_before_2015)-1)], 
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
  csn_daily_after_2015[, 4:(ncol(csn_daily_after_2015)-1)], 
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
mice::md.pattern(csn_daily_before_2015[, 4:(ncol(csn_daily_before_2015)-1)], rotate.names = T) 
mice::md.pattern(csn_daily_after_2015[, 4:(ncol(csn_daily_after_2015)-1)], rotate.names = T) 


##### 2.7 data for interpolation - separate before & after 2015 #####
csn_before_2015 = csn_daily_before_2015
csn_after_2015 = csn_daily_after_2015

# select components with 80% missing 
comp_miss_before = unique(p_miss_with_neg_before2015$CompCode[p_miss_with_neg_before2015$Missing_Rate > 0.8])
comp_miss_after = unique(p_miss_with_neg_after2015$CompCode[p_miss_with_neg_after2015$Missing_Rate > 0.8])

# based on the result, "Cl-" was totally removed due to 100% missing before 2015
# "Accept.PM2.5" was removed due to largely missing from 2016
comp_miss_after = append("Cl-", comp_miss_after) 
comp_miss_after = append("Accept.PM2.5", comp_miss_after) 


# remove rows with largely missing
csn_before_2015[ ,comp_miss_before] <- list(NULL)
csn_after_2015[ ,comp_miss_after] <- list(NULL)

# reorder the files according to site and date
csn_before_2015 = csn_before_2015[with(csn_before_2015, 
                                       order(State, SiteCode, Date)), ]
csn_after_2015 = csn_after_2015[with(csn_after_2015, 
                                     order(State, SiteCode, Date)), ]

# write.csv(csn_before_2015, "CSN_Component_with_missing_Before_2015.csv")
# write.csv(csn_after_2015, "CSN_Component_with_missing_After_2015.csv")

# start the record of each site from the non-NA point
# cause some only start from very late
site.number = length(unique(csn_after_2015$SiteCode))

csn_before = NULL
csn_bef_site_info = NULL
for (i in 1:site.number){ 
  # extract data for a single site
  site.study.bef = unique(csn_before_2015$SiteCode)[i]
  site_single = subset(csn_before_2015, SiteCode == site.study.bef)
  
  # detect the first date when there is no NAs in PM component
  site_bef_1stDate = min(site_single$Date[!is.na(site_single$Al) & 
                                            !is.na(site_single$Cu) & 
                                            !is.na(site_single$SO4) & 
                                            !is.na(site_single$EC.TOR.unadjust.88) &
                                            !is.na(site_single$OC.TOR.unadjusted.88)])
  
  # in case of dataset with only NA
  if(is.na(site_bef_1stDate)){
    row.No = 0
    Date.last.record = NA
  } else{
    site_bef_with_NA = subset(site_single, Date >= site_bef_1stDate)
    row.No = nrow(site_bef_with_NA)
    csn_before = rbind(csn_before, site_bef_with_NA)
    Date.last.record = max(site_bef_with_NA$Date)
  }
    
  site_info = data.frame(SiteCode = site.study.bef, 
                         row.No = row.No,
                         Date.first.record = site_bef_1stDate,
                         Date.last.record = Date.last.record)
  csn_bef_site_info = rbind.fill(csn_bef_site_info, site_info[1, ])
}


csn_after = NULL
csn_aft_site_info = NULL
for (i in 1:site.number){ 
  # extract data for a single site
  site.study.aft = unique(csn_after_2015$SiteCode)[i]
  site_single = subset(csn_after_2015, SiteCode == site.study.aft)
  
  # detect the first date when there is no NAs in PM component
  site_aft_1stDate = min(site_single$Date[!is.na(site_single$Al) & 
                                            !is.na(site_single$Cu) & 
                                            !is.na(site_single$SO4) & 
                                            !is.na(site_single$EC.TOR.88) &
                                            !is.na(site_single$OC.88)])
  
  # in case of dataset with only NA
  if(is.na(site_aft_1stDate)){
    row.No = 0
    Date.last.record = NA
  } else{
    site_aft_with_NA = subset(site_single, Date >= site_aft_1stDate)
    row.No = nrow(site_aft_with_NA)
    csn_after = rbind(csn_after, site_aft_with_NA)
    Date.last.record = max(site_aft_with_NA$Date)
  }
  
  site_info = data.frame(SiteCode = site.study.aft, 
                         row.No = row.No,
                         Date.first.record = site_aft_1stDate,
                         Date.last.record = Date.last.record)
  csn_aft_site_info = rbind.fill(csn_aft_site_info, site_info[1, ])
}

# move state to the first column
csn_before = csn_before %>% relocate(State, .before = SiteCode)
csn_after = csn_after %>% relocate(State, .before = SiteCode)

## The biggest difference compared with previous one is this time, 
## unacceptable Qualifiers were marked Val as NA or 0.00009 before further process
write.csv(csn_before, "CSN_Component_with_missing_Before_2015_2023.02.csv")
write.csv(csn_after, "CSN_Component_with_missing_After_2015_2023.02.csv")

write.csv(csn_bef_site_info, "CSN_site_StartDate_info_Before_2015_2023.02.csv")
write.csv(csn_aft_site_info, "CSN_site_StartDate_info_After_2015_2023.02.csv")

#### extract data with qualifier 6 ####
csn_daily_before = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE/CSN_Component_with_missing_Before_2015.csv")
csn_daily_after = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE/CSN_Component_with_missing_After_2015.csv")
csn_daily_before$X = csn_daily_after$X =NULL

library(ggrepel)
csn_ql_after_MX = subset(csn_daily_after, 
                         grepl("MX", Qualifier, fixed=T))

csn_ql_before_1 = subset(csn_daily_before, 
                         grepl("1", Qualifier, fixed=T))

write.csv(csn_ql_before_1, "CSN_Qualifier_1_before2015.csv")
write.csv(csn_ql_after_MX, "CSN_Qualifier_MX_after2015.csv")

#### not used: CSN OC EC ####
csn_data_org = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original CSN/CSN method collection analysis use 10032022.csv") ## with extracted collection, analysis methods
# csn_data_org = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original CSN/CSN data for analysis 10232022.csv") ## with extracted collection, analysis methods
csn_data_org$X = csn_data_org$X.1 = csn_data_org$X.2 = NULL
csn_data_org$Date = as.Date(csn_data_org$Date, format="%m/%d/%Y")
head(csn_data_org)

csn_ocec = subset(csn_data_org, grepl("OC1", CompName, fixed = T) | 
                    grepl("EC1", CompName, fixed = T))
csn_oc_unadjust = subset(csn_ocec, CompName == "OC1.unadjusted.88")
csn_oc_unadj_urg_tor = subset(csn_oc_unadjust, Collection == "URG3000N" & Analysis == "IMPROVE_TOR")
csn_oc_unadj_urg_imA = subset(csn_oc_unadjust, Collection == "URG3000N" & Analysis == "IMPROVE_A")
min(csn_oc_unadj_urg_tor$Date); max(csn_oc_unadj_urg_tor$Date)
min(csn_oc_unadj_urg_imA$Date); max(csn_oc_unadj_urg_imA$Date)
length(unique(csn_oc_unadj_urg_tor$SiteCode))
length(unique(csn_oc_unadj_urg_imA$SiteCode))


################################################################################################
#### Filling the missing data - in another file ####
################################################################################################



