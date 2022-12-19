##clear environment
# rm(list=ls())

##set working directory
setwd("xxxx")
getwd()
data.dir <- "xxxx"

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

imp_data = read.csv("IMPROVE component only 10092022.csv")
imp_data$X = imp_data$X.1 = NULL
head(imp_data)
imp_data$Date = as.Date(imp_data$Date)
# imp_data$CompName[is.na(imp_data$CompName)] = "Na"
imp_data$Method = imp_data$MethodSimple
imp_data$class = 0
imp_data$class[imp_data$Method == "A-XRF"] = "Element"
imp_data$class[imp_data$Method == "B-IC"] = "Ion"
imp_data$class[imp_data$Method == "C-TOR"] = "OC/EC subgroup"

imp_data$Qualifier = imp_data$Status
imp_data_compare = select(imp_data, Dataset, State, SiteCode, Date, Method, class, CompName, Val, Qualifier, 
                          Unc, MDL, year, month)

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

csn_data$class[csn_data$class == "element.ion.PM"] = "Element"
csn_data$class[csn_data$CompName %in% c("Accept.PM2.5", "PM2.5RC")] = 0
csn_data$class[csn_data$CompName %in% c("Cl-", "K+", "NH4+", "Na+", "NO3", "SO4")] = "Ion"
csn_data$class[csn_data$class == "OC.EC"] = "OC/EC subgroup"
csn_data.2 = csn_data

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

#### overall component distribution check between nearest sites ####
imp_csn_nearest_site_data = read.csv("IMPROVE_CSN_Nearest_site_comparison.csv") 

element.group = unique(imp_data$CompName[imp_data$Method == "A-XRF"])
ion.group = unique(imp_data$CompName[imp_data$Method == "B-IC"])
oc.ec.group = unique(imp_data$CompName[imp_data$Method == "C-TOR"])

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

theme.comp = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                   strip.text.x = element_text(size = 14, colour = "grey25", angle = 0),
                   axis.title.x = element_text(color="grey25", size = 18, vjust=0, margin=margin(0,0,0,300)), 
                   axis.title.y = element_text(color="grey25", size = 18, vjust=1, margin=margin(0,2,0,0)),
                   axis.text.x = element_text(color="grey25", size = 16, angle = 90, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                   axis.text.y = element_text(color="grey25", size = 16, angle = 0, hjust = 0.5))

ggplot(subset(imp_csn_nearest_data, class != 0), 
       aes(IMPROVE, CSN, color = CompName)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ class, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp

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
