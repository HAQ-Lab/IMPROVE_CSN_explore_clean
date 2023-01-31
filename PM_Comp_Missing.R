##clear environment
# rm(list=ls())

##set working directory - for IMPROVE
setwd("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE")
getwd()
data.dir <- "Users/ztttttt/Documents/HEI PMF/R - original IMPROVE"

##packages in need
library(tidyr) # separate{tidyr}, gather{tidyr}, spread{tidyr},  spread is VIP function, str_split_fixed{stringr} is better than separate
library(stats) # aggregate{stats}, VIP function
library(ggplot2)
library(scales) # percent{}
library(stringr) # str_split_fixed, separate one column into multiple
library(dplyr)
library(plyr)
library(lubridate)
library(gridExtra) #grid.arrange{}
library(grid) #textGrob{}
library(imputeTS) #na_ma, na_interpolation ect.
library(mice) # using Markov Chain Monte Carlo simulation to impute the missing entries
library(tibble)
library(missForest) # implementation of random forest algorithm
library(ggsci)

####################################################################################
####### IMPROVE - Fill the Missing #######
####################################################################################
#### generate basic data for filling ####
imp_daily = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/IMPROVE_Component_with_missing.csv")
# imp_meta_para = read.csv("IMPROVE metadata 109 parameters 2010-20.csv")
imp_daily$X = imp_daily$X.1 = NULL
imp_daily$Date = as.Date(imp_daily$Date)
sapply(imp_daily, class)
head(imp_daily)
imp_daily = subset(imp_daily, Date > as.Date("2010-12-31"))
colnames(imp_daily)[37:48] # OPC sub-group
imp_daily_opt_sub = imp_daily[ ,37:48]
imp_miss = cbind(imp_daily[,1:36], imp_daily[,49:ncol(imp_daily)])
# remove those not directly detected, except of RC.PM2.5, ammoniaSO4, ammoniaNO3, OC & EC
imp_miss$SeaSalt = imp_miss$Soil = imp_miss$PM10 = imp_miss$RC.PM10 = imp_miss$site.date.qualifier = NULL
imp_miss$OC_UCD = imp_miss$EC_UCD = imp_miss$OMC = imp_miss$CM_calculated = imp_miss$TC =NULL
imp_miss$OPT = NULL # OPT is deleted cause it is not used to calculated OC & EC
head(imp_miss)
sum(-0.00310, 0.06140, 0.08480, 0.09960) + 0.13030

# imp_miss$DOW = weekdays(imp_miss$Date)
# imp_miss$year = year(imp_miss$Date)
# imp_miss$month = month(imp_miss$Date)
# imp_miss$week = week(imp_miss$Date)
# exclude states not located in the "Continental and mainland USA"
imp_miss = subset(imp_miss, !(State %in% c("AK", "HI", "PR", "VI"))) 
imp_miss = imp_miss[with(imp_miss, order(State, SiteCode, Date)), ]

## pattern of missing data exploration
# p_missing <- unlist(lapply(imp_miss, function(x) sum(is.na(x))))/nrow(imp_miss)
# sort(p_missing[p_missing > 0], decreasing = TRUE)
p_missing <- data.frame(unlist(lapply(imp_miss, function(x) sum(is.na(x))))/nrow(imp_miss))
colnames(p_missing) = "MissingRate"

n.site = length(unique(imp_miss$SiteCode))

theme.comp.1 = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                     strip.text.x = element_text(size = 14, colour = "grey25", angle = 0),
                     axis.title.x = element_text(color="grey25", size = 15, vjust=0, margin=margin(0,0,0,300)), 
                     axis.title.y = element_text(color="grey25", size = 15, vjust=1, margin=margin(0,2,0,0)),
                     axis.text.x = element_text(color="grey25", size = 14, angle = 0, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                     axis.text.y = element_text(color="grey25", size = 14, angle = 0, hjust = 0.5))

theme.comp.2 = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                     strip.text.x = element_text(size = 14, colour = "grey25", angle = 0),
                     axis.title.x = element_text(color="grey25", size = 15, vjust=0, margin=margin(0,0,0,300)), 
                     axis.title.y = element_text(color="grey25", size = 15, vjust=1, margin=margin(0,2,0,0)),
                     axis.text.x = element_text(color="grey25", size = 13, angle = 0, hjust = 0.5, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                     axis.text.y = element_text(color="grey25", size = 13, angle = 0, hjust = 0.5))

#### filling the NAs for each site ####
imp_avg_month_summary = imp_running_avg_summary = imp_linear_summary = imp_sgl_mice_summary = imp_sgl_rf_summary = NULL
interp_org_cmp_summary = NULL

for (i in 1 :n.site){ # i=85, 126~128,  all NA!
  site.study = unique(imp_miss$SiteCode)[i]
  imp_site_single = subset(imp_miss, SiteCode == site.study)
  imp_site_single_1stDate = min(imp_site_single$Date[!is.na(imp_site_single$Al) & 
                                                       !is.na(imp_site_single$Cu) & 
                                                       !is.na(imp_site_single$SO4) & 
                                                       !is.na(imp_site_single$OC1)])
  imp_site_single_use = subset(imp_site_single, Date >= imp_site_single_1stDate)
  row.No = nrow(imp_site_single_use)
  imp_site_single_use$row.No = 1:row.No
  
  ## check the temporal distribution and detect the date window to use for average
  ## choose MONTH as the time window
  # Al.date = ggplot(imp_site_single_use, aes(Date, Al)) + geom_point() + theme.comp.1
  # Al.year = ggplot(imp_site_single_use, aes(as.factor(year), Al)) + geom_boxplot() + theme.comp.2
  # Al.month = ggplot(imp_site_single_use, aes(as.factor(month), Al)) + geom_boxplot() + theme.comp.2
  # Al.week = ggplot(imp_site_single_use, aes(week, Al)) + geom_point() + theme.comp.2
  # Al.DOW = ggplot(imp_site_single_use, aes(DOW, Al)) + geom_boxplot() + theme.comp.2
  # multiplot(Al.year, Al.month, Al.week, Al.DOW, cols = 2)
  
  # SO4.date = ggplot(imp_site_single_use, aes(Date, SO4)) + geom_point() + theme.comp.1
  # SO4.year = ggplot(imp_site_single_use, aes(as.factor(year), SO4)) + geom_boxplot() + theme.comp.2
  # SO4.month = ggplot(imp_site_single_use, aes(as.factor(month), SO4)) + geom_boxplot() + theme.comp.2
  # SO4.week = ggplot(imp_site_single_use, aes(week, SO4)) + geom_point() + theme.comp.2
  # SO4.DOW = ggplot(imp_site_single_use, aes(DOW, SO4)) + geom_boxplot() + theme.comp.2
  # multiplot(SO4.year, SO4.month, SO4.week, SO4.DOW, cols = 2)
  
  # OC2.date = ggplot(subset(imp_site_single_use, OC2 < quantile(imp_site_single_use$OC2, 0.99, na.rm = T)), aes(Date, OC2)) + geom_point() + theme.comp.1
  # OC2.year = ggplot(subset(imp_site_single_use, OC2 < quantile(imp_site_single_use$OC2, 0.99, na.rm = T)), aes(as.factor(year), OC2)) + geom_boxplot() + theme.comp.2
  # OC2.month = ggplot(subset(imp_site_single_use, OC2 < quantile(imp_site_single_use$OC2, 0.99, na.rm = T)), aes(as.factor(month), OC2)) + geom_boxplot() + theme.comp.2
  # OC2.week = ggplot(subset(imp_site_single_use, OC2 < quantile(imp_site_single_use$OC2, 0.99, na.rm = T)), aes(week, OC2)) + geom_point() + theme.comp.2
  # OC2.DOW = ggplot(subset(imp_site_single_use, OC2 < quantile(imp_site_single_use$OC2, 0.99, na.rm = T)), aes(DOW, OC2)) + geom_boxplot() + theme.comp.2
  # multiplot(OC2.year, OC2.month, OC2.week, OC2.DOW, cols = 2)
  
  #### mix error calculation - data preparation
  # imp_site_single_use_all = na.omit(imp_site_single_use)
  # imp_site_single_use = prodNA(imp_site_single_use, noNA = 0.1)

  imp_sgl_avg_month = imp_sgl_running_avg = imp_sgl_linear = imp_site_single_use

  for(j in 5:(ncol(imp_site_single_use) - 5)){
    # interpolation1: linear  
    imp_sgl_linear[, j] = na_interpolation(imp_site_single_use[, j]) 
    # interpolation2: running window average 
    imp_sgl_running_avg[, j] = na_ma(imp_site_single_use[, j], weighting = "simple", k = 6) 
    
    for(l in 1:row.No){
      if(is.na(imp_site_single_use[l, j])){
        month.miss = imp_site_single_use$month[l]
        # interpolation3: average of the month of year
        imp_sgl_avg_month[l, j] = mean(subset(imp_site_single_use, month == month.miss)[, j], na.rm = T) 
      }
    }
  }
  
  # interpolation4: mice, using MCMC, multiple interpolation
  imp_sgl_miss_pattern = cbind(imp_site_single_use[, 8:34], imp_site_single_use[, 38:46])
  mice::md.pattern(imp_sgl_miss_pattern) # using mice to check the missing data pattern

  imp_sgl_mice <- mice(imp_site_single_use[, 5:46], maxit=0, m = 50, remove.collinear = F)
  imp_sgl_mice$loggedEvents
  # collinear exists in the dataset, mice will automatically remove NO3 & SO4 as a result
  # https://stefvanbuuren.name/fimd/sec-toomany.html
  
  imp_sgl_mice_dslist <- complete(imp_sgl_mice, "all")
  ## get the average of data.frame in the array
  imp_sgl_mice_avg = data.frame(aaply(laply(imp_sgl_mice_dslist, as.matrix),  c(2, 3), mean))
  imp_sgl_mice_avg_use = cbind(imp_site_single_use[, 1:4], imp_sgl_mice_avg, imp_site_single_use[, (ncol(imp_site_single_use)-4):ncol(imp_site_single_use)])
  
  #interpolation5: missForest, using random forest
  imp_sgl_rf_mf = missForest(imp_site_single_use[, 5:46])
  imp_sgl_rf = cbind(imp_site_single_use[, 1:4], data.frame(imp_sgl_rf_mf$ximp), imp_site_single_use[, 47:51])
   
  #### time series plotting for BIRM1
  # when taking i = 1, BIRM1 as example, plot time series:
  imp_linear_BIRM1 = imp_sgl_linear; imp_linear_BIRM1$Interporlation = "Linear"
  imp_running_BIRM1 = imp_sgl_running_avg; imp_running_BIRM1$Interporlation = "RunningAvg"
  imp_mice_BIRM1 = imp_sgl_mice_avg_use; imp_mice_BIRM1$Interporlation = "Multiple"
  imp_sf_BIRM1 = imp_sgl_rf; imp_sf_BIRM1$Interporlation = "RandomForest"
  imp_org_BIRM1 = imp_site_single_use; imp_org_BIRM1$Interporlation = "RawData"
  imp_mice_sf_BIRM1 = rbind(imp_org_BIRM1, imp_linear_BIRM1, imp_running_BIRM1, imp_mice_BIRM1, imp_sf_BIRM1)
  
  ggplot(imp_mice_sf_BIRM1, aes(Date, Al, col = Interporlation)) +
    geom_line(linetype="dashed") +
    geom_point(size = 1, alpha = 0.5) + 
    facet_grid(Interporlation~., scales = "free") + 
    scale_color_npg() + 
    theme_bw() +
    ylim(0,2)+
    theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
          strip.text.x = element_text(size = 20, colour = "Orange", angle = 0),
          legend.title = element_text(face="italic", size=16), # family="Times", colour="red", 
          legend.text = element_text(size=14), 
          axis.title.x = element_text(color="grey25", size = 20, vjust=0, margin=margin(0,0,0,300)), 
          axis.title.y = element_text(color="grey25", size = 20, vjust=1, margin=margin(0,2,0,0)),
          axis.text.x = element_text(color="grey25", size = 14, angle = 0, hjust = 0.5, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
          axis.text.y = element_text(color="grey25", size = 14, angle = 0, hjust = 0.5))
  
  
  # when taking i = 1, BIRM1 as example, the mixError from different method is:
  # mixError(imp_sgl_avg_month[, 5:46], imp_site_single_use[, 5:46], imp_site_single_use_all[, 5:46])
  # mixError(imp_sgl_running_avg[, 5:46], imp_site_single_use[, 5:46], imp_site_single_use_all[, 5:46]) # 0.4309416 
  # mixError(imp_sgl_linear[, 5:46], imp_site_single_use[, 5:46], imp_site_single_use_all[, 5:46]) # 0.4946119 
  # mixError(imp_sgl_mice_avg_use[, 5:46], imp_site_single_use[, 5:46], imp_site_single_use_all[, 5:46]) # 0.4798888 
  # mixError(imp_sgl_rf[, 5:46], imp_site_single_use[, 5:46], imp_site_single_use_all[, 5:46]) # 0.1040394 
  
  imp_avg_month_summary = rbind(imp_avg_month_summary, imp_sgl_avg_month)
  imp_running_avg_summary = rbind(imp_running_avg_summary, imp_sgl_running_avg)
  imp_linear_summary = rbind(imp_linear_summary, imp_sgl_linear)
  imp_sgl_mice_summary = rbind(imp_sgl_mice_summary, imp_sgl_mice_avg_use)
  imp_sgl_rf_summary = rbind(imp_sgl_rf_summary, imp_sgl_rf)
  
  single_use_median = sapply(imp_site_single_use[, 5:46], median, na.rm = T)
  single_use_mean = sapply(imp_site_single_use[, 5:46], mean, na.rm = T)
  single_use_sd = sapply(imp_site_single_use[, 5:46], sd, na.rm = T)
  avg_month_median = sapply(imp_sgl_avg_month[, 5:46], median)
  avg_month_mean = sapply(imp_sgl_avg_month[, 5:46], mean)
  avg_month_sd = sapply(imp_sgl_avg_month[, 5:46], sd)
  running_avg_median = sapply(imp_sgl_running_avg[, 5:46], median)
  running_avg_mean = sapply(imp_sgl_running_avg[, 5:46], mean)
  running_avg_sd = sapply(imp_sgl_running_avg[, 5:46], sd)
  linear_median = sapply(imp_sgl_linear[, 5:46], median)
  linear_mean = sapply(imp_sgl_linear[, 5:46], mean)
  linear_sd = sapply(imp_sgl_linear[, 5:46], sd)
  mice_avg_median = sapply(imp_sgl_mice_avg_use[, 5:46], median)
  mice_avg_mean = sapply(imp_sgl_mice_avg_use[, 5:46], mean)
  mice_avg_sd = sapply(imp_sgl_mice_avg_use[, 5:46], sd)
  rf_avg_median = sapply(imp_sgl_rf[, 5:46], median)
  rf_avg_mean = sapply(imp_sgl_rf[, 5:46], mean)
  rf_avg_sd = sapply(imp_sgl_rf[, 5:46], sd)
  
  interp_org_cmp = data.frame(rbind(single_use_median, avg_month_median, running_avg_median, linear_median, mice_avg_median, rf_avg_median, 
                         single_use_mean, avg_month_mean, running_avg_mean, linear_mean, mice_avg_mean,rf_avg_mean, 
                         single_use_sd, avg_month_sd, running_avg_sd, linear_sd, mice_avg_sd, rf_avg_sd))
  interp_org_cmp <- tibble::rownames_to_column(interp_org_cmp, "Comparisons")
  interp_org_cmp$State = imp_site_single_use$State[1]
  interp_org_cmp$SiteCode = imp_site_single_use$SiteCode[1]
  interp_org_cmp_summary = rbind(interp_org_cmp_summary, interp_org_cmp)
  # library(Metrics) # rmse{Metrics}
  # rmse(imp_site_single_use$Al, imp_sgl_avg_month$Al)
  # cor_avg_month = cor(imp_sgl_avg_month[, 5:46], imp_site_single_use[, 5:46], use = "complete.obs")
  
  ### figure out the functions to be used for mice interpolation
  # https://data.library.virginia.edu/getting-started-with-multiple-imputation-in-r/
  # https://datascienceplus.com/imputing-missing-data-with-r-mice-package/
  # sgl_test = cbind(imp_site_single_use[, 10:11], imp_site_single_use[, 26:27], imp_site_single_use[, 40:41])
  # sgl_test = sgl_test[1:100, ]
  # dim(sgl_test)
  # imp_sgl_mice_1 <- mice(sgl_test, maxit=0, m = 5)
  # xyplot(imp_sgl_mice_1,Ca~Cl+SO4+Si+OC,pch=18,cex=1)
  # densityplot(imp_sgl_mice_1)
  # head(complete(imp_sgl_mice_1, c(3, 5)))
  # imp_sgl_mice_1_dslist <- complete(imp_sgl_mice_1, "all")
  ## get the average of data.frame in the array
  # imp_sgl_mice_1_dslist_avg = data.frame(aaply(laply(imp_sgl_mice_1_dslist, as.matrix),  c(2, 3), mean))

  # imp_sgl_mice_2 <- mice(sgl_test, maxit=0, m = 50)
  # xyplot(imp_sgl_mice_2,Ca~Cl+SO4+Si+OC,pch=18,cex=1)
  # densityplot(imp_sgl_mice_2)
  # head(complete(imp_sgl_mice_2, c(3, 5)))
  # xyplot(imp_sgl_mice,Al~Fe+Si,pch=18,cex=1) # compare distribution of original & imputed data
  # densityplot(imp_sgl_mice) # magenta - density of the imputed data (each dataset); blue - the observed data
  # stripplot(imp_sgl_mice, pch = 20, cex = 1.2)
  
}

write.csv(imp_avg_month_summary, "IMPROVE_Impute_monthly-average.csv")
write.csv(imp_running_avg_summary, "IMPROVE_Impute_running-6-average.csv")
write.csv(imp_linear_summary, "IMPROVE_Impute_linear.csv")
write.csv(imp_sgl_mice_summary, "IMPROVE_Impute_multi-mice.csv")
write.csv(imp_sgl_rf_summary, "IMPROVE_Impute_random-forest.csv")
write.csv(interp_org_cmp_summary, "IMPROVE_Impute_method_comparison.csv")


#### 5 sites in MN & MI - prepare concentration & uncertainty data for PMF ####
imp_sgl_rf_summary = read.csv("IMPROVE_Impute_random-forest.csv")
imp_sgl_rf_MN_WI = subset(imp_sgl_rf_summary, State == "MN" | State == "WI")
imp_sgl_rf_MN_WI$X = imp_sgl_rf_MN_WI$Qualifier = imp_sgl_rf_MN_WI$DOW = imp_sgl_rf_MN_WI$year = 
  imp_sgl_rf_MN_WI$month = imp_sgl_rf_MN_WI$week = imp_sgl_rf_MN_WI$row.No = 
  imp_sgl_rf_MN_WI$ammNO3 = imp_sgl_rf_MN_WI$ammSO4 = imp_sgl_rf_MN_WI$PM2.5 = NULL

MN_WI_try = dplyr::select(imp_sgl_rf_MN_WI, SiteCode, Cl, Cr, Cu)[31:nrow(imp_sgl_rf_MN_WI), ]

imp_sgl_rf_MN_WI_conc = cbind(imp_sgl_rf_MN_WI[, 1:3], 
                              imp_sgl_rf_MN_WI %>% 
                                dplyr::select(where(is.numeric)) %>%
                                abs())
col_numeric <- which( sapply(imp_sgl_rf_MN_WI_conc, is.numeric ) )
MN_WI_conc_col_numeric = data.frame(
  sapply( col_numeric, function( y ) {
    quantile( x = unlist( imp_sgl_rf_MN_WI_conc[,  y ] ), 
              c(0.05, .1, 0.25, .5, .75, .8, 0.95),
              na.rm = TRUE )
  }))

imp_sgl_rf_MN_WI_conc_0.8 = MN_WI_conc_col_numeric[6, ]
imp_sgl_rf_MN_WI_conc = cbind(imp_sgl_rf_MN_WI_conc[, 1:3], 
                              imp_sgl_rf_MN_WI_conc %>% 
                                dplyr::select(where()) %>%
                                abs())
imp_sgl_rf_MN_WI_conc_1 = imp_sgl_rf_MN_WI_conc
for(i in 4:ncol(imp_sgl_rf_MN_WI_conc_1)){
  imp_sgl_rf_MN_WI_conc_1[ ,i] = ifelse(imp_sgl_rf_MN_WI_conc[ ,i]>0, 
                                      imp_sgl_rf_MN_WI_conc[ ,i], 
                                      imp_sgl_rf_MN_WI_conc_0.8[ , (i-3)]*0.05)
}
  
imp_sgl_rf_MN_WI_unc = cbind(imp_sgl_rf_MN_WI_conc_1[, 1:3], 
                             imp_sgl_rf_MN_WI_conc_1[, 4:ncol(imp_sgl_rf_MN_WI_conc_1)] * 0.1)
imp_sgl_rf_MN_WI_unc$RC.PM2.5 = imp_sgl_rf_MN_WI_unc$RC.PM2.5 * 30

### below code does not work
# library(purrr) # accumulate{purrr}
# imp_sgl_rf_MN_WI_conc = y <- join %>% 
#   group_by(SiteCode) %>% 
#   mutate(adj_ebit=purrr::accumulate(ebit,~ifelse(.y<0,(.y+.x)/2,.y)))

write.csv(imp_sgl_rf_MN_WI_conc_1, "interpolated_RF_MN_WI_conc.csv")
write.csv(imp_sgl_rf_MN_WI_unc, "interpolated_RF_MN_WI_unc.csv")
write.csv(imp_sgl_rf_MN_WI_conc_1[ ,4:ncol(imp_sgl_rf_MN_WI_conc_1)], "interpolated_RF_MN_WI_conc_use.csv")
write.csv(imp_sgl_rf_MN_WI_unc[ ,4:ncol(imp_sgl_rf_MN_WI_unc)], "interpolated_RF_MN_WI_unc_use.csv")

interp_org_cmp_summary = read.csv("IMPROVE_Impute_method_comparison.csv"); interp_org_cmp_summary$X = NULL

interp_org_cmp_summary$analysis = rep(c("single_use", "avg_month", "running_avg", 
                                        "linear", "mice_avg", "rf_avg"))
interp_org_cmp_summary$value = rep(c("median", "mean", "sd"), each = 6)

interp_org_cmp_summary_plot = interp_org_cmp_summary[, 1:(ncol(interp_org_cmp_summary)-3)]
interp_org_cmp_summary_plot = gather(interp_org_cmp_summary_plot, "CompName", "Concentration", -Comparisons)
interp_org_cmp_summary_plot$analysis = rep(c("single_use", "avg_month", "running_avg", 
                                             "linear", "mice_avg", "rf_avg"))
interp_org_cmp_summary_plot$value = rep(c("median", "mean", "sd"), each = 6)

theme.comp = theme(plot.title = element_text(hjust = 0.05, vjust = 0, size = 16),
                   strip.text.x = element_text(size = 14, colour = "grey25", angle = 0),
                   axis.title.x = element_text(color="grey25", size = 15, vjust=0, margin=margin(0,0,0,300)), 
                   axis.title.y = element_text(color="grey25", size = 15, vjust=1, margin=margin(0,2,0,0)),
                   axis.text.x = element_text(color="grey25", size = 14, angle = 0, hjust = 0, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                   axis.text.y = element_text(color="grey25", size = 14, angle = 0, hjust = 0.5))

ggplot(subset(interp_org_cmp_summary_plot, value == "median"), 
       aes(CompName, Concentration)) + 
  geom_boxplot(size = 1.5, alpha = 0.8) + 
  facet_grid(. ~ analysis, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  # geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp


ggplot(subset(interp_org_cmp_summary, value == "mean"), 
       aes(Al, Ca, color = analysis)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  # facet_grid(. ~ analysis, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  # geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
ggplot(subset(interp_org_cmp_summary, value == "mean"), 
       aes(NO3, SO4, color = analysis)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  # facet_grid(. ~ analysis, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  # geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
ggplot(subset(interp_org_cmp_summary, value == "mean"), 
       aes(OC, EC, color = analysis)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  # facet_grid(. ~ analysis, scales = "free") +  # facet_grid(.~Freq.Flag.Counts) +
  # geom_abline(slope=1, intercept=0, color = "black") +
  theme.comp
dev.off()

#### V.2022.11 - filling the NAs for each site - na.omit & check Mix Error ####
# n.site=85, 126~128,  all NA!, remove these sites
imp_miss_no_NA = subset(imp_miss, !(SiteCode %in% c("DETR1", "RENO1", "RENO2", "RENO3")))
n.site.no.NA = length(unique(imp_miss_no_NA$SiteCode)) # 169 sites to use for IMPROVE

mix_error_intp_summary = NULL

for (i in 1:n.site.no.NA){ 
  site.study = unique(imp_miss_no_NA$SiteCode)[i]
  imp_site_single = subset(imp_miss_no_NA, SiteCode == site.study)
  imp_site_single_1stDate = min(imp_site_single$Date[!is.na(imp_site_single$Al) & 
                                                       !is.na(imp_site_single$Cu) & 
                                                       !is.na(imp_site_single$SO4) & 
                                                       !is.na(imp_site_single$OC1)])
  imp_site_single_use = subset(imp_site_single, Date >= imp_site_single_1stDate)
  row.No = nrow(imp_site_single_use)
  
  # total NA number in components
  na.count = sum(is.na(imp_site_single_use[ ,5:46]))
  na.percent = na.count/(row.No * 42)
  
  # random set a dataframe with same percent of NAs
  imp_site_single_no_NA = na.omit(imp_site_single_use)
  imp_site_single_rdm_NA = cbind(imp_site_single_no_NA[ ,1:4], 
                                 prodNA(imp_site_single_no_NA[ ,5:46], noNA = na.percent),
                                 imp_site_single_no_NA[ ,47:51])
  
  
  ### interpolation
  imp_sgl_avg_month = imp_sgl_running_avg = imp_sgl_linear = imp_site_single_rdm_NA
  
  for(j in 5:(ncol(imp_site_single_rdm_NA) - 5)){
    # interpolation1: linear  
    imp_sgl_linear[, j] = na_interpolation(imp_site_single_rdm_NA[, j]) 
    # interpolation2: running window average 
    imp_sgl_running_avg[, j] = na_ma(imp_site_single_rdm_NA[, j], weighting = "simple", k = 6) 
    
#    for(l in 1:row.No){
#      if(is.na(imp_site_single_rdm_NA[l, j])){
#        month.miss = imp_site_single_rdm_NA$month[l]
#        # interpolation3: average of the month of year
#        imp_sgl_avg_month[l, j] = mean(subset(imp_site_single_rdm_NA, month == month.miss)[, j], na.rm = T) 
#      }
#    }
  }
  
  # interpolation4: mice, using MCMC, multiple interpolation
  imp_sgl_mice <- mice(imp_site_single_rdm_NA[, 5:46], maxit=0, m = 50, remove.collinear = F)
  # imp_sgl_mice$loggedEvents
  # collinear exists in the dataset, mice will automatically remove NO3 & SO4 as a result
  # https://stefvanbuuren.name/fimd/sec-toomany.html
  
  imp_sgl_mice_dslist <- complete(imp_sgl_mice, "all")
  ## get the average of data.frame in the array
  imp_sgl_mice_avg = data.frame(aaply(laply(imp_sgl_mice_dslist, as.matrix),  c(2, 3), mean))
  imp_sgl_mice_avg_use = cbind(imp_site_single_rdm_NA[, 1:4], imp_sgl_mice_avg, imp_site_single_rdm_NA[, (ncol(imp_site_single_rdm_NA)-4):ncol(imp_site_single_rdm_NA)])
  
  #interpolation5: missForest, using random forest
  imp_sgl_rf_mf = missForest(imp_site_single_rdm_NA[, 5:46])
  imp_sgl_rf = cbind(imp_site_single_rdm_NA[, 1:4], data.frame(imp_sgl_rf_mf$ximp), imp_site_single_rdm_NA[, 47:51])
  
  # me.avg.month = mixError(imp_sgl_avg_month[, 5:46], imp_site_single_rdm_NA[, 5:46], imp_site_single_no_NA[, 5:46])
  me.running.avg = mixError(imp_sgl_running_avg[, 5:46], imp_site_single_rdm_NA[, 5:46], imp_site_single_no_NA[, 5:46]) 
  me.linear = mixError(imp_sgl_linear[, 5:46], imp_site_single_rdm_NA[, 5:46], imp_site_single_no_NA[, 5:46]) 
  me.mice = mixError(imp_sgl_mice_avg_use[, 5:46], imp_site_single_rdm_NA[, 5:46], imp_site_single_no_NA[, 5:46]) 
  me.rf = mixError(imp_sgl_rf[, 5:46], imp_site_single_rdm_NA[, 5:46], imp_site_single_no_NA[, 5:46]) 
  
  mix_error_intp = data.frame(State = imp_site_single_rdm_NA$State[1], SiteCode = site.study, 
                              percent.NA = na.percent, row.total = row.No, row.no.NA = nrow(imp_site_single_rdm_NA), 
                              me.running.avg = me.running.avg, 
                              me.linear = me.linear, me.mice = me.mice, me.rf = me.rf)
  rownames(mix_error_intp)[1] = i
  mix_error_intp_summary = rbind(mix_error_intp_summary, mix_error_intp[1, ])
}

write.csv(mix_error_intp_summary, "IMPROVE_interpulation_Mix_Error.csv")

## plotting
colnames(mix_error_intp_summary)[6:9] = c("Running_Avg", "Linear", "Multiple", "Random_Forest")
mix_error_intp_plot = select(mix_error_intp_summary, SiteCode, Running_Avg, Linear, Multiple, Random_Forest)
mix_error_intp_plot = gather(mix_error_intp_plot, "Interpolations", "Mix_Error", -SiteCode)

theme.mix.err = theme(axis.title.y.right = element_blank(),
                      panel.spacing = unit(10, "mm"),   
                      legend.background = element_blank(),
                      strip.text = element_text(face="bold", size=rel(1.5)),
                      strip.background = element_rect(fill="lightblue", colour="grey", size=16),
                      axis.title.x = element_text(color="grey25", size = 24, vjust=-2, margin=margin(0,0,0,300)), 
                      axis.title.y = element_text(color="grey25", size = 24, vjust=2, margin=margin(0,2,0,0)),
                      plot.title=element_text(size=rel(2)), 
                      axis.text.x = element_text(color="grey25", size = 20, angle = 0, hjust = 0.5, vjust = 0.3), plot.margin = unit(c(2,1,2, 2), "lines"),
                      axis.text.y = element_text(color="grey25", size = 20, angle = 0, hjust = 0.5))

ggplot(mix_error_intp_plot, aes(Interpolations, Mix_Error, fill = Interpolations)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.5) +
  scale_fill_npg() + 
  theme_bw() +
  theme.mix.err

#### V.2022.12 - filling the NAs for each site (logged, no negative)- na.omit & check Mix Error ####
# n.site=85, 126~128,  all NA!, remove these sites
imp_miss_no_NA = subset(imp_miss, !(SiteCode %in% c("DETR1", "RENO1", "RENO2", "RENO3")))
n.site.no.NA = length(unique(imp_miss_no_NA$SiteCode)) # 169 sites to use for IMPROVE

running_avg_sum = linear_sum = mice_sum = rf_sum = NULL
mix_error_intp_pstv_summary = NULL

imp_var = data.frame(colnames(imp_miss_no_NA))
colnames(imp_var)[1] = "Variables"
p_neg_summary = p_miss_summary = imp_var

for (i in 1:n.site.no.NA){ 
  # extract data for a single site
  site.study = unique(imp_miss_no_NA$SiteCode)[i]
  site_single = subset(imp_miss_no_NA, SiteCode == site.study)
  
  # detect the first date when there is no NAs in PM component
  site_single_1stDate = min(site_single$Date[!is.na(site_single$Al) & 
                                               !is.na(site_single$Cu) & 
                                               !is.na(site_single$SO4) & 
                                               !is.na(site_single$OC1)])
  site_single_with_neg = subset(site_single, Date >= site_single_1stDate)
  row.No = nrow(site_single_with_neg)

  ## detect the percent of missing concentration values for each component
  p_miss_with_neg <- data.frame(unlist(lapply(site_single_with_neg, 
                                              function(x) 
                                                sum(is.na(x))))/row.No)
  colnames(p_miss_with_neg) = site.study
  
  ## detect the percent of negative concentration values for each component
  p_neg <- data.frame(unlist(lapply(site_single_with_neg, 
                                    function(x) 
                                      sum(x<0, na.rm = T)))/row.No)
  colnames(p_neg) = site.study
  
  ## summarise the missing / negative values for each site
  p_miss_summary = cbind(p_miss_summary, p_miss_with_neg)
  p_neg_summary = cbind(p_neg_summary, p_neg)

  ## remove the rows where all component concentrations are NAs
  col.withAllNa = ncol(site_single_with_neg)
  cols.comp = 5:col.withAllNa # columns for components
  col.component = ncol(site_single_with_neg[, cols.comp]) 
  
  site_single_neg_noAllNA = subset(site_single_with_neg, 
                                   rowSums(is.na(site_single_with_neg[, cols.comp])) != 
                                     col.component)
  
  ## substitute the negative or 0 with 0.000005 before interpolation
  ## 0.000009 was selected cause the present lowest positive value is 0.00001
  ## these value will be set to 1/2 MDL before PMF analysis 
  # question? change the pattern for those with lots of negatives??
  site_single_noAllNA = site_single_neg_noAllNA
  site_single_noAllNA[site_single_noAllNA <= 0] <- 0.000009
  
  ## log all value to avoid negative interpolation 
  site_single_log = cbind(site_single_noAllNA[, 1:4],
                          site_single_noAllNA %>% 
    dplyr::select(where(is.numeric)) %>%
    log())
  

  # total NA number in components
  na.count = sum(is.na(site_single_log[ ,cols.comp]))
  na.percent = na.count/(row.No * 42)
  
  # random set a dataframe with same percent of NAs as the original dataset
  site_single_no_NA = na.omit(site_single_log)
  site_single_rdm_NA = cbind(site_single_no_NA[ ,1:4], 
                             prodNA(site_single_no_NA[ ,cols.comp], 
                                    noNA = na.percent))
  
  
  #### interpolation, below is set for the randomly set NA, for method performance comparison
  sgl_intp_rdNA_running_avg = sgl_intp_rdNA_linear = site_single_rdm_NA

  for(j in cols.comp){
    # interpolation1: linear  
    sgl_intp_rdNA_linear[, j] = na_interpolation(site_single_rdm_NA[, j]) 
    # interpolation2: running window average 
    sgl_intp_rdNA_running_avg[, j] = na_ma(site_single_rdm_NA[, j], weighting = "simple", k = 6) 
  }
  
  # interpolation3: mice, using MCMC, multiple interpolation
  sgl_intp_rdNA_mice <- mice(site_single_rdm_NA[, cols.comp], maxit=0, m = 50, remove.collinear = F)
  # sgl_intp_rdNA_mice$loggedEvents
  # collinear exists in the dataset, mice will automatically remove NO3 & SO4 as a result
  # https://stefvanbuuren.name/fimd/sec-toomany.html
  
  sgl_intp_rdNA_mice_dslist <- complete(sgl_intp_rdNA_mice, "all")
  ## get the average of data.frame in the array
  sgl_intp_rdNA_mice_avg = data.frame(aaply(laply(sgl_intp_rdNA_mice_dslist, as.matrix),  c(2, 3), mean))
  sgl_intp_rdNA_mice_avg_use = cbind(site_single_rdm_NA[, 1:4], sgl_intp_rdNA_mice_avg)
  
  #interpolation4: missForest, using random forest
  sgl_intp_rdNA_rf_mf = missForest(site_single_rdm_NA[, cols.comp])
  sgl_intp_rdNA_rf = cbind(site_single_rdm_NA[, 1:4], data.frame(sgl_intp_rdNA_rf_mf$ximp))
  
  me.running.avg = mixError(sgl_intp_rdNA_running_avg[, cols.comp], site_single_rdm_NA[, cols.comp], site_single_no_NA[, cols.comp]) 
  me.linear = mixError(sgl_intp_rdNA_linear[, cols.comp], site_single_rdm_NA[, cols.comp], site_single_no_NA[, cols.comp]) 
  me.mice = mixError(sgl_intp_rdNA_mice_avg_use[, cols.comp], site_single_rdm_NA[, cols.comp], site_single_no_NA[, cols.comp]) 
  me.rf = mixError(sgl_intp_rdNA_rf[, cols.comp], site_single_rdm_NA[, cols.comp], site_single_no_NA[, cols.comp]) 
  
  mix_error_intp = data.frame(State = site_single_rdm_NA$State[1], SiteCode = site.study, 
                              percent.NA = na.percent, row.total = row.No, row.no.NA = nrow(site_single_rdm_NA), 
                              me.running.avg = me.running.avg, 
                              me.linear = me.linear, me.mice = me.mice, me.rf = me.rf)
  rownames(mix_error_intp)[1] = i
  mix_error_intp_pstv_summary = rbind(mix_error_intp_pstv_summary, mix_error_intp[1, ])
  
  
  ### interpolation, below is set for to get the interpolated data for future analysis
  sgl_intp_running_avg_pro = sgl_intp_linear_pro = site_single_log
  
  for(k in cols.comp){
    # interpolation1: linear  
    sgl_intp_linear_pro[, k] = na_interpolation(site_single_log[, k]) 
    # interpolation2: running window average 
    sgl_intp_running_avg_pro[, k] = na_ma(site_single_log[, k], weighting = "simple", k = 6) 
  }
  # covert logged concentrations to original
  sgl_intp_linear = cbind(sgl_intp_linear_pro[, 1:4], 
                          exp(sgl_intp_linear_pro[, cols.comp]))
  sgl_intp_running_avg = cbind(sgl_intp_running_avg_pro[, 1:4], 
                               exp(sgl_intp_running_avg_pro[, cols.comp]))
  
  # interpolation3: mice, using MCMC, multiple interpolation
  sgl_intp_mice_org <- mice(site_single_log[, cols.comp], maxit=0, m = 50, remove.collinear = F)
  sgl_intp_mice_dslist <- complete(sgl_intp_mice_org, "all")

  ## get the average of data.frame in the array
  sgl_intp_mice_avg = data.frame(aaply(laply(sgl_intp_mice_dslist, as.matrix),  c(2, 3), mean))
  sgl_intp_mice = cbind(site_single_log[, 1:4], 
                        exp(sgl_intp_mice_avg))
  
  #interpolation4: missForest, using random forest
  sgl_intp_rf_mf = missForest(site_single_log[, cols.comp])
  sgl_intp_rf = cbind(site_single_log[, 1:4], 
                      exp(sgl_intp_rf_mf$ximp))
  
  running_avg_sum = rbind(running_avg_sum, sgl_intp_running_avg)
  linear_sum = rbind(linear_sum, sgl_intp_linear)
  mice_sum = rbind(mice_sum, sgl_intp_mice)
  rf_sum = rbind(rf_sum, sgl_intp_rf)
}

rownames(p_miss_summary) = 1:nrow(p_miss_summary)
rownames(p_neg_summary) = 1:nrow(p_neg_summary)
round(rowMeans(p_neg_summary[2:(ncol(p_neg_summary)-1)], ), 3)
write.csv(p_miss_summary, "IMPROVE_Missing_Rate_Site-level.csv")
write.csv(p_neg_summary, "IMPROVE_Negative_Rate_Site-level.csv")

write.csv(mix_error_intp_pstv_summary, "IMPROVE_interpulation_Mix_Error_positive.csv")

write.csv(running_avg_sum, "IMPROVE_interpulation_running-6-average_afterLog.csv")
write.csv(linear_sum, "IMPROVE_interpulation_linear_afterLog.csv")
write.csv(mice_sum, "IMPROVE_interpulation_multi-mice_afterLog.csv")
write.csv(rf_sum, "IMPROVE_interpulation_random-forest_afterLog.csv")

### to find out how was mixError calculated, the dataset was generated from the loop above
rf_intp_example = sgl_intp_rdNA_rf[, cols.comp]
rdm_NA_example = site_single_rdm_NA[, cols.comp]
site_noNa_example = site_single_no_NA[, cols.comp]

mixError(rf_intp_example, rdm_NA_example, site_noNa_example) 
varClass(rf_intp_example)

mis <- is.na(rdm_NA_example)
sqrt(mean((rf_intp_example[mis] - site_noNa_example[mis])^{2}) / stats::var(site_noNa_example[mis]))

# nrmse, a function in package missForest, and internally used by mixError{missForest} for numeric variables 
nrmse <- function(ximp, xmis, xtrue){
  mis <- is.na(xmis)
  sqrt(mean((ximp[mis] - xtrue[mis])^{2}) / stats::var(xtrue[mis]))
}

nrmse(rf_intp_example, rdm_NA_example, site_noNa_example)

## plotting
colnames(mix_error_intp_pstv_summary)[6:9] = c("Running_Avg", "Linear", "Multiple", "Random_Forest")
mix_error_intp_pstv_plot = select(mix_error_intp_pstv_summary, SiteCode, Running_Avg, Linear, Multiple, Random_Forest)
mix_error_intp_pstv_plot = gather(mix_error_intp_pstv_plot, "Interpolations", "Mix_Error", -SiteCode)

theme.mix.err = theme(axis.title.y.right = element_blank(),
                      panel.spacing = unit(10, "mm"),   
                      legend.background = element_blank(),
                      strip.text = element_text(face="bold", size=rel(1.5)),
                      strip.background = element_rect(fill="lightblue", colour="grey", size=16),
                      axis.title.x = element_text(color="grey25", size = 24, vjust=-2, margin=margin(0,0,0,300)), 
                      axis.title.y = element_text(color="grey25", size = 24, vjust=2, margin=margin(0,2,0,0)),
                      plot.title=element_text(size=rel(2)), 
                      axis.text.x = element_text(color="grey25", size = 20, angle = 15, hjust = 0.5, vjust = 0.5), plot.margin = unit(c(2,1,2, 2), "lines"),
                      axis.text.y = element_text(color="grey25", size = 20, angle = 0, hjust = 0.5))

ggplot(mix_error_intp_pstv_plot, aes(Interpolations, Mix_Error, fill = Interpolations)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.5) +
  scale_fill_npg() + 
  theme_bw() +
  theme.mix.err

### heat map for the negative and missing rate distribution for each PM species across sampling sites
## variables not to be included in the heat plot
imp_exclude = c("SiteCode", "Date", "Qualifier", "State", "ammNO3", "ammSO4", "NO2.", "RC.PM2.5")

## transfer dataset for plotting
p_miss_summary_plot = gather(p_miss_summary, "SiteCode", "Missing_Rate", -Variables)
p_miss_summary_plot = subset(p_miss_summary_plot, !(Variables %in% imp_exclude))
colnames(p_miss_summary_plot)[1] = "PM2.5_Species"
p_miss_summary_plot$PM2.5_Species[p_miss_summary_plot$PM2.5_Species == "Cl."] = "Chl"

p_neg_summary_plot = gather(p_neg_summary, "SiteCode", "Negative_Rate", -Variables)
p_neg_summary_plot = subset(p_neg_summary_plot, !(Variables %in% imp_exclude))
colnames(p_neg_summary_plot)[1] = "PM2.5_Species"
p_neg_summary_plot$PM2.5_Species[p_neg_summary_plot$PM2.5_Species == "Cl."] = "Chl"

theme.miss.neg = theme(axis.title.y.right = element_blank(),
                       panel.spacing = unit(10, "mm"),   
                       legend.background = element_blank(),
                       strip.text = element_text(face="bold", size=rel(1.5)),
                       strip.background = element_rect(fill="lightblue", colour="grey", size=16),
                       axis.title.x = element_text(color="grey25", size = 24, vjust=-2, margin=margin(0,0,0,300)), 
                       axis.title.y = element_text(color="grey25", size = 24, vjust=2, margin=margin(0,2,0,0)),
                       plot.title=element_text(size=rel(2)), 
                       axis.text.x = element_text(color="grey25", size = 10, angle = 15, hjust = 0.5, vjust = 0.5), plot.margin = unit(c(2,1,2, 2), "lines"),
                       axis.text.y = element_text(color="grey25", size = 2, angle = 0, hjust = 0.5))

## create a heat map
ggplot(p_miss_summary_plot, 
       aes(x=PM2.5_Species, y=SiteCode, fill=Missing_Rate)) +
  geom_tile() +
  scale_fill_material("teal") +
  theme.miss.neg
  
ggplot(p_neg_summary_plot, 
       aes(x=PM2.5_Species, y=SiteCode, fill=Negative_Rate)) +
  geom_tile() +
  scale_fill_material("deep-orange") +
  theme.miss.neg


####################################################################################
##### CSN - Fill the Missing #####
####################################################################################
#### generate basic data for filling ####
csn_daily = read.csv("CSN_Component_with_missing_Before_2015.csv")
csn_daily = read.csv("CSN_Component_with_missing_After_2015.csv")
csn_daily = read.csv("/Users/ztttttt/Documents/HEI PMF/R - original IMPROVE/CSN_Component_with_missing_After_2015.csv")

csn_daily$X = csn_daily$X.1 = NULL

# convert data.frame to data.table 
# setDT(csn_daily) # cancel this step, multiple functions do not work with data.table

csn_daily$Date = as.Date(csn_daily$Date)
sapply(csn_daily, class)
head(csn_daily)

site_day_NA_count = ddply(csn_daily, 
                          .(SiteCode), 
                          summarise,
                          count = length(Date),
                          Mg.NA = sum(is.na(Mg)),
                          EC.unadjusted.88.NA = sum(is.na(EC.unadjusted.88)),
                          OC.unadjusted.88.NA = sum(is.na(OC.unadjusted.88)),
                          SO4.NA = sum(is.na(SO4)))
site_day_NA_count$Mg.NA.per = round(site_day_NA_count$Mg.NA/site_day_NA_count$count, 3)
site_day_NA_count$EC.NA.per = round(site_day_NA_count$EC.unadjusted.88.NA/site_day_NA_count$count, 3)
site_day_NA_count$OC.NA.per = round(site_day_NA_count$OC.unadjusted.88.NA/site_day_NA_count$count, 3)
site_day_NA_count$SO4.NA.per = round(site_day_NA_count$SO4.NA/site_day_NA_count$count, 3)

# detect two sides with < 100 rows (before 2015)
data_132950002 = subset(csn_daily, SiteCode == 132950002)
data_530530031 = subset(csn_daily, SiteCode == 530530031)

# detect three sides with < 100 rows (after 2015)
data_11130001 = subset(csn_daily, SiteCode == 11130001)
data_60731018 = subset(csn_daily, SiteCode == 60731018)
data_220150008 = subset(csn_daily, SiteCode == 220150008)
# lack Cl- or mostly 9E-06 for EC3

# remove the sites with only NAs
csn_miss = subset(csn_daily, 
                  !(SiteCode %in% 
                      site_day_NA_count$SiteCode[
                        site_day_NA_count$Mg.NA.per > 0.5 & 
                          site_day_NA_count$EC.NA.per > 0.5 & 
                          site_day_NA_count$OC.NA.per > 0.5]))
# 135 out of 136 sites left for data before 2015
# 146 out of 149 sites left for data after 2015

csn_miss = csn_miss[with(csn_miss, order(State, SiteCode, Date)), ]
n.site = length(unique(csn_miss$SiteCode))

# numeric will be calculated later, which we wants to avoid
csn_miss$SiteCode = as.character(csn_miss$SiteCode)

#### V.2022.12 - filling the NAs for each site (logged, no negative)- na.omit & check Mix Error ####

# create data.frame to store results
running_avg_sum = linear_sum = mice_sum = rf_sum = NULL
mix_error_intp_pstv_summary = NULL

# create data.frame with the first column being the colnames of csn_miss for matching data
csn_var = data.frame(colnames(csn_miss))
colnames(csn_var)[1] = "Variables"
p_neg_summary = p_miss_summary = csn_var

for (i in 1:n.site){ 
  # 12, 20, 89, 100(?), 101, 103, 125-6, 129 (before_2015 data)
  # 4-5, 14, 18, 37, 49, 71,75-6, 111-2, 118, 144 (after_2015 data)
  # extract data for a single site
  site.study = unique(csn_miss$SiteCode)[i]
  site_single = subset(csn_miss, SiteCode == site.study)
  
  # detect the first date when there is no NAs in PM component
  site_single_1stDate = min(site_single$Date[!is.na(site_single$Al) & 
                                               !is.na(site_single$Cu) & 
                                               !is.na(site_single$SO4) & 
                                               # !is.na(site_single$OC1)] # (before 2015)
                                               !is.na(site_single$EC1)]) # (after 2015)
  site_single_with_NA = subset(site_single, Date >= site_single_1stDate)
  row.No = nrow(site_single_with_NA)
  
  ## detect the percent of missing concentration values for each component
  p_miss_with_neg <- data.frame(unlist(lapply(site_single_with_NA, 
                                              function(x) 
                                                sum(is.na(x))))/row.No)
  colnames(p_miss_with_neg) = site.study
  
  ## detect the percent of negative concentration values for each component
  p_neg <- data.frame(unlist(lapply(site_single_with_NA, 
                                    function(x) 
                                      sum(x<0, na.rm = T)))/row.No)
  colnames(p_neg) = site.study
  
  ## summarise the missing / negative values for each site
  p_miss_summary = cbind(p_miss_summary, p_miss_with_neg)
  p_neg_summary = cbind(p_neg_summary, p_neg)
  
  ## remove the rows where all component concentrations are NAs
  col.withAllNa = ncol(site_single_with_NA)
  cols.comp = 5:col.withAllNa # columns for PM/components
  col.component = ncol(site_single_with_NA[, cols.comp]) # the code does not work for data.table
  
  site_single_noAllNA = subset(site_single_with_NA, 
                                   rowSums(is.na(site_single_with_NA[, cols.comp])) != 
                                     col.component)
  
  ## substitute the negative or 0 with 0.000005 before interpolation
  ## 0.000009 was selected cause the present lowest positive value is 0.00001
  ## these value will be set to 1/2 MDL before PMF analysis 
  # question? change the pattern for those with lots of negatives??
  # this step has been carried out for CSN data in early steps
  
  ## log all value to avoid negative interpolation 
  site_single_log = cbind(site_single_noAllNA[, 1:4],
                          site_single_noAllNA %>% 
                            dplyr::select(where(is.numeric)) %>%
                            log())
  
  
  # total NA number in components
  na.count = sum(is.na(site_single_log[ ,cols.comp]))
  na.percent = na.count/(row.No * col.component)
  
  # random set a dataframe with same percent of NAs as the original dataset
  site_single_no_NA = na.omit(site_single_log)
  site_single_rdm_NA = cbind(site_single_no_NA[ ,1:4], 
                             prodNA(site_single_no_NA[ ,cols.comp], 
                                    noNA = na.percent))
  
  
  #### interpolation, below is set for the randomly set NA, for method performance comparison
  sgl_intp_rdNA_running_avg = sgl_intp_rdNA_linear = site_single_rdm_NA
  
  for(j in cols.comp){
    # interpolation1: linear  
    sgl_intp_rdNA_linear[, j] = na_interpolation(site_single_rdm_NA[, j]) 
    # interpolation2: running window average 
    sgl_intp_rdNA_running_avg[, j] = na_ma(site_single_rdm_NA[, j], 
                                           weighting = "simple", k = 6) 
  }
  
  # interpolation3: mice, using MCMC, multiple interpolation
  sgl_intp_rdNA_mice <- mice(site_single_rdm_NA[, cols.comp], 
                             maxit=0, m = 50,
                             remove.collinear = F)
  # sgl_intp_rdNA_mice$loggedEvents
  # collinear exists in the dataset, mice will automatically remove NO3 & SO4 as a result
  # https://stefvanbuuren.name/fimd/sec-toomany.html
  
  sgl_intp_rdNA_mice_dslist <- complete(sgl_intp_rdNA_mice, "all")
  ## get the average of data.frame in the array
  sgl_intp_rdNA_mice_avg = data.frame(aaply(
    laply(
      sgl_intp_rdNA_mice_dslist, as.matrix),  
    c(2, 3), 
    mean))
  sgl_intp_rdNA_mice_avg_use = cbind(site_single_rdm_NA[, 1:4], 
                                     sgl_intp_rdNA_mice_avg)
  
  #interpolation4: missForest, using random forest
  sgl_intp_rdNA_rf_mf = missForest(site_single_rdm_NA[, cols.comp])
  sgl_intp_rdNA_rf = cbind(site_single_rdm_NA[, 1:4], 
                           data.frame(sgl_intp_rdNA_rf_mf$ximp))
  
  me.running.avg = mixError(sgl_intp_rdNA_running_avg[, cols.comp], 
                            site_single_rdm_NA[, cols.comp], 
                            site_single_no_NA[, cols.comp]) 
  me.linear = mixError(sgl_intp_rdNA_linear[, cols.comp], 
                       site_single_rdm_NA[, cols.comp], 
                       site_single_no_NA[, cols.comp]) 
  me.mice = mixError(sgl_intp_rdNA_mice_avg_use[, cols.comp], 
                     site_single_rdm_NA[, cols.comp], 
                     site_single_no_NA[, cols.comp]) 
  me.rf = mixError(sgl_intp_rdNA_rf[, cols.comp], 
                   site_single_rdm_NA[, cols.comp], 
                   site_single_no_NA[, cols.comp]) 
  
  mix_error_intp = data.frame(State = site_single_rdm_NA$State[1], 
                              SiteCode = site.study, 
                              percent.NA = na.percent, 
                              row.total = row.No, 
                              row.no.NA = nrow(site_single_rdm_NA), 
                              me.running.avg = me.running.avg, 
                              me.linear = me.linear, 
                              me.mice = me.mice, 
                              me.rf = me.rf)
  rownames(mix_error_intp)[1] = i
  mix_error_intp_pstv_summary = rbind(mix_error_intp_pstv_summary, 
                                      mix_error_intp[1, ])
  
  
  ### interpolation, below is set for to get the interpolated data for future analysis
  sgl_intp_running_avg_pro = sgl_intp_linear_pro = site_single_log
  
  for(k in cols.comp){
    # interpolation1: linear  
    sgl_intp_linear_pro[, k] = na_interpolation(site_single_log[, k]) 
    # interpolation2: running window average 
    sgl_intp_running_avg_pro[, k] = na_ma(site_single_log[, k], 
                                          weighting = "simple", k = 6) 
  }
  # covert logged concentrations to original
  sgl_intp_linear = cbind(sgl_intp_linear_pro[, 1:4], 
                          exp(sgl_intp_linear_pro[, cols.comp]))
  sgl_intp_running_avg = cbind(sgl_intp_running_avg_pro[, 1:4], 
                               exp(sgl_intp_running_avg_pro[, cols.comp]))
  
  # interpolation3: mice, using MCMC, multiple interpolation
  sgl_intp_mice_org <- mice(site_single_log[, cols.comp], 
                            maxit=0, m = 50, 
                            remove.collinear = F)
  sgl_intp_mice_dslist <- complete(sgl_intp_mice_org, "all")
  
  ## get the average of data.frame in the array
  sgl_intp_mice_avg = data.frame(aaply(
    laply(sgl_intp_mice_dslist, as.matrix),  
    c(2, 3), mean))
  sgl_intp_mice = cbind(site_single_log[, 1:4], 
                        exp(sgl_intp_mice_avg))
  
  #interpolation4: missForest, using random forest
  sgl_intp_rf_mf = missForest(site_single_log[, cols.comp])
  sgl_intp_rf = cbind(site_single_log[, 1:4], 
                      exp(sgl_intp_rf_mf$ximp))
  
  running_avg_sum = rbind(running_avg_sum, sgl_intp_running_avg)
  linear_sum = rbind(linear_sum, sgl_intp_linear)
  mice_sum = rbind(mice_sum, sgl_intp_mice)
  rf_sum = rbind(rf_sum, sgl_intp_rf)
}

rownames(p_miss_summary) = 1:nrow(p_miss_summary)
rownames(p_neg_summary) = 1:nrow(p_neg_summary)
round(rowMeans(p_neg_summary[2:(ncol(p_neg_summary)-1)], ), 3)

# out put results for data before 2015
write.csv(p_miss_summary, "CSN_Missing_Rate_Site-level_before_2015.csv")
write.csv(p_neg_summary, "CSN_Negative_Rate_Site-level_before_2015.csv")
write.csv(mix_error_intp_pstv_summary, "CSN_interpulation_Mix_Error_positive_before_2015.csv")

write.csv(running_avg_sum, "CSN_interpulation_running-6-average_afterLog_before_2015.csv")
write.csv(linear_sum, "CSN_interpulation_linear_afterLog_before_2015.csv")
write.csv(mice_sum, "CSN_interpulation_multi-mice_afterLog_before_2015.csv")
write.csv(rf_sum, "CSN_interpulation_random-forest_afterLog_before_2015.csv")

# out put results for data after 2015
write.csv(p_miss_summary, "CSN_Missing_Rate_Site-level_after_2015.csv")
write.csv(p_neg_summary, "CSN_Negative_Rate_Site-level_after_2015.csv")
write.csv(mix_error_intp_pstv_summary, "CSN_interpulation_Mix_Error_positive_after_2015.csv")

write.csv(running_avg_sum, "CSN_interpulation_running-6-average_afterLog_after_2015.csv")
write.csv(linear_sum, "CSN_interpulation_linear_afterLog_after_2015.csv")
write.csv(mice_sum, "CSN_interpulation_multi-mice_afterLog_after_2015.csv")
write.csv(rf_sum, "CSN_interpulation_random-forest_afterLog_after_2015.csv")

### to find out how was mixError calculated, the dataset was generated from the loop above
# based on IMPROVE data
rf_intp_example = sgl_intp_rdNA_rf[, 5:col.withAllNa]
rdm_NA_example = site_single_rdm_NA[, 5:col.withAllNa]
site_noNa_example = site_single_no_NA[, 5:col.withAllNa]

mixError(rf_intp_example, rdm_NA_example, site_noNa_example) 
varClass(rf_intp_example)

mis <- is.na(rdm_NA_example)
sqrt(mean((rf_intp_example[mis] - site_noNa_example[mis])^{2}) / stats::var(site_noNa_example[mis]))

# nrmse, a function in package missForest, and internally used by mixError{missForest} for numeric variables 
nrmse <- function(xcsn, xmis, xtrue){
  mis <- is.na(xmis)
  sqrt(mean((xcsn[mis] - xtrue[mis])^{2}) / stats::var(xtrue[mis]))
}

nrmse(rf_intp_example, rdm_NA_example, site_noNa_example)

## plotting
colnames(mix_error_intp_pstv_summary)[6:9] = c("Running_Avg", "Linear", "Multiple", "Random_Forest")
mix_error_intp_pstv_plot = select(mix_error_intp_pstv_summary, SiteCode, Running_Avg, Linear, Multiple, Random_Forest)
mix_error_intp_pstv_plot = gather(mix_error_intp_pstv_plot, "Interpolations", "Mix_Error", -SiteCode)

theme.mix.err = theme(axis.title.y.right = element_blank(),
                      panel.spacing = unit(10, "mm"),   
                      legend.background = element_blank(),
                      strip.text = element_text(face="bold", size=rel(1.5)),
                      strip.background = element_rect(fill="lightblue", colour="grey", size=16),
                      axis.title.x = element_text(color="grey25", size = 24, vjust=-2, margin=margin(0,0,0,300)), 
                      axis.title.y = element_text(color="grey25", size = 24, vjust=2, margin=margin(0,2,0,0)),
                      plot.title=element_text(size=rel(2)), 
                      axis.text.x = element_text(color="grey25", size = 20, angle = 15, hjust = 0.5, vjust = 0.5), plot.margin = unit(c(2,1,2, 2), "lines"),
                      axis.text.y = element_text(color="grey25", size = 20, angle = 0, hjust = 0.5))

ggplot(mix_error_intp_pstv_plot, aes(Interpolations, Mix_Error, fill = Interpolations)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.5) +
  scale_fill_npg() + 
  theme_bw() +
  theme.mix.err


####################################################################################
#### Test - filling the NAs for whole dataset ####
imp_miss_pattern = cbind(imp_miss[, 8:14],imp_miss[, 16:24], imp_miss[, 26:34], imp_miss[, 38:45])
imp_miss_pattern$OC = NULL
mice::md.pattern(imp_miss_pattern) # using mice to check the missing data pattern

imp_miss_whl = na.omit(imp_miss)
row.No.whl = nrow(imp_miss_whl)
# missing value number from main pattern, 3033*23+10199*26+7117*10+2757*13+9925*39+2275+42 = 831336
# 831336/(nrow(imp_miss)*42) = 0.10208941
imp_rdm_na = prodNA(imp_miss_whl, noNA = 0.1020894)
imp_rdm_na = prodNA(imp_miss_whl, noNA = 0.1020894)

imp_avg_month = imp_running_avg = imp_linear = imp_rdm_na

for(j in 5:(ncol(imp_rdm_na) - 4)){
  # interpolation1: linear  
  imp_linear[, j] = na_interpolation(imp_rdm_na[, j]) 
  # interpolation2: running window average 
  imp_running_avg[, j] = na_ma(imp_rdm_na[, j], weighting = "simple", k = 6) 
  
  for(l in 1:row.No.whl){
    if(is.na(imp_rdm_na[l, j])){
      month.miss.whl = imp_rdm_na$month[l]
      # interpolation3: average of the month of year
      imp_avg_month[l, j] = mean(subset(imp_rdm_na, month == month.miss.whl)[, j], na.rm = T) 
    }
  }
}

# interpolation4: mice, using MCMC, multiple interpolation
imp_mice <- mice(imp_rdm_na[, 5:46], maxit=0, m = 50, remove.collinear = F)
imp_mice$loggedEvents
# collinear exists in the dataset, mice will automatically remove NO3 & SO4 as a result
# https://stefvanbuuren.name/fimd/sec-toomany.html

imp_mice_dslist <- complete(imp_mice, "all")
## get the average of data.frame in the array
imp_mice_avg = data.frame(aaply(laply(imp_mice_dslist, as.matrix),  c(2, 3), mean))
imp_mice_avg_use = cbind(imp_rdm_na[, 1:4], imp_mice_avg, imp_rdm_na[, (ncol(imp_rdm_na)-4):ncol(imp_rdm_na)])

#interpolation5: missForest, using random forest
imp_rf_mf = missForest(imp_rdm_na[, 5:46])
imp_rf = cbind(imp_rdm_na[, 1:4], data.frame(imp_rf_mf$ximp), imp_rdm_na[, 47:51])

# mixError, the normalized root mean squared error for the continuous and the proportion of falsely classified entries for the categorical variables are computed.
mixError(imp_avg_month[, 5:46], imp_rdm_na[, 5:46], imp_miss_whl[, 5:46])
mixError(imp_running_avg[, 5:46], imp_rdm_na[, 5:46], imp_miss_whl[, 5:46]) # 0.5892714
mixError(imp_linear[, 5:46], imp_rdm_na[, 5:46], imp_miss_whl[, 5:46])  # 0.6277518 
mixError(imp_mice_avg_use[, 5:46], imp_rdm_na[, 5:46], imp_miss_whl[, 5:46])  
mixError(imp_rf[, 5:46], imp_rdm_na[, 5:46], imp_miss_whl[, 5:46])  

## not need to carry out the bottom calculation, no need to further normalize the dataset
mixError(scale(imp_avg_month[, 5:46]), scale(imp_rdm_na[, 5:46]), scale(imp_miss_whl[, 5:46]))
mixError(scale(imp_running_avg[, 5:46]), scale(imp_rdm_na[, 5:46]), scale(imp_miss_whl[, 5:46])) # 0.8319173
mixError(scale(imp_linear[, 5:46]), scale(imp_rdm_na[, 5:46]), scale(imp_miss_whl[, 5:46]))  # 0.9494074 
mixError(scale(imp_mice_avg_use[, 5:46]), scale(imp_rdm_na[, 5:46]), scale(imp_miss_whl[, 5:46]))  
mixError(scale(imp_rf[, 5:46]), scale(imp_rdm_na[, 5:46]), scale(imp_miss_whl[, 5:46]))  

####################################################################################
######## methods for interpolation ########
####################################################################################
#### imputeTS, different interpolation ####
library(imputeTS)

#Prerequisite: Create Time series with missing values
x <- ts(c(12,33,NA,25,6,NA,7,48,NA))
y = x
y = na_seadec(y, algorithm = "interpolation")

#Example 1: Replace all NAs by random values that are between min and max of the input time series

a = na_random(x)
b = na_locf(x, option = "locf")   # Last Obs. Carried Forward
c = na_locf(x, option = "nocb")   # Next Obs. Carried Backward
d = na_interpolation(x)           # Linear Interpolation
e = na_seadec(x, algorithm = "interpolation") # Seasonal Adjustment then Linear Interpolation
# "interpolation" - Imputation by Interpolation (default choice)
# "locf" - Imputation by Last Observation Carried Forward
# "mean" - Imputation by Mean Value
# "random" - Imputation by Random Sample
# "kalman" - Imputation by Kalman Smoothing and State Space Models
# "ma" - Imputation by Weighted Moving Average

f = na_mean(x, option = "mean")   # Mean Imputation
g = na_mean(x, option = "median") # Median Imputation
h = na_mean(x, option = "mode")   # Mode Imputation

# Replace all NAs by random values between 1 and 10
na_random(x, lowerBound = 1, upperBound = 10)


#### mice, using MCMC, multiple interpolation ####
mydata = data.frame(x = c(1,2,8,3,7,4,6,5,10),
            y = c(43,NA,28,NA,19,100,67, 83,NA),
            z = c(17,33,NA,6,9,NA,NA, 34,NA))

p_missing <- unlist(lapply(mydata, function(x) sum(is.na(x))))/nrow(mydata)
sort(p_missing[p_missing > 0], decreasing = TRUE)

# Deterministic regression imputation via mice
imp <- mice(mydata, method = "norm.predict", m = 1)
# Store data
data_imp <- complete(imp)
# Multiple Imputation
imp <- mice(mydata, m = 5)
#build predictive model
fit <- with(data = imp, lm(y ~ x + z))
#combine results of all 5 models
combine <- pool(fit)

xyplot(imp,y ~ x + z,pch=18,cex=1) # compare distribution of original & imputed data
densityplot(imp) # magenta - density of the imputed data (each dataset); blue - the observed data
stripplot(imp, pch = 20, cex = 1.2)

imp_1 <- mice(mydata, m = 50)
densityplot(imp_1) 

library(dplyr)
library(mice)
library(foreign) # to import Stata DTA files
library(car)     # for recode

set.seed(145)

## Import ANES dataset
anesimp <- read.dta("anesimputation.dta", 
                    convert.factors = FALSE, missing.type = TRUE)

# Dataset contains values <0. Code all of them as missing 

for(i in 1:ncol(anesimp)){
  anesimp[,i] <- ifelse(anesimp[,i]<0, NA, anesimp[,i]) 
}

## Add occupation variable 

anesocc <- read.csv("anesocc.csv",sep=";",na.strings=c("","NA"))

# Selecting occupation now and industry now variables
anesocc2 <- anesocc %>%
  dplyr::select(caseid, dem_occnow, dem_indnow)

# Coding any text that includes "manu" in it as respondent working in
# manufacturing, excluding manuver
anesocc2 <- anesocc2 %>% 
  mutate(manuf = case_when((grepl("manu",dem_occnow)&!grepl("manuver",dem_occnow)) ~ 1,
                           grepl("manu",anesocc2$dem_indnow) ~ 1,
                           is.na(dem_occnow) ~ NA_real_,
                           is.na(dem_indnow) ~ NA_real_,
                           !is.na(dem_occnow) ~ 0,
                           !is.na(dem_indnow) ~ 0)
  )


anesocc2 <- anesocc2 %>% 
  dplyr::select(manuf)

# combining by columns as they are sorted in the same order
anesimp <- cbind(anesimp,anesocc2)

## Merge M&A data 
maimp <- read.dta("ma.dta")
anesimp <- merge(x=anesimp, y=maimp, by=c("sample_state"))

# Recode variables 
anesimp$patriot_amident <- recode(anesimp$patriot_amident, 
                                  "5=0; 4=1; 3=2; 2=3; 1=4")
anesimp$econ_ecnext_x <- recode(anesimp$econ_ecnext_x, 
                                "1=0; 2=1; 3=2; 4=3; 5=4")
anesimp$pid_x <- recode(anesimp$pid_x, 
                        "1=0; 2=1; 3=2; 4=3; 5=4; 6=5; 7=6")

anesimp$dem_edugroup_x <- recode(anesimp$dem_edugroup_x, 
                                 "1=0; 2=1; 3=2; 4=3; 5=4")

# Treat manuf as a factor 
anesimp$manuf <- as.factor(anesimp$manuf)


# Save the dataframe as another object so that we can use the original dataframe
# for multiple imputation
anesimpor <- anesimp 

## Transform variables for regression
# Treat nationalism as continuous
anesimpor$patriot_amident <- as.numeric(anesimpor$patriot_amident)
# Treat party id as continuous 
anesimpor$pid_x <- as.numeric(anesimpor$pid_x)
# Treat china_econ as dichotomous 
anesimpor$china_econ <- recode(anesimpor$china_econ, "1=0; 3=0; 2=1")
anesimpor$china_econ <- as.factor(anesimpor$china_econ)

# Take the log of Chinese M&A variables - add a small number as variable
# contains 0s
anesimpor$LogMANO <- log(anesimpor$MANo+1.01)
# Treat party id as continuous 

## Estimate an OLS regression

fitols <- lm(ft_hclinton ~ manuf + pid_x + patriot_amident + 
               china_econ + LogMANO, data=anesimpor)

summary(fitols)

#### Mice tutorial
# https://datascienceplus.com/imputing-missing-data-with-r-mice-package/
data <- airquality
data[4:10,3] <- rep(NA,7)
data[1:5,4] <- NA

data <- data[-c(5,6)]
summary(data)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(data,2,pMiss)
apply(data,1,pMiss)

mice::md.pattern(data) # using mice to check the missing data pattern
# 104 samples are complete
# 34 samples miss only the Ozone measurement
# 4 samples miss only the Solar.R value and so on.

library(VIM) # a more helpful visual representation via VIM
aggr_plot <- aggr(data, col=c('navyblue','red'), 
                  numbers=TRUE, sortVars=TRUE, labels=names(data), 
                  cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

marginplot(data[c(1,2)]) # Scatterplot with additional information in the margins
# The red box plot on the left shows the distribution of Solar.R with Ozone missing 
#     while the blue box plot shows the distribution of the remaining datapoints. 

## using mice() to take cares of the imputing
tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
# m=5, number of imputed datasets, 5 is the default value
# meth='pmm' refers to the imputation method, ppm is the predictive mean matching as imputation method
# to check other method, methods(mice)
summary(tempData)
tempData$imp$Ozone
# data$Ozone

completedData <- complete(tempData,1) # replaced with the imputed values in the first of the five datasets
completedData_alll <- (complete(tempData,1) + complete(tempData,2) + complete(tempData,3) + complete(tempData,4) + complete(tempData,5))/5
xyplot(tempData,Ozone ~ Wind+Temp+Solar.R,pch=18,cex=1) # compare distribution of original & imputed data
densityplot(tempData) # magenta - density of the imputed data (each dataset); blue - the observed data
stripplot(tempData, pch = 20, cex = 1.2)

# fit a linear model to the data
modelFit1 <- with(tempData,lm(Temp~ Ozone+Solar.R+Wind)) # fit a a model to each of the imputed dataset and then pool the results together
summary(pool(modelFit1))
fit1_pool = pool(modelFit1)
# modelFit1, containts the results of the fitting performed over the imputed datasets
# pool(), pools them all togethe

# we initialized the mice function with a specific seed
# the results are somewhat dependent on our initial choice
# To reduce this effect, we can impute a higher number of dataset as shown:
tempData2 <- mice(data,m=50,seed=245435)
modelFit2 <- with(tempData2,lm(Temp~ Ozone+Solar.R+Wind))
summary(pool(modelFit2))
complete()

#### missForest ####

# interpolation via Random Forest
# However, according to the link below, it does not make sense; 
# anyway, we have to face the situation that there are NAs in the predictors, when we can not use the prediction
# https://stackoverflow.com/questions/8370455/how-to-use-random-forests-in-r-with-missing-values
# https://www.analyticsvidhya.com/blog/2016/03/tutorial-powerful-packages-imputing-missing-values/
# another solution!
#load data
data("iris")

#seed 10% missing values
iris.mis <- prodNA(iris, noNA = 0.1)
summary(iris.mis)

#impute missing values, using all parameters as default values
iris.imp <- missForest(iris.mis)

#check imputed values
iris.imp$ximp

#check imputation error
iris.imp$OOBerror

#comparing actual data accuracy
iris.err <- mixError(iris.imp$ximp, iris.mis, iris)
iris.err



