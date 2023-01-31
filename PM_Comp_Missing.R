##packages in need
library(tidyr) 
library(stats) 
library(ggplot2)
library(scales) 
library(dplyr)
library(imputeTS) #na_ma, na_interpolation ect.
library(mice) # using Markov Chain Monte Carlo simulation to impute the missing entries
library(tibble)
library(missForest) # implementation of random forest algorithm
library(ggsci)

####################################################################################
####### IMPROVE - Fill the Missing #######
####################################################################################
#### generate basic data for filling ####
imp_daily = read.csv("IMPROVE_Component_with_missing.csv")
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

# exclude states not located in the "Continental and mainland USA"
imp_miss = subset(imp_miss, !(State %in% c("AK", "HI", "PR", "VI"))) 
imp_miss = imp_miss[with(imp_miss, order(State, SiteCode, Date)), ]

## pattern of missing data exploration
p_missing <- data.frame(unlist(lapply(imp_miss, function(x) sum(is.na(x))))/nrow(imp_miss))
colnames(p_missing) = "MissingRate"

n.site = length(unique(imp_miss$SiteCode))

#### filling the NAs for each site ####
imp_avg_month_summary = imp_running_avg_summary = imp_linear_summary = imp_sgl_mice_summary = imp_sgl_rf_summary = NULL
interp_org_cmp_summary = NULL

#### filling the NAs for each site (logged, no negative)- na.omit & check Mix Error ####
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


