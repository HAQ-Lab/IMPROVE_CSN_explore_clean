##set working directory
setwd("/Users/TingZhang/Documents/HEI HAQ PMF/PMF_try/for R")
getwd()
data.dir <- "/Users/TingZhang/Documents/HEI HAQ PMF/PMF_try/for R"

#### CSN EPA-MDL, PMF output linear - clusters ####
csn_pmf_base = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/for R/EPA_CSN_Cluster6_F9.csv")
csn_pmf_base$X = NULL
factor.number = ncol(csn_pmf_base)
cluster.study = "Cluster6"
# cluster.study = "Cluster24"

# there is formatting wrong in date due to influence of excel, so use the date from concentration file
csn_conc = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/EPA MDL/CSN_Cluster6_2011-14_PMF_conc_EPAmdl.csv")
csn_pmf_base$Date = csn_conc$Date
csn_pmf_base$Date = as.Date(csn_pmf_base$Date)

# match the total PM2.5 mass concentration
csn_pmf_base$P25 = csn_conc$PM25

# match the sidecode of sites in given cluster
site_row = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/EPA MDL/CSN_Cluster6_date_number_of_site_EPAmdl.csv")
csn_pmf_base$SiteCode = rep(unique(site_row$SiteCode), site_row$Freq)

csn_pmf_base = csn_pmf_base[!duplicated(csn_pmf_base), ]
csn_pmf_base <- csn_pmf_base[-2, ]

###### CSN EPA-MDL, lm, all sites ######
# use linear regression to get the factor contribution
csn_pmf_base_forLM = csn_pmf_base
csn_pmf_base_forLM$Date = csn_pmf_base_forLM$SiteCode = NULL
csn_pmf_base_lm = lm(P25 ~., data = csn_pmf_base_forLM)
csn_pmf_base_lm_beta = summary(csn_pmf_base_lm)$coefficients[, 1]
csn_pmf_base_lm_beta = data.frame(csn_pmf_base_lm_beta[2:length(csn_pmf_base_lm_beta)])
colnames(csn_pmf_base_lm_beta) = "lm.factor.beta.cluster"
csn_pmf_base_lm_beta$Factor = rownames(csn_pmf_base_lm_beta)
csn_pmf_base_lm_beta$Factor.contribution = (csn_pmf_base_lm_beta$lm.factor.beta.cluster/
                                              sum(csn_pmf_base_lm_beta$lm.factor.beta.cluster))*100

write.csv(csn_pmf_base_lm_beta, 
          paste0("CSN_", cluster.study, "F", factor.number, "2011-14_EPAmdl_PMF_overall.csv"))

###### CSN EPA-MDL, lm, each sites ######
csn_pmf_base_site = subset(csn_pmf_base,
                            SiteCode == site_row$SiteCode[3])
csn_pmf_base_site_lm = lm(P25 ~ Factor.1 + Factor.2 + Factor.3 + Factor.4 + 
                             Factor.5 + Factor.6 + Factor.7 + Factor.8 + Factor.9, 
                           data = csn_pmf_base_site)
csn_pmf_base_site_lm_beta = summary(csn_pmf_base_site_lm)$coefficients[, 1]
csn_pmf_base_site_lm_beta = data.frame(csn_pmf_base_site_lm_beta[2:length(csn_pmf_base_site_lm_beta)])
colnames(csn_pmf_base_site_lm_beta) = "lm.beta.site-ClusterBased"
csn_pmf_base_site_lm_beta$Factor = rownames(csn_pmf_base_site_lm_beta)
csn_pmf_base_site_lm_beta$Factor.contribution = (csn_pmf_base_site_lm_beta$lm.beta.site/
                                                    sum(csn_pmf_base_site_lm_beta$lm.beta.site))*100

write.csv(csn_pmf_base_site_lm_beta, "CSN_Cluster6_2011-14_EPAmdl_PMF_F7_site3.csv")


#### CSN EPA-MDL, PMF output linear - clusters ####
csn_pmf_1site_base = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/for R/EPA_CSN_Cluster6_site3_F9.csv")
csn_pmf_1site_base$X = NULL

csn_pmf_1site_conc = read.csv("/Users/TingZhang/Library/CloudStorage/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/EPA MDL/CSN_Cluster6_site3_2011-14_PMF_conc_EPAmdl.csv")
csn_pmf_1site_base$Date = csn_pmf_1site_conc$Date
csn_pmf_1site_base$PM25 = csn_pmf_1site_conc$PM25

csn_pmf_1site_lm = lm(PM25 ~ Factor.1 + Factor.2 + Factor.3 + Factor.4 + 
                            Factor.5 + Factor.6 + Factor.7 + Factor.8 + Factor.9, 
                          data = csn_pmf_1site_base)
csn_pmf_1site_lm_beta = summary(csn_pmf_1site_lm)$coefficients[, 1]
csn_pmf_1site_lm_beta = data.frame(csn_pmf_1site_lm_beta[2:length(csn_pmf_1site_lm_beta)])
colnames(csn_pmf_1site_lm_beta) = "lm.beta.site-ClusterBased"
csn_pmf_1site_lm_beta$Factor = rownames(csn_pmf_1site_lm_beta)
csn_pmf_1site_lm_beta$Factor.contribution = (csn_pmf_1site_lm_beta$lm.beta.site/
                                                   sum(csn_pmf_base_site_lm_beta$lm.beta.site))*100




#### CSN EPA-MDL, Cluster6 & 7 factors ####
csn_eap_C6f7 = read.csv("EPA_MDL_CSN_Cluster6_F7.csv")
colnames(csn_eap_C6f7)[1] = "Date"

# there is formatting wrong in date due to influence of excel, so use the date from concentration file
csn_epa_C6f7_conc = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/EPA MDL/CSN_Cluster6_2011-14_PMF_conc_EPAmdl.csv")
csn_eap_C6f7$Date = csn_epa_C6f7_conc$Date
csn_eap_C6f7$Date = as.Date(csn_eap_C6f7$Date, format = "%m/%d/%y")

# match the total PM2.5 mass concentration
csn_eap_C6f7$P25 = csn_epa_C6f7_conc$PM25

# match the sidecode
site_row = read.csv("/Users/ztttttt/Dropbox/HEI_PMF_files_Ting/CSN_PMF_Try_data/EPA MDL/CSN_Cluster6_date_number_of_site_EPAmdl.csv")
csn_eap_C6f7$SiteCode = rep(unique(site_row$SiteCode), site_row$Freq)

csn_eap_C6f7 = csn_eap_C6f7[!duplicated(csn_eap_C6f7), ]
csn_eap_C6f7 <- csn_eap_C6f7[-2, ]

###### CSN EPA-MDL C6F7, lm, all sites ######
# use linear regression to get the factor contribution
csn_eap_C6f7_lm = lm(P25 ~ Factor.1 + Factor.2 + Factor.3 + Factor.4 + 
                       Factor.5 + Factor.6 + Factor.7, data = csn_eap_C6f7)
csn_eap_C6f7_lm_beta = summary(csn_eap_C6f7_lm)$coefficients[, 1]
csn_eap_C6f7_lm_beta = data.frame(csn_eap_C6f7_lm_beta[2:length(csn_eap_C6f7_lm_beta)])
colnames(csn_eap_C6f7_lm_beta) = "lm.factor.beta.cluster"
csn_eap_C6f7_lm_beta$Factor = rownames(csn_eap_C6f7_lm_beta)
csn_eap_C6f7_lm_beta$Factor.contribution = (csn_eap_C6f7_lm_beta$lm.factor.beta.cluster/
                                              sum(csn_eap_C6f7_lm_beta$lm.factor.beta.cluster))*100

write.csv(csn_eap_C6f7_lm_beta, "CSN_Cluster6_2011-14_EPAmdl_PMF_F7_overall.csv")

###### CSN EPA-MDL C6F7, lm, each sites ######
csn_eap_C6f7_site3 = subset(csn_eap_C6f7,
                            SiteCode == site_row$SiteCode[3])
csn_eap_C6f7_site3_lm = lm(P25 ~ Factor.1 + Factor.2 + Factor.3 + Factor.4 + 
                             Factor.5 + Factor.6 + Factor.7, 
                           data = csn_eap_C6f7_site3)
csn_eap_C6f7_site3_lm_beta = summary(csn_eap_C6f7_site3_lm)$coefficients[, 1]
csn_eap_C6f7_site3_lm_beta = data.frame(csn_eap_C6f7_site3_lm_beta[2:length(csn_eap_C6f7_site3_lm_beta)])
colnames(csn_eap_C6f7_site3_lm_beta) = "lm.beta.site-ClusterBased"
csn_eap_C6f7_site3_lm_beta$Factor = rownames(csn_eap_C6f7_site3_lm_beta)
csn_eap_C6f7_site3_lm_beta$Factor.contribution = (csn_eap_C6f7_site3_lm_beta$lm.beta.site/
                                                    sum(csn_eap_C6f7_site3_lm_beta$lm.beta.site))*100

write.csv(csn_eap_C6f7_site3_lm_beta, "CSN_Cluster6_2011-14_EPAmdl_PMF_F7_site3.csv")




