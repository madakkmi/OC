#===============================================================================
# Formatting the two-sample MR analysis results
#===============================================================================

rm(list=ls())

library(TwoSampleMR)
library(tidyverse)
library(openxlsx)

setwd("G:/UKBiobank_AnnualReport_urgent!/Iqbal/madakkatel2024large/data/interim")

# load the two-sample MR outputs
load("mr_results_mr-base_all.RData")

# add pleotropy info from MR-Egger
pleot <- mr_pleiotropy_test(dat)

# Filtering based on the latest list of included expoures (id.exposure)
incld <- read.xlsx("filtered_and _ordered_exposures_from_gbdt_final.xlsx")
mr_results <- mr_results %>% filter(id.exposure%in%incld$id.exposure)

# generating confidence interval for IVW and wald methods (from normal distribution)
mr_results.a <- mr_results %>%
  filter(method%in%c("Inverse variance weighted","Wald ratio")) %>%
  mutate(or = round(exp(b),2),
         lci = round(exp(b-1.96*se),2),
         uci = round(exp(b+1.96*se),2),
         pval = sprintf(pval, fmt = '%#.1e'))

# Generating confidence intervals for MR-Egger (t-distribution)
mr_results.b <- mr_results %>%
  filter(method%in%"MR Egger")

mrout_list <- list()
for(idexp in mr_results.b$id.exposure) {
  subset <- mr_results.b[mr_results.b$id.exposure%in%idexp,]
  degfrd <- subset$nsnp[1]
  subset$or = round(exp(subset$b),2)
  subset$lci = round(exp(subset$b-(qt(0.975,df=degfrd)*subset$se)),2)
  subset$uci = round(exp(subset$b+(qt(0.975,df=degfrd)*subset$se)),2)
  subset$pval = sprintf(subset$pval, fmt = '%#.1e')

  mrout_list[[idexp]] <- subset

}
# convert list of dataframes to one big dataframe
mr_results.b2 <- do.call(rbind,mrout_list)
rownames(mr_results.b2) <- NULL

# rowbind the IVW ( and ratio method) data to MR-Egger data
mr_results.combined <- rbind(mr_results.a,mr_results.b2)

# Sub-setting by each outcomes & change to wide format
mr_results.oc.1 <- mr_results.combined %>% filter(id.outcome%in%'ieu-a-1120')
pleot.oc.1 <- pleot %>%
  filter(id.outcome%in%'ieu-a-1120') %>%
  select(id.exposure,pval)
names(pleot.oc.1)[2] <- "p-pleotropy"

mr_results.oc.1 <- merge(mr_results.oc.1, pleot.oc.1, by="id.exposure", all = TRUE)

mr_results.oc.2 <- mr_results.combined %>% filter(id.outcome%in%'ieu-a-1228')
pleot.oc.2 <- pleot %>%
  filter(id.outcome%in%'ieu-a-1228') %>%
  select(id.exposure,pval)
names(pleot.oc.2)[2] <- "p-pleotropy"
mr_results.oc.2 <- merge(mr_results.oc.2, pleot.oc.2, by="id.exposure", all = TRUE)

mr_results.oc.3 <- mr_results.combined %>% filter(id.outcome%in%'ieu-a-1125')
pleot.oc.3 <- pleot %>%
  filter(id.outcome%in%'ieu-a-1125') %>%
  select(id.exposure,pval)
names(pleot.oc.3)[2] <- "p-pleotropy"
mr_results.oc.3 <- merge(mr_results.oc.3, pleot.oc.3, by="id.exposure", all = TRUE)

mr_results.oc.4 <- mr_results.combined %>% filter(id.outcome%in%'ieu-a-1124')
pleot.oc.4 <- pleot %>%
  filter(id.outcome%in%'ieu-a-1124') %>%
  select(id.exposure,pval)
names(pleot.oc.4)[2] <- "p-pleotropy"
mr_results.oc.4 <- merge(mr_results.oc.4, pleot.oc.4, by="id.exposure", all = TRUE)

mr_results.oc.5 <- mr_results.combined %>% filter(id.outcome%in%'ieu-a-1123')
pleot.oc.5 <- pleot %>%
  filter(id.outcome%in%'ieu-a-1123') %>%
  select(id.exposure,pval)
names(pleot.oc.5)[2] <- "p-pleotropy"
mr_results.oc.5 <- merge(mr_results.oc.5, pleot.oc.5, by="id.exposure", all = TRUE)

# Format the estimates from each exposure against each outcome - in wide format
mr_results.oc.wide <- mr_results.oc.1 %>%
  cbind(mr_results.oc.2[,c("outcome","or","lci","uci","pval","p-pleotropy")],
        mr_results.oc.3[,c("outcome","or","lci","uci","pval","p-pleotropy")],
        mr_results.oc.4[,c("outcome","or","lci","uci","pval","p-pleotropy")],
        mr_results.oc.5[,c("outcome","or","lci","uci","pval","p-pleotropy")]
  )

# Add col info
selectedcol <- c("id","trait", "author","consortium","sex","population","sample_size","ncase","ncontrol")
exposure_selected2 <- exposure_selected[,selectedcol]
names(exposure_selected2)[1] <- "id.exposure"

mr_results.oc.wide <- merge(exposure_selected2,mr_results.oc.wide,by="id.exposure",all = TRUE)
mr_results.oc.wide <- mr_results.oc.wide %>% filter(!is.na(b))

#=====================================================================================================
# add more info about the exposures (PLEASE note this can be done in the initial MR-analysis script)
ad.info <-  exposure_selected %>% filter(id%in%mr_results.oc.wide$id.exposure) %>%
  select(id,year,consortium,pmid,unit,note)
names(ad.info)[1] <-  "id.exposure"
mr_results.oc.wide <- mr_results.oc.wide %>%  left_join(ad.info)
#=====================================================================================================


# for ordering - based on "incld"
str(incld)
mr_results.oc.wide <- left_join(incld, mr_results.oc.wide)
mr_results.oc.wide <- mr_results.oc.wide %>%
  mutate(seq2=case_when(method%in%"Inverse variance weighted"~1,
                        method%in%"Wald ratio"~2,
                        method%in%"MR Egger"~3)
  )

# arrange rows and select coloumns
mr_results.oc.wide <- mr_results.oc.wide %>%
  arrange(seq,seq2)      %>%            # sort rows based on seq(from "incld" file and seq2 (method))
  filter(id.exposure%in%list_of_studies) %>%
  select(group,trait,source,id.exposure,author,sex,population,sample_size,nsnp,method,or,lci,uci,pval,`p-pleotropy`,
         or.1,lci.1,uci.1,pval.1,`p-pleotropy.1`,or.2,lci.2,uci.2,pval.2,`p-pleotropy.2`,
         or.3,lci.3,uci.3,pval.3,`p-pleotropy.3`,or.4,lci.4,uci.4,pval.4,`p-pleotropy.4`)

# save the output
#write.xlsx(mr_results.oc.wide,"mr_results.o_v4_allcombined.xlsx")
