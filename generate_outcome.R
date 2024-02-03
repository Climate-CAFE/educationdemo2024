######################################################################################################
## code to generate synthetic daily mortality time series for CA counties, 01/01/2015-12/31/2019    ##
## to be used to in a synthetic study of associations between wildfire pm2.5 exposure and mortality ##
## to demonstrate climate and health analytic methods in 2024 CAFE conference educational demo      ##
## created by: Rachel Nethery                                                                       ##
######################################################################################################

library(eesim)
library(lubridate)
library(sf)
library(dplyr)
library(tidyr)
library(gt)
library(data.table)
library(tidycensus)
library(ggplot2)
library(scales)

## set seed for reproducibility ##
set.seed(567)

##################################################################################################################################################
## county-level wildfire-attributable pm2.5 concentrations for the contiguous US downloaded from https://www.stanfordecholab.com/wildfire_smoke ##
## below the code from associated readme to load data and add non-smoke days into the time series (with zero values of smoke PM2.5)             ##
##################################################################################################################################################

## Load daily county smokePM estimates on smoke days ##
preds = readRDS("county/smokePM2pt5_predictions_daily_county_20060101-20201231.rds")

## Load county shapefile ##
counties = read_sf("county/tl_2019_us_county")

## Load full set of dates ##
dates = seq.Date(ymd("20060101"), ymd("20201231"), by = "day")

## Get full combination of county-days ##
## Warning: this may require a large amount of memory ##
out = expand.grid(GEOID = counties$GEOID, date = dates)

## Match smokePM predictions on smoke days to county-days ##
out = left_join(out, preds, by = c("GEOID", "date"))

## Predict 0 for remaining county-days, which are non-smoke days ##
out = mutate(out, smokePM_pred = replace_na(smokePM_pred, 0))

############################################################################
## subset to only CA counties (counties with fips codes starting with 06) ##
############################################################################

setDT(out)

ca<-out[substr(GEOID,1,2)=='06']

## qa/qc check: should be 58 unique counties in CA ##
length(unique(ca$GEOID))

####################################
## subset to only years 2015-2019 ##
####################################

ca<-ca[date>=ymd('20150101') & date<=ymd('20191231')]

## qa/qc check: should be 5*365 + 2 (leap years)=1826 unique dates ##
length(unique(ca$date))

## rename exposure to 'x' for use in eesims::sim_outcome ##
setnames(ca,'smokePM_pred','x')

#########################################
## load in CA county population counts ##
#########################################

## will use the population estimates from 2015 PEP ##
popest<-get_estimates(geography = 'county',year = 2015,variables = 'POP')
setDT(popest)

## subset to CA only ##
popest<-popest[substr(GEOID,1,2)=='06']

## qa/qc check: should be 58 counties ##
nrow(popest)

## remove unnecessary columns ##
popest[,':='(NAME=NULL,variable=NULL)]

## rename column with populations ##
setnames(popest,'value','pop2015')

## add in baseline counts of the outcome for each county ( rates drawn from unif[1/30000, 1/20000] ) ##
popest[,bline:=pop2015*runif(nrow(popest),min=1/30000,max=1/20000)]

## add in a relative risk for each county ( drawn from unif[1.001, 1.01] ) ##
popest[,rr:=runif(nrow(popest),min=1.001,max=1.01)]

##############################################################################################
## simulate daily synthetic mortality counts using exposures, pop counts, and eesim package ##
##############################################################################################

## loop through counties (have to simulate separately by county) ##
outdf<-NULL
for (i in 1:nrow(popest)){
  testout <- sim_outcome(exposure = ca[GEOID==popest$GEOID[i],.(date,x)], 
                         average_outcome = popest$bline[i],
                         trend = "cos1", 
                         rr = popest$rr[i],
                         amp=.1,
                         start.date = '2015-01-01')
  
  outdf<-rbind(outdf,data.frame('GEOID'=popest$GEOID[i],testout))
}

setDT(outdf)

## make time series plots of smoothed outcome time series by county ##
pdf('synth_ts_plots.pdf',
    height=6,width=9)
for (i in 1:nrow(popest)){
  print(
  ggplot(outdf[GEOID==popest$GEOID[i]], aes(x = date, y = outcome)) +
    geom_line() +
    geom_smooth(method='loess',span=.1,se=F)+
    ggtitle(paste0('County ',popest$GEOID[i],', Population=',popest$pop2015[i]))+
    scale_x_date(breaks=ymd(c('20150101','20160101','20170101','20180101','20190101','20200101')))+
    theme_bw()
  )
}
dev.off()


## export the merged wildfire pm2.5 and synthetic mortality dataset ##
write.csv(outdf,file='full_data.csv',row.names=F)

## remove the wildfire pm2.5 concentrations and export (for dataverse) ##
outdf_sub<-subset(outdf,select=c('GEOID','date','outcome'))
write.csv(outdf_sub,file='synthetic_mortality.csv',row.names = F)