# Padfield et al 2017 Ecology Letters
# Processing of ecosystem-level metabolism data 

# clear workspace
mise::mise(pkgs = TRUE, console = TRUE, figs = TRUE, vars = TRUE)

# load in packages ####
library(chron)
library(caTools)
library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(scales)

# functions used in the analysis ####

# function to calculate multiple AIC scores
all.AIC <- function(x){data.frame(model = x, AIC = MuMIn::AICc(get(x)))}

# path to graphs ####
path_fig <- 'plots'

# load in data ####
d_eco <- readRDS('data/eco_metab_raw.rds')

# make time as.posixct
d_eco <- mutate(d_eco, year = substr(utc, 1, 4), utc = as.POSIXct(strptime(utc, format='%Y-%m-%d %H:%M:%S'))) 

# create average dataframe over 15 minutes ####
d_ave <- 
  # take O2.merge
  d_eco %>%
  # group it by site, day and hourly time intervals
  group_by(site, day, year, time = cut(utc, breaks = '15 min')) %>%
  # calculate the mean DO, light, atmospheric O2 saturation concentration and temp at each hour
  summarise(DO = mean(DO), light = mean(mean.light), atm.sat = mean(atm.sat), temp = mean(temp)) %>%
  # calculate the difference between the T1 and T0 for each of these hourly estimates
  mutate(DO.change = c(0, diff(DO))) %>%
  # do some processing on the data to create required columns
  mutate(day_night = ifelse(light < 5, 'night', 'day'), 
         O2.sur.def = atm.sat - DO, time.length = 15, 
         time = as.POSIXct(strptime(time, format='%Y-%m-%d %H:%M:%S'))) %>%
  # add 7.5 minutes on for the midpoint between the T1 and T0
  mutate(time = time + seconds(60*7.5)) %>%
  mutate(time.of.day = as.POSIXct(strftime(time, format="%H:%M:%S"), format = "%H:%M:%S")) %>%
  data.frame()

# Stream characteristics ####
# load in stream trait data
StreamVars <- read.csv('data/StreamTraits_all.csv')

# units of Stream Vars
# reach - m
# width - m
# depth - m
# velocity - m/s
# reaeration coefficient - m/min

# add columns to Stream Vars ####
StreamVars <- StreamVars %>%
  mutate(area = length*width, 
         volume = length*width*depth,
         discharge = width*depth*velocity) %>%
  data.frame()

# load in chlorophyll data ####
chl <- read.csv('data/chl_stream.csv')

# get averages of chl for each site
av_chl <- chl %>% 
  group_by(site) %>% 
  summarise(., Chla = mean(Chla), Chl_tot = mean(Chl_tot)) %>% 
  # changing units to g Chla m-2
  mutate(., Chla = Chla * 10^-2,
         year = 2016) %>% 
  data.frame()

# merge stream variables and chlorophyll data
StreamVars <- left_join(StreamVars, av_chl, by = c('year', 'site'))

# merge dataframes by site
d <- merge(d_ave, StreamVars, by = c('site', 'year'))

# Calculating reaeration ####

# Temperature correction for reaeration rate, temperature needs to be in ºC
Ktemp <- function(K20, temp){
  return(K20 * 1.025^(temp - 20))
}

# units of equation
# Velocity (V) - m/s
# Mean Depth (D) - m
# K is in m/min

# reaeration equation (from Owens 1974)
K <- function(V, D){
  K <- 50.8*((V*100)^(0.67))*(D*100)^(-0.85) 
  return(K/100/60)
}

d <- mutate(d, reaeration = K(d$velocity, d$depth)*15,
               reaer.coef = reaeration/depth)

# further processing ####
d <- mutate(d, 
               # calculate temperature corrected reaearation
               Ktemp = Ktemp(reaeration, temp),
               reaer.coef_temp = Ktemp(reaer.coef, temp),
               # calculate distance up stream sampled
               dist.samp = 3*velocity*60*15/reaer.coef,
               # work out how much oxygen exchange occurs per timepoint
               DO.change.area = depth*DO.change,
               O2.exch = Ktemp*O2.sur.def,
               # account for the change in rate in terms of depth - into m-2 hr-1
               rate.cor = DO.change.area - O2.exch,
               # assign a flux type to the data
               flux = ifelse(d$day_night == 'night', 'R', 'NPP')) %>%
  data.frame()

# calculating rates dataframe of entire stream metabolism ####

# metabolism functions ####

ER <- function(rate.change, night_or_day, time_interval){
  R <- mean(rate.change[night_or_day=='night'],na.rm=TRUE)*(1440/mean(time_interval))
  return(R)
}

GPP <- function(rate.change, night_or_day, time_interval){
  PP <- (mean(rate.change[night_or_day == 'day'],na.rm=TRUE) - mean(rate.change[night_or_day=='night'],na.rm=TRUE)) * (1440/mean(time_interval))
  return(PP)
}

# change some of the site values to match up with other analysis
d$site <- plyr::mapvalues(d$site, from = c('S11B', 'S5', 'S1', 'S7'), to = c('S11B_high', 'S5_high', 'S1_high', 'S7_high'))

# create final metabolism data set
d_metab <- d %>% group_by(site, year, day) %>%
  summarise(temp = mean(temp, na.rm=TRUE),
            pH = mean(pH),
            conduct = mean(conductivity),
            light = mean(light, na.rm=TRUE),
            tot_light = sum(light, na.rm = T),
            Chla = mean(Chla, na.rm=TRUE),
            width = mean(width,na.rm=TRUE),
            depth = mean(depth,na.rm=TRUE),
            velocity = mean(velocity,na.rm=TRUE),
            reaeration = mean(reaeration,na.rm=TRUE),
            reaer.coef = mean(reaer.coef)/15,
            dist.samp = mean(dist.samp),
            cross.sec.area = mean(width*depth),
            volume = mean(dist.samp * width * depth),
            ER = abs(ER(rate.cor, day_night, time.length)),
            PP = GPP(rate.cor, day_night, time.length),
            NO2 = mean(NO2),
            NO2_NO3 = mean(NO2_NO3),
            NH4 = mean(NH4),
            PO4 = mean(PO4)) %>%
  mutate(., ikt = (1/8.62e-05/(mean(.$temp, na.rm = TRUE) + 273.15)) - 
           (1/8.62e-05/(temp + 273.15))) %>%
  data.frame()

# delete any rows where PP is < 0
d_metab <- d_metab[d_metab$PP > 0,]

d_metab <- mutate(d_metab, NEP = PP - ER) %>%
  gather(., 'flux', 'rate', c(ER, PP, NEP)) %>%
  mutate(., 
         # calculate log rate
         log.rate = log(rate)) %>%
  data.frame()

# delete all fluxes apart from GPP from the analysis
d_metab$flux <- as.factor(d_metab$flux)
d_metab <- d_metab[d_metab$flux != 'NEP',]
d_metab <- d_metab[d_metab$flux != 'ER',]

# save data
saveRDS(d_metab, 'data/eco_level_metab.rds')

# now move onto eco_level_metab_analysis.R ####

# SI plots ####

# Figure S3. temperature vs time
ggplot(d) +
  geom_point(aes(time, temp, col = day_night), shape = 21, show.legend = F) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  scale_color_manual(values = c('orange', 'black')) +
  facet_wrap(~ site + year + day, scales = 'free', labeller = labeller(.multi_line = F)) +
  ylab('Temperature (ºC)') +
  xlab('Time') +
  scale_x_datetime(labels = date_format(format = "%H:%M", tz = 'GMT-1')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Figure S4. light vs time
ggplot(d) +
  geom_point(aes(time, light, col = day_night), shape = 21, show.legend = F, alpha = 0.5) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  scale_color_manual(values = c('orange', 'black')) +
  facet_wrap(~ year + day, scales = 'free', labeller = labeller(.multi_line = F)) +
  ylab(expression(Light~(µmol^-1~m^-2~s^-1))) +
  xlab('Time') +
  scale_x_datetime(labels = date_format(format = "%H:%M", tz = 'GMT-1')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Figure S5. rate vs time
ggplot(d) +
  geom_point(aes(time, rate.cor, col = day_night), shape = 21, show.legend = F) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  scale_color_manual(values = c('orange', 'black')) +
  facet_wrap(~ site + year + day, scales = 'free', labeller = labeller(.multi_line = F)) +
  ylab(expression(Corrected~d_eco~(g~O[2]~m^-2~min^-1))) +
  xlab('Time') +
  scale_x_datetime(labels = date_format(format = "%H:%M", tz = 'GMT-1')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))