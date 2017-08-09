# Padfield et al 2017 Ecology Letters
# Analysis of ecosystem-level metabolism

# clear workspace
mise::mise(pkgs = TRUE, console = TRUE, figs = TRUE, vars = TRUE)

# load in packages ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(lme4)
library(nlme)
library(grid)
library(gridExtra)
library(nlsLoop)
library(broom)

# path to graphs ####
path_fig <- 'plots'

# load in data ####
d_eco <- readRDS('data/eco_level_metab.rds')

# 1. Mixed model of the temperature dependence across years and sites and days ####
MM_all_years_1 <- lme(log.rate ~ ikt, random = ~1|site/year, d_eco, na.action = na.omit, method = 'ML')
MM_all_years_2 <- lme(log.rate ~ 1, random = ~1|site/year, d_eco, na.action = na.omit, method = 'ML')

MuMIn::AICc(MM_all_years_1, MM_all_years_2)
anova(MM_all_years_1, MM_all_years_2)

# best model is model 1, refit with REML
MM_final <- lmer(log.rate ~ ikt + (1|site/year), d_eco, na.action = na.omit)

# Work out r squared and confidence interval around the slope
MuMIn::r.squaredGLMM(MM_final)
lsmeans::lstrends(MM_final, ~ ikt, var = 'ikt')

# create predictions for the model
preds <- data.frame(ikt = seq(min(d_eco$ikt, na.rm = F), max(d_eco$ikt, na.rm = F), length = 20), temp = seq(min(d_eco$temp, na.rm = F), max(d_eco$temp, na.rm = F), length = 20))
# new y
preds$pred <- as.vector(predict(MM_final, re.form = NA, preds))

# compute standard error for predictions
# create design matrix of the model
Designmat <- model.matrix(~ikt, preds)
preds$predvar <- as.vector(diag(Designmat %*% vcov(MM_final) %*% t(Designmat)))

# add confidence around predictions
preds <- mutate(preds, conf.low = pred - (1.96*sqrt(predvar)),
                conf.high = pred + 1.96*sqrt(predvar))

# plot catchment temperature sensitivity ####
d_eco_plot <- ggplot(d_eco) +
  geom_point(aes(temp, log.rate), size = 1.25) +
  geom_line(aes(temp, pred), linetype = 2, preds) +
  #geom_ribbon(aes(x = temp, ymin = pred.int.low, ymax = pred.int.high), alpha=0.1, fill = "green4", preds) +
  geom_ribbon(aes(x = temp, ymin = conf.low, ymax = conf.high), alpha=0.2, fill = "green4", preds) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  xlab(expression(Stream~Temperature~(degree*C))) +
  ylab(expression(ln~Gross~Primary~Productivity~(g~O[2]~m^-2~day^-1)))

# Analysis of just 2016 data ####

# subset data for just 2016
metab_2016 <- filter(d_eco, year == 2016) %>%
  # recalculate ikt for this subset of data to recentre the data
  select(., -ikt) %>%
  mutate(., ikt = (1/8.62e-05/(mean(.$temp, na.rm = TRUE) + 273.15)) - 
           (1/8.62e-05/(temp + 273.15)),
         Chla_area = Chla*width*dist.samp) %>%
  mutate(log.Chla = log(Chla),
         log.Chla_area = log(Chla_area)) %>%
  mutate(standard_rate_area = log.rate - log.Chla_area) %>%
  data.frame()

# 2. Chlorophyll vs temperature ####
fit_chl_temp <- lm(log.Chla ~ ikt, metab_2016)
fit_chl_temp2 <- lm(log.Chla ~ 1, metab_2016)
anova(fit_chl_temp, fit_chl_temp2)
AIC(fit_chl_temp, fit_chl_temp2)
# best model is fit_chl_temp

# create predictions predictions
chl_temp_preds <- data.frame(ikt = seq(min(metab_2016$ikt, na.rm = T), max(metab_2016$ikt, na.rm = T), length.out = 20), temp = seq(min(metab_2016$temp), max(metab_2016$temp), length.out = 20)) %>%
  mutate(., preds = predict(fit_chl_temp, .),
         # get confidence intervals around predictions
         lwr = predict(fit_chl_temp, ., interval = 'confidence')[,2],
         upr = predict(fit_chl_temp, ., interval = 'confidence')[,3])

# plot
plot_chl_temp <- ggplot() +
  geom_point(aes(x = temp, y = log.Chla), size = 1.25, metab_2016) +
  geom_line(aes(temp, preds), linetype = 2, chl_temp_preds) +
  geom_ribbon(aes(x = temp, ymin = lwr, ymax = upr), fill = 'green4', alpha = 0.2, chl_temp_preds) +
  xlab(expression(Stream~Temperature~(degree*C))) + 
  ylab(expression(ln~Biomass~Density~(g~Chl*italic(a)~m^-2))) +
  theme_bw(base_size = 12, base_family = 'Helvetica')

# 3. log rate vs temperature vs chlorophyll of area sampled
MM_2016 <- lme(log.rate ~ ikt + log.Chla_area, random = ~1|site, metab_2016, method = 'ML')
MM_2016_2 <- lme(log.rate ~ log.Chla_area, random = ~1|site, metab_2016, method = 'ML')
MM_2016_3 <- lme(log.rate ~ ikt, random = ~1|site, metab_2016, method = 'ML')
MM_2016_4 <- lme(log.rate ~ 1, random = ~1|site, metab_2016, method = 'ML')

MuMIn::AICc(MM_2016, MM_2016_2, MM_2016_3, MM_2016_4)
anova(MM_2016_2, MM_2016_4)

# best model is model 2
# refit using REML
rate_chloro_mod <- lmer(log.rate ~ log.Chla_area + (1|site), metab_2016, REML = 'T')

# get predictions
preds.chl <- data.frame(log.Chla_area = seq(min(metab_2016$log.Chla_area, na.rm = T), max(metab_2016$log.Chla_area, na.rm = T), length = 20))
# new y
preds.chl$pred <- predict(rate_chloro_mod, re.form = NA, preds.chl)
# create design matrix
Designmat <- model.matrix(~log.Chla_area, preds.chl)

# compute standard error for predictions
preds.chl$predvar <- as.vector(diag(Designmat %*% vcov(rate_chloro_mod) %*% t(Designmat)))
# add confidence interval columns
preds.chl <- mutate(preds.chl, conf.low = pred - (1.96*sqrt(predvar)),
                    conf.high = pred + 1.96*sqrt(predvar))

# get r squared and confidence interval around slope 
MuMIn::r.squaredGLMM(rate_chloro_mod)
lsmeans::lstrends(rate_chloro_mod, ~ log.Chla_area, 'log.Chla_area')
confint(rate_chloro_mod, method = 'boot', nsim = 1000)

# plot
rate_Chla_plot <- ggplot(metab_2016) +
  geom_point(aes(log.Chla_area, log.rate), size = 1.25) +
  geom_line(aes(log.Chla_area, pred), linetype = 2, preds.chl) +
  #geom_ribbon(aes(x = log.Chla_width, ymin = pred.int.low, ymax = pred.int.high), alpha=0.1, fill = "green4", preds.chl) +
  geom_ribbon(aes(x = log.Chla_area, ymin = conf.low, ymax = conf.high), alpha=0.2, fill = "green4", preds.chl) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  xlab(expression(ln~Biomass~(g~Chl*italic(a)~stream^-1))) +
  ylab(expression(ln~Gross~Primary~Productivity~(g~O[2]~m^-2~day^-1)))

# 4. model rate per unit area against temperature ####
MM_chlrate <- lme(standard_rate_area ~ ikt, random = ~1|site/day, metab_2016, method = 'ML')
MM_chlrate2 <- lme(standard_rate_area ~ 1, random = ~1|site/day, metab_2016, method = 'ML')
anova(MM_chlrate, MM_chlrate2)
MuMIn::AICc(MM_chlrate, MM_chlrate2)
# MM_chlrate2 is the best model

# refit model with REML
MM_chlrate2 <- lmer(standard_rate_area ~ 1 + (1|site), metab_2016)

# create predictions
preds.chlrate <- data.frame(ikt = seq(min(metab_2016$ikt, na.rm = F), max(metab_2016$ikt, na.rm = F), length = 20), temp = seq(min(metab_2016$temp, na.rm = F), max(metab_2016$temp, na.rm = F), length = 20))
# new y
preds.chlrate$pred <- as.vector(predict(MM_chlrate2, re.form = NA, preds.chlrate))
# create design matrix
Designmat <- model.matrix(~1, preds.chlrate)

# compute standard error for predictions
preds.chlrate$predvar <- as.vector(diag(Designmat %*% vcov(MM_chlrate2) %*% t(Designmat)))

# add confidence interval columns
preds.chlrate <- mutate(preds.chlrate, conf.low = pred - (1.96*sqrt(predvar)),
                        conf.high = pred + 1.96*sqrt(predvar))

# fourth panel
rate_control_Chla_plot <- ggplot(metab_2016) +
  geom_point(aes(temp, standard_rate_area),size = 1.25) +
  geom_line(aes(temp, pred), linetype = 2, preds.chlrate) +
  #geom_ribbon(aes(x = temp, ymin = pred.int.low, ymax = pred.int.high), alpha=0.1, fill = "green4", preds.chlrate) +
  geom_ribbon(aes(x = temp, ymin = conf.low, ymax = conf.high), alpha=0.2, fill = "green4", preds.chlrate) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  xlab(expression(Stream~Temperature~(degree*C))) +
  ylab(expression(ln~Gross~Primary~Productivity~(g~O[2]~g~Chla^-1~day^-1))) + 
  ylim(-2, 4)

Figure_3 <- grid.arrange(d_eco_plot + theme(legend.position = 'none') + ylim(-4, 3.5), plot_chl_temp + ylim(-7.5, -1.5), rate_Chla_plot + ylim(-1, 2.75), rate_control_Chla_plot + ylim (-5, 4), ncol = 2)

ggsave(file.path(path_fig, 'Figure_3.pdf'), Figure_3, height = 11, width = 11)


# SI analysis ####

# 1. Analysis of nutrients vs temperature ####
EnvTraits <- select(metab_2016, c(temp, pH, conduct, depth, width, velocity, reaeration, NO2, NO2_NO3, NH4, PO4)) %>%
  mutate(NO3 = NO2_NO3 - NO2) %>%
  gather(., 'env.fac', 'value', c(pH, conduct, depth, width, velocity, reaeration, NO2, NO3, NH4, PO4)) %>%
  group_by(env.fac) %>%
  do(tidy(cor.test(.$value, .$temp, method = 'pearson'))) %>%
  data.frame()

# none of the correlations are significant




