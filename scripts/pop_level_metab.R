# Padfield et al 2017 Ecology Letters
# Iceland Biofilm Analysis with NP included ####

# clear workspace
mise::mise(pkgs = TRUE, console = TRUE, figs = TRUE, vars = TRUE)

# load in packages ####
library(ggplot2) 
library(TeamPhytoplankton)
library(dplyr)
library(tidyr)
library(nlme)
library(nlsLoop)
library(grid)
library(gridExtra)
library(lme4)

# path to graphs ####
path_fig <- 'plots'

# functions to use ####

# function for removing legend from ggplot object for Figure 2
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

# function for layouts in grid
vplayout <- function(x, y)
  viewport(layout.pos.row = x, layout.pos.col = y)

# function computing AIC scores
all.AIC <- function(x){data.frame(model = x, AIC = MuMIn::AICc(get(x)))}

# load in data ####
d <- readRDS('data/pop_metab_data.rds') 

# run nlsLoop model ####
# does non-linear regression of every single curve based on est
# shotguns start values between parameter values and picks best model for each est by AIC score
# est = stream:species:GP/NP/R:stream.temp:year
res <- nlsLoop::nlsLoop(log_rate ~ schoolfield.high(ln.c, Ea, Eh, Th, temp = K, Tc = 10),
                      data = d,
                      tries = 1000,
                      id_col = 'est',
                      param_bds = c(-5,10,0.1,2,0.5,5,285,330),
                      r2 = 'Y',
                      supp.errors = 'Y',
                      AICc = 'Y',
                      na.action = na.omit,
                      lower = c(ln.c = -10, Ea=0, Eh=0, Th=0),
                      upper = c(ln.c = 100, Ea=30, Eh=40, Th=350))

# Transformations
# split columns of param_res$params and param_res$predictions
res$params <- separate(res$params, est, c('stream', 'species', 'flux', 'stream_temp', 'year'), sep = ':', remove = F)
res$predictions <- separate(res$predictions, est, c('stream', 'species', 'flux', 'stream_temp', 'year'), sep = ':', remove = F)
res$params$stream_temp <- as.numeric(res$params$stream_temp)
res$predictions$stream_temp <- as.numeric(res$predictions$stream_temp)

# plot all fits
rates <- unite(d, id, c(stream, species, year), sep = ':', remove = F)
res$params <- unite(res$params, id, c(stream, species, year), sep = ':', remove = F)
res$predictions <- unite(res$predictions, id, c(stream, species, year), sep = ':', remove = F)
# plot_all_nlsLoop(paste(path_fig, 'all_pop_level_fits.pdf', sep = '/'), rates, res, id_col = 'id', col_point = 'flux', col_line = 'flux')


# add columns to parameter dataset ####
# calculate Topt
# calculate (1/kTc) - (1/kT)
# calculate rate at stream temperature
res$params <- mutate(res$params, Topt = Topt(Eh, Th, Ea, K = "N"),
                     ikt_stream_temp = inv.temp(stream_temp, 'N', Tc = 10),
                     bTstream = schoolfield.high(ln.c, Ea, Eh, Th, temp = stream_temp + 273.15, Tc = 10))

# Metabolic trait analyses ####
# looking for differences of each metabolic trait across stream temperatures
# mixed effects model with ln_trait ~ E(1/kTc - 1/kT)*flux + intercept

# 1. ln.c ####
# normalisation constant with just flux and ikt.stream temp as factors
lncMM1 <- lmer(ln.c ~ ikt_stream_temp * flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
lncMM2 <- lmer(ln.c ~ ikt_stream_temp + flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
lncMM3 <- lmer(ln.c ~ ikt_stream_temp + (1|id), res$params, na.action = na.fail, REML = FALSE)
lncMM4 <- lmer(ln.c ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
anova(lncMM1, lncMM2)
anova(lncMM2, lncMM3)
anova(lncMM2, lncMM4)
plyr::ldply(ls(pattern = 'lncMM'), all.AIC)
# best model lncMM2
lncMMbest <- update(lncMM2, REML = TRUE)


# 2. Topt ####
# This is only trait not logged as an activation energy of increase does not make sense
ToptMM1 <- lmer(log(Topt) ~ ikt_stream_temp * flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
ToptMM2 <- lmer(log(Topt) ~ ikt_stream_temp + flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
ToptMM3 <- lmer(log(Topt) ~ ikt_stream_temp + (1|id), res$params, na.action = na.fail, REML = FALSE)
ToptMM4 <- lmer(log(Topt) ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
anova(ToptMM1, ToptMM2)
anova(ToptMM2, ToptMM3)
anova(ToptMM2, ToptMM4)
plyr::ldply(ls(pattern = 'ToptMM'), all.AIC)

# best model is ToptMM2
ToptMMbest <- update(ToptMM2, REML = TRUE)

# 3. Ea ####
EaMM1 <- lmer(log(Ea) ~ ikt_stream_temp * flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EaMM2 <- lmer(log(Ea) ~ ikt_stream_temp + flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EaMM3 <- lmer(log(Ea) ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EaMM4 <- lmer(log(Ea) ~ ikt_stream_temp + (1|id), res$params, na.action = na.fail, REML = FALSE)
EaMM5 <- lmer(log(Ea) ~ 1 + (1|id), res$params, na.action = na.fail, REML = FALSE)
anova(EaMM1, EaMM2)
anova(EaMM2, EaMM3)
anova(EaMM2, EaMM4)
plyr::ldply(ls(pattern = 'EaMM'), all.AIC)

# best model is EaMM5 with no predictors

# bootstrap confidence intervals
bE<- boot::boot(res$params$Ea, function(u,i) mean(u[i]), R = 999)
boot::boot.ci(bE, type = c("norm", "basic", "perc"))

# 4. Eh ####
EhMM1 <- lmer(log(Eh) ~ ikt_stream_temp * flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EhMM2 <- lmer(log(Eh) ~ ikt_stream_temp + flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EhMM3 <- lmer(log(Eh) ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
EhMM4 <- lmer(log(Eh) ~ ikt_stream_temp + (1|id), res$params, na.action = na.fail, REML = FALSE)
EhMM5 <- lmer(log(Eh) ~ 1 + (1|id), res$params, na.action = na.fail, REML = FALSE)
anova(EhMM1, EhMM2)
anova(EhMM2, EhMM3)
anova(EhMM2, EhMM4)
plyr::ldply(ls(pattern = 'EhMM'), all.AIC)

# best model is EhMM5 with no predictors

# bootstrapped confidence intervals
bEh <- boot::boot(res$params$Eh, function(u,i) mean(u[i]), R = 999)
boot::boot.ci(bEh, type = c("norm", "basic", "perc"))


# 5. Th ####
ThMM1 <- lmer(log(Th) ~ ikt_stream_temp * flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
ThMM2 <- lmer(log(Th) ~ ikt_stream_temp + flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
ThMM3 <- lmer(log(Th) ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
ThMM4 <- lmer(log(Th) ~ ikt_stream_temp + (1|id), res$params, na.action = na.fail, REML = FALSE)
ThMM5 <- lmer(log(Th) ~ 1 + (1|id), res$params, na.action = na.fail, REML = FALSE)
anova(ThMM1, ThMM2)
anova(ThMM2, ThMM3)
anova(ThMM2, ThMM4)
anova(ThMM4, ThMM5)
plyr::ldply(ls(pattern = 'ThMM'), all.AIC)

# best model is ThMM2
ThMMbest <- update(ThMM2, REML = TRUE)

# 6. rate at stream temperature - bTstream ####
bTstreamMM1 <- lmer(bTstream ~ ikt_stream_temp * flux + (1|id), random = ~ 1|id2, res$params, na.action = na.fail, REML = FALSE)
bTstreamMM2 <- lmer(bTstream ~ ikt_stream_temp + flux + (1|id), random = ~ 1|id2, res$params, na.action = na.fail, REML = FALSE)
bTstreamMM3 <- lmer(bTstream ~ flux + (1|id), res$params, na.action = na.fail, REML = FALSE)
bTstreamMM4 <- lmer(bTstream ~ ikt_stream_temp + (1|id), random = ~ 1|id2, res$params, na.action = na.fail, REML = FALSE)
anova(bTstreamMM1, bTstreamMM2)
anova(bTstreamMM2, bTstreamMM3)
anova(bTstreamMM2, bTstreamMM4)
plyr::ldply(ls(pattern = 'bTstreamMM'), all.AIC)

# best model is bTstreamMM3
bTstreamMMbest <- update(bTstreamMM3, REML = TRUE)

# Figure 2 from manuscript ####

# ln.c plot (Fig 2c) #### 
# create predictions from lncMMfinal predictions
# refit in lme to get confidence intervals around predictions
lncMMbest <- lme(ln.c ~ ikt_stream_temp + flux, random = ~1|id, na.action = na.fail, method = 'REML', res$params)
lnc_preds <- data.frame(expand.grid(ikt_stream_temp = seq(min(res$params$ikt_stream_temp), max(res$params$ikt_stream_temp), length.out = 20), flux = c('GP', 'NP', 'R'))) %>%
  mutate_at(., c('flux'), as.character) %>%
  mutate(., stream_temp = rep(seq(min(res$params$stream_temp), max(res$params$stream_temp), length.out = 20), times = 3),
         mu = as.vector(predict(lncMMbest, level = 0, .)))

# create design matrix of model
Designmat <- model.matrix(formula(lncMMbest)[-2], lnc_preds)

# compute standard error for predictions
lnc_preds$predvar <- diag(Designmat %*% vcov(lncMMbest) %*% t(Designmat))
lnc_preds <- mutate(lnc_preds,
                    lwr_CI = mu - sqrt(predvar),
                    upr_CI = mu + sqrt(predvar),
                    lwr_PI = mu - sqrt(predvar + lncMMbest$sigma^2),
                    upr_PI = mu + sqrt(predvar + lncMMbest$sigma^2))

# plot
plot_lnc <- ggplot() +
  # geom_ribbon(aes(x = stream_temp, ymin = lwr_CI, ymax = upr_CI, fill = flux), alpha = 0.1, lnc_preds) +
  geom_point(aes(x = stream_temp, y = ln.c, col = flux), size = 1.5, res$params) +
  geom_line(aes(x = stream_temp, y = mu, col = flux), lnc_preds) +
  xlab(expression(Stream~Temperature~(degree*C))) +
  ylab(expression(Temperature~Normalised~Rate:~ln~italic(b(T[c])))) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.95,0.7), legend.justification=c(1,0)) +
  scale_colour_manual(values = c('green4', 'blue', 'red4')) +
  scale_fill_manual(values = c('green4', 'blue', 'red4')) +
  ylim(-1, 6)

# bTstream plot (Fig 2d) ####
# refit bTstreamMMbest in lme
bTstreamMMbest <- lme(bTstream ~ flux, random = ~ 1|id, res$params, na.action = na.fail, method = 'REML')

# predictions
bTs_preds <- data.frame(expand.grid(ikt_stream_temp = seq(min(res$params$ikt_stream_temp), max(res$params$ikt_stream_temp), length.out = 20), flux = c('GP', 'NP', 'R'))) %>%
  mutate_at(., c('flux'), as.character) %>%
  mutate(., stream_temp = rep(seq(min(res$params$stream_temp), max(res$params$stream_temp), length.out = 20), times = 3),
         mu = as.vector(predict(bTstreamMMbest, level = 0, .)))

# create design matrix of model
Designmat <- model.matrix(formula(bTstreamMMbest)[-2], bTs_preds)

# compute standard error for predictions
bTs_preds$predvar <- diag(Designmat %*% vcov(bTstreamMMbest) %*% t(Designmat))
bTs_preds <- mutate(bTs_preds,
                    lwr_CI = mu - sqrt(predvar),
                    upr_CI = mu + sqrt(predvar),
                    lwr_PI = mu - sqrt(predvar + bTstreamMMbest$sigma^2),
                    upr_PI = mu + sqrt(predvar + bTstreamMMbest$sigma^2))

plot_bTstream <- ggplot() +
  geom_point(aes(x = stream_temp, y = bTstream), size = 1.5, filter(res$params, flux == "GP")) +
  geom_line(aes(stream_temp, mu), linetype = 2, filter(bTs_preds, flux == 'GP')) +
  geom_ribbon(aes(x = stream_temp, ymin = lwr_CI, ymax = upr_CI), fill = 'green4', alpha = 0.2, filter(bTs_preds, flux == 'GP')) +
  labs(x=expression(Stream~Temperature~(degree*C)), y=expression(GP~at~Stream~Temperature:~ln~italic(gp(T[s])))) +
  theme_bw(base_size = 12, base_family = 'Helvetica')

# plot net photosynthesis (Fig 2a) ####
plot_NP <- ggplot(filter(rates, flux == 'GP')) +
  geom_line(aes(x= K - 273.15, y=log_rate, colour=stream_temp, group = est), filter(res$predictions, flux == 'GP')) +
  geom_point(aes(temp, log_rate, colour=stream_temp, group = est), alpha = 0.25) +
  scale_colour_gradient(breaks = c(8,10, 12, 14, 16, 18, 20, 22, 24, 26), low = 'blue', high = 'red', 'Stream \nTemperature (ºC)') +
  labs(x=expression(Assay~Temperature~(degree*C)), y=expression(ln~Gross~Photosynthesis~(mu*mol~O[2]~mu*g~Chl*italic(a)^{-1}~h^{-1}))) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(strip.background = element_rect(fill = 'white', colour = 'white'),
        legend.key.size = unit(1.25, "cm"))+
  xlim(0,65) +
  ylim(-3, 8)

# plot respiration (Fig 2b) ####
plot_R <- ggplot(filter(rates, flux == 'R')) +
  geom_line(aes(x= K - 273.15, y=log_rate, colour=stream_temp, group = est), filter(res$predictions, flux == 'R')) +
  geom_point(aes(temp, log_rate, colour=stream_temp, group = est), alpha = 0.25) +
  scale_colour_gradient(breaks = c(8,10, 12, 14, 16, 18, 20, 22, 24, 26), low = 'blue', high = 'red', 'Stream \nTemperature (ºC)') +
  labs(x=expression(Assay~Temperature~(degree*C)), y=expression(ln~Respiration~(mu*mol~O[2]~mu*g~Chl*italic(a)^{-1}~h^{-1}))) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  #facet_wrap(~ flux, labeller = labeller(.multi_line = FALSE, flux = facet_labels)) +
  theme(strip.background = element_rect(fill = 'white', colour = 'white'),
        legend.key.size = unit(1.25, "cm"))+
  xlim(0,65) +
  ylim(-3, 8)

# assemble plot ####
# simple method
gridExtra::grid.arrange(plot_NP, plot_R, plot_lnc, plot_bTstream, ncol = 2)

# viewport method ####
grid.Fig1 <- viewport(layout = grid.layout(nrow = 2, ncol = 3, heights = c(0.6, 0.4), widths = c(0.4, 0.4, 0.2)))
grid.NP <- viewport(layout.pos.col = 1, layout.pos.row = 1, name = 'grid.NP')
grid.R <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = 'grid.R')
grid.ln.c <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = 'grid.lnc')
grid.bTs <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = 'grid.bTs')
grid.leg <- viewport(layout.pos.col = 3, layout.pos.row = 1:2, name = 'grid.leg')
splot <- vpTree(grid.Fig1, vpList(grid.NP, grid.R, grid.ln.c, grid.bTs, grid.leg))

# extract legend from plot
leg <- g_legend(plot_NP)

# assemble plot ####
pdf(file.path(path_fig, 'Figure_2.pdf'), width = 11, height = 9)

grid.newpage()
pushViewport(splot)
seekViewport('grid.NP')
print(plot_NP + theme(legend.position = 'none'), vp = vplayout(1,1))
seekViewport('grid.R')
print(plot_R + theme(legend.position = 'none'), vp = vplayout(1,2))
seekViewport('grid.lnc')
print(plot_lnc + theme(legend.position = 'none') + ylim(-2, 6), vp = vplayout(2,1))
seekViewport('grid.bTs')
print(plot_bTstream + ylim(1,5.5), vp = vplayout(2,2))
seekViewport('grid.leg')
grid.draw(leg)

dev.off()
