#' ---
#' title: "Addressing paper reviews"
#' author: "Abhijit Dasgupta"
#' date: "`r format(Sys.Date(), '%B %d, %Y')`"
#' output: html_document
#' ---


#' The paper reviews promoted an idea of modeling with the "complete cohort". We are
#' interpreting that to mean that we want to provide overall p-values showing that the effect
#' of the aorta on syndesmophyte formation (by showing that overall the nadir sector has
#' the lowest height of syndesmophytes) exists, rather than being uniform.
#'
#' A couple of ideas will be tried.
#'
#' 1. Show an upward trend in syndesmophyte heights moving away from the nadir sector, and
#'    that any other ordering of the sectors would not show a generally increasing pattern
#' 2. Model the avg heights and/or the presence of syndesmophytes on location using a
#'    repeated measures approach
#'
#' 2019-01-08
#' ===========
#'
#' We're doing GEE to look at the sector effect, adjusting for the joint and the
#' distance from the aorta of the sector
#'

# Setup ---------------------------------------------------------------------------------------

ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)


munged_data <- readRDS(file.path(datadir,'rda','cleaned_data_2.rds'))
munged_data %<>% mutate( rel_angle = as.factor(rel_angle)) %>%
  mutate(rel_angle = fct_relevel(rel_angle, '0')) %>%
  arrange(ID)


# Recomputing sector distances/labels ----------------------------------------------------------------

munged_data_2 <- munged_data %>%
  mutate(d0 = abs(rel_angle/5), # abs distance from nadir
         d1 = (rel_angle + 90)/5, # counting from left
         d2 = (rel_angle - 90)/5) # counting from right


#+ pattern_plots1, eval = TRUE, echo = FALSE, message = FALSE
# Plot relation between syndesmophyte height and sectors from nadir ---------------------------
library(cowplot)
ggplot(munged_data_2, aes(x = d0, y = SyndHeight)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = 'Smoother'), se = F) +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F) +
  facet_wrap(~Spine) +
  labs(x = 'Sectors from nadir', y = 'Syndesmophyte height')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue','Linear' = 'red'))+
  theme(legend.position = 'bottom')

ggplot(munged_data_2, aes(x = d1, y = SyndHeight)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = 'Smoother'), se = F) +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F) +
  facet_wrap(~Spine) +
  labs(x = 'Sectors from left', y = 'Syndesmophyte height')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue','Linear' = 'red'))+
  theme(legend.position = 'bottom')

ggplot(munged_data_2, aes(x = abs(d2), y = SyndHeight)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = 'Smoother'), se = F) +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F) +
  facet_wrap(~Spine) +
  labs(x = 'Sectors from right', y = 'Syndesmophyte height')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue','Linear' = 'red'))+
  theme(legend.position = 'bottom')

#+ pattern_plots2, echo = FALSE, eval = FALSE
# Plot relation between probability of syndesmophyte and sector location ----------------------

munged_data_2 %>% group_by(Spine, d0) %>%
  summarize(p1 = mean(SyndHeight > 0), p0 = mean(SyndHeight == 0)) %>%
  ggplot(aes(x = d0, y = p0)) +
    geom_point() +
    geom_smooth(aes(color = 'Linear'), method = 'lm', se = F)+
    geom_smooth(aes(color = 'Smoother'), se = F)+
    facet_wrap(~Spine) +
    labs(x = 'Sectors from nadir', y = 'Probability of no syndesmophyte growth')+
    scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue', 'Linear' = 'red'))+
    theme(legend.position = 'bottom')

munged_data_2 %>% group_by(Spine, d1) %>%
  summarize(p1 = mean(SyndHeight > 0), p0 = mean(SyndHeight == 0)) %>%
  ggplot(aes(x = d1, y = p0)) +
  geom_point() +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F)+
  geom_smooth(aes(color = 'Smoother'), se = F)+
  facet_wrap(~Spine) +
  labs(x = 'Sectors from "left"', y = 'Probability of no syndesmophyte growth')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue', 'Linear' = 'red'))+
  theme(legend.position = 'bottom')

munged_data_2 %>% group_by(Spine, d2) %>%
  summarize(p1 = mean(SyndHeight > 0), p0 = mean(SyndHeight == 0)) %>%
  ggplot(aes(x = abs(d2), y = p0)) +
  geom_point() +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F)+
  geom_smooth(aes(color = 'Smoother'), se = F)+
  facet_wrap(~Spine) +
  labs(x = 'Sectors from "right"', y = 'Probability of no syndesmophyte growth')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue', 'Linear' = 'red'))+
  theme(legend.position = 'bottom')

#' It is evident from the plots that the approach Mike was suggesting won't work to distinguish the nadir.
#' We will advance two approaches:
#'
#' 1. Create the permutation distribution for the slope, by permuting the origin point of the sector distances
#' and computing the linear slope of the heights / the log-odds ratio of the trend, and then compare that distribution with when the nadir sector is the origin
#' 1. Use GEE to model the full dataset, first as a main effects model with joint as a covariate, then using the joint x sector interaction, using absense of growth as the outcome. Using the heights leads to issues of distribution and appropriateness of p-values and standard errors due to issues with the non-standard nature of the distribution.
#' 1. If we see signficant interactions, we can then report joint-by-joint models and their p-values, to better understand what is going on.


# GEE modeling --------------------------------------------------------------------------------

library(geepack)

m0 <- geeglm(SyndPresent ~ rel_angle, data = munged_data,
             id = ID,
             family = binomial(link=logit),
             corstr = 'ex',
             std.err = 'san.se')
m1 <- update(m0, . ~ . + Dist2Aorta)
m2 <- update(m0, . ~ . + Spine)
m3 <- update(m0, . ~ . + Spine + Dist2Aorta)
m4 <- update(m0, . ~ . + rel_angle*Spine)
m5 <- update(m0, . ~ . + rel_angle*Dist2Aorta)
m6 <- update(m0, . ~ . + Spine*Dist2Aorta)
m7 <- update(m0, . ~ . + rel_angle*Spine + Dist2Aorta)
m8 <- update(m0, . ~ . + rel_angle*Dist2Aorta + Spine)

anova(m1, m0)
#' Analysis of 'Wald statistic' Table
#'
#' Model 1 SyndPresent ~ rel_angle + Dist2Aorta
#' Model 2 SyndPresent ~ rel_angle
#' Df    X2 P(>|Chi|)
#' 1  1 0.882      0.35

anova(m2, m0)
#' Analysis of 'Wald statistic' Table
#'
#' Model 1 SyndPresent ~ rel_angle + Spine
#' Model 2 SyndPresent ~ rel_angle
#' Df   X2 P(>|Chi|)
#' 1 10 17.7     0.061 .
#' ---
#'   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m3, m0)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Spine + Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle
#` Df   X2 P(>|Chi|)
#` 1 11 25.2    0.0085 **
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m3, m1)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Spine + Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle + Dist2Aorta
#` Df X2 P(>|Chi|)
#` 1 10 24    0.0075 **
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m3, m2)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Spine + Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle + Spine
#` Df   X2 P(>|Chi|)
#` 1  1 7.33    0.0068 **
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m4, m2)
#` Models not nested
#` NULL

anova(m5, m1)
#` Models not nested
#` NULL

anova(m6, m3)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Spine + Dist2Aorta + Spine:Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle + Spine + Dist2Aorta
#` Df   X2 P(>|Chi|)
#` 1 10 28.9    0.0013 **
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m7, m4)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Spine + Dist2Aorta + rel_angle:Spine
#` Model 2 SyndPresent ~ rel_angle + Spine + rel_angle:Spine
#` Df   X2 P(>|Chi|)
#` 1  1 10.5    0.0012 **
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m8, m5)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Dist2Aorta + Spine + rel_angle:Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle + Dist2Aorta + rel_angle:Dist2Aorta
#` Df   X2 P(>|Chi|)
#` 1 10 31.4   0.00051 ***
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m7, m3)

anova(m8, m3)
#` Analysis of 'Wald statistic' Table
#` `
#` Model 1 SyndPresent ~ rel_angle + Dist2Aorta + Spine + rel_angle:Dist2Aorta
#` Model 2 SyndPresent ~ rel_angle + Spine + Dist2Aorta
#` Df  X2 P(>|Chi|)
#` 1 36 404    <2e-16 ***
#` `  ---
#` `  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#' Effect of distance by sector

L <- cbind(matrix(0,36, 48), diag(rep(1,36))) # Contrast matrix
L <- rbind(L[1:18,],rep(0,84),L[19:36,])
L[,38] = 1 # Dist2Aorta location

bl <- doBy::esticon(m8, L)
labs <- str_remove(names(coef(m8)[49:84]),':Dist2Aorta')
labs <- c(labs[1:18], 'rel_angle0', labs[19:36])
bl <- cbind('Sector' =labs, bl) %>%
  mutate(Sector = str_remove(Sector, 'rel_angle')) %>%
  select(Sector, Estimate, Std.Error,`Pr(>|X^2|)`, Lower, Upper)


# GEE with discrete aorta distances -----------------------------------------------------------
#' We will repeat this process with discretized distance from aorta, cut at <2 and <3 mm

munged_data_2 <- munged_data %>%
  mutate(dist2aorta_discrete = cut(round(Dist2Aorta,2),c(0, 2, 3, 100), labels = c('< 2', '\u2265 2 & <3', '\u2265 3'), right = F))

m8.discrete <- geeglm(formula = SyndPresent ~ rel_angle + dist2aorta_discrete + Spine +
                        rel_angle:dist2aorta_discrete, family = binomial(link = logit), data = munged_data_2,
                      id = ID, corstr = "ex", std.err = "san.se")
#' This is problematic since many of the rel_angle:dist2aorta_discrete combos have no data


# Description of sector x height frequencies --------------------------------------------------

blah <- munged_data_2 %>% mutate(rel_angle = as.numeric(as.character(rel_angle))) %>%
  count(rel_angle, dist2aorta_discrete) %>%
  spread(dist2aorta_discrete, n) %>%
  gather(Height, Count, -rel_angle) %>%
  mutate(Height = as.factor(Height))
lattice::cloud(Count ~ rel_angle * Height , data = blah,
               scales = list(arrows=F), zlab = list(z = 'Count', rot=90))

blah2 <- munged_data_2 %>% mutate(rel_angle = as.numeric(as.character(rel_angle))) %>%
  group_by(rel_angle, dist2aorta_discrete) %>% summarize(prob = mean(SyndPresent))
lattice::cloud(prob ~ rel_angle*dist2aorta_discrete, data = blah2,
               xlab = 'Sector', ylab = 'Distance', zlab = list('Probability of growth', rot=90),
               scales = list(arrows = F))
blah3 <- munged_data_2 %>% mutate(rel_angle = as.numeric(as.character(rel_angle))) %>%
  group_by(Spine, rel_angle, dist2aorta_discrete) %>% summarize(prob = mean(SyndPresent))
lattice::cloud(prob ~ rel_angle*dist2aorta_discrete | Spine, data = blah3,
               xlab = '', ylab = '', zlab = '',
               scales = list(arrows = F))
