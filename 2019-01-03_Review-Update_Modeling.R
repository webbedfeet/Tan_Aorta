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
#'    repeated measures approach with sectors as factors. We can evaluate the statistical uncertainty
#'    using cluster bootstrap

#+ setup, include = FALSE
# Setup ---------------------------------------------------------------------------------------

ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)

# system(paste('rsync -azvh ', normalizePath(file.path('data/rda')),' ', datadir))
munged_data <- readRDS(file.path(datadir,'rda/cleaned_data_2.rds'))

# Recomputing sector distances/labels ----------------------------------------------------------------

munged_data_2 <- munged_data %>%
  mutate(d0 = abs(rel_angle/5), # abs distance from nadir
         d1 = (rel_angle + 90)/5, # counting from left
         d2 = (rel_angle - 90)/5) # counting from right


#+ pattern_plots1, eval = TRUE, echo = FALSE, message = FALSE
# Plot relation between syndesmophyte height and sectors from nadir ---------------------------

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
