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


# Plot relation between syndesmophyte height and sectors from nadir ---------------------------

ggplot(munged_data_2, aes(x = d0, y = SyndHeight)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = 'Smoother'), se = F) +
  geom_smooth(aes(color = 'Linear'), method = 'lm', se = F) +
  facet_wrap(~Spine) +
  labs(x = 'Sectors from nadir', y = 'Syndesmophyte height')+
  scale_color_manual(name = 'Legend', values = c('Smoother' = 'blue','Linear' = 'red'))+
  theme(legend.position = 'bottom')
