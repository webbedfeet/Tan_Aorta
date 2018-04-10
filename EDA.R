##%######################################################%##
#                                                          #
####                        EDA                         ####
#                                                          #
##%######################################################%##

ProjTemplate::reload()

cleaned_dat <- readRDS('data/rda/cleaned_data.rds')

ggplot(cleaned_dat, aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha = 0.1, size = 0.5)+
  facet_wrap(~Spine)
