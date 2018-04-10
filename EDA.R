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

ggplot(cleaned_dat %>% filter(Spine=='T5T6', ID=='10'), aes(angle, Dist2Aorta, color=SyndHeight))+
  geom_point()

cleaned_dat %>% group_by(ID, Spine) %>%
  mutate(HeightRank = rank(SyndHeight, ties.method='min')) %>%
  top_n(-5,Dist2Aorta) %>%
  mutate(meanHeightRank = mean(HeightRank)) %>%
  arrange(ID, Spine, Dist2Aorta)-> blah
