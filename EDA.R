##%######################################################%##
#                                                          #
####                        EDA                         ####
#                                                          #
##%######################################################%##

ProjTemplate::reload()

cleaned_dat <- readRDS('data/rda/cleaned_data.rds')

cleaned_dat %>% group_by(Spine) %>% summarize(NumberOfPts = length(unique(ID)))
# Overall association with height and distance ------------------------------------------------

ggplot(cleaned_dat, aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha=.1, size=0.5) + geom_smooth()


ggplot(cleaned_dat, aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha=0.1, size=0.5)+
  geom_smooth(aes(color=Spine), se=F)

ggplot(cleaned_dat, aes(angle, Dist2Aorta))+
  geom_point(alpha=0.1)+geom_smooth(aes(color = Spine))

ggplot(cleaned_dat, aes(angle, SyndHeight))+
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(aes(color=Spine))

ggplot(cleaned_dat %>% filter(Spine=='T5T6'), aes(Dist2Aorta, SyndHeight))+geom_point() +
  facet_wrap(~ID)


cleaned_dat %>% group_by(ID, Spine) %>%
  mutate(height_rank = rank(SyndHeight, ties.method = 'min')) %>%
  top_n(-5, Dist2Aorta) %>% ungroup() %>% count(Spine, angle)
