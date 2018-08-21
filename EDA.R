##%######################################################%##
#                                                          #
####                        EDA                         ####
#                                                          #
##%######################################################%##

ProjTemplate::reload()

cleaned_dat <- readRDS('data/rda/cleaned_data_2.rds')

cleaned_dat %>% group_by(Spine) %>% summarize(NumberOfPts = length(unique(ID)))
# Overall association with height and distance ------------------------------------------------

ggplot(cleaned_dat, aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha=.1, size=0.5) +
  geom_smooth(aes(group = Spine, color = Spine), se=F, size=1.5)

ggplot(filter(cleaned_dat, Spine %in% c('T5T6','T6T7','T7T8','T8T9','T9T10','T10T11')),
       aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(aes(group = Spine, color = Spine), se=F, size=1.5)

ggplot(filter(cleaned_dat, Spine %in% c('T11T12','T12L1','L1L2','L2L3','L3L4')),
       aes(Dist2Aorta, SyndHeight))+
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(aes(group = Spine, color = Spine), se=F, size=1.5)


# Association by spine level ------------------------------------------------------------------

ggplot(cleaned_dat, aes(Dist2Aorta, SyndHeight))+geom_point(alpha = 0.1) + geom_smooth() +
  facet_wrap(~Spine) +
  labs(x = 'Distance from Aorta (mm)',
       y = 'Normalized syndesmophyte height')

sem <- function(x){sqrt(var(x, na.rm=T)/sum(!is.na(x)))}

cleaned_dat_summarized <- cleaned_dat %>%
  group_by(Spine, angle) %>%
  summarize(meanADist = mean(Dist2Aorta, na.rm=T), seADist = sem(Dist2Aorta),
            meanHt = mean(SyndHeight, na.rm=T), seHt = sem(SyndHeight),
            propSynd = mean(SyndPresent, na.rm=T)) %>%
  ungroup()

cleaned_dat_summarized %>% mutate(meanHt = meanHt*50) %>%
  ggplot(aes(x = angle))+
  geom_line(aes(y = meanADist, color='Distance'))+
  geom_line(aes(y = meanHt, color = 'Height'))+
  scale_y_continuous('Mean distance from aorta (mm)', sec.axis = sec_axis(~./50, name = 'Mean normalized height'))+
  facet_wrap(~Spine)+
  scale_color_manual(values = c('blue','red'))+
  labs(color = '', x ='Angle')

cleaned_dat_summarized %>% mutate(propSynd = 50*propSynd) %>%
  ggplot(aes(x = angle)) +
  geom_line(aes(y = meanADist, color='Distance')) +
  geom_line(aes(y = propSynd, color = 'Percent')) +
  scale_y_continuous('Mean distance from aorta (mm)',
                     sec.axis = sec_axis(~./50, name = 'Proportion present')) +
  facet_wrap(~Spine) +
  scale_color_manual(values = c('blue','red')) +
  labs(color = '', x = 'Angle')+
  theme(axis.text = element_text(size = 8, angle = 45))



# Ranking distances and looking for relationships ---------------------------------------------

cleaned_dat <- cleaned_dat %>% group_by(Spine, ID) %>%
  mutate(DistRank = rank(Dist2Aorta, ties.method = 'min')) %>% ungroup()
dat_summarize <- cleaned_dat %>% group_by(Spine, DistRank) %>%
  summarize(PctPos = mean(SyndPresent, na.rm=T),
            meanHt = mean(SyndHeight, na.rm=T)) %>%
  ungroup()
ggplot(dat_summarize, aes(x = DistRank, y = PctPos))+
  geom_point(alpha = 0.5)+
  geom_smooth() +
  facet_wrap(~Spine)+
  labs(x = 'Rank of each sector based on distance from aorta',
       y = 'Percent of individuals with postive growth')

ggplot(dat_summarize, aes(x = DistRank, y = meanHt))+
  geom_point(alpha = 0.5)+
  geom_smooth() +
  facet_wrap(~Spine)+
  labs(x = 'Rank of each sector based on distance from aorta',
       y = 'Avg normalized syndesmophyte height')

ggplot(dat_summarize, aes(x=DistRank, y = meanHt)) + geom_point() + geom_smooth(se = F)+
  facet_wrap(~Spine)

cleaned_dat %>% filter(DistRank==1) %>%
  ggplot(aes(x = Dist2Aorta, y = SyndHeight))+geom_point() + geom_smooth( se = F)
cleaned_dat %>% group_by(Spine) %>% summarize(minDist = mean(Dist2Aorta, na.rm=T)) %>% ungroup() %>%
  ggplot(aes(x = minDist, y = Spine))+geom_point() + geom_line(group = 1) +
  scale_y_discrete(limits = rev(levels(cleaned_dat$Spine)))

##%######################################################%##
#                                                          #
####                 With centered data                 ####
#                                                          #
##%######################################################%##


munged_data <- readRDS('data/rda/cleaned_data_2.rds')
