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



cleaned_dat %>% filter(Spine == 'T10T11') %>%
  ggplot(aes(x = Dist2Aorta,y = SyndHeight, group = ID))+
  geom_smooth(se=F)

cleaned_dat %>% filter(Spine == 'T10T11', ID == '10') %>%
  ggplot(aes(x = Dist2Aorta, y = SyndHeight))+
  geom_line()


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
