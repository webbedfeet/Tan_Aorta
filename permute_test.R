# Permutation testing

# Setup ---------------------------------------------------------------------------------------

ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)

munged_data <- readRDS('data/rda/cleaned_data_2.rds')


# Permutation testing -------------------------------------------------------------------------

x <- seq(-90,80,by=5)
wedges <- as.list(data.frame(rbind(x, x+5, x+10)))

`%nin%` <- Negate(`%in%`)

perm_test_fn_cts <- function(d){
  out <- tibble(wedge_centre = seq(-85,85,by=5))
  out$in_wedge <- lapply(wedges, function(x) d %>% filter(rel_angle %in% x) %>%
                           summarize(M = mean(SyndHeight))) %>% bind_rows() %>% pull(M)
  out$out_wedge <- lapply(wedges, function(x) d %>% filter(rel_angle %nin% x) %>% summarize(M = mean(SyndHeight))) %>% bind_rows() %>% pull(M)
  out <- out %>% mutate(difference = in_wedge-out_wedge)
  return(out)
}

perm_test_fn_bin <- function(d){
  out <- tibble(wedge_centre = seq(-85,85,by=5))
  out$in_wedge <- lapply(wedges, function(x) d %>% filter(rel_angle %in% x) %>%
                           summarize(M = mean(SyndPresent == 1))) %>% bind_rows() %>% pull(M)
  out$out_wedge <- lapply(wedges, function(x) d %>% filter(rel_angle %nin% x) %>%
                            summarize(M = mean(SyndPresent == 1))) %>% bind_rows() %>% pull(M)
  out <- out %>% mutate(difference = in_wedge-out_wedge)
  return(out)
}

outputs <- lapply(split(munged_data, munged_data$Spine), perm_test_fn_cts)
permute_pvals_cts <- map(outputs, ~mean(.$difference <= .$difference[.$wedge_centre==0])) %>%
  bind_rows() %>% gather(Spine, pval_cts)
outputs <- lapply(split(munged_data, munged_data$Spine), perm_test_fn_bin)
permute_pvals_bin <- map(outputs, ~mean(.$difference <= .$difference[.$wedge_centre==0])) %>%
  bind_rows() %>% gather(Spine, pval_bin)
dists <- munged_data %>% group_by(Spine) %>% summarize(Dist2Aorta = min(Dist2Aorta)) %>%
  mutate(Spine = as.character(Spine))
dists %>% left_join(permute_pvals_cts) %>% left_join(permute_pvals_bin)

pdf('graphs/heightCI.pdf')
for(nm in unique(munged_data$Spine)){
  d <- munged_data %>% filter(Spine==nm)
  bl <- d %>% group_by(rel_angle) %>%
    summarize(mn = mean(SyndHeight), sd = sd(SyndHeight), n = length(SyndHeight),
              sem = sd/sqrt(n),
              lcb = mn - 1.96*sem, ucb = mn + 1.96*sem)
  print(ggplot(bl, aes(x = rel_angle, y  = mn, ymin = lcb, ymax = ucb))+
    geom_pointrange()+ggtitle(nm) +
    labs(x = 'Relative angle from nadir (degrees)',
         y = 'Avg syndesmophyte height (mm)'))
}
dev.off()


# Alternatives --------------------------------------------------------------------------------

## Look at a ranking solution

d <- munged_data %>% filter(Spine=='T12L1')
d %>% group_by(ID) %>% mutate(ranking = rank(SyndHeight, ties.method = 'min')) %>% ungroup() -> bl

blah <- munged_data %>%
  group_by(Spine, ID) %>%
  mutate(ranking = rank(SyndHeight, ties.method='min')) %>%
  ungroup()
blah %>% count(Spine, rel_angle, ranking) %>%
  group_by(Spine, rel_angle) %>%
  mutate(perc = n/sum(n)*100) %>%
  ungroup() -> bl2
bl2 %>% filter(rel_angle==0, ranking==1)
bl2 %>% filter(rel_angle==45, ranking==1)
bl2 %>% filter(ranking == 1) %>%
  ggplot(aes(x = rel_angle, y = perc))+geom_line() +
  geom_vline(xintercept = 0, linetype = 2)+facet_wrap(~Spine)+
  xlab('Angle relative to nadir') + ylab('Percentage with growth ranking 1')


# Bootstrapping to evaluate ranking uncertainty -----------------------------------------------

library(rsample)
set.seed(154235)
orig <- bl2 %>% filter(ranking ==1) %>% split(.$Spine)
bs <- munged_data %>% nest(-Spine, -ID) %>% split(.$Spine) %>%
  map(bootstraps, times = 500)

prep_bs_data <- function(d){
  tmp <- map(d$splits, ~as.tibble(.) %>%
               mutate(ranking = map(data, ~rank(.$SyndHeight, ties.method='min')),
                      rel_angle = map(data, 'rel_angle')) %>%
               select(-data) %>% unnest() %>%
               count(rel_angle, ranking) %>%
               group_by(rel_angle) %>% mutate(perc = n/sum(n)*100) %>% ungroup() %>%
               filter(ranking == 1))
  out <- bind_rows(tmp, .id = 'boots')
  return(out)
}
bs_data <- map(bs, prep_bs_data)

pdf('MinRankPerc_Boots.pdf')
for(sp in names(bs_data)){
  plt <- ggplot() +
    geom_line(data = bs_data[[sp]], aes(x = rel_angle, y = perc, group = boots), alpha = 0.02)+
    geom_line(data = orig[[sp]], aes(x = rel_angle, y = perc), size = 2)+
    geom_vline(xintercept = 0, linetype = 2)+
    labs(x = 'Angle relative to nadir', y = 'Percentage with minimum rank')+
    ggtitle(sp)
  print(plt)
}
dev.off()


# Ranking for probabilty of growth, and bootstrapping -----------------------------------------

orig_probs <- munged_data %>% group_by(Spine, rel_angle) %>%
  summarize(prob = mean(SyndPresent)) %>% ungroup()

prep_bs_data_prob <- function(d){
  tmp <- map(d$splits, ~as.tibble(.) %>% unnest() %>%
               group_by(rel_angle) %>%
               summarize(prob = mean(SyndPresent)) %>% ungroup())
  out <- bind_rows(tmp, .id='boots')
  return(out)
}

bs_data_prob <- map(bs, prep_bs_data_prob)

pdf('Prob_Boots.pdf')
for(sp in names(bs_data_prob)){
  plt <- ggplot()+
    geom_point(data = bs_data_prob[[sp]], aes(x = rel_angle, y = prob), alpha = 0.01)+
    geom_line(data = orig_probs %>% filter(Spine==sp), aes(x = rel_angle, y = prob ))+
    geom_vline(xintercept = 0, linetype = 2)+
    labs(x = 'Angle from the nadir', y = 'Probability of having any growth')+
    ggtitle(sp)
  print(plt)
}
dev.off()

pdf('Prob_boots_ci.pdf')
for(sp in names(bs_data_prob)){
  bl <- bs_data_prob[[sp]] %>% group_by(rel_angle) %>% summarize(lcb = quantile(prob, 0.025),
                                                                 ucb = quantile(prob, 0.975)) %>%
    ungroup()
  plt <- ggplot()+
    geom_linerange(data = bl, aes(x = rel_angle, ymin=lcb, ymax = ucb), alpha=0.3)+
    geom_line(data = orig_probs %>% filter(Spine==sp), aes(x = rel_angle, y = prob ))+
    #geom_vline(xintercept = 0, linetype = 2)+
    labs(x = 'Angle from the nadir', y = 'Probability of having any growth')+
    ggtitle(sp)
  print(plt)
}
dev.off()


# Comparing different sectors (Bayesian) -----------------------------------------------------------------

## Let's compare the 0 angle to the +/- 15, 30 and 45 angles. We will evaluate the
## probability that any of these sectors have at least as much growth as the 0th sector.
## We will take a Bayesian approach using `rstanarm`

#library(rstanarm)
ProjTemplate::reload()
munged_data <- readRDS('data/rda/cleaned_data_2.rds')
expit <- function(x) exp(x)/(1+exp(x))

d <- munged_data %>% filter(rel_angle %in% c(0, 15, 30, 45, -15, -30, -45)) %>%
  select(Spine, ID, rel_angle, SyndHeight) %>% spread(rel_angle, SyndHeight) %>%
  gather(rel.angle, SyndHeight, -Spine, -ID, -`0`) %>%
  mutate(inds = ifelse(SyndHeight >= `0`, 1, 0)) %>%
  select(-`0`, -SyndHeight) %>%
  spread(rel.angle, inds)
d %>% group_by(Spine) %>% summarise_at(vars(`-15`:`45`), funs(mean))

# For uniform prior, posterior is Beta(1 + x, 1+n-x) (Conjugate prior)

d1 <- d %>% gather(rel.angle, inds, `-15`:`45`) %>% mutate(inds = as.factor(inds)) %>%
  nest(-Spine, -rel.angle) %>%
  mutate(all_pos = map_int(data, ~all(.$inds==1)),
         posterior.summary = map(data, ~data.frame(matrix(qbeta(c(0.025,0.5, 0.975), 1 + sum(.$inds=='1'), 1+ sum(.$inds=='0')), nrow=1))))

d2 <- d1 %>% select(Spine, rel.angle, posterior.summary) %>% unnest() %>%
  rename(lower.bound = X1, posterior.median = X2, upper.bound = X3)

# d2 <- d1 %>% filter(all_pos == 0) %>%
#   mutate(mods = map(data, ~stan_glm(inds~1, data = ., family=binomial(link='logit'),
#                                     prior_intercept = normal(0,1.75))))
# d3 <- d2 %>%
#   mutate(posterior.median = map_dbl(mods, ~expit(median(as.matrix(.)[,1]))),
#          lower.bound = map_dbl(mods, ~expit(quantile(as.matrix(.)[,1], 0.025))),
#          upper.bound = map_dbl(mods, ~expit(quantile(as.matrix(.)[,1], 0.975)))) %>%
#   select(-data, -mods)

N <- munged_data %>% group_by(Spine) %>% summarize(nobs = length(unique(ID))) %>% ungroup()
final_d <- d2 %>%  left_join(N)


final_d <- final_d %>% arrange(Spine, rel.angle) %>%
  mutate(rel.angle = factor(rel.angle, levels = c('-45','-30','-15','15','30','45')))

library(Cairo) # Need this for including Unicode characters
cairo_pdf('PosteriorProbs.pdf', width = 10, height = 7.5)
ggplot(final_d,
       aes(x = rel.angle, y = posterior.median, ymin = lower.bound, ymax = upper.bound)) +
  geom_pointrange() +
  facet_wrap(~Spine) +
  scale_x_discrete(limits = rev(levels(final_d$rel.angle))) +
  xlab("Relative angle from nadir") +
  ylab('Posterior 95% interval for P(Sector growth \u2265 Sector 0 growth)') +
  coord_flip() +
  theme_bw()
dev.off()


# Update Bayesian analysis, using strict less than criterion ----------------------------------

#library(rstanarm)
munged_data <- readRDS('data/rda/cleaned_data_2.rds')
expit <- function(x) exp(x)/(1+exp(x))

d <- munged_data %>% filter(rel_angle %in% c(0, 15, 30, 45, -15, -30, -45)) %>%
  select(Spine, ID, rel_angle, SyndHeight) %>% spread(rel_angle, SyndHeight) %>%
  gather(rel.angle, SyndHeight, -Spine, -ID, -`0`) %>%
  mutate(inds = ifelse(SyndHeight > `0`, 1, 0)) %>% # Strict inequality
  select(-`0`, -SyndHeight) %>%
  spread(rel.angle, inds)

# For uniform prior, posterior is Beta(1 + x, 1+n-x) (Conjugate prior)

d1 <- d %>% gather(rel.angle, inds, `-15`:`45`) %>% mutate(inds = as.factor(inds)) %>%
  nest(-Spine, -rel.angle) %>%
  mutate(all_pos = map_int(data, ~all(.$inds==1)),
         posterior.summary = map(data, ~data.frame(matrix(qbeta(c(0.025,0.5, 0.975), 1 + sum(.$inds=='1'), 1+ sum(.$inds=='0')), nrow=1))))

d2 <- d1 %>% select(Spine, rel.angle, posterior.summary) %>% unnest() %>%
  rename(lower.bound = X1, posterior.median = X2, upper.bound = X3)

N <- munged_data %>% group_by(Spine) %>% summarize(nobs = length(unique(ID))) %>% ungroup()
final_d <- d2 %>%  left_join(N)


final_d <- final_d %>% arrange(Spine, rel.angle) %>%
  mutate(rel.angle = factor(rel.angle, levels = c('-45','-30','-15','15','30','45')))

ggplot(final_d,
       aes(x = rel.angle, y = posterior.median, ymin = lower.bound, ymax = upper.bound)) +
  geom_pointrange() +
  facet_wrap(~Spine) +
  scale_x_discrete(limits = rev(levels(final_d$rel.angle))) +
  xlab("Relative angle from nadir") +
  ylab('Posterior 95% interval for P(Sector growth > Sector 0 growth)') +
  coord_flip() +
  theme_bw()


# Update Bayesian analysis, restricted to non-fused joints ------------------------------------

library(rstanarm)
munged_data <- readRDS('data/rda/cleaned_data_2.rds')
expit <- function(x) exp(x)/(1+exp(x))

d <- munged_data %>%
  group_by(Spine, ID) %>%
  mutate(not_fused = !(all(SyndHeight == 1))) %>% # Filter out fully fused joints
  ungroup() %>%
  filter(not_fused) %>%
  filter(rel_angle %in% c(0, 15, 30, 45, -15, -30, -45)) %>%
  select(Spine, ID, rel_angle, SyndHeight) %>% spread(rel_angle, SyndHeight) %>%
  gather(rel.angle, SyndHeight, -Spine, -ID, -`0`) %>%
  mutate(inds = ifelse(SyndHeight >= `0`, 1, 0)) %>%
  select(-`0`, -SyndHeight) %>%
  spread(rel.angle, inds)

# For uniform prior, posterior is Beta(1 + x, 1+n-x) (Conjugate prior)

d1 <- d %>% gather(rel.angle, inds, `-15`:`45`) %>% mutate(inds = as.factor(inds)) %>%
  nest(-Spine, -rel.angle) %>%
  mutate(all_pos = map_int(data, ~all(.$inds==1)),
         posterior.summary = map(data, ~data.frame(matrix(qbeta(c(0.025,0.5, 0.975), 1 + sum(.$inds=='1'), 1+ sum(.$inds=='0')), nrow=1))))

d2 <- d1 %>% select(Spine, rel.angle, posterior.summary) %>% unnest() %>%
  rename(lower.bound = X1, posterior.median = X2, upper.bound = X3)

N <- munged_data %>% group_by(Spine) %>% summarize(nobs = length(unique(ID))) %>% ungroup()
final_d <- d2 %>%  left_join(N)


final_d <- final_d %>% arrange(Spine, rel.angle) %>%
  mutate(rel.angle = factor(rel.angle, levels = c('-45','-30','-15','15','30','45')))

ggplot(final_d,
       aes(x = rel.angle, y = posterior.median, ymin = lower.bound, ymax = upper.bound)) +
  geom_pointrange() +
  facet_wrap(~Spine) +
  scale_x_discrete(limits = rev(levels(final_d$rel.angle))) +
  xlab("Relative angle from nadir") +
  ylab('Posterior 95% interval for P(Sector growth > Sector 0 growth)') +
  coord_flip() +
  theme_bw()+
  ggtitle('Among non-fused joints')

