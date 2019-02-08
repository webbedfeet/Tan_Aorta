#' ---
#' title: Adding calcification data
#' author: Abhijit Dasgupta
#' date: 2019-01-19
#' output: word_document
#' ---
#'
#+ setup, include = FALSE
# Setup ---------------------------------------------------------------------------------------

ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)


munged_data <- readRDS(file.path(datadir,'rda','cleaned_data_2.rds'))
munged_data %<>% mutate( rel_angle = as.factor(rel_angle)) %>%
  mutate(rel_angle = fct_relevel(rel_angle, '0')) %>%
  mutate(Spine = fct_relevel(Spine, 'T12L1')) %>%
  arrange(ID)

#+ munging, include = FALSE
# Prepping the calcification data -------------------------------------------------------------

# calcification <- readxl::read_excel('data/raw/calcAnon.xlsx')
calcification <- readxl::read_excel(file.path(datadir, 'calcAnon.xlsx'))
calcification %<>%
  rename(ID = `Patient nb/Level`) %>%
  mutate_all(na_if, 'N/A') %>%
  gather(Spine, calc, -ID) %>%
  mutate(calc = as.numeric(calc))
munged_data <- munged_data %>%
  left_join(calcification, by = c('ID', 'Spine'))

#' The following table shows the results of adding calcification (as a 0-1 variable) to the GEE
#' model that models the relation between presence of syndesmophytes and sector angle, distance from aorta,
#' IDS, and the sector angle x distance from aorta interaction.
#'
#+ GEE, echo = FALSE, message = F, warning = F, include=F, eval=F
# GEE modeling --------------------------------------------------------------------------------

library(geepack)

m0 <- geeglm(SyndPresent ~ rel_angle*Dist2Aorta + Spine + calc,
             data = munged_data %>% mutate(Dist2Aorta = -Dist2Aorta), # Protective as you get closer
             id = ID,
             family = binomial(link = logit),
             corstr = 'ex', # Exchangeable working  correlation
             std.err = 'san.se') # Robust standard errors

broom::tidy(m0) %>%
  mutate(term = str_replace(term, 'rel_angle', 'Sector angle (per 5 degrees)')) %>%
  mutate(term = str_replace(term, 'Spine', 'Spinal joint = ')) %>%
  mutate(term = str_replace(term, 'Dist2Aorta','Aorta distance (mm)')) %>%
  mutate(term = str_replace(term, 'calc', 'Calcification')) %>%
  mutate(term = str_replace(term, ':', ' x ')) %>%
  filter(term != '(Intercept)') %>%
  mutate(OR = round(exp(estimate), 2)) %>%
  mutate(CI = as.character(glue::glue('({round(exp(estimate - 1.96*std.error),2)}, {round(exp(estimate + 1.96 * std.error),2)})'))) %>%
  mutate(`P-value` = str_replace(ProjTemplate::myformat_pvalue(p.value, digits =3), 'e(-\\d{2})',' x 10^\\1^')) %>%
  dplyr::select(term, OR:`P-value`) %>%
  knitr::kable()

#' #### Understanding the science
#'
#' Aortal calcification is being used as a surrogate for inflammation, with the idea that
#' proximity to calcified areas would promote inflammation in the close IDS, leading to
#' the promotion of syndesmophyte growth, since the growth is believed to be inflammation-
#' related. So presence of calcification would be a promotor of growth _biochemically_,
#' while distance of the aorta would be an inhibitor of growth _biomechanically_.
#'
#+ calc-interaction, echo=F, message=F, warning=F, include = F, eval=F

m1 <- update(m0, .~. + calc:Dist2Aorta)

munged_data <- readRDS(file.path(datadir,'rda','cleaned_data_2.rds'))
munged_data %<>% mutate( rel_angle = as.factor(rel_angle)) %>%
  mutate(rel_angle = fct_relevel(rel_angle, '0')) %>%
  mutate(Spine = fct_relevel(Spine, 'T12L1')) %>%
  arrange(ID)
munged_data <- munged_data %>%
  left_join(calcification, by = c('ID', 'Spine'))

m8 <- readRDS(file.path(datadir, 'FinalGEEModel.rds'))

m9 <- update(m8, .~. + calc*Dist2Aorta)

m10 <- update(m8, . ~ . + calc*rel_angle)

m11 <- update(m8, . ~ . +  calc : rel_angle : Dist2Aorta)

#' Let's look to see how calcification affects the effect of distance
#' We'll only look at T11-T12 through L3-L4 since there is no calcification higher on
#' the spine
#+ calc_graph, echo = F, message = F, warning = F
m8 <- readRDS(file.path(datadir, 'FinalGEEModel.rds'))

munged_data2 <- munged_data %>% filter(Spine %in% c('T11T12','T12L1','L1L2','L2L3','L3L4')) %>%
  mutate(Dist2Aorta = -Dist2Aorta)

munged_data3 <- munged_data %>% mutate(dist_cat = cut(Dist2Aorta, c(0,2,3,50), right=F))


L <- cbind(matrix(0,36, 42), diag(rep(1,36))) # Contrast matrix
L <- rbind(L[1:18,],rep(0,78),L[19:36,])
L[,38] = 1 # Dist2Aorta location


m8_1 <- update(m8, data = munged_data2 %>% filter(calc==1))
coef_1 <- doBy::esticon(m8_1, L)
labs <- str_remove(names(coef(m8_1)[43:78]),':Dist2Aorta')
labs <- c(labs[1:18], 'rel_angle0', labs[19:36])
coef_1 <- cbind('Sector' =labs, coef_1) %>%
  mutate(Sector = str_remove(Sector, 'rel_angle')) %>%
  dplyr::select(Sector, Estimate, Std.Error,`Pr(>|X^2|)`, Lower, Upper)

m8_0 <- update(m8, data = munged_data2 %>% filter(calc==0))
coef_0 <- doBy::esticon(m8_0,L)
coef_0 <- cbind('Sector' =labs, coef_0) %>%
  mutate(Sector = str_remove(Sector, 'rel_angle')) %>%
  dplyr::select(Sector, Estimate, Std.Error,`Pr(>|X^2|)`, Lower, Upper)

coefs <- bind_rows(list('Yes' = coef_1, 'No' = coef_0), .id = 'calc') %>%
  as_tibble() %>%
  mutate(x = as.numeric(Sector) - 1 + 2*ifelse(calc=='Yes',1,0)) %>%
  mutate_at(vars(Estimate, Lower, Upper), exp)
ggplot(coefs, aes(x = x, y = Estimate, ymin = Lower, ymax = Upper, color = calc))+
  geom_pointrange() +
  geom_hline(yintercept = 1) +
  scale_color_discrete('Aortic calcification')+
  scale_x_continuous('Sector') +
  scale_y_continuous('Odds ratio (95% CI) for
min AV distance per 1 mm reduction') +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('graphs/calcification.pdf')

