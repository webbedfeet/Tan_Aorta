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
#'    repeated mea
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
