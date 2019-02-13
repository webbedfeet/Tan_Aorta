ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
scientific_10 <- function(x) {
  parse(text=gsub("1e", "10^", scales::scientific_format()(x)))
}


# Figure 4 ------------------------------------------------------------------------------------

m8 <- readRDS(file.path(datadir, 'FinalGEEModel.rds'))
L <- cbind(matrix(0,36, 48), diag(rep(1,36))) # Contrast matrix
L <- rbind(L[1:18,],rep(0,84),L[19:36,])
L[,38] = 1 # Dist2Aorta location

bl <- doBy::esticon(m8, L)
labs <- str_remove(names(coef(m8)[49:84]),':Dist2Aorta')
labs <- c(labs[1:18], 'rel_angle0', labs[19:36])
bl <- cbind('Sector' =labs, bl) %>%
  mutate(Sector = str_remove(Sector, 'rel_angle')) %>%
  select(Sector, Estimate, Std.Error,`Pr(>|X^2|)`, Lower, Upper) %>%
  rename(pval = `Pr(>|X^2|)`)

pval_plt <- bl %>%
  mutate(Sector = as.numeric(Sector)) %>%
  ggplot(aes(Sector, pval))+geom_point()+
  geom_hline(yintercept = 0.05, linetype =3)+
  scale_y_log10('P-value', label=scientific_10,
                breaks = c(0.05, 1e-3, 1e-7, 1e-11))+
  scale_x_continuous('Sector angle') +
  theme_bw()

or_plt <- bl %>%
  mutate(Sector = as.numeric(Sector)) %>%
  mutate(OR = exp(Estimate),
         Lower = exp(Estimate - 1.96*Std.Error),
         Upper = exp(Estimate + 1.96*Std.Error)) %>%
  ggplot(aes(x = Sector, y = OR, ymin = Lower, ymax = Upper))+
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_x_continuous('Sector angle') +
  scale_y_continuous('Odds ratio') +
  theme_bw()
