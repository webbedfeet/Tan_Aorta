# Permutation testing

ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)

munged_data <- readRDS('data/rda/cleaned_data_2.rds')
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

