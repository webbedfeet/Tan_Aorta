ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)

# dat <- readRDS('data/rda/raw_data.rds')
#
# cleanData <- function(d) {
#   angle <- d %>% filter(row >= 4, col == 1) %>% pull(numeric)
#   name <- d %>% filter(row %in% c(1,3), col >= 3) %>%
#     mutate(character = ifelse(is.na(character), as.character(numeric),
#                               character)) %>%
#     select(row, col, character) %>%
#     filter(!is.na(character)) %>%
#     spread(row, character) %>%
#     unite(ID, `1`:`3`, sep = '_') %>%
#     pull(ID)
#   data <- d %>% filter(row >= 4, col >= 3) %>%
#     filter(!is.na(numeric)) %>%
#     select(row, col, numeric) %>%
#     spread(col, numeric) %>%
#     select(-row) %>%
#     set_names(name) %>%
#     cbind(angle) %>%
#     gather(variable, value, -angle) %>%
#     separate(variable, c('ID','Measure'), sep = '_') %>%
#     spread(Measure, value) %>%
#     select(ID, angle, Dist2Aorta, SyndHeight) %>%
#     arrange(ID, angle)
#   return(data)
# }
#
# cleaned_dat <- dat %>% nest(-sheet) %>%
#   mutate(cleaned = map(data, cleanData)) %>%
#   select(sheet, cleaned) %>%
#   unnest() %>% rename(Spine = sheet) %>%
#   mutate(Spine = as.factor(Spine)) %>%
#   mutate(Spine = forcats::lvls_reorder(Spine, c(7:11,4:6,1:3))) %>%
#   mutate(SyndPresent = ifelse(SyndHeight > 0, 1, 0))
#
# saveRDS(cleaned_dat, file = 'data/rda/cleaned_data.rds', compress = T)
# #

dat <- readRDS('data/rda/raw_data_2.rds')
munge_data <- function(d){
  d2 <- enhead(d %>% filter(variable!='angle'), d %>% filter(variable=='angle'), 'W') %>%
    rename(angle = numeric.header) %>%
    select(-ends_with('header')) %>%
    set_names(str_remove(names(.), '.data')) %>%
    select(data_type:angle) %>%
    spatter(variable) %>%
    mutate(rel_angle = seq(-90,90,by=5)) # relative angle from nadir, designated 0
  # d2 %>% mutate(rel_angle = angle - angle[which.min(Dist2Aorta)])
}

munged_data <-  bind_rows(
  lapply(split(dat, dat$sheet), function(d){
    d %>% filter(col > 1) %>%
      behead('N', ID) %>%
      behead('N', Spine) %>%
      behead('N', variable) %>%
      select(row, col, data_type, numeric, variable, ID) %>%
      nest(-ID) %>%
      mutate(munged = map(data, munge_data)) %>%
      select(-data) %>%
      unnest()
  }),
  .id = 'Spine') %>%
  mutate(Spine = as.factor(Spine)) %>%
  mutate(Spine = forcats::lvls_reorder(Spine, c(7:11,4:6,1:3))) %>%
  mutate(SyndPresent = ifelse(SyndHeight > 0, 1, 0))

saveRDS(munged_data, file = 'data/rda/cleaned_data_2.rds', compress = T)

munged_data2 <- munged_data %>% filter(abs(rel_angle) <= 50)
