ProjTemplate::reload()
dat <- readRDS('data/rda/raw_data.rds')

d <- dat %>% filter(sheet == 'T5T6')

cleanData <- function(d) {
  angle <- d %>% filter(row >= 4, col == 1) %>% pull(numeric)
  name <- d %>% filter(row %in% c(1,3), col >= 3) %>%
    mutate(character = ifelse(is.na(character), as.character(numeric),
                              character)) %>%
    select(row, col, character) %>%
    filter(!is.na(character)) %>%
    spread(row, character) %>%
    unite(ID, `1`:`3`, sep = '_') %>%
    pull(ID)
  data <- d %>% filter(row >= 4, col >= 3) %>%
    filter(!is.na(numeric)) %>%
    select(row, col, numeric) %>%
    spread(col, numeric) %>%
    select(-row) %>%
    set_names(name) %>%
    cbind(angle) %>%
    gather(variable, value, -angle) %>%
    separate(variable, c('ID','Measure'), sep = '_') %>%
    spread(Measure, value) %>%
    select(ID, angle, Dist2Aorta, SyndHeight) %>%
    arrange(ID, angle)
  return(data)
}

dat <- tidyxl::xlsx_cells(dir(datadir, full.names = T))

cleaned_dat <- dat %>% nest(-sheet) %>%
  mutate(cleaned = map(data, cleanData)) %>%
  select(sheet, cleaned) %>%
  unnest() %>% rename(Spine = sheet) %>%
  mutate(Spine = as.factor(Spine)) %>%
  mutate(Spine = forcats::lvls_reorder(Spine, c(7:11,4:6,1:3)))

saveRDS(cleaned_dat, file = 'data/rda/cleaned_data.rds', compress = T)
