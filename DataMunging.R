ProjTemplate::reload()
dat <- readRDS('data/rda/raw_data.rds')

d <- dat %>% filter(sheet=='T5T6')

angle <- d %>% filter(row >= 4, col == 1) %>% pull(numeric)
name <- d %>% filter(row < 4, col >= 3) %>%
  mutate(character = ifelse(is.na(character), as.character(numeric), character)) %>%
  select(row, col, character) %>%
  spread(row, character) %>%
  unite(ID, `1`:`3`, sep='_') %>%
  pull(ID)
data <- d %>% filter(row >= 4, col >= 3) %>%
  select(row, col, numeric) %>%
  spread(col, numeric) %>%
  select(-row) %>%
  set_names(name)
