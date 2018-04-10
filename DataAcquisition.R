ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)
# sheet_names <- openxlsx::getSheetNames(dir(datadir, full.names = T))

dat <- tidyxl::xlsx_cells(normalizePath(dir(datadir, full.names = T), winslash = '/'))
saveRDS(dat, 'data/rda/raw_data.rds', compress=T)
