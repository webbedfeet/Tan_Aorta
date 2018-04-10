ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)
# sheet_names <- openxlsx::getSheetNames(dir(datadir, full.names = T))

dat <- tidyxl::xlsx_cells(file.path(dirname(dir(datadir, full.names = T)), basename(dir(datadir))))
saveRDS(dat, 'data/rda/raw_data.rds', compress=T)
