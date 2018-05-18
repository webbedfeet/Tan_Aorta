ProjTemplate::reload()
datadir <- file.path(ProjTemplate::find_dropbox(), 'NIAMS','Ward','Tan_Aorta')
dir.exists(datadir)
# sheet_names <- openxlsx::getSheetNames(dir(datadir, full.names = T))

dat <- tidyxl::xlsx_cells(normalizePath(file.path(datadir, 'syndAortaData.xlsx')))
saveRDS(dat, 'data/rda/raw_data.rds', compress=T)

dat2 <- tidyxl::xlsx_cells(normalizePath(file.path(datadir, 'syndAortaData180.xlsx')))
saveRDS(dat2, 'data/rda/raw_data_2.rds')
