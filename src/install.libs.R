files <- Sys.glob(paste0("*", SHLIB_EXT))
files <- c(files, dir("htslib-1.21/", pattern="a$", full.names=TRUE))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)

