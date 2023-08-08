
# save data as.rds ####
# folder needs to have a '/' as last character

save_data <- function(data, folder){
  base::saveRDS(data, file = paste0(folder, deparse(substitute(data)), ".rds"))
}

# load rds data ####

load_data <- function(filestring, folder){
  data <- base::readRDS(file = paste0(folder, filestring, ".rds"))
  return(data)
}
