## Function to load RDS files ##
Load_My_Rds <- function(dir, file_pattern = "(?i)\\.rds$", slide_regex = "\\.Slide(\\d+)\\.") {
  
  ## List RDS files within a folder ##
  files <- list.files(dir, pattern = file_pattern, full.names = TRUE)
  
  ## In case no files are RDS ##
  if (!length(files)) stop("No matching files: ", dir)
  
  ## Get the names ##
  base_n <- basename(files) 
  ## Get the slide number based on the provided regex, change as necessary ##
  slide_n <- as.integer(ifelse(grepl(slide_regex, base_n, perl = TRUE),
                               sub(paste0(".*", slide_regex, ".*"), "\\1", base_n, perl = TRUE), NA))
  ## Get their order based on the slide number ##
  f_order <- order(is.na(slide_n), slide_n, base_n)
  
  ## Order them ##
  files <- files[f_order]
  base_n <- base_n[f_order]
  slide_n <- slide_n[f_order]
  
  ## Apply the command to load each file ##
  objs <- lapply(files, readRDS)
  
  ## Get a list with all the RDS objects ##
  setNames(objs, make.unique(ifelse(is.na(slide_n), 
                                    paste0("SeO_", seq_along(objs)), 
                                    paste0("SeO_", slide_n)), 
                             sep = "_dup"))
}


## Function to map coordinates on top of the background image ##
