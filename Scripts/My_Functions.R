## Function to shift and fix dataframes ##
Correct_Dataframes <- function(df) {
  ## Only works on dataframes ##
  stopifnot(is.data.frame(df))
  
  ## Get the total number of columns so we know the last's position
  n <- ncol(df)
  ## Make sure there are columns
  if (n == 0) return(df)
  
  ## Copy to work on ##
  out <- df
  
  ## Shift columns to the left (col i becomes original col i + 1) ##
  if (n >= 2) out[1:(n-1)] <- df[2:n]
  
  ## Make last column NAs ##
  out[[n]] <- NA
  
  ## Split the second to last at the comma ##
  ## Left value stays at n - 1 while right value at n ##
  pen <- out[[n-1]]
  has_comma <- !is.na(pen) & grepl(",", pen, fixed = TRUE)
  
  ## Extract left and right values using the final comma ##
  left_vals  <- ifelse(has_comma, sub(",[^,]*$", "", pen), pen)
  right_vals <- ifelse(has_comma, sub("^.*,", "", pen), NA)
  
  ## Assign to the proper columns and remove whitespaces ##
  out[[n-1]] <- trimws(left_vals)
  out[[n]] <- trimws(right_vals)
  
  ## Type convert the values so numeric stays numeric and chr stays chr ##
  out[] <- lapply(out, function(x) type.convert(x, as.is = TRUE))
  
  ## Return the fixed dataframe ##
  out
}


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
