## Function to load paths and report variables ##
Init_project <- function(main_dir) {
  ## Normalize path ##
  main_dir <- normalizePath(main_dir, mustWork = FALSE)
  
  ## Build paths ##
  log_file <- file.path(main_dir, "Results", "Reports", "preprocess_log.txt")
  qc_dir <- file.path(main_dir, "Results", "QC_Images")
  metrics_path <- file.path(main_dir, "Results", "Reports", "preprocess_metrics.csv")
  
  ## Set knitr root if knitr is loaded ##
  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::opts_knit$set(root.dir = main_dir)
  }
  
  ## Load or create metrics ##
  if (file.exists(metrics_path)) {
    metrics <- readr::read_csv(metrics_path, show_col_types = FALSE)
  } else {
    metrics <- tibble::tibble(
      time = as.POSIXct(character()),
      step = character(),
      cells = numeric(),
      genes = numeric(),
      note = character())
  }
  
  ## Put objects in the global environment ##
  assign("main_dir", main_dir, envir = .GlobalEnv)
  assign("log_file", log_file, envir = .GlobalEnv)
  assign("qc_dir", qc_dir, envir = .GlobalEnv)
  assign("metrics", metrics, envir = .GlobalEnv)
  
  invisible(list(main_dir = main_dir,
                 log_file = log_file,
                 qc_dir = qc_dir,
                 metrics = metrics))
}

## Report for the analysis steps ##
log_info <- function(note, fmt, ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- trimws(sprintf(fmt, ...))   # remove leading/trailing spaces
  header <- sprintf("==== %s ====", note)
  date_line <- sprintf("Date: %s", ts)
  
  ## Console ##
  cat("\n", header, "\n", date_line, "\n", msg, "\n", sep = "")
  
  ## File ##
  cat(header, "\n", date_line, "\n", msg, "\n\n",
      file = log_file, append = TRUE, sep = "")
}

## Save info on metrics ##
add_metric <- function(Se_object, step, note = NA_character_) {
  metrics <<- dplyr::bind_rows(metrics,
                               tibble::tibble(time = Sys.time(),
                                              step = step,
                                              cells = ncol(Se_object),
                                              genes = nrow(Se_object),
                                              note = note))
}

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


## Function add percent labels in the quality control figures ##
lbl_cells_percent <- function(b) {
  paste0(format(round(b * ncol(Counts_filtered)), 
                big.mark = ","), " (", scales::percent(b, accuracy = 1), ")")
}
