
# Initiate Project --------------------------------------------------------

## Function to load paths and report variables ##
Init_project <- function(main_dir) {
  ## Normalize path ##
  main_dir <- normalizePath(main_dir, mustWork = FALSE)
  
  ## Set working directory ##
  setwd(dir = main_dir)
  
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


# Loading - Fixing - Generic ----------------------------------------------

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

## Make png backgrounds transparent ##
Transparent <- function(plot) {
  plot <- plot + 
    theme(panel.background  = element_rect(fill = "transparent", colour = NA),
          plot.background   = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
}


# Quality Control ---------------------------------------------------------

## Assess batch effect function ##
AssessClusterBatch <- function(obj,
                               cluster_col,
                               group_col,
                               reduction = "pca",
                               dims_use = 1:10,
                               ## composition plot options ##
                               pct_labels = TRUE,
                               pct_label_min = 0.05, # only label segments >= 5% by default
                               pct_label_digits = 0, # 0 = integer percent
                               show_total_n = TRUE,
                               ## R2 options ##
                               r2_within_cluster = TRUE,
                               min_cells_per_cluster = 200,
                               sample_size = NULL,
                               seed = 1511, 
                               save = FALSE,
                               out_dir = NULL,
                               file_prefix = NULL,
                               save_device = "png",  
                               dpi = 300) {
  
  stopifnot(cluster_col %in% colnames(obj@meta.data))
  stopifnot(group_col %in% colnames(obj@meta.data))
  stopifnot(reduction %in% Reductions(obj))
  
  ## Pull metadata ##
  meta_df <- FetchData(obj, vars = c(cluster_col, group_col)) 
  colnames(meta_df) <- c("Cluster", "Group")
  meta_df$Cluster <- factor(meta_df$Cluster)
  meta_df$Group <- factor(meta_df$Group)
  
  ## Composition table ##
  tab <- table(meta_df$Cluster, meta_df$Group)
  prop_tab <- prop.table(tab, margin = 1)
  
  df_tab <- as.data.frame(tab)
  colnames(df_tab) <- c("Cluster", "Group", "Freq")
  
  ## Totals plus proportions (within cluster) ##
  totals <- df_tab %>%
    group_by(Cluster) %>%
    summarise(Total = sum(Freq), .groups = "drop")
  
  df_tab <- df_tab %>%
    left_join(totals, by = "Cluster") %>%
    mutate(Prop = ifelse(Total > 0, Freq / Total, NA_real_))
  
  ## Percentage labels for each segment ##
  df_tab <- df_tab %>%
    mutate(PctLabel = ifelse(pct_labels & !is.na(Prop) & Prop >= pct_label_min,
                             paste0(round(Prop * 100, pct_label_digits), "%"), ""))
  
  ## Composition bar plots ##
  p_comp <- ggplot(df_tab, aes(x = Cluster, y = Prop, fill = Group)) +
    geom_col(width = 0.85) +
    labs(y = "Proportion within cluster", x = "Cluster") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (pct_labels) {
    p_comp <- p_comp +
      geom_text(aes(label = PctLabel),
                position = position_stack(vjust = 0.5), size = 3)
  }
  
  if (show_total_n) {
    p_comp <- p_comp +
      geom_text(data = totals, aes(x = Cluster, y = 1.02, label = paste0("n=", Total)),
                inherit.aes = FALSE, size = 3) +
      coord_cartesian(ylim = c(0, 1.08))
  }
  
  ## Embeddings ##
  emb <- Embeddings(obj, reduction = reduction)  
  dims_use <- intersect(dims_use, seq_len(ncol(emb)))
  emb <- emb[, dims_use, drop = FALSE]
  
  if (!is.null(sample_size) && nrow(emb) > sample_size) {
    set.seed(seed)
    keep <- sample(seq_len(nrow(emb)), sample_size)
    emb_use <- emb[keep, , drop = FALSE]
    meta_use <- meta_df[keep, , drop = FALSE]
  } else {
    emb_use <- emb
    meta_use <- meta_df
  }
  
  ## Helper: adjusted R2 of x ~ batch ##
  get_r2 <- function(x, batch) {
    batch <- droplevels(factor(batch))
    if (nlevels(batch) < 2) return(NA_real_)
    summary(lm(x ~ batch))$adj.r.squared
  }
  
  ## Global R2 per dimension ##
  r2_global <- apply(emb_use, 2, get_r2, batch = meta_use$Group)
  
  r2_global_df <- tibble(Dim = factor(names(r2_global), levels = names(r2_global)),
                         R2 = as.numeric(r2_global)
  )
  
  p_r2_global <- ggplot(r2_global_df, aes(x = Dim, y = R2)) +
    geom_col() +
    labs(title = paste0("Adjusted R² of ", group_col, " per ", reduction, " dimension"),
         x = "Dimension", 
         y = "Adjusted R²") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## Within cluster R2 ##
  r2_cluster_df <- NULL
  p_r2_cluster  <- NULL
  
  if (isTRUE(r2_within_cluster)) {
    cl_levels <- levels(meta_use$Cluster)
    
    r2_cluster_df <- lapply(cl_levels, function(cl) {
      idx <- which(meta_use$Cluster == cl)
      if (length(idx) < min_cells_per_cluster) return(NULL)
      
      r2 <- apply(emb_use[idx, , drop = FALSE], 2, get_r2, batch = meta_use$Group[idx])
      tibble(Cluster = cl, 
             Dim = factor(names(r2), levels = colnames(emb_use)), 
             R2 = as.numeric(r2), 
             n = length(idx))
    }) %>% bind_rows()
    
    ## Heatmap: clusters x dimensions ##
    p_r2_cluster <- ggplot(r2_cluster_df, aes(x = Dim, y = Cluster, fill = R2)) +
      geom_tile() +
      labs(title = paste0("Within-cluster adjusted R² of ", group_col),
           subtitle = paste0("Clusters with n < ", min_cells_per_cluster, " omitted"),
           x = "Dimension", 
           y = "Cluster") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  ## Save the files if save = TRUE ##
  saved_files <- character(0)
  
  if (isTRUE(save)) {
    if (is.null(out_dir) || !nzchar(out_dir)) {
      stop("When save=TRUE you must provide out_dir (a folder path).")
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (is.null(file_prefix) || !nzchar(file_prefix)) {
      file_prefix <- paste0(cluster_col, "__by__", group_col, "__", reduction)
    }
    
    # helper to save only if plot exists
    save_plot <- function(p, suffix, w, h) {
      if (is.null(p)) return(invisible(NULL))
      f <- file.path(out_dir, paste0(file_prefix, "_", suffix, ".png"))
      ggplot2::ggsave(filename = f, 
                      plot = p,
                      width = w, height = h, dpi = dpi,
                      device = save_device)
      saved_files <<- c(saved_files, f)
      invisible(NULL)
    }
    
    save_plot(p_comp, "composition", w = 12, h = 6)
    save_plot(p_r2_global, "r2_global", w = 8,  h = 4)
    save_plot(p_r2_cluster,"r2_within_cluster_map", w = 10, h = 6)
  }
  
  ## Final list ##
  list(tab = tab,
       prop_tab = prop_tab,
       composition_df = df_tab,
       composition_plot = p_comp,
       r2_global_df = r2_global_df,
       r2_global_plot = p_r2_global,
       r2_within_cluster_df = r2_cluster_df,
       r2_within_cluster_plot = p_r2_cluster, 
       saved_files = saved_files)
}


# Annotation visualisations -----------------------------------------------

## Function to visualise cell annotation results ##
Custom_flighpath <- function(annotation_obj,
                             title_text,
                             plot_seed = 1511,
                             point_size = 0.2,
                             point_alpha = 0.8,
                             label_size = 3) {
  
  stopifnot(!is.null(annotation_obj$logliks),
            !is.null(annotation_obj$profiles),
            !is.null(annotation_obj$clust))
  
  Flight <- InSituType::flightpath_layout(logliks = annotation_obj$logliks,
                                          profiles = annotation_obj$profiles)
  
  cl_levels <- rownames(as.data.frame(Flight$clustpos))
  if (is.null(cl_levels) || !length(cl_levels)) {
    cl_levels <- sort(unique(as.character(annotation_obj$clust)))
  }
  
  mean_conf <- Flight$meanconfidence
  mean_conf_map <- if (!is.null(names(mean_conf)) && all(names(mean_conf) != "")) {
    tibble::tibble(Cluster = names(mean_conf), mean_conf = as.numeric(mean_conf))
  } else {
    tibble::tibble(Cluster = cl_levels,
                   mean_conf = as.numeric(mean_conf)[seq_len(min(length(mean_conf), length(cl_levels)))])
  }
  
  cl_tab <- sort(table(annotation_obj$clust), decreasing = TRUE)
  
  Cluster_stats <- tibble::tibble(Cluster = factor(names(cl_tab), levels = cl_levels),
                                  n_cells = as.integer(cl_tab)) %>%
    dplyr::left_join(mean_conf_map %>% dplyr::mutate(Cluster = factor(Cluster, levels = cl_levels)),
                     by = "Cluster") %>%
    dplyr::mutate(Cluster_lab = sprintf("%s (n=%s, mean conf=%.2f)",
                                        as.character(Cluster),
                                        scales::comma(n_cells),
                                        dplyr::if_else(is.na(mean_conf), 0, mean_conf))) %>%
    dplyr::arrange(dplyr::desc(n_cells))
  
  ## Plot the cells ##
  Celltypes_df <- as.data.frame(Flight$cellpos) %>%
    dplyr::mutate(Cluster = factor(annotation_obj$clust, levels = levels(Cluster_stats$Cluster))) %>%
    dplyr::left_join(Cluster_stats, by = "Cluster") %>%
    dplyr::mutate(Cluster_lab = factor(Cluster_lab, levels = Cluster_stats$Cluster_lab))
  
  ## Cluster label positions ##
  Celltypes_pos_df <- as.data.frame(Flight$clustpos) %>%
    tibble::rownames_to_column("Cluster") %>%
    dplyr::mutate(Cluster = factor(Cluster, levels = levels(Cluster_stats$Cluster))) %>%
    dplyr::left_join(Cluster_stats, by = "Cluster")
  
  set.seed(plot_seed)
  
  ## Scatter plot ##
  ggplot2::ggplot(Celltypes_df, ggplot2::aes(x = x, y = y, colour = Cluster_lab)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggrepel::geom_label_repel(data = Celltypes_pos_df,
                              ggplot2::aes(x = x, y = y, label = Cluster),
                              colour = "black", fill = "white",
                              label.size = 0.2, fontface = "bold",
                              label.r = grid::unit(0.1, "lines"),
                              size = label_size,
                              show.legend = FALSE,
                              max.overlaps = Inf) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title_text, colour = "InSituType Clusters") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 1)) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.title = ggplot2::element_text(face = "bold", size = 11),
                   legend.text = ggplot2::element_text(face = "bold", size = 9))
}

## Helper for the pheatmap visualisations ##
## Scales matrix expression rows and drops genes who dont meet the threshold ##
Scale_rows_pheatmap <- function(mat, min_max = 0.2, min_keep = 0.2) {
  row_max <- matrixStats::rowMaxs(mat)
  scaled <- mat / pmax(row_max, min_max)
  scaled[row_max > min_keep, , drop = FALSE]
}

# InSituType functions ----------------------------------------------------

## Annotation and refinement function ##
## 1. Loads the reference celltype profiles ##
## 2. Filters the genes of the seurat counts based on the reference ##
## 3. Creates the initial annotation ##
## 4. Visualises the initial results ##
## 5. Refines the annotation based on the size of the "unknown" clusters ##
## 6. Re-visualises after refinement ##

Run_insitutype_panel <- function(Seurat_obj,
                                 Ref_matrix = NULL,
                                 ref_tag,
                                 n_clust = 0,
                                 refine = TRUE,
                                 auto_delete_letters = TRUE,
                                 to_delete = NULL,
                                 seed = 1511) {
  
  ## Define the folders ##
  ann_dir <- file.path(main_dir, "Data", "Annotations")
  res_dir <- file.path(main_dir, "Results", "Annotations")
  
  if (!dir.exists(ann_dir)) dir.create(ann_dir, recursive = TRUE)
  if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  
  ## Folder for each panel ##
  ref_tag_dir <- file.path(res_dir, ref_tag)
  if (!dir.exists(ref_tag_dir)) dir.create(ref_tag_dir, recursive = TRUE)
  
  ## Define the paths ##
  ann_rds <- file.path(ann_dir, paste0("Annotation_", ref_tag, ".Rds"))
  ann_ref_rds <- file.path(ann_dir, paste0("Annotation_", ref_tag, "_refined.Rds"))
  flight_png <- file.path(ref_tag_dir, paste0("Flightpath_", ref_tag, ".png"))
  flight_ref_png <- file.path(ref_tag_dir, paste0("Flightpath_", ref_tag, "_refined.png"))
  pheatmap_png <- file.path(ref_tag_dir, paste0("Pheatmap_", ref_tag, ".png"))
  pheatmap_ref_png <- file.path(ref_tag_dir, paste0("Pheatmap_", ref_tag, "_refined.png"))
  
  
  ## Load counts and the reference ##
  Counts <- GetAssayData(Seurat_obj, assay = "RNA", layer = "counts") %>% 
    t()
  
  Neg_probs <- GetAssayData(Seurat_obj, assay = "negprobes", layer = "counts") %>% 
    t()
  
  
  ## Load reference profiles or keep unsupervised ##
  unsupervised <- FALSE
  
  if (is.null(Ref_matrix)) {
    Ref_mat <- NULL
    unsupervised <- TRUE
    
  } else if (is.character(Ref_matrix) && length(Ref_matrix) == 1 && 
             grepl("\\.csv$", Ref_matrix, ignore.case = TRUE)) {
    
    unsupervised <- FALSE
    Ref_mat <- read.csv(Ref_matrix, row.names = 1, header = TRUE, check.names = FALSE)
    
  } else if (inherits(Ref_matrix, c("matrix", "array", "data.frame", "dgCMatrix"))) {
    Ref_mat <- as.matrix(Ref_matrix)
    
  } else {
    stop("Ref_matrix must be either:\n",
         "  - NULL (unsupervised)\n",
         "  - path to a CSV file\n",
         "  - matrix / data.frame / dgCMatrix")
  }
  
  IF_Markers <- Seurat_obj@meta.data %>% 
    select(contains(c("Mean", "Max")))
  
  Neg_Mean_per_total_Count <- mean(rowMeans(Neg_probs)) / mean(rowSums(Counts))
  Per_Cell_BG <- rowSums(Counts) * Neg_Mean_per_total_Count
  
  edited_data_dir <- file.path(main_dir, "Data", "Edited_data")
  
  if (file.exists(file.path(edited_data_dir, "ImmunoFlourescence_Cohorts.Rds"))) {
    IS_Cohort <- readRDS(file.path(edited_data_dir, "ImmunoFlourescence_Cohorts.Rds"))
  } else {
    IF_Markers <- Seurat_obj@meta.data %>% 
      select(contains(c("Mean", "Max")))
    
    IS_Cohort <- fastCohorting(IF_Markers, gaussian_transform = TRUE)
    IS_Cohort <- interaction(IS_Cohort, Seurat_obj$Slide, drop = TRUE)
  }
  
  ## Run InSituType (or load if existing) ##
  if (file.exists(ann_rds)) {
    Annotation <- readRDS(ann_rds)
    
  } else {
    set.seed(seed)
    
    if (is.null(Ref_matrix)) {
      
      ## Unsupervisedvised InSituType ##
      if (n_clust <= 0) {
        stop("In unsupervised mode (Ref_matrix = NULL), n_clust must be > 0.") 
      }
      
      Annotation <- InSituType::insitutype(x = Counts,
                                           neg = Matrix::rowMeans(Neg_probs),
                                           reference_profiles = NULL,
                                           n_clust = n_clust,   
                                           bg = Per_Cell_BG,
                                           cohort = IS_Cohort)
      
    } else {
      
      ## Reference based InSituType ##
      Annotation <- InSituType::insitutype(x = Counts,
                                           neg = Matrix::rowMeans(Neg_probs),
                                           reference_profiles = Ref_mat,
                                           update_reference_profiles = FALSE,
                                           n_clust = n_clust, 
                                           bg = Per_Cell_BG, 
                                           cohort = IS_Cohort)
    } 
    saveRDS(Annotation, ann_rds) 
  }
  
  ## Initial pheatmap ##
  Pheat <- pheatmap::pheatmap(Scale_rows_pheatmap(Annotation$profiles),
                              fontsize_row = 4,
                              col = colorRampPalette(c("white", "darkblue"))(100), 
                              filename = pheatmap_png)
  
  ## Initial flightpath ##
  Flight <- InSituType::flightpath_layout(logliks = Annotation$logliks,
                                          profiles = Annotation$profiles)
  
  Flight_plot <- Custom_flighpath(flightpath_obj = Flight, 
                                  annotation_obj = Annotation, 
                                  unsupervised = unsupervised,
                                  title = ref_tag)
  
  ggplot2::ggsave(plot = Flight_plot, filename = flight_png,
                  width = 16, height = 10, dpi = 600, device = "png")
  
  ## Decide which clusters to refine/delete ##
  # to_delete = TRUE the "unknown" clusters are removed and the cells are supplied to the immediate next cluster ##
  if (is.null(to_delete) && auto_delete_letters) {
    letter_clusters <- intersect(letters, unique(Annotation$clust))
    to_delete <- sort(letter_clusters)
  }
  
  ## Optional refineClusters step ##
  Annotation_ref <- NULL
  Flight_ref <- NULL
  Flight_ref_plot <- NULL
  Pheat_ref <- NULL
  
  if (refine && length(to_delete) > 0) {
    
    if (file.exists(ann_ref_rds)) {
      Annotation_ref <- readRDS(ann_ref_rds)
      
    } else {
      
      Annotation_ref <- InSituType::refineClusters(logliks = Annotation$logliks,
                                                   counts = Counts,
                                                   neg = Matrix::rowMeans(Neg_probs),
                                                   to_delete = to_delete)
      saveRDS(Annotation_ref, ann_ref_rds)
    }
    
    Pheat_ref <- pheatmap::pheatmap(Scale_rows_pheatmap(Annotation_ref$profiles),
                                    fontsize_row = 4,
                                    col = colorRampPalette(c("white", "darkblue"))(100), 
                                    filename = pheatmap_ref_png)
    
    Flight_ref <- InSituType::flightpath_layout(logliks = Annotation_ref$logliks,
                                                profiles = Annotation_ref$profiles)
    
    Flight_ref_plot <- Custom_flighpath(flightpath_obj = Flight_ref, 
                                        annotation_obj = Annotation_ref, 
                                        unsupervised = unsupervised, 
                                        title = ref_tag)
    
    ggplot2::ggsave(plot = Flight_ref_plot, filename = flight_ref_png, 
                    width = 16, height = 10, dpi = 600, device = "png")
    
  }
  
  ## Save the annotations to the common csv file ##
  Clusters_total_path <- file.path(ref_tag_dir, paste0(ref_tag, "_clusters.csv"))
  
  ## Get original annotations ##
  Clusters_current <- tibble::enframe(Annotation$clust, 
                                      name = "cell_ID", 
                                      value = paste0(ref_tag, "_clusters"))
  
  ## Check if there is a refined annotation ##
  Clusters_current_ref <- NULL
  if (file.exists(ann_ref_rds)) {
    Annotation_ref <- readRDS(ann_ref_rds)
    
    Clusters_current_ref <- tibble::enframe(Annotation_ref$clust,
                                            name  = "cell_ID",
                                            value = paste0(ref_tag, "_clusters_ref")) 
  }
  
  ## If there is a preexisting CLuster file load it, otherwise create it ##
  if (file.exists(Clusters_total_path)) {
    Cell_annotations_total <- read_csv(Clusters_total_path)
  } else {
    Cell_annotations_total <- Clusters_current
  }
  
  Cell_annotations_total <- Cell_annotations_total %>% 
    left_join(Clusters_current, by = "cell_ID")
  
  ## If refinement was performed add those clusters too ##
  if (!is.null(Clusters_current_ref)) {
    Cell_annotations_total <- Cell_annotations_total %>% 
      left_join(Clusters_current_ref, by = "cell_ID")
    }
  
  ## Save the updated file ##
  write_csv(x = Cell_annotations_total, file = file.path(ref_tag_dir, paste0(ref_tag, "_clusters.csv")))
  
}

## Function to specify the paths within the main function ##
Build_InSituType_Paths <- function(main_dir,
                                   Ref_family,     
                                   Ref_tag,        
                                   Ref_matrix = NULL,
                                   n_clust = 0,
                                   cohort_by = c("Slide", "Batch", "Group"), 
                                   Auto = FALSE) {
  
  cohort_by <- match.arg(cohort_by)
  base_ann_dir <- file.path(main_dir, "Results/Annotations")
  
  ## In case we want it completely unsupervised ##
  auto_selected <- isTRUE(attr(n_clust, "auto_selected"))
  
  if (is.null(Ref_matrix)) {
    mode <- "unsupervised"
    main_ann_dir <- file.path(base_ann_dir, "Unsupervised")
    
    ## Make the name prettier ##
    Ref_tag2 <- tools::toTitleCase(Ref_tag)
    
    if (auto_selected) {
      ## Completely unsupervised ##
      main_tag <- paste0(Ref_tag2, "_auto_n_clust_", as.integer(n_clust))
    } else {
      ## Unsupervised but with specified clusters ##
      main_tag <- paste0(Ref_tag2, "_n_clust_", as.integer(n_clust))
    }
    
  } else {
    main_ann_dir <- file.path(base_ann_dir, Ref_family)
    
    if (n_clust <= 0) {
      mode <- "supervised"
      main_tag <- paste0(Ref_tag, "_", mode)
    } else {
      mode <- "semi_supervised"
      main_tag <- paste0(Ref_tag, "_", mode, "_n_clust_", n_clust)
    }
  }
  
  ann_dir <- file.path(main_ann_dir, main_tag)
  ann_rds <- file.path(ann_dir, paste0(main_tag, "_annotation.Rds"))
  
  list(base_ann_dir = base_ann_dir,
       main_ann_dir = main_ann_dir,
       ann_dir = ann_dir,
       ann_rds = ann_rds,
       main_tag = main_tag,
       mode = mode,
       cohort_by = cohort_by,
       n_clust = n_clust,
       has_reference = !is.null(Ref_matrix))
}

## Function that chooses to load pre existing Rds annotation or create new ##
Cache_Load_or_Run <- function(rds_path,
                              run_fun,
                              ...,
                              overwrite = FALSE,
                              verbose = TRUE,
                              save_params = TRUE,
                              params = NULL) {
  
  ## Check if the output directory already exists ##
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  
  ## Load the annotation if it exists, unless overwrite is activated ##
  if (!overwrite && file.exists(rds_path)) {
    if (verbose) message("Loading existing annotation: ", rds_path)
    return(readRDS(rds_path))
  }
  
  if (verbose) {
    if (overwrite && file.exists(rds_path)) {
      message("Overwriting existing annotation: ", rds_path)
    } else {
      message("No cached annotation found. Running and saving: ", rds_path)
    }
  }
  
  ##Run the main function later on ##
  obj <- run_fun(...)
  
  ## Save annotation ##
  saveRDS(obj, rds_path)
  
  ## Save run parameters (optional) ##
  if (isTRUE(save_params)) {
    if (is.null(params)) {
      params <- list(...)
    }
    saveRDS(params, file.path(dirname(rds_path), "params.rds"))
  }
  
  obj
}

## Main function to run first InSituType run ##
Run_InSituType_Main <- function(Seurat_obj,
                                main_dir,
                                Ref_family,
                                Ref_tag,
                                Ref_matrix = NULL,
                                n_clust = 0,
                                cohort_by = c("Slide", "Batch", "Group"),
                                overwrite = FALSE,
                                seed = 1511, 
                                visualisation = TRUE,
                                save_plots = TRUE,
                                plot_seed = 1511,
                                heat_min_max = 0.2,
                                heat_min_keep = 0.2) {
  
  ## Build output paths ##
  paths <- Build_InSituType_Paths(main_dir = main_dir,
                                  Ref_family = Ref_family,
                                  Ref_tag = Ref_tag,
                                  Ref_matrix = Ref_matrix,
                                  n_clust = n_clust,
                                  cohort_by = cohort_by)
  
  ## Store TRUE if n_clust is a range ##
  auto_run <- is.null(Ref_matrix) && length(n_clust) > 1L
  
  # This avoids re-running InSituType if you already did an auto-selected run
  if (auto_run && !overwrite) {
    main_ann_dir_auto <- file.path(main_dir, "Results", "Annotations", "Unsupervised")
    if (dir.exists(main_ann_dir_auto)) {
      ## NOTE: pattern must match your naming convention in Build_InSituType_Paths
      Ref_tag2 <- tools::toTitleCase(Ref_tag)
      pat <- paste0("^", Ref_tag2, "_auto_n_clust_\\d+$")
      hit <- list.dirs(main_ann_dir_auto, full.names = TRUE, recursive = FALSE)
      hit <- hit[grepl(pat, basename(hit))]
      
      if (length(hit)) {
        hit <- hit[order(file.info(hit)$mtime, decreasing = TRUE)]
        rds <- list.files(hit[1], pattern = "_annotation\\.Rds$", full.names = TRUE)
        if (length(rds) && file.exists(rds[1])) {
          message("Loading existing auto-selected annotation: ", rds[1])
          return(list(annotation = readRDS(rds[1]),
                      paths = list(main_ann_dir = main_ann_dir_auto,
                                   ann_dir = hit[1],
                                   ann_rds = rds[1])))
        }
      }
    }
  }
  
  ## The actual function ##
  run_insitutype <- function() {
    set.seed(seed)
    
    ## Counts: cells x genes, in sparse format ##
    Counts <- GetAssayData(Seurat_obj, assay = "RNA", layer = "counts") %>% 
      t()
    Neg_probs <- GetAssayData(Seurat_obj, assay = "negprobes", layer = "counts") %>% 
      t()
    
    ## Immunoflourescence markers ##
    IF_Markers <- Seurat_obj@meta.data %>%
      dplyr::select(dplyr::contains(c("Mean", "Max"))) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      dplyr::select(dplyr::where(~ !any(is.na(.))))
    
    ## Cohort ##
    IS_Cohort <- fastCohorting(IF_Markers, gaussian_transform = TRUE)
    IS_Cohort <- interaction(IS_Cohort, Seurat_obj[[cohort_by]][, 1], drop = TRUE)
    
    ## Per cell background noise ##
    neg_means <- Matrix::rowMeans(Neg_probs)
    totals <- Matrix::rowSums(Counts)
    ratio <- neg_means / totals
    ratio <- ratio[is.finite(ratio) & ratio > 0]
    Neg_Mean_per_total_Count <- median(ratio)
    Per_Cell_BG <- totals * Neg_Mean_per_total_Count
    
    ## Unsupervised version, no reference profiles ##
    if (is.null(Ref_matrix)) {
      if (!auto_run && (length(n_clust) != 1L || n_clust <= 0)) {
        stop("Unsupervised mode requires n_clust > 0 (or a range like 5:20).")
      }
      
      return(InSituType::insitutype(x = Counts,
                                    neg = neg_means,
                                    reference_profiles = NULL,
                                    n_clust = n_clust,
                                    bg = Per_Cell_BG,
                                    cohort = IS_Cohort))
    }
    
    ## Semi-/Supervised version, reference profiles but with specified cluster number ##
    Ref_mat <- as.matrix(Ref_matrix)
    
    InSituType::insitutype(x = Counts,
                           neg = neg_means,
                           reference_profiles = Ref_matrix,
                           n_clust = n_clust,            # 0 = supervised; > 0 = semi-supervised
                           bg = Per_Cell_BG,
                           cohort = IS_Cohort)
  }
  
  ## Run first in case we provide a range of n_cluster values ##
  Annotation_tmp <- run_insitutype()
  
  ## Save the optimal number provided by InSituType ##
  k_final <- length(unique(Annotation_tmp$clust))
  
  if (auto_run) {
    k_for_naming <- as.integer(k_final)
    ## Informs Build_InSituType_Paths to use autoK naming ##
    attr(k_for_naming, "auto_selected") <- TRUE  
    
    paths <- Build_InSituType_Paths(main_dir = main_dir,
                                    Ref_family = Ref_family,
                                    Ref_tag = Ref_tag,
                                    Ref_matrix = NULL,
                                    n_clust = k_for_naming,
                                    cohort_by = cohort_by)
  } else {
    paths <- Build_InSituType_Paths(main_dir = main_dir,
                                    Ref_family = Ref_family,
                                    Ref_tag = Ref_tag,
                                    Ref_matrix = Ref_matrix,
                                    n_clust = n_clust,
                                    cohort_by = cohort_by)
  }
  
  ## Ensure directories exist (no overwriting) ##
  dir.create(paths$ann_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Cache/save using the FINAL path ##
  # We pass run_fun that returns the already computed object to avoid re-running
  Annotation <- Cache_Load_or_Run(rds_path = paths$ann_rds,
                                  run_fun = function() Annotation_tmp,
                                  overwrite = overwrite,
                                  params = list(Ref_family = Ref_family,
                                                Ref_tag = Ref_tag,
                                                cohort_by = cohort_by,
                                                n_clust_input = n_clust,
                                                auto_selected = auto_run,
                                                k_final = k_final,
                                                seed = seed))
  
  ## Flightpath Visualisation ##
  flight_plot <- NULL
  
  if (isTRUE(visualisation)) {
    
    title_text <- if (is.null(Ref_matrix)) {
      paste0("InSituType annotation: ", paths$main_tag)
    } else {
      paste0("InSituType annotation using reference: ", paths$main_tag)
    }
    
    flight_plot <- Custom_flighpath(annotation_obj = Annotation,
                                              title = title_text,
                                              plot_seed = plot_seed)
    
    if (isTRUE(save_plots)) {
      ggplot2::ggsave(filename = file.path(paths$ann_dir, paste0(paths$main_tag, "_flightpath.png")),
                      plot = flight_plot, width = 16, height = 10, dpi = 600, device = "png")
    }
    
    pheat_file <- if (isTRUE(save_plots)) {
      file.path(paths$ann_dir, paste0(paths$main_tag, "_profiles_pheatmap.png"))
    } else {
      NULL
    }
    
    pheatmap::pheatmap(Scale_rows_pheatmap(Annotation$profiles, 
                                           min_max = heat_min_max, 
                                           min_keep = heat_min_keep),
                       main = paste0("Pheatmap: ", paths$main_tag), 
                       fontsize_row = 4, 
                       col = grDevices::colorRampPalette(c("white", "darkblue"))(100), 
                       filename = pheat_file)
  }
  
  list(annotation = Annotation, paths = paths)
  
}


# InSituType helper functions ---------------------------------------------

## Function which loads annotations from the output folder ##
Load_Annotation_Rds <- function(main_dir, 
                                query = NULL, 
                                name_from = c("filename", "folder"), 
                                keep_paths_vector = TRUE) {
  
  name_from <- match.arg(name_from)
  base_dir <- file.path(main_dir, "Results", "Annotations")
  stopifnot(dir.exists(base_dir))
  
  rds_files <- list.files(base_dir,
                          pattern = "_annotation\\.Rds$",
                          recursive = TRUE,
                          full.names = TRUE)
  
  if (!length(rds_files)) {
    stop("No *_annotation.Rds files found under: ", base_dir)
  }
  
  if (!is.null(query) && nzchar(query)) {
    rds_files <- rds_files[grepl(query, rds_files, fixed = TRUE)]
    if (!length(rds_files)) {
      stop("No annotation matches query = '", query, "'.")
    }
  }
  
  nm <- if (name_from == "folder") {
    basename(dirname(rds_files))                
  } else {
    sub("_annotation\\.Rds$", "", basename(rds_files))
  }
  nm <- make.unique(gsub("[^A-Za-z0-9._-]+", "_", nm))
  
  ## Load ##
  annotations <- Map(
    f = function(fp, nm_i) {
      obj <- readRDS(fp)

      ## Add object path to the annotation list ##      
      if (is.list(obj)) {
        obj$path <- fp
        obj
      } else {
        ## Add it in a previous step to avoid errors if its not a list ##
        list(value = obj, path = fp, name = nm_i)
      }
    },
    fp = rds_files,
    nm_i = nm)
  
  annotations <- setNames(annotations, nm)
  paths_named <- setNames(rds_files, nm)
  
  out <- list(annotations = annotations, base_dir = base_dir)
  if (isTRUE(keep_paths_vector)) out$paths <- paths_named
  out
  
}

# Reduction visualisations ------------------------------------------------

## Better UMAP visualisation ##
PlotUMAPWithCentroids <- function(Seurat_obj,
                                  reduction = "umap",         
                                  group.by,                   
                                  title = NULL,
                                  legend_ncol = 1,
                                  point_size = 0.1,
                                  point_alpha = 0.5,
                                  label_size = 3, 
                                  save = FALSE, 
                                  width = 16, 
                                  height = 10, 
                                  dpi = 600, 
                                  out_dir = NULL, 
                                  filename = NULL) {
  
  ## Retrieve the reduction's information from the object ##
  emb <- Seurat::Embeddings(Seurat_obj, reduction)
  
  emb_df <- as.data.frame(emb) %>%
    rownames_to_column("cell_ID") %>% 
    rename(UMAP_1 = !!colnames(emb)[1],
           UMAP_2 = !!colnames(emb)[2])
  
  ## Get grouping variable from metadata ##
  meta_df <- Seurat_obj@meta.data %>%
    select(cell_ID, !!group.by)
  
  colnames(meta_df)[2] <- "group_var"  # standard name
  
  plot_df <- emb_df %>%
    left_join(meta_df, by = "cell_ID")
  
  group_sym <- sym("group_var")
  
  ## Get the centroids for the labels ##
  centroids <- plot_df %>%
    group_by(!!group_sym) %>%
    summarise(UMAP_1 = median(UMAP_1),
              UMAP_2 = median(UMAP_2),
              .groups = "drop")
  
  ## Assign the title ##
  if (is.null(title)) {
    tag <- gsub(pattern = "clusters", replacement = "reference", group.by)
    tag <- gsub(pattern = "_", replacement = " ", tag)
    title <- paste0("UMAP projection by ", tag, " annotations")
  } else {
    title <- paste0("UMAP projection by ", title, " annotations")
  }
  
  ## Perform the plotting
  UMAP_plot <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, colour = !!group_sym)) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_label_repel(data = centroids,
                     aes(label = !!group_sym),
                     colour = "black",
                     fill = "white",
                     size = label_size,
                     label.size  = 0.25,
                     label.r = unit(0.15, "lines"),
                     segment.size = 0.2,
                     show.legend = FALSE, 
                     max.overlaps = Inf) +
    # coord_equal() +
    labs(title  = title,
         x = "UMAP 1",
         y = "UMAP 2",
         colour = "Annotation clusters") +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    theme_minimal(base_size = 12, base_family = "sans") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(size = 9))
  
  if (save) {
    
    if (is.null(out_dir)) {
      out_dir <- file.path(main_dir, "Results/UMAP_Plots")
    }
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (is.null(filename)) {
      safe_group <- gsub("[^A-Za-z0-9_]+", "_", group.by)
      safe_red   <- gsub("[^A-Za-z0-9_]+", "_", reduction)
      filename <- paste0("UMAP_", safe_red, "_", safe_group, ".png")
    }
    
    ggsave(plot = UMAP_plot, filename = file.path(out_dir, filename), 
           width = width, height = height, dpi = dpi, device = "png")
  }
  
  return(UMAP_plot)
  
}

PlotInSituTypeUMAP <- function(seurat_obj,
                               insitutype_run,
                               reduction = "umap",
                               meta_col = NULL,
                               label = TRUE,
                               repel = TRUE,
                               pt.size = 1,
                               raster = TRUE,
                               save = FALSE,
                               filename = NULL,
                               width = 8,
                               height = 6,
                               dpi = 300) {
  
  if (is.list(insitutype_run) && !is.null(insitutype_run$annotation)) {
    ann <- insitutype_run$annotation
    paths <- insitutype_run$paths
  } else {
    ann <- insitutype_run
    paths <- NULL
  }
  
  stopifnot(reduction %in% Seurat::Reductions(seurat_obj))
  stopifnot(!is.null(ann$clust))
  
  if (is.null(meta_col) || !nzchar(meta_col)) {
    if (!is.null(paths$main_tag)) {
      meta_col <- paste0(reduction, ": ", paths$main_tag)
    } else {
      meta_col <- "InSituType_clusters"
    }
  }
  
  cl <- ann$clust
  if (is.null(names(cl))) {
    stop("ann$clust has no names. It must be named by cell IDs to align to Seurat_obj.")
  }
  
  shared <- intersect(colnames(seurat_obj), names(cl))
  if (!length(shared)) stop("No overlapping cell IDs between Seurat object and annotation.")
  
  seurat_obj[[meta_col]] <- cl[colnames(seurat_obj)]
  
  # --- NEW: title/subtitle logic ---
  ann_name <- if (!is.null(paths$main_tag) && nzchar(paths$main_tag)) {
    paths$main_tag
  } else {
    meta_col
  }
  
  p <- Seurat::DimPlot(
    seurat_obj,
    reduction = reduction,
    group.by  = meta_col,
    label     = label,
    repel     = repel,
    pt.size   = pt.size,
    raster    = raster
  ) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", title = NULL, subtitle = NULL) +
    patchwork::plot_annotation(
      title = reduction,
      subtitle = ann_name
    )
  
  if (isTRUE(save)) {
    out_dir <- if (!is.null(paths$ann_dir)) paths$ann_dir else getwd()
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (is.null(filename) || !nzchar(filename)) {
      filename <- paste0(meta_col, "_", reduction, ".png")
      filename <- gsub("[^A-Za-z0-9._-]+", "_", filename)
    }
    
    ggplot2::ggsave(
      filename = file.path(out_dir, filename),
      plot = p, width = width, height = height, dpi = dpi, device = "png"
    )
  }
  
  list(obj = seurat_obj, plot = p, meta_col = meta_col)
}

## Reductions visualisation with cell annotations ##
UMAP_with_cell_annotations <- function(seurat_obj,
                                       annotation,                 
                                       reduction_id = "UMAP_RNA",
                                       reduction_title = "UMAP",
                                       facet_by = NULL,           
                                       annotation_name = NULL,    
                                       point_size = 0.1,
                                       point_alpha = 0.6,
                                       label_clusters = FALSE,
                                       label_size = 3,
                                       facet_ncol = 3,
                                       seed = 1511,
                                       save = FALSE,
                                       out_dir = NULL,
                                       filename = NULL,
                                       width = 10,
                                       height = 7,
                                       dpi = 300) {

  ## Get the path from the saved objects ##
  paths <- NULL
  if (is.list(annotation) && !is.null(annotation$annotation)) {
    paths <- annotation$paths
    ann   <- annotation$annotation
  } else {
    ann <- annotation
  }
  
  ## Keep it safe and make sure the requested values exist in the object 
  stopifnot(reduction_id %in% Seurat::Reductions(seurat_obj))
  stopifnot(is.list(ann), !is.null(ann$clust))
  
  ## Get the annotations from the list ##
  cl <- ann$clust
  if (is.null(names(cl))) {
    stop("annotation$clust must be a named vector (names = cell IDs) to align to Seurat cells.")
  }
  
  ## Annotation name in case one is not provided for automation ##
  if (is.null(annotation_name) || !nzchar(annotation_name)) {
    
    if (is.list(ann) && !is.null(ann$path) && nzchar(ann$path)) {
      annotation_name <- basename(dirname(ann$path))  # e.g. Kidney_HCA_supervised
      
    } else if (!is.null(paths) && !is.null(paths$main_tag) && nzchar(paths$main_tag)) {
      annotation_name <- paths$main_tag
      
    } else {
      annotation_name <- "InSituType"
    }
  }
  
  ## Get reduction coordinates ##
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction_id)
  if (ncol(emb) < 2) stop("Reduction '", reduction_id, "' has < 2 dimensions.")
  
  emb_df <- tibble::as_tibble(emb[, 1:2, drop = FALSE], rownames = "cell_ID")
  
  xcol <- paste0(reduction_title, "_1")
  ycol <- paste0(reduction_title, "_2")
  colnames(emb_df)[2:3] <- c(xcol, ycol)
  
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_ID")
  
  ## Attach cell labels ##
  emb_df_ann <- cl %>% 
    tibble(cell_ID = names(.), 
           Annotation = factor(unname(.))) %>% 
    select(cell_ID, Annotation) %>% 
    inner_join(emb_df, by = "cell_ID")
  
  ## Attach facet variable ##
  emb_df_ann_grp <- emb_df_ann
  
  if (!is.null(facet_by) && nzchar(facet_by)) {
    stopifnot(facet_by %in% colnames(meta_df))
    emb_df_ann_grp <- dplyr::left_join(emb_df_ann, 
                                       meta_df[, c("cell_ID", facet_by), drop = FALSE],
                                       by = "cell_ID")
  }
  
  ## Set title for the plot ##
  title_txt <- paste0(reduction_title, " projection: Based on ", 
                      annotation_name, " annotation" )
  
  ## Plot ##
  p <- ggplot2::ggplot(emb_df_ann_grp, ggplot2::aes(x = .data[[xcol]], 
                                                    y = .data[[ycol]], 
                                                    colour = Annotation)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::labs(title = title_txt,
                  x = paste0(reduction_title, " 1"),
                  y = paste0(reduction_title, " 2"),
                  colour = "Cell type") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.position = "right")
  
  p <- p +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1,
                                                   override.aes = list(size = 3, alpha = 1)))
  
  ## Optional facetting by specified group ##
  if (!is.null(facet_by) && nzchar(facet_by)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_by)), ncol = facet_ncol)
  }
  
  # optional cluster labels (centroids)
  if (isTRUE(label_clusters)) {
    centroids <- emb_df_ann_grp %>%
      dplyr::filter(!is.na(Annotation)) %>%
      dplyr::group_by(Annotation) %>%
      dplyr::summarise(x = stats::median(.data[[xcol]], na.rm = TRUE),
                       y = stats::median(.data[[ycol]], na.rm = TRUE),
                       .groups = "drop")
    
    p <- p +
      ggrepel::geom_label_repel(data = centroids,
                                ggplot2::aes(x = x, y = y, label = Annotation),
                                colour = "black",
                                fill = "white",
                                size = label_size,
                                label.size = 0.25,
                                show.legend = FALSE,
                                max.overlaps = Inf, 
                                seed = seed)
  }
  
  ## Save the plots next to the annotation RDS ##
  saved_path <- NULL
  if (isTRUE(save)) {
    
    if (is.null(out_dir) || !nzchar(out_dir)) {
      
      if (is.list(ann) && !is.null(ann$path) && nzchar(ann$path)) {
        out_dir <- dirname(ann$path)
        
      } else if (!is.null(paths) && !is.null(paths$ann_dir) && nzchar(paths$ann_dir)) {
        out_dir <- paths$ann_dir
        
      } else {
        out_dir <- file.path(main_dir, "Results/Annotations/Reduction_plots")
      }
    }
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (is.null(filename) || !nzchar(filename)) {
      filename <- paste0(annotation_name, "_", reduction_id, ".png")
      filename <- gsub("[^A-Za-z0-9._-]+", "_", filename)
    }
    
    saved_path <- file.path(out_dir, filename)
    
    ggplot2::ggsave(filename = saved_path,
                    plot = p,
                    width = width,
                    height = height,
                    dpi = dpi,
                    device = "png")
  }  
  
  list(plot = p,
       data = emb_df_ann_grp,
       annotation_name = annotation_name,
       reduction = reduction_id,
       facet_by = facet_by,
       saved_path = saved_path)
}

# Differential Expression visualisations ----------------------------------
## Rename the two Custom_VLN functions
## Custom vln plot for the markers ##
Custom_VLN <- function(obj, features, group_by = NULL, ncol = 3, pt.size = 0) {
  
  p <- VlnPlot(obj,
               features = features,
               group.by = group_by, 
               pt.size = pt.size,
               cols = NULL,       
               ncol = ncol) 
  
  p & theme_classic(base_size = 12) &
    theme(plot.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(face = "bold"),
          legend.position = "none")
  
}

## Custom vln plot for the markers ##
Custom_VLN <- function(obj, features, group_by = NULL, ncol = 1, pt.size = 0) {
  
  vln_list <- lapply(seq_along(features), function(i) {
    gene <- features[i]
    
    p <- VlnPlot(obj, 
                 features = gene, 
                 group.by = group_by, 
                 pt.size = pt.size,
                 cols = NULL) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            legend.position = "none", 
            plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 30))
    
    # Remove x text for all but the bottom panel
    if (i < length(features)) {
      p <- p +
        theme(axis.text.x  = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 30))
    } else {
      p <- p +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(l = 5)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 65))
    }
    
    p
  })
  
  # Stack: 1 column → N rows
  p_stack <- wrap_plots(vln_list, ncol = 1) 
  
  p_stack <- cowplot::ggdraw(p_stack) +
    cowplot::draw_label(label = "Expression Level",
                        x = 0.02,   
                        y = 0.6,    
                        angle = 90,
                        vjust = 0.5,
                        fontface = "bold") +
    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 65))
  
  p_stack
}

## Plot a number of the top markers by Log2FC ##
Plot_top_markers_VLN <- function(obj,
                                 markers_tbl,
                                 group_by,
                                 top_n = 3,
                                 out_dir = NULL,
                                 prefix = "Kidney_HCA", 
                                 pt_size = 0) {
  
  # 1. keep significant markers
  sig_markers <- markers_tbl %>%
    filter(p_val_adj < 0.05)
  
  # 2. top N genes per "cluster" (this column is created by FindAllMarkers)
  top_markers_by_cluster <- sig_markers %>%
    group_by(cluster) %>%
    slice_max(n = top_n, order_by = avg_log2FC) %>%
    ungroup()
  
  # list of features per cluster
  marker_list <- top_markers_by_cluster %>%
    group_split(cluster) %>%
    set_names(unique(top_markers_by_cluster$cluster)) %>%
    map(~ .x$gene)
  
  # 3. generate one violin plot per cluster
  plots <- imap(marker_list, ~ {
    clust <- .y
    feats <- .x
    
    p <- Custom_VLN(obj = obj,
                    features = feats,
                    group_by = group_by,
                    ncol = length(feats),
                    pt.size = pt_size) +
      patchwork::plot_annotation(title = paste0("Top markers for ", group_by, ": ", clust))
    
    # 4. optionally save to disk
    if (!is.null(out_dir)) {
      fname <- file.path(out_dir, paste0(prefix, "_Cluster_", clust, "_violin.png"))
      ggsave(fname, p, width = 16, height = 10, dpi = 300)
    }
    
    p
  })
  
  return(plots)  # list of ggplot objects
}

## VLN but custom for specific markers and celltypes ##
Plot_manual_markers_VLN <- function(obj,
                                    group_by,
                                    cell_type,
                                    markers,
                                    save_path = NULL,
                                    pt_size = 0) {
  
  # make the stacked violins for the chosen markers
  p <- Custom_VLN(obj = obj,
                  features = markers,
                  group_by = group_by,
                  ncol = 1,
                  pt.size = pt_size) +
    patchwork::plot_annotation(title = paste0("Marker expression for ", group_by, ": ", cell_type)) 
  
  # optionally save to disk
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, p, width = 16, height = 10, dpi = 300)
  }
  
  return(p)
}



## Make into function ##
# B_plots <- FeaturePlot(
#   Seurat_obj,
#   features   = B_Markers,
#   reduction  = "umap",
#   pt.size    = 0.4,
#   min.cutoff = "q05",
#   max.cutoff = "q95",
#   cols       = c("grey90", "midnightblue"),
#   combine    = FALSE           # <- important
# )
# 
# # 2. Clean theme for each plot
# B_plots <- lapply(seq_along(B_plots), function(i) {
#   B_plots[[i]] +
#     theme_void(base_size = 12) +
#     theme(
#       plot.title   = element_text(face = "bold", hjust = 0.5),
#       legend.position = "right",
#       legend.title    = element_text(face = "bold"),
#       legend.text     = element_text(size = 8),
#       plot.margin     = margin(t = 5, r = 5, b = 5, l = 5)
#     ) +
#     labs(title = B_Markers[i], colour = "Expression")
# })
# 
# # 3. Arrange in a grid and add a single global title + single legend
# B_feature_grid <- wrap_plots(B_plots, ncol = 3, guides = "collect") +
#   plot_annotation(
#     title = "B cell markers expression on the UMAP projection"
#   ) &
#   theme(
#     legend.position = "right",
#     legend.title    = element_text(face = "bold"),
#     legend.text     = element_text(size = 8)
#   )
# 
# B_feature_grid
# 
# ggsave(file.path(main_dir, "Results/DE_Plots/FeaturePlots/B_markers_UMAP_grid.png"),
#        B_feature_grid,
#        width  = 16,
#        height = 10,
#        dpi    = 300, device = "png"
# )




Plot_all_sig_genes_VLN <- function(obj,
                                   markers_tbl,
                                   group_by,
                                   annotation_name,
                                   out_dir_base,
                                   p_adj_cutoff = 0.05,
                                   avg_log2FC_cutoff = 1.5,
                                   pt_size = 0,
                                   overwrite = FALSE) {
  
  ## Significant thresholds ##
  sig_tbl <- markers_tbl %>%
    dplyr::filter(p_val_adj <= p_adj_cutoff &
                  avg_log2FC > avg_log2FC_cutoff) %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val_adj) %>%
    dplyr::distinct()
  
  ## Output directory ##
  out_dir_annotation <- file.path(out_dir_base, annotation_name)
  dir.create(out_dir_annotation, showWarnings = FALSE, recursive = TRUE)
  
  ## Choose clusters ##
  clusters <- sort(unique(sig_tbl$cluster))
  
  ## Initiate loop ##
  for (cl in clusters) {
    message("Cluster: ", cl)
    
    ## Create subfolder per cluster ##
    out_dir_cluster <- file.path(out_dir_annotation, cl)
    dir.create(out_dir_cluster, showWarnings = FALSE, recursive = TRUE)
    
    ## Get the significant genes for the current cluster ##
    genes_cl <- sig_tbl %>%
      dplyr::filter(cluster == cl) %>%
      dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
      dplyr::pull(gene) %>%
      unique(gene, avg_log2FC)
    
    ## Visualisation loop ##
    for (i in seq_len(nrow(genes_cl))) {
      
      g <- genes_cl$gene[i]
      fc <- genes_cl$avg_log2FC[i]
      
      ## Make sure the gene names contain acceptable characters ##
      g_safe <- gsub("[^A-Za-z0-9._-]", "_", g)
      
      ## Create the filenames ##
      fname <- paste0(annotation_name, "_", cl, "_", g_safe, "_Vln.png")
      fpath <- file.path(out_dir_cluster, fname)
      
      ## If it already exists and overwrite is set to FALSE skip ##
      if (file.exists(fpath) && !overwrite) next
      
      ## Create the current gene's plot ##
      p <- Custom_VLN(obj = obj,
                      features = c(g),
                      group_by = group_by,
                      ncol = 1,
                      pt.size = pt_size) +
        patchwork::plot_annotation(title = paste0(annotation_name, ": cluster ", cl, " – ", g, 
                                                  " (avg_log2FC=", sprintf("%.3f", fc), ")"))
      
      ggsave(fpath, p, width = 16, height = 9, dpi = 300)
      
    }
  }
}

