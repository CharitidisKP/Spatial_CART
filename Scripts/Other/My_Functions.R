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

## Function to visualise cell annotation results ##
Custom_flighpath <- function(flightpath_obj,
                             annotation_obj,
                             unsupervised = FALSE,
                             title = NULL, 
                             plot_seed = 1511) {
  
  ## Get object name simplify for the title ##
  if (unsupervised == TRUE) {
    title_text <- paste0("InSituType annotation: Unsupervised")
  } else {
    title_text <- paste0("InSituType annotation using the: ", title, " reference profiles")
  }
  
  Cluster_stats <- tibble(Cluster = flightpath_obj$clust) %>%
    group_by(Cluster) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(mean_conf = flightpath_obj$meanconfidence, 
           Cluster_lab = sprintf("%s (n = %s, mean conf = %.2f)", Cluster, 
                                 scales::comma(n_cells), mean_conf), 
           Cluster = factor(Cluster)) %>% 
    arrange(desc(n_cells))
  
  ## Plot the cells ##
  Celltypes_df <- flightpath_obj$cellpos %>% 
    as.data.frame() %>% 
    mutate(Cluster = factor(annotation_obj$clust)) %>% 
    left_join(Cluster_stats, by = "Cluster") %>%
    mutate(Cluster_lab = factor(Cluster_lab, levels = Cluster_stats$Cluster_lab))
  
  ## Cluster label positions ##
  Celltypes_pos_df <- flightpath_obj$clustpos %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Cluster") %>% 
    mutate(Cluster = factor(Cluster)) %>% 
    left_join(Cluster_stats, by = "Cluster") %>%
    mutate(Cluster_lab = factor(Cluster_lab, levels = Cluster_stats$Cluster_lab))
  
  set.seed(plot_seed)
  ## Scatter plot ##
  ggplot(Celltypes_df, aes(x = x, y = y, colour = Cluster_lab)) +
    geom_point(size = 0.2, alpha = 0.8) +
    geom_label_repel(data = Celltypes_pos_df, aes(x = x, y = y, label = Cluster),
                     colour = "black", fill = "white", label.size = 0.2, fontface = "bold",
                     label.r = unit(0.1, "lines"), size = 3, show.legend = FALSE) +
    theme_void() +
    labs(title = title_text, 
         colour = paste0("InSituType Clusters")) +
    guides(colour = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme(legend.position = "right", 
          plot.margin = margin(0, 0, 0, 0), 
          plot.title = element_text(face = "bold", hjust = 0.5), 
          legend.title = element_text(face = "bold", size = 11), 
          legend.text = element_text(face = "bold", size = 9)) 
  
}

## Helper for the pheatmap visualisations ##
## Scales matrix expression rows and drops genes who dont meet the threshold ##
Scale_rows_pheatmap <- function(mat, min_max = 0.2, min_keep = 0.2) {
  row_max <- matrixStats::rowMaxs(mat)
  scaled <- mat / pmax(row_max, min_max)
  scaled[row_max > min_keep, , drop = FALSE]
}

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

## FeaturePlot CUstom ##
B_plots <- FeaturePlot(Seurat_obj,
                       features = B_Markers,
                       reduction = "umap",
                       pt.size = 1,
                       min.cutoff = "q05",
                       max.cutoff = "q95",
                       cols = c("grey30", "firebrick"),
                       combine = FALSE)

# tweak each plot: remove axes, center title, small legend
B_plots <- lapply(seq_along(B_plots), function(i) {
  p <- B_plots[[i]] +
    theme_void(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(face = "bold"),
      legend.text     = element_text(size = 8),
      plot.margin     = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    labs(title = B_Markers[i], colour = "Expression")
  
  p
})
