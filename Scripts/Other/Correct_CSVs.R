## Load the libraries and the sliding function ##
library(tidyverse)
source(file = "My_Functions.R")

## Set the working directory ##
setwd("~/P_lab/Spatial_CART/")

## Load the "raw" files ##
Raw_Counts <- read_csv("Data/ASCLST12502_rawcount_exprMatrix.csv")
Metadata <- read_csv("Data/ASCLST12502_metadata.csv")

## Reorder the columns ##
Metadata <- Correct_Dataframes(Metadata) %>%
  select(cell_id, fov, everything()) %>% 
  dplyr::rename(Cell_ID = cell_id)

## Save the Cell IDs ##
Cell_IDs <- Raw_Counts[1]

## Apply the function ## 
Raw_Counts <- Correct_Dataframes(Raw_Counts)

## Add the names 
Raw_Counts["Cell_ID"] <- Cell_IDs

Raw_Counts <- Raw_Counts %>%
  select(Cell_ID, everything())

# write_csv(x = Raw_Counts, file = "Data/ASCLST12502_rawcount_fixed.csv")
# write_csv(x = Metadata, file = "Data/ASCLST12502_metadata_fixed.csv")
