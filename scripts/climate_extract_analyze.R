### climate extraction, elevation, pca!! ###
# Anri Chomentowska, 2025 #

# set up 
rm(list = ls()) # clear environment
setwd("~/Dropbox/Yale/Research/Map") #change

# libraries
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(gtools)
library(sp)
library(raster)
library(maps)
library(mapproj)
library(ENMTools)
library(dismo)
library(corrr)
library(tidyverse)
library(caret)
library(ggfortify)
library(cluster)
library(elevatr)
#library(rgdal)
#library(rgeos)

### extract elevation ###

# load gbif data
CCM <- read.csv("2025_gbif_clean_cropped.csv")
CCM <- data.frame(CCM[,c(2,3,6,7,15)])

View(CCM)

# get elevation
data <- data.frame(x = CCM$longitude, y = CCM$latitude, names = CCM$species)
prj_dd <- "EPSG:4326"
CCM_elevation <- get_elev_point(data, prj = prj_dd, z = 5, src = "aws")

CCM_elevation_df <- CCM %>%
  mutate(elevation = CCM_elevation$elevation)

View(CCM_elevation_df)
summary(CCM_elevation_df)
table(CCM_elevation_df$species)

# write output
write.csv(CCM_elevation_df, "2025_gbif_elevation.csv", row.names=F)

# add non-gbif locality data from our own specimens
add_CCM <- read.csv("localities_from_specimen.csv")
all_CCM <- bind_rows(CCM_elevation_df, add_CCM)

View(all_CCM)

# map
mymap <- borders("world", xlim = c(-130, -60), ylim = c(-50, 50),
                 colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ mymap +
  geom_point(data = all_CCM, aes(x = longitude,
                                 y = latitude,
                                 color = genus),
             size = 0.5)+
  theme_bw()

# clean the names column
all_CCM <- all_CCM %>%
  mutate(names = trimws(as.character(species)))

# check if names are cleaned properly
print(unique(all_CCM$species))

# summarize elevation data using aggregate()
summary_base <- aggregate(elevation ~ species, data = all_CCM, 
                          FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                              min = min(x, na.rm = TRUE), 
                                              max = max(x, na.rm = TRUE)))

# convert the result to a data frame
summary_base <- do.call(data.frame, summary_base)

# rename columns for clarity
colnames(summary_base) <- c("names", "avg_elevation", "min_elevation", "max_elevation")

# view summarized data
print(summary_base)
View(summary_base)
View(all_CCM)

# save output
write.csv(all_CCM, "2025_CCM_elevation.csv", row.names=F)
write.csv(summary_base, "2025_CCM_elevation_summary.csv", row.names=F)


### now let's get the climate data! ###

# first I need to sort the data frame by species
#all_CCM <- read.csv("2025_CCM_elevation.csv")
data_sorted <- arrange(all_CCM, species)

View(data_sorted)

# define species extent and map
species_extent <- raster::extent(c(min(data_sorted$longitude)-10,
                                   max(data_sorted$longitude)+10),
                                 c(min(data_sorted$latitude)-10,
                                   max(data_sorted$latitude)+10))

mymap <- borders("world", xlim = c(-134, -55), ylim = c(-62, 60),
                 colour=NA, fill="gray30")

ggplot()+ coord_fixed()+ mymap

## download bioclim data at 2.5degree resolution (SKIP IF ALREADY DONE)
Env <- raster::getData("worldclim", var="bio", res=2.5) ##comes in as global by default
Env = raster::crop(Env,species_extent)

# see: https://www.worldclim.org/data/bioclim.html for variable name explanations
rasterVis::splom(Env,varname.cex=0.8)
Env = raster::dropLayer(Env, c("bio6", "bio7", "bio9", "bio10", "bio11", "bio14", "bio16", "bio17", "bio18", "bio19"))

## load bioclim variables (if not downloading for first time)
bio1_l <- raster("wc2-5/bio1.bil")
bio2_l <- raster("wc2-5/bio2.bil")
bio3_l <- raster("wc2-5/bio3.bil")
bio4_l <- raster("wc2-5/bio4.bil")
bio5_l <- raster("wc2-5/bio5.bil")
bio6_l <- raster("wc2-5/bio6.bil")
bio7_l <- raster("wc2-5/bio7.bil")
bio8_l <- raster("wc2-5/bio8.bil")
bio9_l <- raster("wc2-5/bio9.bil")
bio10_l <- raster("wc2-5/bio10.bil")
bio11_l <- raster("wc2-5/bio11.bil")
bio12_l <- raster("wc2-5/bio12.bil")
bio13_l <- raster("wc2-5/bio13.bil")
bio14_l <- raster("wc2-5/bio14.bil")
bio15_l <- raster("wc2-5/bio15.bil")
bio16_l <- raster("wc2-5/bio16.bil")
bio17_l <- raster("wc2-5/bio17.bil")
bio18_l <- raster("wc2-5/bio18.bil")
bio19_l <- raster("wc2-5/bio19.bil")


## crop bioclim layers to the desired extent
# visualize the first layer
plot(bio1_l)

## crop the bioclim layers to extent of the range, then mask the bioclim layer with the range.
bio1 <- raster::crop(bio1_l, species_extent)
bio2 <- raster::crop(bio2_l, species_extent)
bio3 <- raster::crop(bio3_l, species_extent)
bio4 <- raster::crop(bio4_l, species_extent)
bio5 <- raster::crop(bio5_l, species_extent)
bio6 <- raster::crop(bio6_l, species_extent)
bio7 <- raster::crop(bio7_l, species_extent)
bio8 <- raster::crop(bio8_l, species_extent)
bio9 <- raster::crop(bio9_l, species_extent)
bio10 <- raster::crop(bio10_l, species_extent)
bio11 <- raster::crop(bio11_l, species_extent)
bio12 <- raster::crop(bio12_l, species_extent)
bio13 <- raster::crop(bio13_l, species_extent)
bio14 <- raster::crop(bio14_l, species_extent)
bio15 <- raster::crop(bio15_l, species_extent)
bio16 <- raster::crop(bio16_l, species_extent)
bio17 <- raster::crop(bio17_l, species_extent)
bio18 <- raster::crop(bio18_l, species_extent)
bio19 <- raster::crop(bio19_l, species_extent)

# visualize first layer again
plot(bio1)

# save the layers
setwd("~/Dropbox/Yale/Research/Map/bioclim_output") #change
writeRaster(bio1, "bio1.asc", format="ascii", overwrite=TRUE)
writeRaster(bio2, "bio2.asc", format="ascii", overwrite=TRUE)
writeRaster(bio3, "bio3.asc", format="ascii", overwrite=TRUE)
writeRaster(bio4, "bio4.asc", format="ascii", overwrite=TRUE)
writeRaster(bio5, "bio5.asc", format="ascii", overwrite=TRUE)
writeRaster(bio6, "bio6.asc", format="ascii", overwrite=TRUE)
writeRaster(bio7, "bio7.asc", format="ascii", overwrite=TRUE)
writeRaster(bio8, "bio8.asc", format="ascii", overwrite=TRUE)
writeRaster(bio9, "bio9.asc", format="ascii", overwrite=TRUE)
writeRaster(bio10, "bio10.asc", format="ascii", overwrite=TRUE)
writeRaster(bio11, "bio11.asc", format="ascii", overwrite=TRUE)
writeRaster(bio12, "bio12.asc", format="ascii", overwrite=TRUE)
writeRaster(bio13, "bio13.asc", format="ascii", overwrite=TRUE)
writeRaster(bio14, "bio14.asc", format="ascii", overwrite=TRUE)
writeRaster(bio15, "bio15.asc", format="ascii", overwrite=TRUE)
writeRaster(bio16, "bio16.asc", format="ascii", overwrite=TRUE)
writeRaster(bio17, "bio17.asc", format="ascii", overwrite=TRUE)
writeRaster(bio18, "bio18.asc", format="ascii", overwrite=TRUE)
writeRaster(bio19, "bio19.asc", format="ascii", overwrite=TRUE)

## load all of the cropped climate layers
setwd("~/Dropbox/Yale/Research/Map") #change
outputDir <- "bioclim_output"
list <- list.files(outputDir, full.names = T, recursive = FALSE)
list <- mixedsort(sort(list))
envtStack <- stack(list)

View(data_sorted)

# ust the lat long
bioclim <- raster::extract(envtStack, data_sorted[, c("longitude", "latitude")])
View(bioclim)

# add specific points from our latest collections
lenzia <- data.frame(
  longitude = -69.80167,
  latitude = -30.19055)
bioclim_lenzia <- raster::extract(envtStack, df)
pachy <- data.frame(
  longitude = -71.66767,
  latitude = -17.21354)
bioclim_pachy <- raster::extract(envtStack, pachy)

data_combined <- data.frame(data_sorted, bioclim)
View(data_combined)

# filter out NA values
data_combined_clean <- data_combined %>% 
  drop_na(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)

View(data_combined_clean)

# write the output
write.csv(data_combined_clean, "2025_data_elev_bioclim.csv")

## let's summarize again
setwd("~/Dropbox/Yale/Research/Map") #change 
data_combined_clean <- read.csv("2025_data_elev_bioclim.csv")

View(data_combined_clean)

# remove weird negative sea level samples
data_combined_clean <- data_combined_clean %>%
  filter(elevation >= 0)

# list all bioclim variables + elevation to summarize
vars_to_summarize <- c("elevation", paste0("bio", 1:19))

# function to summarize a single variable by species
summarize_by_species <- function(var) {
  summary <- aggregate(data_combined_clean[[var]] ~ data_combined_clean$species, 
                       FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                           min = min(x, na.rm = TRUE), 
                                           max = max(x, na.rm = TRUE)))
  
  summary <- do.call(data.frame, summary)
  colnames(summary) <- c("species", 
                         paste0(var, "_mean"), 
                         paste0(var, "_min"), 
                         paste0(var, "_max"))
  return(summary)
}

# apply to each variable and merge results
summary_list <- lapply(vars_to_summarize, summarize_by_species)

# merge all summaries by species
summary_all <- Reduce(function(x, y) merge(x, y, by = "species"), summary_list)

# let's check
View(summary_all)

#write results
write.csv(summary_all, "2025_bioclim_summary_by_species_FINAL.csv", row.names = FALSE)

## a quick detour to see if there is meaningful difference in using all data vs IQR
# get elevation IQR bounds per species
iqr_bounds <- data_combined_clean %>%
  group_by(species) %>%
  summarize(elev_q1 = quantile(elevation, 0.25, na.rm = TRUE),
            elev_q3 = quantile(elevation, 0.75, na.rm = TRUE),
            groups = "drop")

# keep only rows within each species' elevation IQR
data_in_iqr <- data_combined_clean %>%
  left_join(iqr_bounds, by = "species") %>%
  filter(elevation >= elev_q1, elevation <= elev_q3)

# summarize means for elevation + all bioclim variables
summary_all_iqr <- data_in_iqr %>%
  group_by(species) %>%
  summarize(
    across(all_of(vars_to_summarize),
           ~ mean(.x, na.rm = TRUE),
           .names = "{.col}_mean"),
    .groups = "drop"
    )

View(summary_all_iqr)

write.csv(summary_all_iqr,
          "2025_bioclim_summary_by_species_IQR_means.csv",
          row.names = FALSE)

# load if not continuing from previous steps
summary_all <- read.csv("2025_bioclim_summary_by_species_FINAL.csv")
View(summary_all)

# join the two summary tables
compare_df <- summary_all %>%
  inner_join(summary_all_iqr, by = "species", suffix = c("_full", "_iqr"))

compare_means <- compare_df %>%
  select(species, ends_with("_mean_full"), ends_with("_mean_iqr"))

# pivot
compare_long <- compare_means %>%
  pivot_longer(cols = -species,
               names_to = c("variable", "method"),
               names_pattern = "(.*)_mean_(full|iqr)",
               values_to = "value") %>%
  pivot_wider(id_cols = c(species, variable),
              names_from = method,
              values_from = value)

View(compare_long)

# let's graph
ggplot(compare_long, aes(x = full, y = iqr)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  labs(x = "Mean from full dataset",
       y = "Mean from IQR-filtered dataset",
       title = "Comparison of full vs IQR-trimmed means per species")


### let's summarize and visualize ###

# load data that you produced above
setwd("~/Dropbox/Yale/Research/Map") #change
data_all <- read.csv("2025_data_elev_bioclim.csv") #not yet 'cleaned', i.e. matches phylogeny 
View(data_all)

# load libraries (if starting here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gtools)
library(tidyverse)
library(ggridges)
library(ggrepel)

## clean to match trimmed phylogeny
# do this for the summary data output above in excel too! (weighted average across cymosa / longiscapa)
# and keep just the means (not max min) == FINAL LH and bioclim matrix

# but the following is for visualization purposes
data_all_clean <- data_all %>%
  mutate(
    elevation = as.numeric(elevation),  # Ensure numeric
    species = case_when(
      species %in% c("Cistanthe cymosa", "Cistanthe longiscapa") ~ "Cistanthe aff longiscapa",
      TRUE ~ species
    )
  ) %>%
  filter(!is.na(elevation), elevation >= 0) %>%  # Remove NAs and values < 0
  filter(!species %in% c("Cistanthe thyrsoidea", "Cistanthe arenaria"))

# pivot
bioclim_data_long <- data_all_clean %>%
  dplyr::select(species, 
                elevation,
                bio12
                #starts_with("bio")
  ) %>%
  pivot_longer(cols = -species, names_to = "variable", values_to = "value")

# load life history data and join
life_hist_df <- read.csv("2025_Species_LH.csv", stringsAsFactors = FALSE)
View(life_hist_df)

bioclim_data_annotated <- bioclim_data_long %>%
  left_join(life_hist_df, by = "species")

# just get elevation
elevation_long <- bioclim_data_annotated %>%
  filter(variable == "elevation")

View(elevation_long)

# re-order
View(life_hist_df)
ordered_species <- life_hist_df$species

# make species in elevation_long a factor with order
elevation_long <- elevation_long %>%
  mutate(species = factor(species, levels = ordered_species))

# count number of specimens (occurrence points) per species
n_specimens_df <- elevation_long %>%
  group_by(species) %>%
  summarise(n_specimens = n(), .groups = "drop")

# sanity check
hist(elevation_long$value, breaks = 50)
summary(elevation_long$value) 

# get IQR
iqr_per_species <- elevation_long %>%
  group_by(species, life_history) %>%
  summarize(
    #elevation_IQR = IQR(value, na.rm = TRUE),
    elevation_IQR = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
    .groups = "drop"
  )
View(iqr_per_species)

# plot
ggplot(iqr_per_species, aes(x = life_history, y = elevation_IQR)) +
  geom_boxplot(aes(color = life_history),
               outlier.shape = NA,
               width = 0.5,
               alpha = 0.4) +
  geom_jitter(aes(color = life_history),
              width = 0.1,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Life History Strategy",
       y = "Elevation range (m)",
       title = "Range of Elevation per Species",
       color = "Life History") +
  theme_minimal(base_size = 13)

# species labels to check
geom_text_repel(aes(label = species),
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.3,
                point.padding = 0.2,
                segment.size = 0.2)

## let's look at statistical significance
wilcox.test(elevation_IQR ~ life_history, data = iqr_per_species)

model <- lm(elevation_IQR ~ life_history, data = iqr_per_species)
summary(model)

## more stats
elev_stats <- elevation_long %>%
  group_by(species, life_history) %>%
  summarize(elevation_mean = mean(value, na.rm = TRUE),
            elevation_IQR = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
            .groups = "drop")

elev_stats <- elevation_long %>%
  group_by(species, life_history) %>%
  summarize(elevation_mean = mean(value, na.rm = TRUE),
            elevation_IQR = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(n_specimens_df, by = "species")

View(elev_stats)
summary(elev_stats$n_specimens)

# linear model
sampling_model <- lm(elevation_IQR ~ n_specimens, data = elev_stats)
summary(sampling_model)

# plots
ggplot(elev_stats, aes(x = n_specimens, y = elevation_IQR)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "Number of specimens (GBIF records)",
       y = "Elevation (m)",
       title = "Does sampling effort affect estimated elevational range?") +
  theme_minimal()

ggplot(elev_stats, aes(x = elevation_mean, y = elevation_IQR)) +
  geom_point(aes(color = life_history), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", aes(color = life_history), #color = "black", #formula = y ~ poly(x, 2),
              se = TRUE, linetype = "dotted") +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Mean Elevation (m)",
       y = "Elevation (m)",
       color = "Life History") +
  theme_minimal(base_size = 13)

## MORE stats!
# annual model
annual_model <- lm(elevation_IQR ~ elevation_mean, data = elev_stats, 
                   subset = life_history == "Annual")
summary(annual_model)

# perennial model
perennial_model <- lm(elevation_IQR ~ elevation_mean, data = elev_stats, 
                      subset = life_history == "Perennial")
summary(perennial_model)

#model interaction
model_interaction <- lm(elevation_IQR ~ elevation_mean * life_history, data = elev_stats)
summary(model_interaction)

#model quadratic formula
model_quad <- lm(elevation_IQR ~ poly(elevation_mean, 2), data = elev_stats)
summary(model_quad)

#write outputs
library(broom)
# annual model coefficients
write.csv(tidy(annual_model), "Annual_Model_Coefficients.csv", row.names = FALSE)
# interaction model coefficients
write.csv(tidy(model_interaction), "Interaction_Model_Coefficients.csv", row.names = FALSE)
# quadratic model coefficients
write.csv(tidy(model_quad), "Quadratic_Model_Coefficients.csv", row.names = FALSE)

## PGLS
View(elev_stats)

elev_stats_df <- as.data.frame(elev_stats)
View(elev_stats_df)

elev_stats_df$species <- gsub(" ", "_", elev_stats_df$species)

# load tree
# SCRIPT TO LOAD TREE from phylogenetic_analyses.R

setdiff(elev_stats_df$species, trim_tree_nodeless$tip.label)
setdiff(trim_tree_nodeless$tip.label, elev_stats_df$species)

# with the tree, prep data
comp_data <- comparative.data(phy = trim_tree_nodeless,
                              data = elev_stats_df,
                              names.col = "species",
                              vcv = TRUE,
                              na.omit = FALSE)

# convert life history to numeric
comp_data$data$LH_numeric <- ifelse(comp_data$data$life_history == "Perennial", 1, 0)

# run PGLS
pgls_model <- pgls(elevation_IQR ~ LH_numeric, data = comp_data)

# view summary
summary(pgls_model)

## summary plots
#ridgeplot
ggplot(elevation_long, 
       aes(x = value, 
           y = species, 
           fill = life_history)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01, alpha = 0.4) +
  geom_jitter(aes(color = life_history),
              size = 0.5, 
              height = 0.2, 
              alpha = 0.5,
              width = 0.15) +
  scale_x_continuous(limits = c(0, 5000)) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Elevation (m)", y = "Species",
       fill = "Life History",
       title = "Elevation distributions across species (ridge plot)") +
  theme_minimal(base_size = 13)

#violin
ggplot(elevation_long, 
       aes(x = value, y = species, fill = life_history)) +
  geom_violin(scale = "width",
              trim = TRUE,
              adjust = 1.5,
              draw_quantiles = NULL,
              alpha = 0.4) +
  geom_jitter(aes(color = life_history),
              size = 0.5, 
              height = 0.2, 
              alpha = 0.5,
              width = 0.15) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Elevation (m)", y = "Species",
       fill = "Life History", color = "Life History",
       title = "Elevation distributions across species") +
  theme_minimal(base_size = 13)

#boxplot
ggplot(elevation_long, 
       aes(x = species, 
           y = value,
           fill = life_history)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3, alpha = 0.6) +
  geom_jitter(aes(color = life_history),
              size = 0.5,
              height = 0.2,
              width = 0.15,
              alpha = 0.5) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  coord_flip() +
  labs(x = "Species", y = "Elevation (m)",
       fill = "Life History",
       title = "Elevation by species (phylogenetic order)") +
  theme_minimal(base_size = 13)

#half-violin plots (different attempts)
library(ggdist)

ggplot(elevation_long, 
       aes(x = value, y = species, fill = life_history)) +
  stat_halfeye(adjust = 1.5,
               justification = -0.3,
               .width = 0.8,
               point_alpha = 0,
               alpha = 0.6) +
  geom_jitter(aes(color = life_history),
              size = 0.5,
              height = 0.2,
              width = 0,
              alpha = 0.5) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Elevation (m)", y = "Species",
       fill = "Life History", color = "Life History",
       title = "Elevation distributions (half-violin + raw points)") +
  theme_minimal(base_size = 13)

ggplot(elevation_long, 
       aes(x = value, y = species, fill = life_history)) +
  geom_violin(scale = "width",
              trim = TRUE,
              adjust = 1.5,
              draw_quantiles = NULL,
              alpha = 0.4) +
  geom_jitter(aes(color = life_history),
              size = 0.5, 
              height = 0.2, 
              alpha = 0.8,
              width = 0) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Elevation (m)", y = "Species",
       fill = "Life History", color = "Life History",
       title = "Elevation distributions across species (half-violin style)") +
  theme_minimal(base_size = 13)

library(gghalves)

ggplot(elevation_long, 
       aes(x = value, 
           y = species, 
           fill = life_history)) +
  geom_half_violin(side = "l",
                   trim = FALSE,
                   adjust = 1.5,
                   alpha = 0.6,
                   color = NA) +
  geom_half_point(aes(color = life_history),
                  side = "r",
                  transformation = position_jitter(height = 0.1, width = 0),
                  size = 0.5,
                  alpha = 0.5) +
  scale_fill_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  scale_color_manual(values = c("Perennial" = "#377EB8", "Annual" = "#E41A1C")) +
  labs(x = "Elevation (m)", y = "Species",
       fill = "Life History", color = "Life History",
       title = "Elevation distributions across species (half-violin + raw data)") +
  theme_minimal(base_size = 13)
