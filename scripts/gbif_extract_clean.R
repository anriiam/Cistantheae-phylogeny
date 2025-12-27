### Extracting and cleaning occurence points for the Cistantheae clade! ###
# Anri Chomentowska, 2025 #

# set up 
rm(list = ls()) # clear environment
setwd("~/Dropbox/Yale/Research/Map") #change

# libraries
library(tidyverse)
library(CoordinateCleaner)
library(countrycode)
library(ggplot2)
library(spocc)
library(dismo)


### EXTRACTING ###
# select the genera I want from gbif website and download occ data
df_genera <- read.delim("gbif_download_03262025.csv", header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
View(df_genera)

# check NA species and make sure ok with getting rid of 
no_cn <- df_genera %>%
  filter(is.na(species))
View(no_cn)

# subset to the columns I want
df_genera <- df_genera %>%
  dplyr::select(verbatimScientificName,
                decimalLongitude,
                decimalLatitude,
                basisOfRecord,
                family,
                genus,
                species,
                infraspecificEpithet,
                scientificName,
                stateProvince,
                year,
                month,
                day,
                eventDate, 
                elevation)

View(df_genera)
table(df_genera$species)
barplot(table(df_genera$genus))

### CLEANING ###
# remove NA species
df_genera <- df_genera %>% 
  filter(!is.na(species) & species != "") %>%
  filter(decimalLongitude != is.na(decimalLongitude)) %>%
  filter(decimalLatitude != is.na(decimalLatitude))

View(df_genera)

# flag problems
flags <- CoordinateCleaner::clean_coordinates(x = df_genera,
                                              lon = "decimalLongitude",
                                              lat = "decimalLatitude",
                                              species = "species",
                                              tests = c("centroids", 
                                                        "equal",
                                                        "gbif", 
                                                        "institutions",
                                                        "zeros", 
                                                        "seas"))

# exclude problematic records (317)
df_genera_2 <- df_genera[flags$.summary,]
rm(flags)

View(df_genera_2)

# subset to only preserved specimen 
df_genera_3 <- df_genera_2 %>%
  filter(basisOfRecord == "PRESERVED_SPECIMEN")

df_genera_3 <- df_genera_3 %>%
  rename(
    longitude = decimalLongitude,
    latitude = decimalLatitude)

View(df_genera_3)

## take a look at species distribution
# create data for world coordinates using map_data() function 
world_coordinates <- map_data("world") 
wrl_map <- borders("world", colour=NA, fill="gray30")

base_map <- ggplot()+ coord_fixed()+ wrl_map +
  geom_point(data = df_genera_3,
             aes(x = longitude, y = latitude,
                 color = genus),
             size = 0.74) +
  theme(legend.position.inside = c(.40, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(fill=NA),
        legend.key = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

base_map

# delete the ones in europe/new zealand, etc
df_delete <- filter(df_genera_3,longitude > -50)
df_genera_4 <- setdiff(df_genera_3, df_delete)
rm(df_delete)

table(df_genera_4$species)

# filter out records with no year
df_genera_4 <- df_genera_4 %>%
  filter(year != is.na(year))

View(df_genera_4)
barplot(table(df_genera_4$year))

# filter out records older than 1900
df_genera_4 <- df_genera_4 %>% 
  filter(year > 1900)

view(df_genera_4)
table(df_genera_4$species)

# filter out duplicate data
df_genera_4 <- df_genera_4 %>%
  distinct(longitude, latitude, .keep_all=TRUE)

view(df_genera_4)
table(df_genera_4$species)

## make a plot again to see which points still shouldn't be there
base_map <- ggplot()+ coord_fixed()+ wrl_map +
  geom_point(data = calandrinia_clean,
             aes(x = longitude, y = latitude,
                 color = genus),
             size = 0.74) +
  theme(legend.position.inside = c(.40, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(fill=NA),
        legend.key = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

base_map

## get rid of the point in the middle of the US and in New England
# Filter problematic genera with spatial logic
cistanthe_clean <- df_genera_4 %>%
  filter(genus == "Cistanthe") %>%
  filter(!(longitude > -100 & latitude > 15))  # remove midwest + Caribbean points

montiopsis_clean <- df_genera_4 %>%
  filter(genus == "Montiopsis") %>%
  filter(!(latitude > 0))  # remove points in the Northern Hemisphere

calandrinia_clean <- df_genera_4 %>%
  filter(genus == "Calandrinia") %>%
  filter(!(longitude > -100 & latitude > 25))  # remove east Texas points

# Keep all other genera as-is
others_clean <- df_genera_4 %>%
  filter(!genus %in% c("Cistanthe", "Montiopsis", "Calandrinia"))

# Recombine everything
df_cleaned <- bind_rows(cistanthe_clean,
                        montiopsis_clean,
                        calandrinia_clean,
                        others_clean)

View(df_cleaned)

# make a plot again to see which points still shouldn't be there
base_map <- ggplot()+ coord_fixed()+ wrl_map +
  geom_point(data = df_cleaned,
             aes(x = longitude, y = latitude,
                 color = genus),
             size = 0.74) +
  theme(legend.position.inside = c(.40, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(fill=NA),
        legend.key = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

base_map

## let's take a look at species name and update if need be
df_look <- table(df_cleaned$species)
df_look

df_cleaned <- df_cleaned %>%
  mutate(species = ifelse(species == "Cistanthe pachyphylla",
                          "Philippiamra pachyphylla",
                          species))

spnames <- read.csv("name_list.csv")

df_cleaned_crop <- df_cleaned %>%
  filter(species %in% spnames$Name)

View(df_cleaned_crop)

df_look_again <- table(df_cleaned_crop$species)
df_look_again

setdiff(spnames$Name, df_cleaned$species)

## save final output
write.csv(df_cleaned, "2025_gbif_clean.csv", row.names=F)
write.csv(df_cleaned_crop, "2025_gbif_clean_cropped.csv", row.names=F)

## final map!
library(RColorBrewer)
mymap <- borders("world", xlim = c(-130, -20), ylim = c(-60, 50),
                 colour=NA, fill="gray30")

base_map <- ggplot()+ coord_fixed()+ mymap +
  geom_point(data = df_cleaned_crop %>%
               filter(!genus %in% c("Lewisiopsis", "Calandrinia", "Phemeranthus")),
             aes(x = longitude, y = latitude,
                 color = genus),
             size = 0.74) +
  scale_color_brewer(palette = "Set3") +
  theme(legend.position = c(.40, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(fill=NA),
        legend.key = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

base_map

## visualize genus counts
barplot(table(df_cleaned_crop$genus))


