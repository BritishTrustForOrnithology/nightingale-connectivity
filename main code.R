############################################################
# PLOTS FOR MANUSCRIPT
# M. Kirkland
# Created 02/02/24
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(geojsonsf)
library(ggnewscale) ## For new_scale*()
library(terra)
library(ggforce) ## For geom_mark_hull()
library(lme4)
library(geosphere)
library(ecodist)
library(Hmisc) ## For capitalize()
library(plyr) ## For rbindfill()
library(emmeans)
library(lmerTest) ## To get p-value from lme model 
library(FSA) ## For dunn_test()
library(concaveman)

# Clear working directory
rm(list = ls())

# Read in wintering ground locations for UK birds
setwd(".Github/geolocator data/results/wgr")
filenames <- list.files(getwd(), pattern="*.csv", full.names=TRUE)
uk = tibble(file = filenames) %>%
  tidyr::extract(file, "ID", "(?<=wgr/)(.+)(?=\\_wgr.csv)", remove = FALSE) %>% # Add column with bird tag ID
  mutate(data = lapply(file, read.csv)) %>%
  unnest(data) %>%
  dplyr::select(-file)

countries <- st_read(dsn = "~/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp") # Read in country uk 

# Add breeding lat and long
uk$blat <- ifelse(uk$tag.loc == "Methwold Hythe", 52.5240,
                    ifelse(uk$tag.loc == "Alton Water", 51.9804,
                           ifelse(uk$tag.loc == "Orlestone Forest", 51.0805556,
                                  52.2981)))
uk$blon <- ifelse(uk$tag.loc == "Methwold Hythe", 0.5234,
                    ifelse(uk$tag.loc == "Alton Water", 1.1326,
                           ifelse(uk$tag.loc == "Orlestone Forest", 0.8322222222222,
                                  -0.31504)))


# Ring numbers 
setwd("./geolocator data")
uk.birds <- read.csv("Tagging sites for retrieved data.csv")
uk <- merge(uk, uk.birds, by.x = "ID", by.y = "Tag_no", all = T)
uk$Ring.No[uk$ID == "0AD"] <- "0AD"

# Select all wintering grounds of UK tagged birds estimated using a land mask an
uk.wgr <- uk %>% 
  filter(site == "wgr" & mask == "land") %>% 
  # Some birds have to years of data so split by year
  group_by(ID) %>% 
  # Remove stationary periods after 15th March
  filter(!(Arrival > as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)) + 1,"-03-15"), tz = "GMT")) & 
           # Remove stationary periods before 15th October
          !(Departure < as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)),"-10-15"), tz = "GMT"))) %>% 
  mutate(year = min(as.numeric(substr(Arrival,1,4)))) %>% 
  group_by(ID, year) %>% 
  # Number sites for each bird
  dplyr::mutate(site.num = row_number(),
  # Calculate number of sites
         nsite = length(site.num),
  # Identify mid winter location, using period defined by Hahn et al. 
         overlap.mid = ifelse(Arrival < paste0(substr(Arrival, 1,4),"/12/31") & Departure >= paste0(substr(Arrival, 1,4),"/12/31"), "Y", "N")) %>% 
  dplyr::select(ID, Arrival, Departure, Lon, Lat, Lon.2.5., Lon.97.25., Lat.2.5., Lat.97.25.,blon, blat, tag.loc, site.num, nsite, overlap.mid, Ring.No)

# Look at repeats? 
dups <- unique(uk.wgr$Ring.No[duplicated(uk.wgr[,c("Ring.No","site.num")])])
uk.wgr$site.num[uk.wgr$Ring.No %in% dups] ## Repeated birds all used one winter site 

# Average locations for repeated inds 
uk.avg <- uk.wgr %>% 
  group_by(Ring.No, site.num) %>% 
  dplyr::mutate(Lon = mean(Lon),
                Lat = mean(Lat)) %>% 
  dplyr::distinct(Ring.No, .keep_all=TRUE)

# Read in breeding grounds of birds tagged in Africa
setwd("./geolocator data/results/brgr")
filenames <- list.files(getwd(), pattern="*.csv", full.names=TRUE)
african = tibble(file = filenames) %>%
  tidyr::extract(file, "ID", "(?<=brgr/)(.+)(?=\\_brgr.csv)", remove = FALSE) %>% # Add column with bird tag ID
  mutate(data = lapply(file, read.csv)) %>%
  unnest(data) %>%
  dplyr::select(-file)
# Average repeats 
african.avg <- african
african.avg$indID <- african.avg$ID
african.avg$indID[african.avg$ID == "BV422"] <- african.avg$indID[african.avg$ID == "BR064"]
african.avg <- african.avg %>% 
  group_by(indID, mask, site) %>% 
  dplyr::mutate(Lon = mean(Lon),
                Lat = mean(Lat)) %>% 
  dplyr::distinct(indID, .keep_all=TRUE)

  # Read in ringiafrican.avg# Read in ringi# Read in ringi# Read in ringing data from Kartong Bird Observatory
setwd("./ringing data")
ringdata <- read.csv("KBO_DATA_NIGAL_20231012.csv")

# Select birds recaptured elsewhere
recap <- unique(ringdata[,c('RING', 'PLACE')])
recap <- recap$RING[duplicated(recap$RING)]
recap <- ringdata[ringdata$RING %in% recap,]
recap <- recap[!duplicated(recap[,c('RING', 'PLACE')]),]
recap <- recap %>% 
  group_by(RING) %>% 
  dplyr::mutate(site = row_number()) # Add site number to differentiate between tagging location and recapture site

# Read in GPS tags
setwd("./gps data")
# Breeding locations
gps1 <- st_read("Obs150324_104804_Tag58279.kml")
gps2 <- st_read("Obs130324_154637_Tag58288.kml")
# Wintering locations 
gps3 <- st_read("Obs140524_094210_Tag58392.kml")
gps4 <- st_read("Obs140524_094834_Tag58400.kml")
plot(gps3$geometry)
plot(gps4$geometry)

# Extract coorindates. Last row is the path taken so ignore. 
gps1_coords <- as.data.frame(gps1 %>% st_cast("MULTIPOINT") %>% st_coordinates())
gps2_coords <- as.data.frame(gps2 %>% st_cast("MULTIPOINT") %>% st_coordinates())
gps3_coords <- as.data.frame(gps3 %>% st_cast("MULTIPOINT") %>% st_coordinates())
gps4_coords <- as.data.frame(gps4 %>% st_cast("MULTIPOINT") %>% st_coordinates())

# Read in results of SDM
setwd("./sdm")
sdm <- rast("results-rasters.tif") 
# Crop by extent of records
# Read in dataset of from GEE
files <- list.files(path = "./sdm/GEE results", pattern='.csv', 
                    all.files = T, full.names = T)
recs <- lapply(files, read.csv)
recs_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = F), recs)
# Select records used in SDM
recs_df <- recs_df[complete.cases(recs_df),]
recs_sf <- cbind(st_as_sf(geojson_sf(recs_df$.geo)), recs_df)
# Remove Eastern section from plot
recs_sf <- st_crop(recs_sf, ext(c(-18.0416195374957, 15, 2.99999216327873, 17.5)))
# Create polygon around points
recs_poly <- concaveman(recs_sf, concavity = 1, length_threshold = 5)
sdm_cr <- sdm %>% crop(recs_poly) %>% terra::mask(recs_poly)
plot(sdm_cr)
plot(countries, add = T, col = NA)
plot(recs_sf, add = T, col = "black")

# Data from Hahn et al. https://link.springer.com/article/10.1007/s00442-013-2726-4
setwd(".")
nigal.hahn <- read.csv("LusMeg_F_IT_BG_KernelMedians for Chris_Nov2023.csv") ## Available upon request from authors
# Select wintering sites. Populations are France/Western (1), Italy/Central (2) and Bulgaria/Eastern (3)
# Type = 0: breeding; 1: autumn staging; 2: pre winter site; 3: main winter; 4) spring staging site
nigal.hahn2 <- nigal.hahn[nigal.hahn$type == 3,]
# Make sure columns match Finch data
colnames(nigal.hahn2)[c(3,6,7)] <- c("popunique", "wlon1","wlat1")
nigal.hahn2 <- nigal.hahn2[,c(3,6,7)]
# Add additional info
nigal.hahn2$species <- "common nightingale"
nigal.hahn2$system <- "afro"
nigal.hahn2 <- nigal.hahn2 %>% mutate(country = ifelse(popunique == 1, "France", ifelse(popunique == 2, "Italy", "Bulgaria")))

# PLOTS ~~~~~~~~~~~~~~~~~~~~

# Fig. 2: ####
# Comparison with Finch et al. ~~~~~~~~~~~~~~~~
## Obtained from https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.12635
# Do this first to get coefficient from model for spread corrected for sample size

# Convert dataset to wide format to use only first wintering location
wgr.wide <- pivot_wider(data = uk.avg, 
                        id_cols = c(ID, blon, blat, year), 
                        names_from = site.num, 
                        values_from = c("Lon", "Lat"))

# Rename to match column names from Finch et al. 2017
colnames(wgr.wide)[5:8] <- c("wlon1","wlon2","wlat1","wlat2")

# Custom function to calculate mean inter-individual distance
meanDist <- function(lon, lat){
  mean(lower(distm(cbind(lon, lat), fun = distHaversine)))
}

sdDist <- function(lon, lat){
  sd(lower(distm(cbind(lon, lat), fun = distHaversine)))
}

# Read in data from Finch study
setwd("~/BTO projects/nightingale migration/doi_10.5061_dryad.ss3r7__v1")
other.data <- read.csv("data.csv")

# Replace nightingale data with more accurate data from Hahn et al. 
other.data <- other.data %>% 
  filter(!study == "Hahn et al 2014 Ecol & Evol")

# Add additional UK info to common nightingale data
wgr.wide$species <- "common nightingale"
wgr.wide$popunique <- 4
wgr.wide$country <- "england"
wgr.wide$system <- "afro"

# Remove year as a variable
wgr.wide <- subset(wgr.wide, select = -year)

# Merge datasets
all.data <- bind_rows(other.data, wgr.wide, nigal.hahn2)

# Fix thrush nightingale - set main winter site as final one
all.data[all.data$species == "thrush nightingale",]$wlat1 <- all.data[all.data$species == "thrush nightingale",]$wlat3
all.data[all.data$species == "thrush nightingale",]$wlon1 <- all.data[all.data$species == "thrush nightingale",]$wlon3

# Calculate spread of each population, using land mask
data.pop <- all.data %>% 
  # Create ID for each population of each species
  dplyr::mutate(popid = paste(species, system, popunique, sep = ".")) %>% 
  dplyr::group_by(popid) %>% 
  dplyr::mutate(nind = length(popid),
                # Set unknown age to -1
                age = ifelse(is.na(age), -1, age)) %>% 
  # Exclude juveniles
  filter(age != 0 &
           # Exclude if outside winter range (afro and neo)
           ((system == "afro" & wlat1 <= 20 & wlon1 <= 65 & blon >= -20 & blon <= 65) | (system == "neo" & wlat1 <= 30) | species == "common nightingale") & ## No breeding data for common nightingale so keep all
           # Exclude if tagged overwinter
           (sampledat == 1 | is.na(sampledat)) &
           # Exclude single individual populations and populations with no ID
           nind > 1 & !is.na(popunique) &
           # Exclude northern wheatear from neo (winters in afro)
           !(system == "neo" & species == "northern wheatear")) %>% 
  dplyr::summarise('system' = first(system),
                   'species' = first(species),
                   'age' = first(age),
                   'nind' = first(nind),
                   'ln' = log(nind),
                   'n2' = length(nind) ^ 2,
                   'country' = first(country),
                   # Mean breeding longitude & lat
                   'b.lon' = mean(blon),
                   'b.lat' = mean(blat),  
                   # Longitudinal range of breeding sites 
                   'b.londif' = diff(range(blon)), 
                   # Max distance between breeding sites
                   'b.maxdist' = (max(lower(distm(cbind(blon, blat), fun = distHaversine))))/1000,
                   # Mean distance between breeding sites
                   'b.meandist' = (meanDist(blon, blat))/1000, 
                   # Mean winter longitude & lat
                   'w.lon' = mean(wlon1), 
                   'w.lat' = mean(wlat1),
                   # Longitudinal range of winter sites
                   'w.londif' = diff(range(wlon1)), 
                   # Max distance between winter sites
                   'w.maxdist' = (max(lower(distm(cbind(wlon1, wlat1), fun = distHaversine))))/1000, # First wintering site is always longest 
                   # Mean distance between winter sites
                   'w.meandist' = (meanDist(wlon1, wlat1))/1000,
                   # Standard deviation around the mean
                   'w.sddist' = (sdDist(wlon1, wlat1))/1000)

# Control for population size
# Model population spread by population size
mod1 <- lmer(w.meandist ~ ln + (1|species), data = data.pop, REML = F)
mod2 <- lmer(w.meandist ~ nind + (1|species), data = data.pop, REML = F)
mod3 <- lmer(w.meandist ~ n2 + nind + (1|species), data = data.pop, REML = F)
AIC(mod1,mod2,mod3)

# Check cook's distance and re-run model without outliers 
cooks <- cooks.distance(influence(mod1, obs = TRUE))
mod1b <- lmer(w.meandist ~ ln + (1|species), subset(data.pop, cooks < quantile(cooks, 0.95)))
fit <- data.frame(n = seq(1, 51, 1), ln = log(seq(1, 51, 1)))
fit <- cbind(fit, AICcmodavg::predictSE(mod = mod1b, newdata = fit))

# Predict spread assuming equal sample size 
summary(mod1)
data.pop$pred <- data.pop$w.meandist + 137.42*(max(fit$ln)-data.pop$ln) ## Uses coefficient from model output

# Plot population spread by species
# Lock in factor orders
data.pop <- data.pop %>% arrange(desc(pred)) %>% group_by(species) %>% dplyr::mutate(popunique = row_number())  # Order by increasing spread
data.pop$species <- factor(capitalize(data.pop$species), levels = rev(unique(capitalize(data.pop$species))))

ggplot(data.pop, aes(x = species, y = pred, fill = as.factor(popunique))) + 
  theme_classic() + 
  geom_col(position="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  xlab("") + ylab("Population spread (km)") +
  scale_fill_manual(values = c('#ece7f9','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858','black'))

## Is the difference between nightingale populations significant? 
library("ggpubr")
distsFun <- function(lon, lat){
  lower(distm(cbind(lon, lat), fun = distHaversine))
}

# First explore distances in UK nightingales and wood thrush
wt.dists <- all.data %>% 
  dplyr::filter((species == "common nightingale" & popunique == 4) | (species == "wood thrush" & popunique == 5)) %>% 
  dplyr::group_by(species) %>% 
  dplyr::reframe(dists = distsFun(wlon1, wlat1)/1000)

wilcox.test(dists ~ species, data = wt.dists)

# Now repeat for nightingales, but using mid-winter period
mid.sites <- uk.avg[uk.avg$overlap.mid == "Y",]
# Are there any birds that don't have a mid-winter period?
mid.sites <- rbind(mid.sites,uk.avg[!uk.avg$Ring.No %in% mid.sites$Ring.No & uk.avg$site.num == 1,])
colnames(mid.sites)[5:6] <- c("wlon1","wlat1")
nigal.dists <- bind_rows(mid.sites, nigal.hahn2) %>% 
  mutate(popunique = ifelse(is.na(popunique), 4, popunique),
         n = n()) %>% 
  dplyr::group_by(popunique) %>% 
  dplyr::reframe(dists = distsFun(wlon1, wlat1)/1000)

# Summarise mean and 85% CIs for distances within nigal pops 
nigal.dists.sum <- bind_rows(mid.sites, nigal.hahn2) %>% 
  mutate(popunique = ifelse(is.na(popunique), 4, popunique)) %>% 
  dplyr::group_by(popunique) %>% 
  dplyr::summarize(n = n(), 
                   mean = (meanDist(wlon1, wlat1))/1000,
                   up.ci = mean + ((sdDist(wlon1, wlon1))/1000)/sqrt(n),
                   low.ci = mean - ((sdDist(wlon1, wlon1))/1000)/sqrt(n))

# Fig. 3: ####
ggplot(nigal.dists, aes(x = as.factor(popunique))) +
  geom_violin(mapping = aes(group = as.factor(popunique), y = dists), fill = "grey90") +
  ylab("Inter-individual distances (km)") + xlab("Tagging location") +
  scale_x_discrete(labels = c("France", "Italy", "Bulgaria", "UK")) +
  geom_point(nigal.dists.sum, mapping = aes(y = mean)) +
  geom_errorbar(nigal.dists.sum, mapping = aes(ymax = up.ci, ymin = low.ci, width = 0)) + 
  theme_classic() +
  theme(panel.grid.major.y = element_line())

kruskal.test(dists ~ popunique, data = nigal.dists)

dunnTest(dists ~ as.factor(popunique), data = nigal.dists,
         method = "bh") 

# Does this change if we use the first wintering site? 
nigal.dists2 <- all.data %>% 
  filter(species == "common nightingale") %>% 
  dplyr::group_by(popunique) %>% 
  dplyr::reframe(dists = distsFun(wlon1, wlat1)/1000)

kruskal.test(dists ~ popunique, data = nigal.dists2)

# Remove outlier in Bulgarian population
nigal.dists.out <- bind_rows(mid.sites, nigal.hahn2[nigal.hahn2$wlon1 < 30,]) %>% 
  mutate(popunique = ifelse(is.na(popunique), 4, popunique)) %>% 
  dplyr::group_by(popunique) %>% 
  dplyr::reframe(dists = distsFun(wlon1, wlat1)/1000)

kruskal.test(dists ~ popunique, data = nigal.dists.out)

# Figure 1a: ####
# Comparison of estimated wintering grounds ~~~~~~~~~~~~~~~~
ggplot() +
  geom_sf(data = countries, fill = "grey90", col = "grey50", linewidth = .1) +
  geom_pointrange(data = uk.wgr[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y",], mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .4) +
  geom_errorbarh(data = uk.wgr[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y",], mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .4) +
  geom_pointrange(data = african[african$mask == "land" & african$site == "wgr" & african$tag.loc == "Gambia",], mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .4) +
  geom_errorbarh(data = african[african$mask == "land" & african$site == "wgr" & african$tag.loc == "Gambia",], mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .4) +
  geom_point(data = uk.wgr[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y",], mapping = aes(x = Lon, y = Lat, fill = as.factor(site.num)), pch = 21, size = 4, stroke = .2) + 
  geom_point(data = african[african$mask == "land" & african$site == "wgr" & african$tag.loc == "Gambia",], mapping = aes(x = Lon, y = Lat), pch = 21, size = 4, fill = "red", stroke = .2) +
  scale_fill_manual(values = c("blue", "lightblue")) +
  geom_point(mapping = aes(x = -16.7666672, y = 13.1), size = 2, fill = "yellow", pch = 24, stroke = .2) +
  # Add symbol to show repeats. See tagging sites info for ring numbers.
  geom_point(data = uk.wgr[uk.wgr$ID == "097",], mapping = aes(x = Lon, y = Lat), size = 1, col = "white") +
  geom_point(data = uk.wgr[uk.wgr$ID == "S418",], mapping = aes(x = Lon+.2, y = Lat+.2), size = 1, col = "white") +
  geom_point(data = uk.wgr[uk.wgr$ID == "BE482",], mapping = aes(x = Lon, y = Lat), size = 1, col = "white") + 
  geom_point(data = uk.wgr[uk.wgr$ID == "075",], mapping = aes(x = Lon, y = Lat), size = 1, col = "#fed976") +
  geom_point(data = uk.wgr[uk.wgr$ID == "S412",], mapping = aes(x = Lon, y = Lat), size = 1, col = "#fed976") +
  geom_point(data = uk.wgr[uk.wgr$ID == "098",], mapping = aes(x = Lon, y = Lat), size = 1, col = "pink") +
  geom_point(data = uk.wgr[uk.wgr$ID == "S420",], mapping = aes(x = Lon-.2, y = Lat-.2), size = 1, col = "pink") +
  geom_point(data = uk.wgr[uk.wgr$ID == "BE491",], mapping = aes(x = Lon, y = Lat), size = 1, col = "#addd8e") + 
  geom_point(data = uk.wgr[uk.wgr$ID == "Z411",], mapping = aes(x = Lon, y = Lat), size = 1, col = "#addd8e") +
  geom_point(data = african[african$mask == "land" & african$site == "wgr" & african$ID == "BV422",], mapping = aes(x = Lon, y = Lat), size = 1, col = "white") +
  geom_point(data = african[african$mask == "land" & african$site == "wgr" & african$ID == "BR064",], mapping = aes(x = Lon, y = Lat), size = 1, col = "white") +
  xlim(-19,-3) + 
  ylim(-4,37) + 
  theme_classic() + 
  guides(fill = "none") +
  xlab("") + ylab("") +  
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank())
## 360x756

# Figure 1b: ####
# Results of SDM ~~~~~~~~~~~~~~~~

# Crop to wintering grounds of UK birds
sdm_uk <- sdm_cr %>% crop(ext(c(-18.0416195374957, -10, 6, 17.5)))
sdm_uk_df <- as.data.frame(sdm_uk, xy = T)

ggplot() + 
  geom_sf(data = countries, fill = "white", col = "grey66", linewidth = .1) +
  geom_tile(sdm_uk_df, mapping = aes(x = x, y = y, fill = sum)) +
  geom_sf(data = countries, fill = NA, col = "grey66", linewidth = .1) +
  geom_point(mapping = aes(x = -16.7666672, y = 13.1), size = 4, fill = "yellow", shape = 24, stroke = .2) +
  geom_point(data = uk.wgr[uk.wgr$overlap.mid == "Y" | uk.wgr$site.num == 1,], mapping = aes(x = Lon, y = Lat), shape = 21, size = 4, fill = NA, stroke = .3) +
  geom_point(gps3_coords[which.min(gps3_coords$Y),], mapping = aes(x = X, y = Y), size = 3, shape = 4, stroke = 1) +
  geom_point(gps4_coords[which.min(gps4_coords$Y),], mapping = aes(x = X, y = Y), size = 3, shape = 4, stroke = 1) +
  scale_y_continuous(expand = c(0,0), limits = c(5,17.5)) +
  scale_x_continuous(expand = c(0,0), limits = c(-18,-10), breaks = seq(-18,-10, by = 2)) +
  scale_fill_gradientn(colors = rev(c("#006837", "#00A600", "#1DB000", "#3EBB00", "#63C600", "#8BD000", "#B6DB00", "#E6E600", "#E7CE1D", "#E9BD3A", "#EAB358", "white"))) +
  theme_classic() + 
  guides(shape = "none") +
  xlab("") + ylab("") +
  guides(fill = "none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank()) 
## 300x410

# Figure 1c: ####
## A plot to look at differences between known population spread and estimated population spread from GLS at wintering grounds

uk.wgr.all <- uk %>% 
  filter(site == "wgr") %>% 
  group_by(ID) %>% 
  # Remove stationary periods after 15th March
  filter(!(Arrival > as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)) + 1,"-03-15"), tz = "GMT")) & 
           # Remove stationary periods before 15th October
           !(Departure < as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)),"-10-15"), tz = "GMT"))) %>% 
  mutate(year = min(as.numeric(substr(Arrival,1,4)))) %>% 
  group_by(ID, year, mask) %>% 
  # Number sites for each bird
  dplyr::mutate(site.num = row_number(),
                # Calculate number of sites
                nsite = length(site.num),
                # Identify mid winter location, using period defined by Hahn et al. 
                overlap.mid = ifelse(Arrival < paste0(substr(Arrival, 1,4),"/12/31") & Departure >= paste0(substr(Arrival, 1,4),"/12/31"), "Y", "N")) %>% 
  ungroup() %>% 
  filter(site.num == 1) %>% 
  # Average individuals
  group_by(Ring.No, mask) %>% 
  dplyr::mutate(Lon = mean(Lon),
                Lat = mean(Lat)) %>% 
  dplyr::distinct(Ring.No, .keep_all=TRUE) %>% 
  dplyr::select(ID, Arrival, Departure, Lon, Lat, Lon.2.5., Lon.97.25., Lat.2.5., Lat.97.25.,blon, blat, tag.loc, site.num, nsite, overlap.mid, mask)

gls.error.wgr <- rbind.fill(uk.wgr.all, african.avg[african.avg$site == "wgr",]) %>% 
  mutate(group = ifelse(tag.loc == "Ghana", "Ghana", ifelse(tag.loc == "Gambia", "Gambia", "UK"))) %>% 
  group_by(ID, year, mask) %>% 
  mutate(site.num = row_number()) %>% 
  ungroup() %>% 
  filter(site.num == 1) %>% 
  group_by(group, mask) %>% 
  dplyr::summarise(
    # Sample size
    'n' = n(),
    # Mean distance between winter sites
    'w.meandist' = (meanDist(Lon, Lat))/1000,
    # Standard deviation around the mean
    'w.sddist' = (sdDist(Lon, Lat))/1000) %>% 
  ## Repeated inds so n is 26 for UK birds
  mutate(
    # Logged 
    'ln' = log(n),
    # Corrected for sample size
    'corr' = w.meandist+137.98*(max(fit$ln)-ln),
    # 84% CIs
    'w.ci.up' = corr + 1.405072*(w.sddist/sqrt(n)),
    'w.ci.lw' = corr - 1.405072*(w.sddist/sqrt(n)))

ggplot(gls.error.wgr, aes(x = factor(group, levels = c("UK", "Gambia", "Ghana")), y = corr, col = mask)) +
  geom_point(size = 2, position = position_dodge(width=0.3)) +
  geom_errorbar(mapping = aes(ymax = w.ci.up, ymin = w.ci.lw), width = .2, position = position_dodge(width=0.3)) +
  scale_shape_manual(values = c(16, 16, 4)) +
  scale_colour_manual(values = c("black", "grey70")) +
  ylab("Population spread (km)") + xlab("Tagging location") +
  ylim(-1, 2000) +
  theme_classic() +
  guides(col = "none") +
  theme(axis.title.y = element_text(size = 12,
                           margin = unit(c(0, 3, 0, 0), "mm")),
        axis.title.x = element_text(size = 12,
                                    margin = unit(c(3, 0, 0, 0), "mm"))) +
  geom_point(mapping = aes(x = "Gambia", y = 0), shape = 4, size = 2, col = "black") +
  geom_point(mapping = aes(x = "Ghana", y = 0), shape = 4, size = 2, col = "black")

# Do some stats to see if differences in pop spread at the wintering grounds is significantly different
wgr.first <- rbind.fill(african.avg[african.avg$site == "wgr" & african.avg$mask == "land",],
                        uk.avg[uk.avg$site.num == 1,])
wgr.first$group <- ifelse(wgr.first$tag.loc == "Ghana", "Ghana", ifelse(wgr.first$tag.loc == "Gambia", "Gambia", "England"))

wgr.first.dist <- wgr.first %>% 
  dplyr::group_by(group) %>% 
  dplyr::reframe(dists = distsFun(Lon, Lat)/1000)

kruskal.test(dists ~ group, data = wgr.first.dist)

wgr.mid <- rbind.fill(african.avg[african.avg$site == "wgr" & african.avg$mask == "land",],mid.sites)
wgr.mid$group <- ifelse(wgr.mid$tag.loc == "Ghana", "Ghana", ifelse(wgr.mid$tag.loc == "Gambia", "Gambia", "England"))

wgr.mid.dist <- wgr.mid %>% 
  dplyr::group_by(group) %>% 
  dplyr::reframe(dists = distsFun(Lon, Lat)/1000)

kruskal.test(dists ~ group, data = wgr.mid.dist)

# Fig. 4a: ####
# Comparison of estimated breeding grounds ~~~~~~~~~~~~~~~~
ggplot() + 
  geom_sf(data = countries, fill = "grey90", col = "grey66", linewidth = .1) +
  geom_pointrange(data = uk[uk$mask == "land" & uk$site == "brgr",], mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .6) +
  geom_errorbarh(data = uk[uk$mask == "land" & uk$site == "brgr",], mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .6) +
  geom_pointrange(data = african[african$mask == "land" & african$site == "brgr" & african$tag.loc == "Gambia",], mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .6) +
  geom_errorbarh(data = african[african$mask == "land" & african$site == "brgr" & african$tag.loc == "Gambia",], mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .6) +
  geom_point(data = african[african$mask == "land" & african$site == "brgr" & african$tag.loc == "Gambia",], mapping = aes(x = Lon, y = Lat), size = 3, pch = 21, fill = "red", stroke = .2) +
  geom_point(data = uk, mapping = aes(x = blon, y = blat), size = 3, fill = "yellow", col = "black", pch = 24, stroke = .2) + ## Plot actual breeding locations
  geom_point(data = uk[uk$mask == "land" & uk$site == "brgr",], mapping = aes(x = Lon, y = Lat), size = 3, pch = 21, fill = "blue") +
  # Indicate estimation of second breeding site of same bird
  geom_point(data = african[african$mask == "land" & african$site == "brgr" & african$ID == "BV422" | african$mask == "land" & african$site == "brgr" & african$ID == "BR064",], mapping = aes(x = Lon, y = Lat), size = 1, col = "white") +
  # Add GPS data
  geom_point(gps1_coords[which.max(gps1_coords$Y),], mapping = aes(x = X, y = Y), shape = 22, size = 3, fill = "red") +
  geom_point(gps2_coords[which.max(gps2_coords$Y),], mapping = aes(x = X, y = Y), shape = 22, size = 3, fill = "red") +
  xlim(-8,9) + 
  ylim(44,59) + 
  theme_classic() + 
  xlab("") + ylab("") +
  guides(fill = "none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank()) 
## 420x370

# Is breading ground significantly different for the three pops? 
brgr <- rbind.fill(african.avg[african.avg$site == "brgr" & african.avg$mask == "land",],
                   uk[uk$site == "brgr" & uk$mask == "land",])
brgr$group <- ifelse(brgr$tag.loc == "Ghana", "Ghana", ifelse(brgr$tag.loc == "Gambia", "Gambia", "England"))

brgr.dists <- brgr %>% 
  dplyr::group_by(group) %>% 
  dplyr::reframe(dists = distsFun(Lon, Lat)/1000)

kruskal.test(dists ~ group, data = brgr.dists)
dunnTest(dists ~ as.factor(group), data = brgr.dists,
         method = "bh") 

# Fig. 4b: ####
# Ringing recoveries ~~~~~~~~~~~~~~~~

# Extract national ringing data between 2001 & 2023
setwd("~/BTO projects/nightingale migration/ringing data")
natring <- read.csv("nigal ringing totals.csv")
natring <- natring %>% 
  mutate(species = gsub(" ","",species)) %>% 
  filter(YEAR > 2000,
         species == "Lusciniamegarhynchos") %>% 
  group_by(scheme) %>% 
  dplyr::summarise(total = sum(as.numeric(TOTAL),na.rm = T)) ## Total is in character format so convert to numeric
countries <- merge(countries, natring, by.x = "SOV_A3", by.y = "scheme", all.x = T)

ggplot() +
  geom_sf(data = countries, fill = "grey90", col = "grey66", linewidth = .1) +
  geom_sf(data = countries[!is.na(countries$total),], mapping = aes(fill = total), col = "grey66", linewidth = .1) +
  scale_fill_gradientn(colours = c('#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32')) +
  guides(fill = "none") +
  geom_point(data = recap, mapping = aes(x = LONG, y = LAT), size = 2, shape = 22) +
  geom_point(mapping = aes(y = 51.8253, x = 0.8646), size = 2, shape = 22) +
  xlim(-10,15) + 
  ylim(37,60) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = .5),
        legend.position = "right",
        axis.line = element_blank(),) +
  xlab("") + ylab("") 
## 300x370

# Fig. 4c: ####
## A plot to look at differences between known population spread and estimated population spread from GLS at breeding grounds
gls.error.brgr <- rbind.fill(uk[uk$site == "brgr",], african.avg[african.avg$site == "brgr",]) %>% 
  mutate(group = ifelse(tag.loc == "Ghana", "Ghana", ifelse(tag.loc == "Gambia", "Gambia", "UK"))) %>% 
  group_by(group, mask) %>% 
  dplyr::summarise(
    # Sample size
    'n' = n(),
    # Mean distance between winter sites
    'w.meandist' = (meanDist(Lon, Lat))/1000,
    # Standard deviation around the mean
    'w.sddist' = (sdDist(Lon, Lat))/1000) %>% 
  ## Repeated inds so n is 4 for Gambia birds
  mutate(
    # Logged 
    'ln' = log(n),
    # Corrected for sample size
    'corr' = w.meandist+137.98*(max(fit$ln)-ln),
    # CIs
    'w.ci.up' = corr +  1.405072*(w.sddist/sqrt(n)),
    'w.ci.lw' = corr -  1.405072*(w.sddist/sqrt(n)))

ggplot(gls.error.brgr, aes(x = factor(group, levels = c("UK", "Gambia", "Ghana")), y = corr, col = mask)) +
  geom_point(size = 2, position = position_dodge(width=0.3)) +
  geom_errorbar(mapping = aes(ymax = w.ci.up, ymin = w.ci.lw), width = .2, position = position_dodge(width=0.3)) +
  scale_shape_manual(values = c(16, 16, 4)) +
  scale_colour_manual(values = c("black", "grey70")) +
  ylab("Population spread (km)") + xlab("Tagging location") +
  ylim(-1, 2100) +
  theme_classic() +
  guides(col = "none") +
  theme(axis.title.y = element_text(size = 12,
                                    margin = unit(c(0, 3, 0, 0), "mm")),
        axis.title.x = element_text(size = 12,
                                    margin = unit(c(3, 0, 0, 0), "mm"))) +
  geom_point(mapping = aes(x = "UK", y = 91.2), shape = 4, size = 2, col = "black") 

# Fig. 4d: ####
# Polygons of breeding sites of birds tagged in Gambia and Ghana ~~~~~~~~~~~~~~~~

# Tagging locations from Hahn et al. in text
hahn.tag <- data.frame(blon = c(7.5, 11.8, 27.9, 28.3), blat = c(47.6, 44.6, 42.1, 43.4), country = c( "France", "Italy", "Bulgaria", "Bulgaria"))

# Create minimum convex hulls for each breeding pop.
hull.brgr <- st_as_sf(data.frame(
    x = c(african.avg$Lon[african.avg$mask == "land" &
                        african.avg$site == "brgr" & african.avg$tag.loc == "Ghana"],
          hahn.tag$blon[hahn.tag$country != "Bulgaria"],
          african.avg$Lon[african.avg$mask == "land" &
                        african.avg$site == "brgr" & african.avg$tag.loc == "Gambia"],
          recap$LONG[recap$LAT > 30]),
    y = c(african.avg$Lat[african.avg$mask == "land" &
                        african.avg$site == "brgr" & african.avg$tag.loc == "Ghana"],
          hahn.tag$blat[hahn.tag$country != "Bulgaria"],
          african.avg$Lat[african.avg$mask == "land" &
                        african.avg$site == "brgr" & african.avg$tag.loc == "Gambia"],
          recap$LAT[recap$LAT > 30]),
    country = c(rep("Ghana", 5), rep("Gambia", 10))
  ), coords = c("x", "y"), crs = 4326) %>% 
  group_by(country) %>% 
  dplyr::summarise( geometry = st_combine( geometry ) ) %>%
  sf::st_convex_hull()

hull.wgr <- st_as_sf(data.frame(
  x = c(uk.wgr$Lon[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y"], nigal.hahn2$wlon1[nigal.hahn2$country != "Bulgaria"]),
  y = c(uk.wgr$Lat[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y"], nigal.hahn2$wlat1[nigal.hahn2$country != "Bulgaria"]),
  country = c(rep("Gambia", length(uk.wgr$Lon[uk.wgr$site.num == 1 | uk.wgr$overlap.mid == "Y"])), rep("Ghana", length(nigal.hahn2$wlon1[nigal.hahn2$country != "Bulgaria"]))) 
), coords = c("x", "y"), crs = 4326) %>% 
  group_by(country) %>% 
  dplyr::summarise( geometry = st_combine( geometry ) ) %>%
  sf::st_convex_hull()

# Create minimum convex polygoins
ggplot() + 
  xlim(-18,18) + 
  ylim(2,53) + 
  theme_classic() + 
  xlab("") + ylab("") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank(), 
        legend.position = "none") +
  geom_sf(data = countries, fill = "grey90", col = "grey50", linewidth = .1) +
  # Connection between UK and Gambia
  geom_point(data = african.avg[african.avg$mask == "land" & african.avg$site == "brgr" & african.avg$tag.loc == "Gambia",], mapping = aes(x = Lon, y = Lat), col = "#542788", alpha = .6) +
  geom_point(mapping = aes(x = -16.7666672, y = 13.1), fill = "#542788", pch = 24, col = "#542788") +
  # From ringing data
  geom_point(data = recap[recap$LAT > 30,], mapping = aes(x = LONG, y = LAT), fill = "#542788", col = "#542788", pch = 22, alpha = .6) +
  # Connection between birds tagged in Ghana
  geom_point(data = african.avg[african.avg$mask == "land" & african.avg$site == "brgr" & african.avg$tag.loc == "Ghana",], mapping = aes(x = Lon, y = Lat), col = "#e08214") +
  geom_point(mapping = aes(x = -2.464, y = 7.402), fill = "#e08214", pch = 24, col = "#e08214") +
  # Connection between Italian and French breeding and wintering populations
  geom_point(data = hahn.tag[hahn.tag$country != "Bulgaria",], mapping = aes(x = blon, y = blat), col = "#e08214", fill = "#e08214", shape = 24, size = .9) +
  # Draw polygons around breeding site
  geom_sf(hull.brgr, mapping = aes(fill = country, col = country), alpha = .3) +
  scale_fill_manual(values = c("#542788", "#e08214")) + 
  scale_colour_manual(values = c("#542788", "#e08214"))
## 340x480

# Intersect with land, convert to European equal area reference system and calculate size
hulls.land <- st_intersection(hull.brgr, st_union(countries))
st_area(st_transform(hulls.land, 3035))/1e6

# Supp figures: ####

## Results of SDM ~~~~~~~~~~~~~~~~

## Overall suitability ####

# Convert 'sdm' to a data frame
sdm_df <- as.data.frame(sdm_cr, xy = TRUE)

# Identify points overlapping with regions above the MSS threshold of suitability
nigal.hahn2$ID <- rownames(nigal.hahn2)
all_points <- nigal.hahn2
colnames(all_points)[2:3] <- c("Lon", "Lat")

# Select mid-winter sites
uk_midwinter_pts <- uk.wgr[uk.wgr$overlap.mid == "Y",]

# Include birds without a mid-winter period
uk_midwinter_pts <-
  rbind(uk_midwinter_pts, uk.wgr[!uk.wgr$Ring.No %in% uk_midwinter_pts$Ring.No &
                                   uk.wgr$site.num == 1, ])
uk_midwinter_pts$ID2 <- rownames(uk_midwinter_pts)

all_points_sf <- st_as_sf(rbind.fill(all_points, uk_midwinter_pts), coords = c("Lon", "Lat"), crs = 4326)

# Add ring number to identify repeated individuals
all_points_sf$Ring.No <- ifelse(is.na(all_points_sf$Ring.No), all_points_sf$ID, all_points_sf$Ring.No)
# Add buffer
all_points_sf <- all_points_sf %>%
  st_transform(crs = 32628) %>% # Transform to UTM for Gambia
  st_buffer(100000) %>% ## 100 km buffer
  st_transform(crs = 4326) # Transform back  

# Extract suitability for each individual   
min_suit <- terra::extract(sdm, all_points_sf, fun = max, na.rm = T, bind = T)

# Suitable points
suitable_points <- min_suit[min_suit$sum >= 0.03,]

# Plotting suitability maps
sf_use_s2(FALSE)
p1 <- ggplot() + 
  geom_sf(data = st_intersection(countries,recs_poly), fill = "white", col = "grey66", linewidth = 0.1) +
  geom_tile(sdm_df, mapping = aes(x = x, y = y, fill = sum)) +
  geom_sf(data = st_intersection(countries, recs_poly), fill = NA, col = "grey66", linewidth = 0.1) +
  geom_point(data = uk_midwinter_pts, mapping = aes(x = Lon, y = Lat), shape = 21, size = 3, fill = NA, stroke = 0.4) +
  geom_point(data = nigal.hahn2, mapping = aes(x = wlon1, y = wlat1, shape = as.factor(popunique)), col = "black", size = 3, fill = NA, stroke = 0.4) +
  geom_point(data = uk_midwinter_pts[!(uk_midwinter_pts$ID %in% suitable_points$ID),], mapping = aes(x = Lon, Lat), shape = 21, size = 3, stroke = 0.3, col= "red", fill = "red",  alpha = 0.6) +
  geom_point(data = nigal.hahn2[!(nigal.hahn2$ID %in% suitable_points$ID),], mapping = aes(x = wlon1, wlat1, shape = as.factor(popunique)), size = 3, stroke = 0.3, col= "red", fill = "red", alpha = 0.6) +
  scale_y_continuous(expand = c(0,0), limits = c(min(sdm_df$y), max(sdm_df$y)), labels = seq(4, 16, by = 2)) +
  scale_x_continuous(expand = c(0,0), limits = c(min(sdm_df$x)-1, max(sdm_df$x))) +
  scale_fill_gradientn(colors = rev(c("#006837", "#00A600", "#1DB000", "#3EBB00", "#63C600", "#8BD000", "#B6DB00", "#E6E600", "#E7CE1D", "#E9BD3A", "#EAB358", "white"))) +
  scale_shape_manual(values = c(22, 24, 21)) +
  theme_classic() + 
  xlab("") + ylab("Latitude (°N)") +
  guides(fill = "none", shape = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank()) 

## Changing habitat suitability ####
setwd("~/BTO projects/nightingale migration/sdm")
change_rast <- rast("changing-suitability.tif") 
change_rast_cr <- change_rast %>% crop(recs_poly) %>% terra::mask(recs_poly)
change_df <- as.data.frame(change_rast_cr, xy = T)

p2 <- ggplot() + 
  geom_sf(data = st_intersection(countries,recs_poly), fill = "white", col = "grey66", linewidth = 0.1) +  geom_tile(change_df, mapping = aes(x = x, y = y, fill = sum)) +
  geom_sf(data = st_intersection(countries, recs_poly), fill = NA, col = "grey66", linewidth = 0.1) +
  geom_point(data = uk_midwinter_pts, mapping = aes(x = Lon, y = Lat), shape = 21, size = 3, fill = NA, stroke = .4) +
  geom_point(data = nigal.hahn2, mapping = aes(x = wlon1, y = wlat1, shape = as.factor(popunique)), col = "black", size = 3, fill = NA, stroke = .4) +
  scale_y_continuous(expand = c(0,0), limits = c(min(change_df$y),max(change_df$y)), labels = seq(4, 16, by = 2)) +
  scale_x_continuous(expand = c(0,0), limits = c(min(change_df$x)-1,max(change_df$x)), labels = seq(-15, 10, by = 5)) +
  scale_fill_gradientn(colors = rev(c('#7f3b08','#8c510a','#b35806','#e08214','#fdb863','#fee0b6',
                                      'white','#d8daeb','#8073ac','#542788','#2d004b'))) +
  scale_shape_manual(values = c(22, 24, 21)) +
  theme_classic() + 
  xlab("Longtitude (°E)") + ylab("Latitude (°N)") +
  guides(fill = "none", shape = "none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank()) 

tiff("suitability maps.tif", width = 150, units = "mm", height = 120,
     res = 900)
library(ggpubr)
ggarrange(p1, p2, ncol = 1)
dev.off()

# Summarise past suitability 
past_rast <- rast("past-suitability.tif")
past_df <- as.data.frame(terra::extract(past_rast, all_points_sf, fun = table, bind = T, na.rm = T))
past_df <- gather(past_df, suit, freq, 2:ncol(past_df))
past_sum <- past_df %>% 
  mutate(
    Ring.No = all_points_sf$Ring.No[as.numeric(ID)],
    country = all_points_sf$country[as.numeric(ID)]
  ) %>%
  # Replace missing country values with "UK"
  mutate(country = ifelse(is.na(country), "UK", country)) %>%
  # Filter out rows where country is "Bulgaria" & frequency is 1
  filter(country != "Bulgaria" & freq == 1) %>%
  # Calculate summary statistics
  group_by(country, ID) %>%
  dplyr::summarize(
    mean = mean(as.numeric(suit)),
    prop = sum(freq[as.numeric(suit) >= 0.03]) / sum(freq),
    sd = sd(suit),
    se = sd / sqrt(length(unique(Ring.No))),
    up.ci = mean + se * 1.405072,
    low.ci = mean - se * 1.405072,
    Ring.No = first(Ring.No),
    year = 2001
  )

# Current suitability
curr_rast <- rast("current-suitability.tif")
curr_df <- as.data.frame(terra::extract(curr_rast, all_points_sf, fun = table, bind = T, na.rm = T))
curr_df <- gather(curr_df, suit, freq, 2:ncol(curr_df))
curr_sum <- curr_df %>%
  mutate(
    Ring.No = all_points_sf$Ring.No[as.numeric(ID)],
    country = all_points_sf$country[as.numeric(ID)]
  ) %>%
  # Replace missing country values with "UK"
  mutate(country = ifelse(is.na(country), "UK", country)) %>%
  # Filter out rows where country is "Bulgaria" and frequency is 1
  filter(country != "Bulgaria" & freq == 1) %>%
  # Calculate summary statistics
  group_by(country, ID) %>%
  dplyr::summarize(
    mean = mean(as.numeric(suit)),
    prop = sum(freq[as.numeric(suit) >= 0.03]) / sum(freq),
    sd = sd(suit),
    se = sd / sqrt(length(unique(Ring.No))),
    up.ci = mean + se * 1.405072,
    low.ci = mean - se * 1.405072,
    Ring.No = first(Ring.No),
    year = 2020
  )

# Combine current and past suitability into one dataframe
curr_past <- rbind(curr_sum, past_sum)

emmeans(lmer(mean ~ country*as.factor(year) + (1|Ring.No), data = curr_past), list(pairwise ~ as.factor(year)|country), adjust = "tukey")
## No change in proportion over time 
(mod_res <- emmeans(lmer(mean ~ country*as.factor(year) + (1|Ring.No), data = curr_past), list(pairwise ~ country|as.factor(year)), adjust = "tukey", level = .85))
mod_res_df <- as.data.frame(mod_res$`emmeans of country | year`)

# Plot
rbind(past_df %>% 
        mutate(
          suit = as.numeric(suit),
          country = all_points_sf$country[as.numeric(ID)]
        ) %>%
        # Replace missing country values with "UK"
        mutate(country = ifelse(is.na(country), "UK", country),
               year = 2001) %>%
        # Filter out rows where country is "Bulgaria" & frequency is 1
        filter(country != "Bulgaria" & freq == 1),
        curr_df %>% 
          mutate(
            suit = as.numeric(suit),
            country = all_points_sf$country[as.numeric(ID)]
        ) %>%
        # Replace missing country values with "UK"
        mutate(country = ifelse(is.na(country), "UK", country),
               year = 2020) %>%
        # Filter out rows where country is "Bulgaria" & frequency is 1
        filter(country != "Bulgaria" & freq == 1)) %>% 
ggplot() +
  geom_violin(mapping = aes(x = as.factor(year), y = suit, group = interaction(year,factor(country, levels = c("UK", "France", "Italy"))), col = factor(country, levels = c("UK", "France", "Italy"))), bw = .01, width = .5) +
  geom_point(mod_res_df, mapping = aes(as.factor(year), y = emmean, group = interaction(year,factor(country, levels = c("UK", "France", "Italy")))), position = position_dodge(width = .5), size = .7) +
  geom_errorbar(mod_res_df, mapping = aes(as.factor(year), ymax = upper.CL, ymin = lower.CL, group = interaction(year,factor(country, levels = c("UK", "France", "Italy")))), width = 0, position = position_dodge(width = .5)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line()) +
  xlab("") + ylab("Suitability") +
  geom_hline(yintercept = .03, linetype = "dashed") +
  scale_colour_manual(values = c('#1b9e77','#d95f02','#7570b3'), labels = c("UK", "France", "Italy")) +
  labs(col = NULL) +
  scale_x_discrete(labels = c("2001", "2020")) 

## Wintering locations without land mask #### 
wgr.no.mask <- rbind.fill(
  uk %>%
    filter(site == "wgr" &
             mask == "none") %>%
    # Some birds have to years of data so split by year
    mutate(year = substr(Arrival, 1, 4)) %>%
    group_by(ID) %>%
    # Remove stationary periods after 15th March
    filter(!(Arrival > as.POSIXct(
      paste0(as.numeric(substr(min(
        Arrival
      ), 1, 4)) + 1, "-03-15"), tz = "GMT"
    )) &
      # Remove stationary periods before 15th October
      !(Departure < as.POSIXct(
      paste0(as.numeric(substr(min(
        Arrival
      ), 1, 4)), "-10-15"), tz = "GMT"
    ))) %>%
  mutate(year = min(year)) %>%
  group_by(ID, year) %>%
  # Number sites for each bird
  dplyr::mutate(
    site.num = row_number(),
    # Calculate number of sites
    nsite = length(site.num),
    # Identify mid winter location, using period defined by Hahn et al.
    overlap.mid = ifelse(
      Arrival < paste0(substr(Arrival, 1, 4), "/12/31") &
        Departure >= paste0(substr(Arrival, 1, 4), "/12/31"),
      "Y",
      "N"
    )
  ) %>%
  filter(site.num == 1 |
           overlap.mid == "Y")

,
african %>% filter(site == "brgr" & mask == "none")))

wgr.no.mask$group <- ifelse(
  wgr.no.mask$tag.loc == "Ghana",
  "Ghana",
  ifelse(wgr.no.mask$tag.loc == "Gambia", "Gambia", "UK")
)

ggplot() +
  geom_sf(data = countries, fill = "grey90", col = "grey50", linewidth = .1) +
  geom_pointrange(data = wgr.no.mask, mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .4) +
  geom_errorbarh(data = wgr.no.mask, mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .4) +
  geom_point(data = wgr.no.mask, mapping = aes(x = Lon, y = Lat, fill = group), pch = 21, size = 2, stroke = .2) + 
  scale_fill_manual(values = c("red", "darkgreen", "blue")) +
  geom_point(mapping = aes(x = -16.7666672, y = 13.1), size = 2, fill = "yellow", pch = 24, stroke = .2) +
  geom_point(mapping = aes(x = -2.464, y = 7.402), size = 2, fill = "yellow", pch = 24, stroke = .2) +
  xlim(-25,5) + 
  ylim(-29,50) +
  xlab("") + ylab("") +
  labs(fill = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank()) 

brgr.no.mask <- rbind.fill(
  uk %>%
    filter(site == "brgr" &
             mask == "none"),
  african %>% filter(site == "brgr" & mask == "none")) %>% mutate(group = ifelse(
  tag.loc == "Ghana",
  "Ghana",
  ifelse(tag.loc == "Gambia", "Gambia", "UK")
))

ggplot() +
  geom_sf(data = countries, fill = "grey90", col = "grey50", linewidth = .1) +
  geom_pointrange(data = wgr.no.mask, mapping = aes(x = Lon, y = Lat, ymin = Lat.2.5., ymax = Lat.97.25.), size = .2, alpha = .4) +
  geom_errorbarh(data = wgr.no.mask, mapping = aes(y = Lat, xmin = Lon.2.5., xmax = Lon.97.25.), linewidth = .7, alpha = .4) +
  geom_point(data = wgr.no.mask, mapping = aes(x = Lon, y = Lat, fill = group), pch = 21, size = 2, stroke = .2) + 
  scale_fill_manual(values = c("red", "darkgreen", "blue")) +
  geom_point(mapping = aes(x = -16.7666672, y = 13.1), size = 2, fill = "yellow", pch = 24, stroke = .2) +
  geom_point(mapping = aes(x = -2.464, y = 7.402), size = 2, fill = "yellow", pch = 24, stroke = .2) +
  xlim(-25,5) + 
  ylim(-29,50) +
  xlab("") + ylab("") +
  labs(fill = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line = element_blank()) 
