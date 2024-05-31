############################################################
# MANTELS CORRELATION OF POPULATION CONNECTIVITY
# M. Kirkland
# Created 02/04/24
############################################################

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(geosphere)
library(sf)

# Read in our nightingale data 
setwd("./geolocator data/results/wgr")
filenames <- list.files(getwd(), pattern="*.csv", full.names=TRUE)
data = tibble(file = filenames) %>%
  tidyr::extract(file, "ID", "(?<=wgr/)(.+)(?=\\_wgr.csv)", remove = FALSE) %>% # Add column with bird tag ID
  mutate(data = lapply(file, read.csv)) %>%
  unnest(data) %>%
  dplyr::select(-file)

## Data from Hahn et al. 
setwd(".")
nigal.hahn <- read.csv("LusMeg_F_IT_BG_KernelMedians for Chris_Nov2023.csv")
hahn.wgr <- nigal.hahn %>% 
  dplyr::filter(type == 3) %>% # Select main wintering sites
  dplyr::select(long,lat, population)
colnames(hahn.wgr)[1:3] <- c("Lon", "Lat", "population")

# Create distance matrix for breeding grounds 
brgr.uk <- data %>% 
  dplyr::filter(site == "brgr" & mask == "land") %>% 
  group_by(ID) %>% 
  dplyr::select(Lon, Lat) 

hahn.tag <- data.frame(blon = c(7.5, 11.8, 28.1), blat = c(47.6, 44.6, 42.75), population = c(1,2,3))
hahn.brgr <- merge(hahn.wgr, hahn.tag) 
hahn.brgr <- subset(hahn.brgr, select = c(blon,blat))
colnames(hahn.brgr) <- c("Lon","Lat")

# Select wintering grounds estimated using land mask for birds who also have breeding ground estimates
wgr.uk <- data %>% 
  dplyr::filter(site == "wgr" & mask == "land" & ID %in% brgr.uk$ID) %>% 
  group_by(ID) %>%
  # Remove stationary periods after 15th March
  filter(!(Arrival > as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)) + 1,"-03-15"), tz = "GMT")) & 
           # Remove stationary periods before 15th October
           !(Departure < as.POSIXct(paste0(as.numeric(substr(min(Arrival), 1,4)),"-10-15"), tz = "GMT"))) %>% 
  dplyr::select(Lon, Lat)

# Merge with breeding ground data since multiple wintering sites for individuals
brgr.uk <- merge(subset(wgr.uk, select = c(ID)), brgr.uk)

# Combine datasets
all.wgr <- bind_rows(subset(hahn.wgr, select = -population), subset(wgr.uk, select = -ID))
all.brgr <- bind_rows(hahn.brgr, subset(brgr.uk, select = -ID)) ## Remove ID, but keep in datatset to select only birds with breeding grounds estimates later when calculating distances at wintering groudns

# Distance matrix at the breeding grounds
d.brgr = distm(all.brgr, fun = distHaversine)
dist.brgr = as.dist(d.brgr)
# Distance matrix at the wintering grounds 
d.wgr = distm(all.wgr, fun = distHaversine)
dist.wgr = as.dist(d.wgr)

# Mantels correaltion
mantel(dist.wgr, dist.brgr, method = "spearman", permutations = 9999, na.rm = TRUE)
## This quantifies whether distances between individual breeding sites are maintained during the non-breeding season. Strong positive Mantel coefficients indicate that individuals which breed close together also spend the non-breeding season relatively close together, and vice versa (i.e. low inter-population mixing).

# Read in breeding grounds of birds tagged in Africa
setwd("./geolocator data/results/brgr")
filenames <- list.files(getwd(), pattern="*.csv", full.names=TRUE)
african.brgr = tibble(file = filenames) %>%
  tidyr::extract(file, "ID", "(?<=brgr/)(.+)(?=\\_brgr.csv)", remove = FALSE) %>% # Add column with bird tag ID
  mutate(data = lapply(file, read.csv)) %>%
  unnest(data) %>%
  filter(site == "brgr" & mask == "land") %>% 
  dplyr::select(Lon, Lat, tag.loc, ID)

# Average locations for repeated inds
african.avg <- african.brgr
african.avg$indID <- african.avg$ID
african.avg$indID[african.avg$ID == "BV422"] <- african.avg$indID[african.avg$ID == "BR064"]
african.avg <- african.avg %>% 
  group_by(indID) %>% 
  dplyr::mutate(Lon = mean(Lon),
                Lat = mean(Lat)) %>% 
  dplyr::distinct(indID, .keep_all=TRUE) %>% 
  ungroup %>% 
  select(Lon, Lat, tag.loc)

african.wgr <- african.avg
# Add tagging location in Gambia
african.wgr$Lon <- -16.7666672
african.wgr$Lat <- 13.1
# Add tagging location in Ghana
african.wgr$Lon[african.wgr$tag.loc == "Ghana"] <- -2.464
african.wgr$Lat[african.wgr$tag.loc == "Ghana"] <- 7.402
  
african.avg <- subset(african.avg, select = -tag.loc)
african.wgr <- subset(african.wgr, select = -tag.loc)

# Read in ringing data from Kartong Bird Observatory
setwd("./ringing data")
ringdata <- read.csv("KBO_DATA_NIGAL_20231012.csv")

# Select birds recaptured elsewhere
recap <- unique(ringdata[,c('RING', 'PLACE')])
recap <- recap$RING[duplicated(recap$RING)]
recap <- ringdata[ringdata$RING %in% recap,]
recap <- recap[!duplicated(recap[,c('RING', 'PLACE')]),]
recap <- subset(recap, select = c(LONG, LAT))
colnames(recap) <- c("Lon", "Lat")
recap.brgr <- recap[recap$Lat > 50,]
recap.wgr <- recap[recap$Lat < 50,]

# Read in GPS tags
setwd("./gps data")
gps1 <- st_read("Obs150324_104804_Tag58279.kml")
gps2 <- st_read("Obs130324_154637_Tag58288.kml")
# Extract coorindates. Last row is the path taken so ignore. 
gps1_coords <- as.data.frame(gps1 %>% st_cast("MULTIPOINT") %>% st_coordinates())
gps2_coords <- as.data.frame(gps2 %>% st_cast("MULTIPOINT") %>% st_coordinates())
gps1_coords <- subset(gps1_coords, select = c("X","Y"))
gps2_coords <- subset(gps2_coords, select = c("X","Y"))
colnames(gps1_coords) <- c("Lon","Lat"); colnames(gps2_coords) <- c("Lon","Lat")
gps1_brgr <- gps1_coords[which.max(gps1_coords$Lat),]
gps2_brgr <- gps2_coords[which.max(gps2_coords$Lat),]
gps_wgr <- data.frame(Lon = -16.7666672,
                       Lat = 13.1)

# Combine datasets
all.wgr2 <- bind_rows(all.wgr, recap.wgr, african.wgr, gps_wgr, gps_wgr)
all.brgr2 <- bind_rows(all.brgr, recap.brgr, african.avg, gps1_brgr, gps2_brgr)

# Plot this
countries <- st_read(dsn = "~/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp") # Read in country uk 

ggplot() +
  geom_sf(data = countries, fill = "grey90", col = "grey50", linewidth = .1) +
  geom_point(data = all.brgr2, mapping = aes(x = Lon, y = Lat))
  
# Distance matrix at the breeding grounds
d.brgr2 = distm(all.brgr2, fun = distHaversine)
dist.brgr2 = as.dist(d.brgr2)
# Distance matrix at the wintering grounds 
d.wgr2 = distm(all.wgr2, fun = distHaversine)
dist.wgr2 = as.dist(d.wgr2)

mantel(dist.wgr2, dist.brgr2, method = "spearman", permutations = 9999, na.rm = TRUE)
