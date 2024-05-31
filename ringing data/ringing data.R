############################################################
# PLOT RINGING DATA
# M. Kirkland
# Created 06/11/23
############################################################

# Clear working directory
rm(list = ls())

# Load libraries
library(ggplot2)
library(sf)
library(ecodist)
library(dplyr)

# Read in ringing data from Kartong Bird Observatory
setwd("~/BTO projects/nightingale migration/ringing data")
data <- read.csv("KBO_DATA_NIGAL_20231012.csv")

countries <- st_read(dsn = "~/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp") # Read in country data

# Select birds recaptured elsewhere
recap <- unique(data[,c('RING', 'PLACE')])
recap <- recap$RING[duplicated(recap$RING)]
recap <- data[data$RING %in% recap,]
recap <- recap[!duplicated(recap[,c('RING', 'PLACE')]),]
recap <- recap %>% 
  group_by(RING) %>% 
  mutate(site = row_number()) # Add site number to differentiate between tagging location and recapture site

ggplot() +
  geom_sf(data = countries, fill = "grey88", col = "grey66", linewidth = .5, linetype = "dashed") +
  geom_point(data = recap, mapping = aes(x = LONG, y = LAT, shape = as.factor(site)), size = 2) +
  xlim(-20,10) + 
  ylim(10,55) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  facet_wrap(vars(RING)) +
  xlab("") + ylab("") +
  scale_shape_manual(values = c(4,16))

# Calculate distances between capture sites
meanDist <- function(lon, lat){
  max(lower(distm(cbind(lon, lat), fun = distHaversine)))
}
options(pillar.sigfig = 7) ## tibble ouput is 3 digits so increase this
recap %>% 
  filter(site == 2 & PLACE == "KARTON") %>% 
  group_by(site) %>% 
  summarise('meandist' = (meanDist(LONG, LAT))/1000) 

recap %>% 
  filter(site == 2 & PLACE != "KARTON") %>% 
  group_by(site) %>% 
  summarise('meandist' = (meanDist(LONG, LAT))/1000) 
