############################################################
# HABITAT SUITABILITY MODELLING
# M. Kirkland
# Created 23/01/24
############################################################

# Clear workspace
rm(list = ls())

## devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")

# Load libraries
library(RColorBrewer)
library(corrplot)
library(dismo)
library(ecospat) ## For ecospat.mantel.correlogram()
library(randomForest)
library(terra)
library(sf)
library(PresenceAbsence) ## For calculating thresholds
library(mecofun) ## For crossvalSDM()
library(sdm) ## For fitting ensemble models
library(ggplot2)
library(geojsonsf)
library(ggthemes)
library(car)
library(mgcv)
library(data.table)
library(dplyr)
library(dggridR) ## To create hexagonal grid
library(lubridate) ## For year() and week() functions
library(tidyr)
library(stringr)
library(Boruta)

setwd("./sdm") ## Where data for SDM are shored

# Read in dataset of from GEE
files <- list.files(path = "./sdm/GEE results", pattern='.csv', 
           all.files = T, full.names = T)
gee <- lapply(files, read.csv)
gee.df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = F),
               gee)
lc.df <- read.csv("nigal-lc-prop.csv")
nigal <- merge(gee.df, lc.df)

# DATA PREP #########################################

# Remove missing data
nigal <- nigal[complete.cases(nigal),]
# Extract coordinates
nigal <- cbind(st_coordinates(st_as_sf(geojson_sf(nigal$.geo))),nigal)

# Calculate maximum and minimum monthly NDVI and rainfall
nigal <- transform(nigal, maxNDVI = apply(nigal[,grep("NDVI", colnames(nigal))], 1, max, na.rm = TRUE),
                          minNDVI = apply(nigal[,grep("NDVI", colnames(nigal))], 1, min, na.rm = TRUE),
                          maxRain = apply(nigal[,grep("Rain", colnames(nigal))], 1, max, na.rm = TRUE),
                          minRain = apply(nigal[,grep("Rain", colnames(nigal))], 1, min, na.rm = TRUE),
                          summerRain = apply(subset(nigal, select = c(JuneRain, JulyRain, AugRain, 
                                              SepRain, OctRain)), 1, sum, na.rm = TRUE),
                          winterRain= apply(subset(nigal, select = c(NovRain, DecRain, JanRain, 
                                                  FebRain)), 1, sum, na.rm = TRUE))
                          
# Plot records
countries <- st_read(dsn = "~/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp") # Read in country data 
ggplot() +
  geom_sf(data = countries, fill = "grey88", col = "grey66", linewidth = .5, linetype = "dashed") +
  geom_point(data = nigal %>% arrange(presence), mapping = aes(x = X, y = Y, col = as.factor(presence)), size = 2, alpha = .7) +
  xlim(-26,30) + 
  ylim(3,22) +
  theme(legend.title = element_blank(),
        legend.position = "left") +
  theme_map() +
  labs(col = "") +
  scale_colour_manual(values = c("#6a3d9a", "#ff7f00"), labels = c("Absence", "Presence"))

# Remove everything west of Cameroon
nigal <- nigal[nigal$X < 15,]

# Reduce sample bias in citizen science projects ####
# Create 4km hexagonal grid 
dggs <- dgconstruct(spacing = 4)
nigal_sub <- nigal %>% 
  filter(study == "The African Bird Atlas" | study == "eBird") %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, X, Y)$seqnum,
         year = year(date),
         week = week(date)) %>% 
  group_by(presence, year, week, cell) %>% ## Select one checklist per week per cell. Group by presence to sample them seperately so you don't lose too many detections
  sample_n(size = 1) %>% 
  ungroup() %>% 
  dplyr::select(-c(year, week, cell))
nigal_sub <- rbind(nigal %>% filter(study != "The African Bird Atlas" & study != "eBird"), nigal_sub)

# Plot this
ggplot() +
  geom_sf(data = countries, fill = "grey88", col = "grey66", linewidth = .5, linetype = "dashed") +
  geom_point(data = nigal_sub %>% arrange(presence), mapping = aes(x = X, y = Y, col = as.factor(presence)), size = 2, alpha = .7) +
  xlim(-26,15) + 
  ylim(3,22) +
  theme(legend.title = element_blank(),
        legend.position = "left") +
  theme_map() +
  labs(col = "") +
  scale_colour_manual(values = c("#6a3d9a", "#ff7f00"), labels = c("Absence", "Presence"))

table(nigal_sub$presence)

# Select presence and absence data
nigal_occ <- nigal_sub[nigal_sub$presence == 1,]
backgr <- nigal_sub[nigal_sub$presence == 0,]

# Test for spatial autocorrelation in occurrence dataset
ecospat.mantel.correlogram(dfvar=nigal_occ,
                           colxy=1:2, n=100,
                           colvar=11:40, max=200, nclass=20, nperm=100)

# Check correlation between variables
nums <- unlist(lapply(nigal_sub, is.numeric), use.names = FALSE)  
num.df <- nigal_sub[ , nums]
# Remove coordinates and ID column 
num.df <- num.df[,-(1:3)]
cor <- cor(num.df)
options(max.print=1000000) 
round(cor,3)

# Calcuate p-value
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# Matrix of the p-value of the correlation
p.mat <- cor.mtest(num.df)

tmp = cor # Copy matrix
tmp[ tmp > -0.7 & tmp < 0.7 ] = 0

# Plot correlogram
corrplot(tmp, type="upper", diag = F, # Add coefficient of correlation
         tl.col="black", tl.cex = .7, 
         col=rev(brewer.pal(n=8, name="PuOr"))) # Text label color and rotation

# MODELS #########################################

fold_pres <- kfold(nigal_occ, k = 5) # Add an index that makes five random groups of observations
nigal_test <- nigal_occ[fold_pres == 1, ] 
nigal_train <- nigal_occ[fold_pres != 1, ]
fold_backgr <- kfold(backgr, k = 5) # Add an index that makes five random groups of observations
backgr_test <- backgr[fold_backgr == 1, ] 
backgr_train <- backgr[fold_backgr != 1, ]

# Make a dataset that combines the presence and absence data, the background points with their model environment value, and the occurrence points with their model environment values
sdmdata <- data.frame(rbind(backgr_train, nigal_train))

# This was run on full model to identify important variables for each model
# fs <- Boruta(rmod1, sdmdata, doTrace = 2, maxRuns = 100)
# print(fs)
# plot(fs)

# Run alternative models
rmod1 <-
  presence ~ pr + JanRain + FebRain + OctRain + NovRain + DecRain + JulyRain + JuneRain + AugRain + SepRain + MarchRain + AprRain + MayRain + JanNDVI + FebNDVI + OctNDVI + NovNDVI + DecNDVI + JulyNDVI + JuneNDVI + AugNDVI + SepNDVI + MarchNDVI + AprNDVI + MayNDVI + elevation + pdsi + tmmx + tmmn + Mosaic.natural.vegetation + Shrubland + Urban
rf_test1 <-
  randomForest(rmod1,
               data = sdmdata,
               importance = T,
               ntree = 500)

rmod2 <-
  presence ~ pr + maxRain + minRain + MayRain + AprRain + MarchRain + winterRain + summerRain + maxNDVI + minNDVI + elevation + pdsi + tmmx + tmmn + Mosaic.natural.vegetation + Shrubland + Forest + Urban + Cropland..rainfed 
rf_test2 <-
  randomForest(rmod2,
               data = sdmdata,
               importance = T,
               ntree = 500)

rmod3 <-
  presence ~ pr + JanRain + FebRain + OctRain + NovRain + DecRain + JulyRain + JuneRain + AugRain + SepRain + MarchRain + AprRain + MayRain + maxNDVI + minNDVI + elevation + pdsi + tmmx + tmmn + Mosaic.natural.vegetation + Shrubland 
rf_test3 <-
  randomForest(rmod3,
               data = sdmdata,
               importance = T,
               ntree = 500)

rmod4 <-
  presence ~ pr + winterRain + summerRain + minRain + maxRain + JanNDVI + FebNDVI + OctNDVI + NovNDVI + DecNDVI + JulyNDVI + JuneNDVI + AugNDVI + SepNDVI + MarchNDVI + AprNDVI + MayNDVI + elevation + pdsi + tmmx + tmmn + Mosaic.natural.vegetation + Shrubland + Forest + Urban 
rf_test4 <-
  randomForest(rmod4,
               data = sdmdata,
               importance = T,
               ntree = 500)

auc1 <- evaluate(nigal_test, backgr_test, rf_test1)@auc
auc2 <- evaluate(nigal_test, backgr_test, rf_test2)@auc
auc3 <- evaluate(nigal_test, backgr_test, rf_test3)@auc
auc4 <- evaluate(nigal_test, backgr_test, rf_test4)@auc

# Now run each model 10 times with different training datasets

# Creating the objects  
run = 10 # 10 repeats of the model
thresh_maxSSS <- vector("numeric", run*4) # Multiply by 4 because 4 different models
auc <- vector("numeric", run*4)
preds_rf <- list()
fitted <- list()
tss <- vector("numeric", run*4)
expl_dev <- vector("numeric", run*4)
l_mod <- list()

# Run each model iteratively 10 times
for(i in 1:run){
  fold_pres <- kfold(nigal_occ, k = 5) # Add an index that makes five random groups of observations
  nigal_test <- nigal_occ[fold_pres == 1, ] 
  nigal_train <- nigal_occ[fold_pres != 1, ]
  fold_backgr <- kfold(backgr, k = 5) # Add an index that makes five random groups of observations
  backgr_test <- backgr[fold_backgr == 1, ] 
  backgr_train <- backgr[fold_backgr != 1, ]
  
  # Make a dataset that combines the presence and absence data, the background points with their model environment value, and the occurrence points with their model environment values
  sdmdata <- data.frame(rbind(backgr_train, nigal_train))
  
  # Subset of that dataset, select variables of interest

  # Fit random forests
  rf1 <-
    randomForest(rmod1,
                 data = sdmdata,
                 importance = T,
                 ntree = 500)
  l_mod[[i]] <- rf1
  names(l_mod)[i] <- paste0("1-",i) # Name list
  
  rf2 <-
    randomForest(rmod2,
                 data = sdmdata,
                 importance = T,
                 ntree = 500)
  l_mod[[i+run]] <- rf2
  names(l_mod)[i+run] <- paste0("2-",i) # Name list
  
  rf3 <-
    randomForest(rmod3,
                 data = sdmdata,
                 importance = T,
                 ntree = 500)
  l_mod[[i+run*2]] <- rf3
  names(l_mod)[i+run*2] <- paste0("3-",i) # Name list
  
  rf4 <-
    randomForest(rmod4,
                 data = sdmdata,
                 importance = T,
                 ntree = 500)
  l_mod[[i+run*3]] <- rf4
  names(l_mod)[i+run*3] <- paste0("4-",i) # Name list
  
  # Store fitted values
  fitted[[i]] <- rf1$predicted
  fitted[[i+run]] <- rf2$predicted
  fitted[[i+run*2]] <- rf3$predicted
  fitted[[i+run*3]] <- rf4$predicted
  
  # Model evaluation
  erf1 <- evaluate(nigal_test, backgr_test, rf1)
  erf2 <- evaluate(nigal_test, backgr_test, rf2)
  erf3 <- evaluate(nigal_test, backgr_test, rf3)
  erf4 <- evaluate(nigal_test, backgr_test, rf4)
  
  # Capturing the AUC values
  auc[i] <- erf1@auc
  auc[i+run] <- erf2@auc
  auc[i+run*2] <- erf3@auc
  auc[i+run*3] <- erf4@auc
  
  # Explained deviance 
  expl_dev[i] <- expl_deviance(obs = sdmdata$presence,
                                    pred = rf1$predicted)
  expl_dev[i+run] <- expl_deviance(obs = sdmdata$presence,
                               pred = rf2$predicted)
  expl_dev[i+run*2] <- expl_deviance(obs = sdmdata$presence,
                               pred = rf3$predicted)
  expl_dev[i+run*3] <- expl_deviance(obs = sdmdata$presence,
                               pred = rf4$predicted)
  # Numeric vector of cross-validated prediction
  preds_cv1 <- crossvalSDM(rf1, traindat = sdmdata, colname_species = 'presence', 
                          colname_pred = row.names(rf_test1$importance))
  preds_cv2 <- crossvalSDM(rf2, traindat = sdmdata, colname_species = 'presence', 
                           colname_pred = row.names(rf_test2$importance))
  preds_cv3 <- crossvalSDM(rf3, traindat = sdmdata, colname_species = 'presence', 
                           colname_pred = row.names(rf_test3$importance))
  preds_cv4 <- crossvalSDM(rf4, traindat = sdmdata, colname_species = 'presence', 
                           colname_pred = row.names(rf_test4$importance))
  
  preds_rf[[i]] <- preds_cv1
  preds_rf[[i+run]] <- preds_cv2
  preds_rf[[i+run*2]] <- preds_cv3
  preds_rf[[i+run*3]] <- preds_cv4
  
  # Prepare cross-validation dataset
  thresh_df1 <- data.frame(ID = seq_len(nrow(sdmdata)), obs = sdmdata$presence, pred = preds_cv1)
  thresh_df2 <- data.frame(ID = seq_len(nrow(sdmdata)), obs = sdmdata$presence, pred = preds_cv2)
  thresh_df3 <- data.frame(ID = seq_len(nrow(sdmdata)), obs = sdmdata$presence, pred = preds_cv3)
  thresh_df4 <- data.frame(ID = seq_len(nrow(sdmdata)), obs = sdmdata$presence, pred = preds_cv4)
  # Find the optimal thresholds    
  thresh1 <- optimal.thresholds(DATA = thresh_df1)
  thresh_maxSSS[i] <- thresh1[3,2] 
  cmx_maxSSS1 <- cmx(DATA = thresh_df1, threshold = thresh1[3,2])
  tss[i] <- TSS(cmx_maxSSS1)
  
  thresh2 <- optimal.thresholds(DATA = thresh_df2)
  thresh_maxSSS[i+run] <- thresh2[3,2] 
  cmx_maxSSS2 <- cmx(DATA = thresh_df2, threshold = thresh2[3,2])
  tss[i+run] <- TSS(cmx_maxSSS2)
  
  thresh3 <- optimal.thresholds(DATA = thresh_df3)
  thresh_maxSSS[i+run*2] <- thresh3[3,2] 
  cmx_maxSSS3 <- cmx(DATA = thresh_df3, threshold = thresh3[3,2])
  tss[i+run*2] <- TSS(cmx_maxSSS3)
  
  thresh4 <- optimal.thresholds(DATA = thresh_df4)
  thresh_maxSSS[i+run*3] <- thresh4[3,2] 
  cmx_maxSSS4 <- cmx(DATA = thresh_df4, threshold = thresh4[3,2])
  tss[i+run*3] <- TSS(cmx_maxSSS4)
}

# Deviance explained. How well does the model fit the data? 
mean(expl_dev)
plot(unlist(fitted), unlist(preds_rf), xlab = 'Fitted values', ylab = 'Predicted values from CV', main = "")
abline(0,1,col='red',lwd=2)

# How robust is the model against changes in the input data and, thus, how robust are predictions to other places and times?
## NB The use of AUC in evaluating SDMs has been criticized (Lobo et al. 2008, Jiménez-Valverde 2011).
mean(auc)
mean(tss)

# Weighted average threshold?
w <- (auc-0.5)^2
sum(thresh_maxSSS[1:10]*w[1:10],thresh_maxSSS[11:20]*w[11:20],thresh_maxSSS[21:30]*w[21:30],thresh_maxSSS[31:40]*w[31:40])/sum(w)

# PREDICTIONS #########################################

# Raster for making predictions. Imported from GEE. 
rasts <- rast("prediction-rasters-4km.tif")

# Add in maximum and minimum NDVI and rain variables
maxNDVI <- app(rasts[[grep("NDVI", names(rasts))]], max); names(maxNDVI) <- "maxNDVI"
minNDVI <- app(rasts[[grep("NDVI", names(rasts))]], min); names(minNDVI) <- "minNDVI"
pr <- app(rasts[[grep("Rain", names(rasts))]], sum); names(pr) <- "pr"
maxRain <- app(rasts[[grep("Rain", names(rasts))]], max); names(maxRain) <- "maxRain"
minRain <- app(rasts[[grep("Rain", names(rasts))]], min); names(minRain) <- "minRain"
summerRain <- rasts[["JuneRain"]] + rasts[["JulyRain"]] + rasts[["AugRain"]] + rasts[["SepRain"]] +
  rasts[["OctRain"]]; names(summerRain) <- "summerRain"
winterRain <- rasts[["NovRain"]] + rasts[["DecRain"]] + rasts[["JanRain"]] + rasts[["FebRain"]]; names(winterRain) <- "winterRain"
rasts <- c(rasts, maxNDVI, minNDVI, pr, maxRain, minRain, winterRain, summerRain)

# Crop to Western Africa, as the model is not very good at predicting in Central Africa
rasts <- crop(rasts, ext(c(-18.0416195374957, 15, 2.99999216327873, 17.5)))

# Stack model predictions
pred_stck <- rast()
for (i in 1:length(l_mod)) {
  pred <- predict(rasts, l_mod[[i]], type = "response")
  pred_stck <- c(pred_stck, pred)
}

# Calculate mean predictions for each model weighted by AUC of each model
pred_avg <- weighted.mean(pred_stck, w)

# Plot predictions  
plot(pred_avg, main = "Weighted mean Random Forest predictions")
plot(countries, add = T, col = NA, lwd = .2)  

# Export
writeRaster(pred_avg, "results-rasters.tif", overwrite = T)
save.image("sdm.RData")

# Look at suitability change 
rasts01 <- rast("2001-prediction-rasts.tif")
rasts20 <- rast("2020-prediction-rasts.tif")
rasts01 <- crop(rasts01, ext(c(-18.0416195374957, 15, 2.99999216327873, 17.5)))
rasts20 <- crop(rasts20, ext(c(-18.0416195374957, 15, 2.99999216327873, 17.5)))

# Add in maximum and minimum NDVI and rain variables
maxNDVI <- app(rasts01[[grep("NDVI", names(rasts01))]], max); names(maxNDVI) <- "maxNDVI"
minNDVI <- app(rasts01[[grep("NDVI", names(rasts01))]], min); names(minNDVI) <- "minNDVI"
pr <- app(rasts01[[grep("Rain", names(rasts01))]], mean); names(pr) <- "pr"
maxRain <- app(rasts01[[grep("Rain", names(rasts01))]], max); names(maxRain) <- "maxRain"
minRain <- app(rasts01[[grep("Rain", names(rasts01))]], min); names(minRain) <- "minRain"
summerRain <- rasts01[["JuneRain"]] + rasts01[["JulyRain"]] + rasts01[["AugRain"]] + rasts01[["SepRain"]] +
  rasts01[["OctRain"]]; names(summerRain) <- "summerRain"
winterRain <- rasts01[["NovRain"]] + rasts01[["DecRain"]] + rasts01[["JanRain"]] + rasts01[["FebRain"]]; names(winterRain) <- "winterRain"
rasts01 <- c(rasts01, maxNDVI, minNDVI, pr, minRain, maxRain, winterRain, summerRain)

maxNDVI <- app(rasts20[[grep("NDVI", names(rasts20))]], max); names(maxNDVI) <- "maxNDVI"
minNDVI <- app(rasts20[[grep("NDVI", names(rasts20))]], min); names(minNDVI) <- "minNDVI"
pr <- app(rasts20[[grep("Rain", names(rasts20))]], sum); names(pr) <- "pr"
maxRain <- app(rasts20[[grep("Rain", names(rasts20))]], max); names(maxRain) <- "maxRain"
minRain <- app(rasts20[[grep("Rain", names(rasts20))]], min); names(minRain) <- "minRain"
summerRain <- rasts20[["JuneRain"]] + rasts20[["JulyRain"]] + rasts20[["AugRain"]] + rasts20[["SepRain"]] +
  rasts20[["OctRain"]]; names(summerRain) <- "summerRain"
winterRain <- rasts20[["NovRain"]] + rasts20[["DecRain"]] + rasts20[["JanRain"]] + rasts20[["FebRain"]]; names(winterRain) <- "winterRain"
rasts20 <- c(rasts20, maxNDVI, minNDVI, pr, minRain, maxRain, winterRain, summerRain) 

stck01 <- rast()
for (i in 1:length(l_mod)) { 
  pred01 <- predict(rasts01, l_mod[[i]], type = "response")
  stck01 <- c(stck01, pred01)
}

# Average predictions from random forest
pred_2001 <- weighted.mean(stck01, w)
plot(pred_2001)
# Export 
writeRaster(pred_2001, "past-suitability.tif", overwrite = T)

# Repeat for latest year
stck20 <- rast()
for (i in 1:length(l_mod)) { 
  pred20 <- predict(rasts20, l_mod[[i]], type = "response")
  stck20 <- c(stck20, pred20)
}
pred_2020 <- weighted.mean(stck20, w)
plot(pred_2020)
writeRaster(pred_2020, "current-suitability.tif", overwrite = T)

# Calculate changes and export as raster
change <- pred_2020-pred_2001
plot(change)
writeRaster(change, "changing-suitability.tif", overwrite = T)

# Plot relative change 
rel_change_rast <- round(change, digits = 4)/round(pred_2001, digits = 4)
rel_change_rast[rel_change_rast>0] <- NA
rel_change_rast[is.infinite(rel_change_rast)] <- NA
rel_change_rast <- rel_change_rast * -1
plot(rel_change_rast)
# Export
writeRaster(rel_change_rast, "relative-change.tif", overwrite = T)

# PLOT EFFECTS ###################################
 
library(stringr)

# We want two panels next to each other
par(mfrow=c(1,2))
# Plot the partial responses of one model of each type. Maybe consider model averaging? 
partialPlot(l_mod[[9]], pred.data = nigal_train, x.var = "Mosaic.natural.vegetation", ylab = 'Occurrence probability', xlab = "Forest", main = "")

# Get variable importance measures for random forest
imp_df <- list()
for (d in 1:length(l_mod)) {
  imp_df[[d]] <- data.frame(importance(l_mod[[d]], type = 1))
  imp_df[[d]]$names <- row.names(imp_df[[d]])
}
imp_df <- rbindlist(imp_df)[,lapply(.SD,mean),names] ## Average across models

### How to average response curves across models?
# impvar <- imp_df[order(imp_df[, 2], decreasing=TRUE)]
# op <- par(mfrow=c(2, 3))
# # Plot responses in order of decreasing importance 
# for (i in seq_along(impvar$name)) {
#   partialPlot(rf, nigal, impvar$name[i], xlab=impvar$name[i],
#               main=paste("Partial Dependence on", impvar$name[i]),
#               ylim=c(0,1))
# }
# par(op)

# Plot mean decreased accuracy
imp_df %>% 
  filter(!str_detect(names, "2")) %>% 
  ggplot(aes(x = reorder(names, X.IncMSE),y = X.IncMSE)) +
  geom_col(fill = "grey", col = "black", linewidth = .1) +
  coord_flip() +
  labs(x= "",
       y= "Mean Decrease in Accuracy") +
  theme_classic() +
  scale_x_discrete(labels = rev(c("Minimum NDVI", "Elevation", "PDSI", "Maximum precipitation", "Summer precipitation", "July precipitation", "Maximum temperature", "June precipitation", "August precipitation", "September precipitation", "May precipitation", "June NDVI", "Minimum temperature", "Annual precipitation", "November NDVI", "Winter rain", "April NDVI", "Maximum NDVI", "August NDVI", "March NDVI", "July NDVI", "October precipitation", "September NDVI", "November precipitation", "Mosaic natural vegetation", "May NDVI", "March precipitation", "February NDVI", "April precipitation", "January NDVI", "December NDVI", "Minimum precipitation", "Febraury precipitation", "January precipitation", "October NDVI", "December precipitation", "Shrubland", "Cropland, rainfed", "Forest", "Urban")))
