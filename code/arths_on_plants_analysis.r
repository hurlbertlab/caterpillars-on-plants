# Script for calculating average abundance, biomass, and occurrence per plant species, and for comparing broadly between native and alien plants

library(tidyverse)
library(RCurl)
library(rvest)
library(xml2)
library(glmmTMB)
library(png)
library(lme4)
library(interactions)
library(RColorBrewer)
library(ggpubr)
library(httr)
library(jsonlite)

# load functions
source('code/plant_analysis_functions.r')

# Read in latest CC fullDataset through 2024
cc = read.csv('data/fullDataset_2025-04-18.csv', header = T, quote = '\"')

# Loading plant files from data_repo

# GitHub API URL to list contents of a repo directory
api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-count-data/contents/plantSpecies"

# Send GET request
res <- GET(api_url)

# Parse JSON response
files_info <- fromJSON(content(res, "text"))

# Filter for files with "officialPlantList" in the name
official_data_links <- files_info %>%
  filter(grepl("officialPlantList", name)) %>%
  transmute(
    file_name = name,
    download_url = download_url
  )

inferred_data_links <- files_info %>%
  filter(grepl("inferredPlantNames", name)) %>%
  transmute(
    file_name = name,
    download_url = download_url
  )

plant_origin_links <- files_info %>%
  filter(grepl("plant_origin_status.csv", name)) %>%
  transmute(
    file_name = name,
    download_url = download_url
  )

officialPlantList = read.csv(official_data_links$download_url[nrow(official_data_links)])
inferredPlantNames = read.csv(inferred_data_links$download_url[nrow(inferred_data_links)])

plantOrigin = read.csv(plant_origin_links$download_url) %>%
  select(scientificName, nativeStatus, plantOrigin)

# The dataset for which we have plant species names with NameConfidence >= 2
# NOTE: sciName is the field that includes inferred scientific names in addition to official ones.
# Exclude data from Coweeta sites, where non-caterpillars were not recorded
ccPlants = cc %>%
  left_join(inferredPlantNames[, c('PlantFK', 'InferredSciName', 'NameConfidence')], by = 'PlantFK') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  mutate(sciName = ifelse(Species == "N/A" & NameConfidence >= 2, InferredSciName, sciName)) %>%
  filter(!is.na(sciName),
         !Name %in% c('Coweeta - BB', 'Coweeta - BS', 'Coweeta - RK'),
         Year <= 2024,
         Longitude > -100)

# Some summary statistics for the dataset used in analysis
jdRange = c(152,212)

analysisdata = ccPlants %>%
                 filter(julianday >= jdRange[1], 
                        julianday <= jdRange[2],
                        !WetLeaves)
nSurvs = length(unique(analysisdata$ID)) # 68741
nSites = length(unique(analysisdata$Name)) # 212
range(analysisdata$Latitude) # 32.33, 55.43 (47.78 east of -100W)
range(analysisdata$Longitude) # -97.41, -68.05
range(analysisdata$Year) # 2010, 2024
nBranches = length(unique(analysisdata$Code)) # 5438
nPlantSpecies = length(unique(analysisdata$sciName)) # 363
nNativePlantSpecies = length(unique(analysisdata$sciName[analysisdata$plantOrigin == 'native'])) # 284
nAlienPlantSpecies = length(unique(analysisdata$sciName[analysisdata$plantOrigin == 'alien'])) # 107
analysisdata %>% 
  distinct(ID, ObservationMethod) %>% 
  count(ObservationMethod)                 # 29812 beat sheet surveys, 38929 visual surveys


# Labels and plot colors for arthropod groups
arthropods = data.frame(Group = c('caterpillar', 'spider', 'leafhopper', 'beetle', 'truebugs', 'ant'),
                        GroupLabel = c('caterpillars', 'spiders', 'hoppers', 'beetles', 'true bugs', 'ants'),
                        color = c('limegreen', 'turquoise2', 'dodgerblue', 
                                  'salmon', 'magenta3', 'orange'),
                        color2 = c(rgb(.6, .9, .6), rgb(.6, .96, .98), rgb(.56, .78, 1),
                                   rgb(.99, .75, .72), rgb(.93, 0.6, .93), rgb(1, .85, .4)))




# Color scheme for tree families
treeFams = data.frame(Family = c('Fagaceae', 'Betulaceae', 'Sapindaceae', 'Caprifoliaceae', 'Juglandaceae', 'Rosaceae'),
                      famcolor = c(rgb(230/255, 159/255, 0),
                        rgb(86/255, 180/255, 233/255),
                        rgb(0, 158/255, 115/255),
                        rgb(240/255, 228/255, 66/255),
                        rgb(213/255, 94/255, 0),
                        'salmon'))

# Caterpillar images
caterpillar = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')

#################################################################################
# Site-level analysis by latitude

catPresences = ccPlants %>% 
  filter(julianday >= 152, julianday <= 212, Group == 'caterpillar') %>% 
  group_by(ID, Name, Latitude, sciName, Family, plantOrigin) %>% 
  summarize(presence = ifelse(sum(Quantity) > 0, 1, 0))

catData = ccPlants %>% 
  filter(julianday >= 152, julianday <=212) %>% 
  distinct(ID, Name, Latitude, sciName, Family, plantOrigin) %>%
  left_join(catPresences)

catData$presence[is.na(catData$presence)] = 0

catDataBySiteAll = catData %>%
  group_by(Name, Latitude) %>%
  summarize(nSurvs = n_distinct(ID),
            nSurvsAlien = n_distinct(ID[plantOrigin == 'alien']),
            nSurvsNative = n_distinct(ID[plantOrigin == 'native']),
            nSurvsCatsAlien = n_distinct(ID[plantOrigin == 'alien' & presence]),
            nSurvsCatsNative = n_distinct(ID[plantOrigin == 'native' & presence]),
            pctCatAlien = 100*nSurvsCatsAlien/nSurvsAlien,
            pctCatNative = 100*nSurvsCatsNative/nSurvsNative,
            pctAlienSurveys = 100*nSurvsAlien/(nSurvsAlien + nSurvsNative))

# Based on this plot, sites above 45N latitude basically have no alien plant surveys, therefore they should be excluded from 
# native / alien comparison
par(mfrow = c(1,1))
plot(catDataBySiteAll$Latitude, catDataBySiteAll$pctAlienSurveys, cex = log10(catDataBySiteAll$nSurvs), 
     xlab = "Latitude", ylab = "% of surveys on alien plants", las = 1, cex.lab = 1.5)

catDataBySite = catDataBySiteAll %>%
  filter(nSurvsAlien >= 10, nSurvsNative >= 10, Latitude < 45)

catDataForAnalysis = catData %>%
  filter(Name %in% catDataBySite$Name) %>%
  mutate(scaledLatitude = scale(Latitude))

ccPlantsForAlienNativeComparison = ccPlants %>%
  filter(Latitude < 45)


cols = brewer.pal(6, "YlGnBu")
# Define colour pallete
pal = colorRampPalette(cols)
# Rank variable for colour assignment
catDataBySite$LatIndex = round(99*(catDataBySite$Latitude - min(catDataBySite$Latitude, na.rm = TRUE))/
                                 (max(catDataBySite$Latitude, na.rm = TRUE) - min(catDataBySite$Latitude))) + 1

par(mfrow = c(1,1), mar = c(5, 5, 1, 1), oma = c(0,0,0,0))
plot(catDataBySite$pctCatNative, catDataBySite$pctCatAlien, cex = log10(catDataBySite$nSurvsAlien),
     col = pal(100)[catDataBySite$LatIndex], pch = 16, xlab = "% of native surveys with caterpillars", 
     ylab = "% of alien surveys with caterpillars", cex.lab = 1.5)
abline(a=0, b= 1, col = 'red')
legend("topleft", legend = c(round(min(catDataBySite$Latitude)), round((min(catDataBySite$Latitude)+max(catDataBySite$Latitude))/2),
                             round(max(catDataBySite$Latitude))),
       pch = 16, cex = 1.5, col = pal(100)[c(1, 50, 100)], pt.cex = 2, title = 'Latitude')


# Logistic regression of caterpillar presence as predicted by latitude, plant origin, and their interaction
log.Origin.Latitude = glm(presence ~ plantOrigin + scaledLatitude + plantOrigin*scaledLatitude, 
                                 data = catDataForAnalysis, family = "binomial")

intplotLatOrigin = interact_plot(log.Origin.Latitude, pred = 'scaledLatitude', modx = 'plantOrigin', 
                        interval = TRUE, int.type = 'confidence', int.width = .95,
                        y.label = "Prop. of surveys with caterpillars",
                        legend.main = "Plant origin", line.thickness = 2, cex.lab = 1.5,
                        colors = c('red', 'gray50'))


# Logistic regression of caterpillar presence as predicted by latitude, plant origin, and their interaction,
# with site-level (Name) random effects.
log.Origin.Latitude.Name = glmer(presence ~ plantOrigin + scaledLatitude + plantOrigin*scaledLatitude + (1 | Name), 
                                 data = catDataForAnalysis, family = "binomial")

intplotOriginLatitudeName = interact_plot(log.Origin.Latitude.Name, pred = 'scaledLatitude', modx = 'plantOrigin', 
                        interval = TRUE, int.type = 'confidence', int.width = .95,
                        y.label = "Prop. of surveys with caterpillars",
                        legend.main = "Plant origin", line.thickness = 2, cex.lab = 1.5,
                        colors = c('red', 'gray50'))

plot1 = intplotOriginLatitudeName + 
  theme_bw() +
  annotation_raster(caterpillar, ymin = .041, ymax= .052,xmin = 39, xmax = 43.4) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 20))



# Logistic regression of caterpillar presence as predicted by latitude, plant origin, and their interaction,
# with site-level (Name) and plant species (sciName) random effects.
log.Origin.Latitude.Name.sciName = glmer(presence ~ plantOrigin + scaledLatitude + plantOrigin*scaledLatitude + 
                                           (1 | Name) + (1 | sciName), 
                                        data = catDataForAnalysis, family = "binomial")

intplotsciName = interact_plot(log.Origin.Latitude.Name.sciName, pred = 'scaledLatitude', modx = 'plantOrigin', 
                              interval = TRUE, int.type = 'confidence', int.width = .95,
                              y.label = "Prop. of surveys with caterpillars",
                              legend.main = "Plant origin", line.thickness = 2, cex.lab = 1.5,
                              colors = c('red', 'gray50'))
plot2 = intplotsciName + 
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 20))


# Figure displaying analyses with and without including plant sciName as a random effect
pdf('Figures/nativeStatus_by_latitude.pdf', height = 4, width = 8)
ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()


# Latitudinal range of alien plant species in our dataset
alienRanges = ccPlantsForAlienNativeComparison %>%
  filter(plantOrigin == 'alien') %>%
  group_by(sciName, Family) %>%
  summarize(nSurvs = n_distinct(ID),
            nBranches = n_distinct(PlantFK),
            minLat = min(Latitude),
            maxLat = max(Latitude)) %>%
  ungroup() %>%
  filter(nSurvs >= 10, nBranches >= 5) %>%
  left_join(treeFams) %>%
  arrange(desc(maxLat), desc(minLat))

alienRanges$famcolor[is.na(alienRanges$famcolor)] = 'gray50'

par(mfrow = c(1,1), mar = c(5, 12, 1, 1))
plot(range(c(alienRanges$minLat, alienRanges$maxLat)), c(0, nrow(alienRanges)), type = 'n', xlab = 'Latitudinal Range', yaxt = 'n', ylab = '')
for (i in 1:nrow(alienRanges)) {
  
  lines(c(alienRanges$minLat[i], alienRanges$maxLat[i]), c(i, i), col = alienRanges$famcolor[i], lwd = log(alienRanges$nSurvs[i]))
}
mtext(alienRanges$sciName, 2, at = 1:nrow(alienRanges), col = alienRanges$famcolor, las = 1, line = 1)


# Latitudinal variation within individual host plant species
hostplants = c('Fagus grandifolia', 'Acer negundo', 'Acer rubrum', 
               'Carpinus caroliniana', 'Cercis canadensis', 'Prunus serotina')

# Plotting % of surveys by half degree band

pdf('Figures/latitude_within_species.pdf', height = 6, width = 10)
par(mfrow = c(2, 3), mar = c(6, 6, 3, 1), mgp = c(3, 1, 0))
for (h in hostplants) {
  tmp = catDataForAnalysis %>%
    filter(sciName == h) %>%
    mutate(latBand = round(2*Latitude)/2) %>%  # half degree bands
    group_by(latBand) %>%
    summarize(nSurvs = n(),
              nSurvsWithCats = sum(presence == 1),
              pctCats = 100*nSurvsWithCats/nSurvs)
  
  plot(tmp$latBand, tmp$pctCats, pch = 16, cex = log10(tmp$nSurvs) + 0.5, xlab = 'Latitude band', 
       ylab = '% of surveys', main = paste0(h, " (", formatC(sum(tmp$nSurvs), big.mark = ","), ")"),
       las = 1, cex.lab = 2, cex.axis = 1.5, cex.main = 1.5, col = 'gray50')
  
  tmp.lm = lm(pctCats ~ latBand, data = tmp, weights = log10(nSurvs))
  abline(tmp.lm, lty = 'dashed', lwd = 2, xpd = FALSE)
  
}
dev.off()



# Site-specific alien-native comparisons for the sites with sufficient surveys of each group
nSurvsThreshold = 100

alienNativeSites = catDataForAnalysis %>% 
  count(Name, plantOrigin) %>% 
  filter(n > nSurvsThreshold) %>% 
  count(Name) %>%
  filter(n == 2) 


withinSiteComparisons = data.frame(Name = NULL, Family = NULL, Group = NULL, nAlienSurveys = NULL, 
                                 nNativeSurveys = NULL, nAlienBranches = NULL, nNativeBranches = NULL, 
                                 estimate = NULL, se = NULL, l95 = NULL, u95 = NULL, p = NULL)

for (n in alienNativeSites$Name) {
  for (a in arthropods$Group) {
    
    tmp = comparingNativeAlien(ccPlants[ccPlants$Name == n, ], arthGroup = a, plantFamily = 'All', 
                               jdRange = c(152, 212), # June + July
                               minArths = 0, minBranches = 2)
    
    tmpcomp = data.frame(Name = n,
                         Family = tmp$Family,
                         Group = tmp$Group,
                         nAlienSurveys = tmp$nAlienSurveys,
                         nNativeSurveys = tmp$nNativeSurveys,
                         nAlienBranches = tmp$nAlienBranches,
                         nNativeBranches = tmp$nNativeBranches,
                         nAlienSurvsWithArth = tmp$nAlienSurvsWithArth,
                         nNativeSurvsWithArth = tmp$nNativeSurvsWithArth,
                         propAlienSurvsWithArth = tmp$propAlienSurvsWithArth,
                         propNativeSurvsWithArth = tmp$propNativeSurvsWithArth,
                         errorAlienSurvsWithArth = tmp$errorAlienSurvsWithArth,
                         errorNativeSurvsWithArth = tmp$errorNativeSurvsWithArth,
                         propTestZ = tmp$propTestZ,
                         propTestP = tmp$propTestP,
                         estimate = tmp$model$coefficients$cond[2,1],
                         se = tmp$model$coefficients$cond[2,2],
                         l95 = tmp$confint[1],
                         u95 = tmp$confint[2],
                         p = tmp$model$coefficients$cond[2,4])
    
    if (is.nan(tmpcomp$se)) {
      tmpcomp$estimate = NA
    }
    
    withinSiteComparisons = rbind(withinSiteComparisons, tmpcomp)
    
  }
}

pdf('Figures/withinSite_native_alien_comparisons.pdf', height = 5, width = 8)
par(mfrow = c(2, 3), mar = c(5, 5, 1, .5), mgp = c(3, 1, 0), tck = -0.03)
for (a in arthropods$Group) {
  
  tmp.df = withinSiteComparisons[withinSiteComparisons$Group == a, ]
  
  maxProp = 1.2*100*max(c(tmp.df$propNativeSurvsWithArth, tmp.df$propAlienSurvsWithArth))
  
  plot(100*tmp.df$propNativeSurvsWithArth, 100*tmp.df$propAlienSurvsWithArth, 
       col = arthropods$color[arthropods$Group == a], pch = 16, las = 1,
       cex = log10(tmp.df$nAlienSurveys + tmp.df$nNativeSurveys), cex.lab = 1.5,
       xlab = '% of native surveys', cex.axis = 1.2,
       ylab = '% of alien surveys',
       ylim = c(0, maxProp),
       xlim = c(0, maxProp))
  
  abline(a=0, b = 1, xpd = FALSE)
  
  # 95% CI segments
  # adding 95% CI line segments
  segments(100*tmp.df$propNativeSurvsWithArth - 100*tmp.df$errorNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth,
           100*tmp.df$propNativeSurvsWithArth + 100*tmp.df$errorNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth,
           col = arthropods$color[arthropods$Group == a], lwd = 2)
  
  segments(100*tmp.df$propNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth - 100*tmp.df$errorAlienSurvsWithArth,
           100*tmp.df$propNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth + 100*tmp.df$errorAlienSurvsWithArth,
           col = arthropods$color[arthropods$Group == a], lwd = 2)

  # Add bug icon
  bug = readPNG(paste0('images/', a, '.png'))
  rasterImage(bug, 0.02*maxProp, .75*maxProp, .3*maxProp, maxProp)
  
}
dev.off()




# Native vs alien within-family comparisons
pdf('Figures/withinFamily_native_alien_comparisons.pdf', height = 5, width = 8)
par(mfrow = c(2, 3), mar = c(5, 5, 1, .5), mgp = c(3, 1, 0), tck = -0.03)
for (a in arthropods$Group) {
  
  tmp.df = comparisons[comparisons$Group == a & comparisons$Family != 'All', ]
  #tmp.df$col = ifelse(tmp.df$propTestP <= 0.01, arthropods$color[arthropods$Group == a],
  #                    ifelse(tmp.df$propTestP > 0.01 & tmp.df$propTestP < 0.05,
  #                           arthropods$color2[arthropods$Group == a], 'gray80'))
  tmp.df$col = ifelse(tmp.df$propTestP < 0.05, arthropods$color[arthropods$Group == a], 'gray80')
  tmp.df$pch = ifelse(tmp.df$propTestP < 0.01, 16, 1)
  
  maxProp = 1.2*100*max(c(tmp.df$propNativeSurvsWithArth, tmp.df$propAlienSurvsWithArth))
  
  plot(100*tmp.df$propNativeSurvsWithArth, 100*tmp.df$propAlienSurvsWithArth, 
       pch = tmp.df$pch, las = 1, 
       col = tmp.df$col, #arthropods$color[arthropods$Group == a], 
       cex = log10(tmp.df$nAlienSurveys + tmp.df$nNativeSurveys), cex.lab = 1.8,
       xlab = '% of native surveys', cex.axis = 1.3,
       ylab = '% of alien surveys',
       ylim = c(0, maxProp),
       xlim = c(0, maxProp))
  
  
  abline(a=0, b = 1, xpd = FALSE)
  
  # 95% CI segments
  # adding 95% CI line segments
  segments(100*tmp.df$propNativeSurvsWithArth - 100*tmp.df$errorNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth,
           100*tmp.df$propNativeSurvsWithArth + 100*tmp.df$errorNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth,
           col = tmp.df$col, #arthropods$color[arthropods$Group == a], 
           lwd = 2)
  
  segments(100*tmp.df$propNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth - 100*tmp.df$errorAlienSurvsWithArth,
           100*tmp.df$propNativeSurvsWithArth, 
           100*tmp.df$propAlienSurvsWithArth + 100*tmp.df$errorAlienSurvsWithArth,
           col = tmp.df$col, #arthropods$color[arthropods$Group == a], 
           lwd = 2)
  
  text(100*tmp.df$propNativeSurvsWithArth, 100*tmp.df$propAlienSurvsWithArth, 
       substr(tmp.df$Family, 1, 2), cex = 1.6)
  
  # Add bug icon
  bug = readPNG(paste0('images/', a, '.png'))
  rasterImage(bug, 0.02*maxProp, .75*maxProp, .3*maxProp, maxProp)
  
}
dev.off()




##########################################################################################################
# Figure 1. Comparison of % surveys with caterpillars across tree species

plantList = officialPlantList %>%
  distinct(sciName, rank, Family)


# Tree species with at least 50 surveys, ranked by % of surveys with caterpillars
byTreeSpp = AnalysisBySciName(ccPlants, ordersToInclude = 'caterpillar', 
                              jdRange = c(152, 212)) %>%   # June + July
  left_join(plantList, by = 'sciName') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  filter(rank == 'species',
         nSurveys >= 40,
         nBranches >= 5) %>%
  mutate(color = ifelse(plantOrigin == 'native', 'gray70', 'firebrick2')) %>%
  arrange(desc(fracSurveys)) %>%
  left_join(treeFams, by = 'Family')

byTreeSpp$famcolor[is.na(byTreeSpp$famcolor)] = 'gray50'

numspp = nrow(byTreeSpp)
if (numspp %% 2 != 0) { 
  numspp = numspp + 1 # add 1 to make numspp even if necessary
  byTreeSpp = rbind(byTreeSpp, NA)
}

pdf('Figures/Figure1_ranking_tree_spp_2col.pdf', height = 9, width = 12)
par(mar = c(5, 10, 2, 1), mgp = c(3, 1, 0), mfrow = c(1,2), oma = c(0, 0, 0, 0), xpd = NA)
plot(byTreeSpp$fracSurveys[1:(numspp/2)], (numspp/2):1, yaxt = 'n', ylab = '', xlab = '% of surveys',
     cex.axis = 1.5, cex.lab = 2, pch = 16, col = byTreeSpp$color[1:(numspp/2)], 
     xlim = c(0, 32), ylim = c(1, numspp/2),
     cex = 2*log10(byTreeSpp$nSurveys[1:(numspp/2)])/max(log10(byTreeSpp$nSurveys[1:(numspp/2)])),
     main = paste0("Species rank 1-", numspp/2), cex.main = 1.5)
segments(byTreeSpp$LL95frac[1:(numspp/2)], (numspp/2):1, byTreeSpp$UL95frac[1:(numspp/2)], (numspp/2):1,
         lwd = 2, col = byTreeSpp$color[1:(numspp/2)])
mtext(byTreeSpp$sciName[1:(numspp/2)], 2, at = (numspp/2):1 + .3, line = 1, adj = 1, las = 1, cex = .9)
points(rep(-2, numspp/2), (numspp/2):1, pch = 15, col = byTreeSpp$famcolor[1:(numspp/2)], cex = 1.4)

legend("bottomright", c("Fagaceae", "Betulaceae", "Sapindaceae", "Caprifoliaceae", "Juglandaceae", "Rosaceae", "Other"), 
       col = c(rgb(230/255, 159/255, 0),
               rgb(86/255, 180/255, 233/255),
               rgb(0, 158/255, 115/255),
               rgb(240/255, 228/255, 66/255),
               rgb(213/255, 94/255, 0),
               'salmon',
               'gray50'), 
       pch = 15, cex = 1.3, pt.cex = 2)

plot(byTreeSpp$fracSurveys[(numspp/2 + 1):numspp], (numspp/2):1, yaxt = 'n', ylab = '', xlab = '% of surveys',
     cex.axis = 1.5, cex.lab = 2, pch = 16, col = byTreeSpp$color[(numspp/2 + 1):numspp], 
     xlim = c(0, 32), ylim = c(1, numspp/2),
     cex = 2*log10(byTreeSpp$nSurveys[(numspp/2 + 1):numspp])/max(log10(byTreeSpp$nSurveys[(numspp/2 + 1):numspp]), na.rm = T),
     main = paste0("Species rank ", numspp/2 + 1, "-", numspp), cex.main = 1.5)
segments(byTreeSpp$LL95frac[(numspp/2 + 1):numspp], (numspp/2):1, byTreeSpp$UL95frac[(numspp/2 + 1):numspp], (numspp/2):1,
         lwd = 2, col = byTreeSpp$color[(numspp/2 + 1):numspp])
mtext(byTreeSpp$sciName[(numspp/2 + 1):numspp], 2, at = (numspp/2):1 + .3, line = 1, adj = 1, las = 1, cex = .9)
points(rep(-2, numspp/2), (numspp/2):1, pch = 15, col = byTreeSpp$famcolor[(numspp/2 + 1):numspp], cex = 1.4)    # 

legend("bottomright", c("native", "alien"), 
       col = c('gray70', 'firebrick2'), 
       pch = 16, cex = 1.4, pt.cex = 2, lwd = 2, lty = 'solid')

rasterImage(caterpillar, 16, 35, 32, 45)
dev.off()


# Supplemental Table S1

tableS1 = byTreeSpp %>%
  mutate(fracSurveys = round(fracSurveys, 2),
         LL95frac = round(LL95frac, 2),
         UL95frac = round(UL95frac, 2)) %>%
  select(sciName, nativeStatus, plantOrigin, nSurveys, nBranches, nSites, fracSurveys, LL95frac, UL95frac) %>%
  rename(`Plant species` = sciName,
         `USDA Native Status` = nativeStatus,
         Origin = plantOrigin,
         n_Surveys = nSurveys,
         n_Branches = nBranches,
         n_Sites = nSites,
         `Pct with caterpillars` = fracSurveys,
         `95%CI Lower limit` = LL95frac,
         `95%CI Upper limit` = UL95frac)

write.csv(tableS1, 'data/Table_S1.csv', row.names = F)



# Phylogenetic signal of fracSurveys across plant species
# install.packages('rtrees', repos=c(rtrees='https://daijiang.r-universe.dev', CRAN='https://cloud.r-project.org'))

library(rtrees)
library(ape)
library(phytools)

# Create a dataframe of species, genus, and family for generating a phylogeny
speciesList = speciesList = data.frame(species = byTreeSpp$sciName, 
                                       genus = word(byTreeSpp$sciName, 1), 
                                       family= byTreeSpp$Family, 
                                       fracSurveys = byTreeSpp$fracSurveys)

plantTree = rtrees::get_tree(speciesList[, 1:3], taxon = 'plant')
phyloSpeciesList = gsub("_", " ", plantTree$tip.label)

# Now get the tree species dataframe in the same phylogenetic order as the tip labels
speciesList$species = factor(speciesList$species, levels = phyloSpeciesList)
speciesList = speciesList[order(speciesList$species),]
speciesList$species = as.character(speciesList$species)

# Calculate Blomberg's K and Pagel's lambda
blomK = phylosig(plantTree, speciesList$fracSurveys, method = 'K', nsim = 999, test = TRUE)       # p = 0.62
lambda = phylosig(plantTree, speciesList$fracSurveys, method = 'lambda', nsim = 999, test = TRUE) # p = 1


#################################################################################################
# Figure 2. Analysis of caterpillar species richness by plant species based on iNat identifications

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-count-data/master/"

cc_data_repo <- "https://github.com/hurlbertlab/caterpillars-count-data"
cc_webpage <- read_html(cc_data_repo)
cc_repo_links <- html_attr(html_nodes(cc_webpage, "a"), "href")
cc_data_links <- tibble(link = cc_repo_links[grepl(".csv", cc_repo_links)]) %>%
  mutate(file_name = word(link, 6, 6, sep = "/")) %>%
  distinct()

expert = read.csv(paste(github_raw, filter(cc_data_links, grepl("ExpertIdentification.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

catSpecies = ccPlants %>% 
  select(ID, Name, Region, Year, LocalDate, julianday, Code, sciName, arthID, Group, Quantity) %>%
  right_join(expert, by = c('arthID' = 'ArthropodSightingFK')) %>%
  rename(ID = ID.x, expertID = ID.y) %>%
  filter(julianday >= 152, julianday <= 212) %>% # all of June and July
  group_by(sciName) %>%
  summarize(nCatTaxa = n_distinct(TaxonName[Group == 'caterpillar' & Rank != 'order']),
            nBeetleTaxa = n_distinct(TaxonName[Group == 'beetle' & Rank != 'order']),
            nSpiderTaxa = n_distinct(TaxonName[Group == 'spider' & Rank != 'order']),
            nHopperTaxa = n_distinct(TaxonName[Group == 'leafhopper' & Rank != 'order']),
            nTruebugTaxa = n_distinct(TaxonName[Group == 'truebugs' & Rank != 'order']),
            nAntTaxa = n_distinct(TaxonName[Group == 'ant' & Rank != 'order']),
            nSurveysWithPhotos = n_distinct(ID),
            nBranches = n_distinct(Code),
            Regions = paste(unique(Region), collapse = '-')) %>%
  arrange(desc(nCatTaxa)) %>%
  mutate(label = paste0(substr(sciName, 1, 1), ". ", word(sciName, 2))) %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  mutate(color = ifelse(plantOrigin == 'native', 'gray50', 'firebrick2'),
         native = ifelse(plantOrigin == 'native', 1, 0),
         logCatTaxa = log10(nCatTaxa),
         logPhotos = log10(nSurveysWithPhotos))



# Linear models
lm.native = lm(log10(catSpecies$nCatTaxa[catSpecies$plantOrigin == 'native' & catSpecies$nCatTaxa > 0]) ~ 
                 log10(catSpecies$nSurveysWithPhotos[catSpecies$plantOrigin == 'native' & catSpecies$nCatTaxa > 0]))

lm.alien = lm(log10(catSpecies$nCatTaxa[catSpecies$plantOrigin == 'alien' & catSpecies$nCatTaxa > 0]) ~ 
                log10(catSpecies$nSurveysWithPhotos[catSpecies$plantOrigin == 'alien' & catSpecies$nCatTaxa > 0]))

catSpeciesWithTaxa = catSpecies[catSpecies$nCatTaxa > 0, ]

# slope difference, p = 0.014, beta = 0.40
lm.native.alien.cat = lm(logCatTaxa ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                           data = catSpeciesWithTaxa)

# slope difference, p = 9.6e-5, beta = 0.38
lm.native.alien.beetle = lm(log10(nBeetleTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                     data = catSpecies[catSpecies$nBeetleTaxa > 0, ])

# slope difference, p = 5.0e-5, beta = 0.38
lm.native.alien.spider = lm(log10(nSpiderTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = catSpecies[catSpecies$nSpiderTaxa > 0, ])

# slope difference, p = 0.027, beta = 0.24
lm.native.alien.hopper = lm(log10(nHopperTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = catSpecies[catSpecies$nHopperTaxa > 0, ])

# slope difference, p = 0.49, beta = 0.09
lm.native.alien.truebug = lm(log10(nTruebugTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = catSpecies[catSpecies$nTruebugTaxa > 0, ])

# slope difference, p = 0.048, beta = 0.25
lm.native.alien.ant = lm(log10(nAntTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                             data = catSpecies[catSpecies$nAntTaxa > 0, ])


taxaIntPlotCat = interact_plot(lm.native.alien.cat, pred = 'logPhotos', modx = 'plantOrigin', 
                           interval = TRUE, int.type = 'confidence', int.width = .95,
                           y.label = expression(log[10]~~species~richness),
                           x.label = expression(log[10]~~"#"~surveys~with~photos),
                           legend.main = "Plant origin", line.thickness = 2, cex.lab = 1.5,
                           colors = c('red', 'gray50'), plot.points = TRUE, #vary.lty = FALSE,
                           point.shape = TRUE, point.size = 2.5, show.legend = FALSE)

pdf('Figures/Figure2_caterpillar_taxa.pdf', height = 5, width = 8)
taxaIntPlotCat + 
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  annotation_raster(caterpillar, ymin = 1,ymax= 1.4,xmin = 0,xmax = .8)
dev.off()


catfig = interact_plot(lm.native.alien.cat, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
           theme_bw() +
           theme(axis.title = element_text(size = 14),
           axis.text = element_text(size = 12),
           axis.title.x = element_text(margin = margin(t = 6)), 
           axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
           theme(legend.position = "none") +
  annotation_raster(caterpillar, ymin = 1,ymax= 1.4,xmin = 0,xmax = .8)

beetfig = interact_plot(lm.native.alien.beetle, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  theme(legend.position = "none") +
  annotation_raster(beetleImage, ymin = 1.25, ymax= 1.7,xmin = 0,xmax = .8)


hopfig = interact_plot(lm.native.alien.hopper, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  theme(legend.position = "none") +
  annotation_raster(hopperImage, ymin = 1,ymax= 1.5,xmin = 0,xmax = .8)


spiderfig = interact_plot(lm.native.alien.spider, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  theme(legend.position = "none") +
  annotation_raster(spiderImage, ymin = 1,ymax= 1.6,xmin = 0,xmax = 1)


bugfig = interact_plot(lm.native.alien.truebug, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  theme(legend.position = "none") +
  annotation_raster(truebugImage, ymin = 1,ymax= 1.5,xmin = 0,xmax = .9)


antfig = interact_plot(lm.native.alien.ant, pred = 'logPhotos', modx = 'plantOrigin', 
                       interval = TRUE, int.type = 'confidence', int.width = .95,
                       y.label = expression(log[10]~~species~richness),
                       x.label = expression(log[10]~~surveys~with~photos),
                       line.thickness = 2, cex.lab = 1.5,
                       colors = c('red', 'gray50'), plot.points = TRUE,
                       point.shape = TRUE, point.size = 2.5) + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2)) +
  theme(legend.position = "none") +
  annotation_raster(antImage, ymin = 1,ymax= 1.4,xmin = 0,xmax = .9)

pdf('Figures/Figure4_speciesrichness.pdf', height = 6, width = 10)
ggarrange(catfig, spiderfig, hopfig, beetfig, bugfig, antfig, nrow = 2, ncol = 3,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'))
dev.off()


# Old basic plot
#pdf('Figures/Figure2_caterpillar_taxa.pdf', height = 6, width = 8)
par(mgp = c(4, 1, 0), tck = -0.03, mar = c(5, 7, 1, 1), mfrow = c(1,1), xpd = FALSE)
plot(log10(catSpecies$nSurveysWithPhotos), log10(catSpecies$nCatTaxa), pch = ifelse(catSpecies$plantOrigin == 'native', 16, 17), 
     xlab = expression(log[10] ~ "#" ~ surveys ~ with ~ photos), 
     ylab = expression(log[10] ~ "#" ~ caterpillar ~ taxa), las = 1, 
     cex.axis = 1.5, cex.lab = 2, col = catSpecies$color, cex = ifelse(catSpecies$plantOrigin == 'native', 1.3, 1.5))

abline(lm.native, lwd = 2, col = 'gray30')
abline(lm.alien, col = 'firebrick2', lwd = 2)

legend("topleft", c("native", "alien"), pch = c(16, 17), col = c('gray50', 'firebrick2'), cex = 1.5)
#dev.off()


# Number of taxa of other arthropod groups

SpeciesWithTaxa = catSpecies[catSpecies$nCatTaxa > 0, ]

lm.native.alien = lm(logCatTaxa ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                     data = catSpeciesWithTaxa)


###############################################################################
# Find set of plant families with sufficient data for native-alien comparisons

jdRange = c(152, 212) # June + July

familyStats = ccPlants %>%
  filter(julianday >= jdRange[1], 
         julianday <= jdRange[2],
         !WetLeaves) %>%
  distinct(ID, PlantFK, Family, sciName, plantOrigin, ObservationMethod) %>%
  group_by(Family) %>%
  summarize(nSpeciesA = n_distinct(sciName[plantOrigin == 'alien']),
            nSpeciesN = n_distinct(sciName[plantOrigin == 'native']),
            nBranchesA = n_distinct(PlantFK[plantOrigin == 'alien']),
            nBranchesN = n_distinct(PlantFK[plantOrigin == 'native']),
            nSurveysA = n_distinct(ID[plantOrigin == 'alien']),
            nSurveysN = n_distinct(ID[plantOrigin == 'native'])) %>%
  #arrange(desc(nBranchesA)) %>%
  filter(nBranchesA >= 5,
         nBranchesN >= 5,
         nSurveysA >= 50,
         nSurveysN >= 50)



comparisons = data.frame(Family = NULL, Group = NULL, nAlienSurveys = NULL, nNativeSurveys = NULL,
                         nAlienBranches = NULL, nNativeBranches = NULL, estimate = NULL, se = NULL, 
                         l95 = NULL, u95 = NULL, p = NULL)

for (f in c('All', familyStats$Family)) {
  for (a in arthropods$Group) {
    
    tmp = comparingNativeAlien(ccPlants, arthGroup = a, plantFamily = f, 
                               jdRange = c(152, 212), # June + July
                               minArths = 5)
    
    tmpcomp = data.frame(Family = tmp$Family,
                         Group = tmp$Group,
                         nAlienSurveys = tmp$nAlienSurveys,
                         nNativeSurveys = tmp$nNativeSurveys,
                         nAlienBranches = tmp$nAlienBranches,
                         nNativeBranches = tmp$nNativeBranches,
                         nAlienSurvsWithArth = tmp$nAlienSurvsWithArth,
                         nNativeSurvsWithArth = tmp$nNativeSurvsWithArth,
                         propAlienSurvsWithArth = tmp$propAlienSurvsWithArth,
                         propNativeSurvsWithArth = tmp$propNativeSurvsWithArth,
                         errorAlienSurvsWithArth = tmp$errorAlienSurvsWithArth,
                         errorNativeSurvsWithArth = tmp$errorNativeSurvsWithArth,
                         propTestZ = tmp$propTestZ,
                         propTestP = tmp$propTestP,
                         estimate = tmp$model$coefficients$cond[2,1],
                         se = tmp$model$coefficients$cond[2,2],
                         l95 = tmp$confint[1],
                         u95 = tmp$confint[2],
                         p = tmp$model$coefficients$cond[2,4])
    
    if (is.nan(tmpcomp$se)) {
      tmpcomp$estimate = NA
    }
    
    comparisons = rbind(comparisons, tmpcomp)
    
  }
}

# Add text symbols for reflecting p-values
comparisons = comparisons %>%
  mutate(pText = case_when(propTestP <= 0.001 & propTestZ >= 0 ~ '+++',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ >= 0 ~ '++',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ >= 0 ~ '+',
                           propTestP <= 0.001 & propTestZ < 0 ~ '---',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ < 0 ~ '--',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ < 0 ~ '-',
                           propTestP > 0.05 ~ '')
  )

###############################################################################################
# Figure 3. Comparing % of surveys with arthropod (+- 95% CI) for all plant families and arthropod groups
all = comparisons %>%
  left_join(arthropods, by = 'Group') %>%
  filter(Family == 'All')


pdf('Figures/alien_vs_native.pdf', height = 6, width = 8)
par(mfrow = c(1,1), mgp = c(3, 1, 0), oma = c(0,0,0,0), mar = c(5, 5, 0, 0))
plot(100*all$propNativeSurvsWithArth, 100*all$propAlienSurvsWithArth, col = all$color, cex = 2.5, 
     xlab = "% of native surveys", ylab = "% of alien surveys", cex.lab = 2, cex.axis = 1.75, las = 1, pch = 16, 
     xlim = c(0, 25), ylim = c(0, 25))
abline(a=0, b=1, lty = 'dotted', lwd = 2, col = 'gray30')

# adding 95% CI line segments
segments(100*all$propNativeSurvsWithArth - 100*all$errorNativeSurvsWithArth, 100*all$propAlienSurvsWithArth,
         100*all$propNativeSurvsWithArth + 100*all$errorNativeSurvsWithArth, 100*all$propAlienSurvsWithArth,
         col = all$color, lwd = 3)

segments(100*all$propNativeSurvsWithArth, 100*all$propAlienSurvsWithArth - 100*all$errorAlienSurvsWithArth,
         100*all$propNativeSurvsWithArth, 100*all$propAlienSurvsWithArth + 100*all$errorAlienSurvsWithArth,
         col = all$color, lwd = 3)

# asterisks
asteriskOffset = 1.5
text(100*all$propNativeSurvsWithArth[all$Group %in% c('caterpillar', 'spider', 'beetle')], 
     100*all$propAlienSurvsWithArth[all$Group %in% c('caterpillar', 'spider', 'beetle')] - asteriskOffset,
     '***', cex = 3)
text(100*all$propNativeSurvsWithArth[all$Group == 'leafhopper'], 
     100*all$propAlienSurvsWithArth[all$Group == 'leafhopper'] - asteriskOffset,
     '*', cex = 3)

# bug icons
for(a in arthropods$Group) {
  bug = readPNG(paste0('images/', a, '.png'))
  
  xoffset = ifelse(a == 'beetle', 1, 1.5)
  yheight = ifelse(a %in%  c('beetle', 'caterpillar'), 2, 3)

  rasterImage(bug, 100*all$propNativeSurvsWithArth[all$Group == a] - xoffset, 
              100*all$propAlienSurvsWithArth[all$Group == a] + .8, 
              100*all$propNativeSurvsWithArth[all$Group == a] + xoffset,
              100*all$propAlienSurvsWithArth[all$Group == a] + .8 + yheight)
}

legend("bottomright", legend = c(" *   p < 0.02", "*** p < 0.00001"), cex = 1.5)
dev.off()



pdf('Figures/Figure3_plant_family_comparison.pdf', height = 6, width = 10)
par(mar = c(2, 0, 3, 0), oma = c(3, 24, 0, 1), mfrow = c(1,6), mgp = c(3, .5, 0))

vertOffset = 0.1

# Order of plant families
famOrder = c(1, length(unique(comparisons$Family)):2) # alphabetical, but with "All" at the bottom
#famOrder = 1:length(unique(comparisons$Family)) 

for (a in arthropods$Group) {
  plot(100*comparisons$propNativeSurvsWithArth[comparisons$Group == a], 
       famOrder + vertOffset,
       xlab = "", ylab = "", yaxt = "n", tck = -0.03, xlim = c(0, ifelse(a == 'caterpillar', 22, 44)), 
       ylim = c(1, nrow(familyStats) + 3),
       pch = 16, col = arthropods$color[arthropods$Group == a], cex = 1.8, cex.axis = 1,
       main = arthropods$GroupLabel[arthropods$Group == a], cex.main = 1.7)
  segments(100*comparisons$propNativeSurvsWithArth[comparisons$Group == a] - 
             100*comparisons$errorNativeSurvsWithArth[comparisons$Group == a], 
           famOrder + vertOffset,
           100*comparisons$propNativeSurvsWithArth[comparisons$Group == a] + 
             100*comparisons$errorNativeSurvsWithArth[comparisons$Group == a], 
           famOrder + vertOffset, 
           col = arthropods$color[arthropods$Group == a])
  
  points(100*comparisons$propAlienSurvsWithArth[comparisons$Group == a], 
         famOrder - vertOffset,
         pch = 1, col = arthropods$color[arthropods$Group == a], cex = 1.8)
  segments(100*comparisons$propAlienSurvsWithArth[comparisons$Group == a] - 
             100*comparisons$errorAlienSurvsWithArth[comparisons$Group == a], 
           famOrder - vertOffset,
           100*comparisons$propAlienSurvsWithArth[comparisons$Group == a] + 
             100*comparisons$errorAlienSurvsWithArth[comparisons$Group == a], 
           famOrder - vertOffset, 
           col = arthropods$color[arthropods$Group == a])
  
  abline(h = 1.5, lwd = 2)
  
  text(# Commented out line below makes horizontal placement relative to largest values
    #100*max(c(comparisons$propAlienSurvsWithArth[comparisons$Group == a], 
    #       comparisons$propNativeSurvsWithArth[comparisons$Group == a])) + 2,
    ifelse(a == 'caterpillar', 20, 40),
    famOrder,
    labels = comparisons$pText[comparisons$Group == a], 
    # - symbols are much smaller than +, so making them a larger font size (2 vs 1.5)
    cex = ifelse(comparisons$propTestZ[comparisons$Group == a] < 0, 2, 1.5))
  
  # Add bug icon
  bug = readPNG(paste0('images/', a, '.png'))
  rasterImage(bug, 3, nrow(familyStats) + 1.5, ifelse(a == 'caterpillar', 20, 40), nrow(familyStats) + 3)
  
  # Put plant Family labels along the y-axis for the first plot
  if (a == arthropods$Group[1]) {
    mtext(paste0(comparisons$Family[comparisons$Group == a], 
                 " (", comparisons$nNativeSurveys[comparisons$Group == a], ", ", 
                 comparisons$nAlienSurveys[comparisons$Group == a], ")"),
          2, at = famOrder, 
          las = 1, line = 1, cex = 1.5)
    legend("topleft", inset=c(-1.4,0), legend=c("native","alien"), pch=c(16,1), lty = 'solid', 
           xpd = NA, cex = 2)
  }
}
mtext("% of surveys", 1, outer = TRUE, line = 1.5, cex = 2)

dev.off()



















####################################################################################

# OLD PLOTS

############################################################################################################
# Figure 4. Controlling for geographic variation by focusing on a single family, Sapindaceae, with the most surveys

ccPlants %>% 
  filter(Family=="Sapindaceae", 
         grepl(" ", sciName),    #only return names with spaces, i.e. full scientific names, not just genera
         julianday >= 152, julianday <= 212) %>% # June + July
  distinct(ID, sciName, plantRank, plantOrigin, Region) %>% 
  count(sciName, plantOrigin, plantRank, Region) %>% 
  arrange(Region, desc(n))

# This shows that the following regions have sufficient surveys on Aceraceae (>50 per plantOrigin) for a comparison:
aceraceaeRegions = list('ON',
               c('MA', 'CT','RI'),
               c('DC', 'MD', 'VA'),
               c('NC', 'SC'))

# Similarly, for Rosaceae (although it turns out several of these regions are only represented by 4-8 unique branches):
rosaceaeRegions = list('PA',
                       c('MA', 'CT', 'RI'),
                       c('DC', 'MD', 'VA'),
                       c('NC', 'SC'))

acerComparisons = nativeAlienAcrossRegions(ccPlants, "Sapindaceae", regions = aceraceaeRegions, minSurveys = 50, minBranches = 10) %>%
  mutate(pText = case_when(propTestP <= 0.001 & propTestZ >= 0 ~ '+++',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ >= 0 ~ '++',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ >= 0 ~ '+',
                           propTestP <= 0.001 & propTestZ < 0 ~ '---',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ < 0 ~ '--',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ < 0 ~ '-',
                           propTestP > 0.05 ~ ''))

# Plot - vertical - caterpillars on Aceraceae
pdf('Figures/Figure4_Aceraceae_across_regions.pdf', height = 6, width = 8)
par(mar = c(7, 6, 1, 1), mfrow = c(1,1), mgp = c(3, 1.2, 0), oma = c(0,0,0,0), xpd = FALSE)

horizOffset = 0.05

# Native
plot(4:1 + horizOffset, 100*acerComparisons$propNativeSurvsWithArth, 
     xlab = "", ylab = "% of surveys", xaxt = "n", tck = -0.02, xlim = c(0.5, 4.5), ylim = c(0, 18),
     pch = 16, col = 'gray50', cex = 2, cex.axis = 1.5, las = 1, cex.lab = 2)
segments(4:1 + horizOffset, 100*acerComparisons$propNativeSurvsWithArth - 
           100*acerComparisons$errorNativeSurvsWithArth, 
         4:1 + horizOffset,
         100*acerComparisons$propNativeSurvsWithArth + 
           100*acerComparisons$errorNativeSurvsWithArth, 
         lwd = 2, col = 'gray50')

# Alien
points(4:1 - horizOffset, 100*acerComparisons$propAlienSurvsWithArth, pch = 16, cex = 2, col = 'red')
segments(4:1 - horizOffset, 100*acerComparisons$propAlienSurvsWithArth - 
           100*acerComparisons$errorAlienSurvsWithArth, 
         4:1 - horizOffset,
         100*acerComparisons$propAlienSurvsWithArth + 
           100*acerComparisons$errorAlienSurvsWithArth, 
         lwd = 2, col = 'red')

#abline(h = 1.5, lwd = 2)

text(4:1, 13, labels = acerComparisons$pText, 
     # - symbols are much smaller than +, so making them a larger font size (2 vs 1.5)
     cex = ifelse(acerComparisons$propTestZ < 0, 2, 1.5))

mtext(paste0(acerComparisons$Region, 
             "\n(", acerComparisons$nNativeSurveys, ", ", 
             acerComparisons$nAlienSurveys, ")"),
      1, at = 4:1, las = 1, line = 2, cex = 1.5)

mtext("Region", 1, line = 5, cex = 2)

legend("topleft", legend=c("native","alien"), pch=c(16,16), lty = 'solid', col = c('gray50', 'red'),
       xpd = NA, cex = 1.5)

maple = readPNG("images/maple.png")
rasterImage(maple, 3.9, 14, 4.6, 19)
rasterImage(caterpillar, 3.3, 15, 4.1, 18)

dev.off()


# Supplemental table on the plant species examined, organized by family and plant origin

survsByFamilyOrigin = ccPlants %>% 
  group_by(Family, plantOrigin) %>% 
  summarize(nSurvs = n_distinct(ID))

familyComparisonSpecies = ccPlants %>% 
  filter(Family %in% familyStats$Family) %>%
  group_by(sciName, Family, plantOrigin) %>% 
  summarize(n = n_distinct(ID)) %>% 
  arrange(Family, plantOrigin, desc(n)) %>%
  left_join(survsByFamilyOrigin, by = c('Family', 'plantOrigin')) %>%
  mutate(pct = round(100*n/nSurvs,2))
