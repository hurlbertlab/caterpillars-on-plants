# Script for calculating average abundance, biomass, and occurrence per plant species, and for comparing broadly between native and alien plants

library(tidyverse)
library(png)
library(glmmTMB)
library(interactions)
library(ggpubr)

# Load analysis functions
source('code/plant_analysis_functions.r')

# Read in latest CC fullDataset through 2024
cc = read.csv('data/fullDataset_2025-04-18.csv', header = T, quote = '\"')

# Read in plant files (up-to-date versions in caterpillars-count-data repo)
officialPlantList = read.csv('data/officialPlantList2024-10-04.csv')
inferredPlantNames = read.csv('data/inferredPlantNames_2024-08-20.csv')
plantOrigin = read.csv('data/plant_origin_status.csv') %>%
  select(scientificName, nativeStatus, plantOrigin)

# Expert identifications of photo observations from iNaturalist
expert = read.csv('data/2025-06-17_ExpertIdentification.csv')

# The dataset for which we have plant species names with NameConfidence >= 2
# NOTE: sciName is the field that includes inferred scientific names in addition to official ones.
# Exclude data from Coweeta sites, where non-caterpillars were not recorded

jdRange = c(152,212) # June and July
minLongitude = -100  # restrict data to eastern North America

ccPlants = cc %>%
  left_join(inferredPlantNames[, c('PlantFK', 'InferredSciName', 'NameConfidence')], by = 'PlantFK') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  mutate(sciName = ifelse(Species == "N/A" & NameConfidence >= 2, InferredSciName, sciName)) %>%
  filter(!is.na(sciName),
         !Name %in% c('Coweeta - BB', 'Coweeta - BS', 'Coweeta - RK'),
         julianday >= jdRange[1], 
         julianday <= jdRange[2],
         Longitude >= minLongitude,
         !WetLeaves)

# Some summary statistics for the dataset used in analysis
nSurvs = length(unique(ccPlants$ID)) # 68741
nSites = length(unique(ccPlants$Name)) # 212
range(ccPlants$Latitude) # 32.33, 47.78
range(ccPlants$Year) # 2010, 2024
nBranches = length(unique(ccPlants$Code)) # 5438
nPlantSpecies = length(unique(ccPlants$sciName)) # 363
nNativePlantSpecies = length(unique(ccPlants$sciName[ccPlants$plantOrigin == 'native'])) # 254
nAlienPlantSpecies = length(unique(ccPlants$sciName[ccPlants$plantOrigin == 'alien'])) # 107
ccPlants %>% 
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
                      famcolor = rgb(230/255, 159/255, 0),
                      rgb(86/255, 180/255, 233/255),
                      rgb(0, 158/255, 115/255),
                      rgb(240/255, 228/255, 66/255),
                      rgb(213/255, 94/255, 0),
                      'salmon')

# Arthropod images
caterpillar = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')


#########################################################################################
# Figure 1. Overall arthropod occurrence on native vs alien plants for 6 arthropod groups

allSpeciesComparison = data.frame(Family = NULL, Group = NULL, nAlienSurveys = NULL, nNativeSurveys = NULL,
                         nAlienBranches = NULL, nNativeBranches = NULL, estimate = NULL, se = NULL, 
                         l95 = NULL, u95 = NULL, p = NULL)

for (a in arthropods$Group) {
  
  tmp = comparingNativeAlien(ccPlants, arthGroup = a, plantFamily = 'All', 
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
                                  p = tmp$model$coefficients$cond[2,4]) %>%
    mutate(estimate = ifelse(is.nan(se), NA, estimate))

  allSpeciesComparison = rbind(allSpeciesComparison, tmpcomp)

} #end arthropod loop

# Join in group-specific colors
allSpeciesComp = allSpeciesComparison %>%
  left_join(arthropods, by = 'Group') 


pdf('Figures/Figure1_alien_vs_native_occurrence.pdf', height = 7, width = 8)
par(mfrow = c(1,1), mgp = c(3, 1, 0), oma = c(0,0,0,0), mar = c(5, 5, 1, 1))
plot(100*allSpeciesComp$propNativeSurvsWithArth, 100*allSpeciesComp$propAlienSurvsWithArth, col = allSpeciesComp$color, cex = 2.5, 
     xlab = "% of native surveys", ylab = "% of alien surveys", cex.lab = 2, cex.axis = 1.75, las = 1, pch = 16, 
     xlim = c(0, 23), ylim = c(0, 23))
abline(a=0, b=1, lty = 'dotted', lwd = 2, col = 'gray30')

# adding 95% CI line segments
segments(100*allSpeciesComp$propNativeSurvsWithArth - 100*allSpeciesComp$errorNativeSurvsWithArth, 100*allSpeciesComp$propAlienSurvsWithArth,
         100*allSpeciesComp$propNativeSurvsWithArth + 100*allSpeciesComp$errorNativeSurvsWithArth, 100*allSpeciesComp$propAlienSurvsWithArth,
         col = allSpeciesComp$color, lwd = 3)

segments(100*allSpeciesComp$propNativeSurvsWithArth, 100*allSpeciesComp$propAlienSurvsWithArth - 100*allSpeciesComp$errorAlienSurvsWithArth,
         100*allSpeciesComp$propNativeSurvsWithArth, 100*allSpeciesComp$propAlienSurvsWithArth + 100*allSpeciesComp$errorAlienSurvsWithArth,
         col = allSpeciesComp$color, lwd = 3)

# asterisks
asteriskOffset = 1.5
text(100*allSpeciesComp$propNativeSurvsWithArth[allSpeciesComp$Group %in% c('caterpillar', 'spider', 'beetle')], 
     100*allSpeciesComp$propAlienSurvsWithArth[allSpeciesComp$Group %in% c('caterpillar', 'spider', 'beetle')] - asteriskOffset,
     '***', cex = 3)
text(100*allSpeciesComp$propNativeSurvsWithArth[allSpeciesComp$Group == 'leafhopper'], 
     100*allSpeciesComp$propAlienSurvsWithArth[allSpeciesComp$Group == 'leafhopper'] - asteriskOffset,
     '*', cex = 3)

# bug icons
for(a in arthropods$Group) {
  bug = readPNG(paste0('images/', a, '.png'))
  
  xoffset = ifelse(a == 'beetle', 1, 1.5)
  yheight = ifelse(a %in%  c('beetle', 'caterpillar'), 2, 3)
  
  rasterImage(bug, 100*allSpeciesComp$propNativeSurvsWithArth[allSpeciesComp$Group == a] - xoffset, 
              100*allSpeciesComp$propAlienSurvsWithArth[allSpeciesComp$Group == a] + .8, 
              100*allSpeciesComp$propNativeSurvsWithArth[allSpeciesComp$Group == a] + xoffset,
              100*allSpeciesComp$propAlienSurvsWithArth[allSpeciesComp$Group == a] + .8 + yheight)
}

legend("bottomright", legend = c(" *   p < 0.02", "*** p < 0.0001"), cex = 1.5)
dev.off()



###############################################################################
# Figure 2. Arthropod species richness on native vs alien plants 
#           for 6 arthropod groups, controlling for number of photos

arthSpecies = ccPlants %>% 
  select(ID, Name, Region, Year, LocalDate, julianday, Code, sciName, arthID, Group, Quantity) %>%
  inner_join(expert, by = c('arthID' = 'ArthropodSightingFK')) %>%
  rename(ID = ID.x, expertID = ID.y) %>%
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
         logPhotos = log10(nSurveysWithPhotos))


# Linear models

# slope difference, p = 0.038, beta = 0.35
lm.native.alien.cat = lm(log10(nCatTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                         data = arthSpecies[arthSpecies$nCatTaxa > 0, ])

# slope difference, p = 6.6e-5, beta = 0.40
lm.native.alien.beetle = lm(log10(nBeetleTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = arthSpecies[arthSpecies$nBeetleTaxa > 0, ])

# slope difference, p = 0.001, beta = 0.33
lm.native.alien.spider = lm(log10(nSpiderTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = arthSpecies[arthSpecies$nSpiderTaxa > 0, ])

# slope difference, p = 0.015, beta = 0.28
lm.native.alien.hopper = lm(log10(nHopperTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                            data = arthSpecies[arthSpecies$nHopperTaxa > 0, ])

# slope difference, p = 0.46, beta = 0.10
lm.native.alien.truebug = lm(log10(nTruebugTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                             data = arthSpecies[arthSpecies$nTruebugTaxa > 0, ])

# slope difference, p = 0.08, beta = 0.22
lm.native.alien.ant = lm(log10(nAntTaxa) ~ logPhotos + plantOrigin + plantOrigin*logPhotos, 
                         data = arthSpecies[arthSpecies$nAntTaxa > 0, ])

# Interaction figure panels
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

pdf('Figures/Figure2_speciesrichness.pdf', height = 6, width = 10)
ggarrange(catfig, spiderfig, hopfig, beetfig, bugfig, antfig, nrow = 2, ncol = 3,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'))
dev.off()










familyStats = ccPlants %>%
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
