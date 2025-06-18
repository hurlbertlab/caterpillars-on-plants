# Script for calculating average abundance, biomass, and occurrence per plant species, and for comparing broadly between native and alien plants

library(tidyverse)
library(png)
library(glmmTMB)
library(interactions)
library(ggpubr)
library(rtrees) #install.packages('rtrees', repos=c(rtrees='https://daijiang.r-universe.dev', CRAN='https://cloud.r-project.org'))
library(ape)
library(phytools)

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
                      famcolor = c(rgb(230/255, 159/255, 0),
                      rgb(86/255, 180/255, 233/255),
                      rgb(0, 158/255, 115/255),
                      rgb(240/255, 228/255, 66/255),
                      rgb(213/255, 94/255, 0),
                      'salmon'))

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





###############################################################################
# Figure 3. Arthropod occurrence on native vs alien plants for 6 arthropod groups
#           BY PLANT FAMILY

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

for (f in familyStats$Family) {
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


# Native vs alien within-family comparisons
pdf('Figures/Figure3_withinFamily_native_alien_occurrence.pdf', height = 5, width = 8)
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




#####################################################################################
# Figure 4. Comparison of % surveys with caterpillars across tree species

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

pdf('Figures/Figure4_ranking_tree_spp_2col.pdf', height = 9, width = 12)
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


# Phylogenetic signal of fracSurveys across plant species

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
blomK = phylosig(plantTree, speciesList$fracSurveys, method = 'K', nsim = 999, test = TRUE)       # K = 0.04, p = 0.77
lambda = phylosig(plantTree, speciesList$fracSurveys, method = 'lambda', nsim = 999, test = TRUE) # lambda = 7e-5, p = 1

