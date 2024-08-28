# Script for calculating average abundance, biomass, and occurrence per plant species, and for comparing broadly between native and alien plants

library(tidyverse)
library(RCurl)
library(rvest)
library(xml2)
library(glmmTMB)
library(png)

# load functions
source('code/plant_analysis_functions.r')

# Read in latest CC fullDataset
fd_data_repo <- "https://github.com/hurlbertlab/caterpillars-analysis-public/blob/master/data"
fd_webpage <- read_html(fd_data_repo)
fd_repo_links <- html_attr(html_nodes(fd_webpage, "a"), "href")
fd_data_links <- tibble(link = fd_repo_links[grepl("fullDataset", fd_repo_links)]) %>%
  mutate(file_name = word(link, 7, 7, sep = "/")) %>%
  distinct()

mostRecentFullDataset = fd_data_links$file_name[nrow(fd_data_links)]

cc = read.csv(paste0("https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/", mostRecentFullDataset), header = T, quote = '\"', fill = TRUE)


# Loading plant files from data_repo
data_repo = "https://raw.githubusercontent.com/hurlbertlab/caterpillars-count-data/master/plantSpecies/"

plants_webpage <- read_html("https://github.com/hurlbertlab/caterpillars-count-data/tree/master/plantSpecies")
plants_repo_links <- html_attr(html_nodes(plants_webpage, "a"), "href")
official_data_links <- tibble(link = plants_repo_links[grepl("officialPlantList", plants_repo_links)]) %>%
  mutate(file_name = word(link, 7, 7, sep = "/")) %>%
  distinct()

inferred_data_links <- tibble(link = plants_repo_links[grepl("inferredPlantNames", plants_repo_links)]) %>%
  mutate(file_name = word(link, 7, 7, sep = "/")) %>%
  distinct()

officialPlantList = read.csv(paste0(data_repo, official_data_links$file_name[nrow(official_data_links)]))
inferredPlantNames = read.csv(paste0(data_repo, inferred_data_links$file_name[nrow(inferred_data_links)]))

plantOrigin = read.csv(paste0(data_repo, "plant_origin_status.csv"))

# The dataset for which we have plant species names with NameConfidence >= 2
# NOTE: sciName is the field that includes inferred scientific names in addition to official ones
ccPlants = cc %>%
  left_join(inferredPlantNames[, c('PlantFK', 'InferredSciName', 'NameConfidence')], by = 'PlantFK') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  mutate(sciName = ifelse(Species == "N/A" & NameConfidence >= 2, InferredSciName, sciName)) %>%
  filter(!is.na(sciName))


# Comparison across tree species












# Find set of plant families with sufficient data for comparisons
jdRange = c(152, 194) # 3 weeks before to 3 weeks after summer solstice (173)

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


arthropods = data.frame(Group = c('caterpillar', 'spider', 'leafhopper', 'beetle', 'truebugs', 'ant'),
                        GroupLabel = c('caterpillars', 'spiders', 'hoppers', 'beetles', 'true bugs', 'ants'),
                        color = c('limegreen', 'gray50', 'dodgerblue', 'salmon', 'magenta3', 'orange'),
                        color2 = c('darkgreen', 'black', 'darkblue', 'red', 'purple3', 'orange4'))


comparisons = data.frame(Family = NULL, Group = NULL, nAlienSurveys = NULL, nNativeSurveys = NULL,
                         nAlienBranches = NULL, nNativeBranches = NULL, estimate = NULL, se = NULL, 
                         l95 = NULL, u95 = NULL, p = NULL)

for (f in c('All', familyStats$Family)) {
  for (a in arthropods$Group) {
    
    tmp = comparingNativeAlien(ccPlants, arthGroup = a, plantFamily = f, jdRange = c(152, 194), minArths = 5)
    
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
# Comparing % of surveys with arthropod (+- 95% CI) for all plant families and arthropod groups

par(mar = c(2, 0, 2, 0), oma = c(3, 18, 0, 1), mfrow = c(1,6), mgp = c(3, .5, 0))

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
       main = arthropods$GroupLabel[arthropods$Group == a])
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
          las = 1, line = 1)
    legend("topleft", inset=c(-1,0), legend=c("native","alien"), pch=c(16,1), lty = 'solid', 
           xpd = NA, cex = 1.5)
  }
}
mtext("% of surveys", 1, outer = TRUE, line = 1.5, cex = 1.5)


#######################################################################################################
# Controlling for geographic variation by focusing on a single family, Aceraceae, with the most surveys

ccPlants %>% 
  filter(Family=="Aceraceae", julianday >= 152, julianday <=194) %>% 
  distinct(ID, sciName, plantOrigin, Region) %>% 
  count(sciName, plantOrigin, Region) %>% 
  arrange(Region, desc(n))

# This shows that the following regions have sufficient surveys (>50 per plantOrigin) for a comparison:
regions = list('ON',
               c('MA','RI'),
               c('DC', 'MD', 'VA'),
               'NC')

acerComparisons = data.frame(Region = NULL, Group = NULL, nAlienSurveys = NULL, nNativeSurveys = NULL,
                             nAlienBranches = NULL, nNativeBranches = NULL, estimate = NULL, se = NULL, 
                             l95 = NULL, u95 = NULL, p = NULL)

for (r in 1:length(regions)) {
  
  ccPlantsTmp = ccPlants %>% 
    filter(Region %in% regions[[r]])
  
  acerTmp = comparingNativeAlien(ccPlantsTmp, 
                                 arthGroup = 'caterpillar', 
                                 plantFamily = 'Aceraceae', 
                                 jdRange = c(145, 201), 
                                 minArths = 5)
  
  acerTmpcomp = data.frame(Region = paste(regions[[r]], collapse = "-"),
                           Group = 'caterpillar',
                           nAlienSurveys = acerTmp$nAlienSurveys,
                           nNativeSurveys = acerTmp$nNativeSurveys,
                           nAlienBranches = acerTmp$nAlienBranches,
                           nNativeBranches = acerTmp$nNativeBranches,
                           nAlienSurvsWithArth = acerTmp$nAlienSurvsWithArth,
                           nNativeSurvsWithArth = acerTmp$nNativeSurvsWithArth,
                           estimate = acerTmp$model$coefficients$cond[2,1],
                           se = acerTmp$model$coefficients$cond[2,2],
                           l95 = acerTmp$confint[1],
                           u95 = acerTmp$confint[2],
                           p = acerTmp$model$coefficients$cond[2,4])
  
  if (is.nan(acerTmpcomp$se)) {
    tmpcomp$estimate = NA
  }
  
  acerComparisons = rbind(acerComparisons, acerTmpcomp)
  
}

acerComparisons = acerComparisons %>%
  mutate(propAlienSurvsWithArth = nAlienSurvsWithArth/nAlienSurveys,
         propNativeSurvsWithArth = nNativeSurvsWithArth/nNativeSurveys,
         # Below, multiplying by 1.96 to get half-width of 95% CI
         errorAlienSurvsWithArth = 1.96*((propAlienSurvsWithArth)*(1 - propAlienSurvsWithArth)/
                                           nAlienSurveys)^.5,
         errorNativeSurvsWithArth = 1.96*((propNativeSurvsWithArth)*(1 - propNativeSurvsWithArth)/
                                            nNativeSurveys)^.5,
         
         # z = (p1 - p2)/((p1*(1-p1)/n1) + (p2*(1-p2)/n2))^.5 
         # which I think is preferable when n1 and n2 might differ substantially and/or if p1 or p2 = 0
         
         # As opposed to also commonly used z = (p1 - p2)/((p*(1-p)*(1/n1 + 1/n2)))^.5
         propTestZ = (propNativeSurvsWithArth - propAlienSurvsWithArth)/
           ((propNativeSurvsWithArth*(1-propNativeSurvsWithArth)/nNativeSurveys) + 
              (propAlienSurvsWithArth*(1-propAlienSurvsWithArth)/nAlienSurveys))^.5,
         propTestP = 2*pnorm(q=abs(propTestZ), lower.tail=FALSE),
         pText = case_when(propTestP <= 0.001 & propTestZ >= 0 ~ '+++',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ >= 0 ~ '++',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ >= 0 ~ '+',
                           propTestP <= 0.001 & propTestZ < 0 ~ '---',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ < 0 ~ '--',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ < 0 ~ '-',
                           propTestP > 0.05 ~ ''))

# Plot - vertical - caterpillars on Aceraceae
par(mar = c(6, 6, 1, 1), mfrow = c(1,1), mgp = c(3, 1.2, 0), oma = c(0,0,0,0), xpd = FALSE)

horizOffset = 0.05

plot(4:1 + horizOffset, 100*acerComparisons$propNativeSurvsWithArth, 
     xlab = "", ylab = "% of surveys", xaxt = "n", tck = -0.02, xlim = c(0.5, 4.5), ylim = c(0, 18),
     pch = 16, col = 'black', cex = 2, cex.axis = 1.5, las = 1, cex.lab = 2)
segments(4:1 + horizOffset, 100*acerComparisons$propNativeSurvsWithArth - 
           100*acerComparisons$errorNativeSurvsWithArth, 
         4:1 + horizOffset,
         100*acerComparisons$propNativeSurvsWithArth + 
           100*acerComparisons$errorNativeSurvsWithArth, 
         lwd = 2)

points(4:1 - horizOffset, 100*acerComparisons$propAlienSurvsWithArth, pch = 1, cex = 2)
segments(4:1 - horizOffset, 100*acerComparisons$propAlienSurvsWithArth - 
           100*acerComparisons$errorAlienSurvsWithArth, 
         4:1 - horizOffset,
         100*acerComparisons$propAlienSurvsWithArth + 
           100*acerComparisons$errorAlienSurvsWithArth, 
         lwd = 2)

#abline(h = 1.5, lwd = 2)

text(4:1, 13, labels = acerComparisons$pText, 
     # - symbols are much smaller than +, so making them a larger font size (2 vs 1.5)
     cex = ifelse(acerComparisons$propTestZ < 0, 2, 1.5))

mtext(paste0(acerComparisons$Region, 
             "\n(", acerComparisons$nNativeSurveys, ", ", 
             acerComparisons$nAlienSurveys, ")"),
      1, at = 4:1, las = 1, line = 3, cex = 1.5)

legend("topleft", legend=c("native","alien"), pch=c(16,1), lty = 'solid', 
       xpd = NA, cex = 1.5)

