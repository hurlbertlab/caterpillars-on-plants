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


##########################################################################################################
# Figure 1. Comparison of % surveys with caterpillars across tree species

plantList = officialPlantList %>%
  distinct(sciName, rank)

treeFams = data.frame(Family = c('Fagaceae', 'Betulaceae', 'Aceraceae', 'Caprifoliaceae', 'Juglandaceae'),
                                 #'Ericaceae', 'Hippocastanaceae', 
                                 #'Oleaceae', 'Rosaceae'),
                      famcolor = c(#rgb(0,0,0),
                                rgb(230/255, 159/255, 0),
                                rgb(86/255, 180/255, 233/255),
                                rgb(0, 158/255, 115/255),
                                rgb(240/255, 228/255, 66/255),
                                rgb(213/255, 94/255, 0)))
                                #rgb(0, 114/255, 178/255)))
                                #rgb(204/255, 121/255, 167/255),
                                #'magenta'))

# Tree species with at least 50 surveys, ranked by % of surveys with caterpillars
byTreeSpp = AnalysisBySciName(ccPlants, ordersToInclude = 'caterpillar', jdRange = c(152, 194)) %>%
  left_join(plantList, by = 'sciName') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName')) %>%
  filter(rank == 'species',
         nSurveys >= 50,
         nBranches >= 5) %>%
  mutate(color = ifelse(plantOrigin == 'native', 'gray70', 'firebrick2')) %>%
  arrange(desc(fracSurveys)) %>%
  left_join(treeFams, by = 'Family')

byTreeSpp$famcolor[is.na(byTreeSpp$famcolor)] = 'gray50'

                                

# See alternate plot below
#pdf('Figures/Figure1_ranking_tree_spp.pdf', height = 10, width = 6)
par(mar = c(6, 8, 1, 1), mgp = c(3, 1, 0), mfrow = c(1,1), oma = c(0, 0, 0, 0), xpd = NA)
plot(byTreeSpp$fracSurveys, nrow(byTreeSpp):1, yaxt = 'n', ylab = '', xlab = '% of surveys with caterpillars',
     cex.axis = 1.5, cex.lab = 2, pch = 16, col = byTreeSpp$color, xlim = c(0, 32), ylim = c(5, nrow(byTreeSpp))-2,
     cex = 2*log10(byTreeSpp$nSurveys)/max(log10(byTreeSpp$nSurveys)))
segments(byTreeSpp$LL95frac, nrow(byTreeSpp):1, byTreeSpp$UL95frac, nrow(byTreeSpp):1,
         lwd = 2, col = byTreeSpp$color)
mtext(byTreeSpp$sciName, 2, at = nrow(byTreeSpp):1, line = 1, adj = 1, las = 1, cex = .6)
points(rep(-2, nrow(byTreeSpp)), nrow(byTreeSpp):1, pch = 15, col = byTreeSpp$famcolor)    # color coding by family doesn't add much

legend("bottomright", c("native", "alien", "", "Fagaceae", "Betulaceae", "Aceraceae", "Caprifoliaceae", "Juglandaceae", "Other"), 
       col = c('gray70', 'firebrick2', NA, rgb(230/255, 159/255, 0),
               rgb(86/255, 180/255, 233/255),
               rgb(0, 158/255, 115/255),
               rgb(240/255, 228/255, 66/255),
               rgb(213/255, 94/255, 0),
               'gray50'), 
       pch = c(16, 16, NA, rep(15, 6)), lty = c(rep('solid', 2), rep('blank', 7)), 
       cex = c(rep(1.75, 2), rep(1.3, 7)), pt.cex = c(rep(2, 2), rep(1.8, 7)), lwd = c(rep(2, 2), rep(0, 7)))

caterpillar = readPNG('images/caterpillar.png')
rasterImage(caterpillar, 16, 35, 32, 49)
#dev.off()



# Alternative two column figure

numspp = nrow(byTreeSpp)
if (numspp %% 2 != 0) { numspp = numspp + 1 } # add 1 to make numspp even if necessary

pdf('Figures/Figure1_ranking_tree_spp_2col.pdf', height = 7, width = 12)
par(mar = c(6, 7, 3, 1), mgp = c(3, 1, 0), mfrow = c(1,2), oma = c(0, 0, 0, 0), xpd = NA)
plot(byTreeSpp$fracSurveys[1:(numspp/2)], (numspp/2):1, yaxt = 'n', ylab = '', xlab = '% of surveys',
     cex.axis = 1.5, cex.lab = 2, pch = 16, col = byTreeSpp$color[1:(numspp/2)], 
     xlim = c(0, 32), ylim = c(1, numspp/2),
     cex = 2*log10(byTreeSpp$nSurveys[1:(numspp/2)])/max(log10(byTreeSpp$nSurveys[1:(numspp/2)])),
     main = "Species rank 1-44")
segments(byTreeSpp$LL95frac[1:(numspp/2)], (numspp/2):1, byTreeSpp$UL95frac[1:(numspp/2)], (numspp/2):1,
         lwd = 2, col = byTreeSpp$color[1:(numspp/2)])
mtext(byTreeSpp$sciName[1:(numspp/2)], 2, at = (numspp/2):1, line = 1, adj = 1, las = 1, cex = .6)
points(rep(-2, numspp/2), (numspp/2):1, pch = 15, col = byTreeSpp$famcolor[1:(numspp/2)], cex = 1.2)

legend("bottomright", c("Fagaceae", "Betulaceae", "Aceraceae", "Caprifoliaceae", "Juglandaceae", "Other"), 
       col = c(rgb(230/255, 159/255, 0),
               rgb(86/255, 180/255, 233/255),
               rgb(0, 158/255, 115/255),
               rgb(240/255, 228/255, 66/255),
               rgb(213/255, 94/255, 0),
               'gray50'), 
       pch = 15, cex = 1, pt.cex = 1.5)

plot(byTreeSpp$fracSurveys[(numspp/2 + 1):numspp], (numspp/2):1, yaxt = 'n', ylab = '', xlab = '% of surveys',
     cex.axis = 1.5, cex.lab = 2, pch = 16, col = byTreeSpp$color[(numspp/2 + 1):numspp], 
     xlim = c(0, 32), ylim = c(1, numspp/2),
     cex = 2*log10(byTreeSpp$nSurveys[(numspp/2 + 1):numspp])/max(log10(byTreeSpp$nSurveys[(numspp/2 + 1):numspp])),
     main = "Species rank 45-88")
segments(byTreeSpp$LL95frac[(numspp/2 + 1):numspp], (numspp/2):1, byTreeSpp$UL95frac[(numspp/2 + 1):numspp], (numspp/2):1,
         lwd = 2, col = byTreeSpp$color[(numspp/2 + 1):numspp])
mtext(byTreeSpp$sciName[(numspp/2 + 1):numspp], 2, at = (numspp/2):1, line = 1, adj = 1, las = 1, cex = .6)
points(rep(-2, numspp/2), (numspp/2):1, pch = 15, col = byTreeSpp$famcolor[(numspp/2 + 1):numspp], cex = 1.2)    # 

legend("bottomright", c("native", "alien"), 
       col = c('gray70', 'firebrick2'), 
       pch = 16, cex = 1.25, pt.cex = 1.8, lwd = 2, lty = 'solid')

caterpillar = readPNG('images/caterpillar.png')
rasterImage(caterpillar, 16, 35, 32, 45)
dev.off()





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
  group_by(sciName) %>%
  summarize(nTaxa = n_distinct(TaxonName[Group == 'caterpillar' & Rank != 'order']),
            nSurveysWithPhotos = n_distinct(ID),
            nBranches = n_distinct(Code),
            Regions = paste(unique(Region), collapse = '-')) %>%
  arrange(desc(nTaxa)) %>%
  mutate(label = paste0(substr(sciName, 1, 1), ". ", word(sciName, 2))) %>%
  left_join(plantOrigin[, c('scientificName', 'plantOrigin', 'Family')], by = c('sciName' = 'scientificName')) %>%
  mutate(color = ifelse(plantOrigin == 'native', 'gray50', 'firebrick2'))



# Linear models
lm.native = lm(log10(catSpecies$nTaxa[catSpecies$plantOrigin == 'native' & catSpecies$nTaxa > 0]) ~ 
                 log10(catSpecies$nSurveysWithPhotos[catSpecies$plantOrigin == 'native' & catSpecies$nTaxa > 0]))

lm.alien = lm(log10(catSpecies$nTaxa[catSpecies$plantOrigin == 'alien' & catSpecies$nTaxa > 0]) ~ 
                log10(catSpecies$nSurveysWithPhotos[catSpecies$plantOrigin == 'alien' & catSpecies$nTaxa > 0]))

# Plot
pdf('Figures/Figure2_caterpillar_taxa.pdf', height = 6, width = 8)
par(mgp = c(4, 1, 0), tck = -0.03, mar = c(5, 7, 1, 1))
plot(log10(catSpecies$nSurveysWithPhotos), log10(catSpecies$nTaxa), pch = ifelse(catSpecies$plantOrigin == 'native', 16, 17), 
     xlab = expression(log[10] ~ "#" ~ surveys ~ with ~ photos), 
     ylab = expression(log[10] ~ "#" ~ caterpillar ~ taxa), las = 1, 
     cex.axis = 1.5, cex.lab = 2, col = catSpecies$color, cex = ifelse(catSpecies$plantOrigin == 'native', 1.3, 1.5))

abline(lm.native, lwd = 2, col = 'gray30')
abline(lm.alien, col = 'firebrick2', lwd = 2)

legend("topleft", c("native", "alien"), pch = c(16, 17), col = c('gray50', 'firebrick2'), cex = 1.5)
dev.off()




###############################################################################
# Find set of plant families with sufficient data for native-alien comparisons

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
# Figure 3. Comparing % of surveys with arthropod (+- 95% CI) for all plant families and arthropod groups


pdf('Figures/Figure3_plant_family_comparison.pdf', height = 6, width = 10)
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

dev.off()




############################################################################################################
# Figure 4. Controlling for geographic variation by focusing on a single family, Aceraceae, with the most surveys

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
                             nAlienBranches = NULL, nNativeBranches = NULL, nAlienSurvsWithArth = NULL,
                             nNativeSurvsWithArth = NULL, propAlienSurvsWithArth = NULL,
                             propNativeSurvsWithArth = NULL, errorAlienSurvsWithArth = NULL,
                             errorNativeSurvsWithArth = NULL, propTestZ = NULL,
                             propTestP = NULL,
                             estimate = NULL, se = NULL, 
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
                           propAlienSurvsWithArth = acerTmp$propAlienSurvsWithArth,
                           propNativeSurvsWithArth = acerTmp$propNativeSurvsWithArth,
                           errorAlienSurvsWithArth = acerTmp$errorAlienSurvsWithArth,
                           errorNativeSurvsWithArth = acerTmp$errorNativeSurvsWithArth,
                           propTestZ = acerTmp$propTestZ,
                           propTestP = acerTmp$propTestP,
                           estimate = acerTmp$model$coefficients$cond[2,1],
                           se = acerTmp$model$coefficients$cond[2,2],
                           l95 = acerTmp$confint[1],
                           u95 = acerTmp$confint[2],
                           p = acerTmp$model$coefficients$cond[2,4])
  
  if (is.nan(acerTmpcomp$se)) {
    acerTmpcomp$estimate = NA
  }
  
  acerComparisons = rbind(acerComparisons, acerTmpcomp)
  
}

acerComparisons = acerComparisons %>%
  mutate(pText = case_when(propTestP <= 0.001 & propTestZ >= 0 ~ '+++',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ >= 0 ~ '++',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ >= 0 ~ '+',
                           propTestP <= 0.001 & propTestZ < 0 ~ '---',
                           propTestP > 0.001 & propTestP <= 0.01 & propTestZ < 0 ~ '--',
                           propTestP > 0.01 & propTestP <= 0.05 & propTestZ < 0 ~ '-',
                           propTestP > 0.05 ~ ''))

# Plot - vertical - caterpillars on Aceraceae
pdf('Figures/Figure3_Aceraceae_across_regions.pdf', height = 6, width = 8)
par(mar = c(7, 6, 1, 1), mfrow = c(1,1), mgp = c(3, 1.2, 0), oma = c(0,0,0,0), xpd = FALSE)

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
      1, at = 4:1, las = 1, line = 2, cex = 1.5)

mtext("Region", 1, line = 5, cex = 2)

legend("topleft", legend=c("native","alien"), pch=c(16,1), lty = 'solid', 
       xpd = NA, cex = 1.5)

maple = readPNG("images/maple.png")
rasterImage(maple, 3.9, 14, 4.6, 19)
rasterImage(caterpillar, 3.3, 15, 4.1, 18)

dev.off()


