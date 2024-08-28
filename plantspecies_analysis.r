# Script for calculating average abundance, biomass, and occurrence per plant species, and for comparing broadly between native and alien plants

library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(vioplot)
library(RCurl)
library(rvest)
library(xml2)
library(glmmTMB)


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
  

# FUNCTIONS
# Looking at meanDensity, Biomass, and fracSurveys of caterpillars in different plant families 
AnalysisBySciName = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                             ordersToInclude = 'All',    # or 'caterpillar' like arthGroup
                             minLength = 0,              # minimum arthropod size to include 
                             jdRange = c(132, 232),      # change range of days
                             outlierCount = 10000,       # Outliers 
                             excludeWetLeaves = TRUE)    # exclude surveys with wet leaves                  
  
{ if(length(ordersToInclude)==1 & ordersToInclude[1]=='All') {
  ordersToInclude = unique(surveyData$Group)
}
  
  numUniqueBranches = length(unique(surveyData$PlantFK))
  
  firstFilter = surveyData %>%
    filter(julianday >= jdRange[1], julianday <= jdRange[2]) %>% #subscription notation
    mutate(julianweek = 7*floor(julianday/7) + 4)
  
  if (excludeWetLeaves) {
    firstFilter = filter(firstFilter, !WetLeaves)
  }
  
  effortBySciName = firstFilter %>%
    group_by(sciName) %>% 
    summarize(nSurveys = n_distinct(ID),
              nBranches = n_distinct(PlantFK),
              nSites = n_distinct(Name),
              nRegions = n_distinct(Region))
  
  arthCount = firstFilter %>%
    filter(Length >= minLength, 
           Group %in% ordersToInclude) %>%
    mutate(Quantity2 = ifelse(Quantity > outlierCount, 1, Quantity)) %>% #outlier counts replaced with 1
    group_by(sciName) %>%
    summarize(totalCount = sum(Quantity2, na.rm = TRUE),
              numSurveysGTzero = length(unique(ID[Quantity > 0])),
              totalBiomass = sum(Biomass_mg, na.rm = TRUE)) %>% 
    right_join(effortBySciName, by = 'sciName') %>%
    mutate(meanDensity = totalCount/nSurveys,
           meanBiomass = totalBiomass/nSurveys,
           fracSurveys = 100*numSurveysGTzero/nSurveys,
           LL95frac = fracSurveys - 1.96*(fracSurveys*(100-fracSurveys)/nSurveys)^.5,
           UL95frac = fracSurveys + 1.96*(fracSurveys*(100-fracSurveys)/nSurveys)^.5) %>%
    arrange(sciName) %>%
    data.frame()
  
  arthCount[is.na(arthCount)] = 0
  
  return(arthCount)
}


# Function for comparing arthropod abundance at the branch-scale between native and alien tree species.
# This is done using zero-inflated negative binomial regression due to the large number of 0's in this count data.
# This can be done across all families (in which case plant species is a random effect nested within Family), or
# within a specified plant family (in which case plant species is a random effect).
comparingNativeAlien = function(surveyData, 
                                arthGroup = 'caterpillar',            # 'caterpillar', 'spider', etc.
                                plantFamily = 'All',                  # Plant Family for comparison
                                obsMethod = c('Visual', 'Beat sheet'),# or 'Visual', or 'Beat sheet'
                                jdRange = c(132, 232),                # Range of days over which surveys were done   
                                minSurveys = 10,                      # min # of survey events per group (native/alien)
                                minBranches = 5,                      # min # of unique branches per group (native/alien) 
                                minArths = 10,                        # min # of surveys with at least 1 arthropod
                                excludeWetLeaves = TRUE)              # exclude surveys with wet leaves
{

  require(glmmTMB)
  
  # Warning messages
  if (!arthGroup %in% unique(surveyData$Group)) {
    stop("Not a valid arthropod group name. Should be one of: 'ant', 'aphid', 'bee', 'beetle', 'caterpillar', 
          'daddlylonglegs', 'fly', 'grasshopper', 'leafhopper', 'moths', 'truebugs'.")
  }
  
  if (!plantFamily %in% c('All', unique(surveyData$Family))) {
    stop("Not a valid plant Family name.")
  }
  
  
  if (plantFamily == 'All') {
    fam = unique(surveyData$Family) 
  } else {
    fam = plantFamily
  }
  
  survData = filter(surveyData, 
                    Family %in% fam,
                    plantOrigin %in% c('alien', 'native'),
                    ObservationMethod %in% obsMethod,
                    julianday >= jdRange[1], 
                    julianday <= jdRange[2])
  
  if (excludeWetLeaves) {
    survData = filter(survData, !WetLeaves)
  }
  
  arthData = filter(survData, Group == arthGroup) %>%
    select(ID, Quantity, Length)

  surveyEvents = survData %>%
    distinct(ID, PlantFK, sciName, Family, plantOrigin) %>%
    left_join(arthData, by = 'ID')
  
  surveyEvents[is.na(surveyEvents)] = 0 # fill in 0's where Quantity or Length is NA
  
  nArthRecords = n_distinct(surveyEvents$ID[surveyEvents$Quantity > 0])
  
  nativeSurveyEvents = n_distinct(surveyEvents$ID[surveyEvents$plantOrigin == 'native'])
  alienSurveyEvents = n_distinct(surveyEvents$ID[surveyEvents$plantOrigin == 'alien'])
  
  nativeSurveyBranches = n_distinct(surveyEvents$PlantFK[surveyEvents$plantOrigin == 'native'])
  alienSurveyBranches = n_distinct(surveyEvents$PlantFK[surveyEvents$plantOrigin == 'alien'])
  
  if (nArthRecords < minArths) {
    stop("Not enough surveys with arthropods")
  }
  
  if (nativeSurveyEvents < minSurveys | alienSurveyEvents < minSurveys) {
    stop("Not enough survey events for the comparison.")
  }
  
  if (nativeSurveyBranches < minBranches | alienSurveyBranches < minBranches) {
    stop("Not enough unique survey branches for the comparison.")
  }
  
  
  # Fit zero-inflated negative binomial
  if (plantFamily == 'All') {

    zinfmodel = glmmTMB(Quantity ~ plantOrigin + (1 | Family / sciName), 
                        ziformula = ~1, 
                        family=nbinom2, 
                        data = surveyEvents)
  } else {
    
    zinfmodel = glmmTMB(Quantity ~ plantOrigin + (1 | sciName), 
                        ziformula = ~1, 
                        family=nbinom2, 
                        data = surveyEvents)
  }
  
  ci = confint(zinfmodel)
  
  return(list(Group = arthGroup,
              Family = plantFamily,
              nNativeSurveys = nativeSurveyEvents,
              nAlienSurveys = alienSurveyEvents,
              nNativeBranches = nativeSurveyBranches,
              nAlienBranches = alienSurveyBranches,
              nNativeSurvsWithArth = n_distinct(surveyEvents$ID[surveyEvents$Quantity > 0 &
                                                                  surveyEvents$plantOrigin == 'native']),
              nAlienSurvsWithArth = n_distinct(surveyEvents$ID[surveyEvents$Quantity > 0 &
                                                                 surveyEvents$plantOrigin == 'alien']),
              model = summary(zinfmodel),
              confint = ci[2, 1:2]))
}



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
         nSurveysA >= 10,
         nSurveysN >= 10)
  

arthropods = data.frame(Group = c('caterpillar', 'spider', 'leafhopper', 'beetle', 'truebugs', 'ant'),
                        GroupLabel = c('caterpillars', 'spiders', 'hoppers', 'beetles', 'true bugs', 'ants'),
                        color = c('limegreen', 'gray50', 'dodgerblue', 'salmon', 'magenta', 'orange'),
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

comparisons = comparisons %>%
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
                           propTestP > 0.05 ~ '')
  )

# Plotting comparisons
# --one issue with this plot is that if there are ZERO observations for a group (e.g. caterpillars on alien Cornaceae),
#   then model fails to converge and no estimate or SE are available
par(mar = c(2, 0, 2, 0), oma = c(3, 18, 0, 1), mfrow = c(1,6), mgp = c(3, .5, 0))
for (a in arthropods$Group) {
  plot(comparisons$estimate[comparisons$Group == a], 1:length(unique(comparisons$Family)),
       xlab = "", ylab = "", yaxt = "n", tck = -0.03,
       pch = 16, col = arthropods$color[arthropods$Group == a], cex = 2, xlim = c(-2, 5), 
       main = arthropods$GroupLabel[arthropods$Group == a])
  segments(comparisons$l95[comparisons$Group == a], 1:length(unique(comparisons$Family)),
        comparisons$u95[comparisons$Group == a], 1:length(unique(comparisons$Family)), 
        col = arthropods$color[arthropods$Group == a])
  abline(v = 0, lty = 'dashed')
  
  if (a == arthropods$Group[1]) {
    mtext(paste0(comparisons$Family[comparisons$Group == a], 
                " (", comparisons$nAlienSurveys[comparisons$Group == a], ", ", 
                comparisons$nNativeSurveys[comparisons$Group == a], ")"),
          2, at = 1:length(unique(comparisons$Family)), 
          las = 1, line = 1)
  }
}
mtext("log Native / Alien abundance", 1, outer = TRUE, line = 1.5, cex = 1.5)


# Second comparison plot shows separate alien and native estimates +- 95% CI for % of surveys with arthropod
par(mar = c(2, 0, 2, 0), oma = c(3, 18, 0, 1), mfrow = c(1,6), mgp = c(3, .5, 0))

vertOffset = 0.1

for (a in arthropods$Group) {
  plot(100*comparisons$propNativeSurvsWithArth[comparisons$Group == a], 
       1:length(unique(comparisons$Family)) + vertOffset,
       xlab = "", ylab = "", yaxt = "n", tck = -0.03, xlim = c(0, 44), ylim = c(1, 15),
       pch = 16, col = arthropods$color[arthropods$Group == a], cex = 1.8, 
       main = arthropods$GroupLabel[arthropods$Group == a])
  segments(100*comparisons$propNativeSurvsWithArth[comparisons$Group == a] - 
             100*comparisons$errorNativeSurvsWithArth[comparisons$Group == a], 
           1:length(unique(comparisons$Family)) + vertOffset,
           100*comparisons$propNativeSurvsWithArth[comparisons$Group == a] + 
             100*comparisons$errorNativeSurvsWithArth[comparisons$Group == a], 
           1:length(unique(comparisons$Family)) + vertOffset, 
           col = arthropods$color[arthropods$Group == a])
  
  points(100*comparisons$propAlienSurvsWithArth[comparisons$Group == a], 
         1:length(unique(comparisons$Family)) - vertOffset,
       pch = 1, col = arthropods$color[arthropods$Group == a], cex = 1.8)
  segments(100*comparisons$propAlienSurvsWithArth[comparisons$Group == a] - 
             100*comparisons$errorAlienSurvsWithArth[comparisons$Group == a], 
           1:length(unique(comparisons$Family)) - vertOffset,
           100*comparisons$propAlienSurvsWithArth[comparisons$Group == a] + 
             100*comparisons$errorAlienSurvsWithArth[comparisons$Group == a], 
           1:length(unique(comparisons$Family)) - vertOffset, 
           col = arthropods$color[arthropods$Group == a])
  
  #abline(h = 1:length(unique(comparisons$Family))+0.5)
  
  text(# Commented out line below makes horizontal placement relative to largest values
       #100*max(c(comparisons$propAlienSurvsWithArth[comparisons$Group == a], 
       #       comparisons$propNativeSurvsWithArth[comparisons$Group == a])) + 2,
       40,
       1:length(unique(comparisons$Family)),
       labels = comparisons$pText[comparisons$Group == a], 
       # - symbols are much smaller than +, so making them a larger font size (2 vs 1.5)
       cex = ifelse(comparisons$propTestZ[comparisons$Group == a] < 0, 2, 1.5))
  
  # Put plant Family labels along the y-axis for the first plot
  if (a == arthropods$Group[1]) {
    mtext(paste0(comparisons$Family[comparisons$Group == a], 
                 " (", comparisons$nAlienSurveys[comparisons$Group == a], ", ", 
                 comparisons$nNativeSurveys[comparisons$Group == a], ")"),
          2, at = 1:length(unique(comparisons$Family)), 
          las = 1, line = 1)
    legend("topleft", inset=c(-1,0), legend=c("native","alien"), pch=c(16,1), lty = 'solid', 
           xpd = NA, cex = 1.5)
  }
}
mtext("% of surveys", 1, outer = TRUE, line = 1.5, cex = 1.5)




# START HERE
cc_plus_tallamy = read.csv(file = "data/Plant Analysis/cc_plus_tallamy.csv")

#for every plant family return a graph created by comparingBugs... for an arthropod
comparingBugsonNativeVersusAlienPlants <- function(cc_plus_tallamy,  # Original dataset with native/alien info
                            arthGroup,                               # Arthropod to be analyzed
                            plantFamily,                             # Plant family with both native/alien species
                            jdRange = c(132, 232),                   # Range of days
                            minSurveysPerPlant = 10,                 # Minimum number of surveys done per branch
                            plot = FALSE,                            # Plotting the data
                            comparisonVar = "meanDensity")           # 'meanDensity' or 'fracSurveys' or 'meanBiomass' 
  { familiesWithNativeAndAlienSpecies = cc_plus_tallamy %>%
    group_by(Family) %>%
    summarize(NativeSpp = length(unique(sciName[origin == 'native'])),
             AlienSpp = length(unique(sciName[origin == 'alien']))) %>%
    filter(NativeSpp >= 2 & AlienSpp >= 2)
  
  # Further filtering the data, to pull out the rows with a Family fitting these conditions
  plantCount = cc_plus_tallamy %>%
    filter(Family == plantFamily, 
           julianday >= jdRange[1], 
           julianday <= jdRange[2]) %>%
    count(Family, sciName, origin) 
    
  # Counting how many species for each family: 1 for either just alien or just native, 2 for both 
  if (!plantFamily %in% familiesWithNativeAndAlienSpecies$Family) {
    stop("There are either not enough native or alien species to analyze.")
} else {
    # Could make an error message for each problem; nested ifelse statement
    filteredData = cc_plus_tallamy %>% 
      filter(julianday >= jdRange[1], 
             julianday <= jdRange[2],
             Family == plantFamily,
             sciName %in% plantCount$sciName[plantCount$n >= minSurveysPerPlant]) #%>%
      #filter(ObservationMethod == "BeatSheet")
    
    onlyBugs<-AnalysisBySciName(filteredData, ordersToInclude = arthGroup) %>%
      left_join(plantCount, by = "sciName") 
    
    # Separating data sets
    nativeData = filter(onlyBugs, origin == 'native')
    alienData = filter(onlyBugs, origin == 'alien')
    
    # Completing an analysis and pulling out the means, p-value, etc.
    if(comparisonVar != "fracSurveys") {
      x = log10(nativeData[,comparisonVar] + 0.001)
      y = log10(alienData[,comparisonVar] + 0.001)
    } else {
      x = nativeData[,comparisonVar]
      y = alienData[,comparisonVar]
    }
    
    w = wilcox.test(x, y, exact = FALSE)
    p_value = w$p.value
    
    # Change the column name to meanBiomass, meanDensity, or fracSurveys
    median(nativeData$meanDensity)
    median(alienData$meanDensity)
    #nativeMean = w$estimate[1]
    #alienMean = w$estimate[2]
    
    native_pop_size = nrow(nativeData)
    alien_pop_size = nrow(alienData)
    
    # Creating a pdf for each arthGroup
    # Plotting the analysis
    if(plot == TRUE) {
      y_label = comparisonVar
      
      # Colors associated with each family name changed manually: 
      # Native plants are darker colors (blue4, darkred) and alien plants (steelblue1, firebrick1)
      vioplot(x, y, boxwex = 0.5, xaxt = 'n', las = 1, col = c("darkred", "firebrick1"))
      mtext(paste("p =", round(p_value,3)), adj = 0.9, line = -1.5, cex = 1.15, ylim = 20)
      mtext(c(paste("Native"), paste("Non-native")), 1,
            line = 0.8, at = 1:2, cex = 1.15)
    }
  }
}


## A pdf with graphs depicting density, biomass, % surveyed ##
pdf(file = "Figures/RosaceaeWilCoxTest.pdf",
    width = 11, height = 8.5)
par(mfrow = c(4, 3), mar = c(2, 4, 2, 1), oma = c(0,0,1,0))

# Changing the name of the file and the plant family used, manually
for (group in c("caterpillar", "beetle", "truebugs", "spider")) {
  
  for (plotVar in c("meanDensity", "meanBiomass", "fracSurveys")) {
    
    comparingBugsonNativeVersusAlienPlants(cc_plus_tallamy, plantFamily = "Rosaceae", 
                                           arthGroup = group, comparisonVar = plotVar, plot = TRUE)
  }
}
dev.off()


## A graph of ALL the families in the dataset, still separating the arthropod groups
pdf(file = "Figures/AllFamiliesAllArth.pdf",
    width = 11, height = 8.5)
par(mfrow = c(4, 3), mar = c(3, 4, 3, 1))

for (group in c("caterpillar", "beetle", "truebugs", "spider")) {
  
  for (plotVar in c("meanDensity", "meanBiomass", "fracSurveys")) {
    
    allFamilies = AnalysisBySciName(cc_plus_tallamy, ordersToInclude = group, plotVar = analysis) %>%
      left_join(cc_plus_tallamy, by = 'sciName') %>%
      select(sciName, totalCount, numSurveysGTzero, totalBiomass, nSurveys, meanDensity,
           fracSurveys, meanBiomass, origin) %>%
      distinct(sciName, .keep_all = TRUE) 
  
    allNativeFamilies = filter(allFamilies, origin == 'native')
    allAlienFamilies = filter(allFamilies, origin == 'alien')
    
    Allnative_pop_size = nrow(allNativeFamilies)
    Allalien_pop_size = nrow(allAlienFamilies)
  
    # Change the column name to meanBiomass, meanDensity, or fraSurveys
    median(allNativeFamilies$meanDensity)
    median(allAlienFamilies$meanDensity)
    
    if(plotVar != "fracSurveys") {
      x = log10(allNativeFamilies[,plotVar] + 0.001)
      y = log10(allAlienFamilies[,plotVar] + 0.001)
    } else {
      x = allNativeFamilies[,plotVar]
      y = allAlienFamilies[,plotVar]
    }
    
    vioplot(x, y, boxwex = 0.5, xaxt = 'n', col = c("grey28", "grey76"), las = 1)
    w = wilcox.test(x, y, exact = FALSE)
    p_value = w$p.value
    mtext(paste("p =", round(p_value,3)), adj = 0.91, line = -1.5, cex = 1.15, ylim = 25) 
    mtext(c(paste("Native"), paste("Non-native")), 1,
          line = 0.8, at = 1:2, cex = 0.95)
  }
}  
dev.off()


## Code to calculate the lepS comparisons 
plantCountJuneJuly = cleanDatasetCC %>%
  dplyr::filter(julianday >= 132, julianday <= 232) %>% #change range of days 
  distinct(ID, sciName) %>%
  count(sciName) %>%
  arrange(desc(n)) #%>%
  #filter(sciName %in% plantCountJuneJuly$sciName[plantCountJuneJuly$n >= 10])
# Specifies that only plant species that were surveyed at least 10x in June and July were included
SurveyedCertainAmount = cleanDatasetCC %>%
  filter(sciName %in% plantCountJuneJuly$sciName[plantCountJuneJuly$n >= 10])

# Specifies that only caterpillars (not all arthropods) were analyzed in this analysis
onlyCaterpillars = AnalysisBySciName(SurveyedCertainAmount, ordersToInclude = "caterpillar") %>%
  mutate(Genus = word(sciName, 1)) 

# left_join SurveyWithCaterpillar before
lepSandAllFam <- left_join(onlyCaterpillars, tallamy, by = 'Genus') %>%
  select(sciName:Genus,Family, origin..for.analysis., total.Lep.spp, nSurveys, 
         meanDensity, fracSurveys, meanBiomass) %>%
  rename(origin = origin..for.analysis., lepS = total.Lep.spp) %>%
  filter(nSurveys >= 10) %>%
  mutate(color = case_when(Family == "Rosaceae" & origin == "native" ~ "darkred",
                           Family == "Rosaceae" & origin == "alien" ~ "firebrick1",
                           Family == "Oleaceae" & origin == "native" ~ "blue4",
                           Family == "Oleaceae" & origin == "alien" ~ "steelblue1",
                           TRUE ~ "black")) %>%
  mutate(origin2 = case_when(origin == "native" ~ 1,
                             origin == "alien" ~ 2))


## Figure 6: Comparing the average *caterpillar* density, biomass, and fracSurveys 
## per survey to lepS and conducting a linear regression 
pdf(file = "Figures/LepidopteraAnalysis.pdf", 
    width = 9, height = 6)
par(mfrow = c(2, 2), mar = c(5,5,2,1))


plot(lepSandAllFam$lepS, log10(lepSandAllFam$meanDensity), xlab = "Genus-level Lepidoptera Richness", 
     ylab = "Density per Branch", col = lepSandAllFam$color, las = 1,
     pch = ifelse(lepSandAllFam$origin2 == 1, 16, 17), cex = log10(lepSandAllFam$nSurveys)/2)
lm.density = lm(log10(meanDensity[meanDensity > 0]) ~ lepS[meanDensity > 0], data = lepSandAllFam)
p_value = summary(lm.density)$coefficients[2,4]
text(48, 0.45, bquote(R^2==.(round(summary(lm.density)$r.squared, 2))), cex = 1.05)
mtext(paste("p =", round(p_value,3)), line = -1.15, adj = 0.05)
abline(lm.density)

plot(log10(lepSandAllFam$lepS), log10(lepSandAllFam$meanBiomass), xlab = "Genus-level Lepidoptera Richness", 
     ylab = "Biomass per Branch", col = lepSandAllFam$color, las = 1,
     pch = ifelse(lepSandAllFam$origin2 == 1, 16, 17), cex = log10(lepSandAllFam$nSurveys)/2)
lm.biomass = lm(log10(meanBiomass[meanBiomass > 0]) ~ lepS[meanBiomass > 0], data = lepSandAllFam)
p_value = summary(lm.biomass)$coefficients[2,4]
text(0.25, 2, bquote(R^2==.(round(summary(lm.biomass)$r.squared, 2))), cex = 1.05)
mtext(paste("p =", round(p_value,3)), line = -1.15, adj = 0.05)
abline(lm.biomass)


plot(lepSandAllFam$lepS, lepSandAllFam$fracSurveys, xlab = "Genus-level Lepidoptera Richness", 
     ylab = "% Occurrence per Branch", col = lepSandAllFam$color, las = 1,
     pch = ifelse(lepSandAllFam$origin2 == 1, 16, 17), cex = log10(lepSandAllFam$nSurveys)/2)
lm.surveys = lm(fracSurveys ~ lepS, data = lepSandAllFam)
p_value = summary(lm.surveys)$coefficients[2,4]
text(48, 50, bquote(R^2==.(round(summary(lm.surveys)$r.squared, 2))), cex = 1.05)
mtext(paste("p =", round(p_value,3)), line = -1.15, adj = 0.05)
abline(lm.surveys)

plot(x = 1, type = 'n', ylab = "", yaxt = 'n', xlab = "", xaxt = 'n')

legend("center", legend = c("Rosaceae Natives", "Rosaceae Non-natives", "Oleaceae Natives", "Oleaceae Non-natives", "Native Other Families", "Non-native Other Families"), 
       col = c("darkred","firebrick1", "blue4", "steelblue1", "black", "black"),
       pch= c(16, 17, 16, 17, 16, 17), pt.cex = 1.5, cex = 1.5)

dev.off()


## Summary of the different species of native vs. alien plants in the two families used for analysis ##
plantCountJuneJuly = cleanDatasetCC %>%
  dplyr::filter(julianday >= 132, julianday <= 232) %>% #change range of days 
  distinct(ID, sciName) %>%
  count(sciName) %>%
  arrange(desc(n))  

# Specifies that only plant species that were surveyed at least 10x in June and July were included
SurveyedCertainAmount = cleanDatasetCC %>%
  filter(sciName %in% plantCountJuneJuly$sciName[plantCountJuneJuly$n >= 6])

# Specifies that only caterpillars (not all arthropods) were analyzed in this analysis
AllArths = AnalysisBySciName(SurveyedCertainAmount, ordersToInclude = "All") %>%
  mutate(Genus = word(sciName, 1)) 

ArthsAndTallamy = left_join(AllArths, tallamy, by = 'Genus') %>%
  select(sciName:Genus,Family, origin..for.analysis., total.Lep.spp, nSurveys, meanDensity, fracSurveys, meanBiomass) %>%
  rename(origin = origin..for.analysis., lepS = total.Lep.spp) %>%
  mutate(color = case_when(Family == "Rosaceae" & origin == "native" ~ "darkred",
                           Family == "Rosaceae" & origin == "alien" ~ "firebrick1",
                           Family == "Oleaceae" & origin == "native" ~ "blue4",
                           Family == "Oleaceae" & origin == "alien" ~ "steelblue1",
                           TRUE ~ "black"))

nativeSpecies = filter(ArthsAndTallamy, origin == 'native')  
alienSpecies = filter(ArthsAndTallamy, origin == 'alien')         

# only the top 10 or so
nativeSpecies_desc = nativeSpecies[order(-nativeSpecies$nSurveys),] 
nativeTop20 = nativeSpecies_desc[1:20,]

nativeTop = nativeSpecies_desc %>%
  filter(Family %in% c("Rosaceae", "Oleaceae") | nSurveys >= 500)

alienSpecies_desc = alienSpecies[order(-alienSpecies$nSurveys),]



# Figure 2 code:
pdf(file = "Figures/BranchesVsSpecies.pdf", width = 8, height = 10)
par(mfrow = c(2, 1), mar = c(4,5,2,1), mgp = c(3.5,.75,0))

barplot2 = barplot(log10(nativeTop$nSurveys), ylab = "log10 # of Surveys",
                   pch = 16, main ="(A) Native Species",
                   col = nativeTop$color,
                   xaxt = "n", yaxt="n", ylim = c(0, 4))
axis(2, at = 0:4, labels = c(1, 10, 100, 1000, 10000), las = 1)
legend("topright", legend = c("Rosaceae Natives", "Oleaceae Natives", "Other Families"), 
       col = c("darkred", "blue4", "black"), pch = 15, horiz = FALSE, pt.cex = 1,
       cex = 0.9, bty = "n")
text(barplot2 + 0.3, -.13, labels = nativeTop$sciName, cex = .5, 
     xpd = NA, srt=35, adj=1)

barplot1 = barplot(log10(alienSpecies_desc$nSurveys), 
        ylab = "log 10 # of Surveys",
        pch = 16, main ="(B) Non-native Species", 
        col = alienSpecies_desc$color, las = 1, yaxt = 'n', ylim = c(0,3))
legend("topright", legend = c("Rosaceae Non-natives", "Oleaceae Non-natives", "Other Families"), 
       col = c("firebrick1", "steelblue1", "black"), pch = 15, horiz = FALSE, pt.cex = 1,
       cex = 0.9, bty = "n")
text(barplot1 +.3, -0.05, labels = alienSpecies_desc$sciName, cex = 0.7, 
     xpd = NA, srt=35, adj=1)
axis(2, at = 0:2, labels = c(1, 10, 100), las = 1)


dev.off()

