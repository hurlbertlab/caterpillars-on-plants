# Functions for analyzing the Caterpillars Count! dataset and comparing arthropod abundance or occurrence on different plant species

library(tidyverse)
library(lubridate)

# Looking at meanDensity, Biomass, and fracSurveys of arthropods on different plant species
AnalysisBySciName = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                             ordersToInclude = 'All',    # or 'caterpillar' like arthGroup
                             minLength = 0,              # minimum arthropod size to include 
                             jdRange = c(132, 232),      # change range of days
                             outlierCount = 10000,       # threshold for outliers to be removed
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
# This is done using 2 methods:
# First, a two sample proportion test on proportion of surveys reporting a given arthropod group.
# Second, a zero-inflated negative binomial regression due to the large number of 0's in this count data.
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
  
  nAlienSurvsWithArth = n_distinct(surveyEvents$ID[surveyEvents$Quantity > 0 &
                                                     surveyEvents$plantOrigin == 'alien'])
  nNativeSurvsWithArth = n_distinct(surveyEvents$ID[surveyEvents$Quantity > 0 &
                                                      surveyEvents$plantOrigin == 'native'])
  
  propAlienSurvsWithArth = nAlienSurvsWithArth/alienSurveyEvents
  propNativeSurvsWithArth = nNativeSurvsWithArth/nativeSurveyEvents
  
  errorAlienSurvsWithArth = 1.96*((propAlienSurvsWithArth)*(1 - propAlienSurvsWithArth)/
                                    alienSurveyEvents)^.5
  errorNativeSurvsWithArth = 1.96*((propNativeSurvsWithArth)*(1 - propNativeSurvsWithArth)/
                                     nativeSurveyEvents)^.5
  
  propTestZ = (propNativeSurvsWithArth - propAlienSurvsWithArth)/
    ((propNativeSurvsWithArth*(1-propNativeSurvsWithArth)/nativeSurveyEvents) + 
       (propAlienSurvsWithArth*(1-propAlienSurvsWithArth)/alienSurveyEvents))^.5
  
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
              nNativeSurvsWithArth = nNativeSurvsWithArth,
              nAlienSurvsWithArth = nAlienSurvsWithArth,
              propAlienSurvsWithArth = propAlienSurvsWithArth,
              propNativeSurvsWithArth = propNativeSurvsWithArth,
              errorAlienSurvsWithArth = errorAlienSurvsWithArth,
              errorNativeSurvsWithArth = errorNativeSurvsWithArth,
              propTestZ = propTestZ,
              propTestP = 2*pnorm(q=abs(propTestZ), lower.tail=FALSE),
              model = summary(zinfmodel),
              confint = ci[2, 1:2]))
}
