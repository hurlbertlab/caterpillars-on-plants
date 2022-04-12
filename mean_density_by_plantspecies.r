library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(maps)
library(sp)
library(maptools)

#Added bc "mutate_cond" is not a built-in function
mutate_cond <- function(.data, condition,...,envir=parent.frame()){
  condition <- eval(substitute(condition),.data,envir)
  .data[condition,] <- .data[condition,] %>% mutate(...)
  .data
}

meanDensityBySpecies = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                             ordersToInclude = 'All',       # or 'caterpillar'
                             
                             minLength = 0,         # minimum arthropod size to include 
                             jdRange = c(1,365),
                             outlierCount = 10000,
                             plotVar = 'meanDensity', # 'meanDensity' or 'fracSurveys' or 'meanBiomass'
                             ...)                  

{
  
  if(length(ordersToInclude)==1 & ordersToInclude[1]=='All') {
    ordersToInclude = unique(surveyData$Group)
  }
  
  numUniqueBranches = length(unique(surveyData$PlantFK))
  
  firstFilter = surveyData %>%
    filter(julianday >= jdRange[1], julianday <= jdRange[2]) %>%
    mutate(julianweek = 7*floor(julianday/7) + 4)
  
  effortBySpecies = firstFilter %>%
    group_by(Species) %>%  ## group_by(PlantSpecies)  or group_by(Species)
    summarize(nSurveys = n_distinct(ID))
  
  arthCount = firstFilter %>%
    filter(Length >= minLength, 
           Group %in% ordersToInclude) %>%
    mutate(Quantity2 = ifelse(Quantity > outlierCount, 1, Quantity)) %>% #outlier counts replaced with 1
    group_by(Species) %>%
    summarize(totalCount = sum(Quantity2, na.rm = TRUE),
              numSurveysGTzero = length(unique(ID[Quantity > 0])),
              totalBiomass = sum(Biomass_mg, na.rm = TRUE)) %>% 
    right_join(effortBySpecies, by = 'Species') %>%
    #next line replaces 3 fields with 0 if the totalCount is NA
    mutate_cond(is.na(totalCount), totalCount = 0, numSurveysGTzero = 0, totalBiomass = 0) %>%
    mutate(meanDensity = totalCount/nSurveys,
           fracSurveys = 100*numSurveysGTzero/nSurveys,
           meanBiomass = totalBiomass/nSurveys) %>%
  x  arrange(Species) %>%
    data.frame()
  
  return(arthCount)
}

#Joining so "Species" becomes associated with a sciName
caterpillar_unclean = read.csv('Tasks/caterpillar_plantanalysis.csv')
plants_clean = read.csv('Tasks/plantList_rerun.csv')
cleaned <- left_join(caterpillar_unclean, plants_clean, by= 'Species')

#Created a column with just "Genus" in order to add to "tallamy_shrop...csv"
Genus <- word(cleaned$sciName, 1)
cleanedish <- mutate(cleaned, Genus=Genus)

tallamy = read.csv('data/tallamy_shropshire_2009_plant_genera.csv')
clean_and_tallamy <- left_join(tallamy, cleanedish, by = 'Genus')
clean_and_tallamy <- select(clean_and_tallamy, -"X.1", -"X.2", -"X.x")
