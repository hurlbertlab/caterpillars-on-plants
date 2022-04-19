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
    arrange(Species) %>%
    data.frame()
  
  return(arthCount)
}


fullDataset = read.csv('../caterpillars-analysis-public/data/fullDataset_2022-03-22.csv')
caterpillar_unclean = meanDensityBySpecies(fullDataset, ordersToInclude = "caterpillar")
write.csv(caterpillar_unclean, 'data/Plant Analysis/caterpillar_plantanalysis.csv', row.names = F)


#Joining so the "Species" column becomes associated with a sciName
caterpillar_unclean = read.csv('data/Plant Analysis/caterpillar_plantanalysis.csv') %>%
  select(-X)

plants_clean = read.csv('data/Plant Analysis/plantList_rerun.csv') %>%
  select(-X)

cleaned <- left_join(caterpillar_unclean, plants_clean, by= 'Species') %>%
  mutate(Genus = word(cleaned$sciName, 1)) #Created a column with just "Genus" in order to add to "tallamy_shrop...csv"

tallamy = read.csv('data/tallamy_shropshire_2009_plant_genera.csv') %>%
  mutate(Family = trimws(Family..as.listed.by.USDA.))

alien_families = unique(tallamy$Family[tallamy$origin..for.analysis. == "alien"])

clean_and_tallamy <- left_join(cleaned, tallamy, by = 'Genus') %>%
  select(Species:Genus,Family, origin..for.analysis., total.Lep.spp) %>%
  rename(origin = origin..for.analysis., lepS = total.Lep.spp) %>%
  filter(Family %in% alien_families) %>%
  arrange(Family, origin)



# Compare origin = native to origin = alien
# use log10 transformation because of skew in the distributions
# Adding 0.001 to all values because there are many 0 values for which we can't calculate a log.
# Should maybe revisit the decision to use 0.001 versus some other constant.
## In addition to comparing meanDensity, you can compare meanBiomass and fracSurveys)

nativeData = filter(clean_and_tallamy, origin == 'native')
alienData = filter(clean_and_tallamy, origin == 'alien')

pdf("6families_allspecies.pdf", width = 8.5, height = 11)
par(mfrow = c(6, 3), mar = c(3, 3, 1, 1))

t.test(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001))
boxplot(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "All species - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.34, "p = 0.003")


t.test(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001))
boxplot(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "All species - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 1.25, "p = 0.0001")


t.test(log10(nativeData$fracSurveys + 0.001), log10(alienData$fracSurveys + 0.001))
boxplot(log10(nativeData$fracSurveys + 0.001), log10(alienData$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "All species - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2.3, 0.3, "p = 0.002", cex=0.7)

# Compare origin = native vs alien for each of the 5 plant families with both
# Multi-panel plots for meanDensity, meanBiomass, and fracSurveys for 5 families

#Rosaceae family
rosaceaeNative = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == "native")
rosaceaeAlien = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == 'alien')

t.test(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001))
boxplot(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.3, "p = 0.853")


t.test(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001))
boxplot(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.3, "p = 0.797")


t.test(log10(rosaceaeNative$fracSurveys + 0.001), log10(rosaceaeAlien$fracSurveys + 0.001))
boxplot(log10(rosaceaeNative$fracSurveys + 0.001), log10(rosaceaeAlien$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 1.6, "p = 0.798", cex=0.7)



#Ericaceae family
ericaceaeNative = filter(clean_and_tallamy, Family == "Ericaceae", origin == "native")
ericaceaeAlien = filter(clean_and_tallamy, Family == "Ericaceae", origin == 'alien')

t.test(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001))
boxplot(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")


t.test(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")


t.test(log10(ericaceaeNative$fracSurveys + 0.001), log10(ericaceaeAlien$fracSurveys + 0.001))
boxplot(log10(ericaceaeNative$fracSurveys + 0.001), log10(ericaceaeAlien$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")



#Moraceae family
moraceaeNative = filter(clean_and_tallamy, Family == "Moraceae", origin == "native")
moraceaeAlien = filter(clean_and_tallamy, Family == "Moraceae", origin == 'alien')

t.test(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001))
boxplot(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")


t.test(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001))
boxplot(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")


t.test(log10(moraceaeNative$fracSurveys + 0.001), log10(moraceaeAlien$fracSurveys + 0.001))
boxplot(log10(moraceaeNative$fracSurveys + 0.001), log10(moraceaeAlien$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
#text(2, 0.3, "p = 0.003")



#Oleaceae family
oleaceaeNative = filter(clean_and_tallamy, Family == "Oleaceae", origin == "native")
oleaceaeAlien = filter(clean_and_tallamy, Family == "Oleaceae", origin == 'alien')

t.test(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001))
boxplot(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.1, "p = 0.006")


t.test(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001))
boxplot(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.3, "p = 0.00005")


t.test(log10(oleaceaeNative$fracSurveys + 0.001), log10(oleaceaeAlien$fracSurveys + 0.001))
boxplot(log10(oleaceaeNative$fracSurveys + 0.001), log10(oleaceaeAlien$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2.3, 0.5, "p = 0.002", cex=0.7)



#Ulmaceae family
ulmaceaeNative = filter(clean_and_tallamy, Family == "Ulmaceae", origin == "native")
ulmaceaeAlien = filter(clean_and_tallamy, Family == "Ulmaceae", origin == 'alien')

t.test(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001))
boxplot(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Ulmaceae - meanDensity")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.0099, "p = 0.001")


t.test(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Ulmaceae - meanBiomass")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.3, "p = 0.001")


t.test(log10(ulmaceaeNative$fracSurveys + 0.001), log10(ulmaceaeAlien$fracSurveys + 0.001))
boxplot(log10(ulmaceaeNative$fracSurveys + 0.001), log10(ulmaceaeAlien$fracSurveys + 0.001), 
        xaxt = 'n', las = 1, main = "Ulmaceae - fracSurveys")
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
text(2, 0.3, "p = 0.0008")


dev.off()


# Compare average caterpillar density (and biomass, and fracSurveys) per branch to lepS (the # of species ever recorded for a genus).
# also linear regression using lm.density = lm(meanDensity ~ lepS, data = clean_and_tallamy); summary(lm.density); abline(lm.density)

pdf("Lepidoptera.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 1), mar = c(4,4,1,1))

plot(clean_and_tallamy$lepS, clean_and_tallamy$meanDensity, xlab = "Total Lepidoptera richness", ylab = "Lepidoptera per branch", pch = 16)
text(160, 2.60, "p = 0.001, R2 = 0.065")
lm.density = lm(meanDensity ~ lepS, data = clean_and_tallamy)
summary(lm.density)
abline(lm.density)


plot(clean_and_tallamy$lepS, clean_and_tallamy$meanBiomass, xlab = "Total Lepidoptera richness", ylab = "Lepidoptera per branch", pch = 16)
text(320, 130, "p = 0.476, R2 = 0.003")
lm.biomass = lm(meanBiomass ~ lepS, data = clean_and_tallamy)
summary(lm.biomass)
abline(lm.biomass)


plot(clean_and_tallamy$lepS, clean_and_tallamy$fracSurveys, xlab = "Total Lepidoptera richness", ylab = "Lepidoptera per branch", pch = 16)
text(160, 50, "p = 0.009, R2 = 0.041")
lm.surveys = lm(fracSurveys ~ lepS, data = clean_and_tallamy)
summary(lm.surveys)
abline(lm.surveys)

dev.off()

