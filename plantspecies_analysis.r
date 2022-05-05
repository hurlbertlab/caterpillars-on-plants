library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(gridExtra)
library(maps)
library(sp)
library(maptools)

#Added bc "mutate_cond" is not a built-in function
mutate_cond <- function(.data, condition,...,envir=parent.frame()){
  condition <- eval(substitute(condition),.data,envir)
  .data[condition,] <- .data[condition,] %>% mutate(...)
  .data
}

#Looking at meanDensity, meanBiomass, and fracSurveys of caterpillars in different plant families over the course of the citizen science dataset from June to July

meanDensityBySpecies = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                             ordersToInclude = 'All',       # or 'caterpillar'
                             minLength = 0,         # minimum arthropod size to include 
                             jdRange = c(152, 212), #change range of days
                             outlierCount = 10000,
                             plotVar = 'meanDensity', # 'meanDensity' or 'fracSurveys' or 'meanBiomass'
                             ...)                  

{
  
  if(length(ordersToInclude)==1 & ordersToInclude[1]=='All') {
    ordersToInclude = unique(surveyData$Group)
  }
  
  numUniqueBranches = length(unique(surveyData$PlantFK))
  
  firstFilter = surveyData %>%
    filter(julianday >= jdRange[1], julianday <= jdRange[2]) %>% #subscription notation
    mutate(julianweek = 7*floor(julianday/7) + 4)
  
  effortBySpecies = firstFilter %>%
    group_by(Species) %>% 
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


cleanDataset = read.csv('../caterpillars-count-data/dataCleaning/flagged_dataset_2022-01-27.csv') %>%
  filter(status != "remove")

plantCountJuneJuly = cleanDataset %>%
  filter(julianday >= 152, julianday <= 252) %>% #change range of days 
  distinct(ID, Species) %>%
  count(Species) %>%
  arrange(desc(n))

# Specifies that only plant species that were surveyed at least 10x in June and July were included
filteredData = cleanDataset %>%
  filter(Species %in% plantCountJuneJuly$Species[plantCountJuneJuly$n >= 10])

# Specifics that only caterpillars (not all arthropods) were analyzed in this analysis
caterpillar_unclean = meanDensityBySpecies(filteredData, ordersToInclude = "caterpillar")

write.csv(caterpillar_unclean, 'data/Plant Analysis/caterpillar_plantanalysis.csv', row.names = F)

#Joining so the "Species" column becomes associated with a sciName
caterpillar_unclean = read.csv('data/Plant Analysis/caterpillar_plantanalysis.csv')

plants_clean = read.csv('data/Plant Analysis/plantList_rerun.csv') %>%
  select(-X)

cleaned <- left_join(caterpillar_unclean, plants_clean, by= 'Species')
cleaned.new <- cleaned %>%
  mutate(Genus = word(cleaned$sciName, 1)) #Created a column with just "Genus" in order to add to "tallamy_shrop...csv"

tallamy = read.csv('data/tallamy_shropshire_2009_plant_genera.csv') %>%
  mutate(Family = trimws(Family..as.listed.by.USDA.))

alien_families = unique(tallamy$Family[tallamy$origin..for.analysis. == "alien"])

clean_and_tallamy <- left_join(cleaned.new, tallamy, by = 'Genus') %>%
  select(Species:Genus,Family, origin..for.analysis., total.Lep.spp) %>%
  rename(origin = origin..for.analysis., lepS = total.Lep.spp) %>%
  filter(Family %in% alien_families) %>%
  arrange(Family, origin)



# Compare origin to native and origin to alien species and examining arthropod meanDensity, meanBiomass, and fracSurveys 

nativeData = filter(clean_and_tallamy, origin == 'native')
alienData = filter(clean_and_tallamy, origin == 'alien')

pdf("allSpecies.pdf", width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1))

t.test(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001))
boxplot(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "All species", width = c(0.5, 0.5), ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 138", "N = 27"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[1], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.34, "p = 0.004")


t.test(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001))
boxplot(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "All species", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 138", "N = 27"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[2], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 1.8, "p = 0.0001")


t.test(nativeData$fracSurveys, alienData$fracSurveys)
boxplot(nativeData$fracSurveys, alienData$fracSurveys, 
        xaxt = 'n', las = 1, main = "All species", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 138", "N = 27"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[3], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 50, "p = 0.033")

dev.off()


# Multi-panel plots comparing native vs alien species for 5 plant families for meanDensity, meanBiomass, and fracSurveys

pdf("6families.pdf", width = 8.5, height = 12)
par(mfrow = c(5, 3), mar = c(5, 5, 2, 2))

#Rosaceae family comparison of meanDensity, meanBiomass, and fracSurveys
rosaceaeNative = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == "native")
rosaceaeAlien = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == 'alien')

t.test(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001))
boxplot(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 42", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[1], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.3, "p = 0.853")


t.test(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001))
boxplot(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 42", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = 0.797")


t.test(rosaceaeNative$fracSurveys, rosaceaeAlien$fracSurveys)
boxplot(rosaceaeNative$fracSurveys, rosaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 42", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 45, "p = 0.896")


#Ericaceae family comparison of meanDensity, meanBiomass, and fracSurveys
ericaceaeNative = filter(clean_and_tallamy, Family == "Ericaceae", origin == "native")
ericaceaeAlien = filter(clean_and_tallamy, Family == "Ericaceae", origin == 'alien')

t.test(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001))
boxplot(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 15", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[2], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(2, 0.3, "p = 0.003")


t.test(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 15", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 0.003")


t.test(ericaceaeNative$fracSurveys, ericaceaeAlien$fracSurveys)
boxplot(ericaceaeNative$fracSurveys, ericaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 15", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 0.003")



#Moraceae family comparison of meanDensity, meanBiomass, and fracSurveys
moraceaeNative = filter(clean_and_tallamy, Family == "Moraceae", origin == "native")
moraceaeAlien = filter(clean_and_tallamy, Family == "Moraceae", origin == 'alien')

t.test(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001))
boxplot(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 4", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[3], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(2, 0.3, "p = 0.003")


t.test(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001))
boxplot(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 4", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 0.003")


t.test(moraceaeNative$fracSurveys, moraceaeAlien$fracSurveys)
boxplot(moraceaeNative$fracSurveys, moraceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 4", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 0.003")



#Oleaceae family comparison of meanDensity, meanBiomass, and fracSurveys
oleaceaeNative = filter(clean_and_tallamy, Family == "Oleaceae", origin == "native")
oleaceaeAlien = filter(clean_and_tallamy, Family == "Oleaceae", origin == 'alien')

t.test(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001))
boxplot(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 8"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(c("p = 0.006"), 1, at = 1, line = 0.05, cex = 0.5)
mtext(text=LETTERS[4], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(x=1, y=1, labels="p = 0.006", cex =5)
#it's there but won't move down ####

t.test(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001))
boxplot(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 8"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = 0.00005")


t.test(oleaceaeNative$fracSurveys, oleaceaeAlien$fracSurveys)
boxplot(oleaceaeNative$fracSurveys, oleaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 8"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 15, "p = 0.155")



#Ulmaceae family comparison of meanDensity, meanBiomass, and fracSurveys
ulmaceaeNative = filter(clean_and_tallamy, Family == "Ulmaceae", origin == "native")
ulmaceaeAlien = filter(clean_and_tallamy, Family == "Ulmaceae", origin == 'alien')

t.test(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001))
boxplot(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Ulmaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 9", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[5], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.0099, "p = 0.001")


t.test(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Ulmaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 9", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = 0.001")


t.test(ulmaceaeNative$fracSurveys, ulmaceaeAlien$fracSurveys)
boxplot(ulmaceaeNative$fracSurveys, ulmaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Ulmaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 9", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 5.8, "p = 0.041")


dev.off()

# Comparing the average caterpillar density, biomass, and fracSurveys per survey to lepS and conducting a linear regression
pdf("Lepidoptera.pdf", width = 9, height = 6)
par(mfrow = c(2, 2), mar = c(5,5,2,1))

plot(clean_and_tallamy$lepS, clean_and_tallamy$meanDensity, xlab = "Lepidoptera per survey", ylab = "Density", pch = 16, main ="Mean Density")
text(140, 2.60, "R2 = 0.065, p = 0.001")
#mtext(text=LETTERS[1], xpd=NA, side=3, adj=0, font=2)
lm.density = lm(meanDensity ~ lepS, data = clean_and_tallamy)
summary(lm.density)
abline(lm.density)


plot(log10(clean_and_tallamy$lepS), log10(clean_and_tallamy$meanBiomass), xlab = "log(Lepidoptera per survey)", ylab = "log(Biomass)", pch = 16, main ="Mean Biomass")
text(0.75, 2, "R2 = 0.003, p = 0.469", cex = 0.85)
#mtext(text=LETTERS[2], xpd=NA, side=3, adj=0, font=2)
lm.biomass = lm(meanBiomass ~ lepS, data = clean_and_tallamy)
summary(lm.biomass)
abline(lm.biomass)


plot(clean_and_tallamy$lepS, clean_and_tallamy$fracSurveys, xlab = "Lepidoptera per survey", ylab = "Lepidoptera Richness", pch = 16, main ="% of Surveys")
text(140, 50, "R2 = 0.041, p = 0.009", cex = 0.85)
#mtext(text=LETTERS[3], xpd=NA, side=3, adj=0, font=2)
lm.surveys = lm(fracSurveys ~ lepS, data = clean_and_tallamy)
summary(lm.surveys)
abline(lm.surveys)

dev.off()


# Creating a summary table of # of species, etc.
pdf("Summary.pdf", height=11.5, width=8)


n <- grep("native", clean_and_tallamy$origin)
a <- grep("alien", clean_and_tallamy$origin)
total <- grep("native|alien", clean_and_tallamy$origin)

length(n)
length(a)
length(total)

summary_table <- data.frame("Summary of Species Totals"=c(length(total), length(n), length(a)))
rownames(summary_table) <- c("Alien + Native","Native","Alien")
grid.table(summary_table)

dev.off()
