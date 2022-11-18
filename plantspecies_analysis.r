library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(gridExtra)
library(maps)
library(sp)
library(maptools)

# Added because "mutate_cond" is not a built-in function
mutate_cond <- function(.data, condition,...,envir=parent.frame()){
  condition <- eval(substitute(condition),.data,envir)
  .data[condition,] <- .data[condition,] %>% mutate(...)
  .data
}

# Looking at meanDensity, meanBiomass, and fracSurveys of caterpillars in different plant families 
# over the course of the citizen science dataset from June to July
meanDensityBySciName = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
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
  
  effortBySciName = firstFilter %>%
    group_by(sciName) %>% 
    summarize(nSurveys = n_distinct(ID))
  
  arthCount = firstFilter %>%
    filter(Length >= minLength, 
           Group %in% ordersToInclude) %>%
    mutate(Quantity2 = ifelse(Quantity > outlierCount, 1, Quantity)) %>% #outlier counts replaced with 1
    group_by(sciName) %>%
    summarize(totalCount = sum(Quantity2, na.rm = TRUE),
              numSurveysGTzero = length(unique(ID[Quantity > 0])),
              totalBiomass = sum(Biomass_mg, na.rm = TRUE)) %>% 
    right_join(effortBySciName, by = 'sciName') %>%
    #next line replaces 3 fields with 0 if the totalCount is NA
    mutate_cond(is.na(totalCount), totalCount = 0, numSurveysGTzero = 0, totalBiomass = 0) %>%
    mutate(meanDensity = totalCount/nSurveys,
           fracSurveys = 100*numSurveysGTzero/nSurveys,
           meanBiomass = totalBiomass/nSurveys) %>%
    arrange(sciName) %>%
    data.frame()
  
  return(arthCount)
}

cleanDataset = read.csv('../caterpillars-on-plants/PlantsToIdentify/JoinedPhotoAndOccurrenceToFull.csv', row.names = 1) 
cleanDataset$sciName = gsub("\xa0", " ", cleanDataset$sciName)


# Finding the caterpillars surveyed in the months of June and July  
plantCountJuneJuly = cleanDataset %>%
  dplyr::filter(julianday >= 152, julianday <= 252) %>% #change range of days 
  distinct(ID, sciName) %>%
  count(sciName) %>%
  arrange(desc(n))

# Specifies that only plant species that were surveyed at least 10x in June and July were included
SurveyedCertainAmount = cleanDataset %>%
  filter(sciName %in% plantCountJuneJuly$sciName[plantCountJuneJuly$n >= 10])

# Specifics that only caterpillars (not all arthropods) were analyzed in this analysis
onlyCaterpillars = meanDensityBySciName(SurveyedCertainAmount, ordersToInclude = "caterpillar") %>%
  mutate(Genus = word(sciName, 1)) 

# Joining official full dataset of plants to Tallamy et al. to get native, introduced, etc. data
# This helps obtain the families that should be analyzed (those with native, introduced species)
tallamy = read.csv('data/Plant Analysis/tallamy_shropshire_2009_plant_genera.csv') %>%
  mutate(Family = trimws(Family..as.listed.by.USDA.))

# finding the plant families that are alien
alien_families = unique(tallamy$Family[tallamy$origin..for.analysis. == "alien"])

# left_join SurveyWithCaterpillar before
clean_and_tallamy <- left_join(onlyCaterpillars, tallamy, by = 'Genus') %>%
  select(sciName:Genus,Family, origin..for.analysis., total.Lep.spp, nSurveys, meanDensity, fracSurveys, meanBiomass) %>%
  rename(origin = origin..for.analysis., lepS = total.Lep.spp) %>%
###finding Families that are in both native and alien categories?  
  filter(Family %in% alien_families) %>%
  arrange(Family, origin) #%>%

# Compare origin to native and origin to alien species and examining arthropod meanDensity, meanBiomass, and fracSurveys 
### is giving all the entries.. is this right
nativeData = filter(clean_and_tallamy, origin == 'native') 
alienData = filter(clean_and_tallamy, origin == 'alien')


# A pdf with graphs depicting density, biomass, % surveyed
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/allSpecies.pdf",
    width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1))

t.test(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001))
boxplot(log10(nativeData$meanDensity + 0.001), log10(alienData$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "All species", width = c(0.5, 0.5), ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 96", "N = 18"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[1], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.34, "p = 0.005")


t.test(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001))
boxplot(log10(nativeData$meanBiomass + 0.001), log10(alienData$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "All species", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 96", "N = 18"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[2], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 1.8, "p = 0.148")


t.test(nativeData$fracSurveys, alienData$fracSurveys)
boxplot(nativeData$fracSurveys, alienData$fracSurveys, 
        xaxt = 'n', las = 1, main = "All species", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 96", "N = 18"), 1, at = 1:2, line = 2, cex = 0.75)
#mtext(text=LETTERS[3], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 50, "p = 0.002")

dev.off()


# Multi-panel plots comparing native vs alien species for 5 plant families for meanDensity, meanBiomass, and fracSurveys
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/6families.pdf", 
    width = 8.5, height = 15)
par(mfrow = c(7, 3), mar = c(5, 5, 2, 2))

#Rosaceae family comparison of meanDensity, meanBiomass, and fracSurveys
rosaceaeNative = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == "native")
rosaceaeAlien = dplyr::filter(clean_and_tallamy, Family == "Rosaceae", origin == 'alien')

t.test(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001))
boxplot(log10(rosaceaeNative$meanDensity + 0.001), log10(rosaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 29", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[1], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.3, "p = 0.529")


t.test(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001))
boxplot(log10(rosaceaeNative$meanBiomass + 0.001), log10(rosaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 29", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = 0.817")


t.test(rosaceaeNative$fracSurveys, rosaceaeAlien$fracSurveys)
boxplot(rosaceaeNative$fracSurveys, rosaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Rosaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 29", "N = 3"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 45, "p = 0.024")


#Ericaceae family comparison of meanDensity, meanBiomass, and fracSurveys
ericaceaeNative = filter(clean_and_tallamy, Family == "Ericaceae", origin == "native")
ericaceaeAlien = filter(clean_and_tallamy, Family == "Ericaceae", origin == 'alien')

t.test(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001))
boxplot(log10(ericaceaeNative$meanDensity + 0.001), log10(ericaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 14", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[2], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(2, 0.3, "p = 2.2e-16")


t.test(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ericaceaeNative$meanBiomass + 0.001), log10(ericaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 14", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 2.2e-16")


t.test(ericaceaeNative$fracSurveys, ericaceaeAlien$fracSurveys)
boxplot(ericaceaeNative$fracSurveys, ericaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Ericaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 14", "N = 1"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 2.2e-16")



#Moraceae family comparison of meanDensity, meanBiomass, and fracSurveys
moraceaeNative = filter(clean_and_tallamy, Family == "Moraceae", origin == "native")
moraceaeAlien = filter(clean_and_tallamy, Family == "Moraceae", origin == 'alien')

t.test(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001))
boxplot(log10(moraceaeNative$meanDensity + 0.001), log10(moraceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 3", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[3], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(2, 0.3, "p = 2.2e-16")


t.test(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001))
boxplot(log10(moraceaeNative$meanBiomass + 0.001), log10(moraceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 3", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 2.2e-16")


t.test(moraceaeNative$fracSurveys, moraceaeAlien$fracSurveys)
boxplot(moraceaeNative$fracSurveys, moraceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Moraceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 3", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
#text(2, 0.3, "p = 2.2e-16")



#Oleaceae family comparison of meanDensity, meanBiomass, and fracSurveys
oleaceaeNative = filter(clean_and_tallamy, Family == "Oleaceae", origin == "native")
oleaceaeAlien = filter(clean_and_tallamy, Family == "Oleaceae", origin == 'alien')

t.test(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001))
boxplot(log10(oleaceaeNative$meanDensity + 0.001), log10(oleaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(c("p = 0.050"), 1, at = 1, line = 0.05, cex = 0.5)
mtext(text=LETTERS[4], xpd=NA, side=1, adj=0, font=2, cex=0.75)
#text(x=1, y=1, labels="p = 2.2e-16", cex =5)
#it's there but won't move down ####

t.test(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001))
boxplot(log10(oleaceaeNative$meanBiomass + 0.001), log10(oleaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = 0.358")


t.test(oleaceaeNative$fracSurveys, oleaceaeAlien$fracSurveys)
boxplot(oleaceaeNative$fracSurveys, oleaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Oleaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 6", "N = 2"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 15, "p = 0.161")


####
#Ulmaceae family comparison of meanDensity, meanBiomass, and fracSurveys
ulmaceaeNative = filter(clean_and_tallamy, Family == "Ulmaceae", origin == "native")
ulmaceaeAlien = filter(clean_and_tallamy, Family == "Ulmaceae", origin == 'alien')

t.test(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001))
boxplot(log10(ulmaceaeNative$meanDensity + 0.001), log10(ulmaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Berberidaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[5], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.0099, "p = constant")


t.test(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001))
boxplot(log10(ulmaceaeNative$meanBiomass + 0.001), log10(ulmaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Berberidaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = constant")


t.test(ulmaceaeNative$fracSurveys, ulmaceaeAlien$fracSurveys)
boxplot(ulmaceaeNative$fracSurveys, ulmaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Berberidaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 5.8, "p = constant")

#Theaceae family comparison of meanDensity, meanBiomass, and fracSurveys
theaceaeNative = filter(clean_and_tallamy, Family == "Theaceae", origin == "native")
theaceaeAlien = filter(clean_and_tallamy, Family == "Theaceae", origin == 'alien')

t.test(log10(theaceaeNative$meanDensity + 0.001), log10(theaceaeAlien$meanDensity + 0.001))
boxplot(log10(theaceaeNative$meanDensity + 0.001), log10(theaceaeAlien$meanDensity + 0.001), 
        xaxt = 'n', las = 1, main = "Theaceae", boxwex = 0.5, ylab = "log(Density)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
mtext(text=LETTERS[5], xpd=NA, side=1, adj=0, font=2, cex=0.75)
text(2, 0.0099, "p = constant")


t.test(log10(theaceaeNative$meanBiomass + 0.001), log10(theaceaeAlien$meanBiomass + 0.001))
boxplot(log10(theaceaeNative$meanBiomass + 0.001), log10(theaceaeAlien$meanBiomass + 0.001), 
        xaxt = 'n', las = 1, main = "Theaceae", boxwex = 0.5, ylab = "log(Biomass)", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 0.3, "p = constant")


t.test(theaceaeNative$fracSurveys, theaceaeAlien$fracSurveys)
boxplot(theaceaeNative$fracSurveys, theaceaeAlien$fracSurveys, 
        xaxt = 'n', las = 1, main = "Theaceae", boxwex = 0.5, ylab = "% of Surveys", col = c("burlywood", "rosybrown"))
mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
mtext(c("N = 5", "N = 0"), 1, at = 1:2, line = 2, cex = 0.75)
text(2, 5.8, "p = constant")


dev.off()

# Comparing the average caterpillar density, biomass, and fracSurveys per survey to lepS and conducting a linear regression
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/Lepidoptera.pdf", 
    width = 9, height = 6)
par(mfrow = c(2, 2), mar = c(5,5,2,1))

plot(clean_and_tallamy$lepS, clean_and_tallamy$meanDensity, xlab = "Lepidoptera per survey", ylab = "Density", pch = 16, main ="Mean Density")
text(140, 2.60, "R2 =0.067, p = 0.006")
#mtext(text=LETTERS[1], xpd=NA, side=3, adj=0, font=2)
lm.density = lm(meanDensity ~ lepS, data = clean_and_tallamy)
summary(lm.density)
abline(lm.density)


plot(log10(clean_and_tallamy$lepS), log10(clean_and_tallamy$meanBiomass), xlab = "log(Lepidoptera per survey)", ylab = "log(Biomass)", pch = 16, main ="Mean Biomass")
text(0.75, 2, "R2 = 4.25e-05, p = 0.945", cex = 0.85)
#mtext(text=LETTERS[2], xpd=NA, side=3, adj=0, font=2)
lm.biomass = lm(meanBiomass ~ lepS, data = clean_and_tallamy)
summary(lm.biomass)
abline(lm.biomass)


plot(clean_and_tallamy$lepS, clean_and_tallamy$fracSurveys, xlab = "Lepidoptera per survey", ylab = "Lepidoptera Richness", pch = 16, main ="% of Surveys")
text(140, 50, "R2 = .069, p = 0.005", cex = 0.85)
#mtext(text=LETTERS[3], xpd=NA, side=3, adj=0, font=2)
lm.surveys = lm(fracSurveys ~ lepS, data = clean_and_tallamy)
summary(lm.surveys)
abline(lm.surveys)

dev.off()


# Creating a summary table of # of species, etc.
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/Summary.pdf", 
    height=11.5, width=8)


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

