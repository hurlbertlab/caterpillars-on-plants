library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(gridExtra)
library(maps)
library(sp)
#maptools is going to expire end of 2023; do i need this?
library(maptools)
library(vioplot)

sites = read.csv('2022-08-17_Site.csv', row.names = 1)

cleanDatasetCC = read.csv('PlantsToIdentify/JoinedPhotoAndOccurrenceToFull.csv', row.names = 1) %>%
  mutate(sciName = gsub("\xa0", " ", sciName),
         Genus = word(sciName, 1)) %>%
  filter(sciName != Genus)

# Joining official full dataset of plants to Tallamy et al. to get native, introduced, etc. data
# This helps obtain the families that should be analyzed (those with native, introduced species)
tallamy = read.csv('data/Plant Analysis/tallamy_shropshire_2009_plant_genera.csv') %>%
  mutate(Family = trimws(Family..as.listed.by.USDA.)) %>%
  filter(herbaceous.or.woody %in% c("w"))

# Finding the plant families that are alien
alien_families = unique(tallamy$Family[tallamy$origin..for.analysis. == "alien"])

# Left_join Caterpillars Count! data with Tallamy (alien/native) genus list
cc_plus_tallamy <- left_join(cleanDatasetCC, tallamy, by = 'Genus') %>%
  dplyr::select(ID, PlantFK:ObservationMethod, PlantSpecies:AverageLeafLength, Group:Biomass_mg, 
                sciName:Genus,Family, origin..for.analysis., total.Lep.spp) %>%
  dplyr::rename(origin = origin..for.analysis., lepS = total.Lep.spp)

write.csv(cc_plus_tallamy, 'data/Plant Analysis/cc_plus_tallamy.csv', row.names = F)

# Added because "mutate_cond" is not a built-in function
mutate_cond <- function(.data, condition,...,envir=parent.frame()){
  condition <- eval(substitute(condition),.data,envir)
  .data[condition,] <- .data[condition,] %>% mutate(...)
  .data
}

# Looking at meanDensity, Biomass, and fracSurveys of caterpillars in different plant families 
AnalysisBySciName = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                            ordersToInclude = 'All',    # or 'caterpillar' like arthGroup
                            minLength = 0,              # minimum arthropod size to include 
                            jdRange = c(132, 232),      # change range of days
                            outlierCount = 10000,       # Outliers 
                            plotVar = 'Density',        # 'Density' or 'fracSurveys' or 'Biomass'
                            ...)                  
  
{ if(length(ordersToInclude)==1 & ordersToInclude[1]=='All') {
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

# START HERE
cc_plus_tallamy = read.csv(file = "data/Plant Analysis/cc_plus_tallamy.csv")

#for every plant family return a graph created by comparingBugs... for an arthropod

# Rename the function ; incorporating the other function?
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
    
    # Completing a t.test analysis and pulling out the means and p-value
    if(comparisonVar != "fracSurveys") {
      x = log10(nativeData[,comparisonVar] + 0.001)
      y = log10(alienData[,comparisonVar] + 0.001)
    } else {
      x = nativeData[,comparisonVar]
      y = alienData[,comparisonVar]
    }
    
    w = wilcox.test(x, y, exact = FALSE)
    p_value = w$p.value
    
    #Change the column name to meanBiomass, meanDensity, or fraSurveys
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

#text(1, -0.5, "Density")

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


#comparisonVar = meanDensity, biomass, fracSurveys 
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


## Code to calculate the lepS stuff ## something is still off ##
# need meanDensity, etc. to be calculated but that's all cal in the function currently
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

# ancova.test = lm(meanDensity ~ lepS + origin2 + lepS * origin2, data = lepSandAllFam)
# summary(ancova.test)
# separate lines for each part of the ancova values based on the summary

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

