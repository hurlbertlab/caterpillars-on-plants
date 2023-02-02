library(tidyverse)
library(lubridate)
library(data.table)
library(gsheet)
library(gridExtra)
library(maps)
library(sp)
#maptools is going to expire end of 2023; do i need this?
library(maptools)

cleanDatasetCC = read.csv('../caterpillars-on-plants/PlantsToIdentify/JoinedPhotoAndOccurrenceToFull.csv', row.names = 1) %>%
  mutate(sciName = gsub("\xa0", " ", sciName),
         Genus = word(sciName, 1))

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
                            jdRange = c(152, 212),      # change range of days
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
                            jdRange = c(152, 252),                   # Range of days
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
             sciName %in% plantCount$sciName[plantCount$n >= minSurveysPerPlant])  
    
    onlyBugs = AnalysisBySciName(filteredData, ordersToInclude = arthGroup) %>%
      left_join(plantCount, by = "sciName")
    
    # Separating datasets
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
    
    t = t.test(x, y)
    nativeMean = t$estimate[1]
    alienMean = t$estimate[2]
    p_value = t$p.value
    native_pop_size = nrow(nativeData)
    alien_pop_size = nrow(alienData)
    
    # Creating a pdf for each arthGroup
   # for (j in cc_plus_tallamy$Group) {
    #  for (i in familiesWithNativeAndAlienSpecies$Family) {
        
    # Plotting the analysis
    if(plot == TRUE) {
      plot_title = plantFamily
      y_label = comparisonVar
      #if (comparisonVar = "meanDensity" || comparisonVar = "meanBiomass" || comparisonVar = "fracSurveys") {
       # y_label = "Density" 
        #|| y_label = "Biomass" || y_label = "% Surveyed"
     # }
      
      boxplot(x, y, xaxt = 'n', las = 1, main = paste(plot_title, ", p =", round(p_value,3)),
              boxwex = 0.5, ylab = y_label, col = c("burlywood", "rosybrown"))
      mtext(c("Native", "Alien"), 1, at = 1:2, line = 1)
      mtext(c(paste("N =", native_pop_size), paste("N =", alien_pop_size)), 1, at = 1:2, line = 2, cex = 0.75)
      mtext(text=LETTERS[1], xpd=NA, side=2, adj=0, font=2, cex=0.75)
    }
  }
}

# A pdf with graphs depicting density, biomass, % surveyed
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/ByArthByPlantFam.pdf",
    width = 15, height = 5)
par(mfrow = c(4, 3), mar = c(3, 3, 3, 1))

#creating a vector list for arthGroup and the specific families, can run the familiesWith... group after the function

    
for (group in c("caterpillar", "beetle", "truebugs", "spider")) {
  
  for (plotVar in c("meanDensity", "meanBiomass", "fracSurveys")) {
    
    comparingBugsonNativeVersusAlienPlants(cc_plus_tallamy, plantFamily = "Rosaceae", arthGroup = group, comparisonVar = plotVar, plot = TRUE)
  }
}

    
dev.off()



# Comparing the average caterpillar density, biomass, and fracSurveys per survey to lepS and conducting a linear regression
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/Lepidoptera.pdf", 
    width = 9, height = 6)
par(mfrow = c(2, 2), mar = c(5,5,2,1))

plot(clean_and_tallamy$lepS, clean_and_tallamy$meanDensity, xlab = "Lepidoptera Richness", ylab = "Density", pch = 16, main ="Mean Density")
text(140, 2.60, "R2 =0.067, p = 0.006")
mtext(text=LETTERS[1], xpd=NA, side=3, adj=0, font=2)
lm.density = lm(meanDensity ~ lepS, data = clean_and_tallamy)
summary(lm.density)
abline(lm.density)


plot(log10(clean_and_tallamy$lepS), log10(clean_and_tallamy$meanBiomass), xlab = "Lepidoptera Richness", ylab = "log(Biomass)", pch = 16, main ="Mean Biomass")
text(0.75, 2, "R2 = 4.25e-05, p = 0.945", cex = 0.85)
mtext(text=LETTERS[2], xpd=NA, side=3, adj=0, font=2)
lm.biomass = lm(meanBiomass ~ lepS, data = clean_and_tallamy)
summary(lm.biomass)
abline(lm.biomass)


plot(clean_and_tallamy$lepS, clean_and_tallamy$fracSurveys, xlab = "Lepidoptera Richness", ylab = "Lepidoptera", pch = 16, main ="% of Surveys")
text(140, 50, "R2 = 0.069, p = 0.005", cex = 0.85)
mtext(text=LETTERS[3], xpd=NA, side=3, adj=0, font=2)
lm.surveys = lm(fracSurveys ~ lepS, data = clean_and_tallamy)
summary(lm.surveys)
abline(lm.surveys)

dev.off()


# Creating a summary table of the number of species, etc.
pdf(file = "/Users/colleenwhitener/Documents/2-Junior Year/1-BIOL 395/caterpillars-on-plants/Figures/Summary.pdf", 
    height=11.5, width=8)

summary_table <- matrix(c(114, 94, 18, 3, 29, 6, 7, 1, 1, 1, 13, 1, 4, 1, 0))
colnames(summary_table) <- c("Sum")
rownames(summary_table) <- c("Alien + Native","Native","Alien", 
                             "Rosaceae Alien", "Rosaceae Native",
                             "Oleaceae Alien", "Oleaceae Native",
                             "Berberidaceae Alien", "Berberidaceae Native",
                             "Ericaceae Alien", "Ericaceae Native",
                             "Moraceae Alien", "Moraceae Native", 
                             "Theaceae Alien", "Theaceae Native")
summary_table <- as.table(summary_table)
grid.table(summary_table)
t(summary_table)

dev.off()


