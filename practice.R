# Script for summarizing TRY trait data for plant species

library(readxl)
library(dplyr)

# Read in TRY trait files
try1 = read.csv('TRY/part1-fromTRYfromplantspp.csv')
try2 = read.csv('TRY/part2-fromTRYfromplantspp.csv')

try = rbind(try1, try2)

traitsummary = try %>%
  filter(TraitName != "") %>%
  mutate(NewTraitName = case_when(
    TraitName %in% c("Leaf length", "Leaf length excluding petiole (leaf lamina length)") ~ "LeafLength",
    TraitName %in% c("Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if petiole is in- or excluded)", 
                     "Leaf area (in case of compound leaves: leaf, petiole excluded)", 
                     "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", 
                     "Leaf area (in case of compound leaves: leaflet, petiole excluded)", 
                     "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)") ~ "LeafArea"),
    StdValue = as.numeric(StdValue)) %>%
  group_by(AccSpeciesName, NewTraitName) %>%
  summarize(max = max(StdValue, na.rm = T),
            mean = mean(StdValue, na.rm = T),
            min = min(StdValue, na.rm =T)) 

lengths = filter(traitsummary, NewTraitName == "LeafLength")
areas = filter(traitsummary, NewTraitName == "LeafArea")

traits = left_join(lengths, areas, by = 'AccSpeciesName') %>%
  rename(maxLeafLength = max.x,
         meanLeafLength = mean.x,
         minLeafLength = min.x) %>%
  select(AccSpeciesName, maxLeafLength, meanLeafLength, minLeafLength, NewTraitName.x)
# what does select() do exactly? grabs a subset but why?
# Example incorporating multiple trait names
filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole (leaf lamina length)"))



----------------------------------------------------------------------------------
  
  
# Script for summarizing TRY leaf length(s) trait data for plant species for Part 1 and Part 2
  
library(readxl)
library(dplyr)

part2_fromTRYfromplantspp = read.csv(part2_fromTRYfromplantspp.csv)

leaflengthnopetiolepart2 = part2_fromTRYfromplantspp %>%
  filter(TraitName == "Leaf length excluding petiole (leaf lamina length)") %>%
  #  mutate(lengthCM = case_when(
  #    OrigUnitStr == "cm" ~ OrigValueStr,
  #    OrigUnitStr == "mm" ~ OrigValueStr/10,
  #    OrigUnitStr == "m" ~ OrigValueStr*100
  #  )) %>%
  mutate(lengthCM = as.numeric(OrigValueStr)) %>%
  group_by(AccSpeciesName) %>%
  summarize(maxLength = max(lengthCM, na.rm = T),
            meanLength = mean(lengthCM, na.rm = T),
            minLength = min(lengthCM, na.rm = T))

----------------------------------------------------------------------------------
  
# Summarizing other trait data all leaf length

library(readxl)
library(dplyr)

alltraits_summaries_part1 = part1_fromTRYfromplantspp %>%
  filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole (leaf lamina length)", ") %>%
#      mutate(lengthCM = case_when(
#      OrigUnitStr == "cm" ~ StdValue,
#      OrigUnitStr == "cm2" ~ StdValue
#      OrigUnitStr == "mm" ~ StdValue/10,
#      OrigUnitStr == "mm2" ~ StdValue/100,
#    )) %>%
  mutate(lengthCM = as.numeric(StdValue)) %>%
  group_by(AccSpeciesName) %>%
  summarize(maxLength = max(lengthCM, na.rm = T),
            meanLength = mean(lengthCM, na.rm = T),
            minLength = min(lengthCM, na.rm = T))

# Example incorporating multiple trait names
filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole"))
----------------------------------------------------------------------------------
  
# Summarizing leaf area data from Part 1 and Part 2
library(readxl)
library(dplyr)

alltraits_summaries_part1 = part1_fromTRYfromplantspp %>%
  filter(TraitName %in% c("Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if petiole is in- or excluded)", 
  "Leaf area (in case of compound leaves: leaf, petiole excluded)", 
  "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", 
  "Leaf area (in case of compound leaves: leaflet, petiole excluded)", 
  "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)")) %>%
  mutate(lengthCM = case_when(
    # when the original value (which is in the units in OrigUnitStr) is in mm2 then replace it with the orginal value divided by 100 = cm2
    UnitName == "mm2" ~ StdValue/100,
  )) %>%
  mutate(lengthCM = as.numeric(StdValue)) %>%
  group_by(AccSpeciesName) %>%
  summarize(maxLength = max(lengthCM, na.rm = T),
            meanLength = mean(lengthCM, na.rm = T),
            minLength = min(lengthCM, na.rm = T))

# Example incorporating multiple trait names
filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole"))

----------------------------------------------------------------------------------
# Joining Part 1s and Part 2s
combinedleaflengthtraits <-rbind(leaflengthpart1, leaflengthpart2)

#Running this new combined lists through ITIS?
new = read.csv('Tasks/cleanedListwithNewSpecies.csv')
coniferList = unique(new[, c('sciName', 'isConifer')])

##Comparing the plantList generated by ITIS and the new species list##
plantList = read.csv('Tasks/plantList_rerun.csv')
newSpeciesList = read.csv('Tasks/newSpecieslist.csv')
adddingplantsandnew <- bind_rows(plantList, newSpeciesList)  

#how many time a value occurs in the combined list to see what is species aren't new
n_occur <-data.frame(table(adddingplantsandnew$sciName))

                        
                        
                        