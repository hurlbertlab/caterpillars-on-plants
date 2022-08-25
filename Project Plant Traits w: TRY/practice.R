# Script for summarizing TRY trait data for plant species

library(readxl)
library(dplyr)
library(tidyr)

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
         meanLeafLength = mean.x, ...) %>%
  select(AccSpeciesName, maxLeafLength, ...)

----------------------------------------------------------------------------------
  
  
# Script for summarizing TRY leaf length(s) trait data for plant species for Part 1 and Part 2
  
library(readxl)
library(dplyr)

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
  
# Summarizing other trait data all at once

library(readxl)
library(dplyr)

alltraits_summaries_part1 = part1_fromTRYfromplantspp %>%
  filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole (leaf lamina length)", "Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if petiole is in- or excluded)", "Leaf area (in case of compound leaves: leaf, petiole excluded)", "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", "Leaf area (in case of compound leaves: leaflet, petiole excluded)", "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)")) %>%
      mutate(lengthCM = case_when(
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
  group_by(AccSpeciesName) %>%
  summarize(maxLeafArea = max(StdValue, na.rm = T),
            meanLeafArea = mean(StdValue, na.rm = T),
            minLeafArea = min(StdValue, na.rm = T))

# Example incorporating multiple trait names
filter(TraitName %in% c("Leaf length", "Leaf length excluding petiole"))
----------------------------------------------------------------------------------
#Joining Part 1s and Part 2s
leaflengths = left_join(leaflengthpart1, leaflengthpart2, by = 'AccSpeciesName'