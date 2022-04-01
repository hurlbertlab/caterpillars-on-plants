# Script for summarizing TRY trait data for plant species

library(readxl)
library(dplyr)

# Read in TRY trait files
try1 = read.csv('TRY/part1-fromTRYfromplantspp.csv')
try2 = read.csv('TRY/part2-fromTRYfromplantspp.csv')
masterplants = read.csv('Tasks/plantList_rerun.csv')
traits = read.csv('TRY/traits.csv')

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
         minLeafLength = min.x,
         LeafLength = NewTraitName.x,
         maxLeafArea = max.y,
         meanLeafArea = mean.y,
         minLeafArea = min.y,
         LeafArea = NewTraitName.y)

# Possibly joining the master plant list to the trait list to see what plants are left
joined_try = left_join(
  masterplants, 
  traits, 
  # tells R which columns should match as c('firstframe' = 'secondframe')
  by = c('sciName' = 'AccSpeciesName'))

# check the number of rows where LeafLength is NA
joined_try %>% 
  summarize(n_missing = sum(is.na(LeafLength)))

# check the number of rows where LeafLength is not NA (where we have data)
joined_try %>% 
  summarize(n_present = sum(!is.na(LeafLength)))

#%>% 
#  filter(!AccSpeciesName %in% 8:9)

# index to view sciName where sciName is duplicated (in the list more than one time)
joined_try$sciName[which(duplicated(joined_try$sciName))]

# index to whole rows where sciName is 'Carpinus caroliniana'
joined_try[joined_try$sciName == 'Carpinus caroliniana',]

