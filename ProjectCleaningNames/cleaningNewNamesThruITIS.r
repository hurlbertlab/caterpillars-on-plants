# KEY #
# plantName = the user inputted name for a plant (might be the scientific name, different common names, a genus, etc.)
# cleanedPlantName = a version of the plant name that contains only the scientific name, genus or common name without "spp.", "?", etc.
# sciName = the official ITIS (usually) recognized scientific name
# ITIS_ID = the ID number that matches to a specific plant taxonomy in the Integrated Taxonomic Information System

# NOTE: Still need a way to efficiently deal with ambiguous common names that could refer to different taxonomic concepts
#       depending on where you are, e.g. "scrub oak", "ironwood", etc.


library(dplyr)
library(taxize)
library(stringr)


# This function takes a vector of species names and checks each one with ITIS,
# returning a dataframe with the name, sciName, itis_id, and rank.
cleanNamesThruITIS = function(speciesList) {
  
  # speciesList must have a column called plantName
  plantList = data.frame(Species = speciesList,
                         sciName = NA, 
                         itis_id = NA,    
                         rank = NA)
  
  for (i in 1:nrow(plantList)) {
    
    print(paste(i, "of", nrow(plantList), "/n"))
    
    if (is.na(speciesList[i])) { 
      
      plantList$sciName[i] = NA
      plantList$itis_id[i] = NA
      plantList$rank[i] = NA
    
    } else {
    
      hierarchy = classification(speciesList[i], db = 'itis', accepted = TRUE)[[1]]
      
      # class is logical if taxonomic name does not match any existing names
      if (!is.null(nrow(hierarchy))) {
        plantList$sciName[i] = hierarchy$name[nrow(hierarchy)]
        plantList$itis_id[i] = hierarchy$id[nrow(hierarchy)]
        plantList$rank = hierarchy$rank[nrow(hierarchy)]
      } else {
        plantList$sciName[i] = NA
        plantList$itis_id[i] = NA
        plantList$rank[i] = NA
      }    
    }
  }
  return(plantList)
  
}

## NEW WORKFLOW ##
# 1. read in latest Plants.csv
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)

officialPlantListFiles = list.files('ProjectCleaningNames')[str_detect(list.files('ProjectCleaningNames'), 'officialPlantList_')]
mostRecentOfficialPlantList = officialPlantListFiles[length(officialPlantListFiles)]
officialPlantList = read.csv(paste0('ProjectCleaningNames/', mostRecentOfficialPlantList), header = T)


# 2. Find new names not in plantName of officialPLantList
  # Below isn't a true representation of new species bc the names are not "clean" they still have spp., etc.
new_species <- plants %>% 
  rename(plantName = Species) %>%
  distinct(plantName) %>%
  # select rerun sciName entries that are NOT (!) in sciName from cleaned list
  filter(!plantName %in% officialPlantList$plantName) 

write.csv(new_species, paste("ProjectCleaningNames/newSpecies_", Sys.Date(), ".csv", sep = ""), row.names = F)

if (nrow(new_species) > 0) {
# 3. Run new entries through ITIS / Match new names using taxize
cleanedNewNames = cleanNamesThruITIS(new_species$plantName)

# 3.1  Separate out results that did vs did not match in ITIS
# For results that matched, rename "Species" as "plantName", add notes, isConifer, and cleanedName = plantName

matched_new_species <- filter(cleanedNewNames, !is.na(cleanedNewNames$itis_id)) %>%
  rename(plantName = Species) %>%
  mutate(cleanedName = plantName,
         isConifer = NA,
         notes= NA) %>%
  select(plantName, cleanedName, sciName, itis_id, rank, notes, isConifer)

# 3.3  Append the matched results to officialPlantList and save with date in the filename.

officialPlantList <- rbind(officialPlantList, matched_new_species)

# 3.4  For results that don't match, write to a file and examine manually in Excel (as .csv), and add a new cleanedName if you can figure out what the original name is referring to.
unmatched_new_species <- filter(cleanedNewNames, is.na(cleanedNewNames$itis_id)) %>% 
  rename(plantName = Species) %>%
  mutate(cleanedName = NA,
         isConifer = NA,
         notes= NA) %>%
  select(plantName, cleanedName, sciName, itis_id, rank, notes, isConifer) %>%
write.csv(unmatched_new_species, paste0("ProjectCleaningNames/unmatched_new_species_", Sys.Date(), ".csv"), row.names = F)    

  
# 3.5  Manually go through all rows in unmatched_new_species...csv file and fill in a taxonomically valid cleanedName.
# -- Use resolver.globalnames.org to find potential synonymns or other authorities that recognize the name (put in notes).
# -- After manually fixing all entries, then read in .csv as a dataframe which will have the original plantName and a new cleanedName.
# -- If research does not yield a valid taxonomic name, leave cleanedName as NA.
listOfUnmatchedFiles = list.files('ProjectCleaningNames')[str_detect(list.files('ProjectCleaningNames'), '^unmatched_new_species_with_NA')]
mostRecentUnmatchedFile = listOfUnmatchedFiles[length(listOfUnmatchedFiles)]
  
manually_matched_new_species <- read.csv(paste0('ProjectCleaningNames/', mostRecentUnmatchedFile))

# 3.6  Run the cleanedName column of that dataframe through cleanNamesThruITIS(), rename "Species" as "cleanedName" and join 
#     the results back to the original manually created dataframe that includes both plantName and cleanedName by cleanedName.
cleanedManuallyEnteredNames = cleanNamesThruITIS(manually_matched_new_species$cleanedName) 

manuallyCleanedRecordsWithITIS = left_join(manually_matched_new_species[, c('plantName', 'cleanedName')], 
                                           cleanedManuallyEnteredNames, by = c('cleanedName' = 'Species')) %>%
  mutate(isConifer = NA,
         notes = NA) %>%
  select(plantName, cleanedName, sciName, itis_id, rank, notes, isConifer)


# 4. Manually exmamine names that still don't match
# 5. Add cleaning code to fix typos such that the cleanedPlantName will match in ITIS (the long list of plantName...?)

# 6. Add all new plantName names to officialPlantList (after possibly running through the previous steps a couple of times)("paste" to the bottom of officialPlantList)
officialPlantList = rbind(officialPlantList, cleanedManuallyEnteredNames)

# 7. Writing an updated officialPlantList
write.csv(officialPlantList, paste("ProjectCleaningNames/officialPlantList", Sys.Date(), ".csv", sep = ""), row.names = F)
}