## KEY ##
# userPlantName = the user inputted name for a plant (might be the scientific name, different common names, a genus, etc.)
# cleanedPlantName = a version of the plant name that contains only the scientific name, genus or common name without "spp.", "?", etc.
# sciName = the official ITIS (usually) recognized scientific name, the number of unique species
# ITIS_ID = the ID number that matches to a specific plant taxonomy in the Integrated Taxonomic Information System

# NOTE: Still need a way to efficiently deal with ambiguous common names that could refer to different taxonomic concepts
#       depending on where you are, e.g. "scrub oak", "ironwood", etc.


library(dplyr)
library(taxize)
library(stringr)


# This function takes a vector of species names and checks each one with ITIS,
# returning a dataframe with the name, sciName, itis_id, and rank.
cleanNamesThruITIS = function(speciesList) {
  
  plantList = data.frame(Species = speciesList,
                         sciName = NA, 
                         itis_id = NA,    
                         rank = NA)
  
  for (i in 1:nrow(plantList)) {
    
    print(paste(i, "of", nrow(plantList)))
    
    if (!is.na(speciesList[i]) & nchar(speciesList[i]) >= 3) {  # for names that are at least 3 characters and not NA, try to match
      
      hierarchy = classification(speciesList[i], db = 'itis', accepted = TRUE)[[1]]
      
      # class is logical if taxonomic name does not match any existing names
      if (!is.null(nrow(hierarchy))) {
        plantList$sciName[i] = hierarchy$name[nrow(hierarchy)]
        plantList$itis_id[i] = hierarchy$id[nrow(hierarchy)]
        plantList$rank[i] = hierarchy$rank[nrow(hierarchy)]
        
      } # end if there's a match
    } # end if name should be searched
  } # end for loop
  return(plantList)
}

## NEW WORKFLOW ##
# 1. read in latest Plants.csv
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)

officialPlantListFiles = list.files('ProjectCleaningNames')[str_detect(list.files('ProjectCleaningNames'), 'officialPlantList')]
mostRecentOfficialPlantList = officialPlantListFiles[length(officialPlantListFiles)]
officialPlantList = read.csv(paste0('ProjectCleaningNames/', mostRecentOfficialPlantList), header = T)


# 2. Find new names not in userPlantName of officialPLantList
  # Below isn't a true representation of new species bc the names are not "clean" they still have spp., etc.
new_species <- plants %>% 
  rename(userPlantName = Species) %>%
  distinct(userPlantName) %>%
  # select rerun sciName entries that are NOT (!) in sciName from cleaned list
  filter(!userPlantName %in% officialPlantList$userPlantName) 

# IF THERE ARE NEW SPECIES (nrow(new_species) > 0), MOVE FORWARD. OTHERWISE STOP.

write.csv(new_species, paste("ProjectCleaningNames/newSpecies_", Sys.Date(), ".csv", sep = ""), row.names = F)


# 3. Run new entries through ITIS / Match new names using taxize
cleanedNewNames = cleanNamesThruITIS(new_species$userPlantName)

# 3.1  Separate out results that did vs did not match in ITIS
# For results that matched, rename "Species" as "userPlantName", add notes, isConifer, and cleanedPlantName = userPlantName

matched_new_species <- filter(cleanedNewNames, !is.na(cleanedNewNames$itis_id)) %>%
  rename(userPlantName = Species) %>%
  mutate(cleanedPlantName = userPlantName,
         isConifer = NA,
         notes= NA) %>%
  select(userPlantName, cleanedPlantName, sciName, itis_id, rank, notes, isConifer)

# 3.3  Append the matched results to officialPlantList and save with date in the filename.

officialPlantList <- rbind(officialPlantList, matched_new_species)

# 3.4  For results that don't match, write to a file and examine manually in Excel (as .csv), and add a new cleanedPlantName if you can figure out what the original name is referring to.
unmatched_new_species <- filter(cleanedNewNames, is.na(cleanedNewNames$itis_id)) %>% 
  rename(userPlantName = Species) %>%
  mutate(cleanedPlantName = NA,
         isConifer = NA,
         notes= NA) %>%
  select(userPlantName, cleanedPlantName, sciName, itis_id, rank, notes, isConifer)

write.csv(unmatched_new_species, paste0("ProjectCleaningNames/unmatched_new_species_", Sys.Date(), ".csv"), row.names = F)    

# 3.5  Manually go through all rows in unmatched_new_species...csv file and fill in a taxonomically valid cleanedPlantName.
# -- Use resolver.globalnames.org to find potential synonymns or other authorities that recognize the name (put in notes).
# -- After manually fixing all entries, then read in .csv as a dataframe which will have the original userPlantName and a new cleanedPlantName.
# -- If research does not yield a valid taxonomic name, leave cleanedPlantName as NA.
listOfUnmatchedFiles = list.files('ProjectCleaningNames')[str_detect(list.files('ProjectCleaningNames'), '^unmatched_new_species')]
mostRecentUnmatchedFile = listOfUnmatchedFiles[length(listOfUnmatchedFiles)]
  
manually_matched_new_species <- read.csv(paste0('ProjectCleaningNames/', mostRecentUnmatchedFile))

# 3.6  Run the cleanedPlantName column of that dataframe through cleanNamesThruITIS(), rename "Species" as "cleanedPlantName" and join the results back to the original manually created dataframe that includes both userPlantName and cleanedPlantName by cleanedPlantName.

manually_matched_names_to_clean = manually_matched_new_species$cleanedPlantName[!is.na(manually_matched_new_species$cleanedPlantName)]
cleanedManuallyEnteredNames = cleanNamesThruITIS(manually_matched_names_to_clean) 

manuallyCleanedRecordsWithITIS = left_join(manually_matched_new_species[, c('userPlantName', 'cleanedPlantName')], 
                                           cleanedManuallyEnteredNames, by = c('cleanedPlantName' = 'Species')) %>%
  mutate(isConifer = NA,
         notes = NA) %>%
  select(userPlantName, cleanedPlantName, sciName, itis_id, rank, notes, isConifer)


# 4. Manually examine names that still don't match and correct where possible. Remaining unmatched names will simply be NA.

# 5. Add all new userPlantName names to officialPlantList (after possibly running through the previous steps a couple of times)("paste" to the bottom of officialPlantList)
officialPlantList = rbind(officialPlantList, manuallyCleanedRecordsWithITIS)

# 6. Writing an updated officialPlantList
write.csv(officialPlantList, paste("ProjectCleaningNames/officialPlantList", Sys.Date(), ".csv", sep = ""), row.names = F)
