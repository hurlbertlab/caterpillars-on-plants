library(dplyr)
library(taxize)
library(stringr)

officialPlantList = read.csv('ProjectCleaningNames/cleanedPlantList.csv')
#TEMPORARY
officialPlantList = officialPlantList %>%
  mutate(plantName = cleanedName) %>%
  select(plantName, cleanedName, sciName, itis_id, rank, isConifer, notes)

sites = read.csv(list.files()[grep('Site.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)
plantList_rerun = read.csv('ProjectCleaningNames/plantList_rerun.csv')


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
  
  return(plantList)
  
}





#Obtaining masterList as of Fall 2022
new_species_create_list <- anti_join(plantList_rerun, officialPlantList, by = "cleanedName")  
write.csv(new_species_create_list, paste("ProjectCleaningNames/newSpecies_to_create_list", Sys.Date(), ".csv", sep = ""), row.names = F)    
  ##where are the plantNames in plantList_rerun? is it cleanName like line 7?
  ###run these thru ITIS until everything and then manually until everything has an ID manually (won't run)
  ####once this is added to officialPlantList then that's the starting official list for the new workflow

## NEW WORKFLOW ##
# 1. read in latest Plants.csv
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)

# 2. Find new names not in plantName of officialPLantList
new_species <- anti_join(officialPlantList, plants, by = c("plantName" = "Species"))
write.csv(new_species, paste("ProjectCleaningNames/newSpecies", Sys.Date(), ".csv", sep = ""), row.names = F)

    #OR#
new_species <- plants %>% 
  rename(plantName = Species) %>%
  distinct(plantName) %>%
  # select rerun sciName entries that are NOT (!) in sciName from cleaned list
  filter(!plantName %in% officialPlantList$plantName) 

write.csv(new_species, paste("ProjectCleaningNames/newSpecies", Sys.Date(), ".csv", sep = ""), row.names = F)
##this is just the plantName column?


# 3. Run new entries through ITIS / Match new names using taxize


cleanedNewNames = cleanNamesThruITIS(new_species$plantName)


# 3.1  Separate out results that did vs did not match in ITIS


# 3.2  For results that matched, rename "Species" as "plantName", add notes, isConifer, and cleanedName = plantName


# 3.3  Append the matched results to officialPlantList and save with date in the filename.


# 3.4  For results that don't match, write to a file and examine manually in Excel (as .csv), and add a new cleanedName if you can figure out what the original name is referring to.

# 3.5  Then read in .csv as a dataframe which will have the original plantName and a new cleanedName

# 3.6  Run the cleanedName column of that dataframe through cleanNamesThruITIS(), rename "Species" as "cleanedName" and join the results back to the original manually created dataframe that includes both plantName and cleanedName by cleanedName.






# 4. Manually exmamine names that still don't match
# 5. Add cleaning code to fix typos such that the cleanedPlantName will match in ITIS (the long list of plantName...?)

# 6. Add all new plantName names to officialPlantList (after possibly running through the previous steps a couple of times)("paste" to the bottom of officialPlantList)
officialPlantList <- full_join(new_species, officialPlantList) 
          #OR pasting new species manually to the bottom of the official list

# 7. Writing an updated officialPlantList
write.csv(officialPlantList, paste("ProjectCleaningNames/officialPlantList", Sys.Date(), ".csv", sep = ""), row.names = F)
