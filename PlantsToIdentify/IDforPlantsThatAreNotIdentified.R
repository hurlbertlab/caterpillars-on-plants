# This script is the workflow for identifying individual survey branches that have not been id'ed by site managers where Plant Species name is blank (NA). In those cases, the branch is potentially identifiable to family, genus, or species either because
#  1) users entered a plant species name while conducting a survey, or
#  2) arthropod photos from that survey branch might allow the inference of plant species

library(dplyr)
library(xml2)
library(stringr)
library(rvest)

# Existing table of inferred names (swap out for most recent file)
inferredNames = read.csv("PlantsToIdentify/inferredPlantNames_2024-03-13.csv")

# Specifying data_repo
data_repo <- "https://github.com/hurlbertlab/caterpillars-count-data"
webpage <- read_html(data_repo)
repo_links <- unique(html_attr(html_nodes(webpage, "a"), "href"))
data_links <- tibble(link = repo_links[grepl(".csv", repo_links)]) %>%
  mutate(file_name = word(link, 6, 6, sep = "/"))

# Read data files from data repo links.
# This is necessary because PhotoURL and user-specified plant names are in 'surveys' and 'ArthropodSightings' tables
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-count-data/master/"

sites = read.csv(paste(github_raw, filter(data_links, grepl("Site.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

surveys = read.csv(paste(github_raw, filter(data_links, grepl("Survey.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

plants = read.csv(paste(github_raw, filter(data_links, grepl("Plant.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

ArthropodSighting = read.csv(paste(github_raw, filter(data_links, grepl("ArthropodSighting.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)


# Filtering plants to where Species is N/A. These include species which we may have evaluated previously, but it's possible that new photos or new user-entered names have been added since, so we take all of it.
unidentifiedBranches <- plants %>% 
  filter(Species == "N/A") 


# Filtering surveys to find where PlantSpecies has a name entered by a user at least once
userIdentifiedBranches <- surveys %>%
  filter(!PlantSpecies %in% c("N/A","","Hello","Dvt","Dvz","Dvt","N/a","Tree","Unknown","Unknown, will take picture")) %>%
  select(UserFKOfObserver, PlantSpecies, PlantFK) %>%
  left_join(unidentifiedBranches, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  group_by(PlantFK) %>%
  summarize(PlantSpecies = paste(PlantSpecies, collapse = ", ")) %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  left_join(surveys, by = c('PlantFK', 'PlantSpecies')) %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  filter(Name != "Example Site") %>%
  select(Name, Region, PlantFK, PlantSpecies) %>%
  rename('UserSuggestedName' = 'PlantSpecies')

# All survey branches without a Species name where an arthropod photo has been taken
allBranchesWithPhotos <- surveys %>%
  filter(ObservationMethod == "Visual") %>% 
  left_join(ArthropodSighting, by = c('ID' = 'SurveyFK')) %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  filter(PhotoURL != "",
         Name != "Example Site") %>% 
  group_by(PlantFK, Name, Region) %>%
  summarize(PhotoURL = paste(PhotoURL, collapse = ", ")) %>%
  select(Name, Region, PlantFK, PhotoURL)

# Join to get one dataframe with both user-entered names as well as photos
plantsToIdentify = full_join(userIdentifiedBranches, allBranchesWithPhotos, by = c('Name', 'Region', 'PlantFK')) %>% 
  arrange(Name, PlantFK)

# New branches with either user-entered names or photos that have not been examined before
newPlantsToIdentify = plantsToIdentify %>%
  filter(!PlantFK %in% inferredNames$PlantFK) %>%
  mutate(InferredName = NA,
         NameConfidence = NA,
         Notes = NA,
         New = 'Y')

# Check branches that have been previously examined for which NameConfidence < 3 to see whether there are new user-entered names or photos by comparing the number of characters in the UserSuggestedNames and PhotoURL fields (if new names or photos have been added, the number will be larger)
oldPlantsToIdentify = plantsToIdentify %>%
  filter(PlantFK %in% inferredNames$PlantFK) %>%
  mutate(ncharNamesNew = nchar(UserSuggestedName, keepNA = F),
         ncharPhotoNew = nchar(PhotoURL, keepNA = F)) %>%
  left_join(inferredNames, by = c('Name', 'Region', 'PlantFK')) %>%
  mutate(ncharNamesOld = nchar(UserSuggestedName.y, keepNA = F),
         ncharPhotoOld = nchar(PhotoURL.y, keepNA = F),
         UserSuggestedName = ifelse(ncharNamesNew > ncharNamesOld, UserSuggestedName.x, UserSuggestedName.y),
         PhotoURL = ifelse(ncharPhotoNew > ncharPhotoOld, PhotoURL.x, PhotoURL.y),
         New = ifelse((ncharNamesNew > ncharNamesOld | ncharPhotoNew > ncharPhotoOld) & NameConfidence < 3, 'Y', 'N')) %>%
  select(Name, Region, PlantFK, UserSuggestedName, PhotoURL, InferredName, NameConfidence, Notes, New)

newInferredNames = rbind(oldPlantsToIdentify, newPlantsToIdentify)

# Examine each record where New == 'Y' manually (e.g. in Excel), fill in the inferred name if there's agreement, and assign a confidence rating based on user agreement.
# 1 is the least confident meaning there is disagreement among user-entered names, 
# 2 could mean only one name ever entered, or that it is identifiable to genus but not species from photos
# 3 is the most confident with all entries agreeing multiple times, or photos support id.

write.csv(newInferredNames, paste("PlantsToIdentify/inferredPlantNames_", Sys.Date(), ".csv", sep = ""), row.names = F)
