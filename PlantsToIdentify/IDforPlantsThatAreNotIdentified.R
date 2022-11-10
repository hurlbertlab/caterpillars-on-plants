library(dplyr)
library(xml2)
library(stringr)
library(rvest)

# Loading data_repo
data_repo <- "https://github.com/hurlbertlab/caterpillars-count-data"
webpage <- read_html(data_repo)
repo_links <- html_attr(html_nodes(webpage, "a"), "href")
data_links <- tibble(link = repo_links[grepl(".csv", repo_links)]) %>%
  mutate(file_name = word(link, 6, 6, sep = "/"))

# Read data files from data repo links
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-count-data/master/"

sites = read.csv(paste(github_raw, filter(data_links, grepl("Site.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

surveys = read.csv(paste(github_raw, filter(data_links, grepl("Survey.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

plants = read.csv(paste(github_raw, filter(data_links, grepl("Plant.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

ArthropodSighting = read.csv(paste(github_raw, filter(data_links, grepl("ArthropodSighting.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

fullDataset <- read.csv(list.files('data', full.names = T)[str_detect(list.files('data'), '^fullDataset')])

officialPlantList <- read.csv(list.files('ProjectCleaningNames', full.names = T)[str_detect(list.files('ProjectCleaningNames'), '^officialPlantList_')])


# Filtering plants to where Species is N/A and obtaining that ID
unidentifiedBranches <- filter(plants, plants$Species == "N/A") %>%
                        select(ID, Species)

# Filtering surveys to find where PlantSpecies has a name entered by a user
# and comparing to where Species has not been identified
# (one row per PlantFK)
userIdentifiedBranches <- surveys %>%
  filter(!PlantSpecies %in% c("N/A","","Hello","Dvt","Dvz","Dvt","N/a","Tree","Unknown","Unknown, will take picture")) %>%
  select(UserFKOfObserver, PlantSpecies, PlantFK) %>%
  left_join(unidentifiedBranches, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  group_by(PlantFK) %>%
  summarize(PlantSpecies = paste(PlantSpecies, collapse = ", ")) %>%
  mutate(InferredName = NA, 
         NameConfidence = NA)

write.csv(userIdentifiedBranches, "PlantsToIdentify/userIdentifiedBranches.csv")

# In Excel, fill in the inferred name if there's agreement and confidence rating
# Giving a confidence rating for the most agreed upon name given by users 
# 1 is the least confident meaning disagreement, 2 means only one name ever entered,
# 3 is the most confident with all entries agreeing multiple times)
# Save a new version of the file with InferredName and NameConfidence (and Notes) values where appropriate

#ratedUserIdentifiedBranches <- read.csv(file = 'PlantsToIdentify/userIdentifiedBranches.csv') 

# Trying to ID photos based off of user suggestions (confidence rating) and photo ID
fullUserIdentifiedBranches <- read.csv(file = 'PlantsToIdentify/userIdentifiedBranches.csv') %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  left_join(surveys, by = c('PlantFK', 'PlantSpecies')) %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  select(Name, Region, PlantFK, PlantSpecies, InferredName, NameConfidence, Notes) %>%
  rename('UserSuggestedName' = 'PlantSpecies')
  
# deleted branchesWithPhotos bc the data is already entered and duplicate of some below
allBranchesWithPhotos <- surveys %>%
  filter(ObservationMethod == "Visual") %>% 
  left_join(ArthropodSighting, by = c('ID' = 'SurveyFK')) %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  filter(PhotoURL != "") %>% 
  group_by(PlantFK) %>%
  summarize(PhotoURL = paste(PhotoURL, collapse = ", ")) %>%
  select(PlantFK, PhotoURL)

# In Excel, look at the PhotoURL and paste it after https://caterpillarscount.unc.edu/images/arthropods/
# Joining userIdentifiedBranches.csv and allBranchesWithPhotos.csv to compile a doc with
# the names suggested by users and plants identified by photos and add to fullUserIdentifiedBranches.csv
JoinedDoc <- full_join(fullUserIdentifiedBranches, allBranchesWithPhotos, by = c('PlantFK'))
write.csv(JoinedDoc, "PlantsToIdentify/JoinedDoc.csv")

# Joining the previously unidentified photos to the fullDataset
JoinedPhotoAndOccurrenceToFull <- fullDataset %>%
  left_join(fullUserIdentifiedBranches, by = 'PlantFK') %>%
  mutate(InferredName = ifelse(is.na(InferredName), Species, InferredName)) %>%
  left_join(officialPlantList, by = c("InferredName" = "cleanedPlantName")) %>%
  arrange(desc(NameConfidence))
  
write.csv(JoinedPhotoAndOccurrenceToFull, "PlantsToIdentify/JoinedPhotoAndOccurrenceToFull.csv")

# use this to ensure that the number of rows is not changing while joining the differnet files above
#length(unique(JoinedPhotoAndOccurrenceToFull$PlantFK))
