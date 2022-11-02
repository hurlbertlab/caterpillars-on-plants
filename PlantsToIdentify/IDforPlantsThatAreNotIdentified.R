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

# Filtering plants to where Species is N/A and obtaining that ID
unidentifiedBranches <- filter(plants, plants$Species == "N/A") %>%
                        select(ID, Species)

# Filtering surveys to find where PlantSpecies has a name entered by a user
# and comparing to where Species has not been identified
userIdentifiedBranches <- surveys %>%
  filter(!PlantSpecies %in% c("N/A","","Hello","Dvt","Dvz","Dvt","N/a","Tree","Unknown","Unknown, will take picture")) %>%
  select(UserFKOfObserver, PlantSpecies, PlantFK) %>%
  left_join(unidentifiedBranches, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  group_by(PlantFK) %>%
  summarize(PlantSpecies = paste(PlantSpecies, collapse = ", ")) %>%
  # Giving a confidence rating for the most agreed upon name given by users 
  # 1 is the least confident meaning disagreement, 2 means only one name ever entered,
  # 3 is the most confident with all entries agreeing multiple times)
  mutate(InferredName = NA, 
         ConfidenceInterval = NA, 
         Notes = NA)

write.csv(userIdentifiedBranches, "PlantsToIdentify/userIdentifiedBranches.csv", row.names = F)
#### where did UserofFKofObserver go? if it's added it creates multiple rows with the same UserID

# In Excel, fill in the inferred name if there's agreement and confidence rating
# Obtain the branches with photos as they have an iNat ID and try to ID the plant from the photo
ratedUserIdentifiedBranches <- read.csv(file = 'PlantsToIdentify/userIdentifiedBranches.csv') %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  left_join(surveys, by = c('PlantFK', 'Notes', 'PlantSpecies')) %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  select(Name, Region, PlantFK, PlantSpecies, InferredName, ConfidenceInterval, Notes) %>%
  
branchesWithPhotos <- surveys %>%
  filter(PlantFK %in% ratedUserIdentifiedBranches$PlantFK,
         ObservationMethod == "Visual") %>% 
  left_join(ArthropodSighting, by = c('ID' = 'SurveyFK')) %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  filter(PhotoURL != "") %>% 
  select(Name, Region, PlantFK, PlantSpecies, Notes.x, PhotoURL)

write.csv(branchesWithPhotos, "PlantsToIdentify/branchesWithPhotos.csv", row.names = F)
read.csv('PlantsToIdentify/branchesWithPhotos.csv')
#added in Excel, InferredName and myNotes column to try and use the iNaturalist photos to determine plant species

#think through logic of all the things we want to look at
#don't need two sets of code just the one with stuff you think would be useful
allbranchesWithPhotos <- surveys %>%
  filter(ObservationMethod == "Visual") %>% 
  left_join(ArthropodSighting, by = c('ID' = 'SurveyFK')) %>%
  left_join(plants, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A") %>%
  left_join(sites, by = c('SiteFK' = 'ID')) %>%
  filter(PhotoURL != "") %>% 
  select(Name, Region, PlantFK, PlantSpecies, Notes.x, PhotoURL)

JoinedBranches <- anti_join(allbranchesWithPhotos, branchesWithPhotos, by = c('PlantFK', 'Name', 'Region', 'PlantSpecies', 'Notes.x', 'PhotoURL'))
write.csv(JoinedBranches, "PlantsToIdentify/JoinedBranches.csv", row.names = F)

# In Excel, look at the PhotoURL and paste it after https://caterpillarscount.unc.edu/images/arthropods/
# Try to ID plants based off of photos and add to InferredName in ratedUserIdentifiedBranches.csv
# Joining the two datasets back together (for easier flow should just use allBranches...)
IDBranchesWithPhotos <- read.csv(file = 'PlantsToIdentify/branchesWithPhotos.csv')
IDJoinedBranches <-read.csv(file = 'PlantsToIdentify/JoinedBranches.csv')

PhotosWithAllBranches <- rbind(IDJoinedBranches, IDBranchesWithPhotos)

#Joining PhotosWithAllBranches (photos) to ratedUserIdentifiedBranches (consistency)

ratedUserIdentifiedBranches %>%
  filter(!PlantFK %in% PhotosWithAllBranches$PlantFK)
#probably should've user rbind for code below but have to get the number of cols to match
UserRatedAndPhotoExaminedBranches <- bind_rows(PhotosWithAllBranches, ratedUserIdentifiedBranches)

#number of values definitely changes so that's not good, fix that
JoinedPhotoAndOccurrenceToFull <- fullDataset %>%
  left_join(UserRatedAndPhotoExaminedBranches, by = 'PlantFK') %>%
  mutate(InferredName = ifelse(is.na(InferredName), Species, InferredName)) %>%
  left_join(officialPlantList, by = c("InferredName" = "cleanedPlantName"))
#use this to ensure that the number of rows is not changing while joining the differnet files above
length(unique(JoinedPhotoAndOccurrenceToFull$InferredName))
