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

# Filtering plants to where Species is N/A and obtaining that ID
unidentifiedBranches <- filter(plants, plants$Species == "N/A") %>%
                        select(ID, Species)

# Filtering surveys to find where PlantSpecies has a name entered by a user
# and comparing to where Species has not been identified
userIdentifiedBranches <- surveys %>%
  filter(PlantSpecies != "N/A", PlantSpecies != "") %>%
  select(UserFKOfObserver, PlantSpecies, PlantFK) %>%
  left_join(unidentifiedBranches, by = c('PlantFK' = 'ID')) %>%
  filter(Species == "N/A")
