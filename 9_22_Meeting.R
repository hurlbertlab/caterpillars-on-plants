library(dplyr)
library(xml2)
library(stringr)

data_repo <- "https://github.com/hurlbertlab/caterpillars-count-data"
webpage <- read_html(data_repo)
repo_links <- html_attr(html_nodes(webpage, "a"), "href")
data_links <- tibble(link = repo_links[grepl(".csv", repo_links)]) %>%
  mutate(file_name = word(link, 6, 6, sep = "/"))

## Read data files from data repo links
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-count-data/master/"

sites = read.csv(paste(github_raw, filter(data_links, grepl("Site.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

surveys = read.csv(paste(github_raw, filter(data_links, grepl("Survey.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

plants = read.csv(paste(github_raw, filter(data_links, grepl("Plant.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

# Filtering plants to find NA in species to get their ID
unidentifiedSpecies <- filter(plants, plants$Species == "N/A")

# Filtering surveys to find where plantFK (which is equal to ID in plants) is in the plants$ID
plantFKinplantID <- surveys %>%
  # where PlantFK is in/equal to the ID column 
  filter(PlantFK %in% unidentifiedSpecies$ID) %>%
  #plantFKinplantID[order(plantFKinplantID$PlantSpecies),]

ordered <- plantFKinplantID[order(plantFKinplantID$PlantSpecies),]
  