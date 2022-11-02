library(xml2)
library(stringr)
library(rvest)
library(tidyverse)
library(plyr)

# Loading data_repo
#data_repo <- "https://github.com/hurlbertlab/caterpillars-on-plants/data"
#webpage <- read_html(data_repo)
#repo_links <- html_attr(html_nodes(webpage, "a"), "href")
#data_links <- tibble(link = repo_links[grepl(".csv", repo_links)]) %>%
  mutate(file_name = word(link, 6, 6, sep = "/"))

#fullDataset = read.csv(paste(github_raw, filter(data_links, grepl("fullDataset.csv", file_name))$file_name, sep = ''), header = TRUE, stringsAsFactors = FALSE)

# Creating a matrix by joining BLANK from fullDataset to officialPlantList to get
# the number of times a plant branch is counted at a site
  
fullDataset = read.csv('data/fullDataset_2021-12-01.csv')

officialPlantListFiles = list.files('ProjectCleaningNames')[str_detect(list.files('ProjectCleaningNames'), 'officialPlantList_')]
mostRecentOfficialPlantList = officialPlantListFiles[length(officialPlantListFiles)]
officialPlantList = read.csv(paste0('ProjectCleaningNames/', mostRecentOfficialPlantList), header = T)

MatrixOfSpeciesFoundAtSite<- left_join(fullDataset, officialPlantList, by = c('Species' = 'userPlantName')) %>%
  distinct(sciName, SiteFK) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = SiteFK, values_from = presence)

SiteFKinOfficialPlantList <- left_join(fullDataset, officialPlantList, by = c('Species' = 'userPlantName')) %>%
  d
istinct(sciName, SiteFK) %>%
  dplyr::count(sciName) %>% 
  arrange(desc(n)) %>%
  left_join(MatrixOfSpeciesFoundAtSite, by = 'sciName')




countDuplicateRows <-  ddply(SiteFKinOfficialPlantList,.(sciName, SiteFK), nrow) %>%
  arrange(sciName)
  #rename(countDuplicateRows, sumOfDuplicateRows = V1)
colnames(countDuplicateRows)[3] <- "sumOfDuplicateRows"


#rnames <- officialPlantList$sciName
#cnames <- unique(SiteFKinOfficialPlantList$SiteFK)
#named_matrix <- matrix(SiteFKinOfficialPlantList, nrow = 300, byrow = TRUE, dimnames=list(rnames, cnames))

#is userPlantName accounting for the number of plant branches so is it ok to match with SiteFK??
#want the number of sciNames AT a SiteFK summed
#Red Maple at SiteFK 40 is repeated blank # of times while Flowering Dogwood of this site is repeated # times 
#want to see how many times RedMaple is reported at different sites
#Acer rubrum is sampled10 times at SiteFK 40, 5 at SiteFK 50  
#want to see how many times the sciName repeats 
#numberTimePlantOccurs <- table(SiteFKinOfficialPlantList$sciName) gives the number of times a sciName appears
  
  