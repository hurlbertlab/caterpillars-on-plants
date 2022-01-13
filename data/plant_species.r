# This script automatically matches plant species names provided by users


library(dplyr)
library(taxize)

sites = read.csv(list.files()[grep('Site.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)
plantspp = data.frame(plantName = unique(plants$Species)) %>%
  mutate(plantName = as.character(plantName),
         cleanedPlantName = case_when(
           plantName == "American witch-hazel" ~ "witch-hazel",
           plantName == "American yellowwood" ~ "Kentucky yellowwood",
           plantName == "Arrowwood" ~ "Viburnum",
           plantName == "Arrowwood viburnum" ~ "Viburnum",
           plantName == "Big-leaf dogwood" ~ "Cornus macrophylla",
           plantName == "Blueberry vaccinium sp." ~ "Vaccinium",
           plantName == "Box elder (acer negundo)" ~ "Acer negundo",
           plantName == "Boxelder maple" ~ "Acer negundo",
           plantName == "Burr oak" ~ "Bur oak",
           plantName == "Bush honeysuckle (caprifoliaceae family)" ~ "Bush honeysuckle",
           plantName == "California-laurel" ~ "California laurel",
           plantName == "Crab apple sp" ~ "Malus",
           plantName == "Crepe myrtle" ~ "Crapemyrtle",
           plantName == "Eastern sweet shrub" ~ "Eastern sweetshrub",
           plantName == "Fraser magnolia" ~ "Fraser's magnolia",
           plantName == "Mountain or fraser magnolia" ~ "Fraser's magnolia",
           plantName == "Honey suckle" ~ "Honeysuckle",
           plantName == "Hop hornbeam" ~ "Hophornbeam",
           plantName == "Hop-hornbeam" ~ "Hophornbeam",
           plantName == "Hydrangeas" ~ "Hydrangea",
           plantName == "Pear tree" ~ "Pear",
           plantName == "Red bay" ~ "Redbay",
           plantName == "Red osier dogwood" ~ "Redosier dogwood",
           plantName == "Red-osier dogwood" ~ "Redosier dogwood",
           plantName == "Red osier dogwood (cornus sericea)" ~ "Redosier dogwood",
           plantName == "Sweet gum" ~ "Sweetgum",
           plantName == "Simplocos tinctoria" ~ "Symplocos tinctoria",
           plantName == "Symplocos tinctotria" ~ "Symplocos tinctoria",
           plantName == "Trembling aspen" ~ "Quaking aspen",
           plantName == "Vaccinum corymbosum" ~ "Vaccinium corymbosum",
           plantName == "Virburnum" ~ "Viburnum",
           plantName == "Wild cherry" ~ "Prunus",
           plantName %in% c("", "8", "9", "N/A", "Unknown") ~ "NA",
           TRUE ~ plantName),
         cleanedPlantName = str_replace(cleanedPlantName, " spp.$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " spp$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " sp$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " sp.$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " species$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, "?", "")
         ) 


# Finding matching taxonomic entity in ITIS
plantList = data.frame(cleanedName = unique(plantspp$cleanedPlantName[plantspp$cleanedPlantName != "NA"]), 
                       sciName = NA, 
                       itis_id = NA,    
                       rank = NA)

for (i in 1:nrow(plantList)) {
  
  print(paste(i, "of", nrow(plantList), "/n"))
  
  hierarchy = classification(plantList$cleanedName[i], db = 'itis', accepted = TRUE)[[1]]
  
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





plants2 = left_join(plants, sites[, c('ID', 'Region')], by = c('SiteFK' = 'ID')) %>%
  distinct(Species, Region) %>% arrange(Species, Region)

multiStateSpp = plants2 %>%
  count(Species) %>%
  filter(n > 1)

states = c()
for (s in multiStateSpp$Species) {
  states = c(states, paste(plants2$Region[plants2$Species == s], collapse = ', '))
}

# Get list of regions each species is found in
plantSppRegions = plants2 %>%
  filter(! Species %in% multiStateSpp$Species) %>%
  rbind(data.frame(Species = multiStateSpp$Species, Region = states)) %>%
  arrange(Species)







plantspp2 = full_join(plantspp, plantSppRegions, by = c('plantName' = 'Species'))

write.table(plantspp2, 'z:/projects/caterpillarscount/fia by state/cc_plantlist_20190226.txt', sep = '\t', row.names = F)


allplants = left_join(fia, plantspp, by = 'lower') %>%
  select(commonName, lower, sciName) %>%
  rbind(cc_only)

write.table(allplants, 'z:/projects/caterpillarscount/fia by state/fia_plus_cc_plantlist.txt', sep = '\t', row.names = F)


#fia = read.table('z:/projects/caterpillarscount/FIA by state/tree_ids.txt', header= T, sep = '\t')
#fia$lower = tolower(fia$commonName)

fiacc = read.csv('plantSpecies/FIA_CC_plantList.csv', header = T, quote = '\"')
uniquePlantList = fiacc %>%
  distinct(commonName, sciName)
