library(dplyr)
library(taxize)
library(stringr)

officialPlantList = read.csv('ProjectCleaningNames/cleanedPlantList.csv')


sites = read.csv(list.files()[grep('Site.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)
plants = read.csv(list.files()[grep('Plant.csv', list.files())], header = TRUE, stringsAsFactors = FALSE)

coniferList = unique(plants[, c('Species', 'IsConifer')])

plantspp = data.frame(plantName = unique(plants$Species)) %>%
  mutate(plantName = as.character(plantName),
         cleanedPlantName = case_when(
           plantName == "American witch-hazel" ~ "witch-hazel",
           plantName == "American yellowwood" ~ "Kentucky yellowwood",
           plantName == "Arrowwood" ~ "Viburnum",
           plantName == "Arrowwood viburnum" ~ "Viburnum",
           plantName == "Blue beech" ~ "Carpinus caroliniana",
           plantName == "Big-leaf dogwood" ~ "Cornus macrophylla",
           plantName == "Black haw" ~ "Blackhaw",
           plantName == "Black birch" ~ "Sweet birch",
           plantName == "Blueberry vaccinium sp." ~ "Vaccinium",
           plantName == "Box elder (acer negundo)" ~ "Acer negundo",
           plantName == "Boxelder maple" ~ "Acer negundo",
           plantName == "Broadleaf hawthorn" ~ "Kansas hawthorn",
           plantName == "Burr oak" ~ "Bur oak",
           plantName == "Bush honeysuckle (caprifoliaceae family)" ~ "Bush honeysuckle",
           plantName == "Buttonwood-mangrove" ~ "Button mangrove",
           plantName == "California-laurel" ~ "California laurel",
           plantName == "Carolina lauralcherry" ~ "Carolina laurelcherry",
           plantName == "Carolina cherry laurel" ~ "Carolina laurelcherry",
           plantName == "Carolina laurel cherry" ~ "Carolina laurelcherry",
           plantName == "Cherry and plum" ~ "Prunus",
           plantName == "Chinese tallowtree" ~ "Chinese tallow tree",
           plantName == "Coastal azalea" ~ "Dwarf azalea",
           plantName == "Common spicebush" ~ "Northern spicebush",
           plantName == "Common privet" ~ "Wild privet",
           plantName == "Crab apple sp" ~ "Malus",
           plantName == "Crepe myrtle" ~ "Crapemyrtle",
           plantName == "Cornelian-cherry dogwood" ~ "Cornelian cherry",
           plantName == "Cornus stolonifera" ~ "Cornus sericea ssp. Sericea",
           plantName == "Cucumber magnolia" ~ "Cucumbertree",
           plantName == "Eastern sweet shrub" ~ "Eastern sweetshrub",
           plantName == "Fremont cottonwood" ~ "Fremont's cottonwood",
           plantName == "Fishpole bamboo" ~ "Golden bamboo",
           plantName == "Florida forestiera" ~ "Florida privet",
           plantName == "Fraser magnolia" ~ "Fraser's magnolia",
           plantName == "Great rhododendron" ~ "Great laurel",
           plantName == "Grey dogwood" ~ "Gray dogwood",
           plantName == "Horse chestnut" ~ "Horsechestnut",
           plantName == "Japanese spindle tree" ~ "Japanese spindletree",
           plantName == "Japanese cherry" ~ "Japanese flowering cherry",
           plantName == "Kwanzan cherry" ~ "Japanese flowering cherry",
           plantName == "Leatherleaf mahonia" ~ "Beale's barberry",
           plantName == "Lombardy poplar" ~ "Lombardy's poplar",
           plantName == "Maple leaf viburnum" ~ "Mapleleaf viburnum",
           plantName == "Malus sargentii" ~ "Malus sieboldii",
           plantName == "Malus siberica" ~ "Malus baccata",
           plantName == "Mountain or fraser magnolia" ~ "Fraser's magnolia",
           plantName == "Myrica cerifera" ~ "Morella cerifera",
           plantName == "Northern white-cedar" ~ "Northern white cedar",
           plantName == "Nothoscordum gracile (onionweed)" ~ "Nothoscordum gracile",
           plantName == "Heartleaf birch" ~ "Mountain white birch",
           plantName == "Halesia tetraptera" ~ "Halesia carolina",
           plantName == "Honey suckle" ~ "Honeysuckle",
           plantName == "Hop hornbeam" ~ "Hophornbeam",
           plantName == "Hop-hornbeam" ~ "Hophornbeam",
           plantName == "Hydrangeas" ~ "Hydrangea",
           plantName == "Magnolia figo" ~ "Chinese tulip tree",
           plantName == "Olive quihoui??" ~ "NA",
           plantName == "Ozark witch hazel" ~ "Ozark witchhazel",
           plantName == "Pear tree" ~ "Pear",
           plantName == "Red bay" ~ "Redbay",
           plantName == "Reeves spiraea" ~ "Bridalwreath spirea",
           plantName == "Prunus spp.(apricot)" ~ "Prunus",
           plantName == "Prunus spp.(cherry)" ~ "Prunus",
           plantName == "Prunus spp.(plum)" ~ "Prunus",
           plantName == "Prunus serulata (japanese cherry blossom)" ~ "Prunus serultaa",
           plantName == "Quihoui???" ~ "NA",
           plantName == "Red osier dogwood" ~ "Redosier dogwood",
           plantName == "Red-osier dogwood" ~ "Redosier dogwood",
           plantName == "Red osier dogwood (cornus sericea)" ~ "Redosier dogwood",
           plantName == "Sweet gum" ~ "Sweetgum",
           plantName == "Silver poplar" ~ "White poplar",
           plantName == "Simplocos tinctoria" ~ "Symplocos tinctoria",
           plantName == "Shumard oak" ~ "Shumard's oak",
           plantName == "Speckled japanese aralla" ~ "Paperplant",
           plantName == "Symac?" ~ "NA",
           plantName == "Symplocos tinctotria" ~ "Symplocos tinctoria",
           plantName == "Trembling aspen" ~ "Quaking aspen",
           plantName == "Vaccinum corymbosum" ~ "Vaccinium corymbosum",
           plantName == "Vaccineum stamineum" ~ "Vaccinium stamineum",
           plantName == "Virburnum" ~ "Viburnum",
           plantName == "Vibirnum dentatum" ~ "Viburnum dentatum",
           plantName == "Viburnum cassinoide" ~ "Viburnum cassinoides",
           plantName == "Viburnum trilobum" ~ "	Viburnum opulus var. americanum",
           plantName == "Vaccinium caesariense" ~ "Vaccinium corymbosum",
           plantName == "Wafer ash" ~ "Hoptree",
           plantName == "White elm" ~ "American elm",
           plantName == "Wild cherry" ~ "Prunus",
           plantName == "Winged eonymous" ~ "Winged Euonymus",
           plantName == "Weeping cherry" ~ "Winter-flowering cherry",
           plantName == "Weeping knoss dogwood" ~ "Kousa dogwood",
           plantName == "Yellow buckeye" ~ "Aesculus flava",
           plantName %in% c("", "8", "9", "N/A", "Unknown", "Na") ~ "NA",
           TRUE ~ plantName),
         cleanedPlantName = str_replace(cleanedPlantName, " spp.$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " spp$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " sp$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " sp.$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, " species$", ""),
         cleanedPlantName = str_replace(cleanedPlantName, "\\?", "")
  ) %>%
  left_join(coniferList, by = c('plantName' = 'Species')) %>%
  full_join(officialPlantList, by = c('cleanedPlantName' = 'cleanedName'))

# Run new entries through ITIS
# Join to official master plant list
# Manually correct names, etc
# Re-run creation of plantList


# Example of writing file with date in the filename
write.csv(plantList, paste("ProjectCleaningNames/cleanedPlantList", Sys.Date(), ".csv", sep = ""), row.names = F)