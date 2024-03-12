# Create a master list of native/alien status for 4,877 North American plant species from two original sources.
# Additional species may get added to the list in ad hoc fashion.

library(tidyverse)
library(readxl)

## Expanded USDA plants database
# Per the USDA PLANTS Database website, species were listed as "native" if their L48 native status was listed as N, N?, NI, or NI? and "exotic" if their L48 native status was listed as GP, GP?, I, I?, N?I, W, or W?.

origin = read_csv('data/Plant Analysis/PLANTS_Data_Expanded_2024.csv', locale=locale(encoding="latin1")) %>%
  mutate(scientificName = `Scientific Name`,
         nativeStatus = `Native Status`,
         plantOrigin = case_when(
           grepl("L48(N)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("L48(I,N)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("L48(NI)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("L48(N?)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("L48(NI?)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("NA(N)", nativeStatus, fixed = TRUE) ~ 'native',
           grepl("NA(I)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(I)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(I?)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(GP)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(GP?)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(N?I)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(W)", nativeStatus, fixed = TRUE) ~ 'alien',
           grepl("L48(W?)", nativeStatus, fixed = TRUE) ~ 'alien'
         )) %>%
  select(scientificName, Family, nativeStatus, plantOrigin)

# Manually adding origin information where missing for species in CC database
plantListDir = 'c:/git/caterpillars-analysis-public/data/plants/'
officialPlantListFiles = list.files(plantListDir)[str_detect(list.files(plantListDir), 'officialPlantList')]
mostRecentOfficialPlantList = officialPlantListFiles[length(officialPlantListFiles)]
officialPlantList = read.csv(paste0(plantListDir, 
                                    mostRecentOfficialPlantList), header = T)

# Species in the CC Official Plant List with unknown origin that need to be looked up and manually entered below
unknownorigin = filter(origin, is.na(plantOrigin))

officialPlantList[officialPlantList$sciName %in% unknownorigin$scientificName,]

# Manually specified alien species
origin$plantOrigin[origin$scientificName %in% c("Acer buergerianum",
                                                "Acer diabolicum",
                                                "Cornus macrophylla",
                                                "Fatsia japonica",
                                                "Gardenia",
                                                "Maackia amurensis",
                                                "Parrotia persica",
                                                "Prunus mume",
                                                "Quercus glauca",
                                                "Osmanthus fragrans",
                                                "Gardenia volkensii",
                                                "Camptotheca",
                                                "Viburnum odoratissimum",
                                                "Prunus salicina",
                                                "Galphimia glauca")] = 'alien'

write.csv(origin, 'data/Plant Analysis/usda_plant_origin_status.csv', row.names = F)