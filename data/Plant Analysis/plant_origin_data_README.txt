### README for raw data files for plant origin data.

## PLANTS_Data_Expanded_2024.csv
* provided by Gerry Moore from USDA by email on 12 March 2024 (because website downloads currently broken)
* Did a find-and-replace for all "×" characters to be changed to "X"
--because CC! Official Plant List uses "X", and because "×" gets encoded improperly when writing from R
* Per the USDA PLANTS Database website, species were listed as "native" if their L48 native status was listed as N, N?, NI, or NI? and "exotic" if their L48 native status was listed as GP, GP?, I, I?, N?I, W, or W?.
* Using this as a standardized source of plant origin information, there is no need for the datasets below.


## tallamy_shropshire_2009_plant_genera.csv
* Data associated with Tallamy and Shropshire 2009
* provided by Desiree Narango or Doug Tallamy
* At the genus level, specifies plant origin as well as total number of Lepidoptera species known to use the host genus


## gbif_plant_origin_verbatim.txt
* Citation for "Global Register of Introduced and Invasive Species - United States (Contiguous) (ver.2.0, 2022)"
* GBIF.org (12 March 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.mxhtx8
* Includes 8,527 introduced and invasive species of all kinds (including ~4,000 plants)


## ppp3_10305-sup-0001-dataset s1.xlsx
* Data from Carrero et al. 2022. A standardized checklist of United States native tree species and threat assessments to prioritize and coordinate action.
* Includes a checklist of 882 native tree species of the United States
