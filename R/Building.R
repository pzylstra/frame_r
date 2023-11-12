#' Builds the dataframe site.meta from input tables
#' 
#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10 percent)
#' litter - weight in t/ha of fine dead organic material forming the O horizon
#' diameter - mean diameter of surface litter in m
#' fline - the fireline length in m
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param a The number of the record for the site

siteBuilder <- function(site, Structure, a)
{
  
  # CREATE site.meta
  paramSite <- c('overlapping','overlapping','overlapping','overlapping','overlapping', 'fuelLoad',
                 'meanFuelDiameter','meanFinenessLeaves', 'fireLineLength', 'slope', 'temperature',
                 'deadFuelMoistureProp', 'windSpeed')
  
  
  unitsSite <- c(NA, NA, NA, NA, NA, 't/ha', 'm', 'm', 'm', 'deg', 'degC', NA, 'km/h')
  site.meta <- data.frame(matrix(NA, nrow=13, ncol = 5))
  
  names(site.meta) <- c("stratum", "species", "param", "value", "units")
  site.meta$param <- paramSite
  site.meta$units <- unitsSite
  # ENTER VARIABLES
  site.meta$value[1] <- ifelse(Structure$ns_e[Structure$record==a]=='t',"near surface, elevated, overlapped",
                               ifelse(Structure$ns_e[Structure$record==a]=='f',"near surface, elevated, not overlapped",
                                      "near surface, elevated, automatic"))
  site.meta$value[2] <- ifelse(Structure$ns_m[Structure$record==a]=='t',"near surface, midstorey, overlapped",
                               ifelse(Structure$ns_m[Structure$record==a]=='f',"near surface, midstorey, not overlapped",
                                      "near surface, midstorey, automatic"))
  site.meta$value[3] <- ifelse(Structure$e_m[Structure$record==a]=='t',"elevated, midstorey, overlapped",
                               ifelse(Structure$e_m[Structure$record==a]=='f',"elevated, midstorey, not overlapped",
                                      "elevated, midstorey, automatic"))
  site.meta$value[4] <- ifelse(Structure$e_c[Structure$record==a]=='t',"elevated, canopy, overlapped",
                               ifelse(Structure$e_c[Structure$record==a]=='f',"elevated, canopy, not overlapped",
                                      "elevated, canopy, automatic"))
  site.meta$value[5] <- ifelse(Structure$m_c[Structure$record==a]=='t',"midstorey, canopy, overlapped",
                               ifelse(Structure$m_c[Structure$record==a]=='f',"midstorey, canopy, not overlapped",
                                      "midstorey, canopy, automatic"))
  site.meta$value[6] <- site$litter[site$record==a]
  site.meta$value[7] <- if (is.null(site$diameter[site$record==a])) {
    0.005
  } else {
    site$diameter[site$record==a]
  }
  site.meta$value[8] <- 0.00025
  site.meta$value[9] <- site$fLine[site$record==a]
  site.meta$value[10] <- site$slope[site$record==a]
  site.meta$value[11] <- site$temp[site$record==a]
  site.meta$value[12] <- site$dfmc[site$record==a]
  site.meta$value[13] <- site$wind[site$record==a]
  
  return(site.meta)
}



#' Builds the dataframe strata.meta
#'
#' @param Structure 
#' @param Flora 
#' @param a 
#'
#' @return dataframe
#'

strataBuilder <- function(Structure, Flora, a)
{
  # Collect subsets for record
  st <- Structure[Structure$record==a,]
  fl <- Flora[Flora$record==a,]
  
  # CREATE strata.meta
  strata.meta <- data.frame(matrix(NA, nrow=8, ncol = 5))
  
  names(strata.meta) <- c("stratum", "species", "param", "value", "units")
  
  # Fill levelNames
  strata.meta$param[1] <- "levelName"
  strata.meta$value[1] <- "near surface"
  strata.meta$param[2] <- "plantSeparation"
  strata.meta$value[2] <- st$NS[1]
  strata.meta$units[2] <- "m"
  
  strata.meta$param[3] <- "levelName"
  strata.meta$value[3] <- "elevated"
  strata.meta$param[4] <- "plantSeparation"
  strata.meta$value[4] <- st$El[1]
  strata.meta$units[4] <- "m"
  
  strata.meta$param[5] <- "levelName"
  strata.meta$value[5] <- "midstorey"
  strata.meta$param[6] <- "plantSeparation"
  strata.meta$value[6] <- st$Mid[1]
  strata.meta$units[6] <- "m"
  
  strata.meta$param[7] <- "levelName"
  strata.meta$value[7] <- "canopy"
  strata.meta$param[8] <- "plantSeparation"
  strata.meta$value[8] <- st$Can[1]
  strata.meta$units[8] <- "m"
  
  # Find empty strata and remove them
  deleteVector <- vector()
  for (ro in 2:8) {
    if(is.na(strata.meta$value[ro])){
      deleteVector<- c(deleteVector, ro-1, ro)
    }
  }
  if(length(deleteVector)>0){
    strata.meta <- strata.meta[ -deleteVector, ]}
  
  # Number strata
  rows <- as.numeric(nrow(strata.meta))
  num <- 0.5
  for (ro in 1:rows) {
    numb <- ceiling(num)
    strata.meta$stratum[ro] <- numb
    num = num + 0.5
  }
  
  return(strata.meta)
}



#' Builds the dataframe species.values
#'
#' @param Flora 
#' @param site 
#' @param a 
#'
#' @return dataframe

speciesBuilder <- function(Flora, site, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]
  si <- site[site$record==a,]
  
  # CREATE species.values
  ro <- as.numeric(nrow(fl))
  species.values <- data.frame(matrix(NA, nrow=ro, ncol = 13))
  
  names(species.values) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                             "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")
  
  # Enter values
  species.values$stratum <- as.numeric(fl$stratum)
  species.values$name <- as.character(fl$species)
  species.values$liveLeafMoisture <- fl$moisture
  species.values$hc <- fl$base
  species.values$hp <- fl$top
  species.values$he <- fl$he
  species.values$ht <- fl$ht
  species.values$w <- fl$w
  species.values$clumpDiameter <- pmin((fl$top-fl$base),fl$w)*fl$clump
  species.values$clumpSeparation <- fl$openness*species.values$clumpDiameter
  species.values$composition <- as.numeric(fl$comp)
  species.values$deadLeafMoisture <- si$dfmc
  
  species.values <- species.values[order(species.values$stratum),]
  species.values$species <- as.numeric(c(1:ro))
  
  # Calculate gaps
  for(a in 1:ro) {
    if(is.na(species.values$he[a])){
      species.values$he[a] <- species.values$hc[a]
    }
    if(is.na(species.values$ht[a])){
      species.values$ht[a] <- species.values$hp[a]
    }
    if(is.na(species.values$liveLeafMoisture[a])){
      species.values$liveLeafMoisture[a] <- 1
    }
  }
  
  # Ensure strata are numbered consecutively
  stratCurr <- as.data.frame(unique(species.values$stratum))
  names(stratCurr) <- c("orig")
  stratCurr$cor <- 1:as.numeric(nrow(stratCurr))
  
  #Check values
  for(b in 1:as.numeric(nrow(stratCurr))) {
    species.values$stratum[species.values$stratum == stratCurr$orig[b]] <- stratCurr$cor[b]
  }
  
  return(species.values)
}


#' Formats units for parameter files
#' 
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100 percent ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - Percent composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param a The record number for which to build the table


unitBuilder <- function(Flora, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a,]
  fl <- fl[fl$species != "Litter",] %>%
    arrange(stratum)
  
  # CREATE species.units
  ro <- as.numeric(nrow(fl))
  species.units <- data.frame(matrix(NA, nrow=ro, ncol = 13))
  
  names(species.units) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                            "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")
  
  # Enter values
  species.units$stratum <- fl$stratum
  species.units$species <- c(1:ro)
  species.units$clumpDiameter <- "m"
  species.units$clumpSeparation <- "m"
  species.units$hc <- "m"
  species.units$he <- "m"
  species.units$hp <- "m"
  species.units$ht <- "m"
  species.units$w <- "m"
  
  return(species.units)
}



#' Constructs parameter files from formatted datasets
#' 
#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10 percent)
#' litter - weight in t/ha of fine dead organic material forming the O horizon
#' diameter - mean diameter of surface litter in m
#' fline - the fireline length in m
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100 percent ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - Percent composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param default.species.params Leaf traits database
#' @param a The record number for which to build the table
#' @export
#' @examples
#' record <- 1
#' data(site)
#' data(structure)
#' data(flora)
#' data(traits)
#' base.params <- paramBuilder(site, structure, flora, traits, record)

paramBuilder <- function(site, Structure, Flora, default.species.params, a)
{
  # Construct component tables
  site.meta <- siteBuilder(site, Structure, a)
  strata.meta <- strataBuilder(Structure, Flora, a)
  species.values <- speciesBuilder(Flora, site, a)
  species.units <- unitBuilder(Flora, a)
  
  # Build parameter file
  components <- list(site.meta, strata.meta, species.values, species.units)
  names(components) <- c("site.meta", "strata.meta", "species.values", "species.units")
  param <- ffm_assemble_table(components)
  
  # Find weighted mean of crown widths per stratum
  species.values$weightedW <- species.values$composition * species.values$w
  ww <- species.values%>%
    select(stratum, composition, weightedW)%>%
    group_by(stratum) %>%
    summarize_all(sum) %>%
    mutate(mw = weightedW/composition)
  
  # Ensure separation is greater than weighted mean of crown widths
  for (stNum in 1:max(ww$stratum)) {
    sep <- filter(strata.meta, stratum == stNum)
    param <- ffm_set_stratum_param(param, stNum, "plantSeparation", 
                                   pmax(as.numeric(sep$value[2]),ww$mw[stNum]))
  }
  
    # Change species param leafForm to characters
    default.species.params$leafForm <- as.character(default.species.params$leafForm)
    
    # Fill empty traits
    param <- ffm_complete_params(param,default.species.params)
    
    return(param)
}


#' Builds the dataframe site.meta from input tables
#' Building from MEE inputs requires function siteBuilder
#' 
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - A name for the site
#' species - the name of the species, which will call trait data from 
#' 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - Percent composition of that species in the stratum. 
#'    If absent, all species will be considered equally.
#' base - base height of plant crowns (m)
#' he - he height of plant crowns (m)
#' ht - ht height of plant crowns (m)
#' top - top height of plant crowns (m)
#' w - width of plant crowns (m)
#' Hs - standard deviation of the top height of plant crowns (m)
#' Hr - range of the top height of plant crowns (m)
#' weight - weight in t/ha of fine dead organic material forming 
#'    the surface and suspNS layers
#' diameter - mean diameter of surface and suspNS litter in m
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. 
#'    Acceptable values are TRUE, FALSE, or blank, where the outcome 
#'    will be decided by the relative stratum heights.
#' nsR, eR, mR, cR. Species richness (number of species) in each stratum
#' @param a The number of the record for the site
#' @param fLine The fireline length in m
#' @param slope slope in degrees
#' @param temp Ambient temperature deg. C
#' @param dfmc Moisture content of fine dead fuels in whole numbers (eg 0.1 for 10 percent)
#' @param wind Velocity in km/h

buildSiteMeta <- function(Structure, Flora, a, fLine = 100, slope = 0,
                      temp = 30, dfmc = 0.1, wind = 10)
{
  
  # CREATE site.meta
  paramSite <- c('overlapping','overlapping','overlapping','overlapping','overlapping', 'fuelLoad',
                 'meanFuelDiameter','meanFinenessLeaves', 'fireLineLength', 'slope', 'temperature',
                 'deadFuelMoistureProp', 'windSpeed')
  
  
  unitsSite <- c(NA, NA, NA, NA, NA, 't/ha', 'm', 'm', 'm', 'deg', 'degC', NA, 'km/h')
  site.meta <- data.frame(matrix(NA, nrow=13, ncol = 5))
  
  names(site.meta) <- c("stratum", "species", "param", "value", "units")
  site.meta$param <- paramSite
  site.meta$units <- unitsSite
  # ENTER VARIABLES
  site.meta$value[1] <- ifelse(Structure$ns_e[Structure$record==a]==TRUE,"near surface, elevated, overlapped",
                               ifelse(Structure$ns_e[Structure$record==a]==FALSE,"near surface, elevated, not overlapped",
                                      "near surface, elevated, automatic"))
  site.meta$value[2] <- ifelse(Structure$ns_m[Structure$record==a]==TRUE,"near surface, midstorey, overlapped",
                               ifelse(Structure$ns_m[Structure$record==a]==FALSE,"near surface, midstorey, not overlapped",
                                      "near surface, midstorey, automatic"))
  site.meta$value[3] <- ifelse(Structure$e_m[Structure$record==a]==TRUE,"elevated, midstorey, overlapped",
                               ifelse(Structure$e_m[Structure$record==a]==FALSE,"elevated, midstorey, not overlapped",
                                      "elevated, midstorey, automatic"))
  site.meta$value[4] <- ifelse(Structure$e_c[Structure$record==a]==TRUE,"elevated, canopy, overlapped",
                               ifelse(Structure$e_c[Structure$record==a]==FALSE,"elevated, canopy, not overlapped",
                                      "elevated, canopy, automatic"))
  site.meta$value[5] <- ifelse(Structure$m_c[Structure$record==a]==TRUE,"midstorey, canopy, overlapped",
                               ifelse(Structure$m_c[Structure$record==a]==FALSE,"midstorey, canopy, not overlapped",
                                      "midstorey, canopy, automatic"))
  site.meta$value[6] <- as.numeric(Flora$weight[Flora$record==a & Flora$species=="Litter"])
  site.meta$value[7] <- if (is.null(as.numeric(Flora$diameter[Flora$record==a & Flora$species=="Litter"]))) {
    0.005
  } else {
    as.numeric(Flora$diameter[Flora$record==a & Flora$species=="Litter"])
  }
  site.meta$value[8] <- 0.00025
  site.meta$value[9] <- fLine
  site.meta$value[10] <- slope
  site.meta$value[11] <- temp
  site.meta$value[12] <- dfmc
  site.meta$value[13] <- wind
  
  return(site.meta)
}



#' Builds the dataframe strata.meta from input tables
#' Compatible with MEE inputs
#' 
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. 
#'    Acceptable values are TRUE, FALSE, or blank, where the outcome 
#'    will be decided by the relative stratum heights.
#' nsR, eR, mR, cR. Species richness (number of species) in each stratum
#' @param a The number of the record for the site

buildStrataMeta <- function(Structure, a)
{
  # Collect subsets for record
  st <- Structure[Structure$record==a,]
  
  # CREATE strata.meta
  strata.meta <- data.frame(matrix(NA, nrow=8, ncol = 5))
  
  names(strata.meta) <- c("stratum", "species", "param", "value", "units")
  
  # Fill levelNames
  strata.meta$param[1] <- "levelName"
  strata.meta$value[1] <- "near surface"
  strata.meta$param[2] <- "plantSeparation"
  strata.meta$value[2] <- st$NS[1]
  strata.meta$units[2] <- "m"
  
  strata.meta$param[3] <- "levelName"
  strata.meta$value[3] <- "elevated"
  strata.meta$param[4] <- "plantSeparation"
  strata.meta$value[4] <- st$El[1]
  strata.meta$units[4] <- "m"
  
  strata.meta$param[5] <- "levelName"
  strata.meta$value[5] <- "midstorey"
  strata.meta$param[6] <- "plantSeparation"
  strata.meta$value[6] <- st$Mid[1]
  strata.meta$units[6] <- "m"
  
  strata.meta$param[7] <- "levelName"
  strata.meta$value[7] <- "canopy"
  strata.meta$param[8] <- "plantSeparation"
  strata.meta$value[8] <- st$Can[1]
  strata.meta$units[8] <- "m"
  
  # Find empty strata and remove them
  deleteVector <- vector()
  for (ro in 2:8) {
    if(is.na(strata.meta$value[ro])){
      deleteVector<- c(deleteVector, ro-1, ro)
    }
  }
  if(length(deleteVector)>0){
    strata.meta <- strata.meta[ -deleteVector, ]}
  
  # Number strata
  rows <- as.numeric(nrow(strata.meta))
  num <- 0.5
  for (ro in 1:rows) {
    numb <- ceiling(num)
    strata.meta$stratum[ro] <- numb
    num = num + 0.5
  }
  
  return(strata.meta)
}

#' Builds the dataframe species.values from input tables
#' Building from MEE inputs requires function speciesBuilder
#' 
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - A name for the site
#' species - the name of the species, which will call trait data from 
#' 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - Percent composition of that species in the stratum. 
#'    If absent, all species will be considered equally.
#' base - base height of plant crowns (m)
#' he - he height of plant crowns (m)
#' ht - ht height of plant crowns (m)
#' top - top height of plant crowns (m)
#' w - width of plant crowns (m)
#' Hs - standard deviation of the top height of plant crowns (m)
#' Hr - range of the top height of plant crowns (m)
#' weight - weight in t/ha of fine dead organic material forming 
#'    the surface and suspNS layers
#' diameter - mean diameter of surface and suspNS litter in m
#' @param default.species.params Plant traits database

buildSpeciesValues <- function(Flora, default.species.params, a)
{
  # Collect subsets for site
  fl <- Flora[Flora$record==a & Flora$species != "Litter",] %>%
    mutate("name" = species) %>%
    left_join(default.species.params, by = "name") 
  if (exists("maxW", where = as.environment(fl))) {
    fl <- fl %>%
      select(species, stratum, comp, base, top, he, ht, w, moisture, C.C_rat, G.C_rat, maxW)
  } else {
    fl <- fl %>%
      select(species, stratum, comp, base, top, he, ht, w, moisture, C.C_rat, G.C_rat)
  }
  
  # CREATE species.values
  ro <- as.numeric(nrow(fl))
  species.values <- data.frame(matrix(NA, nrow=ro, ncol = 13))
  
  names(species.values) <- c("stratum", "species", "name", "clumpDiameter", "clumpSeparation", "composition",
                             "deadLeafMoisture", "hc", "he", "hp", "ht", "w", "liveLeafMoisture")
  
  # Enter values
  species.values$stratum <- as.numeric(fl$stratum)
  species.values$name <- as.character(fl$species)
  species.values$liveLeafMoisture <- fl$moisture
  species.values$hc <- fl$base
  species.values$hp <- fl$top
  species.values$he <- fl$he
  species.values$ht <- fl$ht
  species.values$w <- fl$w
  if (exists("maxW", where = as.environment(fl))) {
    species.values$clumpDiameter <- pmin(pmin((as.numeric(fl$top)-as.numeric(fl$base)), as.numeric(fl$w))*fl$C.C_rat, as.numeric(fl$maxW))
  } else {
    species.values$clumpDiameter <- pmin((as.numeric(fl$top)-as.numeric(fl$base)), as.numeric(fl$w))*fl$C.C_rat
  }
  
  
  species.values$clumpSeparation <- fl$G.C_rat*species.values$clumpDiameter
  species.values$composition <- as.numeric(fl$comp)
  
  species.values <- species.values[order(species.values$stratum),]
  species.values$species <- as.numeric(c(1:ro))
  
  # Calculate gaps
  for(a in 1:ro) {
    if(is.na(species.values$he[a])){
      species.values$he[a] <- species.values$hc[a]
    }
    if(is.na(species.values$ht[a])){
      species.values$ht[a] <- species.values$hp[a]
    }
    if(is.na(species.values$liveLeafMoisture[a])){
      species.values$liveLeafMoisture[a] <- 1
    }
  }
  
  # Ensure strata are numbered consecutively
  stratCurr <- as.data.frame(unique(species.values$stratum))
  names(stratCurr) <- c("orig")
  stratCurr$cor <- 1:as.numeric(nrow(stratCurr))
  
  #Check values
  for(b in 1:as.numeric(nrow(stratCurr))) {
    species.values$stratum[species.values$stratum == stratCurr$orig[b]] <- stratCurr$cor[b]
  }
  
  return(species.values)
}

#' Constructs parameter files from imported tables
#'
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. 
#'    Acceptable values are TRUE, FALSE, or blank, where the outcome 
#'    will be decided by the relative stratum heights.
#' nsR, eR, mR, cR. Species richness (number of species) in each stratum
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - A name for the site
#' species - the name of the species, which will call trait data from 
#' 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - Percent composition of that species in the stratum. 
#'    If absent, all species will be considered equally.
#' base - base height of plant crowns (m)
#' he - he height of plant crowns (m)
#' ht - ht height of plant crowns (m)
#' top - top height of plant crowns (m)
#' w - width of plant crowns (m)
#' Hs - standard deviation of the top height of plant crowns (m)
#' Hr - range of the top height of plant crowns (m)
#' weight - weight in t/ha of fine dead organic material forming 
#'    the surface and suspNS layers
#' diameter - mean diameter of surface and suspNS litter in m
#' @param a The number of the record for the site
#' @param fLine The fireline length in m
#' @param slope slope in degrees
#' @param temp Ambient temperature deg. C
#' @param dfmc Moisture content of fine dead fuels in whole numbers (eg 0.1 for 10 percent)
#' @param wind Velocity in km/h
#' @param default.species.params Plant traits database
#' @param a The record number for which to build the table
#' @export

buildParams <- function(Structure, Flora, default.species.params, a,
                        fLine = 100, slope = 0, temp = 30, dfmc = 0.1, wind = 10)
{
  # Remove empty species, leaving litter
  Flora <- Flora[which(Flora$species == "Litter" | Flora$comp != 0), ]
  
  # Construct component tables
  site.meta <- buildSiteMeta(Structure = Structure, Flora = Flora, a = a, fLine = fLine, slope = slope,
                             temp = temp, dfmc = dfmc, wind = wind)
  strata.meta <- buildStrataMeta(Structure, a)
  species.values <- buildSpeciesValues(Flora, default.species.params, a)
  species.values$deadLeafMoisture <- dfmc
  species.units <- unitBuilder(Flora, a)
  
  # Build parameter file
  components <- list(site.meta, strata.meta, species.values, species.units)
  names(components) <- c("site.meta", "strata.meta", "species.values", "species.units")
  param <- ffm_assemble_table(components)
  
  # Find weighted mean of crown widths per stratum
  species.values$weightedW <- species.values$composition * as.numeric(species.values$w)
  ww <- species.values%>%
    select(stratum, composition, weightedW)%>%
    group_by(stratum) %>%
    summarize_all(sum) %>%
    mutate(mw = weightedW/composition)
  
  # Ensure separation is greater than weighted mean of crown widths
  for (stNum in 1:max(ww$stratum)) {
    sep <- filter(strata.meta, stratum == stNum)
    param <- ffm_set_stratum_param(param, stNum, "plantSeparation", 
                                   pmax(as.numeric(sep$value[2]),ww$mw[stNum]))
  }
  
  # Change species param leafForm to characters
  default.species.params$leafForm <- as.character(default.species.params$leafForm)
  
  # Fill empty traits
  param <- ffm_complete_params(param,default.species.params)
  
  return(param)
}

#' Constructs a default.species.params table using traits available in the austraits database
#' Collects traits from ausTraits for use in shade and fire effects modelling
#'
#' @param version Version of austraits
#' @param path 
#' @param shadeTolerance The Whole Plant Light Compensation Point below which 
#' plants are considered shade tolerant
#'
#' @return list of tables
#' @export
#' @examples
#' Traits <- ausTraitTable(version = "3.0.2", path = "data/austraits")
#'

ausTraitTable <- function(version = "3.0.2", path = "data/austraits", shadeTolerance = 3) {
  
  Austraits <- austraits::load_austraits(version = version, path = path)
  
  # Create tables of each trait in FRaME units
  leaf_length  <- (austraits::extract_trait(Austraits, "leaf_length"))$traits %>%
    group_by(taxon_name)%>%
    summarise_if(is.numeric, median) %>%
    mutate(leafLength = value/1000) %>%
    select(taxon_name, leafLength)
  
  leaf_width  <- (austraits::extract_trait(Austraits, "leaf_width"))$traits %>%
    group_by(taxon_name)%>%
    summarise_if(is.numeric, median) %>%
    mutate(leafWidth = value/1000) %>%
    select(taxon_name, leafWidth)
  
  leaf_shape  <- (austraits::extract_trait(Austraits, "leaf_shape"))$traits 
  leaf_shape$leafForm <- case_when(leaf_shape$value == "needle terete" ~ "Round",
                                   leaf_shape$value == "terete" ~ "Round",
                                   leaf_shape$value == "subulate terete" ~ "Round",
                                   leaf_shape$value == "needle" ~ "Round",
                                   leaf_shape$value == "linear terete" ~ "Round",
                                   leaf_shape$value == "clavate linear terete" ~ "Round",
                                   leaf_shape$value == "terete_compressed" ~ "Round",
                                   leaf_shape$value == "scale terete" ~ "Round",
                                   leaf_shape$value == "trigonous_terete" ~ "Round",
                                   leaf_shape$value == "half_terete" ~ "Round",
                                   leaf_shape$value == "terete_compressed" ~ "Round",
                                   TRUE ~ "Flat")
  leaf_shape  <- leaf_shape %>% select(taxon_name, leafForm)
  
  leaf_thickness  <- (austraits::extract_trait(Austraits, "leaf_thickness"))$traits %>%
    group_by(taxon_name)%>%
    summarise_if(is.numeric, median) %>%
    mutate(leafThickness = value/1000) %>%
    select(taxon_name, leafThickness)
  
  # List species
  tNames <- unique(Austraits[[1]]$taxon_name)
  T_names <- tNames[!grepl("\\d+", tNames)]
  
  
  traits <- data.frame(taxon_name = T_names) 
  list_of_dfs <- list(traits, leaf_length, leaf_width, leaf_shape, leaf_thickness)
  result <- list_of_dfs %>%
    purrr::reduce(left_join, by = "taxon_name") %>%
    mutate(name = taxon_name,
           propDead = 0) %>%
    select(name, propDead, leafForm, leafThickness, leafWidth, leafLength) %>%
    mutate(leafSeparation = NA,
           stemOrder = NA,
           ignitionTemp = NA,
           moisture = 1,
           G.C_rat = NA,
           C.C_rat = NA)
  
  result <- arrange(result, name)
  
  # Secondary traits
  traitNames <- data.frame(taxon_name = T_names) 
  
  #___________________________________
  # Shade tolerance
  Dark_respiration  <- (austraits::extract_trait(Austraits, "leaf_dark_respiration_per_dry_mass"))$traits %>%
    group_by(taxon_name)%>%
    summarise_if(is.numeric, median) %>%
    mutate(Dark_respiration = value,
           WPLCP = (390*value)+0.78,
           Shade = WPLCP < shadeTolerance) %>%
    select(taxon_name, WPLCP, Shade)
  # Approximate from area measurement
  respArea  <- (austraits::extract_trait(Austraits, "leaf_dark_respiration_per_area"))$traits %>%
    group_by(taxon_name)%>%
    summarise_if(is.numeric, median)
  list_of_dfs <- list(traitNames, Dark_respiration, respArea)
  sh <- list_of_dfs %>%
    purrr::reduce(left_join, by = "taxon_name")
  sh <- sh[complete.cases(sh),]
  LM<-lm(sh$value ~ sh$WPLCP)
  LMSum <- base::summary(LM)
  eqThres <- LMSum$coefficients[1] + LMSum$coefficients[2] * shadeTolerance
  respArea <- respArea %>%
    mutate(Shade = value < eqThres)
  ShadeTolerance <- respArea[!(respArea$taxon_name %in% Dark_respiration$taxon_name), ]%>%
    right_join(Dark_respiration, by = c("taxon_name", "Shade")) %>%
    select(taxon_name, Shade) 
  
  list_of_dfs <- list(traitNames, ShadeTolerance)
  resultII <- list_of_dfs %>%
    purrr::reduce(left_join, by = "taxon_name") %>%
    mutate(Genus = word(taxon_name, 1),
           name = taxon_name)
  
  # Generalise genera
  genera <- resultII[complete.cases(resultII), ] %>%
    group_by(Genus)%>%
    summarise(Shade = as.logical(round(mean(Shade),0))) %>%
    mutate(name = Genus) %>%
    select(name, Shade)
  
  resultII <- resultII %>%
    select(name, Shade)%>%
    full_join(genera, by = c("name", "Shade"))
  resultII <- arrange(resultII, name)
  
  remList <- vector()
  for (r in 2:nrow(resultII)) {
    if (tolower(resultII$name[r]) == tolower(resultII$name[r-1])) {
      remList[length(remList)+1] <- r-1
    }
  }
  if (length(remList)>0) {
    resultII <- resultII[-remList,]
  }
  
  return(list(result, resultII))
}


#' Checks traits table for column names and adds data from a new table
#'
#' @param traits Existing trait table
#' @param traitsNew Table containing new trait data
#' @param deleteReplicates If set to TRUE, retains the first mention of a species, then removes all following mentions
#' @param fill If set to TRUE, fills empty numeric values with mean values of the genus
#'
#' @return table
#' @export
#' @examples 
#' Traits <- ausTraitTable(version = "3.0.2", path = "data/austraits")
#' TraitsNew <- read.csv("Traits.csv")
#' Tr <- updateTraits(traits = Traits, traitsNew = TraitsNew, deleteReplicates = TRUE, fill = TRUE)
#'

updateTraits <- function(traits, traitsNew, deleteReplicates = TRUE, printReplicates = TRUE, fill = TRUE) {
  # Ensure traits table has all necessary columns
  cols_to_check <- c("name", "propDead", "leafForm", "leafThickness", "leafWidth", 
                     "leafLength", "leafSeparation", "stemOrder", "ignitionTemp", "moisture", "G.C_rat", "C.C_rat")
  
  for (colName in cols_to_check) {
    if(!colName %in% names(traits)) {
      traits$colName == ""
    }
  }
  
  # Add any missing species to existing list
  missing_species <- unique(traitsNew$name)[!(unique(traitsNew$name) %in% unique(traits$name))]
  newSp <- traitsNew[traitsNew$name %in% missing_species, ]
  traits <- rbind(traits,newSp)
  traits <- arrange(traits, name)
  
  newCols <- names(traitsNew)
  
  for (sp in unique(traitsNew$name)) {
    for (colName in newCols) {
      if(!is.na(traitsNew[traitsNew$name == sp,colName])){traits[traits$name == sp,colName] <- traitsNew[traitsNew$name == sp,colName]}
    }
  }
  # Generalise genera
  traits <- traits %>% 
    mutate(Genus = word(name, 1),
           leafThickness = as.numeric(leafThickness),
           leafWidth = as.numeric(leafWidth),
           leafLength = as.numeric(leafLength),
           leafSeparation = as.numeric(leafSeparation),
           stemOrder = as.numeric(stemOrder),
           ignitionTemp = as.numeric(ignitionTemp),
           moisture = as.numeric(moisture),
           G.C_rat = as.numeric(G.C_rat),
           C.C_rat = as.numeric(C.C_rat))
  genera <- traits  %>%
    group_by(Genus)%>%
    summarise_if(is.numeric, ~mean(., na.rm=TRUE)) %>%
    mutate(leafForm = "Flat")
  
  # Find most common leaf form
  for (g in 1:nrow(genera)) {
    genSub <- traits[traits$Genus == genera$Genus[g],]
    Round <- nrow(genSub[genSub$leafForm == "Round",])
    Flat <- nrow(genSub[genSub$leafForm == "Flat",])
    if (Round > Flat) {
      genera$leafForm[g] <- "Round"
    }
  }
  
  # Fill empty traits with means if chosen
  if (fill) {
    traits <- left_join(traits, genera, by = "Genus")%>%
      mutate(propDead = ifelse(!is.na(propDead.x), propDead.x, propDead.y),
             leafThickness = ifelse(!is.na(leafThickness.x), leafThickness.x, leafThickness.y),
             leafForm = ifelse(!is.na(leafForm.x), leafForm.x, leafForm.y),
             leafWidth = ifelse(!is.na(leafWidth.x), leafWidth.x,leafWidth.y),
             leafLength = ifelse(!is.na(leafLength.x), leafLength.x, leafLength.y),
             leafSeparation= ifelse(!is.na(leafSeparation.x), leafSeparation.x,leafSeparation.y),
             stemOrder = ifelse(!is.na(stemOrder.x), stemOrder.x, stemOrder.y),
             ignitionTemp = ifelse(!is.na(ignitionTemp.x), ignitionTemp.x, ignitionTemp.y),
             moisture = ifelse(!is.na(moisture.x), moisture.x, moisture.y),
             G.C_rat = ifelse(!is.na(G.C_rat.x), G.C_rat.x, G.C_rat.y),
             C.C_rat = ifelse(!is.na(C.C_rat.x), C.C_rat.x, C.C_rat.y))
  }
  
  traits <- traits %>%
    dplyr::select(name, propDead, leafForm, leafThickness, leafWidth, leafLength, leafSeparation, stemOrder, ignitionTemp, moisture, G.C_rat, C.C_rat)
  genera <- genera %>%
    mutate(name = Genus) %>%
    dplyr::select(name, propDead, leafForm, leafThickness, leafWidth, leafLength, leafSeparation, stemOrder, ignitionTemp, moisture, G.C_rat, C.C_rat)
  traits <- genera  %>%
    rbind(traits)
  
  traits$propDead[which(is.nan(traits$propDead))] <- NA
  traits$leafThickness[which(is.nan(traits$leafThickness))] <- NA
  traits$leafWidth[which(is.nan(traits$leafWidth))] <- NA
  traits$leafLength[which(is.nan(traits$leafLength))] <- NA
  traits$leafSeparation[which(is.nan(traits$leafSeparation))] <- NA
  traits$stemOrder[which(is.nan(traits$stemOrder))] <- NA
  traits$ignitionTemp[which(is.nan(traits$ignitionTemp))] <- NA
  traits$moisture[which(is.nan(traits$moisture))] <- NA
  traits$G.C_rat[which(is.nan(traits$G.C_rat))] <- NA
  traits$C.C_rat[which(is.nan(traits$C.C_rat))] <- NA
  traits <- arrange(traits,name)
  
  if (deleteReplicates) {
    traits <- remReps(traits, printReplicates)
  }
  
  return(traits)
}


#' Removes replicates in a formatted trait table
#'
#' @param traits the table default.species.params
#' @param printReplicates Set TRUE if names of species removed are printed
#'
#' @return
#' @export
#'
#' 
remReps <- function(traits, printReplicates = TRUE){
  
  if(any(table(tolower(traits$name)) > 1)) {
    repNames <- names(which(table(tolower(traits$name)) > 1))
    remList <- vector()
    n <- 1
    for (s in repNames) {
      recs <- length(which(tolower(traits$name) == s))-2
      remList[n:(n+recs)] <- which(tolower(traits$name) == s)[-1]
      n<-n+1+recs
    }
    if (length(remList)>0) {
      traits <- traits[-remList,]
      cat(length(remList), "rows of duplicate species were removed.\n")
      if (printReplicates) {
        print(repNames)
      }
    }
  }
  return(traits)
}