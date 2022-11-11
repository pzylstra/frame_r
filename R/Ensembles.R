#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#'
#' @param base.params Parameter input table
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export
#' @examples
#' 
#' SPECIFY INPUTS
#' record <- 1
#' data(site)
#' data(structure)
#' data(flora)
#' data(traits)
#' data(Weather)
#' base.params <- paramBuilder(site, structure, flora, traits, record)
#' 
#' RUN WEATHERSET
#' weatherSet(base.params, weather, db.path = "out.db", jitters = 5, l = 0.1, 
#' Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41, updateProgress = NULL)
#' 
#' LOAD AND ORGANISE RESULTS
#' res<-ffm_db_load("out.db")

weatherSet <- function(base.params, weather, db.path = "out_mc.db", jitters = 10, l = 0.1,
                       Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41,updateProgress = NULL)
{
  # Collect initial veg descriptors
  
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  # Run the model, updating the base parameter table
  # with MC values at each iteration
  
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  for (i in 1:max(weather$tm)) {
    
    # Read weather values from the table
    
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Update parameter table
    tbl <- base.params %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        ## Create database and delete the last part
        db.recreate <- (i * j) == 1
        
        tblMod <- plantVar(tbl, Strata, Species, l = l,
                        Ms = Ms, Pm = Pm, Mr = Mr, Hs = Hs, Hr = Hr)
        # Run the model
        ffm_run(tblMod, db.path, db.recreate = db.recreate)
      }
    }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ",max(weather$tm) - i )
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
  
  
  cat("Finished.  Output written to", db.path)
}

#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#' using plantVarS to modify individual species with their own measured variability
#'
#' @param base.params Parameter input table
#' @param a A unique identifier for the record being run
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

weatherSetS <- function(base.params, weather, Variation, Structure, a, db.path = "out_mc.db", jitters = 10, l = 0.1,
                        Ms = 0.01, Pm = 1, Mr = 1.001, updateProgress = NULL)
{
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  
  # Read weather values from the table
  for (i in 1:max(weather$tm)) {
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Select species for a random point in the forest and import the weather parameters
    tbl <- specPoint(base.params, Structure, a) %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    Strata <- strata(tbl)
    Species <- species(tbl)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        # Recreate database on first run, then add following runs to this
        db.recreate <- (i * j) == 1
        
        # Vary plant traits for each species within their range
        tbl <- plantVarS(tbl, Strata, Species, Variation, a, l = l,
                         Ms = Ms, Pm = Pm, Mr = Mr)
        # Run the model
        ffm_run(tbl, db.path, db.recreate = db.recreate)
      }
    }
    
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ",max(weather$tm) - i )
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
  
  cat("Finished.  Output written to", db.path)
}


#' Models fires using sites constructed from imported tables
#' 
#' Private function in development

#' @param site A dataframe with the six fields:
#' record - a unique, consecutively numbered identifier per site
#' site - the name of each record
#' slope - slope in degrees
#' wind - velocity in km/h
#' temp - ambient temperature deg. C
#' dfmc - moisture content of fine dead fuels in whole numbers (eg 0.1 for 10%)
#' litter - weight in t/ha of fine dead organic material forming the O horizon
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
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param traits A dataframe of plant traits, with the fields:
#' name - the name of the species
#' propDead - the proportion (0-1) of the foliage that is dead
#' leafForm - allowable values are 'Flat' or 'Round' for terete or needle leaves
#' leafThickness, leafWidth, leafLength, leafSeparation - leaf dimensions (m)
#' stemOrder - stem ramification
#' ignitionTemp - minimum temperature at which the leaf can ignite (deg C)
#' @param default.species.params Leaf traits database
#' @return dataframe

fireSet <- function(site, Structure, Flora, traits = default.species.params)
{
  # Number of sites to model
  sn <- max(site$record)
  
  # Enter trait database
  default.species.params <- traits
  
  # Initial record
  # Build parameter files
  param <- paramBuilder(site, Structure, Flora, 1)
  #Run the model
  ffm_run(param, "Behav.db")
  
  #Load and display the outputs
  res<-ffm_db_load("Behav.db")
  
  #Build tables
  surf <- surf(res$SurfaceResults)%>%
    mutate(fireNo = 1)
  x <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)%>%
    mutate(fireNo = 1)
  runs <- frameSummary(x, surf)%>%
    mutate(cat = round(wind_kph,0),
           fireNo = 1)
  IP <- repFlame(res$IgnitionPaths)%>%
    mutate(fireNo = 1)
  
  if(sn >1){
    
    # Loop through records
    for (a in 2:sn) {
      
      # Build parameter files
      param <- paramBuilder(site, Structure, Flora, a)
      #Run the model
      ffm_run(param, "Behav.db")
      
      # Build a table with a row for each output. Consider a script for different levels of analysis (eg simple to comprehensive),
      # Then list these in If Then scenarios that can be identified in the initial function.
      # Build wider level functions for weather, growth etc that can call this one.
      
      #Load and display the outputs
      res<-ffm_db_load("Behav.db")
      
      #Build tables
      surfa <- surf(res$SurfaceResults)%>%
        mutate(fireNo = a)
      xa <- stratum(res$FlameSummaries, res$Sites, res$ROS, surf)%>%
        mutate(fireNo = a)
      runsa <- frameSummary(xa, surfa)%>%
        mutate(cat = round(wind_kph,0),
               fireNo = a)
      IPa <- repFlame(res$IgnitionPaths)%>%
        mutate(fireNo = a)
      
      # Append to existing table
      surf <- rbind(surf, surfa)
      x <- rbind(x, xa)
      runs <- rbind(runs, runsa)
      IP <- rbind(IP, IPa)
    }
  }
  # Export csv files
  write.csv(surf,"Surface.csv")
  write.csv(IP,"IP.csv")
  write.csv(x,"All.csv")
  return(runs)
}



#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#' @param base.params Parameter input table
#' @param Strata Strata descriptor table output by the function 'strata'
#' @param Species Species descriptor table output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Hs Standard deviation of plant height variations
#' @param Hr Truncates plant height variability by +/- Hr * height
#' @return dataframe
#' @export

plantVar <- function (base.params, Strata, Species,
                      l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001, Hs = 0.2, Hr = 1.41)
{
  test <- FALSE
  while (test == FALSE) {
    tbl <- base.params
    
    # 1. VARY LEAF TRAITS
    tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
    
    # 2. VARY STRUCTURE
    SpeciesN <- 1
    SpeciesP <- 1
    StN <- as.numeric(count(Strata))
    for (si in 1:StN) {
      
      # 2a. Determine which strata are present
      if (runif(1) <= Strata$cover[si]) {
        
        # 2b. Vary moisture in species present
        for (t in 1:Strata$speciesN[si]) {
          Mrand <- Pm * rtnorm(n = 1, mean = Species$lfmc[SpeciesN],
                               sd = Ms, a = Species$lfmc[SpeciesN]/Mr, b = Species$lfmc[SpeciesN] *
                                 Mr)
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", Mrand)
          (SpeciesN = SpeciesN + 1)
        }
      }
      
      # 2c. Mark species for removal if stratum is not present
      else {
        for (f in 1:Strata$speciesN[si]) {
          tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                               "liveLeafMoisture", 10000)
          SpeciesN = SpeciesN + 1
        }
      }
      
      # 2d. Vary plant crown dimensions
      for (p in 1:Strata$speciesN[si]) {
        peak <- rtnorm(n = 1, mean = Species$hp[SpeciesP],
                       sd = Hs, a = Species$hp[SpeciesP]/Hr, b = Species$hp[SpeciesP] *
                         Hr)
        tbl <- tbl %>%
          ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
          ffm_set_species_param(si, SpeciesP, "ht", peak * Species$htR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "he", peak * Species$heR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "hc", peak *Species$hcR[SpeciesP]) %>%
          ffm_set_species_param(si, SpeciesP, "w", peak * Species$wR[SpeciesP])
        SpeciesP = SpeciesP + 1
      }
    }
    x <- tbl$species[tbl$value==10000]
    x <- x[!is.na(x)]
    test <- length(x) < max(Species$sp)
  }
  
  # 2e. Remove species and strata not present
  x <- tbl$species[tbl$value==10000]
  x <- x[!is.na(x)]
  if (length(x)>0) {
    for (sp in 1:length(x)) {
      tbl <- subset(tbl,species!=x[sp] | is.na(species))
    }
    
    # Remove the empty strata
    Strat <- strata(tbl)
    s <- unique(Strat$stratum)
    sequence <- unique(tbl$stratum)
    sequence <- sequence[!is.na(sequence)]
    x <- sequence[!sequence %in% s]
    for (st in 1:length(x)) {
      tbl <- subset(tbl,stratum!=x[st] | is.na(stratum))
    }
    
    # Renumber strata
    sList<- data.frame('stratum' = s, 'n' = 1:length(s))
    tbl <- left_join(tbl,sList, by = "stratum") %>%
      mutate(stratum = n) %>%
      select(!n)
    
    # Renumber species
    s <- unique(tbl$species)
    s <- s[!is.na(s)]
    sList<- data.frame('species' = s, 'n' = 1:length(s))
    tbl <- left_join(tbl,sList, by = "species") %>%
      mutate(species = n) %>%
      select(!n)
    
    # Rename strata
    s <- unique(tbl$stratum)
    s <- s[!is.na(s)]
    tbl$value[tbl$stratum==s[1]&tbl$param=="levelName"] = "near surface"
    tbl$value[tbl$stratum==max(s)&tbl$param=="levelName"] = "canopy"
  }
  # Remove marked species and tidy params file
  tbl <- speciesDrop(tbl)
  return(tbl)
}

#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#' Differs from plantVar by modifying individual species by their own rules
#' @param base.params Parameter input table
#' @param a A unique identifier for the record being run
#' @param Strata Strata descriptor table output by the function 'strata'
#' @param Species Species descriptor table output by the function 'species'
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @return dataframe
#' @export

plantVarS <- function (base.params, Strata, Species, Variation, a, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001)
{
  
  tbl <- base.params
  # Vary leaf traits
  tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
  SpeciesN <- 1
  SpeciesP <- 1
  
  #Filter variation table to record
  varRec <- Variation[Variation$record ==a, ]
  varRec <- varRec[order(varRec$stratum),]
  
  # Loop through plant strata
  
  StN <- as.numeric(count(Strata))
  
  for (si in 1:StN) {
    
    # Vary leaf moisture to randomly place points in the community and decide which plants will be present.
    # Where the random number is > stratum cover, all species are made too moist to burn and excluded from modelling.
    # Otherwise, species are varied by Ms & Mr, and multiplied by Pm
    
    if (runif(1) <= Strata$cover[si]) {
      for (t in 1:Strata$speciesN[si]) {
        Mrand <- Pm * rtnorm(n = 1, mean = Species$lfmc[SpeciesN],
                             sd = Ms, a = Species$lfmc[SpeciesN]/Mr, b = Species$lfmc[SpeciesN] * Mr)
        tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                             "liveLeafMoisture", Mrand)
        (SpeciesN = SpeciesN + 1)
      }
    } else {
      for (f in 1:Strata$speciesN[si]) {
        tbl <- tbl %>% ffm_set_species_param(si, SpeciesN,
                                             "liveLeafMoisture", 100)
        SpeciesN = SpeciesN + 1
      }
    }
    
    # Modify plant dimensions for each species within the stratum
    
    for (p in 1:Strata$speciesN[si]) {
      Hr <- varRec$Hr[SpeciesP]
      peak <- rtnorm(n = 1, mean = Species$hp[SpeciesP],
                     sd = varRec$Hs[SpeciesP], a = Species$hp[SpeciesP]/Hr, b = Species$hp[SpeciesP] *
                       Hr)
      tbl <- tbl %>%
        ffm_set_species_param(si, SpeciesP, "hp", peak) %>%
        ffm_set_species_param(si, SpeciesP, "ht", peak * Species$htR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "he", peak * Species$heR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "hc", peak *Species$hcR[SpeciesP]) %>%
        ffm_set_species_param(si, SpeciesP, "w", peak * Species$wR[SpeciesP])
      SpeciesP = SpeciesP + 1
    }
  }
  # Remove marked species and tidy params file
  tbl <- speciesDrop(tbl)
  
  return(tbl)
}


#' Randomly modifies plant traits within defined ranges for non-deterministic predictions
#' 
#' Differs from plantVarS by using FRaME-updated tables
#' 
#' @param base.params Parameter input table
#' @param a A unique identifier for the record being run
#' @param Strata Strata descriptor table output by the function 'strata'
#' @param Species Species descriptor table output by the function 'species'
#' @param Flora  A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @return dataframe
#' @export
#' 
plantVarFrame <- function (base.params, Strata, Species, Flora, a, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001)
{
  
  tbl <- base.params
  # Vary leaf traits
  tbl <- ffm_param_variance(tbl, max.prop = l, method = "uniform")
  SpeciesN <- 1
  SpeciesP <- 1
  
  #Filter variation table to record
  varRec <- Flora[Flora$record==a & Flora$species!="Litter",]
  varRec$Hs[is.na(varRec$Hs)]<-0.001
#  varRec$Hr[varRec$Hr==0]<-0.001 
  varRec$Hr<-pmin(2,varRec$Hr+1) 
  varRec$Hs[varRec$Hs==0]<-0.001 
#  varRec$Hr<-as.numeric(varRec$Hr)+1 
  varRec <- varRec %>% 
    mutate(name = species,
           st = as.double(stratum)) %>%
    select(st, name, Hs, Hr)
  varRec <- left_join(Species, varRec, by = c("st", "name"))
  
  # Loop through plant strata
  
  StN <- as.numeric(count(Strata))
  dCount <- 0
  
  for (st in 1:StN) {
    
    # Vary leaf moisture to randomly place points in the community and decide which plants will be present.
    # Where the random number is > stratum cover, all species are marked with liveLeafMoisture == 100.
    # Otherwise, species are varied by Ms & Mr, and multiplied by Pm
    
    if (runif(1) <= Strata$cover[st] | dCount == StN - 1) {
      for (t in 1:Strata$speciesN[st]) {
        Mrand <- Pm * rtnorm(n = 1, mean = Species$lfmc[SpeciesN],
                             sd = Ms, a = Species$lfmc[SpeciesN]/Mr, b = Species$lfmc[SpeciesN] * Mr)
        tbl <- ffm_set_species_param(tbl, st, SpeciesN,
                                     "liveLeafMoisture", Mrand)
        (SpeciesN <- SpeciesN + 1)
      }
    } else {
      dCount <- dCount+1
      for (f in 1:Strata$speciesN[st]) {
        tbl <- tbl %>% ffm_set_species_param(st, SpeciesN,
                                             "liveLeafMoisture", 100)
        SpeciesN <- SpeciesN + 1
      }
    }
    
    # Modify plant dimensions for each species within the stratum
    
    for (p in 1:Strata$speciesN[st]) {
      Hr <- as.numeric(varRec$Hr[SpeciesP])
      peak <- rtnorm(n = 1, mean = Species$hp[SpeciesP],
                     sd = as.numeric(varRec$Hs[SpeciesP]), a = Species$hp[SpeciesP]/Hr, b = Species$hp[SpeciesP]*Hr)
      
      # APPROXIMATE FIX ADDED TO DEAL WITH SCALA CODE FAILING
      # https://github.com/pzylstra/frame_scala/blob/master/forest/src/main/scala/ffm/forest/DefaultVegetationWindModel.scala#L50
      if (st == max(Species$st)) {
        peak <- max(peak, 1)
      }
      tbl <- tbl %>%
        ffm_set_species_param(st, SpeciesP, "hp", peak) %>%
        ffm_set_species_param(st, SpeciesP, "ht", peak * Species$htR[SpeciesP]) %>%
        ffm_set_species_param(st, SpeciesP, "he", min((peak * Species$heR[SpeciesP]),(peak * Species$htR[SpeciesP]))) %>%
        ffm_set_species_param(st, SpeciesP, "hc", min(0.9*peak, peak * Species$hcR[SpeciesP])) %>%
        ffm_set_species_param(st, SpeciesP, "w", peak * Species$wR[SpeciesP])
      SpeciesP = SpeciesP + 1
    }
  }
  
  # Remove marked species and tidy params file
  tbl <- speciesDrop(tbl)
  
  # Check and correct for stratum overlaps
  Strata <- strata(tbl)
  Species <- species(tbl)
  for (row in as.numeric(Strata$stratum)) {
    if (row > 1) {
      # Check for overlap between base of stratum and top of lower stratum
      if (Strata$base[row] < Strata$top[row-1]) {
        # Check that the lower stratum is not as tall as the upper stratum
        if (Strata$top[row-1] < Strata$top[row]) {
          # Correct overlapping strata to average of heights
          newV <- (Strata$base[row]+Strata$top[row-1])/2
          # Make sure the new heights don't make the heights of the lower stratum conflict
          if (newV > max(Species$hc[Species$st==row-1])) {
            Species$hc[Species$st == row] <- newV
            Species$he[Species$st == row] <- newV
            Species$hp[Species$st == row-1] <- newV
            Species$ht[Species$st == row-1] <- newV
          } else {
            cat("The base of stratum", row, "is too low to correct.", "\n")
          }
          # Loop through and update species in each stratum
          for (sp in unique(Species$sp)) {
            tbl$value[tbl$param == "hc" & tbl$species == sp] <- Species$hc[sp]
            tbl$value[tbl$param == "he" & tbl$species == sp] <- Species$he[sp]
            tbl$value[tbl$param == "ht" & tbl$species == sp] <- Species$ht[sp]
            tbl$value[tbl$param == "hp" & tbl$species == sp] <- Species$hp[sp]
          }
        } else {
          cat("Stratum", row, "is as tall as the next stratum")
        }
      }
    }
  }
  
  return(tbl)
}


#' Selects random species from the available list, weighted by their frequency
#' Modifies a parameter table to the shortened list
#'
#' @param base.params A parameter file
#' @param a A unique identifier for the record being run
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @return dataframe
#' @export

specPoint <- function(base.params, Structure, a)
{
  
  Species <- species(base.params)
  
  # Add temporary fields
  Species$wComp = 0
  Species$include = 1
  
  # Count strata
  StN <- as.numeric(max(base.params$stratum[!is.na(base.params$stratum)]))
  
  
  # For each stratum, identify the species being considered, then choose how many of these
  # will be modelled from a range set by the recorded maximum point richness of each stratum
  
  richList <- c(as.numeric(Structure[Structure$record == a, ]$nsR), as.numeric(Structure[Structure$record == a, ]$eR),
                as.numeric(Structure[Structure$record == a, ]$mR), as.numeric(Structure[Structure$record == a, ]$cR))
  richList <- richList[!is.na(richList)]
  
  for (StratNo in 1:StN) {
    
    SpL <- as.numeric(min(Species[Species$st == StratNo, ]$sp))
    SpU <- as.numeric(max(Species[Species$st == StratNo, ]$sp))
    SpN <- SpU-SpL+1
    
    # Species richness for point in the stratum
    R <- richList[StratNo]
    choose <- round(runif(n=1)*(min(R,SpN)-1),0)+1
    
    # Limit stratum species list to a random selection weighted by species occurrence
    for (a in SpL:SpU) {
      Species$wComp[a] = runif(n=1)*Species$comp[a]
    }
    # Identify unneeded records
    low <- Rfast::nth(Species[Species$st == StratNo, ]$wComp, choose, descending = TRUE)
    
    for (sp in SpL:SpU) {
      Species$include[sp] = if (Species$wComp[sp] < low) { 0 }
      else {  1 }
    }
  }
  Species <- Species%>%
    mutate(species = as.character(sp))
  Species[,"new"] <- cumsum(Species$include)
  base.params <- base.params %>%
    mutate(species = as.character(species))
  
  param <- left_join(base.params, Species, by="species")%>%
    subset(include != 0 | is.na(include))%>%
    mutate(species = new)%>%
    select(stratum, species, param, value, units)
  rownames(param) <- seq(length=nrow(param))
  
  return(param)
}



#' Models fire behaviour in each species of a parameter table,
#' summarising the combustibility from the length of flame divided
#' by the length of segment ignited
#'
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
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param default.species.params Leaf traits database
#' @param a The record number for which to build the table
#' @return dataframe
#' @export

spComb <- function(Flora, Structure, default.species.params, a)
{
  
  #Build standardised tables and parameter file
  Si <- as.data.frame(list('record' = a, 'site' = "test", 'slope'=0, 'wind'=5, 'temp'=30, 'dfmc'=0.05, 'litter'=5, 'fLine'=1))
  
  f <- Flora[Flora$record==a, ]
  f$base = 0.25
  f$he = 0.25
  f$ht = 1.25
  f$top = 1.25
  f$w = 1
  
  s <- Structure[Structure$record==a, ]
  s$NS = 100
  s$El = 100
  s$Mid = 100
  s$Can = 100
  
  base.params <- paramBuilder(site = Si, Structure = s, Flora = f, default.species.params, a) 
  
  specflam <- function(base.params, st, sp) {
    
    #Subset tab to relevant species
    SP <- base.params[base.params$stratum==st & base.params$species==sp, ]
    SP <- rbind(SP, c(st, NA, "levelName", "canopy", NA))
    SP <- rbind(SP, c(st, NA, "plantSeparation", "100", "m"))
    SP <- SP[!is.na(SP$stratum),]
    end <- base.params %>% filter(is.na(stratum))
    tab <- rbind(SP, end)
    
    for (i in 5:15) {
      db.recreate <- i == 5
      tab$value[tab$param == "fuelLoad"] = i
      ffm_run(params = tab, db.path = "out.db", db.recreate = db.recreate)
    }
    res <- ffm_db_load("out.db")
    IP <- repFlame(res$IgnitionPaths) %>% mutate(rat = flameLength/length)
    name <- (tab[which(tab$param == "name"), ])$value[1]
    flam <- round(mean(IP$rat), 2)
    out <- as.data.frame(list(Stratum = st, Species = sp, 
                              name = name, Combustibility = flam))
    return(out)
  }
  
  #Create list of species & strata
  can <- base.params %>%
    filter(!is.na(stratum)) %>%
    select(stratum, species)
  candidates <- distinct(can) %>%
    filter(!is.na(species))
  
  n <- as.numeric(nrow(candidates))
  combustibility <- specflam(base.params, candidates$stratum[1], candidates$species[1])
  cat("Measured 1 of", n, "species")
  
  for (j in 2:n){
    summ <- specflam(base.params, candidates$stratum[j], candidates$species[j])
    combustibility <- rbind(combustibility, summ)
    cat("Measured", j, "of", n, "species")
  }
  return(combustibility)
}



#' Models fire behaviour across defined slopes and DFMCs, 
#' with varied plants and random winds and within a defined range
#' 
#' Private function under development
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param replicates Number of wind replicates per slope/DFMC combination
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windMax Maximum wind velocity tested (km/h)
#' @param slopes A string of slope values to be tested
#' @param DFMCs A string of dead fuel moisture values to be tested
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param heightSD Standard deviation of plant height
#' @param heightRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param updateProgress Progress bar for use in the dashboard
 

conDrivers <- function (base.params, db.path = "out_mc.db", replicates, windMin, 
windMax, slopes= c(0, 10), DFMCs = c(0.05, 0.1, 0.15), moistureMultiplier, moistureSD, moistureRange, 
heightSD, heightRange, leafVar, updateProgress = NULL) 
{
  Strata <- strata(base.params)
  Species <- species(base.params)
  dat <- expand.grid(slope = slopes, DFMC = DFMCs)
  Niter <- nrow(dat) * replicates
  pbar <- txtProgressBar(max = Niter, style = 3)
  
  count <- 1
  for (w in 1:replicates) {
    for (i in 1:nrow(dat)) {
      db.recreate <- count == 1
      s <- dat[i, "slope"]
      d <- dat[i, "DFMC"]
      w <- runif(1) * (windMax - windMin) + windMin
      base.params <- base.params %>% 
        ffm_set_site_param("slope",s, "deg") %>% 
        ffm_set_site_param("temperature", 30) %>% 
        ffm_set_site_param("deadFuelMoistureProp", d) %>% 
        ffm_set_site_param("windSpeed", w)
      base.params <- plantVar(base.params, Strata, Species, 
                              l = leafVar, Ms = moistureSD, Pm = moistureMultiplier, 
                              Mr = moistureRange, Hs = heightSD, Hr = heightRange)
      ffm_run(base.params, db.path, db.recreate = db.recreate)
      setTxtProgressBar(pbar, count)
    }
  }
}

#' Removes identified species from a param file and tidies the file
#' 
#' Species where liveLeafMoisture == 100 are removed, then
#' empty strata removed and species & strata renumbered consecutively
#' 
#' @param base.params The params file to be adjusted

speciesDrop <- function(base.params) {
  
  # 1. List table components
  cutList <- as.vector(base.params[base.params$param == "liveLeafMoisture" & base.params$value == 100,])$species
  
  if (length(cutList) > 0) {
    
    end <- base.params[is.na(base.params$stratum),]
    strDesc <- base.params[is.na(base.params$species) & !is.na(base.params$stratum),]
    
    # 2. Remove identified species
    for (sp in cutList) {
      base.params <- base.params[which(base.params$species != sp),]
    }
    
    # 3. Rebuild table
    stList <- unique(base.params$stratum[!is.na(base.params$stratum)])
    for (st in stList) {
      Rest <- strDesc[which(strDesc$stratum == st),]
      base.params <- rbind(base.params, Rest)
    }
    
    # Renumber strata
    stNew <- as.data.frame(stList)
    if (nrow(stNew) > 0) {
      stNew$new = 1:nrow(stNew)
      base.params <- base.params %>%
        left_join(stNew, by = c("stratum" = "stList")) %>%
        mutate(stratum = new) %>%
        select(!"new")
    }
    
    # Renumber species
    spList <- unique(base.params$species)
    spList <- spList[!is.na(spList)]
    spNew <- as.data.frame(spList)
    if (nrow(spNew) > 0) {
      spNew$new = 1:nrow(spNew)
      base.params <- base.params %>%
        left_join(spNew, by = c("species" = "spList")) %>%
        mutate(species = new) %>%
        select(!"new")
    }
    
    # Set stratum names
    M <- as.integer(max(base.params$stratum, na.rm = TRUE))
    base.params[which(base.params$stratum == M & base.params$param == "levelName"),]$value <- "canopy"
    
    if (nrow(stNew) > 1) {
      m <- as.integer(min(base.params$stratum, na.rm = TRUE))
      base.params[which(base.params$stratum == m & base.params$param == "levelName"),]$value <- "near surface"
    }
    
    if (nrow(stNew) > 2) {
      m <- sort(as.integer(unique(base.params$stratum[!is.na(base.params$stratum)])))[2]
      base.params[which(base.params$stratum == m & base.params$param == "levelName"),]$value <- "elevated"
    }
    
    if (nrow(stNew) > 3) {
      m <- sort(as.integer(unique(base.params$stratum[!is.na(base.params$stratum)])))[3]
      base.params[which(base.params$stratum == m & base.params$param == "levelName"),]$value <- "midstorey"
    }
    
    # Restore site data and reorder
    base.params <- rbind(base.params, end)%>%
      arrange(as.integer(stratum), as.integer(species)) %>%
      
      mutate_all(dplyr::funs(as.character))
  }
  
  return(base.params)
}



#' Models probabilistic fire behaviour
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param slope Mean slope (deg)
#' @param slopeSD Standard deviation of slope
#' @param slopeRange Truncates variability by +/- mean * range
#' @param temp Mean ambient temperature (deg.C)
#' @param tempSD Standard deviation of temp
#' @param tempRange Truncates variability by +/- mean * range
#' @param DFMC Mean DFMC (Percent ODW)
#' @param DFMCSD Standard deviation of DFMC
#' @param DFMCRange Truncates variability by +/- mean * range
#' @param wind Mean wind velocity (km/h)
#' @param windSD Standard deviation of wind velocity
#' @param windRange Truncates variability by +/- mean * range
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param heightSD Standard deviation of plant height
#' @param heightRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export
#' @examples
#' 
#' SPECIFY INPUTS
#' record <- 1
#' data(site)
#' data(structure)
#' data(flora)
#' data(traits)
#' base.params <- paramBuilder(site, structure, flora, traits, record)
#' 
#' MODEL PROBABILISTIC FIRE BEHAVIOUR
#' probFire(base.params, db.path = "out.db", jitters = 50,
#'          slope = 0, slopeSD = 2, slopeRange = 5, 
#'          temp = 15, tempSD = 5, tempRange = 3,
#'          DFMC = 0.1, DFMCSD = 0.01, DFMCRange = 2, 
#'          wind = 5, windSD = 1, windRange = 2,
#'          moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.5,
#'          heightSD = 2, heightRange = 1.41, 
#'          leafVar = 0.1,
#'          updateProgress = NULL)
#' 
#' LOAD AND ORGANISE RESULTS
#' res<-ffm_db_load("out.db")

probFire <- function(base.params, db.path = "out_mc.db", jitters,
                     slope, slopeSD, slopeRange, temp, tempSD, tempRange,
                     DFMC, DFMCSD, DFMCRange, wind, windSD, windRange,
                     moistureMultiplier, moistureSD, moistureRange,
                     heightSD, heightRange, leafVar,updateProgress = TRUE) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Set input limits
  DFMCRange <- pmax(1.0001, DFMCRange)
  
  if (jitters > 2) {
    for (j in 1:jitters) {
      db.recreate <- j == 1
      # Modify environmental parameters
      s <- rtnorm(n = 1, mean = slope, sd = slopeSD,
                  a = slope-(slopeRange/2), b = slope+(slopeRange/2))
      t <- rtnorm(n = 1, mean = temp, sd = tempSD,
                  a = temp-(tempRange/2), b = temp+(tempRange/2)) 
      d <- rtnorm(n = 1, mean = DFMC, sd = DFMCSD,
                  a = pmax(0.02, DFMC-(DFMCRange/2)), b = pmin(0.199,DFMC+(DFMCRange/2))) 
      w <- rtnorm(n = 1, mean = wind, sd = windSD,
                  a = wind-(windRange/2), b = wind+(windRange/2))   
      
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", t) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
      
      # Modify vegetation
      tbl <- plantVar(base.params, Strata, Species, 
                      l = leafVar, Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange, 
                      Hs = heightSD, Hr = heightRange)
      ffm_run(tbl, db.path, db.recreate = db.recreate)
      Sys.sleep(0.25)
      ####UpdateProgress
      if (updateProgress == TRUE) {
        cat("Remaining replicates:", jitters - j, '\n')
      }
    }
  }  else  {
    print("Probabilistic analysis requires a minimum of 3 replicates")
  }
}


#' Models probabilistic fire behaviour using species-specific variability from FRaME tables
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param slope Mean slope (deg)
#' @param slopeSD Standard deviation of slope
#' @param slopeRange Truncates variability by +/- mean * range
#' @param temp Mean ambient temperature (deg.C)
#' @param tempSD Standard deviation of temp
#' @param tempRange Truncates variability by +/- mean * range
#' @param DFMC Mean DFMC (Percent ODW)
#' @param DFMCSD Standard deviation of DFMC
#' @param DFMCRange Truncates variability by +/- mean * range
#' @param wind Mean wind velocity (km/h)
#' @param windSD Standard deviation of wind velocity
#' @param windRange Truncates variability by +/- mean * range
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @param Flora A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export


probFire_Frame <- function(base.params, Structure, Flora, a, db.path = "out_mc.db",
                           slope, slopeSD, slopeRange, temp, tempSD, tempRange,
                           DFMC, DFMCSD, DFMCRange, wind, windSD, windRange, 
                           jitters = 50, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.001, 
                           updateProgress = NULL, testN = 10) {
  
  pbar <- txtProgressBar(max = jitters, style = 3)
  
  #Set input limits
  DFMCRange <- pmax(1.0001, DFMCRange)
  
  if (jitters > 2) {
    for (j in 1:jitters) {
      db.recreate <- j == 1
      # Modify environmental parameters
      s <- rtnorm(n = 1, mean = slope, sd = slopeSD,
                  a = slope-(slopeRange/2), b = slope+(slopeRange/2))
      t <- rtnorm(n = 1, mean = temp, sd = tempSD,
                  a = temp-(tempRange/2), b = temp+(tempRange/2)) 
      d <- rtnorm(n = 1, mean = DFMC, sd = DFMCSD,
                  a = pmax(0.02, DFMC-(DFMCRange/2)), b = pmin(0.199,DFMC+(DFMCRange/2))) 
      w <- rtnorm(n = 1, mean = wind, sd = windSD,
                  a = wind-(windRange/2), b = wind+(windRange/2))   
      
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", t) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
      
      # Select species for a random point in the forest and import the weather parameters
#      tbl <- specPoint(base.params, Structure, a) %>%
#        ffm_set_site_param("windSpeed", w, "km/h") %>%
#        ffm_set_site_param("temperature", t, "degc") %>%
#        ffm_set_site_param("deadFuelMoistureProp", d)
      
#      Strata <- strata(tbl)
#      Species <- species(tbl)
      
      # Recreate database on first run, then add following runs to this
#      db.recreate <- j == 1
      
#      TBL <- plantVarFrame(tbl, Strata, Species, Flora, a, l = l,
#                           Ms = Ms, Pm = Pm, Mr = Mr)
      
      # Choose random point and species, and vary plant traits for each species within their range
      TBL <- canopyCheck(base.params, testN = testN, Flora, Structure, a = a, l = l, Ms = Ms, Pm = Pm, Mr = Mr)
      # Run the model
      ffm_run(TBL, db.path = db.path, db.recreate = db.recreate)
      
      Sys.sleep(0.25)
      ####UpdateProgress
      if (is.function(updateProgress)) {
        text <- paste0("Number of remaining steps is ",max(weather$tm) - j )
        updateProgress(detail = text)
      }
      setTxtProgressBar(pbar, j)
    }
  }  else  {
    print("Probabilistic analysis requires a minimum of 3 replicates")
  }
}



#' Models fire behaviour across ranged variables
#' @param base.params Input parameter file
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windReps Number of wind speeds tested
#' @param windStep Gap (km/h) between wind steps
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param heightSD Standard deviation of plant height
#' @param heightRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param updateProgress Progress bar for use in the dashboard
#' @export
#' @examples
#' 
#' SPECIFY INPUTS
#' record <- 1
#' data(site)
#' data(structure)
#' data(flora)
#' data(traits)
#' base.params <- paramBuilder(site, structure, flora, traits, record)
#' 
#' MODEL PROBABILISTIC FIRE BEHAVIOUR
#' drivers(base.params, db.path = "out.db", jitters = REPLICATES, windMin = 0, windReps = 30,
#' windStep = 2, moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.001, 
#' heightSD = 0.2, heightRange = 1.41, leafVar = 0.1, updateProgress = NULL)
#' 
#' LOAD AND ORGANISE RESULTS
#' res<-ffm_db_load("out.db")

drivers <- function(base.params, db.path = "out_mc.db", jitters, windMin, windReps, windStep,
                    moistureMultiplier, moistureSD, moistureRange, heightSD, heightRange, 
                    leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  slopes <- c(-10,0,10,20,30)
  DFMCs <- c(0.05, 0.1, 0.15)
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  Niter <- nrow(dat) * jitters
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / jitters)
    db.recreate <- i == 1
    s <- dat[set, "slope"]
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / jitters)) {
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", 30) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #  for (j in 1:jitters) {
    #Randomise plant parameters
    tbl <- plantVar(base.params, Strata, Species, l = leafVar, Ms = moistureSD, 
                            Pm = moistureMultiplier, Mr = moistureRange, 
                            Hs = heightSD, Hr = heightRange)
    ffm_run(tbl, db.path, db.recreate = db.recreate)
    
    #  }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}




#' Models fire behaviour across ranged variables using species specific details
#' Discontinued function
#' 
#' @param base.params Input parameter file
#' @param a A unique identifier for the record being run
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windReps Number of wind speeds tested
#' @param windStep Gap (km/h) between wind steps
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard


driversS <- function(base.params, a, db.path = "out_mc.db", jitters, windMin, windReps, windStep,
                     moistureMultiplier, moistureSD, moistureRange, Variation,
                     leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  slopes <- c(-10,0,10,20,30)
  DFMCs <- c(0.05, 0.1, 0.15)
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  Niter <- nrow(dat) * jitters
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / jitters)
    db.recreate <- i == 1
    s <- dat[set, "slope"]
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / jitters)) {
      base.params <- base.params %>%
        ffm_set_site_param("slope", s, "deg") %>%
        ffm_set_site_param("temperature", 30) %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #  for (j in 1:jitters) {
    
    #Randomise plant parameters
    tbl <- plantVarS(base.params, Strata, Species, Variation, a, l = leafVar, 
                             Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange)
    ffm_run(tbl, db.path, db.recreate = db.recreate)
    
    #  }
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}




#' Models fire behaviour across ranged variables using species specific details
#' @param base.params Input parameter file
#' @param a A unique identifier for the record being run
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions for each row in the weather table
#' @param windMin Lowest wind velocity tested (km/h)
#' @param windReps Number of wind speeds tested
#' @param windStep Gap (km/h) between wind steps
#' @param slopes List of slope values for testing
#' @param DFMCs List of DFMC values for testing
#' @param temperature Standardised air temperature for the test
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param Variation A database of plant variability in traits, with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' Hs - Standard deviation of plant height variations
#' Hr - Truncates plant height variability by +/- Hr * height
#' @param updateProgress Progress bar for use in the dashboard
#' @export


driversRand <- function(base.params, a, db.path = "out_mc.db", replicates, windMin, windReps, windStep,
                        slopes, DFMCs, moistureMultiplier, moistureSD, moistureRange, temperature, Variation,
                        leafVar,updateProgress = NULL) {
  
  # Collect original descriptors
  Strata <- strata(base.params)
  Species <- species(base.params)
  
  #Range environmental values
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  if (is.null(slopes)){
    dat <- expand.grid(DFMC = DFMCs, wind = winds)
  } else {
    dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  }
  Niter <- nrow(dat) * replicates
  
  #Set test temperature
  base.params <- base.params %>%
    ffm_set_site_param("temperature", temperature, "degC")
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / replicates)
    db.recreate <- i == 1
    if (!is.null(slopes)){
      s <- dat[set, "slope"]}
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / replicates)) {
      if (!is.null(slopes)){
        base.params <- base.params %>%
          ffm_set_site_param("slope", s, "deg")}
      base.params <- base.params %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #Randomise plant parameters
    tbl <- plantVarS(base.params, Strata, Species, Variation, a, l = leafVar, 
                             Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange)
    ffm_run(tbl, db.path, db.recreate = db.recreate)
    Sys.sleep(0.25)
    
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}

#' Models fire behaviour across ranged variables using species specific details
#'
#' @param a A unique identifier for the record being run
#' @param db.path Name of the exported database
#' @param slopes List of slope values for testing
#' @param DFMCs List of DFMC values for testing
#' @param moistureMultiplier Multiplies all LFMC values by this number
#' @param moistureSD Standard deviation of moisture
#' @param moistureRange Truncates variability by +/- mean * range
#' @param leafVar Variation around input leaf dimensions, equivalent to l
#' @param Flora  A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param replicates 
#' @param testN 
#' @param updateProgress Progress bar for use in the dashboard
#' @param default.species.params 
#' @param tAir 
#' @param winds 
#'
#' @export

driversFrame <- function(Flora, Structure, default.species.params, a, db.path = "out_mc.db", replicates = 3, winds = seq(0, 60, 5),
                         slopes, DFMCs, tAir, moistureMultiplier, moistureSD, moistureRange, leafVar, testN = 10, updateProgress = NULL) {
  
  # Collect original descriptors
  base.params <- suppressWarnings(frame::buildParams(Flora, Structure, default.species.params, a,
                                                     fLine = 100, slope = 0, temp = 30, dfmc = 0.1, wind = 10))
  
  #Range environmental values
  winds <- seq(windMin, (windReps*windStep+windMin), windStep)
  
  #Dataframe of orthogonal combinations
  if (is.null(slopes)){
    dat <- expand.grid(DFMC = DFMCs, wind = winds)
  } else {
    dat <- expand.grid(slope = slopes, DFMC = DFMCs, wind = winds)
  }
  Niter <- nrow(dat) * replicates
  
  #Set test temperature
  base.params <- base.params %>%
    ffm_set_site_param("temperature", temperature, "degC")
  
  #Loop through combinations
  pbar <- txtProgressBar(max = Niter, style = 3)
  for (i in 1:Niter) {
    set <- ceiling(i / replicates)
    db.recreate <- i == 1
    if (!is.null(slopes)){
      s <- dat[set, "slope"]}
    d <- dat[set, "DFMC"]
    w <- dat[set, "wind"]
    
    #Update environmental parameters if on a new row
    if (set > ceiling((i-1) / replicates)) {
      if (!is.null(slopes)){
        base.params <- base.params %>%
          ffm_set_site_param("slope", s, "deg")}
      base.params <- base.params %>%
        ffm_set_site_param("deadFuelMoistureProp", d) %>%
        ffm_set_site_param("windSpeed", w)
    }
    
    #Randomise plant parameters
#    base.params <- frame::specPoint(base.params, Structure, a)
#    Strata <- strata(base.params)
#    Species <- species(base.params)
#    tbl <- plantVarFrame(base.params, Strata, Species, Flora, a, l = leafVar,
#                         Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange)
#    ffm_run_robust(base.params=tbl, db.path, db.recreate, testN,
#                   Flora = Flora, Structure = Structure, a = a, l = leafVar,
#                   Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange)
    
    # Choose random point and species, and vary plant traits for each species within their range
    TBL <- canopyCheck(base.params, testN = testN, Flora, Structure, a = a, l = l, Ms = Ms, Pm = Pm, Mr = Mr)
    # Run the model
    ffm_run(TBL, db.path = db.path, db.recreate = db.recreate)
    Sys.sleep(0.25)
    
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", Niter - i)
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
}


#' Updates parameter files with weather from a dataset,
#' then models fire from non-deterministic plant parameters
#' using plantVarFrame to modify individual species with their own measured variability
#' Differs from weatherSetS by using  FRaME formatted tables
#'
#' @param base.params Parameter input table
#' @param a A unique identifier for the record being run
#' @param weather A dataframe with the four fields:
#' tm - Sequential numbering of the records
#' T - Air temperature (deg C)
#' W - Wind velocity (km/h)
#' DFMC - Dead fuel moisture content (proportion ODW)
#' @param db.path Name of the exported database
#' @param jitters Number of repetitions
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' nsR, eR, mR, cR - maximum species richness recorded for each stratum
#' @param Flora  A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param updateProgress Progress bar for use in the dashboard
#'
#' @return dataframe
#' @export

weatherSet_Frame <- function(base.params, weather, Structure, Flora, a, db.path = "out_mc.db", jitters = 10, l = 0.1,
                             Ms = 0.01, Pm = 1, Mr = 1.001, updateProgress = NULL)
{
  pbar <- txtProgressBar(max = max(weather$tm), style = 3)
  
  # Read weather values from the table
  for (i in 1:max(weather$tm)) {
    cat("\n", " Modelling time-step number", i, "\n")
    w <- weather$W[[i]]
    t <- weather$T[[i]]
    d <- max(0.01,min(0.199,weather$DFMC[[i]]))
    
    # Select species for a random point in the forest and import the weather parameters
    tbl <- specPoint(base.params, Structure, a) %>%
      ffm_set_site_param("windSpeed", w, "km/h") %>%
      ffm_set_site_param("temperature", t, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", d)
    
    Strata <- strata(tbl)
    Species <- species(tbl)
    
    if (jitters > 0) {
      for (j in 1:jitters) {
        # Recreate database on first run, then add following runs to this
        db.recreate <- (i * j) == 1
        
        # Vary plant traits for each species within their range
        TBL <- plantVarFrame(tbl, Strata, Species, Flora, a, l = l,
                             Ms = Ms, Pm = Pm, Mr = Mr)
        # Run the model
        ffm_run(TBL, db.path, db.recreate = db.recreate)
      }
    }
    
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ",max(weather$tm) - i )
      updateProgress(detail = text)
    }
    setTxtProgressBar(pbar, i)
  }
  
  cat("Finished.  Output written to", db.path)
}


#' Modifies inputs until model runs
#' Discontinued
#'
#' @param base.params 
#' @param db.path 
#' @param db.recreate 
#' @param testN Number of replicates to allow
#' @param Flora  A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param a A unique identifier for the record being run
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#'
#' @return \code{TRUE} if the run completed and results were written to
#'   the output database successfully; \code{FALSE} otherwise.
#'

ffm_run_robust <- function(base.params, db.path, db.recreate = TRUE, testN = 5,
                           Flora, Structure, a = 1, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5){
  rep<-1
  while (rep < testN) {
    if (!ffm_run(base.params, db.path = db.path, db.recreate = db.recreate)) {
      base.params <- canopyCheck(base.params, testN = testN, Flora, Structure, a = a, l = l, Ms = Ms, Pm = Pm, Mr = Mr)
      rep <- rep+1
      if (rep >= testN) {
        print("Parameter file is faulty")
        test <- FALSE
      }
    } else {
      rep <- testN
      test <- TRUE
    }
  }
  return(test)
}


#' Modifies inputs until params file fits criteria
#'
#' @param base.params 
#' @param testN Number of replicates to allow
#' @param Flora  A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' hc, he, ht, hp & w - canopy dimensions for that species (m)
#' clump - mean ratio of clump diameter to crown diameter
#' openness - proportion of plant canopy occupied by gaps between clumps
#' @param Structure A dataframe with the fields:
#' record - a unique, consecutively numbered identifier per site
#' site - a unique identifier per site
#' NS, El, Mid & Can - the mean separation between plants (m) per stratum
#' ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
#' on the left grow directly beneath those in the stratum on the right. Acceptable values
#' are t, f, or blank, where the outcome will be decided by the relative stratum heights.
#' @param a A unique identifier for the record being run
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#'
#' @return Dataframe
#' @export
#'
canopyCheck <- function(base.params, testN = 5, Flora, Structure, a = 1, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5) {
  
  rep<-1
  while (rep < testN) {
    
    canopyStratum <- base.params$stratum[base.params$value == "canopy"][which(!is.na(base.params$stratum[base.params$value == "canopy"]))]
    
    if (mean(as.numeric(base.params$value[which((base.params$param == "hp" | base.params$param == "ht") & base.params$stratum == canopyStratum)])) < 1) {
      base.params <- frame::specPoint(base.params, Structure, a)
      Strata <- strata(base.params)
      Species <- species(base.params)
      base.params <- frame::plantVarFrame(base.params, Strata, Species, Flora, a, l,
                                          Ms = Ms, Pm = Pm, Mr = Mr)
      rep <- rep+1  
      if (rep >= testN) {
        stop("Vegetation is too low to model accurately")
      }
    } else {
      rep <- testN
    }
  }
  return(base.params)
}

