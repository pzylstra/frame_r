#' Creates an 'F_structure' table for FRaME from plant modelling
#'
#' @param dat The output from stratify_community
#' @param age Years since disturbance
#' @param rec Number of the record
#' 
#' @export
#'

buildStructureP <- function(dat, age, rec = 1) {
  # Summarise strata
  strata <- dat %>%
    group_by(Stratum) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(di = sqrt(2/d),
           cr = sqrt(1/d),
           sep = (di+cr)/2) %>%
    select(Stratum, sep)
  
  # Create structure table
  struct <- data.frame(matrix(ncol = 15, nrow = 1))
  colnames(struct) <- c("record", "site", "NS", "El", "Mid", "Can", "ns_e", "ns_m", "e_m", "e_c", "m_c", "nsR", "eR", "mR", "cR")
  struct$record <- rec
  struct$site <- age
  struct$NS <- signif(strata$sep[1],digits = 3)
  ns <- filter(dat, dat$Stratum == 1)
  struct$nsR <- length(unique(ns$species))
  if (nrow(strata) == 4) {
    E <- filter(dat, dat$Stratum == 2)
    struct$eR <- length(unique(E$species))
    struct$El <- signif(strata$sep[2],digits = 3)
    M <- filter(dat, dat$Stratum == 3)
    struct$mR <- length(unique(M$species))
    struct$Mid <- signif(strata$sep[3],digits = 3)
    C <- filter(dat, dat$Stratum == 4)
    struct$cR <- length(unique(C$species))
    struct$Can <- signif(strata$sep[4],digits = 3)
    
  } else if (nrow(strata) == 3) {
    E <- filter(dat, dat$Stratum == 2)
    struct$eR <- length(unique(E$species))
    struct$El <- signif(strata$sep[2],digits = 3)
    C <- filter(dat, dat$Stratum == 3)
    struct$cR <- length(unique(C$species))
    struct$Can <- signif(strata$sep[3],digits = 3)
  } else {
    C <- filter(dat, dat$Stratum == 2)
    struct$cR <- length(unique(C$species))
    struct$Can <- signif(strata$sep[2],digits = 3)
  }
  return(struct)
}


#' Creates an 'F_flora' table for FRaME from plant modelling
#'
#' @param comm The output from stratify_community
#' @param tr An optional table of input traits
#' @param age Years since disturbance
#' @param rec Number of the record
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param litter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' 
#' @export
#'

buildFloraP <- function(comm, tr, age, rec = 1, moist = 1, litter = 15, diameter = 0.005) {
  
  # Summarise strata
  strata <- comm %>%
    group_by(Stratum, species) %>%
    summarise_if(is.numeric, mean) %>%
    select(species, Stratum, height, base, he, ht, w)
  sdStrat <- comm %>%
    group_by(Stratum, species) %>%
    summarise_if(is.numeric, sd) %>%
    select(species, Stratum, height)
  MStrat <- comm %>%
    group_by(Stratum, species) %>%
    summarise_if(is.numeric, max) %>%
    select(species, Stratum, height)
  mStrat <- comm %>%
    group_by(Stratum, species) %>%
    summarise_if(is.numeric, min) %>%
    select(species, Stratum, height)
  
  # Create flora table
  flo <- data.frame(matrix(ncol = 14, nrow = nrow(strata)+1))
  colnames(flo) <- c("record", "site", "species", "stratum", "comp", "base", "he", "ht", "top", "w", "Hs", "Hr", "weight", "diameter")
  flo$record <- rec
  flo$site <- age
  n <- 1
  for (st in as.numeric(unique(comm$Stratum))) {
    setSt <- comm[comm$Stratum == st,]
    for (sp in unique(setSt$species)) {
      set <- setSt[setSt$species == sp,]
      flo$species[n] <- as.character(sp)
      flo$stratum[n] <- as.character(st)
      flo$comp[n] <- max(sum(set$d),0.00001)
      flo$base[n] <- max(strata$base[n], 0)
      flo$he[n] <- max(strata$he[n], 0)
      flo$ht[n] <- max(strata$ht[n], 0)
      flo$top[n] <- max(strata$height[n], 0)
      flo$w[n] <- max(strata$w[n], 0.00001)
      flo$Hs[n] <- if(is.na(sdStrat$height[n])) {
        1.00001
      } else {
        sdStrat$height[n]
      }
      flo$Hr[n] <- if(is.na(sdStrat$height[n])) {
        1.00001
      } else {
        MStrat$height[n] - mStrat$height[n]
      }
      n <- n+1
    }
  }
  flo$species[n] <- "Litter"
  flo$weight[n] <- litter
  flo$diameter[n] <- diameter
  return(flo)
}


#' Creates traits table for FRaME from plant modelling
#' 
#' Leaf dimensions default to ausTraits mean leaf size
#'
#' @param comm The output from stratify_community
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param lwRatio Ratio of leaf length to width
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param deltaL Leaf density (g/cm3)
#' @param lma List of LMA values per species (kgm−2)
#' 
#' @export
#'

buildTraitsP <- function(comm, propDead = 0, leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                         ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46, lma) {
  
  # Summarise strata
  data <- comm %>%
    mutate(c_vol = (pi*(w/2)^2*(height-ht)/3)+(pi*(w/2)^2*(ht-he))+(pi*(w/2)^2*(he-base)),
           nLeaves = area_leaf/leafA,
           branchA = frontalArea * C.C_rat,
           branchV = 4/3*sqrt(branchA/pi)*branchA,
           nClumps = (c_vol/(1+G.C_rat))/branchV,
           leavesClump = nLeaves/nClumps,
           clumpD = sqrt(branchA/pi),
           lSep = (branchV*ram)/((leavesClump/0.88)^(1/1.18)))
  
  summ <- data %>%
    group_by(species) %>%
    summarise_if(is.numeric, mean) %>%
    select(species, lSep)
  spec <- as.numeric(unique(summ$species))
  
  # Create flora table
  flo <- data.frame(matrix(ncol = 12, nrow = length(spec)))
  colnames(flo) <- c("name", "propDead", "leafForm", "leafThickness", "leafWidth", "leafLength", 
                     "leafSeparation", "stemOrder", "ignitionTemp", "moisture", "G.C_rat", "C.C_rat")
  n <- 1
  for (sp in spec) {
    flo$name[sp] <- as.character(sp)
    flo$propDead[sp]  <- 0
    flo$leafForm[sp]  <- "Flat"
    flo$leafThickness[sp]  <- lma[sp]/(deltaL*1000)
    flo$leafLength[sp]  <- sqrt(2*leafA*lwRat)
    flo$leafWidth[sp]  <- flo$leafLength[sp] / lwRat
    flo$leafSeparation[sp]  <- summ$lSep[sp]
    flo$stemOrder[sp]  <- ram
    flo$ignitionTemp[sp]  <- ignitionTemp
    flo$moisture[sp]  <- moist
    flo$G.C_rat[sp] <- G.C_rat
    flo$C.C_rat[sp] <- C.C_rat
    n <- n+1
  }
  return(flo)
}


#' Creates traits table for FRaME from plant modelling
#' 
#' Reads traits from an input table
#'
#' @param comm The output from stratify_community
#' @param tr A table of input traits
#' @param deltaL Leaf density (g/cm3)
#' 
#' @export
#'

collectTraitsP <- function(comm, deltaL = 0.46) {
  
  # Join to plant traits and clean empty data
  comm["leafForm"][is.na(comm["leafForm"])] <- "Flat"
  comm["moisture"][is.na(comm["moisture"])] <- 1
  comm["propDead"][is.na(comm["propDead"])] <- 0
  comm["ignitionTemp"][is.na(comm["ignitionTemp"])] <- 260
  comm["lwRat"][is.na(comm["lwRat"])] <- 3
  comm["leaf_area"][is.na(comm["leaf_area"])] <- 0.002547
  comm["bark_density"][is.na(comm["bark_density"])] <- 200
  comm["G.C_rat"][is.na(comm["G.C_rat"])] <- 3
  comm["C.C_rat"][is.na(comm["C.C_rat"])] <- 0.1
  comm["stemOrder"][is.na(comm["stemOrder"])] <- 5
  
  # Summarise strata
  comm <- comm %>%
    mutate(c_vol = (pi*(w/2)^2*(height-ht)/3)+(pi*(w/2)^2*(ht-he))+(pi*(w/2)^2*(he-base)),
           leaf_thickness = case_when(is.na(leaf_thickness) ~ lma/(deltaL*1000),
                                      TRUE ~ leaf_thickness),
           name = case_when(is.na(name) ~ species,
                            TRUE ~ name),
           nLeaves = area_leaf/leaf_area,
           branchA = frontalArea * C.C_rat,
           branchV = 4/3*sqrt(branchA/pi)*branchA,
           nClumps = (c_vol/(1+G.C_rat))/branchV,
           leavesClump = nLeaves/nClumps,
           clumpD = sqrt(branchA/pi),
           lSep = (branchV*stemOrder)/((leavesClump/0.88)^(1/1.18)))
  
  sp_names <- unique(comm$name)
  
  form <- comm %>%
    select(species, leafForm)
  form <- form[!duplicated(form), ]
  
  summ <- comm %>%
    group_by(species) %>%
    summarise_if(is.numeric, mean)
  spec <- unique(summ$species)
  
  # Create flora table
  flo <- data.frame(matrix(ncol = 12, nrow = length(spec)))
  colnames(flo) <- c("name", "propDead", "leafForm", "leafThickness", "leafWidth", "leafLength", 
                     "leafSeparation", "stemOrder", "ignitionTemp", "moisture", "G.C_rat", "C.C_rat")
  n <- 1
  for (sp in 1: length(spec)) {
    flo$name[sp] <- sp_names[sp]
    flo$propDead[sp]  <- summ$propDead[sp]
    flo$leafForm[sp]  <- form$leafForm[sp]
    flo$leafThickness[sp]  <- summ$leaf_thickness[sp]
    flo$leafLength[sp]  <- sqrt(2*summ$leaf_area[sp]*summ$lwRat[sp])
    flo$leafWidth[sp]  <- flo$leafLength[sp] / summ$lwRat[sp]
    flo$leafSeparation[sp]  <- summ$lSep[sp]
    flo$stemOrder[sp]  <- summ$stemOrder[sp]
    flo$ignitionTemp[sp]  <- summ$ignitionTemp[sp]
    flo$moisture[sp]  <- summ$moisture[sp]
    flo$G.C_rat[sp] <- summ$G.C_rat[sp]
    flo$C.C_rat[sp] <- summ$C.C_rat[sp]
    n <- n+1
  }
  return(flo)
}


#' Constructs parameter files for FRaME using 
#' outputs from plant
#'
#' @param dat The results of run_scm_collect
#' @param tr An optional table of input traits
#' @param age Years since disturbance
#' @param rec Number of the record
#' @param sample Proportion of cohorts to test (0-1)
#' @param transects Number of repeats for each sample
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param openness Ratio of gaps to clumps of leaves
#' @param clump Ratio of clump frontal area diameter to crown frontal area diameter
#' @param litter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param lwRatio Ratio of leaf length to width
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' @param deltaL Leaf density (g/cm3)
#' @param lma List of LMA values per species (kgm−2)
#'
#' @export
#'

frameTables <- function(dat, tr, age, rec = 1, sample = 0.5, transects = 10, propDead = 0, 
                        leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                        ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, 
                        deltaL = 0.46, lat = -35, map = 1000, mat = 20, lma) {
  
  comm <- stratify_community(dat, tr, age, lat, map, mat, sample, transects)
  Structure <- buildStructureP(comm, age, rec)
  comm <- frame:::updateSpecies(comm, tr)
  Flora <- buildFloraP(comm, age, rec, moist)
  Traits <- if(!missing(tr)) {
    collectTraitsP(comm, deltaL)
  } else {
    buildTraitsP(comm, propDead, leafForm, lwRat, leafA, ram,
                 ignitionTemp, moist, G.C_rat, C.C_rat, deltaL, lma)
  }
  
  return(list(Flora, Structure, Traits))
}

#' Updates species names from an optional table
#'
#' @param comm The output from stratify_community
#' @param tr An optional table of input traits
#'  

updateSpecies <- function(comm, tr){
  
  # Update species names if tr table used
  if(!missing(tr)) {
    comm <- left_join(comm,tr, by = c("species" = "Species")) %>%
      mutate(species = name)
  } 
}


#' Constructs an age sequence of parameter files 
#' built from plant modelling
#' 
#'
#' @param dat The results of run_scm_collect
#' @param tr An optional table of input traits
#' @param max Maximum years since disturbance
#' @param interval Time interval for sampling
#' @param rec Number of the record
#' @param sample Proportion of cohorts to test (0-1)
#' @param transects Number of repeats for each sample
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param openness Ratio of gaps to clumps of leaves
#' @param clump Ratio of clump frontal area diameter to crown frontal area diameter
#' @param litter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param lwRatio Ratio of leaf length to width
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' @param deltaL Leaf density (g/cm3)
#' @param lma List of LMA values per species (kgm−2)
#'
#' @export

frameDynTab <- function(dat, tr, upper, interval, sample = 0.5, transects = 10, propDead = 0, leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                        ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46, lat = -35, map = 1000, mat = 20, lma) {
  Flora <- data.frame()
  Structure <- data.frame()
  T <- data.frame()
  rec <- 1
  for (age in seq(interval,upper,by = interval)) {
    tabs <- frameTables(dat, tr, age, rec, sample, transects, propDead, leafForm, lwRat, leafA, ram,
                        ignitionTemp, moist, G.C_rat, C.C_rat, deltaL, lat, map, mat, lma)
    Flora <- rbind(Flora,tabs[[1]])
    Structure <- rbind(Structure,tabs[[2]])
    T <- rbind(T,tabs[[3]])
    rec <- rec+1
  }
  Traits <- T %>%
    group_by(name) %>%
    summarise_if(is.numeric, mean) %>%
    mutate(leafForm = "Flat")
  
  return(list(Flora, Structure, Traits)) 
}



#' Models fire behaviour across a range of ages
#' from plant modelling
#'
#' @param dat The results of frameDynTab
#' @param max Maximum years since disturbance
#' @param interval Time interval for sampling
#' @param rec Number of the record
#' @param sample Proportion of cohorts to test (0-1)
#' @param transects Number of repeats for each sample
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param openness Ratio of gaps to clumps of leaves
#' @param clump Ratio of clump frontal area diameter to crown frontal area diameter
#' @param litter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param lwRatio Ratio of leaf length to width
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' @param deltaL Leaf density (g/cm3)
#' @param lma List of LMA values per species (kgm−2)
#'
#' @export


firePlant <- function(dat, db.path = "out.plant.db", reps = 5,
                      slope = 0, slopeSD = 2, slopeRange = 5, 
                      temp = 30, tempSD = 5, tempRange = 3,
                      DFMC = 0.1, DFMCSD = 0.01, DFMCRange = 2, 
                      wind = 5, windSD = 1, windRange = 2,
                      moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.5,
                      fLine = 100, leafVar = 0.1, updateProgress = NULL) {
  
  Flora <- dat[[1]]
  Structure <- dat[[2]]
  default.species.params <- dat[[3]]
  r <- unique(Flora$record)
  out <-data.frame()
  
  for (n in r) {
    print(n)
    f <- filter(Flora, record == n)
    s <- filter(Structure, record == n)
    base.params <- suppressWarnings(frame::buildParams(Structure = s, Flora = f, default.species.params, a = n,
                                                       fLine = fLine, slope = slope, temp = temp, dfmc = DFMC, wind = wind))
    heightSD <- mean(f$Hs, na.rm = TRUE)
    heightRange <- mean(f$Hr, na.rm = TRUE)
    
    #MODEL PROBABILISTIC FIRE BEHAVIOUR
#    probFire(base.params, db.path = db.path, jitters = reps,
#             slope = slope, slopeSD = slopeSD, slopeRange = slopeRange, 
#             temp = temp, tempSD = tempSD, tempRange = tempRange,
#             DFMC = DFMC, DFMCSD = DFMCSD, DFMCRange = DFMCRange, 
#             wind = wind, windSD = windSD, windRange = windRange,
#             moistureMultiplier = moistureMultiplier, moistureSD = moistureSD, moistureRange = moistureRange,
#             heightSD = heightSD, heightRange = heightRange, 
#             leafVar = leafVar, updateProgress = TRUE)
    
    probFire_Frame(base.params, Structure = s, Flora = f, a = n, db.path = db.path,
                   slope = slope, slopeSD = slopeSD, slopeRange = slopeRange, 
                   temp = temp, tempSD = tempSD, tempRange = tempRange,
                   DFMC = DFMC, DFMCSD = DFMCSD, DFMCRange = DFMCRange, 
                   wind = wind, windSD = windSD, windRange = windRange,
                   jitters = reps, l = leafVar, 
                   Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange, 
                   updateProgress = updateProgress)
    
    #SUMMARISE BEHAVIOUR
    res<-ffm_db_load("out.plant.db")
    runs <- frame::frameSummary(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults)
    IP <- frame::repFlame(res$IgnitionPaths)
    scorch <- frame::flora(runs, IP, Param = base.params, Test = 80)
    outa <- left_join(runs,scorch, by = "repId") %>%
      mutate(step = n)
    out <- rbind(out, outa)
    
  }
  return(out)
}