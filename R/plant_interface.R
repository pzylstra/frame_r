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
  datW <- dat %>%
    group_by(Stratum) %>%
    summarise_if(is.numeric, mean) %>%
    select(Stratum, w)
  strata <- dat %>%
    group_by(Stratum) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(di = sqrt(2/d),
           cr = sqrt(1/d),
           sep = (di+cr)/2) %>%
    select(Stratum, sep)
  strata <- left_join(strata, datW) %>%
    mutate(sep = pmax(sep,w)) %>%
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
#' @param sLitter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' 
#' @export
#'

buildFloraP <- function(comm, tr, age, rec = 1, moist = 1, sLitter = 15, diameter = 0.005) {
  
  comm <- left_join(as.data.frame(comm), tr, by = c("species" = "Species"))%>%
    mutate(species = name)
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
  flo$weight[n] <- sLitter
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
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param lwRat Ratio of leaf length to width
#' @param ram Stem ramification
#' @param deltaL Leaf density (g/cm3)
#'
#' @export
#'

buildTraitsP <- function(comm, propDead = 0, leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                         ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46) {
  
  # Summarise strata
  data <- comm %>%
    mutate(c_vol = (pi*(w/2)^2*(height-ht)/3)+(pi*(w/2)^2*(ht-he))+(pi*(w/2)^2*(he-base)),
           nLeaves = area_leaf/leafA,
           branchA = frontalArea * C.C_rat,
           branchV = 4/3*sqrt(branchA/pi)*branchA,
           nClumps = (c_vol/(1+G.C_rat))/branchV,
           leavesClump = nLeaves/nClumps,
           clumpD = 2*sqrt(branchA/pi),
           lSep = ((clumpD*ram)/((leavesClump/0.88)^(1/1.18)))/100)
  
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
    flo$leafThickness[sp]  <- lma[sp]/(deltaL*1000) # For LMA kg/m2
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
#' @param tr Table of input traits
#' @param deltaL Default leaf density (g/cm3)
#' @param propDead Default proportion of foliage dead
#' @param leafForm Default leaf shape: Flat or Round
#' @param lwRat Default ratio of leaf length to width
#' @param leafA Default area of a leaf in m2
#' @param ignitionTemp Default temperature of the endotherm (degC)
#' @param moist Default leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Default ratio of gaps to clumps of leaves
#' @param C.C_rat Default ratio of clump to canopy size
#' @param ram Default branch ramification (stem order)
#' @param bark_density Default bark density (kg)
#' 
#' @export
#'

collectTraitsP <- function(comm, tr, 
                           leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5, propDead = 0,
                           ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46,
                           bark_density = 200) {
  comm <- left_join(as.data.frame(comm), tr, by = c("species" = "Species"))%>%
    mutate(species = name)
  
  # Join to plant traits and clean empty data
  comm["leafForm"][is.na(comm["leafForm"])] <- leafForm
  comm["moisture"][is.na(comm["moisture"])] <- moist
  comm["propDead"][is.na(comm["propDead"])] <- propDead
  comm["ignitionTemp"][is.na(comm["ignitionTemp"])] <- ignitionTemp
  comm["lwRat"][is.na(comm["lwRat"])] <- lwRat
  comm["leaf_area"][is.na(comm["leaf_area"])] <- leafA
  comm["bark_density"][is.na(comm["bark_density"])] <- bark_density
  comm["G.C_rat"][is.na(comm["G.C_rat"])] <- G.C_rat
  comm["C.C_rat"][is.na(comm["C.C_rat"])] <- C.C_rat
  comm["stemOrder"][is.na(comm["stemOrder"])] <- ram
  comm[["leaf_thickness"]][is.na(comm[["leaf_thickness"]])] <- comm$lma/(deltaL*1000)
  comm[["name"]][is.na(comm[["name"]])] <- comm$species
  
  comm <- comm %>%
    mutate(c_vol = (pi*(w/2)^2*(height-ht)/3)+(pi*(w/2)^2*(ht-he))+(pi*(w/2)^2*(he-base)),
           nLeaves = area_leaf/leaf_area,
           branchA = frontalArea * C.C_rat,
           branchV = 4/3*sqrt(branchA/pi)*branchA,
           nClumps = (c_vol/(1+G.C_rat))/branchV,
           leavesClump = max(nLeaves/nClumps,1),
           clumpD = 2*sqrt(branchA/pi),
           lSep = ((clumpD*stemOrder)/((leavesClump/0.88)^(1/1.18)))/100)
  
  # Summarise strata
  sp_names <- unique(comm$name)
  
  form <- comm %>%
    select(species, leafForm)
  form <- form[!duplicated(form), ]
  
  summ <- comm %>%
    group_by(species) %>%
    summarise_if(is.numeric, mean)
  spec <- unique(summ$species)
  
  # Create trait table
  t <- data.frame(matrix(ncol = 12, nrow = length(spec)))
  colnames(t) <- c("name", "propDead", "leafForm", "leafThickness", "leafWidth", "leafLength", 
                   "leafSeparation", "stemOrder", "ignitionTemp", "moisture", "G.C_rat", "C.C_rat")
  n <- 1
  for (sp in 1: length(spec)) {
    t$name[sp] <- sp_names[sp]
    t$propDead[sp]  <- summ$propDead[sp]
    t$leafForm[sp]  <- form$leafForm[sp]
    t$leafThickness[sp]  <- summ$leaf_thickness[sp]
    t$leafLength[sp]  <- sqrt(2*summ$leaf_area[sp]*summ$lwRat[sp])
    t$leafWidth[sp]  <- t$leafLength[sp] / summ$lwRat[sp]
    t$leafSeparation[sp]  <- summ$lSep[sp]
    t$stemOrder[sp]  <- summ$stemOrder[sp]
    t$ignitionTemp[sp]  <- summ$ignitionTemp[sp]
    t$moisture[sp]  <- summ$moisture[sp]
    t$G.C_rat[sp] <- summ$G.C_rat[sp]
    t$C.C_rat[sp] <- summ$C.C_rat[sp]
    n <- n+1
  }
  return(t)
}


#' Constructs parameter files for FRaME using 
#' outputs from plant
#'
#' @param dat The results of run_scm_collect
#' @param tr Table of input traits
#' @param age Years since disturbance
#' @param rec Number of the record
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param deltaL Leaf density (g/cm3)
#' @param sLitter Weight of surface litter (t/ha)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param lwRat Ratio of leaf length to width
#' @param ram Stem ramification
#' @param hw Environmental difference in plant height: width ratio
#'
#' @export
#'

frameTables <- function(dat, tr, age, propSamp = 0.75, transects = 10, sepSig = 0.1, rec = 1, propDead = 0, 
                        leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                        ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, 
                        deltaL = 0.46, hw = 0, sLitter = 15, diameter = 0.005, minCov = 0.0001) {
  
  comm <- suppressMessages(stratify_community(dat, tr, age, hw, propSamp = propSamp, transects = transects, sepSig = sepSig, minCov = minCov))
  Structure <- suppressMessages(buildStructureP(comm, age, rec))
  Flora <- buildFloraP(comm, tr, age, rec, moist, sLitter, diameter)
  Traits <- collectTraitsP(comm, tr)
  
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
#' @param interval List of sampling intervals, length 3.
#' @param breaks List of time periods for each sampling interval, length 3
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param propDead Proportion of foliage dead
#' @param leafForm Flat or Round
#' @param leafA Area of a leaf in m2
#' @param ignitionTemp Temperature of the endotherm (degC)
#' @param moist Leaf moisture (ratio moisture weight to dry weight)
#' @param G.C_rat Ratio of gaps to clumps of leaves
#' @param C.C_rat Ratio of clump to canopy size
#' @param deltaL Leaf density (g/cm3)
#' @param diameter Mean diameter of surface litter pieces (m)
#' @param lwRat Ratio of leaf length to width
#' @param ram Stem ramification
#' @param hw Environmental difference in plant height: width ratio
#' @param mat Mean annual temperature (deg C)
#' @param propSamp Values closer to 0 have more accurate ratios of components but miss some cohorts
#' @param transects More transects ensure more cohorts 
#' @param sepSig 
#'
#' @export

frameDynTab <- function(dat, tr, breaks = c(20,50,200), interval = c(2,5,10), propDead = 0, leafForm = "Flat", lwRat = 3, leafA = 0.002547, ram = 5,
                        ignitionTemp = 260, moist = 1, G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46, hw = 0, mat = 17, diameter = 0.005, propSamp = 0.75, 
                        transects = 10, sepSig = 0.1, minCov = 0.0001) {
  
  cat("Collecting FRaME parameters", "\n")
  
  Flora <- data.frame()
  Structure <- data.frame()
  Tr <- data.frame()
  stepsA <- seq(interval[1],breaks[1], by = interval[1])
  stepsB <- seq(breaks[1]+interval[2],breaks[2], by=interval[2])
  stepsC <- seq(breaks[2]+interval[3],breaks[3], by=interval[3])
  steps <- append(stepsA,stepsB)
  steps <- append(steps,stepsC)
  
  # Uses partial plant outputs for litter
  result <-  dat%>% 
    plant::tidy_patch() %>% 
    FF16_expand_state() 
  tab <- result$species%>%
    drop_na()
  mHt <- max(tab$height)
  max <- 7.9*log(mHt)-12.64
  rate <- 0.12775*exp(0.109*mat)
  rec <- 1
  for (age in steps) {
    sLitter <- frame::litter(negEx = 1, max, rate, a = 1, b = 1, age)
    tabs <- frameTables(dat, tr, age, propSamp, transects, sepSig, rec, propDead, leafForm, lwRat, leafA, ram,
                        ignitionTemp, moist, G.C_rat, C.C_rat, deltaL, hw, sLitter, diameter, minCov = minCov)
    Flora <- rbind(Flora,tabs[[1]])
    Structure <- rbind(Structure,tabs[[2]])
    Tr <- rbind(Tr,tabs[[3]])
    rec <- rec+1
  }
  Traits <- Tr %>%
    group_by(name) %>%
    summarise_if(is.numeric, mean) %>%
    mutate(leafForm = "Flat")
  
  return(list(Flora, Structure, Traits)) 
}



#' Models fire behaviour across a range of ages
#' from plant modelling
#'
#' @param db.path 
#' @param reps 
#' @param slope 
#' @param slopeSD 
#' @param slopeRange 
#' @param temp 
#' @param tempSD 
#' @param tempRange 
#' @param DFMC 
#' @param DFMCSD 
#' @param DFMCRange 
#' @param wind 
#' @param windSD 
#' @param windRange 
#' @param moistureMultiplier 
#' @param moistureSD 
#' @param moistureRange 
#' @param fLine Length of the active fire front (m)
#' @param leafVar 
#' @param updateProgress 
#' @param dat The results of frameDynTab
#'
#' @export

firePlant <- function(dat, db.path = "out.plant.db", reps = 5,
                      slope = 0, slopeSD = 2, slopeRange = 5, 
                      temp = 30, tempSD = 5, tempRange = 3,
                      DFMC = 0.1, DFMCSD = 0.01, DFMCRange = 2, 
                      wind = 5, windSD = 1, windRange = 2, vAir500 = 2,
                      moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.5,
                      fLine = 1000, leafVar = 0.1, updateProgress = NULL) {
  
  Flora <- dat[[1]]
  Structure <- dat[[2]]
  default.species.params <- dat[[3]]
  r <- unique(Flora$record)
  RUNS <-data.frame()
  PATHS <- data.frame()
  Ages <- unique(Flora$site)
  
  for (n in r) {
    cat("Modelling age", Ages[n], "\n")
    f <- filter(Flora, record == n)
    s <- filter(Structure, record == n)
    base.params <- suppressWarnings(frame::buildParams(Structure = s, Flora = f, default.species.params, a = n,
                                                       fLine = fLine, slope = slope, temp = temp, dfmc = DFMC, wind = wind))
    
    Strata <- strata(base.params)
    
    # Find foliage projective cover
    for (st in 1:nrow(Strata)) {
      if (st == 1) {
        FPC <- Strata$cover[st]
      } else {
        FPC <- FPC + ((1-FPC) * Strata$cover[st])
      }
    }
    
    probFire_Frame(base.params, Structure = s, Flora = f, a = n, db.path = db.path,
                   slope = slope, slopeSD = slopeSD, slopeRange = slopeRange, 
                   temp = temp, tempSD = tempSD, tempRange = tempRange,
                   DFMC = DFMC, DFMCSD = DFMCSD, DFMCRange = DFMCRange, 
                   wind = wind, windSD = windSD, windRange = windRange,
                   jitters = reps, l = leafVar, 
                   Ms = moistureSD, Pm = moistureMultiplier, Mr = moistureRange, 
                   updateProgress = updateProgress)
    
    #SUMMARISE BEHAVIOUR
    res<-ffm_db_load(db.path)
    runs <- suppressMessages(frame::frameSummary(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults)%>%
                               mutate(site = f$site[1],
                                      FPC = FPC,
                                      Spotting = spotFire(flameHeight = fh, slope = slope_degrees, FPC, windExposure = 1, vAir = wind_kph, vAir500, fireArea = ros_kph*(fLine/1000)),
                                      fReach = max(lengthPlant * cos(flameAngle), lengthSurface * cos(angleSurface), Spotting)))
    IP <- frame::repFlame(res$IgnitionPaths) %>%
      mutate(site = f$site[1])
    scorch <- suppressMessages(frame::flora(runs, IP, Param = base.params, Test = 80)) %>%
      select(!wind_kph)
    outa <- suppressMessages(left_join(runs,scorch, by = "repId") )
    RUNS <- rbind(RUNS, outa)
    PATHS <- rbind(PATHS, IP)
    
  }
  return(list(RUNS, PATHS))
}