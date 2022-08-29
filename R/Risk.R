#' Finds randomised lightning ignition time and location
#'
#' @param weather Hourly weather dataset from function frameWeather
#' @param smoulder Time for which an ignition may smoulder (hours)
#' @param Extinction Extinction moisture threshold
#'
#' @return dataframe
#' @export
#'
 
ignitions <- function(weather, smoulder = 24, Extinction = 0.15){
  
  # Location of strike in landscape
  landscapeLoc <- extraDistr::rtnorm(n = 1, mean = 0, sd = 0.3, a = -1, b = 1)
  landform <- case_when(landscapeLoc >= 0.95 ~ "moistureD",
                        landscapeLoc >= 0.68 ~ "moistureC",
                        landscapeLoc >= 0 ~ "moistureB",
                        TRUE ~ "moistureA")
  
  # Add random likelihood to lightning, then find the time of the strike
  weather$r <- runif(nrow(weather))
  weather$probStrike <- weather$cgStrikes*weather$r
  strike <- weather %>% arrange(desc(probStrike))
  strikeTime <- strike$Hour[1]
  
  # Determine whether ignition occurs either immediately or within smoulder time
  if (weather[strikeTime,landform]<= Extinction) {
    ignition <- TRUE
    ignitionTime <- strikeTime
  } else {
    weatherB <- weather %>%
      mutate(Hour = Hour + max(Hour))
    weatherLong <- rbind(weather,weatherB)
    moistSet <- (weatherLong[c(strikeTime:(strikeTime+smoulder-1)),landform])
    ignitionTime <- which(moistSet<0.15)[1] + strikeTime
    if (!is.na(ignitionTime)) {
      ignition <- TRUE
    } else {ignition <- FALSE}
  }
  
  out <- data.frame(strikeTime = strikeTime, ignition = ignition, ignitionTime = ignitionTime, landscapeLoc = landscapeLoc)
  return(out)
}

#' Support function to randomise river width
#' 
#' @param lRiver Likelihood of encountering a large river
#' @param mRiver Likelihood of encountering a medium river

river <- function(lRiver, mRiver){
  test <- runif(1)
  case_when(test <= lRiver ~ 150,
            test <= mRiver ~ 50,
            TRUE ~ 5
  )
}

#' Finds distance of fire spread for a part of a landscape
#'
#' @param base.params Parameter file
#' @param weather Formatted weather dataset
#' @param t Record number in the weather set for which tio model
#' @param hourStep Number of hours for which to model
#' @param Extinction Extinction litter moisture content (Percent ODW)
#' @param fReach Initial flame reach (m)
#' @param fEdge Fire edge: Head = 1, Flank = 0, Tail = -1
#' @param wHeading Wind direction in relation to slope exposure: Up sun slope = 1, cross slope = 0, down sun slope = -1
#' @param slopeM Mean slope
#' @param slopeSD Slope standard deviation
#' @param slopeLength Distance from ridgeline to river
#' @param mRiver Likelihood of encountering a medium river
#' @param lRiver Likelihood of encountering a large river
#' @param igLoc Ignition location between -1 (base of sunny slope), 0 (ridgeline) & 1 (base of shade slope)
#' @param firelineW 
#' @param Test 
#' @param sensitive 
#' @param girdleH 
#' @param woodDensity 
#' @param barkDensity 
#' @param bark 
#' @param comBark 
#' @param resBark 
#' @param phloem 
#' @param moisture 
#' @param bMoisture 
#' @param distance 
#' @param trail 
#' @param var 
#' @param Altitude 
#' @param necT 
#' @param surfDecl 
#' @param minR 
#' @param a 
#' @param fireArea 
#' @param FPC 
#' @param vAir500 
#' @param targSp 
#' @param testN 
#' @param Strata 
#' @param Species 
#' @param Flora 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param Structure 
#' @param strikeTime 
#'
#' @return dataframe
#' @export
#'

frameSpread <- function (base.params, weather, a, igLoc = -1, t = 1, hourStep = 3, strikeTime = 1,
                         Extinction = 0.15, fReach = 0, fEdge = 1, fireArea = 100, 
                         wHeading = 0, firelineW = 1, slopeM = 6.7, slopeSD = 4.9, 
                         slopeLength = 391, mRiver = 0.38, lRiver = 0.06, Test = 80, 
                         sensitive = 1, girdleH, woodDensity = 1000, barkDensity = 300, 
                         bark = 0.04, comBark = 700, resBark = 0, phloem = 0.01, 
                         moisture = 0.6, bMoisture = 0.5, FPC = 0.5, distance = 2, 
                         trail = 1000, var = 10, Altitude = 200, vAir500 = 2, necT = 60, 
                         surfDecl = 10, minR = 1, targSp = "c", testN = 5,
                         Strata, Species, Flora, Structure, l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5) 
{
  aTxt <- case_when(fEdge == 1 ~ "Head", 
                    fEdge == -1 ~ "Tail", 
                    TRUE ~ "Flank")
  dbH <- paste(aTxt, hourStep, "_", a, "_", round(runif(1, min = 1, max = 1e+07), 3), ".db", sep = "")
  Reach <- vector()
  Reach[1] <- fReach
  travDist <- vector()
  travDist[1] <- 0
  landscapeLoc <- vector()
  landStep <- 2
  sMort <- vector()
  sMort[1] <- 0
  gMort <- vector()
  gMort[1] <- NA
  R <- vector()
  R[1] <- 0
  aspList <- vector()
  stemTests <- data.frame()
  targSp <- base.params$species[max(which(base.params$param == "name" & base.params$value == targSp))]
  remTime <- hourStep
  landscapeLoc[landStep - 1] <- igLoc
  while (remTime > 0) {
    recreate <- remTime == hourStep
    landformDistance <- dplyr::case_when(landscapeLoc[landStep - 1] == 1 ~ frame:::river(lRiver, mRiver), 
                                         landscapeLoc[landStep - 1] >= 0.95 ~ slopeLength - (landscapeLoc[landStep - 1] * slopeLength), 
                                         landscapeLoc[landStep - 1] >= 0.68 ~ 0.95 * slopeLength - (landscapeLoc[landStep - 1] * slopeLength), 
                                         landscapeLoc[landStep - 1] >= 0 ~ 0.68 * slopeLength - (landscapeLoc[landStep - 1] * slopeLength), 
                                         TRUE ~ -(landscapeLoc[landStep - 1] * slopeLength))
    asp <- dplyr::case_when(landscapeLoc[landStep - 1] == 1 ~ 100, 
                            landscapeLoc[landStep - 1] >= 0.95 ~ weather[t, "moistureD"], 
                            landscapeLoc[landStep - 1] >= 0.68 ~ weather[t, "moistureC"], 
                            landscapeLoc[landStep - 1] >= 0 ~ weather[t, "moistureB"], 
                            TRUE ~ weather[t, "moistureA"])
    slopeLand <- dplyr::case_when(landscapeLoc[landStep - 1] >= 0.95 ~ slopeM + slopeSD + 2 * slopeSD, 
                                  landscapeLoc[landStep - 1] >= 0.68 ~ slopeM + slopeSD, 
                                  TRUE ~ slopeM)
    if (asp < Extinction) {
      if (asp %in% aspList) {
        R[landStep] <- R[which(asp == aspList)[1]]
        Rs <- NA
        Reach[landStep] <- Reach[which(asp == aspList)[1]]
        sMort[landStep] <- sMort[which(asp == aspList)[1]]
        gMort[landStep] <- gMort[which(asp == aspList)[1]]
        runTest <- TRUE
      } else {
        TBL <- base.params %>% 
          ffm_set_site_param("windSpeed", weather$Wind[t] * (0.0361 * wHeading * slopeM - 5e-04 * wHeading * slopeM^2 + 1) * fEdge, "km/h") %>% 
          ffm_set_site_param("temperature", weather$Temp[t], "degc") %>% 
          ffm_set_site_param("deadFuelMoistureProp", asp) %>% 
          ffm_set_site_param("fireLineLength", as.numeric(firelineW)) %>% 
          ffm_set_site_param("slope", slopeLand * wHeading * fEdge)
        runTest <- ffm_run_robust(base.params=TBL, db.path = dbH, db.recreate = recreate, testN = testN,
                               Strata = Strata, Species = Species, Flora = Flora, Structure = Structure, a = a, l = l, Ms = Ms, Pm = Pm, Mr = Mr)
        if (runTest) {
          res <- ffm_db_load(dbH)
          R[landStep] <- max(res$ROS$ros) * 1000
          Rs <- spotFire(flameHeight = max(res$FlameSummaries$flameHeight), slope = (slopeLand * wHeading * fEdge), FPC, 
                         windExposure = wHeading * landscapeLoc[landStep - 1], 
                         vAir = weather$Wind[t] * (0.0361 * wHeading * slopeM - 5e-04 * wHeading * slopeM^2 + 1) * fEdge, vAir500, fireArea)
          Reach[landStep] <- max(max(res$FlameSummaries$reach <- res$FlameSummaries$flameLength * cos(res$FlameSummaries$flameAngle)), 
                                 max(res$SurfaceResults$reach <- res$SurfaceResults$flameLength * cos(res$SurfaceResults$flameAngle)), Rs)
          litterDepth <- as.numeric(TBL$value[TBL$param == "fuelLoad"])/0.5
          if (sensitive) {
            runs <- suppressMessages(frame::frameSummary(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults))
            IP <- frame::repFlame(res$IgnitionPaths)
            scorch <- suppressMessages(frame::flora(runs, IP, Param = TBL, targSp = targSp, Test = Test))
            sMort[landStep] <- as.numeric(scorch$targSp[1] >= 50) * sensitive
          }
          else {
            sMort[landStep] <- 0
          }
          stems <- data.frame(matrix(ncol = 3, nrow = length(girdleH)))
          colnames(stems) <- c("aspect", "testHeight", "girdle")
          if (is.null(girdleH)) {
            gMort[landStep] <- NA
          }
          else {
            for (Height in girdleH) {
              if (R[landStep] <= minR) {
                stems$girdle <- 0
              }
              else {
                stem <- frame::girdle(Surf = runs, Plant = IP, 
                                      Height = Height, woodDensity = woodDensity, 
                                      phloem = phloem, barkDensity = barkDensity, 
                                      bark = bark, comBark = comBark, resBark = resBark, 
                                      Altitude = Altitude, moisture = moisture, 
                                      bMoisture = bMoisture, distance = distance, 
                                      trail = trail, var = var, diameter = litterDepth, 
                                      RH = weather$RH[t], Pressure = weather$MSLP[t], 
                                      startTemp = weather$Temp[t], necT = necT, 
                                      surfDecl = surfDecl, updateProgress = NULL)
                stems$girdle <- max(stem$girdling, na.rm = TRUE)
              }
              stems$aspect <- asp
              stems$testHeight <- Height
              stemTests <- rbind(stemTests, stems)
              gMort[landStep] <- max(stems$girdle, na.rm = TRUE)
            }
          }
          aspList[landStep] <- asp
        } else {
          R[landStep] <- 0
          Rs <- 0
          Reach[landStep] <- 0
          sMort[landStep] <- 0
          gMort[landStep] <- 0
        }
      }
    } else {
      R[landStep] <- 0
      Rs <- 0
      Reach[landStep] <- 0
      sMort[landStep] <- 0
      gMort[landStep] <- 0
      runTest <- TRUE
    }
    Barrier <- 0
    if (R[landStep] < minR) {
      ll <- vector()
      seg <- 1
      ll[seg] <- landscapeLoc[landStep - 1]
      Barrier <- landformDistance
      wider <- TRUE
      while (wider == TRUE) {
        seg <- seg + 1
        ll[seg] <- dplyr::case_when(ll[seg - 1] == 1 && weather$moistureA[t] >= 0.15 ~ -1, 
                                    (ll[seg - 1] >= 0.95 && weather$moistureD[t] >= 0.15) ~ 1, 
                                    (ll[seg - 1] >= 0.68 && weather$moistureC[t] >= 0.15) ~ 0.95, 
                                    (ll[seg - 1] >= 0 && weather$moistureB[t] >= 0.15) ~ 0.68, 
                                    (ll[seg - 1] < 0 && weather$moistureA[t] >= 0.15) ~ 0, 
                                    TRUE ~ ll[seg - 1])
        wider <- ll[seg] != ll[seg - 1] && Barrier <= slopeLength * 2
        if (wider) {
          nextSeg <- dplyr::case_when(ll[seg] == 1 ~ frame:::river(lRiver, mRiver), 
                                      ll[seg] >= 0.95 ~ slopeLength - (0.95 * slopeLength), 
                                      ll[seg] >= 0.68 ~ 0.95 * slopeLength - (0.68 * slopeLength), 
                                      ll[seg] >= 0 ~ 0.68 * slopeLength, 
                                      TRUE ~ -(ll[seg] * slopeLength))
          Barrier <- Barrier + nextSeg
        }
      }
      if (Barrier <= slopeLength * 2) {
        dfmc <- dplyr::case_when(ll[seg] >= 0.95 ~ weather[t, "moistureA"], 
                                 ll[seg] == 0.68 ~ weather[t, "moistureD"], 
                                 ll[seg] == 0 ~ weather[t, "moistureC"], 
                                 TRUE ~ weather[t, "moistureB"])
      }
      else {
        dfmc <- 100
      }
      spotI <- as.numeric(runif(1) <= min(1, max(0, -2.75 + 20.95 * runif(1, min = 0, max = 0.5) - 32 * dfmc)))
      if (Reach[landStep - 1] * spotI > Barrier) {
        travDist[landStep] <- Barrier
      }
      else {
        travDist[landStep] <- 0
      }
    }
    else {
      travDist[landStep] <- min((R[landStep] * hourStep), landformDistance)
    }
    if (round(landscapeLoc[landStep - 1] + travDist[landStep]/slopeLength, 5) > 1) {
      newL <- round(landscapeLoc[landStep - 1] + travDist[landStep]/slopeLength, 5) - 2
    }
    else {
      newL <- round(landscapeLoc[landStep - 1] + travDist[landStep]/slopeLength, 5)
    }
    landStep <- landStep + 1
    landscapeLoc[landStep - 1] <- newL
    remTime <- if (travDist[landStep - 1] > 0) {
      remTime - min(1, (landformDistance/(R * hourStep)))
    }
    else {
      0
    }
  }
  Step <- as.data.frame(travDist)
  if(!runTest) {
    Step$travDist <- NA
  }
  Step$t <- t
  Step$strikeTime <- strikeTime
  Step$ignitionTime <- weather$Hour[1]
  Step$landscapeLoc <- landscapeLoc
  if(runTest) {
    Step$ROS <- R 
  } else {
    Step$ROS <- NA
  }
  if(runTest) {
    Step$Reach <- Reach
  } else {
    Step$Reach <- NA
  }
  Step$Edge <- aTxt
  if(runTest) {
    Step$ScorchDeath <- sMort
  } else {
    Step$ScorchDeath <- NA
  }
  if(runTest) {
    Step$Girdling <- gMort
  } else {
    Step$Girdling <- NA
  }
  out <- list(Step, stemTests)
  return(out)
}



#' Finds the area burned and effects of a fire
#'
#' @param Flora 
#' @param Structure 
#' @param default.species.params 
#' @param a 
#' @param weather 
#' @param smoulder 
#' @param Extinction 
#' @param hourStep 
#' @param tArea 
#' @param slopeM 
#' @param slopeSD 
#' @param slopeLength 
#' @param rangeDir 
#' @param mRiver 
#' @param lRiver 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param Test 
#' @param sensitive 
#' @param fireN 
#' @param girdleH 
#' @param woodDensity 
#' @param phloem 
#' @param moisture 
#' @param barkDensity 
#' @param bark 
#' @param comBark 
#' @param resBark 
#' @param bMoisture 
#' @param distance 
#' @param trail 
#' @param var 
#' @param Altitude 
#' @param necT 
#' @param surfDecl 
#' @param minR Minimum ROS (m/h)
#' @param vAir500 Multiplier for wind speed at 500hPa
#' @param targSp 
#' @param testN 
#'
#' @return list
#' @export
#'
#' 

burnPrint <- function(Flora, Structure, default.species.params, a, weather, smoulder = 24, Extinction = 0.15, hourStep = 3, tArea = 10,
                      slopeM = 6.7, slopeSD = 4.9, slopeLength = 391, rangeDir = 270, mRiver = 0.38, lRiver = 0.06, 
                      l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, Test = 80, sensitive = 1, fireN = 1,
                      girdleH = c(0.1, 1), woodDensity = 1000, phloem = 0.01, moisture = 0.6,  
                      barkDensity = 300, bark = 0.04, comBark = 700, resBark = 0,  bMoisture = 0.5, vAir500 = 2,
                      distance = 0.1,  trail = 10,  var = 10, Altitude = 200, necT = 60, surfDecl = 10, minR = 1, targSp = "a", testN = 5) {
  
  
  # Starting parameters
  base.params <- suppressWarnings(frame::buildParams(Structure, Flora, default.species.params, a,
                                                     fLine = 1, slopeM, temp = 30, dfmc = 0.05, wind = 10))
  
  # Create empty output vectors
  A <- vector()
  headReach <- vector()
  flankReach <- vector()
  rosH <- vector()
  rosF <- vector()
  H <- vector()
  Fa <- vector()
  Fb <- vector()
  Ta <- vector()
  hour <- vector()
  L <- vector()
  W <- vector()
  mS <- vector()
  mSA <- vector()
  distNS <- vector()
  distWE <- vector()
  distSN <- vector()
  fGrowth <- data.frame()
  gTest <- data.frame()
  
  # Starting values
  L[1] <- 1
  W[1] <- 1
  Stop <- 0
  stopNS <- 0
  stopWE <- 0
  stopSN <- 0
  t <- 0
  fReachA <- 0
  fReachB <- 0
  fReachC <- 0
  fireArea <- 1e-04 # Actively burning area in ha
  
  # Lightning strike, ignition, and following weather conditions
  ig <- frame::ignitions(weather, smoulder = smoulder, Extinction = Extinction)
    igLoc <- ig$landscapeLoc[1]
  # Run fire if ignition occurs, otherwise return empty dataframe
  if (ig$ignition[1]) {
    weatherB <- weather %>%
      mutate(Hour = Hour + max(weather$Hour))
    wEvent <- filter(weather, Hour >= ig$ignitionTime) %>%
      rbind(weatherB)
    
    # Spread fire from ignition point, over time
    
    while(Stop == 0) {
      t <- t+1
      hour[t] <- t * hourStep
      
      wHeading <- case_when(abs(wEvent$Direction[t]-rangeDir) == 270  ~ 1,
                            (wEvent$Direction[t]-rangeDir == 0) || (abs(wEvent$Direction[t]-rangeDir) == 180) ~ 0,
                            TRUE ~ -1) # 1 ~ upslope, 0 ~ along, -1 ~ down
      
      cat("Time", hour[t], "\n")
      
      # Build new pseudo-transect
      tbl <- frame::specPoint(base.params, Structure, a)
      Strata <- frame::strata(tbl)
      Species <- frame::species(tbl)
      
      # Find foliage projective cover
      for (st in 1:nrow(Strata)) {
        if (st == 1) {
          FPC <- Strata$cover[st]
        } else {
          FPC <- FPC + ((1-FPC) * Strata$cover[st])
        }
      }
      
      TBL <- frame::plantVarFrame(tbl, Strata, Species, Flora, a, l,
                                  Ms = Ms, Pm = Pm, Mr = Mr)
      
      # Model edges
      
      #N-S
      if (stopNS < smoulder) {
        fEdge <- case_when(wEvent$Direction[t] == 0 ~ 1,  #1~Head, 0~Flank, -1~Tail
                           wEvent$Direction[t] == 180 ~ -1,
                           TRUE ~ 0)
        edgeA <- frameSpread(base.params = TBL, weather = wEvent, 
                             a, igLoc = igLoc, t = t, hourStep = hourStep, strikeTime = ig$strikeTime[1], 
                             vAir500, FPC, fireArea, slopeM = slopeM, slopeSD = slopeSD, 
                             slopeLength = slopeLength, mRiver = mRiver, 
                             lRiver = lRiver, wHeading = wHeading, firelineW = W[t], 
                             fReach = fReachA, fEdge = fEdge, Extinction = Extinction, 
                             Test = Test, sensitive = sensitive, girdleH = girdleH, 
                             woodDensity = woodDensity, phloem = phloem, 
                             moisture = moisture, barkDensity = barkDensity, 
                             bark = bark, comBark = comBark, resBark = resBark, 
                             bMoisture = bMoisture, distance = distance, 
                             trail = trail, var = var, Altitude = Altitude, 
                             necT = necT, surfDecl = surfDecl, minR = minR, 
                             targSp = targSp, testN = testN, Strata = Strata, Species = Species, 
                             Flora = Flora, Structure = Structure, l = l, Ms = Ms, Pm = Pm, Mr = Mr)
        tabA <- edgeA[[1]] %>%
          mutate(fAge = Structure$site[1],
                 fNo = fireN,
                 tStep = t,
                 fAspect = "N")
        fGrowth <- rbind(fGrowth,tabA)
        if (length(girdleH)>0) {
          tabB <- edgeA[[2]] %>%
            mutate(fAge = Structure$site[1],
                   fNo = fireN,
                   tStep = t,
                   fAspect = "N")
          gTest <- rbind(gTest, tabB)
        } 
        
        finA <- filter(tabA, row_number()==n())
        igLoc <- finA$landscapeLoc[1]
        fReachA <- finA$Reach[1]
        distNS[t] <- sum(tabA$travDist)
      } else {
        distNS[t] <- 0
      }
      
      #Smoulder counter
      if (distNS[t] == 0) {
        stopNS <- stopNS+1
      } else {
        stopNS <- 0
      }
      
      #W-E
      if (stopWE < smoulder) {
        fEdge <- case_when(wEvent$Direction[t] == 270 ~ 1,  #1~Head, 0~Flank, -1~Tail
                           wEvent$Direction[t] == 90 ~ -1,
                           TRUE ~ 0)
        edgeB <- frameSpread(base.params = TBL, weather = wEvent, 
                             a, igLoc = igLoc, t = t, hourStep = hourStep, strikeTime = ig$strikeTime[1],
                             vAir500, FPC, fireArea, slopeM = slopeM, slopeSD = slopeSD, 
                             slopeLength = slopeLength, mRiver = mRiver, 
                             lRiver = lRiver, wHeading = wHeading, firelineW = L[t], 
                             fReach = fReachB, fEdge = fEdge, Extinction = Extinction, 
                             Test = Test, sensitive = sensitive, girdleH = girdleH, 
                             woodDensity = woodDensity, phloem = phloem, 
                             moisture = moisture, barkDensity = barkDensity, 
                             bark = bark, comBark = comBark, resBark = resBark, 
                             bMoisture = bMoisture, distance = distance, 
                             trail = trail, var = var, Altitude = Altitude, 
                             necT = necT, surfDecl = surfDecl, minR = minR, 
                             targSp = targSp, testN = testN, Strata = Strata, Species = Species, 
                             Flora = Flora, Structure = Structure, l = l, Ms = Ms, Pm = Pm, Mr = Mr) 
        tabA <- edgeB[[1]] %>%
          mutate(fAge = Structure$site[1],
                 fNo = fireN,
                 tStep = t,
                 fAspect = "W")
        fGrowth <- rbind(fGrowth,tabA)
        if (length(girdleH)>0) {
          tabB <- edgeB[[2]] %>%
            mutate(fAge = Structure$site[1],
                   fNo = fireN,
                   tStep = t,
                   fAspect = "W")
          gTest <- rbind(gTest, tabB)
        }
        
        finB <- filter(tabA, row_number()==n())
        igLoc <- finB$landscapeLoc[1]
        fReachB <- finB$Reach[1]
        distWE[t] <- sum(tabA$travDist)
      } else {
        distWE[t] <- 0
      }
      
      #Smoulder counter
      if (distWE[t] == 0) {
        stopWE <- stopWE+1
      } else {
        stopWE <- 0
      }
      
      #S-N
      if (stopSN < smoulder) {
        fEdge <- case_when(wEvent$Direction[t] == 180 ~ 1,  #1~Head, 0~Flank, -1~Tail
                           wEvent$Direction[t] == 0 ~ -1,
                           TRUE ~ 0)
        edgeC <- frameSpread(base.params = TBL, weather = wEvent, 
                             a, igLoc = igLoc, t = t, hourStep = hourStep, strikeTime = ig$strikeTime[1],
                             vAir500, FPC, fireArea, slopeM = slopeM, slopeSD = slopeSD, 
                             slopeLength = slopeLength, mRiver = mRiver, 
                             lRiver = lRiver, wHeading = wHeading, firelineW = W[t], 
                             fReach = fReachC, fEdge = fEdge, Extinction = Extinction, 
                             Test = Test, sensitive = sensitive, girdleH = girdleH, 
                             woodDensity = woodDensity, phloem = phloem, 
                             moisture = moisture, barkDensity = barkDensity, 
                             bark = bark, comBark = comBark, resBark = resBark, 
                             bMoisture = bMoisture, distance = distance, 
                             trail = trail, var = var, Altitude = Altitude, 
                             necT = necT, surfDecl = surfDecl, minR = minR, 
                             targSp = targSp, testN = testN, Strata = Strata, Species = Species, 
                             Flora = Flora, Structure = Structure, l = l, Ms = Ms, Pm = Pm, Mr = Mr) 
        tabA <- edgeC[[1]] %>%
          mutate(fAge = Structure$site[1],
                 fNo = fireN,
                 tStep = t,
                 fAspect = "W")
        fGrowth <- rbind(fGrowth,tabA)
        if (length(girdleH)>0) {
          tabB <- edgeC[[2]] %>%
            mutate(fAge = Structure$site[1],
                   fNo = fireN,
                   tStep = t,
                   fAspect = "W")
          gTest <- rbind(gTest, tabB)
        }
        
        finC <- filter(tabA, row_number()==n())
        igLoc <- finC$landscapeLoc[1]
        fReachC <- finC$Reach[1]
        distSN[t] <- sum(tabA$travDist)
      } else {
        distSN[t] <- 0
      }
      
      #Smoulder counter
      if (distSN[t] == 0) {
        stopSN <- stopSN+1
      } else {
        stopSN <- 0
      }
      
      # Calculate burnt area using elliptical spread
      L[t+1] <- L[t]+distNS[t] + distSN[t]
      W[t+1] <- W[t]+2*distWE[t]
      A[t] <- min((pi*0.5*L[t+1]*0.5*W[t+1]), tArea)
      tArea <- tArea - (A[t]-A[t-1])
      fireArea <- ((A[t]-A[t-1])/hourStep)/10000
      
      
      # End loop if fire has stopped spreading
      if (t > 1) {
        Stop <- as.numeric(A[t] == A[t-1])   
      }
    }
    
    event <- as.data.frame(hour)
    event$Area <- round(A,0)
    event$fire <- fireN
  } else {
    # Values if fire fails to spread
    hour <- 0
    event <- as.data.frame(hour)
    event$Area <- 0
    event$fire <- fireN
    gTest <- as.data.frame(matrix(nrow = 1, ncol = 3))
    colnames(gTest) <- c('girdle', 'aspect', 'testHeight')
    gTest$girdle <- 0
    fGrowth <- as.data.frame(matrix(nrow = 1, ncol = 14))
    colnames(fGrowth) <- c('travDist', 't', 'strikeTime', 'ignitionTime', 'landscapeLoc', 'ROS', 'Reach', 'Edge', 'ScorchDeath', 'Girdling', 'fAge', 'fNo', 'tStep', 'fAspect')
    fGrowth$travDist <- 0
    fGrowth$t <- 1
    fGrowth$strikeTime <- ig$strikeTime[1]
    fGrowth$ignitionTime <- ig$ignitionTime[1]
    fGrowth$landscapeLoc <- igLoc
    fGrowth$ROS <- 0
    fGrowth$Reach <- 0
    fGrowth$ScorchDeath <- 0
    fGrowth$Girdling <- 0
    fGrowth$fAge <- Structure$site[1]
    fGrowth$fNo <- fireN
  }
    # Invalidate values if one edge failed to spread
    if (is.na(mean(event$Area))) {
      event$Area <- NA
    }
    
    if (nrow(fGrowth)>0) {
      if (is.na(mean(fGrowth$travDist))) {
        fGrowth$travDist <- NA
      }
      if (is.na(mean(fGrowth$ROS))) {
        fGrowth$ROS <- NA
      }
      if (is.na(mean(fGrowth$Reach))) {
        fGrowth$Reach <- NA
      }
      if (is.na(mean(fGrowth$ScorchDeath))) {
        fGrowth$ScorchDeath <- NA
      }
      if (is.na(mean(fGrowth$Girdling))) {
        fGrowth$Girdling <- NA
      }
    }
    
    if (nrow(gTest)>0) {
      if (is.na(mean(gTest$girdle))) {
        gTest$girdle <- NA
      }
    }
  
  out <- list(event, fGrowth, gTest)
  return(out)
}


#' Find likelihood and consequence of death for canopy trees
#'
#' @param Flora 
#' @param Structure 
#' @param default.species.params 
#' @param a 
#' @param weather 
#' @param lightning 
#' @param smoulder 
#' @param Extinction 
#' @param hourStep 
#' @param slopeM 
#' @param slopeSD 
#' @param slopeLength 
#' @param rangeDir 
#' @param mRiver 
#' @param lRiver 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param Test 
#' @param sensitive 
#' @param girdleH 
#' @param woodDensity 
#' @param phloem 
#' @param moisture 
#' @param barkDensity 
#' @param bark 
#' @param comBark 
#' @param resBark 
#' @param bMoisture 
#' @param distance 
#' @param trail 
#' @param var 
#' @param Altitude 
#' @param necT 
#' @param surfDecl 
#' @param vAir500 
#' @param minR Minimum ROS (m/h)
#' @param targSp 
#' @param testN 
#' @param fires Number of fires to model for each age
#'
#' @return list
#' @export
#'

frameRisk <- function(Flora, Structure, default.species.params, a = a, weather, lightning = 0.05, fires = 5,
                      smoulder = 24, Extinction = 0.15, hourStep = 3, 
                      slopeM = 6.7, slopeSD = 4.9, slopeLength = 391, rangeDir = 270, mRiver = 0.38, lRiver = 0.06, 
                      l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, Test = 80, sensitive = 1, 
                      girdleH = c(), woodDensity = 1000, phloem = 0.01, moisture = 0.6, vAir500 = 2,
                      barkDensity = 300, bark = 0.04, comBark = 700, resBark = 0,  bMoisture = 0.5,
                      distance = 0.1,  trail = 10,  var = 10, Altitude = 200, necT = 60, surfDecl = 10, minR = 1, 
                      targSp = "a", testN = 5) {
  
  # Starting values
  studyArea <- fires/lightning
  tArea <- studyArea * 1000000
  event <- data.frame()
  fGrowth <- data.frame()
  riskTab <- data.frame(matrix(nrow = 1, ncol = 4))
  colnames(riskTab) <- c('Age', 'Likelihood', 'Consequence', 'Risk')
  gird <- data.frame()
  
  fireN <- 1
  while (fireN <= fires) {
    cat("Running fire", fireN, "of", fires,"\r")
    incident <- burnPrint(Flora, Structure, default.species.params, a = a, weather, smoulder = smoulder, Extinction = Extinction, hourStep = hourStep, tArea = tArea,
                          slopeM = slopeM, slopeSD = slopeSD, slopeLength = slopeLength, rangeDir = rangeDir, mRiver = mRiver, lRiver = lRiver, 
                          l = l, Ms = Ms, Pm = Pm, Mr = Mr, Test = Test, sensitive = sensitive, fireN = fireN,
                          girdleH = girdleH, woodDensity = woodDensity, phloem = phloem, moisture = moisture, vAir500 = vAir500,
                          barkDensity = barkDensity, bark = bark, comBark = comBark, resBark = resBark,  bMoisture = bMoisture,
                          distance = distance,  trail = trail,  var = var, Altitude = Altitude, necT = necT, surfDecl = surfDecl, minR = minR, 
                          targSp = targSp, testN = testN)
    # Cancel the function if there is an error
    if (is.error(incident)) {
      out <- list(riskTab, event, fGrowth, gird)
      return(out)
      stop()
    } 
    fireN <- fireN+1
    tArea <- max(tArea - incident[[1]][nrow(incident[[1]]), 2],0) # Partial solution. Lack of spatial relationship assumes that each fire can pick out available unburned veg
    
    # End loop if study area all burnt
    if (tArea == 0) {
      fireN <- fires
    }
    
    # Collate outputs
    event <- rbind(event,incident[[1]]) 
    fGrowth <- rbind(fGrowth,incident[[2]])
    gird  <- rbind(gird,incident[[3]])
  }
  
  fGrowth$Mortality <- fGrowth$travDist*pmax(fGrowth$ScorchDeath,fGrowth$Girdling, na.rm = TRUE)
  fGrowth$wMortality <- fGrowth$Mortality * fGrowth$travDist
  event$fAge <- fGrowth$fAge[1]
  for (step in 1:nrow(event)) {
    if (step == 1) {
      event$growth[step] <- event$Area[step]
    } else {
      if (event$fire[step] == event$fire[step-1]) {
        event$growth[step] <- event$Area[step] - event$Area[step-1]
      } else {
        event$growth[step] <- event$Area[step]
      }
    }
  }
  
  event$cumHa <- round(cumsum(event$growth)/10000,1)
  riskTab$Age <- fGrowth$fAge[1]
  riskTab$Likelihood <- event$cumHa[nrow(event)]/(studyArea*100)
  riskTab$Consequence <- round((sum(fGrowth$wMortality, na.rm = TRUE)/sum(fGrowth$travDist, na.rm = TRUE))/100,4)
  riskTab$Risk <- round(riskTab$Likelihood * riskTab$Consequence,4)
  
  out <- list(riskTab, event, fGrowth, gird)
  return(out)
}


#' Internal function for riskDynamics
#'
#' @param a 
#'
#' @return list

parRisk <- function(a) {
  FloraA <- filter(Flora, record == a)
  StructureA <- filter(Structure, record == a)
  base.params <- suppressWarnings(frame::buildParams(StructureA, FloraA, default.species.params, a,
                                                     fLine = 1, slopeM, temp = 30, dfmc = 0.05, wind = 10))
  
  hCan <- max(FloraA$top, na.rm = TRUE)
  LAI <- LAIcomm(base.params, yu = hCan, yl = 0)
  WRF <- windReduction(base.params, test = 1.2)
  litterW <- max(FloraA$weight, na.rm = TRUE)
  
  weather <- frame::frameWeather(clim = clim, m, LAI, WRF, hCan, rholitter, litterW,
                                 lat, slope = slopeM, slopeSD, rangeDir, cardinal)
  
  res <- frame::frameRisk(Flora = FloraA, Structure = StructureA, default.species.params, a, weather, lightning, fires,
                          smoulder, Extinction, hourStep, slopeM, slopeSD, slopeLength, rangeDir, mRiver, lRiver, 
                          l, Ms, Pm, Mr, Test, sensitive, girdleH, woodDensity, phloem, moisture, vAir500,  
                          barkDensity, bark, comBark, resBark,  bMoisture,
                          distance,  trail,  var, Altitude, necT, surfDecl, minR, targSp, testN)
  
  nA <- paste0("ARa",a)
  nB <- paste0("ARb",a)
  nC <- paste0("ARc",a)
  nD <- paste0("ARd",a)
  
  nA <- res[[1]]
  nB <- res[[2]]
  nC <- res[[3]]
  nD <- res[[4]]
  
  out <- list(nA, nB, nC, nD)
  return(out)
}


#' Calculates forest risk parameters over a growth series
#' using parallel ports
#'
#' @param fireDat List of three dataframes: Flora, Structure & default.species.params
#' @param clim Formatted AM & PM weather conditions over a time sequence
#' @param lightning Number of lightning strikes per km2
#' @param lat Latitude (degrees)
#' @param Altitude Altitude (m.a.s.l.)
#' @param slopeM Mean slope (degrees)
#' @param slopeSD Standard deviation of slope (degrees)
#' @param slopeLength Distance from ridge to gully (m)
#' @param rangeDir 
#' @param mRiver 
#' @param lRiver 
#' @param cardinal 
#' @param smoulder 
#' @param Extinction 
#' @param m 
#' @param rholitter 
#' @param hourStep 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param Test 
#' @param sensitive 
#' @param girdleH 
#' @param phloem 
#' @param moisture 
#' @param barkDensity 
#' @param bark 
#' @param comBark 
#' @param resBark 
#' @param bMoisture 
#' @param vAir500 
#' @param testN 
#' @param distance 
#' @param trail 
#' @param var 
#' @param necT 
#' @param surfDecl 
#' @param minR 
#' @param freeCores 
#' @param tr 
#' @param fires Number of fires to model per age
#' @param woodDensity 
#' @param targSp 
#'
#' @return list
#' @export
#'

riskDynamics <- function(fireDat, tr, clim, lightning = 0.05, fires = 5, lat = -42.48, Altitude = 200, 
                         slopeM = 6.7, slopeSD = 4.9, slopeLength = 391, rangeDir = 270, mRiver = 0.38, lRiver = 0.06,   
                         cardinal = TRUE, smoulder = 24, Extinction = 0.15, m = 0.15, rholitter = 550, hourStep = 3, 
                         l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, Test = 80, sensitive = 1, girdleH = c(), phloem = 0.01, moisture = 0.6,  
                         barkDensity = 300, bark = 0.04, comBark = 700, resBark = 0,  bMoisture = 0.5, vAir500 = 2, testN = 5,
                         distance = 2,  trail = 1000,  var = 10, necT = 60, surfDecl = 10, minR = 0.001, freeCores = 1){
  
  ages <- max(fireDat[[2]]$record)
  Flora <- fireDat[[1]]
  Structure <- fireDat[[2]]
  default.species.params <- fireDat[[3]]
  
  # Collect inputs
  woodDensity <- tr$rho[tr$hmat == max(tr$hmat)]
  targSp <- Flora$species[which(Flora$top == max(Flora$top, na.rm = TRUE))]
  
  # Create a cluster of cores with replicated R on each
  nCores <- max(parallel::detectCores() - freeCores,1)
  cl <- parallel::makeCluster(nCores)
  # Load the packages
  parallel::clusterEvalQ(cl,
                         { library(dplyr)
                           library(tidyr)
                           library(frame)
                           library(assertthat)
                           library(extraDistr)})
  # Load the inputs
  parallel::clusterExport(cl,varlist=c('Flora', 'Structure', 'default.species.params', 'clim', 'lightning', 'fires', 'lat', 'Altitude',
                                       'slopeM', 'slopeSD','slopeLength', 'rangeDir', 'mRiver', 'lRiver', 'cardinal',
                                       'smoulder', 'Extinction', 'm', 'rholitter', 'hourStep',
                                       'l', 'Ms', 'Pm', 'Mr', 'Test', 'sensitive', 'girdleH', 'woodDensity', 
                                       'phloem', 'moisture', 'barkDensity', 'bark', 'comBark', 'resBark',  'bMoisture', 'vAir500', 'testN',
                                       'distance',  'trail',  'var', 'necT', 'surfDecl', 'minR', 'targSp'), environment())
  
  system.time(parT <- parallel::parLapply(cl, 1:ages, parRisk))
  parallel::stopCluster(cl)
  
  # Risk table
  riskDyn <- data.frame()
  event <- data.frame()
  fGrowth <- data.frame()
  gird <- data.frame()
  for (n in 1:ages) {
    Na <- as.data.frame(parT[[n]][1])
    Nb <- as.data.frame(parT[[n]][2])
    Nc <- as.data.frame(parT[[n]][3])
    Nd <- as.data.frame(parT[[n]][4])
    riskDyn <- rbind(riskDyn,Na)
    event <- rbind(event,Nb)
    fGrowth <- rbind(fGrowth,Nc)
    gird <- rbind(gird,Nd)
  }
  out <- list(riskDyn, event, fGrowth, gird)
  return(out)
}


#' Calculates risk (likelihood & consequence) for plants in a given climate & terrain
#'
#' @param dynDat 
#' @param a 
#' @param ignitions 
#' @param hourStep 
#' @param Area 
#' @param weather 
#' @param vAir 
#' @param wStDev 
#' @param tAir 
#' @param dfmc 
#' @param slope 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param Test 
#' @param fireN 
#'
#' @return list
#'

plantRisk <- function(dynDat, a, ignitions = 1, hourStep = 3, Area = 10, weather,
                       vAir = 25, wStDev = 5, tAir = 30, dfmc = 0.05, slope = 0,
                       l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, Test = 80, fireN = 1) {
  
  Flora <- filter(dynDat[[1]], record == a)
  Structure <- filter(dynDat[[2]], record == a)
  default.species.params <- dynDat[[3]]
  events <- data.frame()
  
  tArea <- Area
  f <- 1
  while (f <= (ignitions*Area)) {
    cat("Fire", f, "of", (ignitions*Area), "fires", "\n")
    fire <- burnPrint(Flora, Structure, default.species.params, a, hourStep, tArea,
                      vAir, wStDev, tAir, dfmc, slope,
                      l, Ms, Pm, Mr, Test, fireN = f)
    tArea <- max(0, tArea - max(fire$Area))
    cat("Remaining area is", tArea, "of", Area, "ha", "\n")
    if (tArea == 0) {
      f <- (ignitions*Area)
    }
    f <- f+1
    events <- rbind(events, fire)
  }
  
  risk <- data.frame(matrix(ncol = 3, nrow = 1))
  colnames(risk) <- c("Likelihood", "Girdle", "Scorch")
  risk$Likelihood <- (Area - tArea)/Area
  risk$Scorch <- sum(events$scorchA)/Area
  
  return(list(events, risk))
}



#' Calculates spotting distance
#'
#' @param flameHeight Flame height (m)
#' @param slope Degrees
#' @param FPC Foliage Projective Cover (0 - 100 percent)
#' @param windExposure <1 protected, >=1 Exposed
#' @param vAir Wind speed (km/h)
#' @param vAir500 Multiplier of wind speed to estimate 500 hPA wind
#' @param fireArea Area of active fire (ha)
#'
#' @return value
#' @export
#'

spotFire <- function(flameHeight, slope, FPC, windExposure, vAir, vAir500 = 2, fireArea){
  
  # Max distance from Storey et al (2020)
  spMax <- exp(-2.762602963 + 0.013923896*slope + 0.014596962*FPC + 3.439632761*windExposure +
                   0.032052944*vAir + 0.005534678*vAir500*vAir + 0.423010927*log(fireArea))
  
  # Approx distance from Gould et al (2007)
  spG <- 11.98*flameHeight^2.19
  
  return(min(spMax, spG))
}


#' Internal function for fireDynamics
#'
#' @param a 
#'
#' @return list
#' @export
#'

parBurn <- function(n) {
  
  f <- filter(Flora, record == n)
  s <- filter(Structure, record == n)
  base.params <- suppressWarnings(frame::buildParams(Structure = s, Flora = f, default.species.params, a = n,
                                                     fLine = fLine, slope = slope, temp = temp, dfmc = DFMC, wind = wind))
  
  Strata <- strata(base.params)
  db.path <- paste("age",n,".db", sep = "")
  
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
  
  return(list(outa, IP))
  
}


#' Models probabilistic fire behaviour across multiple ages on parallel cores
#'
#' @param fireDat 
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
#' @param fLine 
#' @param updateProgress 
#' @param freeCores 
#' @param l 
#' @param Ms 
#' @param Pm 
#' @param Mr 
#' @param vAir500 
#'
#' @return list
#' @export
#'

fireDynamics <- function(fireDat, slope = 5, slopeSD = 2, slopeRange = 5, 
                         temp = 30, tempSD = 5, tempRange = 3,
                         DFMC = 0.1, DFMCSD = 0.01, DFMCRange = 2, 
                         wind = 10, windSD = 5, windRange = 5, fLine = 1000,
                         moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.5,
                         reps = 5, leafVar = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, vAir500 = 2, 
                         testN = 5, updateProgress = NULL, freeCores = 1){
  
  cat("Modelling risk", "\n", "\n")
  
  # 1. Compile inputs
  Flora <- fireDat[[1]]
  Structure <- fireDat[[2]]
  default.species.params <- fireDat[[3]]
  r <- unique(Flora$record)
  
  # 2. Create a cluster of cores with replicated R on each
  nCores <- max(parallel::detectCores() - freeCores,1)
  cl <- parallel::makeCluster(nCores)
  # 3. Load the packages
  parallel::clusterEvalQ(cl,
                         { library(dplyr)
                           library(tidyr)
                           library(frame)
                           library(assertthat)
                           library(extraDistr)})
  # 4. Load the inputs
  parallel::clusterExport(cl,varlist=c('Flora', 'Structure', 'default.species.params', 
                                       'slope', 'slopeSD', 'slopeRange', 
                                       'temp', 'tempSD', 'tempRange',
                                       'DFMC', 'DFMCSD', 'DFMCRange', 
                                       'wind', 'windSD', 'windRange', 'fLine',
                                       'moistureMultiplier', 'moistureSD', 'moistureRange',
                                       'reps', 'leafVar', 'Ms', 'Pm', 'Mr', 'vAir500', 
                                       'updateProgress'), environment())
  
  # 5. Send each rep to a different core to be processed
  system.time(parT <- parallel::parLapply(cl, r, parBurn))
  parallel::stopCluster(cl)
  
  # 6. Summarise and export results
  runs <- data.frame()
  IP <- data.frame()
  for (n in 1:length(r)) {
    Na <- as.data.frame(parT[[n]][1])
    Nb <- as.data.frame(parT[[n]][2])
    runs <- rbind(runs,Na)
    IP <- rbind(IP,Nb)
  }
  out <- list(runs, IP)
  return(out)
}