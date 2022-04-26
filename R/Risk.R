#' Calculate the expected area burned
#'
#' Ember ignition likelihood calculated from 
#' Ellis P F 2011 Fuelbed ignition potential and bark morphology explain the 
#' notoriety of the eucalypt messmate “stringybark” for intense spotting Int. J. Wildl. Fire 20 897–907
#' 
#' @param a A unique identifier for the record being run
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
#' site - Name of the site
#' species - the name of the species, which will call trait data from 'default.species.params'
#' moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
#' stratum - numeric value from 1 to 4, counting from lowest stratum
#' comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
#' base, he, ht, top & w - canopy dimensions for that species (m)
#' Hs - Standard deviation of plant heights
#' Hr - Range of plant heights
#' weight - weight of litter (t/h)
#' diameter - mean diameter of litter (m)
#' @param default.species.params Leaf traits database
#' @param hourStep Number of hours in each calculation timestep
#' @param Area Area of the study area (sq km)
#' @param vAir
#' @param wStDev
#' @param tAir
#' @param dfmc
#' @param slope
#' @param fireN Number of the fire being studied
#' @param l Variation around input leaf dimensions
#' @param Ms Standard deviation of LFMC
#' @param Pm Multiplier of mean LFMC
#' @param Mr Truncates LFMC variability by +/- Mr * LFMC
#' @return dataframe
#' @export


burnPrint <- function(Flora, Structure, default.species.params, a, hourStep = 3, Area = 10,
                      vAir = 25, wStDev = 5, tAir = 30, dfmc = 0.05, slope = 0,
                      l = 0.1, Ms = 0.01, Pm = 1, Mr = 1.5, fireN = 1) {
  
  
  # Starting parameters
  base.params <- suppressWarnings(frame::buildParams(Structure, Flora, default.species.params, a,
                                                     fLine = 1, slope, temp = tAir, dfmc = dfmc, wind = vAir))
  
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
  
  # Starting values
  L[1] <- 0.001
  W[1] <- 0.001
  Stop <- 0
  t <- 0
  
  while(Stop == 0) {
    t <- t+1
    hour[t] <- t * hourStep
    cat("Step", t, "\n")
    # Model behaviour for a random point
    tbl <- frame::specPoint(base.params, Structure, a) %>%
      ffm_set_site_param("windSpeed", vAir, "km/h") %>%
      ffm_set_site_param("temperature", tAir, "degc") %>%
      ffm_set_site_param("deadFuelMoistureProp", dfmc) %>%
      ffm_set_site_param("fireLineLength", as.numeric(W[t]))
    Strata <- frame::strata(tbl)
    Species <- frame::species(tbl)
    
    # Choose random point and vary plant traits for each species within their range
    TBL <- frame::plantVarFrame(tbl, Strata, Species, Flora, a, l,
                                Ms = Ms, Pm = Pm, Mr = Mr)
    # Head fire
    ffm_run(TBL, db.path = "out.head.db", db.recreate = TRUE)
    resH<-ffm_db_load("out.head.db")
    rosH[t] <- max(resH$ROS$ros)
    
    # Spotting distance from function fit to Vesta data
    spotH <- 11.98*max(resH$FlameSummaries$flameHeight)^2.19
    
    # Likelihood of ignition by spotting
    spotI <- as.numeric(runif(1) <= min(1,max(0,-2.75 + 20.95*runif(1, min = 0, max = 0.5) - 32*dfmc)))
    headReach[t] <- max(max(resH$FlameSummaries$reach <- resH$FlameSummaries$flameLength * cos(resH$FlameSummaries$flameAngle)),
                        max(resH$SurfaceResults$reach <- resH$SurfaceResults$flameLength * cos(resH$SurfaceResults$flameAngle)),spotH * spotI, 0.5)
    
    # Flanks
    fWind <- max(0,rnorm(1, mean = 0, sd = wStDev))
    TBL <- TBL %>%
      ffm_set_site_param("windSpeed", fWind, "km/h")
    ffm_run(TBL, db.path = "out.flank.db", db.recreate = TRUE)
    resF<-ffm_db_load("out.flank.db")
    rosF[t] <- max(resF$ROS$ros)
    flankReach[t] <- max(max(resF$FlameSummaries$reach <- resF$FlameSummaries$flameLength * cos(resF$FlameSummaries$flameAngle)),
                         max(resF$SurfaceResults$reach <- resF$SurfaceResults$flameLength * cos(resF$SurfaceResults$flameAngle)),0.5)
    
    # Determine whether edges can cross breaks
    H[t] <- rosH[t]*as.numeric(runif(1) > (1/(1+headReach[t])))
    Fa[t] <- rosF[t]*as.numeric(runif(1) > (1/(1+flankReach[t])))
    Fb[t] <- rosF[t]*as.numeric(runif(1) > (1/(1+flankReach[t])))
    Ta[t] <- rosF[t]*as.numeric(runif(1) > (1/(1+flankReach[t])))
    
    # Prevent spread on edges already blocked
    # Prop flank blocked
    fBlock <- (hourStep*H[t])/L[t]
    # Prop front blocked
    frBlock <- (hourStep*Fa[t] + hourStep*Fb[t])/W[t]
    if (t > 1) {
      if (H[t-1] == 0) {
        H[t] <- frBlock*H[t]
      }
      if (Fa[t-1] == 0) {
        Fa[t] <- fBlock*Fa[t]
      }
      if (Fb[t-1] == 0) {
        Fb[t] <- fBlock*Fb[t]
      }
      if (Ta[t-1] == 0) {
        Ta[t] <- frBlock*Ta[t]
      }
    }
    
    # Calculate burnt area using elliptical spread
    L[t+1] <- L[t]+hourStep*H[t]+(hourStep*Ta[t])
    W[t+1] <- W[t]+(hourStep*Fa[t])+(hourStep*Fb[t])
    A[t] <- min((pi*0.5*L[t+1]*0.5*W[t+1]), Area)
    
    # End loop if fire has stopped spreading
    if (t > 1) {
      Stop <- as.numeric(A[t] == A[t-1])   
    }
    
    # Find the width of the new fire front
    #    fFront[t+1] <- fFront[t] + W[t+1]*1000
    
  }
  event <- as.data.frame(hour)
  event$Area <- A
  event$ROSh <- rosH
  event$ROSF <- rosF
  event$Length <- tail ( L, -1)
  event$Width <- tail ( W, -1)
  event$headReach <- headReach
  event$flankReach <- flankReach
  event$fire <- fireN
  
  return(event)
}