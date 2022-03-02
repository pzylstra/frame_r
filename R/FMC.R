
#' Finds solar radiation reaching the ground surface
#' 
#' @param lat Latitude in degrees
#' @param hr Hour of the day, 0-24, decimal format
#' @param month Name of the month, capitalised (e.g. "January")
#' @param aspect Aspect in degrees
#' @param slope Slope in degrees
#' @param LAI Leaf Area Index
#' @param cloud Proportion cloud cover (0-1)
#' @export

insol <- function(lat = -32, hr = 12, month = "December", aspect = 0, slope = 0, LAI = 3, cloud = 0){
  radians <- function(deg) {(deg * pi) / (180)}
  hAngle <- ifelse(((15*(12-hr))<0), (360+(15*(12-hr))),(15*(12-hr)))
  decl <- dplyr::case_when(month == "March" ~ 0,
                           month == "April" ~ -11.75,
                           month == "May" ~ -20.35,
                           month == "June" ~ -23.50,
                           month == "July" ~ -20.35,
                           month == "August" ~ -11.75,
                           month == "September" ~ 0,
                           month == "October" ~ 11.75,
                           month == "November" ~ 20.35,
                           month == "December" ~ 23.50,
                           month == "January" ~ 20.35,
                           TRUE ~ 11.75)
  terrainE <- (radians(cos(radians(hAngle-aspect))*slope))
  zenith <- acos(sin(radians(abs(lat)))*sin(radians(decl))+cos(radians(abs(lat)))*cos(radians(decl))*cos(radians(hAngle)))-terrainE
  shade <- 1-exp(-(2/(pi*tan(max(0.01,1.570796-zenith))))*LAI)
  out <- max(0,((1-shade)*(1-cloud)*1000*cos(zenith)))
  return(out)
}

# Finds saturation specific humidity from Stull (1988) Eqn 7.5.2c and d

QSat <- function(tAir = 293.46, pAir = 101180){
  esat <- 611.2 * exp(17.67 * (tAir - 273.15) / (tAir - 29.66))
  Q <- 0.622 * (611.2 * exp(17.67 * (tAir - 273.15) / (tAir - 29.66))) / pAir
  return(Q)
}

# Finds water vapour flux
# Calculates RH at the litter surface from inversion of Nelson (1984) EMC equation

Ema <- function(m, nelsonA = 5.2, nelsonB = -19, conLitter = 0.0006, tAir = 293.46, vAir = 3.52, pAir = 101180, insolation = 0, rhAir = 58){
  SpecHum <- rhAir/100 *QSat(tAir, pAir) # Specific humidity
  tLeaf <- tAir + (8.15 - 2.25 * exp(-0.6 * vAir) - 0.0312 * tAir + (0.021 + (-0.04 + 0.0006 * tAir - 0.00000125 * tAir ^ 2) * exp(-0.6 * vAir)) * insolation) # temperature of the leaf
  fA <- -((-nelsonA / nelsonB) - 0.001 * (tLeaf - 293.15)) * nelsonB   
  rhLeaf <- exp((-18.0153 / (1.9872 * tLeaf)) * exp(m * nelsonB + fA))
  out <- 1.14 * conLitter * (rhLeaf * QSat(tLeaf, pAir) - SpecHum) # Water vapour flux from litter
  return(out)}


#' Finds the moisture content of dead leaf litter at a given time-step
#' 
#' Function uses the 'Single differential equation model' of Matthews et al (2010)
#' 
#' @param m The starting moisture (proportion ODW)
#' @param nelsonA Constant from Nelson FMC model
#' @param nelsonB Constant from Nelson FMC model
#' @param tAir Atmospheric temperature (K)
#' @param vAirWind speed (m/s)
#' @param pAir Atmospheric pressure (Pa)
#' @param rhAir Relative humidity (%)
#' @param sigma Surface area: volume ratio of litter particles (m2/m3)
#' @param rhoLitter Density of leaves (kg/m3)
#' @param conLitter Litter conductance (m/s)
#' @param dt Seconds per time step
#' @param EPS Sensitivity of iteration
#' @param insolation Solar energy at the soil surface (kW/m2)
#' @export

simplefmc <- function(m, tAir = 293.46, vAir = 3.52, pAir = 101180, rhAir = 58, dt = 3600, insolation = 150, rhoLitter = 550,
                      nelsonA = 5.2, nelsonB = -19, conLitter = 0.0003, sigma = 3000, EPS = 0.01)
{
  SpecHum <- rhAir/100 *QSat(tAir, pAir) # Specific humidity
  xStart <- m
  xEnd <- m
  Ftest <- 0 - sigma * Ema(m, nelsonA, nelsonB, conLitter, tAir, vAir, pAir, insolation, rhAir)
  
  while(abs(Ftest) > 0.0001){
    Ftest <- rhoLitter * (xStart - xEnd) / dt
    Ftest <- Ftest - sigma * Ema(m, nelsonA, nelsonB, conLitter, tAir, vAir, pAir, insolation, rhAir)
    
    if(abs(Ftest) > 0.0001){
      temp <- xEnd
      h <- EPS * abs(temp)
      if(h == 0) {h <- EPS}
      xEnd <- temp + h
      m <- xEnd
      h <- xEnd - temp
      F2 <- rhoLitter * (xStart - xEnd) / dt - sigma * Ema(m, nelsonA, nelsonB, conLitter, tAir, vAir, pAir, insolation, SpecHum)
      xEnd <- temp
      m <- xEnd
      
      Ftest <- -Ftest * h / (F2 - Ftest)
      
      xEnd <- xEnd + Ftest
      m <- xEnd
    }
  }
  
  fmc <- m
}