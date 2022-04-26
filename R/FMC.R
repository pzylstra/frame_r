
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
#  esat <- 611.2 * exp(17.67 * (tAir - 273.15) / (tAir - 29.66))
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
#' @param vAir Wind speed (m/s)
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


tidyWeather <- function(clim) {
  clim <- clim %>%
    mutate(tAM = tAM + 273.15,
           tPM = tPM + 273.15,
           tMin = tMin + 273.15,
           tMax = tMax + 273.15,
           wAM = wAM/3.6,
           wPM = wPM/3.6,
           cAM = cAM/8,
           cPM = cPM/8,
           MSLP = MSLP*100,
           specHumAM = (rhAM/100)*frame:::QSat(tAM, MSLP),
           specHumPM = (rhPM/100)*frame:::QSat(tPM, MSLP))
  
  # Create sequence
  Temp  <- clim[,c('tAM','tPM','tMin','tMax','Day')] %>% pivot_longer(c(tMin, tAM, tMax, tPM), names_to = "Time", values_to = "Temp")  %>%
    mutate(Time = case_when(Time == "tMin" ~ 3,
                            Time == "tAM" ~ 9,
                            Time == "tMax" ~ 12,
                            Time == "tPM" ~ 15))
  Temp$Hour <- (Temp$Day * 24 - 24) + Temp$Time
  RH  <- clim[,c('rhAM','rhPM','Day')] %>% pivot_longer(c(rhAM, rhPM), names_to = "Time", values_to = "RH") %>%
    mutate(Time = case_when(Time == "rhAM" ~ 9,
                            Time == "rhPM" ~ 15))
  specHum  <- clim[,c('specHumAM','specHumPM','Day')] %>% pivot_longer(c(specHumAM, specHumPM), names_to = "Time", values_to = "specHum") %>%
    mutate(Time = case_when(Time == "specHumAM" ~ 9,
                            Time == "specHumPM" ~ 15))
  Wind  <- clim[,c('wAM','wPM','Day')] %>% pivot_longer(c(wAM, wPM), names_to = "Time", values_to = "Wind") %>%
    mutate(Time = case_when(Time == "wAM" ~ 9,
                            Time == "wPM" ~ 15))
  cloud  <- clim[,c('cAM','cPM','Day')] %>% pivot_longer(c(cAM, cPM), names_to = "Time", values_to = "Cloud") %>%
    mutate(Time = case_when(Time == "cAM" ~ 9,
                            Time == "cPM" ~ 15))
  Rain  <- clim[,c('Rain', 'MSLP', 'Day')] %>% 
    mutate(Time = 9)
  out <- left_join(Temp, cloud, by = c("Day", "Time"))
  out <- left_join(out, RH, by = c("Day", "Time"))
  out <- left_join(out, specHum, by = c("Day", "Time"))
  out <- left_join(out, Wind, by = c("Day", "Time"))
  out <- left_join(out, Rain, by = c("Day", "Time"))
  out <- out %>%
    select(Day, Hour, Rain, MSLP, Cloud, Temp, RH, specHum, Wind)
  out <- out[order(out$Hour),]
  
  return(out)
}


#' Formats hourly weather data for frame, with DFMC
#' 
#' Function uses the 'Single differential equation model' of Matthews et al (2010)
#' 
#' @param clim A dataset with the fields:
#' tAM (9am temp, degC)
#' tPM (3pm temp, degC)
#' tMin (Minimum daily temp, degC)
#' tMax (Maximum daily temp, degC)
#' rhAM (9am Relative humidity, %)
#' rhPM (3pm Relative humidity, %)
#' wAM (9am wind, km/h)
#' wPM (3pm wind, km/h)
#' cAM (Morning cloud, oktas)
#' cPM (Afternoon cloud, oktas)
#' MSLP (Mean Sea Level Pressure, hPa)
#' Rain (Daily rainfall, mm)
#' @param m The starting moisture (proportion ODW)
#' @param nelsonA Constant from Nelson FMC model
#' @param nelsonB Constant from Nelson FMC model
#' @param tAir Atmospheric temperature (K)
#' @param vAir Wind speed (m/s)
#' @param pAir Atmospheric pressure (Pa)
#' @param rhAir Relative humidity (%)
#' @param sigma Surface area: volume ratio of litter particles (m2/m3)
#' @param rhoLitter Density of leaves (kg/m3)
#' @param conLitter Litter conductance (m/s)
#' @param dt Seconds per time step
#' @param EPS Sensitivity of iteration
#' @param insolation Solar energy at the soil surface (kW/m2)
#' @export

frameWeather <- function(clim, m = 15, LAI = 3, WRF = 3, rholitter = 550, litterW = 10,
                         lat = -31.89, altitude = 8, slope = 0, rangeDir = 0,  dt = 3600,
                         nelsonA = 5.2, nelsonB = -19, conLitter = 0.0003, sigma = 3000, EPS = 0.01) {
  climDay <- tidyWeather(clim)
  Rain <- climDay[,c('Hour', 'Rain')]
  
  # Interpolate
  Hour <- seq(from = 1, to = max(climDay$Hour), by = 1)
  
  TempI <- splines::interpSpline( climDay$Hour, climDay$Temp )
  MSLPI <- splines::interpSpline( climDay$Hour, climDay$MSLP, na.action = na.omit )
  CloudI <- splines::interpSpline( climDay$Hour, climDay$Cloud, na.action = na.omit )
  sRHI <- splines::interpSpline( climDay$Hour, climDay$specHum, na.action = na.omit )
  WindI <- splines::interpSpline( climDay$Hour, climDay$Wind, na.action = na.omit )
  
  Temp <- predict( TempI, Hour )$y
  MSLP <- predict( MSLPI, Hour )$y
  Cloud <- predict( CloudI, Hour )$y
  sRH <- predict( sRHI, Hour )$y
  Wind <- predict( WindI, Hour )$y
  
  # Combine data
  hemisphere <- if (lat<0) {
    -1
  } else {1}
  
  out <- as.data.frame(Hour)
  out$Temp <- Temp
  out$sRH <- sRH
  out$MSLP <- MSLP
  out$RH <- (sRH / frame:::QSat(Temp, MSLP))*100
  out$Wind <- Wind
  out$Cloud <- Cloud
  
  # Add rain
  out <- out %>%
    left_join(Rain, by = "Hour")
  out$Rain[is.na(out$Rain)] <- 0
  
  # Calculate inputs & moisture
  out <- out %>%
    mutate(RH = pmax(0, pmin(100, RH)),
           Wind = pmax(0, Wind),
           Cloud = pmax(0, pmin(1, Cloud)),
           Declination = hemisphere*23.45*cos((360*((Hour/24)-172)/365)*pi/180),
           hourAngle = abs(12-((Hour/24)-floor(Hour/24))*24)*15,
           solarAltitude = asin(cos(lat*pi/180)*cos(Declination*pi/180)*cos(hourAngle*pi/180)+sin(lat*pi/180)*sin(Declination*pi/180)),
           TerrainA = (cos((hourAngle-(rangeDir+90))*pi/180)*slope)*pi/180,
           TerrainB = (cos((hourAngle-(rangeDir-90))*pi/180)*slope)*pi/180,
           ZenithA = acos(sin(lat*pi/180)*sin(Declination*pi/180)+cos(lat*pi/180)*cos(Declination*pi/180)*cos(hourAngle*pi/180))-TerrainA,
           ZenithB = acos(sin(lat*pi/180)*sin(Declination*pi/180)+cos(lat*pi/180)*cos(Declination*pi/180)*cos(hourAngle*pi/180))-TerrainB,
           ShadeA = 1-exp(-(2/(pi*tan(pmax(0.01,1.570796-ZenithA))))*LAI),
           ShadeB = 1-exp(-(2/(pi*tan(pmax(0.01,1.570796-ZenithB))))*LAI),
           InsolationA = pmax(0,((1-ShadeA)*Cloud*1000*cos(ZenithA))),
           InsolationB = pmax(0,((1-ShadeB)*Cloud*1000*cos(ZenithB))),
           RainAdj = pmax(0,Rain-(0.001*(Rain/2)+0.08*LAI)),
           Wetting =((10*pmin(RainAdj,(0.01*RainAdj+0.36*(litterW/5))))/litterW))
  
  for (t in out$Hour) {
    if (t == 1) {
      mA <- m/100
      mB <- m/100
    } else {
      mA <- out$moistureA[t-1]
      mB <- out$moistureB[t-1]
    } 
    out$moistureA[t] <- pmax(0.01,pmin(mA,(frame::simplefmc(m = mA,
                                                            tAir = out$Temp[t],
                                                            vAir = out$Wind[t],
                                                            pAir = out$MSLP[t],
                                                            rhAir = out$RH[t],
                                                            insolation = out$InsolationA[t],
                                                            rhoLitter = rholitter,
                                                            dt = dt,
                                                            nelsonA = nelsonA,
                                                            nelsonB = nelsonB,
                                                            conLitter = conLitter,
                                                            sigma = sigma,
                                                            EPS = EPS)))) - max(0,0.5*(mA-1)) + out$Wetting[t]
    out$moistureB[t] <- pmax(0.01,pmin(mB,(frame::simplefmc(m = mB,
                                                            tAir = out$Temp[t],
                                                            vAir = out$Wind[t],
                                                            pAir = out$MSLP[t],
                                                            rhAir = out$RH[t],
                                                            insolation = out$InsolationB[t],
                                                            rhoLitter = rholitter,
                                                            dt = dt,
                                                            nelsonA = nelsonA,
                                                            nelsonB = nelsonB,
                                                            conLitter = conLitter,
                                                            sigma = sigma,
                                                            EPS = EPS)))) - max(0,0.5*(mB-1)) + out$Wetting[t]
  }
  out <- out %>%
    mutate(Temp = Temp - 273.15,
           Wind = Wind * 3.6,
           MSLP = MSLP/100,
           wetBulb = Temp*atan(0.151977*(RH+8.313659)^0.5)+atan(Temp+RH)-atan(RH-1.676331)+0.00391838*RH^(3/2)*atan(0.023101*RH)-4.686035,
           cgStrikes = 2.708*10^-46*exp(3.863*wetBulb))
  return(out)
}