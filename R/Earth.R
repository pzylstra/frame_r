#' Soil heating
#'
#' Calculates the dynamic heating of a soil at 1cm increments, to 5cm depth
#'
#' Assumes all to be A horizon, with constant, uncompacted density
#'
#' Utilises the output tables from 'threat' and adds to these
#' the Reynolds Number, heat transfer coefficients, Newton's convective energy transfer coefficient,
#' and the temperature of the object each second.
#'
#' Reynolds Number utilises a standard formulation (e.g. Gordon, N. T., McMahon, T. A. & Finlayson, B. L.
#' Stream hydrology: an introduction for ecologists. (Wiley, 1992))
#'
#' Convective heat transfer coefficients use the widely adopted formulations of
#' Williams, F. A. Urban and wildland fire phenomenology. Prog. Energy Combust. Sci. 8, 317–354 (1982),
#' and Drysdale, D. An introduction to fire dynamics. (John Wiley and Sons, 1985)
#' utilising a Prandtl number of 0.7.
#'
#' Heat is transferred into the earth using Fourier's Law. Spread continues for a period after the
#' passage of the fire front, equal to the duration of the surface flame, as determined using
#' Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
#' Int. J. Wildl. Fire 10, 137–143 (2001).
#'
#' Default temperature of the resident flame is the average of the surface maximums in
#' Cawson, J. G., Nyman, P., Smith, H. G., Lane, P. N. J. & Sheridan, G. J.
#' How soil temperatures during prescribed burning affect soil water repellency,
#' infiltration and erosion. Geoderma 278, 12–22 (2016).
#'
#' Heating area is set to 1m2, flat, with a characteristic length of 1m
#'
#' Broad germination and seed death temperatures are based on
#' Auld, T. D. & O’Connel, M. A. Predicting patterns of post‐fire germination in
#' 35 eastern Australian Fabaceae. Aust. J. Ecol. 16, 53–70 (1991).
#' 
#' Predicts death of fine roots at 60C, but does not yet include a time component.
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param step The increment of soil depth at which each calculation will be modelled (m)
#' @param diameter Diameter of the surface fuels burning (mm)
#' @param surface Temperature at the surface of the soil, under burning fuels
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the bark and wood
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail Number of seconds to continue modelling after the front has passed
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param texture Soil texture. Allowable values are: "sand", "loamy sand", "sandy loam", "sandy clay loam",
#' "sand clay", "loam", "clay loam", "silt loam", "clay", "silty clay", "silty clay loam", "silt"
#' @param peat Organic proportion of the soil
#' @param grain Allowable values are "fine" or "coarse"
#' @param unfrozen Proportion of soil unfrozen, between 0 and 1
#' @param soilTemp The starting temperature under the ground (deg C)
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export


soil <- function(Surf, Plant, step = 0.01, diameter = 6, surface = 677, RH = 0.2,
                 moisture = 0.1, distance = 50, trail = 600, var = 10, Pressure = 1013.25,
                 Altitude = 0, texture = "sand", peat = 0, grain = "fine", unfrozen = 1, 
                 startTemp = 25,updateProgress = NULL)
  
{
  # Collect step distance, time, and total distance
  residence <- 0.871*diameter^1.875
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  TIME <- Ta + trail
  Horiz <- distance
  
  # Description of the protection
  densityD <- denSoil(texture)
  mass <- step * densityD
  R <- sqrt(1/pi)
  startM <- moisture
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height=0, var, Pressure, Altitude) %>%
    summarise_all(mean)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           #Incoming heat from surface
           tempS = ifelse(Horiz <0, max(tempAir, surface), tempAir),
           qc = h * (tempS - startTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = max(0, qc)+qr,
           
           #A________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drainA = ifelse(startTemp>99,
                           ifelse(moisture>0,mWater*2256400,0),0),
           #Thermal values
           cpA = cpSoil((startTemp+273.15), texture, peat, moisture),
           saturationA = satSoil(texture, moisture),
           kA = kSoil(texture, saturationA),
           # Conduction from above and below, less latent heat of evaporation
           fAD = ((kA * (tempS - startTemp)) / step),
           fAU = 0,
           fourierA = fAD + fAU - max(0,min((fAD + fAU), drainA)),
           tempA = (fourierA / (mass * cpA) + startTemp),
           # Change in proportion water this step
           moistureA = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((Qi/2256400)/mWater)),
                                                  moisture), moisture),
           
           #B________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drainB = ifelse(startTemp>99,
                           ifelse(moisture>0,mWater*2256400,0),0),
           #Thermal values
           cpB = cpSoil((startTemp+273.15), texture, peat, moistureA),
           saturationB = satSoil(texture, moistureA),
           kB = kSoil(texture, saturationB),
           # Conduction from above and below, less latent heat of evaporation
           fBD = ((kB * (tempA - startTemp)) / step),
           fBU = 0,
           fourierB = fBD + fBU - max(0,min((fBD + fBU), drainB)),
           tempB = (fourierB / (mass * cpB) + startTemp),
           # Change in proportion water this step
           moistureB = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((fourierA/2256400)/mWater)),
                                                  moisture), moisture),
           
           #C________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drainC = ifelse(startTemp>99,
                           ifelse(moisture>0,mWater*2256400,0),0),
           #Thermal values
           cpC = cpSoil((startTemp+273.15), texture, peat, moistureB),
           saturationC = satSoil(texture, moistureB),
           kC = kSoil(texture, saturationC),
           # Conduction from above and below, less latent heat of evaporation
           fCD = ((kC * (tempB - startTemp)) / step),
           fCU = 0,
           fourierC = fCD + fCU - max(0,min((fCD + fCU), drainC)),
           tempC = (fourierC / (mass * cpC) + startTemp),
           # Change in proportion water this step
           moistureC = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((fourierB/2256400)/mWater)),
                                                  moisture),moisture),
           
           #D________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drainD = ifelse(startTemp>99,
                           ifelse(moisture>0,mWater*2256400,0),0),
           #Thermal values
           cpD = cpSoil((startTemp+273.15), texture, peat, moistureC),
           saturationD = satSoil(texture, moistureC),
           kD = kSoil(texture, saturationD),
           # Conduction from above and below, less latent heat of evaporation
           fDD = ((kD * (tempC - startTemp)) / step),
           fDU = 0,
           fourierD = fDD + fDU - max(0,min((fDD + fDU), drainD)),
           tempD = (fourierD / (mass * cpD) + startTemp),
           # Change in proportion water this step
           moistureD = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((fourierC/2256400)/mWater)),
                                                  moisture),moisture),
           
           #E________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drainE = ifelse(startTemp>99,
                           ifelse(moisture>0,mWater*2256400,0),0),
           #Thermal values
           cpE = cpSoil((startTemp+273.15), texture, peat, moistureD),
           saturationE = satSoil(texture, moistureD),
           kE = kSoil(texture, saturationE),
           # Conduction from above and below, less latent heat of evaporation
           fED = ((kE * (tempD - startTemp)) / step),
           fEU = 0,
           fourierE = fED + fEU - max(0,min((fED + fEU), drainE)),
           tempE = (fourierE / (mass * cpE) + startTemp),
           # Change in proportion water this step
           moistureE = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((fourierD/2256400)/mWater)),
                                                  moisture),moisture))
  #Collect values for the next step
  tempA <- Ca$tempA
  moistureA <- Ca$moistureA
  kA <- Ca$kA
  tempB <- Ca$tempB
  moistureB <- Ca$moistureB
  kB <- Ca$kB
  tempC <- Ca$tempC
  moistureC <- Ca$moistureC
  kC <- Ca$kC
  tempD <- Ca$tempD
  moistureD <- Ca$moistureD
  kD <- Ca$kD
  tempE <- Ca$tempE
  moistureE <- Ca$moistureE
  kE <- Ca$kE
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  pbar <- txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height=0, var, Pressure, Altitude) %>%
      summarise_all(mean) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             #Incoming heat from surface
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <0, max(tempAir, surface), tempAir)),
             qc = h * (tempS - tempA),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = max(0, qc)+qr,
             
             #A________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureA*mass,
             # Energy removed by current water quantity
             drainA = ifelse(tempA>99,
                             ifelse(moistureA>0,mWater*2256400,0),0),
             #Thermal values
             cpA = cpSoil((tempA+273.15), texture, peat, moistureA),
             saturationA = satSoil(texture, moistureA),
             kA = kSoil(texture, saturationA),
             # Conduction from above and below, less latent heat of evaporation
             fAD = ((kA * (tempS - tempA)) / step),
             fAU = ((kB * (tempB - tempA)) / step),
             fourierA = fAD + fAU - max(0,min((fAD + fAU), drainA)),
             tempA = (fourierA / (mass * cpA) + tempA),
             # Change in proportion water this step
             moistureA = ifelse(tempA>99,ifelse(moistureA>0,max(0,moistureA-((Qi/2256400)/mWater)),
                                                moistureA),moistureA),
             
             #B________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureB*mass,
             # Energy removed by current water quantity
             drainB = ifelse(tempB>99,
                             ifelse(moistureB>0,mWater*2256400,0),0),
             #Thermal values
             cpB = cpSoil((tempB+273.15), texture, peat, moistureB),
             saturationB = satSoil(texture, moistureB),
             kB = kSoil(texture, saturationB),
             # Conduction from above and below, less latent heat of evaporation
             fBD = ((kB * (tempA - tempB)) / step),
             fBU = ((kC * (tempC - tempB)) / step),
             fourierB = fBD + fBU - max(0,min((fBD + fBU), drainB)),
             tempB = (fourierB / (mass * cpB) + tempB),
             # Change in proportion water this step
             moistureB = ifelse(tempB>99,ifelse(moistureB>0,max(0,moistureB-((fourierA/2256400)/mWater)),
                                                moistureB),moistureB),
             
             #C________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureC*mass,
             # Energy removed by current water quantity
             drainC = ifelse(tempC>99,
                             ifelse(moistureC>0,mWater*2256400,0),0),
             #Thermal values
             cpC = cpSoil((tempC+273.15), texture, peat, moistureC),
             saturationC = satSoil(texture, moistureC),
             kC = kSoil(texture, saturationC),
             # Conduction from above and below, less latent heat of evaporation
             fCD = ((kC * (tempB - tempC)) / step),
             fCU = ((kC * (tempD - tempC)) / step),
             fourierC = fCD + fCU - max(0,min((fCD + fCU), drainC)),
             tempC = (fourierC / (mass * cpC) + tempC),
             # Change in proportion water this step
             moistureC = ifelse(tempC>99,ifelse(moistureC>0,max(0,moistureC-((fourierB/2256400)/mWater)),
                                                moistureC),moistureC),
             
             #D________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureD*mass,
             # Energy removed by current water quantity
             drainD = ifelse(tempD>99,
                             ifelse(moistureD>0,mWater*2256400,0),0),
             #Thermal values
             cpD = cpSoil((tempD+273.15), texture, peat, moistureD),
             saturationD = satSoil(texture, moistureD),
             kD = kSoil(texture, saturationD),
             # Conduction from above and below, less latent heat of evaporation
             fDD = ((kD * (tempC - tempD)) / step),
             fDU = ((kD * (tempE - tempD)) / step),
             fourierD = fDD + fDU - max(0,min((fDD + fDU), drainD)),
             tempD = (fourierD / (mass * cpD) + tempD),
             # Change in proportion water this step
             moistureD = ifelse(tempD>99,ifelse(moistureD>0,max(0,moistureD-((fourierC/2256400)/mWater)),
                                                moistureD),moistureD),
             
             #E________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureE*mass,
             # Energy removed by current water quantity
             drainE = ifelse(tempE>99,
                             ifelse(moistureE>0,mWater*2256400,0),0),
             #Thermal values
             cpE = cpSoil((tempE+273.15), texture, peat, moistureE),
             saturationE = satSoil(texture, moistureE),
             kE = kSoil(texture, saturationE),
             # Conduction from above and below, less latent heat of evaporation. Below unknown.
             fED = ((kE * (tempD - tempE)) / step),
             fEU = 0,
             fourierE = fED + fEU - max(0,min((fED + fEU), drainE)),
             tempE = (fourierE / (mass * cpE) + tempE),
             # Change in proportion water this step
             moistureE = ifelse(tempE>99,ifelse(moistureE>0,max(0,moistureE-((fourierD/2256400)/mWater)),
                                                moistureE),moistureE))
    
    Ca <- suppressMessages(rbind(Ca, Cb))
    
    #Collect values for the next step
    tempA <- Cb$tempA
    moistureA <- Cb$moistureA
    kA <- Cb$kA
    tempB <- Cb$tempB
    moistureB <- Cb$moistureB
    kB <- Cb$kB
    tempC <- Cb$tempC
    moistureC <- Cb$moistureC
    kC <- Cb$kB
    tempD <- Cb$tempD
    moistureD <- Cb$moistureD
    kD <- Cb$kB
    tempE <- Cb$tempE
    moistureE <- Cb$moistureE
    kE <- Cb$kE
    
    setTxtProgressBar(pbar,t)
    ##  progress bar
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", TIME - t)
      updateProgress(detail = text)
    }
    t = t + 1
    Horiz = Horiz - ROS
  }
  
  # Create table
  Ca <- Ca %>%
    select(t, repId, tempS, tempA, tempB, tempC, tempD, tempE,
           moistureA, moistureB, moistureC, moistureD, moistureE)%>%
    mutate(startM = startM,
           seedDa = ifelse(tempA>100, 1, 0),
           seedDb = ifelse(tempB>100, 1, 0),
           seedDc = ifelse(tempC>100, 1, 0),
           seedDd = ifelse(tempD>100, 1, 0),
           seedDe = ifelse(tempE>100, 1, 0),
           seedga = ifelse(tempA>60, 1, 0)-seedDa,
           seedgb = ifelse(tempB>60, 1, 0)-seedDb,
           seedgc = ifelse(tempC>60, 1, 0)-seedDc,
           seedgd = ifelse(tempD>60, 1, 0)-seedDd,
           seedge = ifelse(tempE>60, 1, 0)-seedDe,
           rootDa = ifelse(tempA>60, 1, 0),
           rootDb = ifelse(tempB>60, 1, 0),
           rootDc = ifelse(tempC>60, 1, 0),
           rootDd = ifelse(tempD>60, 1, 0),
           rootDe = ifelse(tempE>60, 1, 0),
           orgA = max(0,min(1,0.0038*tempA-0.7692)),
           orgB = max(0,min(1,0.0038*tempB-0.7692)),
           orgC = max(0,min(1,0.0038*tempC-0.7692)),
           orgD = max(0,min(1,0.0038*tempD-0.7692)),
           orgE = max(0,min(1,0.0038*tempE-0.7692)),
           repelA = ifelse(tempA>175, 1, 0),
           repelB = ifelse(tempB>175, 1, 0),
           repelC = ifelse(tempC>175, 1, 0),
           repelD = ifelse(tempD>175, 1, 0),
           repelE = ifelse(tempE>175, 1, 0)
    )
  
  
  return(Ca)
}

#####################################################################

# Soil density
#
# Finds soil density from texture
# 
# Porosity of soils taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
# Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).
# 
# Density equation is from Peters-Lidard, C. D., Blackburn, E., Liang, X. & Wood, E. F.
# The effect of soil thermal conductivity parameterization on surface energy fluxes and temperatures.
# J. Atmos. Sci. 55, 1209–1224 (1998).


denSoil <- function(texture="loam")
{
  #Dry density
  porosity <- ifelse(texture=="sand",0.437,
                     ifelse(texture=="loamy sand",0.437,
                            ifelse(texture=="sandy loam", 0.453,
                                   ifelse(texture=="sandy clay loam",0.398,
                                          ifelse(texture=="sand clay",0.43,
                                                 ifelse(texture=="loam",0.463,
                                                        ifelse(texture=="clay loam",0.464,
                                                               ifelse(texture=="silt loam",0.501,
                                                                      ifelse(texture=="clay",0.475,
                                                                             ifelse(texture=="silty clay", 0.479,
                                                                                    ifelse(texture=="silty clay loam",0.471,0.5)))))))))))
  return((1-porosity)*2700)
}

#####################################################################

# Thermal conductivity of soil
#
# Model drawn from Johansen, O.
# Thermal conductivity of soils. PhD Thesis (University of Trondheim, 1971)
# Modified by Farouki, O. Thermal properties of soils. (Trans Tech, 1986)
# 
# Porosity taken from Rawls, W. J., Brakensiek, D. L. & Saxton, K. E.
# Estimation of soil water properties. Transactions of the ASAE 25, 1316–1320 & 1328 (1982).

kSoil <- function(texture="loam", saturation=0.3, grain="fine", unfrozen=1)
{
  saturation <-max(saturation, 0.1)
  porosity <- ifelse(texture=="sand",0.437,
                     ifelse(texture=="loamy sand",0.437,
                            ifelse(texture=="sandy loam", 0.453,
                                   ifelse(texture=="sandy clay loam",0.398,
                                          ifelse(texture=="sand clay",0.43,
                                                 ifelse(texture=="loam",0.463,
                                                        ifelse(texture=="clay loam",0.464,
                                                               ifelse(texture=="silt loam",0.501,
                                                                      ifelse(texture=="clay",0.475,
                                                                             ifelse(texture=="silty clay", 0.479,
                                                                                    ifelse(texture=="silty clay loam",0.471,0.5)))))))))))
  quartz <- ifelse(texture=="sand",0.92,
                   ifelse(texture=="loamy sand",0.82,
                          ifelse(texture=="sandy loam", 0.6,
                                 ifelse(texture=="sandy clay loam",0.6,
                                        ifelse(texture=="sand clay",0.52,
                                               ifelse(texture=="loam",0.4,
                                                      ifelse(texture=="clay loam",0.35,
                                                             ifelse(texture=="silt loam",0.25,
                                                                    ifelse(texture=="clay",0.25,
                                                                           ifelse(texture=="silty clay", 0.1,
                                                                                  ifelse(texture=="silty clay loam",0.1,
                                                                                         ifelse(texture=="silt",0.1,0))))))))))))
  minerals <- ifelse(texture=="sand",2,
                     ifelse(texture=="loamy sand",2,
                            ifelse(texture=="sandy loam", 2,
                                   ifelse(texture=="sandy clay loam",2,
                                          ifelse(texture=="sand clay",2,
                                                 ifelse(texture=="loam",2,
                                                        ifelse(texture=="clay loam",2,
                                                               ifelse(texture=="silt loam",2,
                                                                      ifelse(texture=="clay",2,3)))))))))
  
  kersten <- ifelse(grain == "fine",log10(saturation)+1,0.7*log10(saturation)+1)
  
  #Dry density (kg/m3)
  densityD <- (1-porosity)*2700
  # Dry thermal conductivity
  kDry <- (0.137*densityD+64.7)/(2700-0.947*densityD)
  # Solids thermal conductivity
  kS <- 7.7^quartz*minerals^(1-quartz)
  # Saturated thermal conductivity
  kSat <- kS^(1-porosity)*2.2^(porosity-unfrozen)*0.57^unfrozen
  
  return(kersten*(kSat-kDry)+kDry)
}


#####################################################################

# Specific heat of soil
#
# Finds volumetric specific heat from the mineral, organic and water components of the soil
#
# Specific heats of soil components estimated from figures 108 & 111 in
# Farouki, O. Thermal properties of soils. (Trans Tech, 1981).
#
# Water specific heat 4185 J/kg.K


cpSoil <- function(temp = 300, texture="loam", peat = 0.2, moisture=0.3)
{
  cpClay <- 3.66*temp-188
  cpSand <- 2.305*temp+67
  cpSilt <- 5.58*temp-854
  cpPeat <- 5.025*temp-151
  
  clay <- ifelse(texture=="sand",0.05,
                 ifelse(texture=="loamy sand",0.05,
                        ifelse(texture=="sandy loam", 0.15,
                               ifelse(texture=="sandy clay loam",0.25,
                                      ifelse(texture=="sand clay",0.4,
                                             ifelse(texture=="loam",0.2,
                                                    ifelse(texture=="clay loam",0.35,
                                                           ifelse(texture=="silt loam",0.15,
                                                                  ifelse(texture=="clay",0.7,
                                                                         ifelse(texture=="silty clay", 0.5,
                                                                                ifelse(texture=="silty clay loam",0.35,
                                                                                       ifelse(texture=="silt",0.05,0))))))))))))
  sand <- ifelse(texture=="sand",0.9,
                 ifelse(texture=="loamy sand",0.85,
                        ifelse(texture=="sandy loam", 0.65,
                               ifelse(texture=="sandy clay loam",0.6,
                                      ifelse(texture=="sand clay",0.5,
                                             ifelse(texture=="loam",0.4,
                                                    ifelse(texture=="clay loam",0.3,
                                                           ifelse(texture=="silt loam",0.2,
                                                                  ifelse(texture=="clay",0.15,
                                                                         ifelse(texture=="silty clay", 0.05,
                                                                                ifelse(texture=="silty clay loam",0.1,
                                                                                       ifelse(texture=="silt",0.05,0))))))))))))
  silt <- ifelse(texture=="sand",0.05,
                 ifelse(texture=="loamy sand",0.1,
                        ifelse(texture=="sandy loam", 0.2,
                               ifelse(texture=="sandy clay loam",0.15,
                                      ifelse(texture=="sand clay",0.1,
                                             ifelse(texture=="loam",0.4,
                                                    ifelse(texture=="clay loam",0.35,
                                                           ifelse(texture=="silt loam",0.65,
                                                                  ifelse(texture=="clay",0.15,
                                                                         ifelse(texture=="silty clay", 0.45,
                                                                                ifelse(texture=="silty clay loam",0.55,
                                                                                       ifelse(texture=="silt",0.9,0))))))))))))
  
  cpMineral <- clay*cpClay + sand*cpSand + silt*cpSilt
  cpDry <- peat*cpPeat+(1-peat)*cpMineral
  
  return(moisture*4185+(1-moisture)*cpDry)
}

#####################################################################

# Soil saturation
#
# Finds saturation from ODW moisture and texture
#
# Field capacity of soils taken from Salter, P. J. & Williams, J. B.
# The influence of texture on the moisture characteristics of soil.
# V. Relationships between particle-size composition and moisturecontents
# at the upper and lower limits of available-water. J. Soil Sci. 20, 126–131 (1969).

satSoil <- function(texture="loam", moisture=0.3)
{
  #Field capacity
  fcWW <- ifelse(texture=="sand",0.14,
                 ifelse(texture=="loamy sand",0.18,
                        ifelse(texture=="sandy loam", 0.26,
                               ifelse(texture=="sandy clay loam",0.26,
                                      ifelse(texture=="sand clay",0.29,
                                             ifelse(texture=="loam",0.30,
                                                    ifelse(texture=="clay loam",0.34,
                                                           ifelse(texture=="silt loam",0.39,
                                                                  ifelse(texture=="clay",0.42,
                                                                         ifelse(texture=="silty clay", 0.47,
                                                                                ifelse(texture=="silty clay loam",0.43,
                                                                                       ifelse(texture=="silt",0.45,0.45))))))))))))
  fc <- 1/((1/fcWW)-1)
  
  return(min(1,moisture/fc))
}


