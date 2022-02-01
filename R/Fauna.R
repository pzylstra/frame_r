#' Fire risk for an exposed arboreal animal
#'
#' Calculates the degree of injury or likelihood of mortality to an animal
#' caused by an approaching fire front
#'
#' Utilises the output tables from 'threat' and 'radiation', and adds to these
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
#' Mammal pelage is given a standardised emissivity of 0.86, based on:
#' McGowan, N. E., Scantlebury, D. M., Maule, A. G. & Marks, N. J.
#' Measuring the emissivity of mammal pelage. Quant. Infrared Thermogr. J. 6733, 1–9 (2018).
#'
#' Skin Cp is set to 3.5, averaged from Duck, F. A.
#' Physical properties of tissues: a comprehensive reference book. (Elsevier Science, 1990).
#'
#' Finds furDensity - the density of fur or feathers on the animal using standard values of 66 fibres per mm2,
#' fibre diameter of 0.01mm (10μm), and α-keratin density of 1300 kg.m-3. Fibre count averaged from
#' Liwanag, H. E. M., Berta, A., Costa, D. P., Abney, M. & Williams, T. M.
#' Morphological and thermal properties of mammalian insulation: the evolution of fur for aquatic living.
#' Biol. J. Linn. Soc. 106, 926–939 (2012)
#'
#' Finds thermal dose from all heat inputs, using the formula from
#' Ciesielski, M., Mochnacki, B. & Szopa, R. Numerical modeling of biological tissue heating.
#' Admissible thermal dose. Sci. Res. Inst. Math. Comput. Sci. 10, 11–20 (2011).
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height directly ove ground (m) at which the species is expected to shelter from a fire.
#' @param low The closest horizontal distance between the flame origin and the point (m)
#' @param high The furthest horizontal distance between the flame origin and the point (m)
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param RH The relative humidity (0-1)
#' @param Class Class of animal. Allowable values are "mammalia", "aves", "amphibia" and "reptilia"
#' @param Dimension The "Characteristic length" of the animal (m)
#' @param Area The surface area of the animal (m^2)
#' @param protection The thickness of fur, feather or scales covering the animal (m)
#' @param count The number of fibres per square mm
#' @param fibre The mean fibre diameter of hairs in mm
#' @param Specific_heat The specific heat of the fur, feather or scale material only (kJ/kg/deg C)
#' @param skinCp The specific heat of the animal skin (kJ/kg/C)
#' @param skinK Thermal conductivity of the animal skin (W/m/C)
#' @param objectTemp The body temperature of the animal (deg C)
#' @param Shape The approximate shape of the animal - either "Flat", "Sphere", or "Cylinder"
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

arboreal <- function(Surf, Plant, percentile = 0.5, Height = 1, low = 1, high = 50, var = 10, Pressure = 1013.25,
                     Altitude = 0, RH = 0.51, Class = "mammalia", Dimension = 0.1, Area = 0.2,
                     protection = 0.02, count = 66, fibre = 0.01, Specific_heat = 2.5, skinCp = 3.5,
                     skinK = 2, objectTemp = 38, Shape = "Flat",updateProgress = NULL)
{
  # Collect testing stats
  ROS <- mean(Surf$ros_kph)/3.6
  TIME <- round((high - low)/ROS)
  Horiz <- high
  dens <- (count*pi*(fibre/2)^2)
  Volume <- Area * protection
  R <- sqrt(Area/pi)
  cl <- ifelse(Class == "amphibia", 1, 0)
  tempR <- objectTemp
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude)%>%
    mutate(t = 1,
           furDensity = (1300*dens)+(1-dens)*Density,
           furCp = (Specific_heat*dens)+(1-dens)*cpAir,
           specificHeat = (cl*skinCp)+abs(cl-1)*furCp,
           Mass = Volume * furDensity,
           Re = (Plume_velocity*Density*Dimension)/viscosity,
           hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           hSphere = 2 + 0.6*Re^(0.5)*0.888,
           hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
           h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
           #Incoming
           qc = h * Area *(tempAir - tempR),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           # Thermal conductivity
           kAir = 0.00028683*(tempAir+273.15)^0.7919,
           kFur = kAir+0.004853*protection,
           k = (cl*skinK)+abs(cl-1)*kFur,
           # Fourier conduction through pelage
           fourier = pmin(Qi,(Area * k * (tempAir - tempR)) / protection),
           tempFourier = pmax(tempR, (0.001 * fourier / (Mass * specificHeat) + tempR)),
           # Thermal dose
           thermalDose = max((0.69 * fourier * R),0),
           fatalityLikelihood = pmax(0,pmin(100, 75.405 * log(thermalDose) - 518.12)),
           #Outgoing
           fourierO = pmin(0,(Area * k * (tempAir - tempFourier)) / protection),
           tempR = pmax(objectTemp, tempFourier + (0.001 * fourierO / (Mass * specificHeat))),
           qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempR+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR,
           # Conservation of energy
           deltaT = 0.001*Q / (Mass * specificHeat),
           tempConservation = pmax(objectTemp, (tempR + deltaT)))
  
  # Average the object temperature. These become starting values for each step of the loop
  #  SiteM <- mean(Ca$tempConservation)
  #  SiteMCond <- mean(Ca$tempR)
  #  TDM <- mean(Ca$thermalDose)
  SiteM <- quantile(Ca$tempConservation, percentile)
  SiteMCond <- quantile(Ca$tempR, percentile)
  TDM <- quantile(Ca$thermalDose, percentile)
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  pbar <-  txtProgressBar(max = TIME,style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      mutate(t = t,
             furDensity = (1300*dens)+(1-dens)*Density,
             furCp = (Specific_heat*dens)+(1-dens)*cpAir,
             specificHeat = (cl*skinCp)+abs(cl-1)*furCp,
             Mass = Volume * furDensity,
             tempConservation = SiteM,
             tempR = SiteMCond,
             thermalDose = TDM,
             Re = (Plume_velocity*Density*Dimension)/viscosity,
             hFlat = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             hSphere = 2 + 0.6*Re^(0.5)*0.888,
             hCylinder = 0.35 + 0.47*Re^(1/2)*0.837,
             h = ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)),
             #Incoming
             qc = h * Area *(tempAir - tempR),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             # Thermal conductivity
             kAir = 0.00028683*(tempAir+273.15)^0.7919,
             kFur = kAir+0.004853*protection,
             k = (cl*skinK)+abs(cl-1)*kFur,
             # Fourier conduction through pelage
             fourier = pmin(Qi,(Area * k * (tempAir - tempR)) / protection),
             tempFourier = pmax(tempR, (0.001 * fourier / (Mass * specificHeat) + tempR)),
             # Thermal dose
             thermalDose = thermalDose + 0.69 * pmax(fourier, 0) * R,
             fatalityLikelihood = pmax(0,pmin(100, 75.405 * log(thermalDose) - 518.12)),
             #Outgoing
             fourierO = pmin(0,(Area * k * (tempAir - tempFourier)) / protection),
             tempR = pmax(objectTemp, tempFourier + (0.001 * fourierO / (Mass * specificHeat))),
             qrO = pmin(0,Area*0.86*0.0000000000567*((tempAir+273.15)^4 - (tempR+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR,
             # Conservation of energy
             deltaT = 0.001*Q / (Mass * specificHeat),
             tempConservation = pmax(objectTemp, (tempConservation + deltaT)))
    Ca <- rbind(Ca, Cb)
    
    # Average the object temperature
    #    SiteM <- mean(Cb$tempConservation)
    #    SiteMCond <- mean(Cb$tempR)
    #    TDM <- mean(Cb$thermalDose)
    SiteM <- quantile(Cb$tempConservation, percentile)
    SiteMCond <- quantile(Cb$tempR, percentile)
    TDM <- quantile(Cb$thermalDose, percentile)
    setTxtProgressBar(pbar,t)
    ## progress bar
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", TIME - t)
      updateProgress(detail = text)
    }
    t = t + 1
    Horiz = Horiz - ROS
  }
  
  # Create and export table
  Ca <- Ca %>%
    mutate(plume_kph = Plume_velocity * 3.6) %>%
    select(t, repId, ros_kph, wind_kph, plume_kph, presAtm, tempR,
           tempAir, viscosity, pAlpha, Re, h,
           Mass, furDensity, specificHeat,
           flameTemp, epsilon, E, phi, att,
           qc, qr, qrO, qR, Q,
           deltaT, tempConservation,
           kAir, kFur, k, fourier, tempFourier, tempR, fourierO,
           thermalDose, fatalityLikelihood)
  write.csv(Ca, "arboreal.csv")
  return(Ca)
}

#####################################################################

#' Fire risk for an animal sheltering in a wooden hollow
#'
#' Calculates the likelihood of mortality to an animal
#' caused by an approaching fire front
#'
#' Utilises the output tables from 'threat' and 'radiation', and adds to these
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
#' Finds animal mortality within a hollow based on the maximum tolerable temperature for a given
#' vapour pressure deficit, based on data from Lawrence, G. E.
#' Ecology of vertebrate animals in relation to chaparral fire in the Sierra Nevada foothills.
#' Ecology 47, 278–291 (1966)
#'
#' Heat is transferred into the hollow using Fourier's Law
#'
#' Thermal conductivity of bark is modelled as per Martin, R. E.
#' Thermal properties of bark. For. Prod. J. 13, 419–426 (1963)
#'
#' Specific heat of bark is modelled using Kain, G., Barbu, M. C., Hinterreiter, S., Richter, K. & Petutschnigg, A.
#' Using bark as a heat insulation material. BioResources 8, 3718–3731 (2013)
#'
#' Thermal conductivity of wood is modelled using an approach from Kollmann, F. F. P. & Cote, W. A.
#' Principles of wood science and technology I. Solid wood. (Springer-Verlag, 1968)
#'
#' Evaporates water at 100 degrees C
#'
#' Specific heat of wood is derived from an established empirical relationship in Volbehr, B.
#' Swelling of wood fiber. PhD Thesis. (University of Kiel, 1896)
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height directly over ground (m) at which the species is expected to shelter from a fire.
#' @param woodDensity The density of wood in the tree or log housing the hollow (kg/m3)
#' @param barkDensity The density of bark in the tree or log housing the hollow (kg/m3)
#' @param wood The thickness of wood on the thinnest side of the hollow (m)
#' @param bark The thickness of bark on the thinnest side of the hollow (m)
#' @param comBark Temperature directly under the burning bark (C)
#' @param resBark Flame residence in the plant bark (s)
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the wood
#' @param bMoisture The proportion oven-dry weight of moisture in the bark
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail The number of seconds to continue modelling after all flames have extinguished
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param Dimension The "Characteristic length" of the hollow (m)
#' @param Area The surface area of the thinnest side of the hollow (m^2)
#' @param diameter depth of the litter layer (mm)
#' @param hollowTemp The starting temperature inside the hollow (deg C)
#' @param Shape The approximate shape of the hollow exterior - either "Flat", "Sphere", or "Cylinder"
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

hollow <- function(Surf, Plant, percentile = 0.5, Height = 1, woodDensity = 700, barkDensity = 500,
                   wood = 0.1, bark = 0.02, comBark = 700, resBark = 45, RH = 0.2, moisture = 0.2, bMoisture = 0.5, distance = 5, trail = 360, var = 10, Pressure = 1013.25,
                   Altitude = 0, Dimension = 0.3, Area = 0.03, diameter = 0.005, surfDecl = 10,
                   startTemp = 25, Shape = "Flat",updateProgress = NULL)
{
  
  # Post-front surface flame heating
  lengthSurface <- mean(Surf$lengthSurface)
  residence <- 0.871*diameter^1.875
  depth <- diameter/1000
  
  # Collect step distance, time, and total distance
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  Tb <- round(distance/ROS)
  TIME <- Ta + trail
  Horiz <- distance
  
  # Description of the protection
  if (bark > 0) {
    Material <- "bark"
    step <- bark/4
  } else {
    bMoisture <- moisture
    barkDensity <- woodDensity
    step <- 0.8* wood
    wood <- 0.2 * wood
  }
  mass <- step * barkDensity
  massW <- wood * woodDensity
  R <- sqrt(0.01/pi)
  startM <- moisture
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
    summarise_all(mean)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density*Dimension)/viscosity,
           h = 0.35 + 0.47*Re^(1/2)*0.837,
           #Incoming heat from surface
           pt = pmax(0, t-Tb),
           comBark = ifelse(pt <= resBark, comBark, 0),
           postS = frame:::bole(lengthSurface,residence, depth, h = Height,
                                surfDecl = 10, t = pt),
           tempS = ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir),
           qc = h * Area * (tempS - startTemp),
           att = frame:::tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           
           #STEP A ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterA = bMoisture*mass,
           # Energy removed by current water quantity
           drainA = ifelse(startTemp>99,
                           ifelse(bMoisture>0,mWaterA*2256400,0),0),
           #Thermal values - (cp: J/kg/deg, k: W/m/deg)
           cpA = frame:::cp(Material = Material, temp = startTemp, moist = bMoisture),
           kA = frame:::k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
           # Conduction from above and below, less latent heat of evaporation
           fAD = ((Area * kA * (tempS - startTemp)) / step),
           fAU = 0,
           fourierA = fAD + fAU - max(0, min((fAD + fAU), drainA)),
           tempA = (fourierA / (mass * cpA) + startTemp),
           # Change in proportion water this step
           moistureA = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((Qi/2256400)/mWaterA)),
                                                  bMoisture),bMoisture),
           
           #STEP B ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water (kg)
           mWaterB = bMoisture*mass,
           # Energy removed by current water quantity
           drainB = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterB*2256400,0),0),
           #Thermal values - (cp: J/kg/deg, k: W/m/deg)
           cpB = frame:::cp(Material = Material, temp = startTemp, moist = bMoisture),
           kB = frame:::k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
           # Conduction from above and below, less latent heat of evaporation
           fBD = ((Area * kB * (tempA - startTemp)) / step),
           fBU = 0,
           fourierB = fBD + fBU - max(0, min((fBD + fBU), drainB)),
           tempB = (fourierB / (mass * cpB) + startTemp),
           # Change in proportion water this step
           moistureB = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((fourierA/2256400)/mWaterB)),
                                                  bMoisture), bMoisture),
           
           #STEP C ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterC = bMoisture*mass,
           # Energy removed by current water quantity
           drainC = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterC*2256400,0),0),
           #Thermal values - (cp: J/kg/deg, k: W/m/deg)
           cpC = frame:::cp(Material = Material, temp = startTemp, moist = bMoisture),
           kC = frame:::k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
           # Conduction from above and below, less latent heat of evaporation
           fCD = ((Area * kC * (tempB - startTemp)) / step),
           fCU = 0,
           fourierC = fCD + fCU - max(0, min((fCD + fCU), drainC)),
           tempC = (fourierC / (mass * cpC) + startTemp),
           # Change in proportion water this step
           moistureC = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((fourierB/2256400)/mWaterC)),
                                                  bMoisture),bMoisture),
           
           #STEP D ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterD = bMoisture*mass,
           # Energy removed by current water quantity
           drainD = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterD*2256400,0),0),
           #Thermal values - (cp: J/kg/deg, k: W/m/deg)
           cpD = frame:::cp(Material = Material, temp = startTemp, moist = bMoisture),
           kD = frame:::k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
           # Conduction from above and below, less latent heat of evaporation
           fDD = ((Area * kD * (tempC - startTemp)) / step),
           fDU = 0,
           fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
           tempD = (fourierD / (mass * cpD) + startTemp),
           # Change in proportion water this step
           moistureD = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((fourierC/2256400)/mWaterD)),
                                                  bMoisture),bMoisture),
           
           # Wood ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterE = moisture*massW,
           # Energy removed by current water quantity
           drainE = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterE*2256400,0),0),
           #Thermal values - (cp: J/g/deg, k: W/m/deg)
           cpE = frame:::cp(temp = startTemp, moist = moisture),
           kE = frame:::k(temp = startTemp, moist = moisture, density = woodDensity),
           # Conduction from above and below, less latent heat of evaporation
           fED = ((Area * kE * (tempD - startTemp)) / step),
           fEU = 0,
           fourierE = fED + fEU - max(0, min((fED + fEU), drainE)),
           tempE = (fourierE / (massW * cpE) + startTemp),
           # Change in proportion water this step
           moistureE = ifelse(startTemp>99,ifelse(moisture>0,max(0,moisture-((fourierD/2256400)/mWaterE)),
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
  Horiz <- Horiz - ROS
  pbar <-  txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      summarise_all(mean)%>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density*Dimension)/viscosity,
             h = 0.35 + 0.47*Re^(1/2)*0.837,
             #Incoming heat from surface
             pt = pmax(0, t-Tb),
             comBark = ifelse(pt <= resBark, comBark, 0),
             postS = frame:::bole(lengthSurface,residence, depth, h = Height,
                                  surfDecl = 10, t = pt),
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir)),
             qc = h * Area * (tempS - tempA),
             att = frame:::tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             
             
             #STEP A ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterA = moistureA*mass,
             # Energy removed by current water quantity
             drainA = ifelse(tempA>99,
                             ifelse(moistureA>0,mWaterA*2256400,0),0),
             #Bark thermal values - (cp: J/kg/deg, k: W/m/deg)
             cpA = frame:::cp(Material = Material, temp = tempA, moist = moistureA),
             kA = frame:::k(Material = Material, temp = tempA, moist = moistureA, density = barkDensity),
             # Conduction from above and below, less latent heat of evaporation
             fAD = ((Area * kA * (tempS - tempA)) / step),
             fAU = ((Area * kB * (tempB - tempA)) / step),
             fourierA = fAD + fAU - max(0, min((fAD + fAU), drainA)),
             tempA = (fourierA / (mass * cpA) + tempA),  
             # Change in proportion water this step
             moistureA = ifelse(tempA>99,ifelse(moistureA>0,max(0,moistureA-((Qi/2256400)/mWaterA)),
                                                moistureA),moistureA),
             
             
             #STEP B ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterB = moistureB*mass,
             # Energy removed by current water quantity
             drainB = ifelse(tempB>99,
                             ifelse(moistureB>0,mWaterB*2256400,0),0),
             #Thermal values - (cp: J/kg/deg, k: W/m/deg)
             cpB = frame:::cp(Material = Material, temp = tempB, moist = moistureB),
             kB = frame:::k(Material = Material, temp = tempB, moist = moistureB, density = barkDensity),
             # Conduction from above and below, less latent heat of evaporation
             fBD = ((Area * kB * (tempA - tempB)) / step),
             fBU = ((Area * kC * (tempC - tempB)) / step),
             fourierB = fBD + fBU - max(0, min((fBD + fBU), drainB)),
             tempB = (fourierB / (mass * cpB) + tempB),
             # Change in proportion water this step
             moistureB = ifelse(tempB>99,ifelse(moistureB>0,max(0,moistureB-((fourierA/2256400)/mWaterB)),
                                                moistureB),moistureB),
             
             #STEP C ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterC = moistureC*mass,
             # Energy removed by current water quantity
             drainC = ifelse(tempC>99,
                             ifelse(moistureC>0,mWaterC*2256400,0),0),
             #Thermal values - (cp: J/kg/deg, k: W/m/deg)
             cpC = frame:::cp(Material = Material, temp = tempC, moist = moistureC),
             kC = frame:::k(Material = Material, temp = tempC, moist = moistureC, density = barkDensity),
             # Conduction from above and below, less latent heat of evaporation
             fCD = ((Area * kC * (tempB - tempC)) / step),
             fCU = ((Area * kC * (tempD - tempC)) / step),
             fourierC = fCD + fCU - max(0, min((fCD + fCU), drainC)),
             tempC = (fourierC / (mass * cpC) + tempC),
             # Change in proportion water this step
             moistureC = ifelse(tempC>99,ifelse(moistureC>0,max(0,moistureC-((fourierB/2256400)/mWaterC)),
                                                moistureC),moistureC),
             
             #STEP D ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterD = moistureD*mass,
             # Energy removed by current water quantity
             drainD = ifelse(tempD>99,
                             ifelse(moistureD>0,mWaterD*2256400,0),0),
             #Thermal values - (cp: J/kg/deg, k: W/m/deg)
             cpD = frame:::cp(Material = Material, temp = tempD, moist = moistureD),
             kD = frame:::k(Material = Material, temp = tempD, moist = moistureD, density = barkDensity),
             # Conduction from above and below, less latent heat of evaporation
             fDD = ((Area * kD * (tempC - tempD)) / step),
             fDU = ((Area * kD * (tempE - tempD)) / step),
             fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
             tempD = (fourierD / (mass * cpD) + tempD),
             # Change in proportion water this step
             moistureD = ifelse(tempD>99,ifelse(moistureD>0,max(0,moistureD-((fourierC/2256400)/mWaterD)),
                                                moistureD),moistureD),
             
             #wood ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterE = moistureE*massW,
             # Energy removed by current water quantity
             drainE = ifelse(tempE>99,
                             ifelse(moistureE>0,mWaterE*2256400,0),0),
             #Thermal values - (cp: J/g/deg, k: W/m/deg)
             cpE = frame:::cp(temp = tempE, moist = moistureE),
             kE = frame:::k(temp = tempE, moist = moistureE, density = barkDensity),
             # Conduction from above and below, less latent heat of evaporation. Below unknown.
             fED = ((Area * kE * (tempD - tempE)) / wood),
             fEU = 0,
             fourierE = fED + fEU - max(0, min((fED + fEU), drainE)),
             tempE = (fourierE / (massW * cpE) + tempE),
             # Change in proportion water this step
             moistureE = ifelse(tempE>99,ifelse(moistureE>0,max(0,moistureE-((fourierD/2256400)/mWaterE)),
                                                moistureE),moistureE))
    
    Ca <- rbind(Ca, Cb)
    
    #Collect values for the next step
    tempA <- Cb$tempA
    moistureA <- Cb$moistureA
    kA <- Cb$kA
    tempB <- Cb$tempB
    moistureB <- Cb$moistureB
    kB <- Cb$kB
    tempC <- Cb$tempC
    moistureC <- Cb$moistureC
    kC <- Cb$kC
    tempD <- Cb$tempD
    moistureD <- Cb$moistureD
    kD <- Cb$kD
    tempE <- Cb$tempE
    moistureE <- Cb$moistureE
    kE <- Cb$kE
    
    setTxtProgressBar(pbar,t)
    ##  progress bar
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0(TIME - t, "Remaining steps.")
      updateProgress(detail = text)
    }
    t <- t + 1
    Horiz <- Horiz - ROS
  }
  
  # Create table
  Ca <- Ca %>%
    select(t, tempS, tempA, tempB, tempC, tempD, tempE,
           moistureA, moistureB, moistureC, moistureD, moistureE) %>%
    mutate(mortality = ifelse(tempE < 67.5-0.3017*30.17*RH, 0, 1))
  
  return(Ca)
}

#####################################################################

#' Fire risk for an animal sheltering underground
#'
#' Calculates the likelihood of mortality to an animal
#' caused by an approaching fire front
#'
#' Utilises the output tables from 'threat' and 'radiation', and adds to these
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
#' Finds animal mortality within a hollow based on the maximum tolerable temperature for a given
#' vapour pressure deficit, based on data from Lawrence, G. E.
#' Ecology of vertebrate animals in relation to chaparral fire in the Sierra Nevada foothills.
#' Ecology 47, 278–291 (1966)
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
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param Plant The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param diameter Diameter of the surface fuels burning (mm)
#' @param surface Temperature at the surface of the soil, under burning fuels
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
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
#' @param depth The depth at which the animal shelters beneath the soil
#' @param soilTemp The starting temperature under the ground (deg C)
#' @param updateProgress Progress bar for use in the dashboard
#' @return dataframe
#' @export

underground <- function(Surf, Plant, diameter = 6, surface = 677, percentile = 0.5, RH = 0.2,
                        moisture = 0.2, distance = 50, trail = 300, var = 10, Pressure = 1013.25,
                        Altitude = 0, texture = "clay", peat = 0.1, grain = "fine", 
                        unfrozen = 1, depth = 0.1, soilTemp = 25,updateProgress = NULL)
  
{
  # Collect step distance, time, and total distance
  residence <- 0.871*diameter^1.875
  ROS <- mean(Surf$ros_kph)/3.6
  Ta <- round(distance/ROS+residence)
  TIME <- Ta + trail
  Horiz <- distance
  HA <- abs(Horiz)
  
  # Description of the protection
  densityD <- denSoil(texture)
  mass <- depth * densityD
  R <- sqrt(1/pi)
  soilTemp <- soilTemp
  
  #Starting values
  Ca <- threat(Surf, Plant, HA, Height=0, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
           #Incoming heat
           tempS = ifelse(Horiz <0, pmax(tempAir, surface), tempAir),
           qc = h * (tempS - soilTemp),
           att = tau(D=HA, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = qr*att,
           Qi = pmax(0, qc)+qr,
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*mass,
           # Energy removed by current water quantity
           drain = ifelse(soilTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           Qi = max(Qi-drain,0),
           #Thermal values
           cpSoil = cpSoil((soilTemp+273.15), texture, peat, moisture),
           saturation = satSoil(texture, moisture),
           kSoil = kSoil(texture, saturation, grain = "grain", unfrozen = unfrozen),
           # Fourier conduction
           fourier = ifelse(Horiz>0, pmin(Qi,(kSoil * pmax(0,tempS - soilTemp)) / depth),
                            (kSoil * pmax(0,tempS - soilTemp)) / depth),
           tempSoil = pmax(soilTemp, (fourier / (mass * cpSoil) + soilTemp)),
           # Change in proportion wood water this step
           moisture = ifelse(moisture>0,max(0,moisture-((fourier/2256400)/mWater)),
                             moisture),
           # Mortality
           mortality = ifelse(tempSoil < 67.5-30.17*RH, 0, 1),
           #Outgoing
           fourierO = pmin(0,(kSoil * (tempS - tempSoil)) / depth),
           soilTemp = tempSoil + (fourierO / (mass * cpSoil)),
           qrO = pmin(0,0.0000000000567*((tempS+273.15)^4 - (soilTemp+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR)
  
  soilTemp <- quantile(Ca$soilTemp, percentile)
  moisture <- quantile(Ca$moisture, percentile)
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  HA <- abs(Horiz)
  
  pbar <-  txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, HA, Height=0, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888),
             #Incoming
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <0, pmax(tempAir, surface), tempAir)),
             qc = h * (tempS - soilTemp),
             att = tau(D=HA, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = qr*att,
             Qi = pmax(0, qc)+qr,
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moisture*mass,
             # Energy removed by current water quantity
             drain = ifelse(soilTemp>95,
                            ifelse(moisture>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             Qi = max(Qi-drain,0),
             #Thermal values
             cpSoil = cpSoil((soilTemp+273.15), texture, peat, moisture),
             saturation = satSoil(texture, moisture),
             kSoil = kSoil(texture, saturation, grain = "grain", unfrozen = unfrozen),
             # Fourier conduction
             fourier = ifelse(Horiz>0, pmin(Qi,(kSoil * pmax(0,tempS - soilTemp)) / depth),
                              (kSoil * pmax(0,tempS - soilTemp)) / depth),
             tempSoil = pmax(soilTemp, (fourier / (mass * cpSoil) + soilTemp)),
             # Change in proportion wood water this step
             moisture = ifelse(moisture>0,max(0,moisture-((fourier/2256400)/mWater)),
                               moisture),
             # Mortality
             mortality = ifelse(tempSoil < 67.5-0.3017*30.17*RH, 0, 1),
             #Outgoing
             fourierO = pmin(0,(kSoil * (tempS - tempSoil)) / depth),
             soilTemp = tempSoil + (fourierO / (mass * cpSoil)),
             qrO = pmin(0,0.0000000000567*((tempS+273.15)^4 - (soilTemp+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR)
    Ca <- rbind(Ca, Cb)
    
    soilTemp <- quantile(Cb$soilTemp, percentile)
    moisture <- quantile(Cb$moisture, percentile)
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
    HA <- abs(Horiz)
  }
  
  # Create table
  Ca <- Ca %>%
    select(t, repId, ros_kph, tempAir, tempS, tempSoil, moisture,
           cpSoil, kSoil, att, qc, qr, qrO, qR, Q, fourier, fourierO, mortality)
  write.csv(Ca, "underground.csv")
  return(Ca)
}