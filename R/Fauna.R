#' Fire risk for an exposed arboreal mammal
#'
#' Calculates the degree of injury or likelihood of mortality 
#' to an exposed mammal caused by an approaching fire front
#'
#'
#' @param Surf The dataframe 'runs' exported from Monte Carlos as 'Summary.csv'
#' @param IP The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
#' @param Height The height directly over ground (m) at which the species is expected to shelter from a fire.
#' @param distance The starting horizontal distance between the flame origin and the point (m)
#' @param trail Number of seconds to continue modelling after the fire has passed
#' @param diameter Diameter of the surface fuels burning (mm)
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param RH Relative humidity (0-1)
#' @param bodyLength The "Characteristic length" of the animal (m)
#' @param surfaceArea The surface area of the animal (m^2)
#' @param bodyMass The mass of the animal (kg)
#' @param protection The thickness of fur covering the animal (m)
#' @param fibreCount The number of fibres per square mm
#' @param fibreDiameter The mean fibre diameter of hairs (mm)
#' @param fibreCp Specific heat of fibres (kJ/kg/C)
#' @param fiberSolid The proportion of the fibre (0-1)
#' @param skinCp Specific heat of the animal skin (kJ/kg/C)
#' @param skinDensity Density of the animal skin (kg/m3)
#' @param skinK Thermal conductivity of the animal skin (W/m/C)
#' @param bodyTemp The body temperature of the animal (deg C)
#' @param Shape The approximate shape of the animal - either "Flat", "Sphere", or "Cylinder"
#' @param xRate Horizontal speed of animal movement (m/s)
#' @param yRate Vertical speed of animal movement (m/s)
#' @param yMax Maximum height animal can reach (m)
#' @param surfDecl The slope of the surface (degrees)
#'
#' @return dataframe
#' @export

mammal <- function(Surf, IP, Height = 1, distance = 5, xRate = 0, trail = 360, diameter = 6, surfDecl = 10, var = 10, Pressure = 1013.25,
                   Altitude = 0, yRate = 0.65, yMax = 10, RH = 0.51, bodyLength = 0.1, surfaceArea = 0.2, bodyMass = 1,
                   protection = 0.0017, fibreCount = 100, fibreDiameter = 0.01, fibreCp = 2.5, fiberSolid = 0.5, skinCp = 3.5, skinDensity = 1020,
                   skinK = 0.187, bodyTemp = 37, Shape = "Cylinder")
{
  # Collect testing stats
  ROS <- mean(Surf$ros_kph)/3.6
  residence <- 0.871*diameter^1.875
  Ta <- round(distance/ROS+residence)
  Tb <- round(distance/ROS)
  TIME <- Ta + trail
  Horiz <- distance
  dens <- fiberSolid * fibreCount*pi*(fibreDiameter/2)^2
  Volume <- surfaceArea * protection
  R <- sqrt(surfaceArea/pi)
  tPelage <- bodyTemp
  tEpiderm <- bodyTemp
  tDermP <- bodyTemp
  tDermR <- bodyTemp
  skinK <- skinK / 1000 # Convert to kW.m.K
  epidermisT <- (10.01*bodyMass^0.143 + 47.7*bodyMass^0.202) / 1000000
  epiMass <- surfaceArea * epidermisT * skinDensity
  dermisT <- 756*bodyMass^0.187 / 1000000
  dermMass <- surfaceArea * dermisT * skinDensity
  
  # Starting values
  Ca <- threat(Surf, IP, Horiz, Height, var, Pressure, Altitude, residence, surfDecl)%>%
    summarise_all(mean)%>%
    mutate(t = 1,
           compression = ((1/(0.1705*Plume_velocity+1.0332))+0.00012*fibreCount-0.1593),
           furDensity = (1000*dens)+(1-dens)*Density,
           furCp = (fibreCp*dens)+(1-dens)*cpAir,
           pelMass = Volume * furDensity,
           Re = (Plume_velocity*Density*bodyLength)/viscosity,
           h = hFauna(Shape = Shape, Re = Re),
           #Incoming
           qc = h * surfaceArea *(tempAir - tPelage),
           att = tau(D = Horiz, flameTemp = flameTemp, temperature = (temperature+273.15), rh = RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           
           # PELAGE ________________________________________________________
           kAir = 0.00028683*(tempAir+273.15)^0.7919,
           kWind = 0.4349*Plume_velocity-0.016*Plume_velocity^2-0.4703*log(fibreCount)+3.63,
           kFur = kWind*(kAir+0.004853*protection)/1000, # Convert to kW/m.K
           fAD = ((kFur * (tempAir - tPelage)) / (protection*compression)),
           fAU = ((kFur * (tEpiderm - tPelage)) / (protection*compression)),
           fourierA = fAD + fAU,
           tPelage = 0.001 * fourierA / (pelMass * furCp) + tPelage,
           
           # EPIDERMIS ________________________________________________________
           fBD = ((skinK * (tPelage - tEpiderm)) / epidermisT),
           fBU = ((skinK * (tDermP - tEpiderm)) / epidermisT),
           fourierB = fBD + fBU,
           tEpiderm = max(bodyTemp, 0.001 * fourierB / (epiMass * skinCp) + tEpiderm),
           B1 = ifelse(tEpiderm >= 60, 1, 0),
           
           # PAPILLARY DERMIS ________________________________________________________
           fCD = ((skinK * (tEpiderm - tDermP)) / (0.2*dermisT)),
           fCU = ((skinK * (tDermR - tDermP)) / (0.2*dermisT)),
           fourierC = fCD + fCU,
           tDermP = max(bodyTemp, 0.001 * fourierC / ((0.2*dermMass) * skinCp) + tDermP),
           B2 = ifelse(tDermP >= 60, 1, 0),
           
           # RETICULAR DERMIS ________________________________________________________
           fDD = ((skinK * (tDermP - tDermR)) / (0.8*dermisT)),
           fDU = ((skinK * (bodyTemp - tDermR)) / (0.8*dermisT)),
           fourierD = fDD + fDU,
           tDermR = max(bodyTemp, 0.001 * fourierD / ((0.8*dermMass) * skinCp) + tDermR),
           B3 = ifelse(tDermR >= 60, 1, 0),
           
           # VP MORTALITY ________________________________________________________
           VPmortality = ifelse(tempAir < 67.5-0.3017*30.17*RH, 0, 1))
  
  # Preallocate list to store all time steps
  out_list <- vector("list", TIME)
  out_list[[1]] <- Ca
  
  # Collect values for the next step
  tPelage <- Ca$tPelage
  tEpiderm <- Ca$tEpiderm
  tDermP <- Ca$tDermP
  tDermR <- Ca$tDermR
  
  # Advance one second's travel
  Horiz <- Horiz - ROS + xRate
  Altitude <- min(Altitude + yRate, yMax)
  
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <- threat(Surf, IP, Horiz, Height, var, Pressure, Altitude, residence, surfDecl) %>%
      summarise_all(mean)%>%
      mutate(t = t,
             compression = ((1/(0.1705*Plume_velocity+1.0332))+0.00012*fibreCount-0.1593),
             furDensity = (1300*dens)+(1-dens)*Density,
             furCp = (fibreCp*dens)+(1-dens)*cpAir,
             pelMass = Volume * furDensity,
             Re = (Plume_velocity*Density*bodyLength)/viscosity,
             h = hFauna(Shape = Shape, Re = Re),
             #Incoming
             qc = h * surfaceArea *(tempAir - tPelage),
             att = tau(D = Horiz, flameTemp = flameTemp, temperature = (temperature+273.15), rh = RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             
             # PELAGE ________________________________________________________
             kAir = 0.00028683*(tempAir+273.15)^0.7919,
             kWind = 0.4349*Plume_velocity-0.016*Plume_velocity^2-0.4703*log(fibreCount)+3.63,
             kFur = kWind*(kAir+0.004853*protection)/1000, # Convert to kW/m.K
             fAD = ((kFur * (tempAir - tPelage)) / (protection*compression)),
             fAU = ((kFur * (tEpiderm - tPelage)) / (protection*compression)),
             fourierA = fAD + fAU,
             tPelage = 0.001 * fourierA / (pelMass * furCp) + tPelage,
             
             # EPIDERMIS ________________________________________________________
             fBD = ((skinK * (tPelage - tEpiderm)) / epidermisT),
             fBU = ((skinK * (tDermP - tEpiderm)) / epidermisT),
             fourierB = fBD + fBU,
             tEpiderm = max(bodyTemp, 0.001 * fourierB / (epiMass * skinCp) + tEpiderm),
             B1 = ifelse(tEpiderm >= 60, 1, 0),
             
             # PAPILLARY DERMIS ________________________________________________________
             fCD = ((skinK * (tEpiderm - tDermP)) / (0.2*dermisT)),
             fCU = ((skinK * (tDermR - tDermP)) / (0.2*dermisT)),
             fourierC = fCD + fCU,
             tDermP = max(bodyTemp, 0.001 * fourierC / ((0.2*dermMass) * skinCp) + tDermP),
             B2 = ifelse(tDermP >= 60, 1, 0),
             
             # RETICULAR DERMIS ________________________________________________________
             fDD = ((skinK * (tDermP - tDermR)) / (0.8*dermisT)),
             fDU = ((skinK * (bodyTemp - tDermR)) / (0.8*dermisT)),
             fourierD = fDD + fDU,
             tDermR = max(bodyTemp, 0.001 * fourierD / ((0.8*dermMass) * skinCp) + tDermR),
             B3 = ifelse(tDermR >= 60, 1, 0),
             
             # VP MORTALITY ________________________________________________________
             VPmortality = ifelse(tempAir < 67.5-0.3017*30.17*RH, 0, 1))
    
    # store this time step
    out_list[[t]] <- Cb
    
    # Stop when full-thickness burn or VP mortality occurs
    if (Cb$B3 == 1 || Cb$VPmortality == 1) {
      break
    }
    
    # Collect values for the next step
    tPelage <- Cb$tPelage
    tEpiderm <- Cb$tEpiderm
    tDermP <- Cb$tDermP
    tDermR <- Cb$tDermR
    
    Horiz <- Horiz - ROS + xRate
    Altitude <- min(Altitude + yRate, yMax)
  }
  
  # Bind all time steps into one data frame
  Ca <- dplyr::bind_rows(out_list)
  
  return(Ca)
}



#' Randomise tree hollows
#'
#'Uses values from Hofman et al () to randomise hollow parameters for Monte Carlo simulations
#'
#' @param wallThickness_min Minimum hollow wall thickness (m)
#' @param wallThickness_max Maximum hollow wall thickness (m)
#' @param lengthHollow_min Minimum hollow length (m)
#' @param lengthHollow_max Maximum hollow length (m)
#' @param diameterHollow_min Minimum hollow diameter (m)
#' @param diameterHollow_max Maximum hollow diameter (m)
#' @param hollow_tilt_deg Mean hollow tilt from horizontal (degrees)
#' @param HTs Standard deviation of hollow tilt (degrees)
#' @param HTmin Minimum hollow tilt (degrees)
#' @param tempDepression Mean temperature depression inside hollow (degrees C)
#' @param TDs Standard deviation of temperature depression (degrees C)
#' @param branchTrunkRat Ratio of hollows located in branches vs trunks (0-1)
#'
#' @returns Dataframe
#' @export
#'

randomiseHollow <- function(branchTrunkRat = 0.5,
                            wallThickness_min = 0.001, wallThickness_max = 0.13,
                            lengthHollow_min = 0.3, lengthHollow_max = 5,
                            diameterHollow_min = 0.08, diameterHollow_max = 0.35,
                            hollow_tilt_deg = 45, HTs = 10, HTmin = -10,
                            tempDepression = 3.9, TDs = 2.4) {
  
  locationHollow <- if(runif(1) < branchTrunkRat) "branch" else "trunk"
  
  tilt_angle <- if(locationHollow == "trunk") 90 else 
    round(extraDistr::rtnorm(1, mean = hollow_tilt_deg, sd = HTs, a = HTmin, b = Inf),0)
  
  out <- data.frame(
    locationHollow = locationHollow,
    wallThickness = round(runif(1, min = wallThickness_min, max = wallThickness_max),3),
    lengthHollow  = round(runif(1, min = lengthHollow_min, max = lengthHollow_max),2),
    diameterHollow= round(runif(1, min = diameterHollow_min, max = diameterHollow_max),2),
    hollow_az_deg = round(runif(1, min = 0, max = 360),0),
    hollow_tilt_deg = tilt_angle,
    tempDep = round(extraDistr::rtnorm(1, mean = tempDepression, sd = TDs),1)
  )
  return(out)
}


#####################################################################
# FUNCTIONS CALCULATING HOLLOW HEATING VIA AIR ENTRAINMENT


#' Calculates atmospheric pressure at a given altitude
#'
#' @param Altitude Altitude (m.a.s.l.)
#' @param Pressure Atmospheric pressure at sea level (Pa)
#' @param ambientTemp Ambient temperature (K)
#' @param lapseRate K/m (standard tropospheric lapse rate)
#' @param g Acceleration due to gravity (m/s^2)
#' @param Rgas J/(kg·K) for dry air
#'
#' @returns Value
#' 
siteAP <- function(Altitude,         
                   Pressure = 101325,  
                   ambientTemp = 281.65,  
                   lapseRate = 0.0065, 
                   g = 9.80665,      
                   Rgas = 287.058) { 
  
  #  theta <- 1 - lapseRate * h / T0_sea
  theta <- 1 - lapseRate * Altitude / (ambientTemp +(lapseRate * Altitude))
  Pressure * (theta^(g/(Rgas*lapseRate)))
}

#' Finds air density
#'
#' @param P Atmospheric pressure
#' @param T_K Temperature (K)
#' @param Rgas Specific gas constant for dry air
#'
#' @returns Value
#' 
air_density <- function(P, T_K, Rgas = 287.058) {
  P / (Rgas * T_K)
}


#' Finds entrainment depth at a given angle
#'
#' @param Plume_velocity Velocity of airflow (m/s) 
#' @param lengthHollow Length of the hollow (m) 
#' @param beta_deg Angle of the hollow to the airflow (degrees). 0 = facing the flow, 90 = perpendicular to the flow 
#' @param Cp_cross Pressure coefficient for suction in crossflow (typically 0.4–1.0) 
#' @param Cp_stag Pressure coefficient for stagnation at head-on flow 
#' @param Altitude Height above sea level (m) 
#' @param hollowTemp Starting temperature inside the hollow (C) 
#' @param tempAir Temperature of the convective plume heating the hollow (C) 
#' @param mode Considers expansion law for trapped gas
#' @param gamma Ratio of specific heats for air
#' @param Pressure Atmospheric pressure at sea level (Pa) 
#' @param ambientTemp Ambient temperature at the site (K) 
#' @param lapseRate K/m (standard tropospheric lapse rate) 
#' @param g Acceleration due to gravity (m/s^2) 
#' @param Rgas Specific gas constant for dry air J/(kg·K) 
#'
#' @return Value
#' 

entrainment_depth_angle <- function(Plume_velocity, lengthHollow,
                                    beta_deg = 90,
                                    Cp_cross = 0.9, Cp_stag = 1.0,
                                    Altitude = 0,
                                    hollowTemp = 20, tempAir = 100, #Air temp (C) inside and outside of the hollow
                                    mode = c("adiabatic","isothermal"),
                                    gamma = 1.4,
                                    Pressure = 101325, 
                                    ambientTemp = 281.65,
                                    lapseRate = 0.0065, g = 9.80665, Rgas = 287.058) {
  mode <- match.arg(mode)
  
  # Ambient pressure
  P_amb   <- siteAP(Altitude, Pressure, ambientTemp, lapseRate, g, Rgas)
  T_out_K <- tempAir + 273.15
  rho_out <- air_density(P_amb, T_out_K, Rgas)
  
  q    <- 0.5 * rho_out * Plume_velocity^2
  beta <- beta_deg * pi/180
  
  # Angle dependence: stagnation vs crossflow
  dP <- q * (Cp_stag * cos(beta)^2 - Cp_cross * sin(beta)^2)
  
  nfac <- if (mode == "adiabatic") gamma else 1.0
  x0   <- lengthHollow * (-dP) / (nfac * P_amb)   # suction (dP<0) => x0>0; stagnation => x0<0
  
  list(x0 = x0, dP = dP, P_amb = P_amb)
}

# Angle-dependent effective diffusivity ________________________________________
effective_diffusivity <- function(Plume_velocity, diameterHollow, beta_deg, D_mol = 2e-5, C_t = 0.02) {
  beta <- beta_deg * pi/180
  D_mol + C_t * Plume_velocity * diameterHollow * sin(beta)^2
}

mixing_fraction_over_time <- function(t, lengthHollow, x0, D_eff, k = 2) {
  delta <- k * sqrt(D_eff * pmax(t, 0))
  x_mix <- if (x0 > 0) pmin(lengthHollow, x0 + delta) else pmin(lengthHollow, delta)
  x_mix / lengthHollow
}



#' Finds the angle of the convective plume at a given height 
#'
#' @param base.params A parameter table used by the FRaME workflow
#' @param Plume_velocity Velocity of convective airflow (m/s)
#' @param Height Height above ground (m)
#'
#' @return value
#' @export
#'

plumeAngle <- function(base.params, Plume_velocity, Height) {
  x <- as.numeric(base.params$value[base.params$param == "windSpeed"])/windReduction(base.params, test = Height) / 3.6
  y <- Plume_velocity
  angle <- atan2(y, x) * 180 / pi
  angle
}


#' Calculate the angle between the wind direction and hollow axis
#'
#' Computes the relative angle (`beta_deg`) between a wind vector and a hollow axis,
#' accounting for both horizontal (azimuth) and vertical (elevation/tilt) angles
#' of each. A value of 0° means the wind flows directly *into* the hollow mouth,
#' while 90° indicates crossflow across the opening.
#'
#' @param wind_az_deg Numeric. Wind azimuth in degrees, measured clockwise from
#'   North (0° = North, 90° = East, 180° = South, 270° = West).
#' @param wind_elev_deg Numeric. Wind elevation angle in degrees above the horizontal
#'   (positive = upward, negative = downward). Default is 0 (horizontal flow).
#' @param hollow_az_deg Numeric. Hollow azimuth in degrees, measured clockwise from
#'   North. Defines the horizontal orientation of the hollow’s outward axis.
#' @param hollow_tilt_deg Numeric. Hollow tilt angle in degrees above the horizontal
#'   (positive = tilted upward, negative = downward). Default is 0 (horizontal).
#' @param wind_is_from Logical. If `TRUE` (default), the `wind_az_deg` value is the
#'   direction from which the wind blows (meteorological convention). If `FALSE`,
#'   it is interpreted as the direction towards which the wind is moving.
#'
#' @return Numeric scalar giving the relative angle between the wind vector and
#'   the hollow axis, in degrees. Values range from 0° (direct inflow) to 180°
#'   (direct outflow).
#'
#' @examples
#' # East-facing hollow, horizontal wind from West (into mouth)
#' betaHollow(wind_az_deg = 270, hollow_az_deg = 90)
#'
#' # Same hollow, wind from North (crossflow)
#' betaHollow(wind_az_deg = 0, hollow_az_deg = 90)
#'
#' # Hollow tilted 30° downward, wind along same direction and tilt
#' betaHollow(wind_az_deg = 270, wind_elev_deg = -30,
#'                     hollow_az_deg = 90, hollow_tilt_deg = -30)
#'
#' @export

betaHollow <- function(wind_az_deg, wind_elev_deg = 0, hollow_az_deg, hollow_tilt_deg = 0,
                       wind_is_from = TRUE) {
  deg2rad <- function(d) d * pi / 180
  clamp   <- function(x) pmax(-1, pmin(1, x))
  
  wind_az_towards <- if (wind_is_from) (wind_az_deg + 180) %% 360 else wind_az_deg
  
  unit_vec <- function(az_deg, el_deg) {
    az <- deg2rad(az_deg); el <- deg2rad(el_deg)
    c(x = cos(el) * sin(az),   # East
      y = cos(el) * cos(az),   # North
      z = sin(el))             # Up
  }
  
  w <- unit_vec(wind_az_towards, wind_elev_deg)
  p_out <- unit_vec(hollow_az_deg, hollow_tilt_deg)
  cos_beta <- clamp(sum(w * p_out))
  beta_rad <- acos(cos_beta)
  beta_deg <- beta_rad * 180 / pi
  as.numeric(beta_deg)
}


#' Heating of a hollow through air entrainment
#'
#' @param Plume_velocity Velocity of airflow (m/s)
#' @param lengthHollow Length of the hollow (m)
#' @param diameterHollow Diameter of the hollow (m)
#' @param beta_deg Angle of the hollow to the airflow (degrees). 0 = facing the flow, 90 = perpendicular to the flow
#' @param Cp_cross Pressure coefficient for suction in crossflow (typically 0.4–1.0)
#' @param Cp_stag Pressure coefficient for stagnation at head-on flow
#' @param Altitude Height above sea level (m)
#' @param hollowTemp Starting temperature inside the hollow (C)
#' @param tempAir Temperature of the convective plume heating the hollow (C)
#' @param mode Considers expansion law for trapped gas
#' @param times Time sequence to evaluate (s)
#' @param D_mol Molecular diffusivity of air (m^2/s)
#' @param C_t Coefficient for shear-enhanced turbulent diffusivity
#' @param k Coefficient for depth of mixing layer in the hollow
#' @param Pressure Atmospheric pressure at sea level (Pa)
#' @param ambientTemp Ambient temperature at the site (K)
#' @param lapseRate K/m (standard tropospheric lapse rate)
#' @param g Acceleration due to gravity (m/s^2)
#' @param Rgas Specific gas constant for dry air J/(kg·K)
#'
#' @returns value
#' @export

entrainmentHeating <- function(Plume_velocity, lengthHollow, diameterHollow,
                               Altitude = 0, hollowTemp = 30, tempAir = 100,
                               Pressure = 101325, ambientTemp = 281.65,
                               times = seq(0, 60, by = 1),
                               beta_deg = 90, Cp_cross = 0.9, Cp_stag = 1.0,
                               mode = c("adiabatic","isothermal"),
                               D_mol = 2e-5, C_t = 0.02, k = 2,
                               lapseRate = 0.0065, g = 9.80665, Rgas = 287.058) {
  
  mode <- match.arg(mode)
  ent <- entrainment_depth_angle(Plume_velocity, lengthHollow, beta_deg, Cp_cross, Cp_stag,
                                 Altitude, hollowTemp, tempAir, mode, gamma = 1.4,
                                 Pressure = Pressure, ambientTemp = ambientTemp,
                                 lapseRate = lapseRate, g = g, Rgas = Rgas)
  
  P_amb   <- ent$P_amb
  T_in0_K <- hollowTemp + 273.15
  T_out_K <- tempAir + 273.15
  rho_in  <- air_density(P_amb, T_in0_K, Rgas)
  rho_out <- air_density(P_amb, T_out_K, Rgas)
  
  # Angle-dependent diffusivity
  D_eff <- effective_diffusivity(Plume_velocity, diameterHollow, beta_deg, D_mol = D_mol, C_t = C_t)
  
  # Mixing & temperature evolution
  frac <- mixing_fraction_over_time(times, lengthHollow, ent$x0, D_eff = D_eff, k = k)
  T_in_time <- (1 - frac) * T_in0_K + frac * T_out_K
  
  out <- data.frame(
    t_s = times,
    beta_deg = beta_deg,
    fraction_mixed = frac,
    T_in_C = T_in_time - 273.15,
    x0_m = ent$x0,
    dP_Pa = ent$dP,
    D_eff = D_eff,
    P_amb_Pa = P_amb,
    rho_in = rho_in,
    rho_out = rho_out,
    Altitude = Altitude,
    Pressure = Pressure, ambientTemp = ambientTemp,
    lapseRate = lapseRate, g = g, Rgas = Rgas
  )
  out
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
#' @param IP The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
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
#' @param Shape The approximate shape of the hollow exterior - either "Flat", "Sphere", or "Cylinder"
#' @param surfDecl The slope of the surface (degrees)
#' @param base.params A parameter table used by the FRaME workflow
#' Numeric. Hollow azimuth in degrees, measured clockwise from
#'   North. Defines the horizontal orientation of the hollow’s outward axis.
#' @param hollow_tilt_deg Numeric. Hollow tilt angle in degrees above the horizontal
#'   (positive = tilted upward, negative = downward). Default is 0 (horizontal).
#' @param wind_is_from Logical. If `TRUE` (default), the `wind_az_deg` value is the
#'   direction from which the wind blows (meteorological convention). If `FALSE`,
#'   it is interpreted as the direction towards which the wind is moving.
#' @param wind_az_deg Numeric. Wind azimuth in degrees, measured clockwise from
#'   North (0° = North, 90° = East, 180° = South, 270° = West).
#' @param startTemp The initial temperature of the hollow (C)
#' @param hollow_az_deg Numeric. Hollow azimuth in degrees, measured clockwise from
#'  North (0° = North, 90° = East, 180° = South, 270° = West).
#' @param diameterHollow Diameter of the hollow (m)
#' @param mortalityEnd Stops loop if mortality occurs
#'
#' @return dataframe
#' @export

hollow <- function(Surf, IP, base.params, Height = 1, woodDensity = 700, barkDensity = 500, wood = 0.1, bark = 0.02,
                   comBark = 700, resBark = 45, RH = 0.5, moisture = 0.2, bMoisture = 0.5,
                   distance = 5, trail = 360, var = 10, Pressure = 1013.25, Altitude = 0,
                   Dimension = 0.3, diameterHollow = 0.2, Area = 0.03, diameter = 6, surfDecl = 2, startTemp = 21, Shape = "Cylinder", 
                   wind_is_from = TRUE, hollow_az_deg = 90, hollow_tilt_deg = 0, wind_az_deg = 270, mortalityEnd = TRUE
) {
  # --- repId to carry through results ---
  repId_val <- if ("repId" %in% names(Surf)) {
    u <- unique(Surf$repId)
    if (length(u) == 1) u else u[1]
  } else NA
  
  # --- Precompute constants (unchanged) ---
  lengthSurface <- mean(Surf$lengthSurface)
  residence     <- 0.871 * diameter^1.875
  depth         <- diameter / 1000
  ROS   <- mean(Surf$ros_kph) / 3.6
  Ta    <- round(distance/ROS + residence)
  Tb    <- round(distance/ROS)
  TIME  <- Ta + trail
  Horiz <- distance
  
  if (bark > 0) {
    Material <- "bark"; step_thk <- bark/4
  } else {
    Material <- "wood"; comBark <- 0; bMoisture <- moisture
    barkDensity <- woodDensity; step_thk <- 0.8 * wood; wood <- 0.2 * wood
  }
  mass   <- step_thk * barkDensity
  massW  <- wood * woodDensity
  
  keep_from_threat <- c(
    "tempAir","Plume_velocity","Density","viscosity","flameTemp",
    "temperature","epsilon","qr","lengthSurface","pAlpha"
  )
  
  # ---- t = 1
  Ca <- threat(Surf, IP, Horiz, Height, var, Pressure, Altitude, residence, surfDecl) %>%
    dplyr::select(dplyr::any_of(keep_from_threat)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
    dplyr::mutate(
      repId = repId_val,                
      t  = 1,
      Re = (Plume_velocity * Density * Dimension) / viscosity,
      h  = 0.35 + 0.47 * Re^(1/2) * 0.837,
      pt      = pmax(0, t - Tb),
      comBark = ifelse(pt <= resBark, comBark, 0),
      postS   = bole(lengthSurface, residence, depth, h = Height, surfDecl = 10, t = pt),
      tempS   = ifelse(Horiz <= 0, pmax(tempAir, postS, comBark), tempAir),
      qc      = h * Area * (tempS - startTemp),
      att     = tau(D = Horiz, flameTemp = flameTemp, temperature = (temperature + 273.15), rh = RH),
      qr      = 0.86 * qr * att,
      Qi      = pmax(0, qc) + qr,
      
      # STEP A
      mWaterA = bMoisture * mass,
      drainA  = ifelse(startTemp > 99, ifelse(bMoisture > 0, mWaterA * 2256400, 0), 0),
      cpA = cp(Material = Material, temp = startTemp, moist = bMoisture),
      kA  = k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
      fAD = ((Area * kA * (tempS - startTemp)) / step_thk),
      fAU = 0,
      fourierA  = fAD + fAU - max(0, min((fAD + fAU), drainA)),
      tempA     = (fourierA / (mass * cpA) + startTemp),
      moistureA = ifelse(startTemp > 99,
                         ifelse(bMoisture > 0, max(0, bMoisture - ((Qi/2256400)/mWaterA)), bMoisture),
                         bMoisture),
      
      # STEP B
      mWaterB = bMoisture * mass,
      drainB  = ifelse(startTemp > 99, ifelse(moisture > 0, mWaterB * 2256400, 0), 0),
      cpB = cp(Material = Material, temp = startTemp, moist = bMoisture),
      kB  = k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
      fBD = ((Area * kB * (tempA - startTemp)) / step_thk),
      fBU = 0,
      fourierB  = fBD + fBU - max(0, min((fBD + fBU), drainB)),
      tempB     = (fourierB / (mass * cpB) + startTemp),
      moistureB = ifelse(startTemp > 99,
                         ifelse(bMoisture > 0, max(0, bMoisture - ((fourierA/2256400)/mWaterB)), bMoisture),
                         bMoisture),
      
      # STEP C
      mWaterC = bMoisture * mass,
      drainC  = ifelse(startTemp > 99, ifelse(moisture > 0, mWaterC * 2256400, 0), 0),
      cpC = cp(Material = Material, temp = startTemp, moist = bMoisture),
      kC  = k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
      fCD = ((Area * kC * (tempB - startTemp)) / step_thk),
      fCU = 0,
      fourierC  = fCD + fCU - max(0, min((fCD + fCU), drainC)),
      tempC     = (fourierC / (mass * cpC) + startTemp),
      moistureC = ifelse(startTemp > 99,
                         ifelse(bMoisture > 0, max(0, bMoisture - ((fourierB/2256400)/mWaterC)), bMoisture),
                         bMoisture),
      
      # STEP D
      mWaterD = bMoisture * mass,
      drainD  = ifelse(startTemp > 99, ifelse(moisture > 0, mWaterD * 2256400, 0), 0),
      cpD = cp(Material = Material, temp = startTemp, moist = bMoisture),
      kD  = k(Material = Material, temp = startTemp, moist = bMoisture, density = barkDensity),
      fDD = ((Area * kD * (tempC - startTemp)) / step_thk),
      fDU = 0,
      fourierD  = fDD + fDU - max(0, min((fDD + fDU), drainD)),
      tempD     = (fourierD / (mass * cpD) + startTemp),
      moistureD = ifelse(startTemp > 99,
                         ifelse(bMoisture > 0, max(0, bMoisture - ((fourierC/2256400)/mWaterD)), bMoisture),
                         bMoisture),
      
      # WOOD (E)
      mWaterE = moisture * massW,
      drainE  = ifelse(startTemp > 99, ifelse(moisture > 0, mWaterE * 2256400, 0), 0),
      cpE = cp(temp = startTemp, moist = moisture),
      kE  = k(temp = startTemp, moist = moisture, density = woodDensity),
      fED = ((Area * kE * (tempD - startTemp)) / step_thk),
      fEU = 0,
      fourierE  = fED + fEU - max(0, min((fED + fEU), drainE)),
      tempE     = (fourierE / (massW * cpE) + startTemp),
      moistureE = ifelse(startTemp > 99,
                         ifelse(moisture > 0, max(0, moisture - ((fourierD/2256400)/mWaterE)), moisture),
                         moisture)
    )
  
  # Find internal heating by air entrainment
  Plume_velocity <- mean(Ca$Plume_velocity) 
  tempAir <- mean(Ca$tempAir)
  wind_elev_deg <- plumeAngle(base.params, Plume_velocity, Height) 
  Ca$wind_elev_deg <- wind_elev_deg
  
  # Until plume depth is modelled, intersection is taken from threat()
  b <- betaHollow(wind_az_deg, wind_elev_deg, hollow_az_deg, hollow_tilt_deg, wind_is_from)
  Ca$plume_hollow_deg <- b
  Ca$tempEntrainment <- entrainmentHeating(Plume_velocity, lengthHollow = Dimension, diameterHollow = diameterHollow, 
                                           Altitude = Altitude, beta_deg = b, Pressure = Pressure * 100, 
                                           ambientTemp = Surf$temperature[1] + 273.15, hollowTemp = Ca$tempE, tempAir = tempAir,
                                           times = 1)$T_in_C
    # Test for vapour pressure mortality
  Ca$VPmortality = ifelse(max(Ca$tempE, Ca$tempEntrainment) < 67.5 - 0.3017 * 30.17 * RH, 0, 1)
  
  # seed for loop
  tempA <- Ca$tempA; moistureA <- Ca$moistureA; kA <- Ca$kA
  tempB <- Ca$tempB; moistureB <- Ca$moistureB; kB <- Ca$kB
  tempC <- Ca$tempC; moistureC <- Ca$moistureC; kC <- Ca$kC
  tempD <- Ca$tempD; moistureD <- Ca$moistureD; kD <- Ca$kD
  tempE <- max(Ca$tempE, Ca$tempEntrainment); moistureE <- Ca$moistureE; kE <- Ca$kE
  Horiz <- Horiz - ROS
  
  out <- vector("list", TIME)
  out[[1]] <- Ca
  
  # ---- loop
  for (t in 2:TIME) {
    Cb <- threat(Surf, IP, Horiz, Height, var, Pressure, Altitude, residence, surfDecl) %>%
      dplyr::select(dplyr::any_of(keep_from_threat)) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
      dplyr::mutate(
        repId = repId_val,
        t  = t,
        Re = (Plume_velocity * Density * Dimension) / viscosity,
        h  = 0.35 + 0.47 * Re^(1/2) * 0.837,
        pt      = pmax(0, t - Tb),
        comBark = ifelse(pt <= resBark, comBark, 0),
        postS   = bole(lengthSurface, residence, depth, h = Height, surfDecl = 10, t = pt),
        tempS   = ifelse(t > Ta, tempAir, ifelse(Horiz <= 0, pmax(tempAir, postS, comBark), tempAir)),
        qc      = h * Area * (tempS - tempA),
        att     = tau(D = Horiz, flameTemp = flameTemp, temperature = (temperature + 273.15), rh = RH),
        qr      = 0.86 * qr * att,
        Qi      = pmax(0, qc) + qr,
        
        # A
        mWaterA = moistureA * mass,
        drainA  = ifelse(tempA > 99, ifelse(moistureA > 0, mWaterA * 2256400, 0), 0),
        cpA = cp(Material = Material, temp = tempA, moist = moistureA),
        kA  = k(Material = Material, temp = tempA, moist = moistureA, density = barkDensity),
        fAD = ((Area * kA * (tempS - tempA)) / step_thk),
        fAU = ((Area * kB * (tempB - tempA)) / step_thk),
        fourierA  = fAD + fAU - max(0, min((fAD + fAU), drainA)),
        tempA     = (fourierA / (mass * cpA) + tempA),
        moistureA = ifelse(tempA > 99,
                           ifelse(moistureA > 0, max(0, moistureA - ((Qi/2256400)/mWaterA)), moistureA),
                           moistureA),
        
        # B
        mWaterB = moistureB * mass,
        drainB  = ifelse(tempB > 99, ifelse(moistureB > 0, mWaterB * 2256400, 0), 0),
        cpB = cp(Material = Material, temp = tempB, moist = moistureB),
        kB  = k(Material = Material, temp = tempB, moist = moistureB, density = barkDensity),
        fBD = ((Area * kB * (tempA - tempB)) / step_thk),
        fBU = ((Area * kC * (tempC - tempB)) / step_thk),
        fourierB  = fBD + fBU - max(0, min((fBD + fBU), drainB)),
        tempB     = (fourierB / (mass * cpB) + tempB),
        moistureB = ifelse(tempB > 99,
                           ifelse(moistureB > 0, max(0, moistureB - ((fourierA/2256400)/mWaterB)), moistureB),
                           moistureB),
        
        # C
        mWaterC = moistureC * mass,
        drainC  = ifelse(tempC > 99, ifelse(moistureC > 0, mWaterC * 2256400, 0), 0),
        cpC = cp(Material = Material, temp = tempC, moist = moistureC),
        kC  = k(Material = Material, temp = tempC, moist = moistureC, density = barkDensity),
        fCD = ((Area * kC * (tempB - tempC)) / step_thk),
        fCU = ((Area * kC * (tempD - tempC)) / step_thk),
        fourierC  = fCD + fCU - max(0, min((fCD + fCU), drainC)),
        tempC     = (fourierC / (mass * cpC) + tempC),
        moistureC = ifelse(tempC > 99,
                           ifelse(moistureC > 0, max(0, moistureC - ((fourierB/2256400)/mWaterC)), moistureC),
                           moistureC),
        
        # D
        mWaterD = moistureD * mass,
        drainD  = ifelse(tempD > 99, ifelse(moistureD > 0, mWaterD * 2256400, 0), 0),
        cpD = cp(Material = Material, temp = tempD, moist = moistureD),
        kD  = k(Material = Material, temp = tempD, moist = moistureD, density = barkDensity),
        fDD = ((Area * kD * (tempC - tempD)) / step_thk),
        fDU = ((Area * kD * (tempE - tempD)) / step_thk),
        fourierD  = fDD + fDU - max(0, min((fDD + fDU), drainD)),
        tempD     = (fourierD / (mass * cpD) + tempD),
        moistureD = ifelse(tempD > 99,
                           ifelse(moistureD > 0, max(0, moistureD - ((fourierC/2256400)/mWaterD)), moistureD),
                           moistureD),
        
        # WOOD
        mWaterE = moistureE * massW,
        drainE  = ifelse(tempE > 99, ifelse(moistureE > 0, mWaterE * 2256400, 0), 0),
        cpE = cp(temp = tempE, moist = moistureE),
        kE  = k(temp = tempE, moist = moistureE, density = woodDensity),
        fED = ((Area * kE * (tempD - tempE)) / wood),
        fEU = 0,
        fourierE  = fED + fEU - max(0, min((fED + fEU), drainE)),
        tempE     = (fourierE / (massW * cpE) + tempE),
        moistureE = ifelse(tempE > 99,
                           ifelse(moistureE > 0, max(0, moistureE - ((fourierD/2256400)/mWaterE)), moistureE),
                           moistureE),
        wind_elev_deg = wind_elev_deg,
        plume_hollow_deg = b,
        
        # Find internal heating by air entrainment
        tempEntrainment = entrainmentHeating(Plume_velocity, lengthHollow = Dimension, diameterHollow = diameterHollow, 
                                             Altitude = Altitude, beta_deg = b, Pressure = Pressure * 100, 
                                             ambientTemp = Surf$temperature[1] + 273.15, hollowTemp = tempE, tempAir = tempAir,
                                             times = 1)$T_in_C,
        
        # Test for vapour pressure mortality
        VPmortality = ifelse(max(tempE, tempEntrainment) < 67.5 - 0.3017 * 30.17 * RH, 0, 1)
      )
    
    out[[t]] <- Cb
    
    if (mortalityEnd && Cb$VPmortality == 1) {
      break
    }
    
    
    # update state
    tempA <- Cb$tempA; moistureA <- Cb$moistureA; kA <- Cb$kA
    tempB <- Cb$tempB; moistureB <- Cb$moistureB; kB <- Cb$kB
    tempC <- Cb$tempC; moistureC <- Cb$moistureC; kC <- Cb$kC
    tempD <- Cb$tempD; moistureD <- Cb$moistureD; kD <- Cb$kD
    tempE <- max(Cb$tempE, Cb$tempEntrainment); 
    moistureE <- Cb$moistureE; kE <- Cb$kE
    Horiz <- Horiz - ROS
    
  }
  
  Ca <- dplyr::bind_rows(out)
  
  Ca
}

#____________________________________________________________________

#' Calculate probabilistic first-order fire impacts on wildlife
#'
#' @param Conditions Table of Weather and terrain conditions to be tested
#' @param Forests List of Structure and Flora tables for the forest stands to be tested
#' @param Trees Table of tree species in the analyses, listing physiological parameters
#' @param fLine Length of the fireline (m)
#' @param edgeDistance Distance from the edge of the fireground (m)
#' @param Altitude Altitude above sea level (m)
#' @param wind_az_deg Numeric. Wind azimuth in degrees, measured clockwise from
#'   North (0° = North, 90° = East, 180° = South, 270° = West).
#' @param hollowHprop Proportion of the tree height where hollows occur
#' @param stags Proportion of trees with hollow centres
#' @param wallThickness_min Minimum hollow wall thickness (m)
#' @param wallThickness_max Maximum hollow wall thickness (m)
#' @param lengthHollow_min Minimum hollow length (m)
#' @param lengthHollow_max Maximum hollow length (m)
#' @param diameterHollow_min Minimum hollow diameter (m)
#' @param diameterHollow_max Maximum hollow diameter (m)
#' @param hollow_tilt_deg Numeric. Hollow tilt angle in degrees above the horizontal
#'   (positive = tilted upward, negative = downward). Default is 0 (horizontal).
#' @param HTs Standard deviation of hollow tilt (degrees)
#' @param HTmin Minimum hollow tilt (degrees)
#' @param tempDepression Mean temperature depression inside hollow (degrees C)
#' @param TDs Standard deviation of temperature depression (degrees C)
#' @param climbProp Proportion of the tree height to which an animal can climb
#' @param glideAngle Angle of glide descent (degrees)
#' @param climbingSpeed Speed that the animal can climb (m/s)
#' @param xRate Horizontal speed of animal movement (m/s)
#' @param yRate Vertical speed of animal movement (m/s)
#' @param yMax Maximum height animal can reach (m)
#' @param shelter Proportion of animals likely to remain in hollows for shelter
#' @param glide Proportion of animals that leave hollows that are likely to glide to an adjacent tree
#' @param bodyLength Length of the animal body (m)
#' @param surfaceArea Surface area of the animal body (m2)
#' @param bodyMass Mass of the animal body (kg)
#' @param fibreLength The length of fur fibres covering the animal (m)
#' @param fibreCount The number of fibres per square mm
#' @param fibreDiameter The mean fibre diameter of hairs (mm)
#' @param fibreCp Specific heat of fibres (kJ/kg/C)
#' @param fiberSolid The proportion of the fibre (0-1)
#' @param skinCp Specific heat of the animal skin (kJ/kg/C)
#' @param skinDensity Density of the animal skin (kg/m3)
#' @param skinK Thermal conductivity of the animal skin (W/m/C)
#' @param bodyTemp The body temperature of the animal (deg C)
#' @param diameter Max diameter of surface litter particles (mm)
#' @param surfDecl Exponent describing the rate of post-front flame decay in surface litter fires
#' @param var The angle in degrees that the plume spreads above/below a central vector;defaults to 10
#' @param reps Number of probabilistic replicates per condition
#' @param distance Starting distance between the animal and the fire front (m)
#' @param trail Number of seconds to continue modelling after all flames have extinguished
#' @param freeCores Number of CPU cores to leave free for other processes
#' @param minHeight Minimum allowable height of hollow bearing trees (m)
#' @param deadHollow Proportion of hollows that are dead
#' @param keepCanopy Logical. If TRUE, randomised sites are only chosen where the tallest stratum is present
#' @param testN Number or tries to attempt in canopyCheck
#' @param branchTrunkRat Ratio of hollows located in branches vs trunks (0-1)
#'
#' @returns List of dataframes: 
#' 1) Runs - detailed results for each replicate; 
#' 2) IP - summary of ignition paths in each plant; 
#' 3) Impact - summary of mortalities associated with each run
#' @export
#'
#' 
frameWildlife <- function(Conditions, Forests, Trees, # Datasets
                          fLine = 100, edgeDistance = 100, Altitude = 550, wind_az_deg = 270,  # Fire conditions
                          branchTrunkRat = 0.5, minHeight = 15, hollowHprop = 0.65, stags = 0.36, # Hollow specifics (0.72 had basal scars; assume half of these are hollow)
                          deadHollow = 0.074, wallThickness_min = 0.001, wallThickness_max = 0.13, lengthHollow_min = 0.3, lengthHollow_max = 5,
                          diameterHollow_min = 0.08, diameterHollow_max = 0.35, hollow_tilt_deg = 45, HTs = 10, HTmin = -10,
                          tempDepression = 3.9, TDs = 2.4, keepCanopy = TRUE,                 
                          climbProp = 0.9, glideAngle = 45, climbingSpeed = 0.65,              # Animal behaviour
                          xRate = 0, yRate = 0.65, yMax = 10, shelter = 0.9, glide = 0.5,      # Strategy
                          bodyLength = 0.4, surfaceArea = 0.2, bodyMass = 1.4,                 # Animal physiology
                          fibreLength = 0.063, fibreCount = 20, fibreDiameter = 0.017, fibreCp = 2.5, fiberSolid = 0.5, 
                          skinCp = 3.5, skinDensity = 1020, skinK = 0.187, bodyTemp = 37,
                          diameter = 6, surfDecl = 10, var = 10,                               # Surface fire
                          reps = reps, distance = 10, trail = 600, testN = 100, freeCores = 2)
{
  nCores <- max(parallel::detectCores() - freeCores, 1)
  cat("Setting up parallel cluster with", nCores, "cores\n")
  cl <- parallel::makeCluster(nCores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # 2) Load packages on each core
  parallel::clusterEvalQ(cl, {
    library(frame); library(frameAnalyses); library(dplyr); library(extraDistr)
  })
  
  # 3) Split Conditions
  rows <- split(Conditions, seq_len(nrow(Conditions)))
  
  # 4) Hold per-stand results
  runs_all     <- vector("list", length(Forests))
  ip_all       <- vector("list", length(Forests))
  shelter_all  <- vector("list", length(Forests))
  exposed_all  <- vector("list", length(Forests))
  
  # Iterate stands
  for (i in seq_along(Forests)) {
    forest <- if (!is.null(names(Forests))) names(Forests)[i] else paste0("stand_", i)
    stand  <- Forests[[i]]
    Flora  <- stand[[1]]
    Structure <- stand[[2]]
    default.species.params <- stand[[3]]
    
    cat("Analysing forest", forest, "\n")
    cat("Running", reps * length(rows), "probabilistic replicates with load-balanced parallel processing\n")
    
    # Export stand-specific objects
    parallel::clusterExport(
      cl,
      varlist = c("Flora","Structure","default.species.params","distance","trail","reps", "testN",
                  "forest","fLine","Trees","edgeDistance","Altitude","wind_az_deg", "branchTrunkRat", "minHeight",
                  "hollowHprop","stags", "deadHollow", "wallThickness_min","wallThickness_max","lengthHollow_min","lengthHollow_max",
                  "diameterHollow_min","diameterHollow_max","hollow_tilt_deg","HTs","HTmin","tempDepression","TDs",
                  "climbProp","glideAngle","climbingSpeed","shelter","glide", "keepCanopy",
                  "bodyLength","surfaceArea","bodyMass","fibreLength","fibreCount","fibreDiameter","fibreCp","fiberSolid",
                  "skinCp","skinDensity","skinK","bodyTemp","xRate","yRate","yMax","diameter","surfDecl","var"),
      envir = environment()
    )
    
    parT <- parallel::parLapplyLB(cl, rows, parArboreal)
    
    # Collapse this stand’s results
    runs_i <- data.table::rbindlist(lapply(parT, `[[`, 1), use.names = TRUE, fill = TRUE)
    ip_i   <- data.table::rbindlist(lapply(parT, `[[`, 2), use.names = TRUE, fill = TRUE)
    
    she_i <- data.table::rbindlist(
      Filter(Negate(is.null), lapply(parT, `[[`, 3)),
      use.names = TRUE, fill = TRUE
    )
    exp_i <- data.table::rbindlist(
      Filter(Negate(is.null), lapply(parT, `[[`, 4)),
      use.names = TRUE, fill = TRUE
    )
    
    # Tag with forest and store
    runs_all[[i]]     <- dplyr::mutate(runs_i, forest = forest)
    ip_all[[i]]       <- dplyr::mutate(ip_i,   forest = forest)
    shelter_all[[i]]  <- if (nrow(she_i)) dplyr::mutate(she_i, forest = forest) else NULL
    exposed_all[[i]]  <- if (nrow(exp_i)) dplyr::mutate(exp_i, forest = forest) else NULL
  }
  
  # 5) Final bind across stands (NULL-safe)
  Runs    <- data.table::rbindlist(runs_all,   use.names = TRUE, fill = TRUE)
  IP      <- data.table::rbindlist(ip_all,     use.names = TRUE, fill = TRUE)
  Shelter <- data.table::rbindlist(Filter(Negate(is.null), shelter_all), use.names = TRUE, fill = TRUE)
  Exposed <- data.table::rbindlist(Filter(Negate(is.null), exposed_all), use.names = TRUE, fill = TRUE)
  
  list(Runs = Runs, IP = IP, Shelter = Shelter, Exposed = Exposed)
}


#________________________________________________________________________________
#' Helper function to run fire and impact model for arboreal fauna
#'
#' @param r A single row from the Conditions table
#'
#' @returns List of dataframes:
#' 1) Runs - detailed results for each replicate;
#' 2) IP - summary of ignition paths in each plant;
#' 3) Shelter_list - summary of mortalities associated with sheltering animals
#' 4) Exposed_list - summary of mortalities associated with exposed animals
#' @export
#'
#' 
parArboreal <- function(r) { 
  
  Runs_list   <- vector("list", reps)
  IP_list     <- vector("list", reps)
  Shelter_list <- vector("list", reps)
  Exposed_list <- vector("list", reps)
  
  base.params <- frame::buildParams(Structure, Flora, default.species.params, a = 1,
                                    fLine = fLine, slope = r$slope[1], temp = r$temp[1], dfmc = r$DFMC[1], wind = r$wind[1])
  
  # Adjust LFMC between 80% and 120% of input based on Drought Factor
  Pm <- 1.6-0.08*r$DF[1]
  
  #MODEL EACH FIRE
  for (j in seq_len(reps)) {
    db.path <- paste(forest, "_rep", r$Conditions[1], "var", j, ".db", sep = "")
    
    # Vary plant traits for each species within their range
    test  <- 0
    tries <- 0
    while (test == 0 && tries < 100) {
      TBL <- frameAnalyses::canopyCheck(
        base.params, testN = testN, Flora, Structure,
        a = 1, footprint = 1, l = 0.1, Ms = 0.01, Pm = Pm, Mr = 1.5, threshold = minHeight/hollowHprop,
        keepCanopy = keepCanopy
      )
      Species <- frame::species(TBL)
      
      top_st <- max(Species$st, na.rm = TRUE)
      if (!is.finite(top_st)) stop("All Species$st are NA; cannot determine top stratum.")
      
      # get name/comp for species in the top stratum that are in Trees
      treeList <- Species[Species$st == top_st & Species$name %in% Trees$spName, c("name","comp")]
      
      # proceed only if we actually have at least one candidate row
      if (nrow(treeList) > 0 && sum(treeList$comp, na.rm = TRUE) > 0) {
        test <- 1
      }
      
      tries <- tries + 1
    }
    
    if (test == 0) {
      stop("Failed to generate a non-empty treeList after 100 tries.")
    }
    
    
    # Collect behaviour and environmental attributes
    Strategy <- if(shelter > runif(1) & stags <= runif(1)) {"Shelter"} else if (glide > runif(1)) {"glide"} else {"climb"}
    
    probs <- treeList$comp
    probs[is.na(probs)] <- 0
    if (!length(probs) || sum(probs) <= 0) probs <- rep(1, length(treeList$name))
    
    first_or_na <- function(x) if (length(x)) x[1] else NA
    
    hollowTree <- sample(treeList$name, 1, prob = probs)
    htN_idx <- which(TBL$value == hollowTree & TBL$stratum == top_st)
    htN     <- as.numeric(first_or_na(TBL$species[htN_idx]))
    escapeTree <- sample(treeList$name, 1, prob = probs)
    etN_idx <- which(TBL$value == escapeTree & TBL$stratum == top_st)
    etN     <- first_or_na(TBL$species[etN_idx])
    
    siteTrees <- r %>%
      mutate("Forest" = forest,
             "Strategy" = Strategy,
             "hollowTreeSp" = hollowTree, "escapeTreeSp" = escapeTree,
             "hollowTreeH" = as.numeric(first_or_na(TBL$value[TBL$species == htN & TBL$param == "hp"])),  
             "treeSpacing" = Structure$Can[1], 
             "escapeTreeH" = as.numeric(first_or_na(TBL$value[TBL$species == etN & TBL$param == "hp"]))) %>%
      mutate("hollowH" = hollowTreeH * hollowHprop, "climbH_hollow" = hollowTreeH * climbProp) %>%
      mutate("climbH_escape" = max(escapeTreeH * climbProp, 1),
             "landingH" = pmax(pmin(climbH_hollow - treeSpacing * tan((glideAngle * pi/180)), climbH_escape - 1),0))
    hollowTreeDat <- Trees[Trees$spName == siteTrees$hollowTreeSp[1], ]
    if (nrow(hollowTreeDat) == 0) stop("No trait row for hollowTreeSp: ", siteTrees$hollowTreeSp[1])
    
    
    # 1. FIRE BEHAVIOUR
    ffm_run(TBL, db.path = db.path, db.recreate = TRUE)
    res<-ffm_db_load(db.path)# Collect results
    outa <- frameSummaryBeta(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults, res$IgnitionPaths) %>%
      mutate(tm = j,
             repId = r$Conditions[1])
    outb <- repFlame(res$IgnitionPaths) %>%
      mutate(tm = j,
             repId = r$Conditions[1])
    
    scorch <- suppressMessages(frame::frameSeverity(outa, outb, Param = TBL, Test = 80)) %>%
      select(!wind_kph)
    
    # 2. IMPACTS ON GLIDERS
    hollowCharacteristics <- randomiseHollow(branchTrunkRat = branchTrunkRat, wallThickness_min = wallThickness_min, wallThickness_max = wallThickness_max,
                                             lengthHollow_min = lengthHollow_min, lengthHollow_max = lengthHollow_max,
                                             diameterHollow_min = diameterHollow_min, diameterHollow_max = diameterHollow_max,
                                             hollow_tilt_deg = hollow_tilt_deg, HTs = HTs, HTmin = HTmin)
    
    outc <- NULL
    outd <- NULL 
    
    # Modify hollow descriptions based on whether hollow is living or dead 
    hollowLiving <- runif(1) > deadHollow
    woodM <- if(hollowLiving) {Trees$moisture[Trees$spName == siteTrees$hollowTreeSp[1]]} else {r$DFMC[1]}
    bark <- if(hollowLiving) {Trees$bark[Trees$spName == siteTrees$hollowTreeSp[1]]} else {0}
    resBark <- if(hollowLiving) {Trees$resBark[Trees$spName == siteTrees$hollowTreeSp[1]]} else {0}
    bMoisture <- if(hollowLiving) {Trees$bMoisture[Trees$spName == siteTrees$hollowTreeSp[1]]} else {r$DFMC[1]}
    
    if (Strategy == "Shelter") {
      Effect <- hollow(Surf = outa, IP = outb, TBL,
                       Height = siteTrees$hollowH[1], woodDensity = hollowTreeDat$woodDensity[1], barkDensity = hollowTreeDat$barkDensity[1],
                       wood = hollowCharacteristics$wallThickness[1], bark = bark,
                       comBark = Trees$comBark[Trees$spName == siteTrees$hollowTreeSp[1]], resBark = resBark,
                       RH = r$RH[1], moisture = Trees$moisture[Trees$spName == siteTrees$hollowTreeSp[1]],
                       bMoisture = bMoisture, distance = distance, trail = trail,
                       var = var, Pressure = r$Pressure[1], Altitude = Altitude, Dimension = hollowCharacteristics$lengthHollow[1],
                       wind_az_deg = 270, Area = hollowCharacteristics$lengthHollow[1] * (hollowCharacteristics$diameterHollow[1] +
                                                                                            hollowCharacteristics$wallThickness[1] * 2),
                       diameterHollow = hollowCharacteristics$diameterHollow[1], diameter = diameter, surfDecl = surfDecl,
                       startTemp = r$temp[1] - hollowCharacteristics$tempDep[1], hollow_az_deg = hollowCharacteristics$hollow_az_deg[1],
                       hollow_tilt_deg = hollowCharacteristics$hollow_tilt_deg[1]
      )
      
      impRes <- siteTrees %>%
        mutate(Altitude = Altitude,
               hollowLocation = hollowCharacteristics$locationHollow[1],
               hollowLiving = hollowLiving,
               woodDensity = hollowTreeDat$woodDensity[1],
               wallThickness = hollowCharacteristics$wallThickness[1],
               wallMoisture = woodM,
               barkDensity = ifelse(hollowLiving, hollowTreeDat$barkDensity[1], NA),
               barkThickness = bark,
               barkMoisture = ifelse(hollowLiving, bMoisture, NA),
               lengthHollow = hollowCharacteristics$lengthHollow[1],
               diameterHollow = hollowCharacteristics$diameterHollow[1],
               hollow_tilt_deg = hollowCharacteristics$hollow_tilt_deg[1],
               Burns = NA,
               Mortality = if (max(Effect$VPmortality, na.rm = TRUE) == 1) "Asphyxiation" else "Survive"
        ) %>%
        mutate(Stat = if_else(Mortality != "Survive", 1, 0))
      
      outc <- Effect %>%
        mutate(Strategy = "Shelter",
               repId = r$Conditions[1],
               tm = j)
      
    } else if (Strategy == "climb") {
      Effect <- mammal(Surf = outa, IP = outb,
                       Height = siteTrees$hollowH[1],distance = distance, trail = trail, diameter = diameter,
                       surfDecl = surfDecl, var = var, Pressure = r$Pressure[1], Altitude = Altitude, RH = r$RH[1],
                       bodyLength = bodyLength, surfaceArea = surfaceArea, bodyMass = bodyMass, protection = fibreLength,
                       fibreCount = fibreCount, fibreDiameter = fibreDiameter, fibreCp = fibreCp, fiberSolid = fiberSolid, 
                       skinCp = skinCp, skinDensity = skinDensity, skinK = skinK, bodyTemp = bodyTemp, xRate = xRate, yRate = yRate, 
                       yMax = siteTrees$climbH_hollow[1], 
      )
      impRes <- siteTrees %>%
        mutate(Altitude = Altitude,
               hollowLocation = hollowCharacteristics$locationHollow[1],
               hollowLiving = hollowLiving,
               woodDensity = hollowTreeDat$woodDensity[1],
               wallThickness = hollowCharacteristics$wallThickness[1],
               wallMoisture = woodM,
               barkDensity = ifelse(hollowLiving, hollowTreeDat$barkDensity[1], NA),
               barkThickness = bark,
               barkMoisture = ifelse(hollowLiving, bMoisture, NA),
               lengthHollow = hollowCharacteristics$lengthHollow[1],
               diameterHollow = hollowCharacteristics$diameterHollow[1],
               hollow_tilt_deg = hollowCharacteristics$hollow_tilt_deg[1],
               Burns = case_when(
                 max(Effect$B3, na.rm = TRUE) == 1 ~ "Full thickness",
                 max(Effect$B2, na.rm = TRUE) == 1 ~ "Partial thickness",
                 max(Effect$B1, na.rm = TRUE) == 1 ~ "Superficial",
                 TRUE ~ "NA"
               ),
               Mortality = if (max(Effect$B2, na.rm = TRUE) == 1) "Burnt" else 
                 if (max(Effect$VPmortality, na.rm = TRUE) == 1) "Asphyxiation" else "Survive"
        ) %>%
        mutate(Stat = if_else(Mortality != "Survive", 1, 0))
      
      outd <- Effect %>%
        mutate(Strategy = "Climb",
               repId = r$Conditions[1],
               tm = j)
      
    } else {
      Effect <- mammal(Surf = outa, IP = outb,
                       Height = siteTrees$landingH[1],distance = distance, trail = trail, diameter = diameter,
                       surfDecl = surfDecl, var = var, Pressure = r$Pressure[1], Altitude = Altitude, RH = r$RH[1],
                       bodyLength = bodyLength, surfaceArea = surfaceArea, bodyMass = bodyMass, protection = fibreLength,
                       fibreCount = fibreCount, fibreDiameter = fibreDiameter, fibreCp = fibreCp, fiberSolid = fiberSolid, 
                       skinCp = skinCp, skinDensity = skinDensity, skinK = skinK, bodyTemp = bodyTemp, xRate = xRate, yRate = yRate, 
                       yMax = siteTrees$climbH_escape[1], 
      )
      impRes <- siteTrees %>%
        mutate(Altitude = Altitude,
               hollowLocation = hollowCharacteristics$locationHollow[1],
               hollowLiving = hollowLiving,
               woodDensity = hollowTreeDat$woodDensity[1],
               wallThickness = hollowCharacteristics$wallThickness[1],
               wallMoisture = woodM,
               barkDensity = ifelse(hollowLiving, hollowTreeDat$barkDensity[1], NA),
               barkThickness = bark,
               barkMoisture = ifelse(hollowLiving, bMoisture, NA),
               lengthHollow = hollowCharacteristics$lengthHollow[1],
               diameterHollow = hollowCharacteristics$diameterHollow[1],
               hollow_tilt_deg = hollowCharacteristics$hollow_tilt_deg[1],
               Burns = case_when(
                 max(Effect$B3, na.rm = TRUE) == 1 ~ "Full thickness",
                 max(Effect$B2, na.rm = TRUE) == 1 ~ "Partial thickness",
                 max(Effect$B1, na.rm = TRUE) == 1 ~ "Superficial",
                 TRUE ~ "NA"
               ),
               Mortality = if (max(Effect$B2, na.rm = TRUE) == 1) "Burnt" else 
                 if (max(Effect$VPmortality, na.rm = TRUE) == 1) "Asphyxiation" else "Survive"
        ) %>%
        mutate(Stat = if_else(Mortality != "Survive", 1, 0))
      
      outd <- Effect %>%
        mutate(Strategy = "glide",
               repId = r$Conditions[1],
               tm = j)
      
    }
    outa <- suppressMessages(left_join(outa, scorch, by = "repId")) %>%
      left_join(impRes, by = c("repId" = "Conditions"))
    
    Runs_list[[j]]   <- outa
    IP_list[[j]]     <- outb
    Shelter_list[[j]] <- outc
    Exposed_list[[j]] <- outd
  }
  
  Runs   <- dplyr::bind_rows(Runs_list)
  IP     <- dplyr::bind_rows(IP_list)
  She <- dplyr::bind_rows(Shelter_list)
  Exp <- dplyr::bind_rows(Exposed_list)
  
  return(list(Runs, IP, She, Exp))
  
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
#' @param IP The dataframe 'IP' exported from Monte Carlos as 'IP.csv'.
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

underground <- function(Surf, IP, diameter = 6, surface = 677, percentile = 0.5, RH = 0.2,
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
  Ca <- threat(Surf, IP, HA, Height=0, var, Pressure, Altitude)%>%
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
    Cb <-threat(Surf, IP, HA, Height=0, var, Pressure, Altitude) %>%
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
  return(Ca)
}