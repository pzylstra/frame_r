#' @title Calculate phi, the view factor
#' @description Configuration or view factor
#'
#' Calculates the configuration factor for radiative heat transfer as used
#' in AS3959 for bushfire risk assessment of built structures
#'
#' @param Fl Flame length (m)
#' @param Fa Flame angle (degrees)
#' @param S Slope angle (degrees)
#' @param D Horizontal distance from the flame origin to the point (m)
#' @param H Vertical distance from the flame origin to the point (m)


phi_fn <- function(Fl, Fa, S, D, H) {
  Fx  <- (Fl * 0.5) * cos(Fa)
  Sep <- abs(D - Fx)
  Sl  <- (S * pi) / 180
  tSl <- tan(Sl)
  
  n <- length(Fl)
  out <- numeric(n)
  
  # Masks
  m_smallFl <- (Fl < 0.1)
  m_badSep  <- (!m_smallFl) & (Sep <= 0)
  m_eval    <- (!m_smallFl) & (Sep > 0)
  
  # Assign trivial cases
  out[m_smallFl] <- 0
  out[m_badSep]  <- 1
  
  if (any(m_eval)) {
    i <- which(m_eval)
    
    Sep_i <- Sep[i]
    Fl_i  <- Fl[i]
    Fa_i  <- Fa[i]
    
    X1 <- (Fl_i * sin(Fa_i) - (Fx[i] + D) * tSl - H) / Sep_i
    X2 <- (H + Sep_i * tSl) / Sep_i
    Y1 <- 50 / Sep_i
    
    X1a <- sqrt(1 + X1^2)
    X2a <- sqrt(1 + X2^2)
    Y1a <- sqrt(1 + Y1^2)
    
    # Use atan2 for stability: atan(y/x) == atan2(y, x)
    term <- (1/pi) * (
      (X1 / X1a) * atan2(Y1, X1a) + (Y1 / Y1a) * atan2(X1, Y1a) +
        (X2 / X2a) * atan2(Y1, X2a) + (Y1 / Y1a) * atan2(X2, Y1a)
    )
    
    out[i] <- term
  }
  
  out
}


#' @title tau
#' @description Atmospheric attenuation
#'
#' Calculates the atmospheric transmissivity of radiation as used
#' in AS3959 for bushfire risk assessment of built structures. Equations from
#' Fuss & Hamins
#' An estimate of the correction applied to radiant flame measurements due to
#' attenuation by atmospheric CO2and H2O.
#' Fire Saf. J. 37, 181–190 (2002)
#'
#' @param D Horizontal distance from the flame origin to the point (m)
#' @param flameTemp Flame temperature 
#' @param temperature Air temperature 
#' @param rh Relative humidity


tau <- function(D = 200, flameTemp = 1300, temperature = 288, rh = 0.51)
{
  
  a0 <-1.486-0.002003*temperature+0.0000468*flameTemp-0.06052*rh
  a1 <-0.01225-0.000059*temperature+0.00000166*flameTemp-0.001759*rh
  a2 <--0.0001489+0.0000006893*temperature-0.00000001922*flameTemp+0.00002092*rh
  a3 <-0.0000008381-0.000000003823*temperature+0.00000000010511*flameTemp-0.0000001166*rh
  a4 <--0.000000001685+0.000000000007637*temperature-0.0000000000002085*flameTemp+0.0000000002350*rh
  return(a0+a1*D+a2*D^2+a3*D^3+a4*D^4)
}



#' @title threat
#' @description Calculates the convective and radiative heat from a flame that is 
#' incident upon a designated point
#'
#' Finds the temperature, velocity, dynamic viscosity, atmospheric pressure, and density of a plume at a point
#' Gives the mean flame temperature (K) and emissive power of each flame
#'
#' tempAir: Air temperature is modelled from dynamic flame segments using
#' Weber R.O., Gill A.M., Lyons P.R.A., Moore P.H.R., Bradstock R.A., Mercer G.N. (1995)
#' Modelling wildland fire temperatures. CALMScience Supplement, 4, 23–26.
#'
#' cpAir: Specific heat of air is found from an empirical function fit to data from
#' Hilsenrath, J. et al.
#' Circular of the Bureau of Standards no. 564: tables of thermal properties of gases comprising tables of
#' thermodynamic and transport properties of air, argon, carbon dioxide, carbon monoxide hydrogen, nitrogen,
#' oxygen, and steam. (Department of Commerce, 1955), and
#' Kyle, B. G. Chemical and process thermodynamics. (Englewood Cliffs / Prentice Hall, 1984)
#'
#' Density: Air density (kg/m^3) is found using the Ideal Gas Law.
#'
#' presAtm: Air pressure (hPa) is calculated using the Barometric Formula with standard Values.
#'
#' viscosity: Dynamic viscosity (10^-6 Pa.s) is calculated using an empirical relationship fit to a subset of
#' air pressures within an order of magnitude of sea level (R2 = 0.999), using data from
#' Kadoya, K., N, M. & Nagashima, A. Viscosity and thermal conductivity of dry air in the gaseous phase.
#' J. Phys. Chem. Ref. Data 14, 947–970 (1985)
#'
#' Plume_velocity: Plume velocity (m/s) is calculated allowing for crosswinds, using
#' Oka Y., Sugawa O., Imamura T. (2008)
#' Correlation of temperature rise and velocity along an inclined fire plume axis in crosswinds.
#' Fire Safety Journal, 43, 391–400.
#'
#' The Froude number is set to 1.5 as per Ma T.G., Quintiere J.G. (2003)
#' Numerical simulation of axi-symmetric fire plumes: Accuracy and limitations.
#' Fire Safety Journal, 38, 467–492.
#'
#' flameTemp: Average flame temperature is integrated from the modelled temperature over the length of the flame
#'
#' epsilon: The emissivity of a flame is calculated from mean flame thickness using
#' Àgueda, A., Pastor, E., Pérez, Y. & Planas, E.
#' Experimental study of the emissivity of flames resulting from the combustion of forest fuels.
#' Int. J. Therm. Sci. 49, 543–554 (2010)
#'
#' phi: The configuration factor is calculated using the function phi_fn, from
#' Tan, Z., Midgley, S. & Douglas, G.
#' A computerised model for bushfire attack assessment and its applications in bushfire protection planning.
#' Int. Congr. Model. Simul. Adv. Appl. Manag. Decis. Making, MODSIM05 538–545 (2005)
#'
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param Horizontal The horizontal distance in metres from the flame origin to the point
#' @param Height The vertical distance in metres from the flame origin to the point
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param residence Surface flame residence time
#' @param IP The dataframe produced by the function 'repFlame'.
#' @param surfDecl Exponent describing the rate of post-front flame decay in surface litter 
#'
#' @return dataframe
#' @export

threat <- function (Surf, IP, Horizontal = 10, Height = 10, var = 10,
                    Pressure = 1013.25, Altitude = 0, residence = 30, surfDecl = 2)
{
  # Shared terms
  slope0         <- Surf$slope_degrees[1]
  El             <- Horizontal * tan((slope0 * pi)/180)
  interceptAngle <- atan((Height + El) / Horizontal)
  tan_var        <- tan((var * pi)/180)
  cos_intercept  <- cos(interceptAngle)
  pAlpha_common  <- abs(Horizontal / cos_intercept)
  
  # Calculate surface flame heating
  s <- Surf %>%
    dplyr::mutate(
      hor         = Horizontal,
      Intercept   = ifelse(((interceptAngle - tan_var) > angleSurface) |
                             ((interceptAngle + tan_var) < angleSurface), 0, 1),
      Alpha       = 1/(2 * lengthSurface^2),
      C           = 950 * lengthSurface * exp(-Alpha * lengthSurface^2),
      pAlphas     = pAlpha_common,
      temp_pointS = ifelse(pAlphas < lengthSurface,
                           950 + exp(-Alpha * pAlphas^2),
                           C/pAlphas) * Intercept + temperature
    ) %>%
    dplyr::select(
      repId, hor, ros_kph, wind_kph, angleSurface, lengthSurface, pAlphas,
      temperature, slope_degrees, temp_pointS, epsilon
    )
  
  # Find post-front heating
  post <- s %>%
    dplyr::mutate(
      t               = max(0, -hor) / (ros_kph/3.6),
      tail            = min(1, max(ceiling(t), 0)),
      surfPost        = lengthSurface * exp(-(surfDecl / residence) * t),
      Alpha           = 1/(2 * surfPost^2),
      C               = 950 * surfPost * exp(-Alpha * surfPost^2),
      pAlphaPost      = pmax(0, Height),
      temp_point_post = tail * ifelse(pAlphaPost < surfPost,
                                      950 + exp(-Alpha * pAlphaPost^2),
                                      C/pAlphaPost) + temperature
    )
  
  # Calculate heating from burning plants
  if (nrow(IP) == 0) {
    # No plant flame: carry 'post' forward and give neutral plant terms.
    p <- post %>%
      dplyr::mutate(
        # No plant plume contribution; NA so it's ignored by pmax(…, na.rm = TRUE)
        temp_pointP = NA_real_,
        flameTempP  = NA_real_,
        # Use surface geometry as the default "representative" length/angle
        repLength   = lengthSurface,
        repAngle    = angleSurface,
        # Keep a defined distance for downstream pAlpha logic
        pAlphap     = pAlpha_common
      )
  } else {
    # Plant flame present: do the join and compute plant plume temps, but guard zero/NaN cases
    p <- suppressMessages(
      IP %>%
        dplyr::left_join(post, by = "repId") %>%
        dplyr::mutate(
          Angle = dplyr::if_else(
            (x1 == x0) & (y1 == y0),
            0,
            atan2(y1 - y0, x1 - x0)
          ),
          InterceptP = ifelse(
            ((interceptAngle - tan_var) > Angle) |
              ((interceptAngle + tan_var) < Angle), 0, 1
          ),
          # guard denominators
          denom      = 2 * pmax(1e-9, abs(flameLength)) * pmax(1e-9, abs(flameLength - length)),
          AlphaP     = 1 / denom,
          CP         = 950 * pmax(1e-9, flameLength) * exp(-AlphaP * (flameLength - length)^2),
          pAlphap    = pAlpha_common,
          temp_pointP = ifelse(
            pAlphap < pmax(1e-9, flameLength),
            950 + exp(-AlphaP * (pAlphap - length)^2),
            CP / pmax(1e-9, pAlphap)
          ) * InterceptP + temperature,
          flameTempP = dplyr::if_else(
            flameLength <= 0 | is.na(flameLength),
            NA_real_,
            1045 * exp(0.155 * (length / flameLength))
          ),
          # provide defaults for downstream geometry use
          repLength = dplyr::coalesce(repLength, lengthSurface),
          repAngle  = dplyr::coalesce(repAngle,  angleSurface)
        ) %>%
        dplyr::select(-denom)
    )
  }
  
  # Calculate heat transfer inputs
  Con <- p %>%
    dplyr::mutate(
      tempAir = pmax(temp_pointS, temp_pointP, temp_point_post, na.rm = TRUE),
      
      # if all sources were NA (extremely unlikely), fall back to ambient
      tempAir = dplyr::coalesce(tempAir, temperature),
      
      cpAir   = 0.00019 * (tempAir + 273.13) + 0.9405,
      pAlpha  = pmax(0.05, dplyr::case_when(
        !is.na(temp_pointS)     & tempAir == temp_pointS     ~ pAlphas,
        !is.na(temp_point_post) & tempAir == temp_point_post ~ pAlphaPost,
        TRUE                                             ~ pAlphap
      )),
      viscosity      = 0.000000388 * (tempAir + 273.15)^0.6818,
      presAtm        = (Pressure/10) * exp(-0.00012 * Altitude),
      Density        = (presAtm * 1000) / (287.05 * (tempAir + 273.15)),
      Plume_velocity = sqrt(
        ((1.24 * (pAlpha/(pAlpha + 2*hor))) *
           sqrt(((tempAir - temperature)/(temperature + 273.15)) * 9.81 * pAlpha))^2 +
          (((wind_kph/3.6) * 0.47 * pmax((pAlpha + hor)/pAlpha, 0)^0.25))^2
      ),
      flameTemp = dplyr::if_else(is.na(flameTempP) | flameTempP <= 0, 1045, flameTempP),
      E              = pmax(0, epsilon * 0.0000000567 * (flameTemp^4 - (tempAir + 273.15)^4)),
      # Use plant angle/length only if plants exceed surface length
      repAngle  = dplyr::coalesce(repAngle,  angleSurface),
      repLength = dplyr::coalesce(repLength, lengthSurface),
      phi       = phi_fn(repLength, repAngle, slope_degrees, Horizontal, Height),
      qr        = E * phi
    ) %>%
    dplyr::select(repId, ros_kph, wind_kph, temperature, lengthSurface, pAlpha, tempAir, cpAir,
                  viscosity, presAtm, Density, Plume_velocity, flameTemp, epsilon, E, phi, qr
    )
  
  Con
}






#' @title cp
#' @description Specific heat capacity of a material
#'
#' @param Material Type of material, either 'bark' or 'wood'
#' @param temp Temperature (C)
#' @param moist Moisture content (proportion oven dry weight)
#'
#' @return value
#' @export
#'

cp <- function(Material = "wood", temp, moist){
  cp <- if (Material == "bark") {
    (1105+4.85*temp)*(1-moist)+moist*4185+1276*moist
  } else {
    1080+408*moist+2.53*temp+6.28*moist*temp
  }
  
  return(cp)
}

k <- function(Material = "bark", temp, moist, density){
  kAir <- 0.00028683*(temp+273.15)^0.7919
  cp <- if (Material == "bark") {
    rhoM <- (moist+moist^2)*density
    (2.104*density+5.544*rhoM+3.266*temp-166.216)*10^-4
  } else {
    kWood(temp, density, kAir)
  }
  
  return(cp)
}


#' @title hFauna
#'
#' @param Shape Shape of the object, either 'Flat', 'Sphere' or 'Cylinder'
#' @param Re Reynolds number
#'
#' @return value
#' @export
#'

hFauna <- function(Shape = "Cylinder", Re) {
  hFlat <- ifelse(Re > 300000,0.037*Re^(4/5)*0.888, 0.66*Re^0.5*0.888)
  hSphere <- 2 + 0.6*Re^(0.5)*0.888
  hCylinder <- 0.35 + 0.47*Re^(1/2)*0.837
  h <- ifelse(Shape == "Flat", hFlat, ifelse(Shape == "Sphere", hSphere, hCylinder)) 
  return(h)
}