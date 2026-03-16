#' @importFrom dplyr %>% arrange mutate count filter
#' 
#' @title floraTarget
#' @description  Scorch height using an isotherm, including target species
#'
#' Calculates the height to which vegetation will be consumed,
#' and the height to a designated temperature isotherm reached for one second.
#'
#' Output fields are:
#' Height - overall scorch height (m)
#' ns, e, m, c - height of consumption in each stratum (m)
#' b1, b2, b3, b4 - percentage of strata 1 (ns) to strata 4 (c) consumed
#' sc1, sc2, sc3, sc4 - percentage of strata 1 (ns) to strata 4 (c) scorched
#' d1, d2, d3, d4 - death (1) or survival (0) of standing foliage per stratum (50% or more scorch)
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Plant The dataframe produced by the function 'repFlame'.
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @param Test The temperature of the flora, default 70 degC
#' @param targSp Name of species to be watched for scorch
#'
#' @return dataframe
#' @export

floraTarget <- function(Surf, Plant, Param = Param, targSp = "a", Test = 70)
{
  S <- strata(Param)
  n <- max(S$stratum)
  
  #Max burn height
  c <- Plant[Plant$level == "Canopy", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(c = y1)%>%
    select(repId, c)%>%
    right_join(Surf)
  c[is.na(c)] <- 0
  m <- Plant[Plant$level == "MidStorey", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(m = y1)%>%
    select(repId, m)%>%
    right_join(c)%>%
    mutate(m = ifelse(c>0, S$top[3], m))
  m[is.na(m)] <- 0
  e <- Plant[Plant$level == "Elevated", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(e = y1)%>%
    select(repId, e)%>%
    right_join(m)%>%
    mutate(e = ifelse(m>0, S$top[2], e))
  e[is.na(e)] <- 0
  ns <- Plant[Plant$level == "NearSurface", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(ns = y1)%>%
    select(repId, ns)%>%
    right_join(e)%>%
    mutate(ns = ifelse(e>0, S$top[1], ns))
  ns[is.na(ns)] <- 0
  
  #Heating from surface fire
  ground <- ns%>%
    mutate(Alpha = 1/(2*lengthSurface^2),
           C = 950*lengthSurface*exp(-Alpha*lengthSurface^2),
           pAlpha = abs(C/(Test-temperature)),
           R = pAlpha*cos(angleSurface),
           El = R * tan((slope_degrees * pi)/180),
           ht = (sin(angleSurface)*pAlpha - El) * extinct,
           ns = ns * extinct,
           e = e * extinct,
           m = m * extinct,
           c = c * extinct
    )%>%
    select(repId, extinct, temperature, slope_degrees, wind_kph, ht, ns, e, m, c)
  
  #Heating from plants
  plant<-Plant%>%
    left_join(ground)%>%
    mutate(Angle = atan((y1-y0)/(x1-x0)),
           Alpha = 1/(2*flameLength*(flameLength-length)),
           C = 950*flameLength*exp(-Alpha*(flameLength-length)^2),
           pAlpha = abs(C/(Test-temperature)),
           Reach = pAlpha*cos(Angle),
           El = Reach * tan((slope_degrees * pi)/180),
           htP = (sin(Angle)*pAlpha + y0 - El) * extinct)%>%
    select(repId, htP)%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    right_join(ground)
  plant[is.na(plant)] <- 0
  
  #Find the corresponding number of each stratum
  NS <- S[S$name == "near surface", ]
  nNS <-ifelse(count(NS)==0, 0,as.numeric(as.character(NS$stratum)))
  E <- S[S$name == "elevated", ]
  nE <- ifelse(count(E)==0, 0,as.numeric(as.character(E$stratum)))
  M <- S[S$name == "midstorey", ]
  nM <-ifelse(count(M)==0, 0,as.numeric(as.character(M$stratum)))
  C <- S[S$name == "canopy", ]
  nC <- ifelse(count(C)==0, 0,as.numeric(as.character(C$stratum)))
  nSp <- targSp %in% unique(Param$value[which(Param$param=="name")])
  
  
  #Final heating
  Iso <- plant%>%
    mutate(Height = pmax(ht, htP),
           nsBase = ifelse(count(NS)==0, 0, S$base[which(S$name == "near surface")]),
           eBase = ifelse(count(E)==0, 0, S$base[which(S$name == "elevated")]),
           mBase = ifelse(count(M)==0, 0, S$base[which(S$name == "midstorey")]),
           cBase = ifelse(count(C)==0, 0, S$base[which(S$name == "canopy")]),
           nsTop = ifelse(count(NS)==0, 0, S$top[which(S$name == "near surface")]),
           eTop = ifelse(count(E)==0, 0, S$top[which(S$name == "elevated")]),
           mTop = ifelse(count(M)==0, 0, S$top[which(S$name == "midstorey")]),
           cTop = ifelse(count(C)==0, 0, S$top[which(S$name == "canopy")]),
           spBase = ifelse(nSp, as.numeric(filter(Param, species==targSp & param == "hc")$value), 0),
           spTop = ifelse(nSp, as.numeric(filter(Param, species==targSp & param == "hp")$value), 0),
           #Calculate scorch and burn
           b1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(ns-nsBase)/(nsTop-nsBase))),0)),
           b2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(e-eBase)/(eTop-eBase))),0)),
           b3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(m-mBase)/(mTop-mBase))),0)),
           b4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(c-cBase)/(cTop-cBase))),0)),
           sc1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(Height-nsBase)/(nsTop-nsBase))),0)),
           sc2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(Height-eBase)/(eTop-eBase))),0)),
           sc3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(Height-mBase)/(mTop-mBase))),0)),
           sc4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(Height-cBase)/(cTop-cBase))),0)),
           severity = case_when(b4 > 50 ~ 5,
                                b4 > 0 ~ 4,
                                sc4 > 0 ~ 3,
                                b1+b2 > 50 ~ 2,
                                b1+b2 > 0 ~ 1,
                                TRUE ~ 0),
           targSp = ifelse(spTop == 0, 0, round(pmin(100,pmax(0,100*(Height-spBase)/(spTop-spBase))),0)))%>%
    select(repId, wind_kph, Height, ns, e, m, c, b1, b2, b3, b4, sc1, sc2, sc3, sc4, severity, targSp)
  
  return(Iso)
}

#' @title flora
#' @description Scorch height using an isotherm
#'
#' Calculates the height to which vegetation will be consumed,
#' and the height to a designated temperature isotherm reached for one second.
#'
#' Output fields are:
#' Height - overall scorch height (m)
#' ns, e, m, c - height of consumption in each stratum (m)
#' b1, b2, b3, b4 - percentage of strata 1 (ns) to strata 4 (c) consumed
#' sc1, sc2, sc3, sc4 - percentage of strata 1 (ns) to strata 4 (c) scorched
#' d1, d2, d3, d4 - death (1) or survival (0) of standing foliage per stratum (50% or more scorch)
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Plant The dataframe produced by the function 'repFlame'.
#' @param Param A parameter dataframe used for FRaME,
#' such as produced using readLegacyParamFile
#' @param Test The temperature of the flora, default 70 degC
#'
#' @return dataframe
#' @export

flora <- function(Surf, Plant, Param = Param, Test = 70)
{
  S <- strata(Param)
  n <- max(S$stratum)
  
  #Max burn height
  c <- Plant[Plant$level == "Canopy", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(c = y1)%>%
    select(repId, c)%>%
    right_join(Surf)
  c[is.na(c)] <- 0
  m <- Plant[Plant$level == "MidStorey", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(m = y1)%>%
    select(repId, m)%>%
    right_join(c)%>%
    mutate(m = ifelse(c>0, S$top[3], m))
  m[is.na(m)] <- 0
  e <- Plant[Plant$level == "Elevated", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(e = y1)%>%
    select(repId, e)%>%
    right_join(m)%>%
    mutate(e = ifelse(m>0, S$top[2], e))
  e[is.na(e)] <- 0
  ns <- Plant[Plant$level == "NearSurface", ]%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    mutate(ns = y1)%>%
    select(repId, ns)%>%
    right_join(e)%>%
    mutate(ns = ifelse(e>0, S$top[1], ns))
  ns[is.na(ns)] <- 0
  
  #Heating from surface fire
  ground <- ns%>%
    mutate(Alpha = 1/(2*lengthSurface^2),
           C = 950*lengthSurface*exp(-Alpha*lengthSurface^2),
           pAlpha = abs(C/(Test-temperature)),
           R = pAlpha*cos(angleSurface),
           El = R * tan((slope_degrees * pi)/180),
           ht = (sin(angleSurface)*pAlpha - El) * extinct,
           ns = ns * extinct,
           e = e * extinct,
           m = m * extinct,
           c = c * extinct
    )%>%
    select(repId, extinct, temperature, slope_degrees, wind_kph, ht, ns, e, m, c)
  
  #Heating from plants
  plant<-Plant%>%
    left_join(ground)%>%
    mutate(Angle = atan((y1-y0)/(x1-x0)),
           Alpha = 1/(2*flameLength*(flameLength-length)),
           C = 950*flameLength*exp(-Alpha*(flameLength-length)^2),
           pAlpha = abs(C/(Test-temperature)),
           Reach = pAlpha*cos(Angle),
           El = Reach * tan((slope_degrees * pi)/180),
           htP = (sin(Angle)*pAlpha + y0 - El) * extinct)%>%
    select(repId, htP)%>%
    group_by(repId) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    right_join(ground)
  plant[is.na(plant)] <- 0
  
  #Find the corresponding number of each stratum
  NS <- S[S$name == "near surface", ]
  nNS <-ifelse(count(NS)==0, 0,as.numeric(as.character(NS$stratum)))
  E <- S[S$name == "elevated", ]
  nE <- ifelse(count(E)==0, 0,as.numeric(as.character(E$stratum)))
  M <- S[S$name == "midstorey", ]
  nM <-ifelse(count(M)==0, 0,as.numeric(as.character(M$stratum)))
  C <- S[S$name == "canopy", ]
  nC <- ifelse(count(C)==0, 0,as.numeric(as.character(C$stratum)))
  
  
  #Final heating
  Iso <- plant%>%
    mutate(Height = pmax(ht, htP),
           nsBase = if(count(NS)==0){0} else {S$base[which(S$name == "near surface")]},
           eBase = if(count(E)==0) {0} else {S$base[which(S$name == "elevated")]},
           mBase = if(count(M)==0) {0} else {S$base[which(S$name == "midstorey")]},
           cBase = if(count(C)==0) {0} else {S$base[which(S$name == "canopy")]},
           nsTop = if(count(NS)==0) {0} else {S$top[which(S$name == "near surface")]},
           eTop = if(count(E)==0) {0} else {S$top[which(S$name == "elevated")]},
           mTop = if(count(M)==0) {0} else {S$top[which(S$name == "midstorey")]},
           cTop = if(count(C)==0) {0} else {S$top[which(S$name == "canopy")]},
           #Calculate scorch and burn
           b1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(ns-nsBase)/(nsTop-nsBase))),0)),
           b2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(e-eBase)/(eTop-eBase))),0)),
           b3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(m-mBase)/(mTop-mBase))),0)),
           b4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(c-cBase)/(cTop-cBase))),0)),
           sc1 = ifelse(nsTop == 0, 0, round(pmin(100,pmax(0,100*(Height-nsBase)/(nsTop-nsBase))),0)),
           sc2 = ifelse(eTop == 0, 0, round(pmin(100,pmax(0,100*(Height-eBase)/(eTop-eBase))),0)),
           sc3 = ifelse(mTop == 0, 0, round(pmin(100,pmax(0,100*(Height-mBase)/(mTop-mBase))),0)),
           sc4 = ifelse(cTop == 0, 0, round(pmin(100,pmax(0,100*(Height-cBase)/(cTop-cBase))),0)),
           severity = case_when(b4 > 50 ~ 7,
                                b4 > 0 ~ 6,
                                sc4 > 80 ~ 5,
                                sc4 > 20 ~ 4,
                                sc4 > 0 ~ 3,
                                b1+b2 > 50 ~ 2,
                                b1+b2 > 0 ~ 1,
                                TRUE ~ 0))%>%
    select(repId, wind_kph, Height, ns, e, m, c, b1, b2, b3, b4, sc1, sc2, sc3, sc4, severity)
  
  return(Iso)
}


#' Compute fire severity classes and metrics by replicate
#'
#' @description
#' Calculates per-replicate burn and scorch metrics for each vegetation stratum
#' (near surface, elevated, midstorey, canopy), combining heating from the
#' surface fire and from plant flames, then assigns an integer severity class
#' and a short description.
#'
#' @details
#' The function:
#' \itemize{
#'   \item extracts per-stratum maximum flame heights from \code{IP}, skipping strata with no usable values;
#'   \item applies a cascade rule: if an upper stratum burns, lower strata are set to their top height;
#'   \item computes surface-fire heating from \code{runs} and plant-flame heating from \code{IP};
#'   \item maps heating to stratum bases and tops from \code{strata(Param)};
#'   \item derives burn and scorch fractions for each stratum and a 0–7 severity class plus text label.
#'   }
#'
#' @param runs A data frame (one or more rows per \code{repId}) describing the
#'   surface fire and site conditions. Required columns include:
#'   \itemize{
#'     \item \code{repId} (integer or character id),
#'     \item \code{slope_degrees} (numeric, degrees),
#'     \item \code{temperature} (numeric),
#'     \item \code{wind_kph} (numeric),
#'     \item \code{extinct} (numeric attenuation factor; \code{NA} is treated as 1),
#'     \item \code{lengthSurface} (numeric),
#'     \item \code{angleSurface} (numeric, radians).
#'   }
#'
#' @param IP A data frame of plant flame segments (zero or more rows per
#'   \code{repId}) with at least:
#'   \itemize{
#'     \item \code{repId},
#'     \item \code{level} (character; one of "NearSurface","Elevated","MidStorey","Canopy"),
#'     \item \code{x0}, \code{y0}, \code{x1}, \code{y1} (numeric segment endpoints),
#'     \item \code{length} (numeric segment length),
#'     \item \code{flameLength} (numeric).
#'   }
#'   Rows where all numeric fields are missing are ignored.
#'
#' @param Param Parameter table used by \code{strata(Param)} to derive stratum
#'   bounds. \code{strata(Param)} must return a data frame with at least:
#'   \itemize{
#'     \item \code{name} (character; "near surface","elevated","midstorey","canopy"),
#'     \item \code{base} (numeric base height),
#'     \item \code{top} (numeric top height).
#'   }
#'
#' @param Test Numeric “target” temperature used in the heating reach
#'   calculations. Default is \code{70}.
#'
#' @return A data frame with one row per \code{repId} containing:
#' \itemize{
#'   \item \code{repId}, \code{wind_kph},
#'   \item \code{Height} (final heating height = max of surface and plant heating),
#'   \item attenuated burn heights: \code{ns}, \code{e}, \code{m}, \code{c},
#'   \item burn fractions as rounded integers: \code{b1} (near surface), \code{b2} (elevated),
#'         \code{b3} (midstorey), \code{b4} (canopy),
#'   \item scorch fractions as rounded integers: \code{sc1}, \code{sc2}, \code{sc3}, \code{sc4},
#'   \item \code{severity} (integer in 0–7),
#'   \item \code{sevDesc} (character label for the severity class).
#' }
#'
#' @section Severity classes:
#' \itemize{
#'   \item 7: canopy burn greater than 0.5,
#'   \item 6: canopy burn greater than 0 up to 0.5,
#'   \item 5/4/3: canopy scorch high/moderate/low,
#'   \item 2/1: combined near-surface plus elevated burn greater than 0.5 or greater than 0,
#'   \item 0: no damage.
#' }
#'
#' @export
#'
#' @importFrom dplyr as_tibble mutate select group_by summarise transmute right_join left_join
#' @importFrom dplyr case_when across if_any n_distinct pull coalesce


frameSeverity <- function(runs, IP, Param = Param, Test = 70)
{
  # ---- Helpers ----
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  # tiny eps to guard against zero denominators; tolerance to avoid "numerical dust"
  eps <- .Machine$double.eps
  tol <- 1e-9
  
  # Precompute per-replicate trig and safe fields
  runs <- dplyr::as_tibble(runs) %>%
    dplyr::mutate(
      slope_rad = slope_degrees * (pi / 180),
      slope_tan = tan(slope_rad),
      cos_angleSurface = cos(angleSurface),
      sin_angleSurface = sin(angleSurface),
      extinct = dplyr::coalesce(extinct, 1)
    )
  
  # Strata parameters
  S  <- strata(Param)
  NS <- S[S$name == "near surface", , drop = FALSE]
  E  <- S[S$name == "elevated",     , drop = FALSE]
  M  <- S[S$name == "midstorey",    , drop = FALSE]
  C  <- S[S$name == "canopy",       , drop = FALSE]
  
  nsTop <- if (nrow(NS) == 0) 0 else NS$top[1]; nsBase <- if (nrow(NS) == 0) 0 else NS$base[1]
  eTop  <- if (nrow(E ) == 0) 0 else E$top[1];  eBase  <- if (nrow(E ) == 0) 0 else E$base[1]
  mTop  <- if (nrow(M ) == 0) 0 else M$top[1];  mBase  <- if (nrow(M ) == 0) 0 else M$base[1]
  cTop  <- if (nrow(C ) == 0) 0 else C$top[1];  cBase  <- if (nrow(C ) == 0) 0 else C$base[1]
  
  # Helper: get per-level max y1, aligned to Surf$repId, skipping empty/all-NA rows
  level_height <- function(IP, level_name, runs) {
    IP %>%
      dplyr::filter(level == level_name) %>%
      dplyr::filter(dplyr::if_any(dplyr::where(is.numeric), ~ !is.na(.x))) %>%
      dplyr::group_by(repId) %>%
      dplyr::summarise(y1 = safe_max(y1), .groups = "drop") %>%
      dplyr::transmute(repId, h = y1) %>%
      dplyr::right_join(dplyr::select(runs, repId), by = "repId") %>%
      dplyr::mutate(h = dplyr::coalesce(h, 0)) %>%
      dplyr::pull(h)
  }
  
  # ---- Per-level heights ----
  c  <- level_height(IP, "Canopy",      runs)
  m  <- level_height(IP, "MidStorey",   runs)
  e  <- level_height(IP, "Elevated",    runs)
  ns <- level_height(IP, "NearSurface", runs)
  
  # Cascade: if an upper stratum burns, force lower to full height (with tolerance)
  m  <- ifelse(c  > tol, mTop, m)
  e  <- ifelse(m  > tol, eTop, e)
  ns <- ifelse(e  > tol, nsTop, ns)
  
  # ---- Surface heating ----
  ground <- dplyr::tibble(repId = runs$repId) %>%
    dplyr::left_join(runs, by = "repId") %>%
    dplyr::mutate(
      Alpha   = 1 / (2 * pmax(lengthSurface^2, eps)),
      C       = 950 * lengthSurface * exp(-Alpha * lengthSurface^2),
      pAlpha  = abs(C / pmax(abs(Test - temperature), eps)),
      R       = pAlpha * cos_angleSurface,
      El      = R * slope_tan,
      ht      = (sin_angleSurface * pAlpha - El) * extinct,
      ns      = ns * extinct,
      e       = e  * extinct,
      m       = m  * extinct,
      c       = c  * extinct
    ) %>%
    dplyr::select(repId, extinct, temperature, slope_degrees, wind_kph, ht, ns, e, m, c)
  
  # ---- Plant heating ----
  IP <- IP %>%
    dplyr::filter(dplyr::if_any(dplyr::where(is.numeric), ~ !is.na(.x))) %>%
    dplyr::left_join(ground, by = "repId") %>%
    dplyr::mutate(
      Angle = atan2((y1 - y0), (x1 - x0)),
      den   = pmax(abs(flameLength * (flameLength - length)), eps),
      Alpha = 1 / (2 * den),
      C     = 950 * flameLength * exp(-Alpha * (flameLength - length)^2),
      pAlpha = abs(C / pmax(abs(Test - temperature), eps)),
      pAlpha = pmin(pAlpha, 1e6),  # clamp extreme values (optional, stabilises outliers)
      Reach = pAlpha * cos(Angle),
      El    = Reach * (tan(slope_degrees * (pi/180))),  # uses slope_degrees directly here (equivalent to slope_tan)
      htP   = (sin(Angle) * pAlpha + y0 - El) * extinct
    ) %>%
    dplyr::group_by(repId) %>%
    dplyr::summarise(htP = safe_max(htP), .groups = "drop") %>%
    dplyr::right_join(ground, by = "repId") %>%
    dplyr::mutate(
      htP = dplyr::coalesce(htP, 0)
    )
  
  # Sanitize any remaining non-finite numbers before severity
  IP <- IP %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ dplyr::if_else(is.finite(.x), .x, 0)))
  
  # ---- Severity mapping (use span guards to avoid 0/0) ----
  nsSpan <- max(nsTop - nsBase, 0)
  eSpan  <- max(eTop  - eBase,  0)
  mSpan  <- max(mTop  - mBase,  0)
  cSpan  <- max(cTop  - cBase,  0)
  
  Iso <- IP %>%
    dplyr::mutate(
      Height = pmax(ht, htP),
      
      b1 = ifelse(nsSpan == 0, 0, round(pmin(100, pmax(0, 100 * (ns     - nsBase) / nsSpan)), 0)),
      b2 = ifelse(eSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (e      - eBase)  / eSpan )), 0)),
      b3 = ifelse(mSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (m      - mBase)  / mSpan )), 0)),
      b4 = ifelse(cSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (c      - cBase)  / cSpan )), 0)),
      
      sc1 = ifelse(nsSpan == 0, 0, round(pmin(100, pmax(0, 100 * (Height - nsBase) / nsSpan)), 0)),
      sc2 = ifelse(eSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (Height - eBase)  / eSpan )), 0)),
      sc3 = ifelse(mSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (Height - mBase)  / mSpan )), 0)),
      sc4 = ifelse(cSpan  == 0, 0, round(pmin(100, pmax(0, 100 * (Height - cBase)  / cSpan )), 0)),
      
      severity = dplyr::case_when(
        b4 > 50 ~ 7L,
        b4 > 0  ~ 6L,
        sc4 > 80 ~ 5L,
        sc4 > 20 ~ 4L,
        sc4 > 0  ~ 3L,
        b1 + b2 > 50 ~ 2L,
        b1 + b2 > 0  ~ 1L,
        TRUE ~ 0L
      ),
      
      sevDesc = dplyr::case_when(
        b4 > 50 ~ "Crown fire (>0.5)",
        b4 > 0  ~ "Crown fire (0-0.5)",
        sc4 > 80 ~ "High scorch (>0.8)",
        sc4 > 20 ~ "Moderate scorch (0.2-0.8)",
        sc4 > 0  ~ "Low scorch (0-0.2)",
        b1 + b2 > 50 ~ "Understorey fire (>0.5)",
        b1 + b2 > 0  ~ "Understorey fire (0-0.5)",
        TRUE ~ "No damage"
      )
    ) %>%
    dplyr::select(repId, wind_kph, Height, ns, e, m, c,
                  b1, b2, b3, b4, sc1, sc2, sc3, sc4, severity, sevDesc)
  
  # Final sanity: ensure one row per repId
  stopifnot(dplyr::n_distinct(Iso$repId) == nrow(Iso))
  
  return(Iso)
}



#' @title cambium
#' @description Finds radial bole necrosis depth
#'
#' Depth to which the cambium of a tree recieves lethal heating
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
#' Heat is transferred into the bark and timber using Fourier's Law
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
#' Continues heating of the bole in the wake of the front for the duration of the surface fire for a period determined using
#' Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
#' Int. J. Wildl. Fire 10, 137–143 (2001). Flame lengths are decreased exponentially over this period
#'
#' Bark is assumed to ignite, and burn for an average of resBark seconds, with the default value of 45s used as a mean for
#' figure 6c in Penman, T. D., Cawson, J. G., Murphy, S. & Duff, T. J.
#' Messmate Stringybark: bark ignitability and burning sustainability in relation to fragment dimensions,
#' hazard and time since fire. Int. J. Wildl. Fire 26, 866–876 (2017).
#'
#' Heats an area of 0.01m2
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Plant The dataframe produced by the function 'repFlame'.
#' @param percentile defines which heating statistics are used for each second, from 0 (min) to 1 (max)
#' @param Height The height on the bole directly over ground (m)
#' @param woodDensity The density of wood in the tree or log housing the hollow (kg/m3)
#' @param barkDensity The density of bark in the tree or log housing the hollow (kg/m3)
#' @param comBark Temperature directly under the burning bark (C)
#' @param bark The thickness of bark on the thinnest side of the hollow (m)
#' @param resBark Flame residence in the tree bark (s)
#' @param phloem Thickness of the tree phloem (depth to cambium, m)
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the wood
#' @param bMoisture The proportion oven-dry weight of moisture in the bark
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail The number of seconds to continue modelling after all flames have extinguished
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param diameter depth of the litter layer (mm)
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param startTemp The starting temperature of wood and bark (deg C)
#' @param necT Temperature of necrosis (deg C)
#' @param surfDecl adjusts the rate at which surface flame length declines after the front
#' @param updateProgress Progress bar for use in the dashboard
#'
#' @return dataframe
#' @export

cambium <- function(Surf, Plant, percentile = 0.95, Height = 0.1, woodDensity = 700, barkDensity = 500,
                    bark = 0.04, comBark = 700, resBark = 45, phloem = 0.01, RH = 0.2,
                    moisture = 1, bMoisture = 0.2, distance = 50, trail = 100, var = 10, diameter = 20, Pressure = 1013.25,
                    Altitude = 0, startTemp = 25, necT = 60, surfDecl = 10,updateProgress=NULL)
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
  # Convert units from kg/m^3 to kg
  massB <- 0.01 * bark * barkDensity
  massW <- 0.0001 * woodDensity
  R <- sqrt(0.01/pi)
  startM <- moisture
  barkTemp <- startTemp
  woodTemp <- startTemp
  
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = 0.35 + 0.47*Re^(1/2)*0.837,
           #Incoming heat from surface
           pt = pmax(0, t-Tb),
           comBark = ifelse(pt <= resBark, comBark, 0),
           postS = bole(lengthSurface,residence, depth, h = Height,
                        surfDecl = 10, t = pt),
           tempS = ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir),
           qc = h * (tempS - barkTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
           qr = 0.86*qr*att,
           Qi = pmax(0, qc)+qr,
           
           #BARK________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = bMoisture*massB,
           # Energy removed by current water quantity
           drain = ifelse(barkTemp>95,
                          ifelse(bMoisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiA = max(Qi-drain,0),
           #Bark thermal values
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpBark = ((1105+4.85*barkTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture)/1000,
           kBark = (2.104*barkDensity+5.544*rhoM+3.266*barkTemp-166.216)*10^-4,
           # Fourier conduction
           fourierA = ifelse(Horiz>0, pmin(QiA,(kBark * pmax(0,tempS - barkTemp)) / bark),
                             (kBark * pmax(0,tempS - barkTemp)) / bark),
           barkTemp = pmax(barkTemp, (0.001 * fourierA / (massB * cpBark) + barkTemp)),
           # Change in proportion water this step
           moistureA = ifelse(bMoisture>0,max(0,bMoisture-((fourierA/2256400)/mWater)),
                              bMoisture),
           
           #1ST CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water (kg)
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiB = max(fourierA-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodB = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodB = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierB = pmin(QiB,(kWoodB * pmax(0,barkTemp - woodTemp)) / 0.01),
           woodTempB = pmax(woodTemp, (0.001 * fourierB / (massW * cpWoodB) + woodTemp)),
           # Change in proportion water this step
           moistureB = ifelse(moisture>0,max(0,moisture-((fourierB/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           barkTemp = pmax(barkTemp, barkTemp + (0.001 * fourierB / (massB * cpBark))),
           
           #2ND CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiC = max(fourierB-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodC = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodC = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierC = pmin(QiC,(kWoodC * pmax(0,woodTempB - woodTemp)) / 0.01),
           woodTempC = pmax(woodTemp, (0.001 * fourierC / (massW * cpWoodC) + woodTemp)),
           # Change in proportion water this step
           moistureC = ifelse(moisture>0,max(0,moisture-((fourierC/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempB = pmax(woodTempB, woodTempB + (0.001 * fourierC / (massW * cpWoodB))),
           
           #3RD CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiD = max(fourierC-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodD = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodD = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierD = pmin(QiD,(kWoodD * pmax(0,woodTempC - woodTemp)) / 0.01),
           woodTempD = pmax(woodTemp, (0.001 * fourierD / (massW * cpWoodD) + woodTemp)),
           # Change in proportion water this step
           moistureD = ifelse(moisture>0,max(0,moisture-((fourierD/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempC = pmax(woodTempC, woodTempC + (0.001 * fourierD / (massW * cpWoodC))),
           
           #4TH CM________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWater = moisture*massW,
           # Energy removed by current water quantity
           drain = ifelse(woodTemp>95,
                          ifelse(moisture>0,mWater*2256400,0),0),
           # Adjusted incoming energy
           QiE = max(fourierD-drain,0),
           #Thermal values
           kAir = 0.00028683*(woodTemp+273.15)^0.7919,
           cpWoodE = 1.08+0.00408*(100*moisture)+0.00253*woodTemp+0.0000628*(100*moisture)*woodTemp,
           kWoodE = kWood(woodTemp, woodDensity, kAir),
           # Fourier conduction
           fourierE = pmin(QiE,(kWoodE * pmax(0,woodTempD - woodTemp)) / 0.01),
           woodTempE = pmax(woodTemp, (0.001 * fourierE / (massW * cpWoodE) + woodTemp)),
           # Change in proportion water this step
           moistureE = ifelse(moisture>0,max(0,moisture-((fourierE/2256400)/mWater)),
                              moisture),
           # Heat draw-down from above
           woodTempD = pmax(woodTempD, woodTempD + (0.001 * fourierE / (massW * cpWoodD))),
           
           #Outgoing BARK________________________________________________________
           fourierOA = pmin(0,(kBark * (tempS - barkTemp)) / bark),
           barkTemp = pmin(barkTemp, barkTemp + (0.001 * fourierOA / (massB * cpBark))),
           qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (barkTemp+273.15)^4)),
           qR = qr + qrO,
           Q = qc + qR,
           
           
           #Outgoing 1ST________________________________________________________
           fourierOB = pmin(0,(kWoodB * (barkTemp - woodTempB)) / 0.01),
           woodTempB = pmin(woodTempB, woodTempB + (0.001 * fourierOB / (massW * cpWoodB))),
           #Outgoing 2ND________________________________________________________
           fourierOC = pmin(0,(kWoodC * (woodTempB - woodTempC)) / 0.01),
           woodTempC = pmin(woodTempC, woodTempC + (0.001 * fourierOC / (massW * cpWoodC))),
           #Outgoing 3RD________________________________________________________
           fourierOD = pmin(0,(kWoodD * (woodTempC - woodTempD)) / 0.01),
           woodTempD = pmin(woodTempD, woodTempD + (0.001 * fourierOD / (massW * cpWoodD))),
           #Outgoing 4TH________________________________________________________
           fourierOE = pmin(0,(kWoodE * (woodTempD - woodTempE)) / 0.01),
           woodTempE = pmin(woodTempE, woodTempE + (0.001 * fourierOE / (massW * cpWoodE))))
  
  
  barkTemp <- quantile(Ca$barkTemp, percentile)
  moistureA <- quantile(Ca$moistureA, percentile)
  woodTempB <- quantile(Ca$woodTempB, percentile)
  moistureB <- quantile(Ca$moistureB, percentile)
  woodTempC <- quantile(Ca$woodTempC, percentile)
  moistureC <- quantile(Ca$moistureC, percentile)
  woodTempD <- quantile(Ca$woodTempD, percentile)
  moistureD <- quantile(Ca$moistureD, percentile)
  woodTempE <- quantile(Ca$woodTempE, percentile)
  moistureE <- quantile(Ca$moistureE, percentile)
  
  # Advance one second's travel
  Horiz = Horiz - ROS
  pbar <-  txtProgressBar(max = TIME, style = 3)
  # Loop through each time step and collect outputs
  for(t in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      mutate(t = t,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = 0.35 + 0.47*Re^(1/2)*0.837,
             #Incoming heat from surface
             pt = pmax(0, t-Tb),
             comBark = ifelse(pt <= resBark, comBark, 0),
             postS = bole(lengthSurface,residence, depth, h = Height,
                          surfDecl = 10, t = pt),
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir)),
             qc = h * (tempS - barkTemp),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
             qr = 0.86*qr*att,
             Qi = pmax(0, qc)+qr,
             
             
             #BARK________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureA*massB,
             # Energy removed by current water quantity
             drain = ifelse(barkTemp>95,
                            ifelse(moistureA>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiA = max(Qi-drain,0),
             #Bark thermal values
             rhoM = (moistureA+moistureA^2)*barkDensity,
             cpBark = ((1105+4.85*barkTemp)*(1-moistureA)+moistureA*4185+1276*moistureA)/1000,
             kBark = (2.104*barkDensity+5.544*rhoM+3.266*barkTemp-166.216)*10^-4,
             # Fourier conduction
             fourierA = ifelse(Horiz>0, pmin(QiA,(kBark * pmax(0,tempS - barkTemp)) / bark),
                               (kBark * pmax(0,tempS - barkTemp)) / bark),
             barkTemp = pmax(barkTemp, (0.001 * fourierA / (massB * cpBark) + barkTemp)),
             # Change in proportion water this step
             moistureA = ifelse(moistureA>0,max(0,moistureA-((fourierA/2256400)/mWater)),
                                moistureA),
             
             
             #1ST CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureB*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempB>95,
                            ifelse(moistureB>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiB = max(fourierA-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempB+273.15)^0.7919,
             cpWoodB = 1.08+0.00408*(100*moistureB)+0.00253*woodTempB+0.0000628*(100*moistureB)*woodTempB,
             kWoodB = kWood(woodTempB, woodDensity, kAir),
             # Fourier conduction
             fourierB = pmin(QiB,(kWoodB * pmax(0,barkTemp - woodTempB)) / 0.01),
             woodTempB = pmax(woodTempB, (0.001 * fourierB / (massW * cpWoodB) + woodTempB)),
             # Change in proportion water this step
             moistureB = ifelse(moistureB>0,max(0,moistureB-((fourierB/2256400)/mWater)),
                                moistureB),
             # Heat draw-down from above
             barkTemp = pmin(barkTemp, barkTemp - (0.001 * fourierB / (massB * cpBark))),
             
             #2ND CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureC*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempC>95,
                            ifelse(moistureC>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiC = max(fourierB-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempC+273.15)^0.7919,
             cpWoodC = 1.08+0.00408*(100*moistureC)+0.00253*woodTempC+0.0000628*(100*moistureC)*woodTempC,
             kWoodC = kWood(woodTempC, woodDensity, kAir),
             # Fourier conduction
             fourierC = pmin(QiC,(kWoodC * pmax(0,woodTempB - woodTempC)) / 0.01),
             woodTempC = pmax(woodTempC, (0.001 * fourierC / (massW * cpWoodC) + woodTempC)),
             # Change in proportion water this step
             moistureC = ifelse(moistureC>0,max(0,moistureC-((fourierC/2256400)/mWater)),
                                moistureC),
             # Heat draw-down from above
             woodTempB = pmin(woodTempB, woodTempB - (0.001 * fourierC / (massW * cpWoodB))),
             
             #3RD CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureD*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempD>95,
                            ifelse(moistureD>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiD = max(fourierC-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempD+273.15)^0.7919,
             cpWoodD = 1.08+0.00408*(100*moistureD)+0.00253*woodTempD+0.0000628*(100*moistureD)*woodTempD,
             kWoodD = kWood(woodTempD, woodDensity, kAir),
             # Fourier conduction
             fourierD = pmin(QiD,(kWoodD * pmax(0,woodTempC - woodTempD)) / 0.01),
             woodTempD = pmax(woodTempD, (0.001 * fourierD / (massW * cpWoodD) + woodTempD)),
             # Change in proportion water this step
             moistureD = ifelse(moistureD>0,max(0,moistureD-((fourierD/2256400)/mWater)),
                                moistureD),
             # Heat draw-down from above
             woodTempC = pmin(woodTempC, woodTempC - (0.001 * fourierD / (massW * cpWoodC))),
             
             #4TH CM________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWater = moistureE*massW,
             # Energy removed by current water quantity
             drain = ifelse(woodTempE>95,
                            ifelse(moistureE>0,mWater*2256400,0),0),
             # Adjusted incoming energy
             QiE = max(fourierD-drain,0),
             #Thermal values
             kAir = 0.00028683*(woodTempE+273.15)^0.7919,
             cpWoodE = 1.08+0.00408*(100*moistureE)+0.00253*woodTempE+0.0000628*(100*moistureE)*woodTempE,
             kWoodE = kWood(woodTempE, woodDensity, kAir),
             # Fourier conduction
             fourierE = pmin(QiE,(kWoodE * pmax(0,woodTempD - woodTempE)) / 0.01),
             woodTempE = pmax(woodTempE, (0.001 * fourierE / (massW * cpWoodE) + woodTempE)),
             # Change in proportion water this step
             moistureE = ifelse(moistureE>0,max(0,moistureE-((fourierE/2256400)/mWater)),
                                moistureE),
             # Heat draw-down from above
             woodTempD = pmin(woodTempD, woodTempD - (0.001 * fourierE / (massW * cpWoodD))),
             
             #Outgoing 1st________________________________________________________
             fourierOA = pmin(0,(kBark * (tempS - barkTemp)) / bark),
             barkTemp = pmin(barkTemp, barkTemp + (0.001 * fourierOA / (massB * cpBark))),
             qrO = pmin(0,0.86*0.0000000000567*((tempS+273.15)^4 - (barkTemp+273.15)^4)),
             qR = qr + qrO,
             Q = qc + qR,
             
             
             #Outgoing 1ST________________________________________________________
             fourierOB = pmin(0,(kWoodB * (barkTemp - woodTempB)) / 0.01),
             woodTempB = pmin(woodTempB, woodTempB + (0.001 * fourierOB / (massW * cpWoodB))),
             #Outgoing 2ND________________________________________________________
             fourierOC = pmin(0,(kWoodC * (woodTempB - woodTempC)) / 0.01),
             woodTempC = pmin(woodTempC, woodTempC + (0.001 * fourierOC / (massW * cpWoodC))),
             #Outgoing 3RD________________________________________________________
             fourierOD = pmin(0,(kWoodD * (woodTempC - woodTempD)) / 0.01),
             woodTempD = pmin(woodTempD, woodTempD + (0.001 * fourierOD / (massW * cpWoodD))),
             #Outgoing 4TH________________________________________________________
             fourierOE = pmin(0,(kWoodE * (woodTempC - woodTempE)) / 0.01),
             woodTempE = pmin(woodTempE, woodTempD + (0.001 * fourierOE / (massW * cpWoodE))))
    
    Ca <- rbind(Ca, Cb)
    
    barkTemp <- quantile(Cb$barkTemp, percentile)
    moistureA <- quantile(Cb$moistureA, percentile)
    woodTempB <- quantile(Cb$woodTempB, percentile)
    moistureB <- quantile(Cb$moistureB, percentile)
    woodTempC <- quantile(Cb$woodTempC, percentile)
    moistureC <- quantile(Cb$moistureC, percentile)
    woodTempD <- quantile(Cb$woodTempD, percentile)
    moistureD <- quantile(Cb$moistureD, percentile)
    woodTempE <- quantile(Cb$woodTempE, percentile)
    moistureE <- quantile(Cb$moistureE, percentile)
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
    select(t, repId, tempS, barkTemp, woodTempB, woodTempC, woodTempD, woodTempE,
           moistureA, moistureB, moistureC, moistureD, moistureE)%>%
    mutate(necrosis = ifelse(woodTempE>=necT, 4,
                             ifelse(woodTempD>=necT, 3,
                                    ifelse(woodTempC>=necT, 2,
                                           ifelse(woodTempB>=necT, 1, 0)))),
           girdling = ifelse(necrosis>=phloem, 1, 0))
  
  return(Ca)
}



#' @title Girdling
#' @description Finds radial bole necrosis depth, dividing bark into four steps
#'
#' Depth to which the stem of a plant recieves lethal heating
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
#' Heat is transferred into the bark and timber using Fourier's Law
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
#' Continues heating of the bole in the wake of the front for the duration of the surface fire for a period determined using
#' Burrows, N. D. Flame residence times and rates of weight loss of eucalypt forest fuel particles.
#' Int. J. Wildl. Fire 10, 137–143 (2001). Flame lengths are decreased exponentially over this period
#'
#' Bark is assumed to ignite, and burn for an average of resBark seconds, with the default value of 45s used as a mean for
#' figure 6c in Penman, T. D., Cawson, J. G., Murphy, S. & Duff, T. J.
#' Messmate Stringybark: bark ignitability and burning sustainability in relation to fragment dimensions,
#' hazard and time since fire. Int. J. Wildl. Fire 26, 866–876 (2017).
#'
#' Heats an area of 0.01m2
#'
#'
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Plant The dataframe produced by the function 'repFlame'
#' @param Height The height on the bole directly over ground (m)
#' @param woodDensity The density of wood in the plant or log housing the hollow (kg/m3)
#' @param barkDensity The density of bark in the plant or log housing the hollow (kg/m3)
#' @param comBark Temperature directly under the burning bark (C)
#' @param bark The thickness of bark on the thinnest side of the hollow (m)
#' @param resBark Flame residence in the plant bark (s)
#' @param phloem Thickness of the plant phloem (depth to cambium, m)
#' @param RH The relative humidity (0-1)
#' @param moisture The proportion oven-dry weight of moisture in the wood
#' @param bMoisture The proportion oven-dry weight of moisture in the bark
#' @param distance The furthest horizontal distance between the flame origin and the point (m)
#' @param trail The number of seconds to continue modelling after all flames have extinguished
#' @param var The angle in degrees that the plume spreads above/below a central vector
#' @param diameter depth of the litter layer (mm)
#' @param Pressure Sea level atmospheric pressure (hPa)
#' @param Altitude Height above sea level (m)
#' @param startTemp The starting temperature of wood and bark (deg C)
#' @param necT Temperature of necrosis (deg C)
#' @param surfDecl adusts the rate at which surface flame length declines after the front
#' @param updateProgress Progress bar for use in the dashboard
#'
#' @return dataframe
#' @export
#'

girdle <- function(Surf, Plant, Height = 0.1, woodDensity = 700, barkDensity = 500,
                    bark = 0.04, comBark = 700, resBark = 45, phloem = 0.01, RH = 0.2,
                    moisture = 0.6, bMoisture = 0.5, distance = 5, trail = 100, var = 10, 
                    diameter = 20, Pressure = 1013.25,Altitude = 0, startTemp = 25, 
                    necT = 60, surfDecl = 10,updateProgress=NULL)
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
  step <- bark/4
  # Convert units from kg/m^3 to kg
  mass <- step * barkDensity
  massW <- phloem * woodDensity
  R <- sqrt(0.01/pi)
  startM <- moisture
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
    summarise_all(mean)%>%
    mutate(t = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = 0.35 + 0.47*Re^(1/2)*0.837,
           #Incoming heat from surface
           pt = pmax(0, t-Tb),
           comBark = ifelse(pt <= resBark, comBark, 0),
           postS = bole(lengthSurface,residence, depth, h = Height,
                        surfDecl = 10, t = pt),
           tempS = ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir),
           qc = h * (tempS - startTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpA = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kA = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fAD = ((kA * (tempS - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpB = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kB = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fBD = ((kB * (tempA - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpC = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kC = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fCD = ((kC * (tempB - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpD = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kD = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fDD = ((kD * (tempC - startTemp)) / step),
           fDU = 0,
           fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
           tempD = (fourierD / (mass * cpD) + startTemp),
           # Change in proportion water this step
           moistureD = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((fourierC/2256400)/mWaterD)),
                                                  bMoisture),bMoisture),
           
           #Phloem ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterE = moisture*massW,
           # Energy removed by current water quantity
           drainE = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterE*2256400,0),0),
           #Thermal values - (cp: J/g/deg, k: W/m/deg)
           kAir = 0.00028683*(startTemp+273.15)^0.7919,
           cpE = 1080+408*moisture+2.53*startTemp+6.28*moisture*startTemp,
           kE = kWood(startTemp, woodDensity, kAir),
           # Conduction from above and below, less latent heat of evaporation
           fED = ((kE * (tempD - startTemp)) / step),
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
             Re = (Plume_velocity*Density)/viscosity,
             h = 0.35 + 0.47*Re^(1/2)*0.837,
             #Incoming heat from surface
             pt = pmax(0, t-Tb),
             comBark = ifelse(pt <= resBark, comBark, 0),
             postS = bole(lengthSurface,residence, depth, h = Height,
                          surfDecl = 10, t = pt),
             tempS = ifelse(t>Ta, tempAir, ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir)),
             qc = h * (tempS - tempA),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
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
             rhoM = (moistureA+moistureA^2)*barkDensity,
             cpA = ((1105+4.85*tempA)*(1-moistureA)+moistureA*4185+1276*moistureA),
             kA = (2.104*barkDensity+5.544*rhoM+3.266*tempA-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fAD = ((kA * (tempS - tempA)) / step),
             fAU = ((kB * (tempB - tempA)) / step),
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
             rhoM = (moistureB+moistureB^2)*barkDensity,
             cpB = ((1105+4.85*tempB)*(1-moistureB)+moistureB*4185+1276*moistureB),
             kB = (2.104*barkDensity+5.544*rhoM+3.266*tempB-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fBD = ((kB * (tempA - tempB)) / step),
             fBU = ((kC * (tempC - tempB)) / step),
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
             rhoM = (moistureC+moistureC^2)*barkDensity,
             cpC = ((1105+4.85*tempC)*(1-moistureC)+moistureC*4185+1276*moistureC),
             kC = (2.104*barkDensity+5.544*rhoM+3.266*tempC-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fCD = ((kC * (tempB - tempC)) / step),
             fCU = ((kC * (tempD - tempC)) / step),
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
             rhoM = (moistureD+moistureD^2)*barkDensity,
             cpD = ((1105+4.85*tempD)*(1-moistureD)+moistureD*4185+1276*moistureD),
             kD = (2.104*barkDensity+5.544*rhoM+3.266*tempD-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fDD = ((kD * (tempC - tempD)) / step),
             fDU = ((kD * (tempE - tempD)) / step),
             fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
             tempD = (fourierD / (mass * cpD) + tempD),
             # Change in proportion water this step
             moistureD = ifelse(tempD>99,ifelse(moistureD>0,max(0,moistureD-((fourierC/2256400)/mWaterD)),
                                                moistureD),moistureD),
             
             #phloem ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterE = moistureE*massW,
             # Energy removed by current water quantity
             drainE = ifelse(tempE>99,
                             ifelse(moistureE>0,mWaterE*2256400,0),0),
             #Thermal values - (cp: J/g/deg, k: W/m/deg)
             kAir = 0.00028683*(tempE+273.15)^0.7919,
             cpE = 1080+408*moistureE+2.53*tempE+6.28*moistureE*tempE,
             kE = kWood(tempE, woodDensity, kAir),
             # Conduction from above and below, less latent heat of evaporation. Below unknown.
             fED = ((kE * (tempD - tempE)) / phloem),
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
      text <- paste0("Number of remaining steps is ", TIME - t)
      updateProgress(detail = text)
    }
    t <- t + 1
    Horiz <- Horiz - ROS
  }
  
  # Create table
  Ca <- Ca %>%
    select(t, repId, tempS, tempA, tempB, tempC, tempD, tempE,
           moistureA, moistureB, moistureC, moistureD, moistureE)%>%
    mutate(girdling = ifelse(tempE>=necT, 1, 0))
  
  return(Ca)
}



#####################################################################

#' @title bole
#' @description Air temperature above ambient at the tree bole behind the flame front
#'
#' Dynamic air temperature at bole height, declining flame length exponentially
#'
#' Air temperature is modelled from dynamic flame segments using
#' Weber R.O., Gill A.M., Lyons P.R.A., Moore P.H.R., Bradstock R.A., Mercer G.N. (1995)
#' Modelling wildland fire temperatures. CALMScience Supplement, 4, 23–26.
#'
#' pAlphas is set to bole height - depth of surface litter
#'
#' @param lengthSurface Length of the surface flame (m)
#' @param residence Residence time of the flame (s)
#' @param depth Depth of the surface litter (m)
#' @param h Height of the bole above ground (m)
#' @param surfDecl Rate of decline of the surface flame length (m/s)
#' @param t Time since the flame front passed (s)
#'
#' @return value
#' @export
#'

bole <- function(lengthSurface = 2, residence = 300, depth = 0.05, h = 0.1, surfDecl = 10, t = 1)
{
  surfPost <- lengthSurface*exp(-(surfDecl/residence)*t)
  Alpha <- 1/(2 * surfPost^2)
  C <- 950 * surfPost * exp(-Alpha * surfPost^2)
  pAlphas <- pmax(0,h-depth)
  return(ifelse(pAlphas < surfPost,
                950 + exp(-Alpha * pAlphas^2),
                C/pAlphas))
}



#' @title kWood
#' @description Thermal conductivity of dry wood
#'
#' Model drawn from Kollmann, F. F. P. & Cote, W. A.
#' Principles of wood science and technology I. Solid wood. (Springer-Verlag, 1968)
#'
#' @param rhoW Density of the wood (kg/m3)
#' @param kAir Thermal conductivity of air (W/m/deg)
#' @param T Temperature of the wood (deg C)
#'
#' @return value
#' @export
#'

kWood <- function(T=100, rhoW=700, kAir = 0.026)
{
  T <- T+273.15
  bridge <- 0.0019*T+0.0503
  r <- 1-(rhoW/1500)
  kA <- (1-r)*0.766+r*kAir
  kB <- 1/((1-r)/0.43+(r/kAir))
  return(bridge*kA+(1-bridge)*kB)
}


#' @title plantFlame
#' @description Collects flame segments for a specified plant in a stratum
#' For each segment, points A-C mark the origin, end of burning foliage, 
#' and tip of the flame segment, respectively.
#'
#' @param paths Output table from the function repFlame
#' @param Stratum Name of the stratum being studied
#' @param Species Name of the species being studied
#' @param repId Number of the repId being studied
#' @param site Name of the site being studied
#'
#' @return dataframe
#' @export

plantFlame <- function(paths, Stratum, Species, repId, site) {
  
  IP <- paths[paths$level == Stratum & paths$species == Species & paths$repId == repId & paths$site == site,] %>%
    mutate(angle = atan((y1-y0)/(x1-x0)),
           x2 = flameLength * cos(angle) + x0,
           y2 = flameLength * sin(angle) + y0)
  
  x0 <- IP %>%
    dplyr::select(segIndex, x0)
  x1 <- IP %>%
    dplyr::select(segIndex, x1)
  IPX <- left_join(x0, x1)
  x2 <- IP %>%
    dplyr::select(segIndex, x2)
  IPX <- left_join(IPX,x2)
  X <- reshape2::melt(IPX, id.vars="segIndex") %>%
    mutate(x = value,
           point = case_when(variable == "x0" ~ "A",
                             variable == "x1" ~ "B",
                             variable == "x2" ~ "C")) %>%
    dplyr::select(segIndex, point, x)
  
  
  y0 <- IP %>%
    dplyr::select(segIndex, y0)
  y1 <- IP %>%
    dplyr::select(segIndex, y1)
  IPy <- left_join(y0, y1)
  y2 <- IP %>%
    dplyr::select(segIndex, y2)
  IPy <- left_join(IPy,y2)
  Y <- reshape2::melt(IPy, id.vars="segIndex") %>%
    mutate(y = value,
           point = case_when(variable == "y0" ~ "A",
                             variable == "y1" ~ "B",
                             variable == "y2" ~ "C")) %>%
    dplyr::select(segIndex, point, y)
  
  XY <- left_join(X,Y, by = c("segIndex", "point"))
  out <- XY[order(XY$segIndex),]
  
  return(out)
}


#' @title LAIp
#' @description Calculates LAI for a slice of a plant
#'
#' @param base.params Parameter input file 
#' @param sp Number of the species
#' @param yu Height at the top of the slice (m)
#' @param yl Height at the base of the slice (m)
#'
#' @return list
#' @export 
#'

LAIp <- function(base.params, sp = 1, yu = 100, yl = 0)
{
  
  # Subset species
  spPar <- subset(base.params, species == sp)
  
  # Collect parameters
  Cd <- as.numeric(spPar$value[spPar$param == "clumpDiameter"])
  Cs <- as.numeric(spPar$value[spPar$param == "clumpSeparation"])
  ram <- as.numeric(spPar$value[spPar$param == "stemOrder"])
  hc <- as.numeric(spPar$value[spPar$param == "hc"])
  he <- as.numeric(spPar$value[spPar$param == "he"])
  ht <- as.numeric(spPar$value[spPar$param == "ht"])
  hp <- as.numeric(spPar$value[spPar$param == "hp"])
  W <- as.numeric(spPar$value[spPar$param == "w"])
  ll <- as.numeric(spPar$value[spPar$param == "leafLength"])
  lw <- as.numeric(spPar$value[spPar$param == "leafWidth"])
  ls <- as.numeric(spPar$value[spPar$param == "leafSeparation"])
  
  if (yu <= min(hc,he)) {
    l <- 0
    LA <- 0
  } else if (yl >= max(ht,hp)) {
    l <- 0
    LA <- 0
  } else {
    
    # Find w at cutoff cones
    topRise <- abs(hp-ht)
    baseFall <- abs(he-hc)
    # Top of slice, plant top
    wa <- ifelse(topRise==0,
                 W,
                 (W/topRise)*max(0,(hp-max(ht,yu))))
    # Base of slice, plant top
    wb <- ifelse(topRise==0,
                 W,
                 (W/topRise)*abs(hp-max(yl,ht)))
    # Top of slice, plant base
    wc <- ifelse(baseFall==0,
                 W,
                 (W/baseFall)*abs(min(he,yu)-hc))
    # Base of slice, plant base
    wd <- ifelse(baseFall==0,
                 W,
                 (W/baseFall)*max(0,(min(he,yl)-hc)))
    
    # Leaf area per clump
    clumpLeaves <- 0.88*(Cd*ram/ls)^1.18
    laClump <- (ll*lw/2) * clumpLeaves
    
    # Plant volume in slice
    upperAboveSlice <- max(0,(hp-max(ht,yu))*(pi*(wa/2)^2))/3
    upperFull <- ((hp-max(yl,ht))*(pi*(wb/2)^2))/3
    upper <- max(0,upperFull-upperAboveSlice)
    centre <- max(0,(min(yu,ht)-max(yl,he)))*(pi*(W/2)^2)
    lowerFull <- ((min(yu,he)-hc)*(pi*(wc/2)^2))/3
    lowerBelowSlice <- (max(0,(min(he,yl)-hc))*(pi*(wd/2)^2))/3
    lower <- max(0,lowerFull-lowerBelowSlice)
    vol <- upper+centre+lower
    
    # Clumps in slice
    volC <- ((4/3)*pi*((Cd+Cs)/2)^3)
    nClumps <- vol/volC
    
    #LAI per plant in slice
    LA <- laClump * nClumps
    l<- ifelse(max(wa,wb,wc,wd)==0,
               0,
               LA/(pi*(max(wa,wb,wc,wd)/2)^2))
  }
  return(list(l,LA))
}

#' Community Leaf Area Index
#'
#' Computes the community-level leaf area index (LAI) for a vertical slice,
#' using vectorized lookups from `species(base.params)` and `strata(base.params)`
#' and per-species `LAIp()` values. 
#'
#' @param base.params A parameter table used by the FRaME workflow. Must be
#'   compatible with `species()`, `strata()`, and `LAIp()`. It should contain
#'   per-species rows and a `param` and `value` column for species parameters.
#' @param yu Numeric. Upper height of the slice (same units as species heights).
#' @param yl Numeric. Lower height of the slice (same units as species heights).
#'
#' @return Numeric scalar. The total community LAI for the slice.
#'
#' @details
#' For each species, `LAIp()` is evaluated for the slice `[yl, yu]`. Species
#' weights are taken from `species(base.params)$comp`, and cover is taken from
#' `strata(base.params)$cover[species$st]`. The returned value is the sum of
#' `LAIp * Cover * Weight` across species.
#'
#' @examples
#' \dontrun{
#' # Example with user project structures:
#' LAIcomm(base.params, yu = 10, yl = 8)
#' }
#'
#' @export

LAIcomm <- function(base.params, yu = 100, yl = 0) {
  # Extract common tables
  spec <- species(base.params)
  str  <- strata(base.params)
  
  N <- nrow(spec)
  ids <- seq_len(N)
  
  # Compute per-species LAIp
  LAIp_vec <- vapply(ids, function(n) {
    (LAIp(base.params, sp = n, yu = yu, yl = yl))[[1]]
  }, numeric(1))
  
  # Lookup vectors
  Weight <- as.numeric(spec$comp)
  Cover  <- as.numeric(str$cover[spec$st])
  Stratum <- as.numeric(spec$st)
  Sep <- as.numeric(str$separation[Stratum])
  w <- vapply(ids, function(n) {
    spPar <- base.params[base.params$species == n, ]
    as.numeric(spPar$value[spPar$param == "w"][1])
  }, numeric(1))
  
  Weight[LAIp_vec <= 0] <- 0
  LAIw <- LAIp_vec * Cover * Weight
  LAI <- sum(LAIw, na.rm = TRUE)
  if (is.nan(LAI)) LAI <- 0
  LAI
}


#' Vertical wind attenuation profile by canopy slices
#'
#' Builds a height profile by slicing the canopy and accumulating a
#' multiplicative wind factor across slices.
#'
#' @param base.params A parameter table used by the FRaME workflow. Must be
#'   compatible with `species()`, `strata()`, and `LAIcomm()`.
#' @param slices Integer. Number of vertical slices between ground and
#'   the top height.
#'
#' @return A data frame with one row per slice plus an initial baseline row,
#'   containing:
#' \itemize{
#'   \item \code{l}: LAI per slice (first row is zero).
#'   \item \code{gam}: Slice-level attenuation coefficient (first row is zero).
#'   \item \code{w}: Cumulative wind multiplier starting at 1.
#'   \item \code{Slice}: Slice index starting at 1.
#'   \item \code{z}: Relative height from 1 (top) to 0 (bottom).
#'   \item \code{hm}: Absolute height corresponding to \code{z}.
#' }
#'
#' @examples
#' \dontrun{
#' prof <- profileDet_fast(base.params, slices = 20)
#' head(prof)
#' }
#'
#' @export

profileDet <- function(base.params, slices = 10) {
  stopifnot(slices >= 1L)
  
  # One-time values
  top   <- max(species(base.params)$hp)
  step  <- top / slices
  
  yu <- top - (0:(slices - 1)) * step
  yl <- yu - step
  
  # LAI for each slice
  LAIslice <- vapply(seq_len(slices),
                     function(i) LAIcomm(base.params = base.params, yl = yl[i], yu = yu[i]),
                     numeric(1))
  
  g <- 1.785 * (LAIslice ^ 0.372)
  mult <- exp(g * ((yl / yu) - 1))
  W    <- c(1, cumprod(mult))
  l_vec   <- c(0, LAIslice)
  gam_vec <- c(0, g)
  
  Slice <- seq_len(slices + 1L)
  wind <- data.frame(
    l     = l_vec,
    gam   = gam_vec,
    w     = W,
    Slice = Slice
  )
  
  wind$z  <- 1 - ((wind$Slice - 1) * (1 / slices))
  wind$hm <- wind$z * top
  
  wind
}


#' @title windReduction
#' @description Finds the Wind Reduction Factor for a param file
#'
#' @param base.params Parameter input file 
#' @param test Height at which the WRF is calculated (m)
#'
#' @return value
#' @export
#'

windReduction <- function(base.params, test = 1.2)
{
  s <- strata(base.params)
  t <- max(s$top, na.rm = TRUE)
  
  # Short-circuit: vegetation too low / undefined → no wind reduction
  if (!is.finite(t) || t <= 0 || test > t) {
    return(1)
  }
  
  slice <- max(1L, round(t / test))
  det   <- profileDet(base.params, slices = slice)
  w     <- det[nrow(det) - 1, ]
  wrf   <- as.numeric(round(1 / w$w[1], 1))
  
  return(wrf)
}



#' @title dryside
#' @description Calculates likelihood of basal scarring on a tree
#' 
#' Assumes lee-side vortex at bole causes flame to lean backward to a distance 2.5*DBH
#' Estimate taken from Malcolm Gill A 1974 
#' Toward an understanding of fire scar formation: field observation and laboratory simulation 
#' For. Sci. 20 198–205
#'
#' @param Surf The dataframe produced by the function 'summary',
#' @param Plant The dataframe produced by the function 'repFlame'.
#' @param DBH Diameter at breast height (m)
#' @param Height Height of the tree (m)
#' @param woodDensity Density of the wood (kg/m^3)
#' @param barkDensity Density of the bark (kg/m^3)
#' @param bark Thickness of the bark (m)
#' @param comBark Conductivity of the bark (W/m/deg)
#' @param resBark Resistance of the bark (deg)
#' @param phloem Thickness of the phloem (m)
#' @param RH Relative humidity
#' @param moisture Moisture content of the wood
#' @param bMoisture Moisture content of the bark
#' @param distance Distance from the flame front (m)
#' @param trail Distance from the flame front to the tree (m)
#' @param var Variability of the flame front
#' @param diameter Diameter of the flame front (m)
#' @param Pressure Atmospheric pressure (hPa)
#' @param Altitude Altitude of the tree (m)
#' @param startTemp Starting temperature (deg C)
#' @param necT Temperature at which the tree is necrotic (deg C)
#' @param surfDecl Decline of the surface temperature (deg C)
#' @param updateProgress Function to update progress
#'
#' @return dataframe
#'

dryside <- function(Surf, Plant, DBH = 1, Height = 0.1, woodDensity = 700, barkDensity = 500,
                    bark = 0.04, comBark = 700, resBark = 45, phloem = 0.01, RH = 0.2,
                    moisture = 0.6, bMoisture = 0.5, distance = 5, trail = 100, var = 10, 
                    diameter = 20, Pressure = 1013.25,Altitude = 0, startTemp = 25, 
                    necT = 60, surfDecl = 10,updateProgress=NULL)
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
  
  # Duration of post-front vortex
  vortex <- round(2.5*DBH/ROS,0)
  
  # Description of the protection
  step <- bark/4
  # Convert units from kg/m^3 to kg
  mass <- step * barkDensity
  massW <- phloem * woodDensity
  R <- sqrt(0.01/pi)
  startM <- moisture
  
  #Starting values
  Ca <- threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
    summarise_all(mean)%>%
    mutate(tS = 1,
           #Convective transfer
           Re = (Plume_velocity*Density)/viscosity,
           h = 0.35 + 0.47*Re^(1/2)*0.837,
           #Incoming heat from surface
           Pt = pmax(0, tS - Tb),
           comBark = ifelse(Pt <= resBark, comBark, 0),
           postS = if(Pt < vortex) {950} else {bole(lengthSurface,residence, depth, h = Height,
                                                    surfDecl = 10, t = Pt)},
           tempS = ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir),
           qc = h * (tempS - startTemp),
           att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpA = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kA = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fAD = ((kA * (tempS - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpB = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kB = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fBD = ((kB * (tempA - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpC = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kC = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fCD = ((kC * (tempB - startTemp)) / step),
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
           rhoM = (bMoisture+bMoisture^2)*barkDensity,
           cpD = ((1105+4.85*startTemp)*(1-bMoisture)+bMoisture*4185+1276*bMoisture),
           kD = (2.104*barkDensity+5.544*rhoM+3.266*startTemp-166.216)*10^-4,
           # Conduction from above and below, less latent heat of evaporation
           fDD = ((kD * (tempC - startTemp)) / step),
           fDU = 0,
           fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
           tempD = (fourierD / (mass * cpD) + startTemp),
           # Change in proportion water this step
           moistureD = ifelse(startTemp>99,ifelse(bMoisture>0,max(0,bMoisture-((fourierC/2256400)/mWaterD)),
                                                  bMoisture),bMoisture),
           
           #Phloem ________________________________________________________
           ### Water effects: evaporation and energy drain
           # Mass of water
           mWaterE = moisture*massW,
           # Energy removed by current water quantity
           drainE = ifelse(startTemp>99,
                           ifelse(moisture>0,mWaterE*2256400,0),0),
           #Thermal values - (cp: J/g/deg, k: W/m/deg)
           kAir = 0.00028683*(startTemp+273.15)^0.7919,
           cpE = 1080+408*moisture+2.53*startTemp+6.28*moisture*startTemp,
           kE = kWood(startTemp, woodDensity, kAir),
           # Conduction from above and below, less latent heat of evaporation
           fED = ((kE * (tempD - startTemp)) / step),
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
  for(tS in 2:TIME){
    Cb <-threat(Surf, Plant, Horiz, Height, var, Pressure, Altitude) %>%
      summarise_all(mean)%>%
      mutate(tS = tS,
             #Convective transfer
             Re = (Plume_velocity*Density)/viscosity,
             h = 0.35 + 0.47*Re^(1/2)*0.837,
             #Incoming heat from surface
             Pt = pmax(0, tS - Tb),
             comBark = ifelse(Pt <= resBark, comBark, 0),
             postS = if(Pt < vortex) {950} else {bole(lengthSurface,residence, depth, h = Height,
                                                      surfDecl = 10, t = Pt)},
             tempS = ifelse(tS > Ta, tempAir, ifelse(Horiz <=0, pmax(tempAir, postS, comBark), tempAir)),
             qc = h * (tempS - tempA),
             att = tau(D=Horiz, flameTemp=flameTemp, temperature=(temperature+273.15), rh=RH),
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
             rhoM = (moistureA+moistureA^2)*barkDensity,
             cpA = ((1105+4.85*tempA)*(1-moistureA)+moistureA*4185+1276*moistureA),
             kA = (2.104*barkDensity+5.544*rhoM+3.266*tempA-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fAD = ((kA * (tempS - tempA)) / step),
             fAU = ((kB * (tempB - tempA)) / step),
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
             rhoM = (moistureB+moistureB^2)*barkDensity,
             cpB = ((1105+4.85*tempB)*(1-moistureB)+moistureB*4185+1276*moistureB),
             kB = (2.104*barkDensity+5.544*rhoM+3.266*tempB-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fBD = ((kB * (tempA - tempB)) / step),
             fBU = ((kC * (tempC - tempB)) / step),
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
             rhoM = (moistureC+moistureC^2)*barkDensity,
             cpC = ((1105+4.85*tempC)*(1-moistureC)+moistureC*4185+1276*moistureC),
             kC = (2.104*barkDensity+5.544*rhoM+3.266*tempC-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fCD = ((kC * (tempB - tempC)) / step),
             fCU = ((kC * (tempD - tempC)) / step),
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
             rhoM = (moistureD+moistureD^2)*barkDensity,
             cpD = ((1105+4.85*tempD)*(1-moistureD)+moistureD*4185+1276*moistureD),
             kD = (2.104*barkDensity+5.544*rhoM+3.266*tempD-166.216)*10^-4,
             # Conduction from above and below, less latent heat of evaporation
             fDD = ((kD * (tempC - tempD)) / step),
             fDU = ((kD * (tempE - tempD)) / step),
             fourierD = fDD + fDU - max(0, min((fDD + fDU), drainD)),
             tempD = (fourierD / (mass * cpD) + tempD),
             # Change in proportion water this step
             moistureD = ifelse(tempD>99,ifelse(moistureD>0,max(0,moistureD-((fourierC/2256400)/mWaterD)),
                                                moistureD),moistureD),
             
             #phloem ________________________________________________________
             ### Water effects: evaporation and energy drain
             # Mass of water
             mWaterE = moistureE*massW,
             # Energy removed by current water quantity
             drainE = ifelse(tempE>99,
                             ifelse(moistureE>0,mWaterE*2256400,0),0),
             #Thermal values - (cp: J/g/deg, k: W/m/deg)
             kAir = 0.00028683*(tempE+273.15)^0.7919,
             cpE = 1080+408*moistureE+2.53*tempE+6.28*moistureE*tempE,
             kE = kWood(tempE, woodDensity, kAir),
             # Conduction from above and below, less latent heat of evaporation. Below unknown.
             fED = ((kE * (tempD - tempE)) / phloem),
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
    
    setTxtProgressBar(pbar,tS)
    ##  progress bar
    Sys.sleep(0.25)
    ####UpdateProgress
    if (is.function(updateProgress)) {
      text <- paste0("Number of remaining steps is ", TIME - tS)
      updateProgress(detail = text)
    }
    tS <- tS + 1
    Horiz <- Horiz - ROS
  }
  
  # Create table
  Ca <- Ca %>%
    select(tS, repId, tempS, tempA, tempB, tempC, tempD, tempE,
           moistureA, moistureB, moistureC, moistureD, moistureE)%>%
    mutate(scar = ifelse(tempE>=necT, 1, 0))
  
  return(Ca)
}