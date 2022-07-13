#' Valid parameters and their default units.
#' 
#' This data frame contains all valid parameter labels with their corresponding 
#' section (site, stratum or species) and default units.
#' A value of \code{NA} for units indicates that a parameter is dimensionless 
#' (e.g. a proportional quantity or a text value).
#'
#' @format A data frame with columns: section (values site,
#' stratum, species); param (standard parameter label);
#' units (unit abbreviation or \code{NA} if dimensionless).
#'
"ParamInfo"


#' @format a data frame with 13 columns, and a row for each species per stratum, per site
#'\describe{
#' \item{record}{Consecutive integer to reference the row}
#' \item{site}{A string to name the row}
#' \item{species}{The name of the species, e.g. Eucalyptus rossii}
#' \item{moisture}{Moisture of live leaves (proportion oven dry weight) e.g 0.8}
#' \item{stratum}{Numbered from 1 (near surface) to 4 (canopy)}     
#' \item{comp}{Integer giving the count of this species}
#' \item{base}{Height (m) of the base of the plant crown}
#' \item{he}{Height (m) of the lower edge of the plant crown}
#' \item{ht}{Height (m) of the upper edge of the plant crown}
#' \item{top}{Height (m) of the top of the plant crown}
#' \item{w}{Width (m) of the plant crown}
#' \item{openness}{Ratio of gaps between leaf clumps to clump diameter}
#' \item{clump}{Ratio of individual clump volume to crown volume}
#' 
#'}
#'
"f_flora"



#'@format a dataframe with 8 columns, and a row for each site
#' \describe{
#' \item{record}{Consecutive integer to reference the row}
#' \item{site}{A string to name the row}
#' \item{slope}{Slope (degrees)}
#' \item{wind}{Wind velocity (km/h)}
#' \item{temp}{Ambient air temperature (deg.C)}
#' \item{dfmc}{Dead fuel moisture content (proportion oven dry weight)}
#' \item{litter}{Weight (t/ha) of fine surface litter}
#' \item{fLine}{Length of the active fire line (m)}
#' \item{diameter}{mean diameter of litter particles (m)}
#'}   
#'     
"f_site"




#'@format a dataframe with 11 columns and 1 row for each site
#'\describe{  
#' \item{record}{Consecutive integer to reference the row} 
#' \item{site}{A string to name the row}
#' \item{NS}{Separation between plants in the near surface stratum (m from centre to centre)}   
#' \item{El}{Separation between plants in the elevated stratum (m from centre to centre)}      
#' \item{Mid}{Separation between plants in the midstorey stratum (m from centre to centre)}       
#' \item{Can}{Separation between plants in the canopy stratum (m from centre to centre)}        
#' \item{ns_e}{Near surface plants grow directly beneath elevated plants. 
#' Acceptable values are t, f, or blank, 
#' where the outcome will be decided by the relative stratum heights.}         
#' \item{ns_m}{Near surface plants grow directly beneath midstorey plants. 
#' Acceptable values are t, f, or blank, where the outcome will be decided by the 
#' relative stratum heights.}          
#' \item{e_m}{Elevated plants grow directly beneath midstorey plants.
#'  Acceptable values are t, f, or blank, where the outcome will be 
#'  decided by the relative stratum heights.}          
#' \item{e_c}{Elevated plants grow directly beneath canopy plants. 
#' Acceptable values are t, f, or blank, where the outcome will be decided by 
#' the relative stratum heights.}        
#' \item{m_c}{Midstorey plants grow directly beneath canopy plants. 
#' Acceptable values are t, f, or blank, where the outcome will be 
#' decided by the relative stratum heights.}            
#' \item{nsR}{Optional field used only by the function frame::specPoint().
#'  Gives the maximum recorded species richness for the near-surface stratum.}         
#' \item{eR}{Optional field used only by the function frame::specPoint(). 
#' Gives the maximum recorded species richness for the elevated stratum.}          
#' \item{mR}{Optional field used only by the function frame::specPoint(). 
#' Gives the maximum recorded species richness for the midstorey stratum.}               
#' \item{cR}{Optional field used only by the function frame::specPoint().
#'  Gives the maximum recorded species richness for the canopy stratum.}        
#'}                    
#'                     
#'                       
"f_structure"







#'@format a data frame with 9 columns, and a row for each species
#'\describe{
#' \item{name}{The name of the species}
#' \tem{propDead}{Proportion of the foliage that is dead}
#' \item{leafForm}{Allowable values are flat or round (e.g. terete)}
#' \item{leafThickness}{Thickness of the leaf (m)}
#' \item{leafWidth}{Width at the widest axis of the leaf (m)}
#' \item{leafLength}{Length of the leaf (m)}
#' \item{leafSeparation}{Separation between leaves along a stem}
#' \item{stemOrder}{Stem ramification}
#' \item{ignitionTemp}{Piloted ignition temperature of the leaf (deg.C)}
#' \item{Field}{Description}
#'}
"f_traits"



#'
#'@format a data frame with 4 columns and a row for each time step
#'\describe{
#' \item{tm}{Consecutive integer to reference the time period}
#' \item{T}{Ambient air temperature (deg. C)}
#' \item{W}{Wind velocity (km/h)}
#' \item{DFMC}{Dead fuel moisture content (proportion oven dry weight)}
#'
#'}
#'
"f_weather"



#'
#'@format a data frame with 9 columns and a row for each point of measurement
#'\describe{
#' \item{Site}{A numeric record of the site}
#' \item{SiteName}{A string name for the site}
#' \item{Point}{The point number in the transect, often the metre mark}
#' \item{Species}{The name of the species, consistent with other tables}
#' \item{base}{Crown base height hc (m)}
#' \item{top}{Crown top height hp (m)}
#' \item{he}{Crown he height (m)}
#' \item{ht}{Crown ht height (m)}
#' \item{width}{Crown width w (m)}
#'
#'}
#'
#'@source {Raw survey data collected using a 3d transect at Warrungup Spring reserve, April 2019}
#'
"tuart_survey"



#'
#'@format a data frame with 12 columns and a row for each species
#'\describe{
#' \item{name}{The name of the species, consistent with other tables}
#' \item{propDead}{Value between 0 & 1 giving the proportion of crown foliage that is dead}
#' \item{leafForm}{Either Flat or Round}
#' \item{leafThickness}{Thickness of the leaf (m)}
#' \item{leafWidth}{Width of the leaf at the widest point (m)}
#' \item{leafLength}{Length of the leaf (m)}
#' \item{leafSeparation}{Mean distance between leaf petioles on a 1st-order stem (m)}
#' \item{stemOrder}{Stem ramification or order}
#' \item{ignitionTemp}{Tempereature of the endotherm or piloted ignition temp (C)}
#' \item{moisture}{Ratio of moisture weight to oven-dry weight}
#' \item{G.C_rat}{Ratio of gaps to clumps, or crown space to foliage clumps}
#' \item{C.C_rat}{Ratio of individual clump area to crown area}
#'
#'}
#'
#'@source {Measurements from plants surveyed at Warrungup Spring reserve, April 2019}
#'
"tuart_traits"



#'
#'@format a data frame with 5 columns and a row for each time step
#'\describe{
#' \item{tm}{Consecutive integer to reference the time period}
#' \item{T}{Ambient air temperature (deg. C)}
#' \item{W}{Wind velocity (km/h)}
#' \item{DFMC}{Dead fuel moisture content (proportion oven dry weight)}
#' \item{time}{Hour of the day)}
#'
#'}
#'
#'@source {BOM AWS data for 16th May 2018, Mandurah weather station 009977. DFMC calculated using Matthews simple model}
#'
"tuart_weather"














