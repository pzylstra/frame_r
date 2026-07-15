#Original
#.onLoad <- function(libname, pkgname) {
#  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
#  assign("._ffm_settings", list(), pos = 1)
#}

.onLoad <- function(libname, pkgname) {
  # initialise settings in the namespace
  pkg_env <- getNamespace(pkgname)
  assign(".ffm_settings", list(), envir = pkg_env)
}

.onAttach <- function(libname, pkgname) {
  # only print the startup message
  packageStartupMessage("FRaME: Fire Research And Modelling Environment")
}



# Original
#.onUnload <- function(libpath) {
#  if (exists(._ffm_settings, where = 1)) rm(._ffm_settings, pos = 1)
#}

.onUnload <- function(libpath) {
  # Check if ._ffm_settings exists in the global environment
  if (exists("._ffm_settings", envir = .GlobalEnv)) {
    # Remove ._ffm_settings from the global environment
    rm("._ffm_settings", envir = .GlobalEnv)
  }
}



# Declare global variables and functions
utils::globalVariables(c("Age", "Alpha", "AlphaP", "Altitude", "Angle", "C", "C.C_rat", "C.C_rat.x", "C.C_rat.y", "Cloud", "Co", "CP", "Cover",
                         "Day", "Declination", "Density", "E", "El", "Extinction", "FF16_expand_state", "flameTempP", "Flora",
                         "G.C_rat", "G.C_rat.x", "G.C_rat.y", "Genus", "Hc", "He", "Height", "Hour", "Hp", "Ht", "InsolationA",
                         "InsolationB", "InsolationC", "InsolationD", "Intercept", "InterceptP", "LegacyParamLookup", "MSLP",
                         "Mr", "Ms", "ParamInfo", "Plume_velocity", "Pm", "Point", "Pt", "Q", "Qi", "QiA", "QiB", "QiC", "QiD", "QiE", "R",
                         "RH", "RHA", "RainAdj", "Reach", "Richness", "Separation", "Shade", "ShadeA", "ShadeB", "ShadeC",
                         "ShadeD", "Site", "Slice", "Species", "Spotting", "Stratum", "Structure", "TempA", "TerrainA",
                         "TerrainB", "TerrainC", "TerrainD", "Test", "WPLCP", "Wc", "We", "Weight", "Width", "Wp", "Wt", "Ww",
                         "ZenithA", "ZenithB", "ZenithC", "ZenithD", "a", "across", "angle", "angleSurface",
                         "angle_degrees", "aov", "area_leaf", "att", "b", "b1", "b2", "b3", "b4", "bMoisture", "bark",
                         "barkDensity", "base", "baseM", "branchA", "branchV", "cAM", "cBase", "cPM", "cTop", "c_vol",
                         "cardinal", "case_when", "clim", "clumpD", "cluster", "comBark", "comp", "comp.x", "comp.y",
                         "complete.cases", "composition", "compression", "cor", "count", "cover", "cpA", "cpAir", "cpB",
                         "cpBark", "cpC", "cpD", "cpE", "cpWoodB", "cpWoodC", "cpWoodD", "cpWoodE", "cr", "d", "dAM", "dPM", "dbeta",
                         "deadFuelMoistureProp", "default.species.params", "delim", "density", "desc", "di",
                         "distance", "drain", "drainA", "drainB", "drainC", "drainD", "drainE", "drop_na", "eBase", "eTop",
                         "epsilon", "everything", "extinct", "fAD", "fAU", "fBD", "fBU", "fCD", "fCU", "fDD", "fDU", "fED", "fEU",
                         "ffm_import_legacy_params", "fh", "filter", "fires", "flameAngle", "flameHeight",
                         "flameLength", "flameTemp", "fourier", "fourierA", "fourierB", "fourierC", "fourierD",
                         "fourierE", "fourierO", "fourierOA", "fourierOB", "fourierOC", "fourierOD", "fourierOE",
                         "frontalArea", "fuelLoad", "full_join", "furCp", "furDensity", "girdleH", "h", "hc", "hcR", "he",
                         "heM", "heR", "height", "heightM", "heightPlant", "heightSurface", "hor", "hourAngle",
                         "hourStep", "hp", "ht", "htM", "htP", "htR", "ignitionTemp", "ignitionTemp.x", "ignitionTemp.y",
                         "is", "is.error", "kAir", "kBark", "kFur", "kWind", "kWoodB", "kWoodC", "kWoodD", "kWoodE", "kmeans",
                         "l", "lAngleAccounting", "lRiver", "lSep", "lat", "leafForm", "leafForm.x", "leafForm.y",
                         "leafLength", "leafLength.x", "leafLength.y", "leafSeparation", "leafSeparation.x",
                         "leafSeparation.y", "leafThickness", "leafThickness.x", "leafThickness.y",
                         "leafWidth", "leafWidth.x", "leafWidth.y", "leaf_area", "leavesClump", "lengthPlant",
                         "lengthSurface", "level", "lightning", "line", "lm", "lma", "m", "mBase", "mRiver", "mTop", "mWater",
                         "mWaterA", "mWaterB", "mWaterC", "mWaterD", "mWaterE", "maxW", "median", "minR", "moisture",
                         "moisture.x", "moisture.y", "mortality", "nClumps", "nLeaves", "n_distinct", "na.omit",
                         "name", "necT", "necrosis", "new.value", "nls", "nls.control", "nsBase", "nsTop", "optimize",
                         "pAlpha", "pAlphaPost", "pAlphap", "pAlphas", "pN", "param", "pelMass", "phloem",
                         "pivot_longer", "plantVarFrame", "point", "postS", "predict", "presAtm", "probFire_Frame", "probStrike",
                         "propDead", "propDead.x", "propDead.y", "propSilicaFreeAsh", "pt", "qR", "qc", "qrO",
                         "quantile", "rangeDir", "rbeta", "record", "repAngle", "repHeight", "repId", "repLength",
                         "resBark", "rhAM", "rhPM", "rhoM", "rholitter", "right_join", "rnorm", "ros_kph", "row_number",
                         "rtnorm", "runIndex", "runif", "sN", "saturation", "saturationA", "saturationB",
                         "saturationC", "saturationD", "saturationE", "sc1", "sc2", "sc3", "sc4", "sd", "section", "seedDa",
                         "seedDb", "seedDc", "seedDd", "seedDe", "segIndex", "select", "sensitive", "sep", "separation",
                         "setTxtProgressBar", "severity", "slope", "slopeLength", "slopeM", "slopeSD",
                         "slope_degrees", "smoulder", "spBase", "spName", "spTop", "specHumAM", "specHumPM", "specPoint", "spread",
                         "st", "stemOrder", "stemOrder.x", "stemOrder.y", "stratify_community", "summarise",
                         "summarise_all", "summarise_if", "summarize_all", "summarize_if", "surfDecl",
                         "surfPost", "tAM", "tMax", "tMin", "tPM", "tail", "targSp", "taxon_name", "tempAir", "tempS",
                         "tempSoil", "temp_pointP", "temp_pointS", "temp_point_post", "temperature", "testN",
                         "top", "topM", "trail", "txtProgressBar", "vAir500", "value", "var", "viscosity", "w", "wAM", "wM",
                         "wPM", "wR", "weighted.mean", "weightedW", "wetBulb", "where", "wid", "windSpeed", "wind_kph",
                         "woodDensity", "word", "write.csv", "x", "x0", "x1", "y", "y0", "y1", "z", "zeta", "phi"))
