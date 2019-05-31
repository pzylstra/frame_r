.onLoad <- function(libname, pkgname) {
  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
  assign("._ffm_settings", list(), pos = 1)
}

.onUnload <- function(libpath) {
  if (exists(._ffm_settings, where = 1)) rm(._ffm_settings, pos = 1)
}
