.onLoad <- function(libname, pkgname) {
  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
  info <- rscala::scalaInfo(verbose = FALSE)
  if (is.null(info))
    stop("Cannot find Scala on this system\n",
         "Please install manually or run rscala::scalaInstall()\n")
  
  packageStartupMessage("Establishing connection to Scala...", appendLF = FALSE)
  
  rscala::scalaPackage(pkgname, assign.name = "frame_scala__")
  
  packageStartupMessage("  ready")
}

.onUnload <- function(libpath) {
  scalaPackageUnload()
}
