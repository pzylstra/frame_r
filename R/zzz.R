.onLoad <- function(libname, pkgname) {
  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
  info <- rscala::scalaInfo()
  if (is.null(info))
    stop("Cannot find Scala on this system\n",
         "Please install manually or run rscala::scalaInstall()\n")
  
  rscalaPackage(pkgname)
  
  packageStartupMessage("Establishing connection to Scala...", appendLF = FALSE)
  rscalaLoad()
  packageStartupMessage("  ready")
}
