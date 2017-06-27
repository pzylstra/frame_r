.cacheEnv <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
  info <- rscala::scalaInfo()
  if (is.null(info))
    stop("Cannot find Scala on this system\n",
         "Please install manually or run rscala::scalaInstall()\n")
  
  packageStartupMessage("Establishing connection to Scala...", appendLF = FALSE)
  .rscalaPackage(pkgname)

  # Give the default interpreter a better name and put
  # it in an environment so that we can update the interpreter
  # object later if required.
  assign("interp", s, .cacheEnv)
  
  packageStartupMessage("  ready")
}
