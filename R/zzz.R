.cacheEnv <- new.env(parent = emptyenv())

# This flag variable is a work-around for building package docs. When 
# devtools builds the package it calls roxygen2::roxygenize to build
# the docs. This function needs to load the package, but attempting to
# do so results in rscala complaining that it cannot find the jar files 
# required for the package because the package has not been built - Catch 22.
# I think we can fix this by updating to rscala 3.x. Until then, do the 
# following when you need to update the package docs:
#
# 1. Ensure that "Generate documentation with roxygen" is unchecked in
#    the RStudio tools->Project options->Build tools dialog.
#
# 2. Delete all files in the 'man' folder.
# 
# 3. Set the DUMMY_LOAD flag (below) to TRUE.
#
# 4. Call devtools::document()
# 
# 5. Set the DUMMY_LOAD flag (below) to FALSE.
#
# 6. Build the package from the RStudio 'Build' tab.

DUMMY_LOAD <- FALSE

.onLoad <- function(libname, pkgname) {
  if (!DUMMY_LOAD) {
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
}
