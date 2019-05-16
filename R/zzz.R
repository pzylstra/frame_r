.onLoad <- function(libname, pkgname) {
  packageStartupMessage("frame: Fire Research And Modelling Environment")
  
  ver <- .check_java_version()
  if (is.na(ver)) {
    warning("Java was not found on the local system path\n",
            "Java version 1.8 or higher is required")
    
  } else if (ver < 1.8) {
    warning("Java version ", ver, " found\n",
            "Version 1.8 or higher is required")
    
  } else {
    packageStartupMessage("Compatible Java version found (", ver, ")")
  }
}


.check_java_version <- function() {
  x <- system2("java", args = "-version", stdout = TRUE)

  found <- any( grepl("java version", x, ignore.case = TRUE) )
  
  if (!found) NA
  else as.numeric(stringr::str_extract(x[1], "\\d\\.\\d+"))
}
