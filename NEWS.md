# NEWS

## Version 0.2.2
 * Function readLegacyParamFile has been deprecated. New code should now use
 `ffm_read_legacy_params`.
 
 * Added function `ffm_write_params` to write parameters to file as a CSV-format
 table. This is just a simple wrapper that calls `write.csv` with appropriate
 default arguments and optionally checks parameters before writing.

## Version 0.2.1
 * Functions `ffm_run` and `ffm_run_command` now attempt to deal with embedded
 spaces in paths and file names.

## Version 0.2.0

 * Removed dependency on the `rScala` package. Function `ffm_run` is now a
 simple function (rather than an S3 generic). It passes simulation parameters
 and details of the output database to the external Scala application via a
 temporary file and system call. When the function is run with an existing
 output file specified, users can select between adding simulation results to
 the existing database as new replicates, or discarding previous results.
 
 * Updated Java JAR files for version 0.2 of the Scala application and its
 dependencies.
   
 * Added function `ffm_check_java` to check that Java version 1.8 (aka Java 8)
 or higher is installed and accessible from the command line. This check is run
 the first time that `ffm_run` is called during a session.
 
 * Added function `ffm_check_params` to test whether all required parameters
 have been provided for the simulation.
 
 * The data frame `DefaultSpeciesParams` has been removed from the package to
 avoid clashes with default data provided by the user. A data frame of default
 parameter values can now be passed to function `ffm_run` which will call
 `ffm_complete_params`. Alternatively, the user can call `ffm_complete_params`
 directly before running the simulation.
 
 * Function `ffm_complete_params` now allows for species to appear in more than
 one stratum and, where this is the case, each species occurrence is considered
 separately. Default values for a given parameter may now be specified for just
 a subset of the species included in the defaults data frame, with others having
 `NA` values.
 
 * Functions to work with the copy of `DefaultSpeciesParams` previously included
 in the package, such as `ffm_find_species` and `ffm_species_known` have been
 removed.
   
