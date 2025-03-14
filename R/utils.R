#' check package version and load
#' @param pkg character with the name of the package
#' @param minimum_version character with the minimum required version
#' @param quietly logical indicating if messages shall be printed
#' @note
#' this is the same code as in the inlabru package
check_package_version_and_load <-
  function(pkg, minimum_version, quietly = FALSE) {
    version <- tryCatch(utils::packageVersion(pkg),
                        error = function(e) NA_character_
    )
    if (is.na(version)) {
      if (!quietly) {
        message(paste0("Package '", pkg, "' is not installed."))
      }
      return(NA_character_)
    }
    if (version < minimum_version) {
      if (!quietly) {
        message(paste0(
          "Installed '", pkg, "' version is ", version, " but ",
          "version >= ", minimum_version, " is required."
        ))
      }
      return(NA_character_)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!quietly) {
        message("Package '", pkg, "' not loaded safely.")
      }
      return(NA_character_)
    }
    return(version)
  }
