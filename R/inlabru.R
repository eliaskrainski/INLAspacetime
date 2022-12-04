#' @title Mapper object for automatic inlabru interface
#' @aliases bru_get_mapper
#'
#' @description
#' Return an `inlabru` `bru_mapper` object that can be used for computing
#' model matrices for the space-time model components.  The `bru_get_mapper()`
#' function is called by the `inlabru` methods to automatically obtain the
#' needed mapper object (from `inlabru` `2.7.0.9001`; before that, use
#' `mapper = bru_get_mapper(model)` explicitly).
#' @param model The model object (of class `stModel_cgeneric`, from
#' `stModel.define`)
#' @param \dots Unused.
#' @return A `bru_mapper` object of class `bru_mapper_multi` with
#' sub-mappers `space` and `time` based on the model `smesh` and `tmesh`
#' objects.
#' @seealso [inlabru::bru_get_mapper()]
#'
#' @rawNamespace S3method(inlabru::bru_get_mapper, stModel_cgeneric)
bru_get_mapper.stModel_cgeneric <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper_multi(list(
    space = inlabru::bru_mapper(model[["smesh"]]),
    time = inlabru::bru_mapper(model[["tmesh"]], indexed = TRUE)
  ))
}
