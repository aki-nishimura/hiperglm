#' @export
print.hglm <- function(hglm_out) {
  cat("`hiper_glm` output\n")
}

#' @export
coef.hglm <- function(hglm_out) {
  return(hglm_out$coef)
}

#' @export
vcov.hglm <- function(hglm_out) {
  return(hglm_out$cov)
}
