#' Operating_depth computes the maximum (default, limit = "bottom") and minimum (when limit = "top") operating depth for one or more gases.
#' O2: tibble, vector or scalar. Proportion of O2 in the mix in percent.
#' PO2: double. PO2 for which the limit is to be calculated (typically, 0.18 or 1.61)
#' limit: character. One of "bottom" or "top". Should be set to "bottom" (default) for Maximum OD, and to "top" for Minimum OD.
Operating_depth <- function(O2, PO2, limit = c("bottom", "top")) {
    ppO2 <- O2 / 100
    pressure <- PO2 / ppO2
    depth <- (pressure - 1) * 10
    depth[depth < 0] <- 0
    if (limit == "top") {
        return(ceiling(depth))
    } else {
        return(floor(depth))
    }
}
