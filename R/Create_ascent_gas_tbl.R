#' Create the ascent gas tibble
#'
#' @param Gas_list 3-colum tibble with proportions of O2, N2 and He in percent.
#' @param Max_depth double. Maximum depth of the dive
#'
#' @return A tibble containing one row for each segment of the ascent
#' @export
#'
#' @examples
Create_ascent_gas_tbl <- function(Gas_list, Max_depth){
    # First, we create a vector of all potential stops from Max_depth to the surface
    Pit_stops <- rev(seq(0, Max_depth, 3))
    # We repeat each element 4 times in order to create deco stops at each depth
    Pit_stops <- rep(Pit_stops, each = 4)
    # We then adjust the initial and final depths depending on wether Max_depth is included in the vector or not
    if (Pit_stops[1] != Max_depth) {
        Pit_stops <- c(Max_depth, Pit_stops)
    } else {
        Pit_stops <- Pit_stops[-(1:3)]
    }
    Pit_stops <- Pit_stops[-length(Pit_stops)]

    # Then, we switch to the tibble format
    Ascent_gas <- tibble::as.tibble(matrix(Pit_stops, ncol=2, byrow = TRUE))
    Ascent_gas <- Ascent_gas[-nrow(Ascent_gas),]
    colnames(Ascent_gas) <- c("Start_depth", "End_depth")

    # Finally, we pick the appropriate gas for each segment
    pick_gas <- function(Start_depth, Gas_list) {
        index <- max(which(Gas_list$Max_OD >= Start_depth))
        return(Gas_list[index,])
    }

    # And we concatenate the 2 tibbles to produce the final table
    tmp <- Ascent_gas$Start_depth %>% purrr::map_df(~pick_gas(.x, Gas_list = Gas_list))
    Ascent_gas <- dplyr::bind_cols(tmp, Ascent_gas)
    return(Ascent_gas)
}
