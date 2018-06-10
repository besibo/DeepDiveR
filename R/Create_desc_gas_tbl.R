#' This function produces a tibble of gases to breathe during the descent
#' Gas_list: a tibble of 5 columns: O2, N2, He, Min_OD and Max_OD
#' Returns a tibble with two additional columns (Start_depth and End_depth) and only relevant gases sorted in the correct order for the descent
Create_desc_gas_tbl <- function(Gas_list, Max_depth) {
    # Since gases are ordered by increasing O2 concentration, the first gaz to breath is the first gas in the tibble that has a Min_OD of 0
    First_gas_index <- which.min(Gas_list$Min_OD)

    # Then, eventual subsequent gases are the rows of the tibble in reverse order
    Desc_gas_index <- First_gas_index:1

    # Here, we create a new tibble that will contain all relevant information about gases we'll breath during the descent
    Desc_gas <- Gas_list[Desc_gas_index, ]

    # We then compute the start and end depth for each gas and add these to the Desc_gaz tibble
    Start_depth <- c(0, Desc_gas$Max_OD)
    Start_depth[length(Start_depth)] <- Max_depth
    End_depth <- Start_depth[-1]
    Start_depth <- Start_depth[-length(Start_depth)]
    Desc_gas <- Desc_gas %>%
        mutate(Start_depth, End_depth)

    return(Desc_gas)
}
