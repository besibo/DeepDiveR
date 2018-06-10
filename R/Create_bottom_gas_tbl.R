#' This function creates a tibble containing all relevant informations about the bottom gas
#' Desc_gas: a tibble as produced by Create_desc_gas_tbl()
Create_bottom_gas_tbl <- function(Desc_gas){
    # The bottom gas is the last descent gas
    Bottom_gas <- Desc_gas[nrow(Desc_gas), ]
    # It is used to stay at Max_depth
    Bottom_gas$Start_depth <- Bottom_gas$End_depth
    return(Bottom_gas)
}
