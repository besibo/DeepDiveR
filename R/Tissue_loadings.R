#' Computes loadings in N2 and He for the 16 tissues defined by Bühlmann (table ZH-L16 C)
#'
#' @param Segment a one row tibble containing (at least) the following variables: N2, He, Start_depth, End_depth, Duration, N2_Load_Start and He_Load_Start
#' @param Penalty either 1 (most permissive decompression model), 2 or 3 (default, most conservative decompression model)
Tissue_loadings <- function(Segment, Penalty = 3) {

    Period <- ZHL16_C %>% select(Molecule, Periode)
    Ambiant_pressure <- Segment$Start_depth + 10

    Quot	<- c(0.627,0.567,0.493)
    P_H2O	<- Quot[Penalty]
    Pi_N2	<- (Ambiant_pressure - P_H2O) * Segment$N2 / 100
    Pi_He	<- (Ambiant_pressure - P_H2O) * Segment$He / 100

    if (Segment$Start_depth != Segment$End_depth) {
        Speed <- round((Segment$End_depth-Segment$Start_depth)/Segment$Duration)
        R_N2 <- Speed * Segment$N2 / 100
        R_He <- Speed * Segment$He / 100
    } else {
        R_N2 <- 0
        R_He <- 0
    }

    K <- log(2) / Period$Periode
    K_N2 <- K[1:16]
    K_He <- K[17:32]

    Po_N2 <- Segment$N2_Load_Start[[1]]
    Po_He <- Segment$He_Load_Start[[1]]

    t <- Segment$Duration

    P_N2 <- Pi_N2 + R_N2 * (t - 1/K_N2) - (Pi_N2 - Po_N2 - R_N2/K_N2) * exp(-K_N2 * t)
    P_He <- Pi_He + R_He * (t - 1/K_He) - (Pi_He - Po_He - R_He/K_He) * exp(-K_He * t)
    list(P_N2, P_He)
}
