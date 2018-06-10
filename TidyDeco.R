library(tidyverse) ; library(DeepDiveR)

# Dive parameters
Max_depth <- 98
Bottom_time <- 17
Speed_desc <- 30     ; Speed_asc <- 10
Low_gradient <- 0.20 ; High_gradient <- 0.75
PO2_min <- 0.18 ; PO2_max <- 1.61
Gas_list <- rbind(c(35, 0), c(65,0), c(100,0), c(10,55), c(18, 34))
Penalty <- 3
Steps <- 0.5
Manual_gaz_change_down <- FALSE ; Gaz_change_down <- c(20,40)
Manual_gaz_change_up <- FALSE   ; Gaz_change_up <- c(60,30,6)

# ------------------------------------------------------------------------------------

# 1. Create a full list of gases
colnames(Gas_list) <- c("O2", "He")
Gas_list <- as.tibble(Gas_list) %>%
    mutate(N2 = 100 - O2 - He) %>%
    select(O2, N2, He) %>%
    arrange(O2)

# 2. Compute Minimum and Maximum Operating Depths for each gaz
Gas_list <- Gas_list %>%
    mutate(Min_OD = Operating_depth(O2, PO2_min, "top"),
           Max_OD = Operating_depth(O2, PO2_max, "bottom"))

# 3. Create a tibble containing gases to used during the descent
Desc_gas <- Create_desc_gas_tbl(Gas_list = Gas_list, Max_depth = Max_depth) %>%
    mutate(Type = "Descent")

# 4. Create a tibble containing the gas to use as the bottom of the dive
Bottom_gas <- Create_bottom_gas_tbl(Desc_gas)  %>%
    mutate(Type = "Bottom")

# 5. Create a tibble containing the gas to use during the ascent
Ascent_gas <- Create_ascent_gas_tbl(Gas_list = Gas_list, Max_depth = Max_depth) %>%
    mutate(Type = "Ascent")

# 6. Concatenate all 3 tables, compute PO2s and EADs, and add empty list columns for loadings
All_segments <- bind_rows(Desc_gas, Bottom_gas, Ascent_gas) %>%
    mutate(EAD_min = EAD(N2, Depth = Start_depth),
           EAD_max = EAD(N2, Depth = End_depth),
           PO2_start = (O2 / 100) * (Start_depth/10 + 1),
           PO2_end   = (O2 / 100) * (End_depth/10 + 1)) %>%
    mutate(Time_start = NA,
           Time_end = NA) %>%
    mutate(N2_Load_Start = list(rep(0,16)),
           N2_Load_End = list(rep(0,16)),
           He_Load_Start = list(rep(0,16)),
           He_Load_End = list(rep(0,16)))

# 7. Compute times for the descent and bottom segments
Desc_segments <- All_segments %>%
    filter(Type %in% c("Descent", "Bottom"))

## Computations for the first segment
Desc_segments$Time_start[1] <- 0

Depth_drop <- Desc_segments %>%
    transmute(Drop = End_depth - Start_depth) %>%
    mutate(Duration = Drop / Speed_desc)

Desc_segments$Time_end[1] <- Desc_segments$Time_start[1] + Depth_drop$Duration[1]

## Computations for all other segments
for (i in 2:nrow(Desc_segments)) {
    Desc_segments$Time_start[i] <- Desc_segments$Time_end[i-1]
    Desc_segments$Time_end[i] <- Desc_segments$Time_start[i] + Depth_drop$Duration[i]
}

## Adjusting the last value (bottom time)
Desc_segments$Time_end[nrow(Desc_segments)] <- Bottom_time

## Add a "Duration" variable
Desc_segments <- Desc_segments %>%
    mutate(Duration = Time_end - Time_start)


# 8. Compute tissue loadings for both N2 and He
## Manually compute the first segment
Quot	<- c(0.627,0.567,0.493)
P_H2O	<- Quot[Penalty]
Desc_segments$N2_Load_Start[[1]] <- rep((10 - P_H2O)*0.795, 16)
tmp <- Tissue_loadings(Desc_segments[1,], Penalty = Penalty)
Desc_segments$N2_Load_End[[1]] <- tmp[[1]]
Desc_segments$He_Load_End[[1]] <- tmp[[2]]

## Compute loadings for all other segments
for (i in 2:nrow(Desc_segments)) {

    Desc_segments$N2_Load_Start[[i]] <- Desc_segments$N2_Load_End[[i-1]]
    Desc_segments$He_Load_Start[[i]] <- Desc_segments$He_Load_End[[i-1]]
    tmp <- Tissue_loadings(Desc_segments[i,], Penalty = Penalty)
    Desc_segments$N2_Load_End[[i]] <- tmp[[1]]
    Desc_segments$He_Load_End[[i]] <- tmp[[2]]
}




# ------------------------------------------------------------------------------------

# !!! Check if there is at least one gas breathable from the surface !!!
# # 5. Determine which gas can be used from the surface
# Gas_list <- Gas_list %>%
#     mutate(From_surface = (Min_OD == 0))

# ------------------------------------------------------------------------------------


both <- function(half.t,prof.init,prof.fin,vitesse,gaz,molec=c("azote","helium"),dur.qr,t=NULL,Pprec=NULL) {
    vitesse <- if (prof.fin>prof.init) vitesse else -1*vitesse
    t		<- if (prof.init != prof.fin) (prof.fin-prof.init) / vitesse else t
    pp.He <- pp.comp(half.t[,1],prof.init=prof.init,prof.fin=prof.fin,vitesse=vitesse,gaz=gaz[2]/100,molec=molec[2],dur.qr=dur.qr,t=t,Pprec=Pprec[,1])
    pp.N2 <- pp.comp(half.t[,2],prof.init=prof.init,prof.fin=prof.fin,vitesse=vitesse,gaz=gaz[3]/100,molec=molec[1],dur.qr=dur.qr,t=t,Pprec=Pprec[,2])
    out <- cbind(pp.He,pp.N2)
    row.names(out) <- c("1b",2:16)
    out
}

# Pi : pression de gaz inerte inspirée (alvéolaire) initiale
# Po : pression de gaz inerte initiale dans le compartiment
# R : rythme du changement de pression du gaz inspiré avec le changement de pression ambiante
# t : durée
# k : constante de temps = ln 2 / période du compartiment


M.val <- function(a=a,b=b,ppart,prof) {
    coef.a 	<- (a[,2]*ppart[,2] + a[,1]*ppart[,1]) / apply(ppart,1,sum)
    coef.b 	<- (b[,2]*ppart[,2] + b[,1]*ppart[,1]) / apply(ppart,1,sum)
    M.val 	<- coef.a+(prof+10)/coef.b
    M.val
}

GF.paliers <- function(a,b,ppart,GF,prof.max) {

    coef.a <- (a[,2]*ppart[,2] + a[,1]*ppart[,1]) / apply(ppart,1,sum)
    coef.b <- (b[,2]*ppart[,2] + b[,1]*ppart[,1]) / apply(ppart,1,sum)

    charge.inerte <- apply(ppart,1,sum)

    pression.toleree <- (charge.inerte - coef.a * GF) / (GF/coef.b - GF + 1)
    if (any(pression.toleree < 0)) pression.toleree[which(pression.toleree<0)] <- 0

    stops <- seq(0,prof.max,3)
    out <- list(min(stops[stops>max(pression.toleree-10)]),pression.toleree)
    out
}


# ============================================================================
# = Fonction calculant la profondeur d'entrée dans la zone de décompression =
# ============================================================================
# La fonction procéde par dychotomie : on se place é la profondeur moyenne et on compare la pression é la tension du tissus directeur.
# Puis on recommence en redivisant l'intervalle en 2 jusqu'é ce que la différence entre pression ambiante et tension soit trés faible.
# é chaque fois, on calcule la charge aditionnelle des tissus entre leur charge au départ du fond, et leur charge é l'arrivée au niveau du palier.
# La profondeur d'entrée dans la zone de déco est alors arrondie au métre suppérieur (éa va dans le sens de la sécurité !).
deco.zone <- function(ppart.fond,prof.max,vit.rem,melange,periode) {

    pression.max <- prof.max + 10
    pression.min <- 10

    for (i in 1:100) {
        pression.milieu <- (pression.max + pression.min) / 2

        ppart.remontee <- both(half.t=periode,prof.init=prof.max,prof.fin=(pression.milieu-10),vitesse=vit.rem,gaz=melange,molec=c("azote","helium"),dur.qr=1,t=NULL,Pprec=ppart.fond)
        tension.tissus.max <- max(apply(ppart.remontee,1,sum))

        diff <- pression.milieu - tension.tissus.max

        if (diff <= 0.001) start.deco.zone <- pression.milieu - 10
        if (diff > 0) pression.max <- pression.milieu else pression.min <- pression.milieu

    }
    first.possible.stop <- start.deco.zone - start.deco.zone%%3
    out <- c(round(start.deco.zone,1),first.possible.stop)
    names(out) <- c("Debut zone deco","Premier palier possible")
    out
}

# ===========================
# = Pourcentage du gradient =
# ===========================
# Connaissant les M-values, la profondeur et les pressions partielles en gaz inerte dans chaque tissus, cette fonction calcul le pourcentage du gradient atteint par la saturation.
# Autrement dit, cette fonction calcul é quel point les tensions sont proches des M-values en pourcentage de gradient, c'est é dire par rapport é la ligne de pression ambiante
pc.grad <- function(ppart,prof,Mval) {
    Pamb <- prof+10
    tension <- apply(ppart,1,sum)
    out <- (tension - Pamb) * 100 / (Mval - Pamb)
}


# ==================================================
# = C'est LA fonction qui calcule la décompression =
# ==================================================
BIGONE <- function (prof.max=90,temps=20,vit.desc=30,vit.rem=10,gradient=c(0.20,0.75),gaz=rbind(c(15,53),c(40,0)),durete=3,PO2.min=0.18,PO2.max=1.61,pas=0.5,manu.down=F,gaz.change.down=c(20,40),manu.up=F,gaz.change.up=c(60,30,6)){

    # ================================================================================================
    # = ===================================== CREATION DES VARIABLES UTILES ======================== =
    # ================================================================================================

    # Data frame des M-Values pour l'azote et l'hélium pour les 16 compartiments (on utilise le 1b et non le 1)
    ZHL16.C <- data.frame(
        comp = 1:16,
        N2.periode = c(5,8,12.5,18.5,27,38.3,54.3,77,109,146,187,239,305,390,498,635),
        a.N2 = c(11.696,10.000,8.618,7.562,6.667,5.600,4.947,4.500,4.187,3.798,3.497,3.223,2.850,2.737,2.523,2.327),
        b.N2 = c(0.5577867,0.6513809,0.7221781,0.7824726,0.8126117,0.8433837,0.8692629,0.8910274,0.9091736,0.9221689,0.9318796,0.9402915,0.9476876,0.9543806,0.9602458,0.9653441),
        He.periode = c(1.88,3.02,4.72,6.99,10.21,14.48,20.53,29.11,41.20,55.19,70.69,90.34,115.29,147.42,188.24,240.03),
        a.He = c(16.189,13.830,11.919,10.458,9.22,8.205,7.305,6.502,5.95,5.545,5.333,5.189,5.181,5.176,5.172,5.119),
        b.He = c(0.4770082,0.5747126,0.6526989,0.7222824,0.7582076,0.7956715,0.8278831,0.8552857,0.8757334,0.8903134,0.8996851,0.9072764,0.9121591,0.9170946,0.9216590,0.9266982)
    )

    # Recupération des données de periode et des coef a et b pour les 2 gaz
    periode <- cbind(ZHL16.C$He.periode,ZHL16.C$N2.periode)
    a <- cbind(ZHL16.C$a.He,ZHL16.C$a.N2)
    b <- cbind(ZHL16.C$b.He,ZHL16.C$b.N2)

    # Initialisation des paramétres de gradient et de profondeur
    GF.Lo <- gradient[1]
    GF.Hi <- gradient[2]

    # Création de la matrice des mélanges gazeux
    if (is.vector(gaz)) gaz <- t(as.matrix(gaz))
    mel		 	<- prof.melange(prof.max=prof.max,melange=gaz,heliair=F,PO2.min=PO2.min,PO2.max=PO2.max)

    if (manu.down)	{
        mel.desc 	<- melange.matrix.desc.manuel(mel,p.max=prof.max,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.down)
    } else {	mel.desc 	<- melange.matrix.desc(mel,p.max=prof.max,PO2.min=PO2.min,PO2.max=PO2.max) }

    nb.mel.desc	<- nrow(mel.desc)
    mel.fond	<- mel.desc[nb.mel.desc,]
    if (is.vector(mel.desc)) mel.desc <- t(as.matrix(mel.desc))
    if (is.vector(mel.fond)) mel.fond <- t(as.matrix(mel.fond))

    # Création et initialisation des objets où seront stockées les valeurs de chaque segment.
    ppart.tot 	<- NULL
    Mval.tot  	<- NULL
    pc.grad.tot	<- NULL
    tension.tot	<- NULL
    temps.tot	<- NULL
    nb.segment <- 0


    # ==========================================================================================
    # = =================================== DEBUT DE LA PLONGEE ============================== =
    # ==========================================================================================
    # Saturation é la descente
    ppart.desc <- NULL
    for (i in 1:(nb.mel.desc-1)) {
        prof.init 		<- mel.desc[i,4]
        prof.fin 		<- mel.desc[i,5]
        gaz 			<- mel.desc[i,1:3]
        tps.desc		<- (prof.fin-prof.init) / vit.desc

        ppart.desc		<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.desc,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.desc,Pprec=ppart.desc)
        Mval.desc		<- M.val(a,b,ppart.desc,prof.fin)
        pc.grad.desc	<- pc.grad(ppart.desc,prof.fin,Mval.desc)
        tension.desc	<- apply(ppart.desc,1,sum)

        ppart.tot 	<- cbind(ppart.tot,ppart.desc)
        Mval.tot	<- cbind(Mval.tot,Mval.desc)
        pc.grad.tot	<- cbind(pc.grad.tot,pc.grad.desc)
        tension.tot	<- cbind(tension.tot,tension.desc)
        temps.tot	<- c(temps.tot,tps.desc)
        nb.segment <- nb.segment + 1
    }

    #  Saturation au fond
    prof.init		<- mel.desc[nb.mel.desc,4]
    prof.fin		<- mel.desc[nb.mel.desc,5]
    gaz				<- as.matrix(mel.fond[1:3])
    tps.fond		<- temps - sum(temps.tot)

    ppart.fond		<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.desc,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.fond,Pprec=ppart.desc)
    Mval.fond		<- M.val(a,b,ppart.desc,prof.max)
    pc.grad.fond	<- pc.grad(ppart.fond,prof.max,Mval.fond)
    tension.fond	<- apply(ppart.fond,1,sum)

    ppart.tot		<- cbind(ppart.tot,ppart.fond)
    Mval.tot 		<- cbind(Mval.tot,Mval.fond)
    pc.grad.tot		<- cbind(pc.grad.tot,pc.grad.fond)
    tension.tot		<- cbind(tension.tot,tension.fond)
    temps.tot		<- c(temps.tot,tps.fond)
    nb.segment		<- nb.segment + 1


    # ========================================================================================
    # = ====== DETERMINATION DE LA ZONE DE DECO ET DE LA PROFONDEUR DU PREMIER PALIER ====== =
    # ========================================================================================
    # Au départ du fond, détermination de la zone oé commence la décompression. é partir de cette zone, ralentir la remontée é au moins 10 m/min
    dec.zone <- deco.zone(ppart.fond,prof.max,vit.rem,gaz,periode)

    #  Au départ du fond, détermination de la profondeur du premier palier é effectuer.
    #  Attention, pendant la remontée jusqu'au premier palier, les tissus les plus rapides désaturent. Le premier palier réel sera donc probablement moins profond et sera recalculé.
    first.stop <- GF.paliers(a,b,ppart.fond,GF.Lo,prof.max)[[1]]
    k=0
    repeat {
        k = k+1
        #  Détermination des mélanges é utiliser jusqu'é la profondeur du SECOND palier
        prof.stops <- seq(first.stop-3,0,-3)

        if (manu.up) {
            mel.rem		<- melange.matrix.paliers.2(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.up)
        } else {
            mel.rem		<- melange.matrix.paliers(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max)
        }
        mel.pal		<- mel.rem[[2]]
        mel.rem		<- mel.rem[[1]]
        nb.mel.rem 	<- nrow(mel.rem)
        ppart.rem	<- ppart.fond
        # IL EST LA LE PB !!! QUAND IL Y A PLUSIEURS MELANGES, POUR LE SECOND, LE TROISIEME ETC, IL NE FAUT PAS REPARTIR DE PPART.FOND !!!

        #  Calcul de la décompression jusqu'au SECOND palier
        for (i in 1:nb.mel.rem) {
            prof.init	<- mel.rem[i,4]
            prof.fin 	<- mel.rem[i,5]
            gaz			<- mel.rem[i,1:3]
            tps.rem		<- (prof.init-prof.fin) / vit.rem

            ppart.rem	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.rem,Pprec=ppart.rem)
            Mval.rem	<- M.val(a,b,ppart.rem,prof.fin)
            pc.grad.rem	<- pc.grad(ppart.rem,prof.fin,Mval.rem)
            tension.rem	<- apply(ppart.rem,1,sum)
        }

        #  Si la tension d'un des tissus dépasse GF.Lo, alors on est allé trop loin et il faut faire un palier é la profondeur du premier palier.
        #  Sinon, on recommence toute la procédure en montant de 3 m
        if ((max(pc.grad.rem)/100) < GF.Lo) {
            first.stop <- first.stop - 3
        } else {
            break
        }

    }

    prof.stops 	<- seq(first.stop,0,-3)

    if (manu.up) {
        mel.rem		<- melange.matrix.paliers.2(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.up)
    } else {
        mel.rem		<- melange.matrix.paliers(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max)
    }

    mel.pal		<- mel.rem[[2]]
    mel.rem		<- mel.rem[[1]]
    nb.mel.rem 	<- nrow(mel.rem)

    # =========================================================================================
    # = ================ CALCUL DE LA DECOMPRESSION JUSQU'AU PREMIER PALIER ================= =
    # =========================================================================================

    for (i in 1:nb.mel.rem) {

        prof.init	<- mel.rem[i,4]
        prof.fin 	<- mel.rem[i,5]
        gaz			<- mel.rem[i,1:3]
        tps.rem		<- (prof.init-prof.fin) / vit.rem

        ppart.rem	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.rem,Pprec=ppart.fond)
        Mval.rem	<- M.val(a,b,ppart.rem,prof.fin)
        pc.grad.rem	<- pc.grad(ppart.rem,prof.fin,Mval.rem)
        tension.rem	<- apply(ppart.rem,1,sum)

        ppart.tot 	<- cbind(ppart.tot,ppart.rem)
        Mval.tot	<- cbind(Mval.tot,Mval.rem)
        pc.grad.tot	<- cbind(pc.grad.tot,pc.grad.rem)
        tension.tot	<- cbind(tension.tot,tension.rem)
        temps.tot	<- c(temps.tot,tps.rem)
        nb.segment <- nb.segment + 1
    }

    # ==========================================================================================
    # = =================================== DEBUT DES PALIERS ================================ =
    # ==========================================================================================

    #  Determination de la profondeur de tous les paliers.
    #  La durée de chaque palier est initialement fixée é 0 minute.
    #  Les incréments de temps seront d'une minute
    minute 		<- rep(0,length(prof.stops))
    names(minute) <- prof.stops

    #  Calcul de tous les GF. Ces valeurs sont les pourcentages de gradient é ne pas dépasser pour chaque palier.
    all.GF <- seq(GF.Lo,GF.Hi,length.out=length(prof.stops))
    # if (length(all.GF)>1) all.GF <- c(all.GF[-2],GF.Hi)

    #  Détermination de la durée de chaque palier
    if (first.stop != 0) {
        ppart.up <- ppart.rem
        tps.rem	<- 3 / vit.rem
        for (i in 1:(length(prof.stops)-1)) {

            prof.init 	<- mel.pal[(2*i - 1),4]
            prof.fin	<- mel.pal[(2*i - 1),5]
            gaz 		<- mel.pal[(2*i - 1),1:3]

            ppart.palier 	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=minute[i],Pprec=ppart.up)
            Mval.palier		<- M.val(a,b,ppart.palier,prof.fin)
            pc.grad.palier	<- pc.grad(ppart.palier,prof.fin,Mval.palier)
            tension.palier	<- apply(ppart.palier,1,sum)

            prof.init 	<- mel.pal[(2*i),4]
            prof.fin	<- mel.pal[(2*i),5]
            gaz 		<- mel.pal[(2*i),1:3]

            ppart.up.new	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=(3/vit.rem),Pprec=ppart.palier)
            Mval.up			<- M.val(a,b,ppart.up.new,prof.fin)
            pc.grad.up		<- pc.grad(ppart.up.new,prof.fin,Mval.up)
            tension.up		<- apply(ppart.up.new,1,sum)

            while ((max(pc.grad.up)/100) > all.GF[i]) {
                minute[i] 	<- minute[i] + pas
                prof.init 	<- mel.pal[(2*i - 1),4]
                prof.fin	<- mel.pal[(2*i - 1),5]
                gaz 		<- mel.pal[(2*i - 1),1:3]

                ppart.palier 	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=minute[i],Pprec=ppart.up)
                Mval.palier		<- M.val(a,b,ppart.palier,prof.fin)
                pc.grad.palier	<- pc.grad(ppart.palier,prof.fin,Mval.palier)
                tension.palier	<- apply(ppart.palier,1,sum)

                prof.init 	<- mel.pal[(2*i),4]
                prof.fin	<- mel.pal[(2*i),5]
                gaz 		<- mel.pal[(2*i),1:3]

                ppart.up.new	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=(3/vit.rem),Pprec=ppart.palier)
                Mval.up			<- M.val(a,b,ppart.up.new,prof.fin)
                pc.grad.up		<- pc.grad(ppart.up.new,prof.fin,Mval.up)
                tension.up		<- apply(ppart.up.new,1,sum)
            }

            ppart.tot 		<- cbind(ppart.tot,ppart.palier,ppart.up.new)
            Mval.tot		<- cbind(Mval.tot,Mval.palier,Mval.up)
            pc.grad.tot		<- cbind(pc.grad.tot,pc.grad.palier,pc.grad.up)
            tension.tot		<- cbind(tension.tot,tension.palier,tension.up)
            nb.segment		<- nb.segment + 2
            temps.tot		<- c(temps.tot,minute[i],tps.rem)
            mess			<- "Profile de decompression"
            ppart.up		<- ppart.up.new
        }


    }

    out <- list(pressions.parielles=ppart.tot,Valeurs.Gradient=all.GF,Premier.palier=first.stop,Profondeur.paliers=prof.stops,Zone.decompression=dec.zone,Duree.paliers=minute,M.values=Mval.tot,PC.Gradient=pc.grad.tot,tensions=tension.tot,Duree.segments=temps.tot,Nombre.segments=nb.segment,Mel.descente=mel.desc,Mel.fond=mel.fond,Mel.rem=mel.rem,Mel.pal=mel.pal)
    out
}



# =================================================
# = Mélanges, profondeurs et pressions partielles =
# =================================================
# Pour un ou des mélanges, cette fonction calcule les pressions partielles d'O2 en surface et au fond, indique les profondeurs d'utilisation possibles et l'équivalent narcose.
# Le mélange entré peut étre une matrice de n x 2 gaz (pc.O2,pcHe), un vecteur (O2,He) ou un vecteur de plusieurs fractions He ou O2 pour un héliair.
prof.melange <- function(prof.max=90,melange=mel,heliair=F,gaz="O2",PO2.min=0.18,PO2.max=1.61) {

    if (is.matrix(melange) & heliair==F) {
        mel <- cbind(melange,100 - apply(melange,1,sum))
        colnames(mel) <- c("O2","He","N2")
    }

    if (is.vector(melange) & heliair==F) {
        mel <- c(melange,100-sum(melange))
        mel <- t(as.matrix(mel))
        colnames(mel) <- c("O2","He","N2")
    }


    if (heliair == T & gaz != "O2") mel <- mel.heliair(pc.He=melange)
    if (heliair == T & gaz == "O2") mel <- mel.heliair(pc.O2=melange)

    nb.mel <- dim(mel)[1]

    PO2.surf		<- round(mel[,1]/100,2)
    PO2.fond		<- round(mel[,1]*(prof.max/10+1)/100,2)
    PN2.fond		<- round(mel[,3]*(prof.max/10+1)/100,2)
    prof.narcose.eq	<- round(PN2.fond *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    prof.plafond	<- ceiling(10*(PO2.min/(mel[,1]/100)-1))
    prof.plafond[prof.plafond<0] <- 0
    prof.sol	<- floor(10*(PO2.max/(mel[,1]/100)-1))
    prof.util.max <- apply(as.matrix(prof.sol),1,min,prof.max)

    prof.maxi <- rep(prof.max,nb.mel)
    mess <- apply(as.matrix(prof.sol),1,alert.mess,prof.max)

    out <- data.frame(mel,PO2.surf,PO2.fond,prof.narcose.eq,prof.plafond,prof.sol,mess,prof.util.max)
    colnames(out) <- c("O2","He","N2","PO2.surface",paste("PO2 a",prof.max,"m"),"P.eq.narc","Prof.plafond","Prof.sol","Message","Prof.util.max")
    out
}

# ====================
# = Message d'alerte =
# ====================
# Cette fonction renvoie un message d'alerte lorsque que la profondeur max d'une plongée dépasse la profondeur possible d'utilisation d'un mélange
alert.mess <- function(prof.utilisation.max,prof.maxi) {
    mess <- if (prof.utilisation.max<prof.maxi) "Danger : prof max trop elevee !" else "OK"
}

# ===============================================
# = Calcul des proportions d'un mélange héliair =
# ===============================================
# L'argument pc.He correspond au pourcentage d'hélium dans le mélange final. Il peut s'agir d'un vecteur.
# Méme chose pour pc.O2. Seul l'un des deux arguments doit étre renseigné. Si les 2 le sont, seul He est pris en compte
mel.heliair <- function(pc.He = NULL, pc.O2 = NULL) {

    if (!is.null(pc.He)) {
        remaining 	<- 100 - pc.He
        pc.O2		<- 0.21*remaining
        pc.N2		<- 0.79*remaining
    } else {
        pc.N2		<- pc.O2*0.79/0.21
        pc.He		<- 100 - pc.N2 - pc.O2
    }

    mel				<- round(cbind(pc.O2,pc.He,pc.N2),1)
    colnames(mel)	<- c("O2","He","N2")
    mel
}


# mel est un objet fourni par la fonction prof.melange() (i.e. un data.frame contenant 9 colonnes)
melange.matrix.desc <- function(mel,p.max=90,PO2.min=0.18,PO2.max=1.61) {

    mel <- mel[order(mel[,1]),]
    mel <- mel[,-(9:10)]
    mel <- cbind(mel,p.min=mel$Prof.plafond,p.max=mel$Prof.sol)
    if (is.vector(mel)) mel <- t(as.matrix(mel))
    if (is.data.frame(mel)) mel <- as.matrix(mel)
    mel[,10] <- apply(as.matrix(mel[,10]),1,min,p.max)

    nb.gaz <- nrow(mel)

    if (nb.gaz == 1) gaz.descente <- mel

    if (nb.gaz > 1) {
        premier.gaz <- which.min(mel[,7])
        gaz.descente <- mel[premier.gaz:1,]
        if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))

        if (nrow(gaz.descente)>=2) {
            for (i in 1:(premier.gaz-1)) {
                gaz.descente[i+1,9] <- gaz.descente[i,10]
            }
        }
    }

    gaz.descente[nrow(gaz.descente),10] <- min(gaz.descente[nrow(gaz.descente),10],p.max)

    gaz.descente <- gaz.descente[,c(1:3,9,10)]
    if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))
    gaz.descente <- gaz.descente[c(1:nrow(gaz.descente),nrow(gaz.descente)),]
    gaz.descente[nrow(gaz.descente),4] <- p.max

    # On calcule ensuite la PO2 au début et é la fin de chaque segment
    PO2.init <- round((as.matrix(gaz.descente[,4])/10 + 1)*gaz.descente[,1]/100,2)
    PO2.fin  <- round((as.matrix(gaz.descente[,5])/10 + 1)*gaz.descente[,1]/100,2)
    gaz.descente <- cbind(gaz.descente,PO2.init,PO2.fin)

    # Et enfin, on calcule les profondeur narcose équivalentes au début et é la fin de chaque segment
    PN2		<- round(gaz.descente[,3]*(gaz.descente[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(gaz.descente[,3]*(gaz.descente[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    gaz.descente <- cbind(gaz.descente,prof.narcose.eq,prof.narcose.eq.f)

    for (i in 1:(nrow(gaz.descente)-1)) {
        if (identical(gaz.descente[i,],gaz.descente[i+1,])) {
            gaz.descente <- gaz.descente[-(i+1),]
        }
    }

    if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))

    # On nomme toutes les colonnes
    colnames(gaz.descente) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")
    rownames(gaz.descente) <- NULL

    gaz.descente
}

melange.matrix.desc.manuel <- function(mel,p.max=90,PO2.min=0.18,PO2.max=1.61,prof.chang=c(20,60)) {

    mel <- mel[order(mel[,1]),]
    mel <- mel[,-(9:10)]
    mel <- cbind(mel,p.min=mel$Prof.plafond,p.max=mel$Prof.sol)
    if (is.vector(mel)) mel <- t(as.matrix(mel))
    if (is.data.frame(mel)) mel <- as.matrix(mel)
    mel[,10] <- apply(as.matrix(mel[,10]),1,min,p.max)

    nb.gaz <- nrow(mel)

    if (nb.gaz == 1) gaz.descente <- mel

    if (nb.gaz > 1) {
        premier.gaz <- which.min(mel[,7])
        gaz.descente <- mel[premier.gaz:1,]
        if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))

        if (nrow(gaz.descente)>=2) {
            gaz.descente[-(nrow(gaz.descente)),10] <- prof.chang
            for (i in 1:(premier.gaz-1)) {
                gaz.descente[i+1,9] <- gaz.descente[i,10]
            }
        }
    }

    gaz.descente[nrow(gaz.descente),10] <- min(gaz.descente[nrow(gaz.descente),10],p.max)

    gaz.descente <- gaz.descente[,c(1:3,9,10)]
    if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))
    gaz.descente <- gaz.descente[c(1:nrow(gaz.descente),nrow(gaz.descente)),]
    gaz.descente[nrow(gaz.descente),4] <- p.max

    # On calcule ensuite la PO2 au début et é la fin de chaque segment
    PO2.init <- round((as.matrix(gaz.descente[,4])/10 + 1)*gaz.descente[,1]/100,2)
    PO2.fin  <- round((as.matrix(gaz.descente[,5])/10 + 1)*gaz.descente[,1]/100,2)
    gaz.descente <- cbind(gaz.descente,PO2.init,PO2.fin)

    # Et enfin, on calcule les profondeur narcose équivalentes au début et é la fin de chaque segment
    PN2		<- round(gaz.descente[,3]*(gaz.descente[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(gaz.descente[,3]*(gaz.descente[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    gaz.descente <- cbind(gaz.descente,prof.narcose.eq,prof.narcose.eq.f)

    for (i in 1:(nrow(gaz.descente)-1)) {
        if (identical(gaz.descente[i,],gaz.descente[i+1,])) {
            gaz.descente <- gaz.descente[-(i+1),]
        }
    }

    if (is.vector(gaz.descente)) gaz.descente <- t(as.matrix(gaz.descente))

    # On nomme toutes les colonnes
    colnames(gaz.descente) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")
    rownames(gaz.descente) <- NULL

    gaz.descente
}

melange.matrix.paliers <- function(mel,p.max=90,prof.paliers,PO2.min=0.18,PO2.max=1.61) {

    # On trie les gaz par ordre croissant de fraction d'oxygéne et on supprime la colonne des erreurs
    mel <- mel[order(mel[,1]),]
    mel <- mel[,-9]
    mel <- as.matrix(mel)

    # Mise en forme de variables utiles
    prof.paliers <- rev(prof.paliers)
    nb.gaz <- nrow(mel)
    mel.rem <- NULL
    mel.mat <- NULL

    # Création de la matrice des gaz é respirer é la remontée jusqu'au premier palier et de la matrice des gaz qui seront respirés lors des paliers
    # Il est possible d'effectuer jusqu'é deux changements de gaz entre le fond et le premier palier
    if (nb.gaz > 2) {
        if(max(prof.paliers) <= mel[3,8]) {
            mel.rem <- mel[1:3,]
            mel.pal	<- mel[3:nb.gaz,]
        } else {
            if (max(prof.paliers) <= mel[2,8]) {
                mel.rem <- mel[1:2,]
                mel.pal	<- mel[2:nb.gaz,]
            } else {
                mel.rem <- mel[1,]
                mel.pal <- mel
            }
        }
    } else {
        if (nb.gaz > 1) {
            if (max(prof.paliers) <= mel[2,8]) {
                mel.rem <- mel[1:2,]
                mel.pal	<- mel[2:nb.gaz,]
            } else {
                mel.rem <- mel[1,]
                mel.pal <- mel
            }
        } else {
            if (nb.gaz == 1) {
                mel.rem <- mel
                mel.pal <- mel
            }
        }
    }

    if (is.vector(mel.rem)) mel.rem <- t(as.matrix(mel.rem))
    if (is.vector(mel.pal)) mel.pal <- t(as.matrix(mel.pal))
    nb.gaz.pal <- nrow(mel.pal)
    nb.gaz.rem <- nrow(mel.rem)

    if (nb.gaz.pal == 1) 	mel.mat <- matrix(as.vector(rep(mel.pal,each=(length(prof.paliers)))),ncol=ncol(mel))

    # Création de la matrice des gaz é respirer pour chaque palier et entre les paliers
    prof.palier.max <- apply(as.matrix(mel.pal[,8]),1,min,max(prof.paliers))
    prof.palier.max <- c(prof.palier.max,0)
    prof.palier.max <- prof.palier.max - prof.palier.max%%3

    if (nb.gaz.pal > 1) {
        index.limite <- numeric(nb.gaz.pal)
        for (i in 1:length(prof.palier.max)) index.limite[i] <- which(prof.paliers==prof.palier.max[i])

        for (i in nb.gaz.pal:1) {
            nb.rep.mel	<- if (i == nb.gaz.pal) index.limite[i] else index.limite[i] - index.limite[i+1]
            mel.matrice	<- matrix(mel.pal[i,],nrow=nb.rep.mel,ncol=ncol(mel),byrow=T)
            mel.mat		<- rbind(mel.matrice,mel.mat)
        }
    }

    # On supprime les colonnes inutiles de mel.mat, on y ajoute 2 colonnes contenant les profondeurs des paliers et on copie chaque ligne 2 fois
    prof.paliers <- rev(prof.paliers)
    mel.mat <- cbind(mel.mat[,1:3],as.matrix(prof.paliers),as.matrix(prof.paliers))
    mel.rem <- mel.rem[,c(1:3,9,9)]
    if (is.vector(mel.rem)) mel.rem <- t(as.matrix(mel.rem))
    if (is.vector(mel.pal)) mel.pal <- t(as.matrix(mel.pal))

    if (nb.gaz.rem == 1) {
        mel.rem[1,5] <- max(prof.paliers)
    } else {
        for (i in 1:(nb.gaz.rem-1)) {
            mel.rem[i,5] <- mel.rem[i+1,4]
        }
        mel.rem[nb.gaz.rem,5] <- max(prof.paliers)
    }

    mel.mat <- mel.mat[rep(1:nrow(mel.mat),each=2),]
    nb.seg <- nrow(mel.mat)

    # On modifie la colonne 5 pour qu'elle contienne les profondeurs de fin de chaque segment de la remontée
    # La colonne 4 reste inchangée et contient les profondeurs initiales de chaque segment de remontée
    p.init <- mel.mat[,5]
    p.init <- c(p.init[-1],0)
    mel.mat[,5] <- p.init

    # On ajoute une ligne correspondant au segment de la remontée et on supprime les 2 derniéres lignes inutiles
    mel.mat <- mel.mat[-c(nrow(mel.mat)-1,nrow(mel.mat)),]
    # mel.mat <- mel.mat[c(1,1:nrow(mel.mat)),]
    # mel.mat[1,4] <- min(mel.rem[,5])

    # On calcule ensuite la PO2 au début et é la fin de chaque segment
    PO2.init <- round((as.matrix(mel.mat[,4])/10 + 1)*mel.mat[,1]/100,2)
    PO2.fin  <- round((as.matrix(mel.mat[,5])/10 + 1)*mel.mat[,1]/100,2)
    mel.mat <- cbind(mel.mat,PO2.init,PO2.fin)

    PO2.init <- round((as.matrix(mel.rem[,4])/10 + 1)*mel.rem[,1]/100,2)
    PO2.fin  <- round((as.matrix(mel.rem[,5])/10 + 1)*mel.rem[,1]/100,2)
    mel.rem <- cbind(mel.rem,PO2.init,PO2.fin)


    # # Et enfin, on calcule les profondeur narcose équivalentes au début et é la fin de chaque segment
    PN2		<- round(mel.mat[,3]*(mel.mat[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(mel.mat[,3]*(mel.mat[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    mel.mat <- cbind(mel.mat,prof.narcose.eq,prof.narcose.eq.f)

    PN2		<- round(mel.rem[,3]*(mel.rem[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(mel.rem[,3]*(mel.rem[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    mel.rem <- cbind(mel.rem,prof.narcose.eq,prof.narcose.eq.f)

    # # On nomme toutes les colonnes
    colnames(mel.mat) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")
    colnames(mel.rem) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")

    out <- list(remontee=mel.rem,paliers=mel.mat)
    out

}

melange.matrix.paliers.2 <- function(mel,p.max=90,prof.paliers,PO2.min=0.18,PO2.max=1.61,prof.chang=c(60,25,6)) {

    # On trie les gaz par ordre croissant de fraction d'oxygéne et on supprime la colonne des erreurs
    mel <- mel[order(mel[,1]),]
    mel <- mel[,-9]
    mel <- as.matrix(mel)

    if (nrow(mel)>1) {
        mel[2:nrow(mel),8] <- prof.chang
        mel[2:nrow(mel),9] <- prof.chang
    }

    # Mise en forme de variables utiles
    prof.paliers <- rev(prof.paliers)
    nb.gaz <- nrow(mel)
    mel.rem <- NULL
    mel.mat <- NULL

    # Création de la matrice des gaz é respirer é la remontée jusqu'au premier palier et de la matrice des gaz qui seront respirés lors des paliers
    # Il est possible d'effectuer jusqu'é deux changements de gaz entre le fond et le premier palier
    if (nb.gaz > 2) {
        if(max(prof.paliers) <= mel[3,8]) {
            mel.rem <- mel[1:3,]
            mel.pal	<- mel[3:nb.gaz,]
        } else {
            if (max(prof.paliers) <= mel[2,8]) {
                mel.rem <- mel[1:2,]
                mel.pal	<- mel[2:nb.gaz,]
            } else {
                mel.rem <- mel[1,]
                mel.pal <- mel
            }
        }
    } else {
        if (nb.gaz > 1) {
            if (max(prof.paliers) <= mel[2,8]) {
                mel.rem <- mel[1:2,]
                mel.pal	<- mel[2:nb.gaz,]
            } else {
                mel.rem <- mel[1,]
                mel.pal <- mel
            }
        } else {
            if (nb.gaz == 1) {
                mel.rem <- mel
                mel.pal <- mel
            }
        }
    }

    if (is.vector(mel.rem)) mel.rem <- t(as.matrix(mel.rem))
    if (is.vector(mel.pal)) mel.pal <- t(as.matrix(mel.pal))
    nb.gaz.pal <- nrow(mel.pal)
    nb.gaz.rem <- nrow(mel.rem)

    if (nb.gaz.pal == 1) 	mel.mat <- matrix(as.vector(rep(mel.pal,each=(length(prof.paliers)))),ncol=ncol(mel))

    # Création de la matrice des gaz é respirer pour chaque palier et entre les paliers
    prof.palier.max <- apply(as.matrix(mel.pal[,8]),1,min,max(prof.paliers))
    prof.palier.max <- c(prof.palier.max,0)
    prof.palier.max <- prof.palier.max - prof.palier.max%%3

    if (nb.gaz.pal > 1) {
        index.limite <- numeric(nb.gaz.pal)
        for (i in 1:length(prof.palier.max)) index.limite[i] <- which(prof.paliers==prof.palier.max[i])

        for (i in nb.gaz.pal:1) {
            nb.rep.mel	<- if (i == nb.gaz.pal) index.limite[i] else index.limite[i] - index.limite[i+1]
            mel.matrice	<- matrix(mel.pal[i,],nrow=nb.rep.mel,ncol=ncol(mel),byrow=T)
            mel.mat		<- rbind(mel.matrice,mel.mat)
        }
    }

    # On supprime les colonnes inutiles de mel.mat, on y ajoute 2 colonnes contenant les profondeurs des paliers et on copie chaque ligne 2 fois
    prof.paliers <- rev(prof.paliers)
    mel.mat <- cbind(mel.mat[,1:3],as.matrix(prof.paliers),as.matrix(prof.paliers))
    mel.rem <- mel.rem[,c(1:3,9,9)]
    if (is.vector(mel.rem)) mel.rem <- t(as.matrix(mel.rem))
    if (is.vector(mel.pal)) mel.pal <- t(as.matrix(mel.pal))

    if (nb.gaz.rem == 1) {
        mel.rem[1,5] <- max(prof.paliers)
    } else {
        for (i in 1:(nb.gaz.rem-1)) {
            mel.rem[i,5] <- mel.rem[i+1,4]
        }
        mel.rem[nb.gaz.rem,5] <- max(prof.paliers)
    }

    mel.mat <- mel.mat[rep(1:nrow(mel.mat),each=2),]
    nb.seg <- nrow(mel.mat)

    # On modifie la colonne 5 pour qu'elle contienne les profondeurs de fin de chaque segment de la remontée
    # La colonne 4 reste inchangée et contient les profondeurs initiales de chaque segment de remontée
    p.init <- mel.mat[,5]
    p.init <- c(p.init[-1],0)
    mel.mat[,5] <- p.init

    # On ajoute une ligne correspondant au segment de la remontée et on supprime les 2 derniéres lignes inutiles
    mel.mat <- mel.mat[-c(nrow(mel.mat)-1,nrow(mel.mat)),]
    # mel.mat <- mel.mat[c(1,1:nrow(mel.mat)),]
    # mel.mat[1,4] <- min(mel.rem[,5])

    # On calcule ensuite la PO2 au début et é la fin de chaque segment
    PO2.init <- round((as.matrix(mel.mat[,4])/10 + 1)*mel.mat[,1]/100,2)
    PO2.fin  <- round((as.matrix(mel.mat[,5])/10 + 1)*mel.mat[,1]/100,2)
    mel.mat <- cbind(mel.mat,PO2.init,PO2.fin)

    PO2.init <- round((as.matrix(mel.rem[,4])/10 + 1)*mel.rem[,1]/100,2)
    PO2.fin  <- round((as.matrix(mel.rem[,5])/10 + 1)*mel.rem[,1]/100,2)
    mel.rem <- cbind(mel.rem,PO2.init,PO2.fin)


    # # Et enfin, on calcule les profondeur narcose équivalentes au début et é la fin de chaque segment
    PN2		<- round(mel.mat[,3]*(mel.mat[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(mel.mat[,3]*(mel.mat[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    mel.mat <- cbind(mel.mat,prof.narcose.eq,prof.narcose.eq.f)

    PN2		<- round(mel.rem[,3]*(mel.rem[,4]/10+1)/100,2)
    prof.narcose.eq	<- round(PN2 *10 / 0.79 - 10,1)
    prof.narcose.eq[prof.narcose.eq<0] <- 0

    PN2.f		<- round(mel.rem[,3]*(mel.rem[,5]/10+1)/100,2)
    prof.narcose.eq.f	<- round(PN2.f *10 / 0.79 - 10,1)
    prof.narcose.eq.f[prof.narcose.eq.f<0] <- 0

    mel.rem <- cbind(mel.rem,prof.narcose.eq,prof.narcose.eq.f)

    # # On nomme toutes les colonnes
    colnames(mel.mat) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")
    colnames(mel.rem) <- c("O2","He","N2","Prof.init","Prof.fin","PO2.init","PO2.fin","Prof.narc.eq.init","Prof.narc.eq.fin")

    out <- list(remontee=mel.rem,paliers=mel.mat)
    out

}

# ====================================================================================
# = Fonction qui renvoie uniquement les infos concernant les segments de la descente =
# ====================================================================================
descente <- function(ex) {
    out <- ex$Mel.descente[-nrow(ex$Mel.descente),]
    out
}

# =====================================================================================================================
# = Fonction qui produit un tableau synthétique reprenant des infos essentielles pour tous les segments d'une plongée =
# =====================================================================================================================
profil <- function(ex) {
    out <- rbind(ex$Mel.descente,ex$Mel.rem,ex$Mel.pal)
    out <- cbind(out[,1:5],Duree=round(ex$Duree.segments,1),out[,6:9])
    Run.Time <- floor(cumsum(out[,6]))
    out <- cbind(out[,1:6],RT.i=c(0,Run.Time[-length(Run.Time)]),RT.f=Run.Time,out[,7:10])
    rownames(out) <- NULL
    out
}

# ====================================================================================
# = Fonction qui calcule, pour chaque gaz, le volume consommé au cours de la plongée =
# ====================================================================================
conso <- function(pro,conso.fond,conso.pal,nb.segment,nb.seg.pal) {
    pres.abs <- (abs(pro[,4] - pro[,5])/2 + pmin(pro[,4],pro[,5]))/10 +1
    conso <- rep(c(conso.fond,conso.pal),c(nb.segment-nb.seg.pal-1,nb.seg.pal+1))
    conso.seg <- ceiling(pres.abs * pro[,6] * conso)
    mel <- unique(pro[,1:3])
    mel <- mel[order(mel[,1]),]
    n.mel <- nrow(mel)
    vol.mel <- numeric(n.mel)
    for (i in 1:n.mel) {
        vol.mel[i] <- sum(conso.seg[pro[,1]==mel[i,1]])

    }
    noms.mel <- paste(mel[,1],"-",mel[,2],"-",mel[,3])
    names(vol.mel) <- noms.mel
    vol.mel
}

# =============================================================================================================================
# = Fonction qui renvoie M-val, tensions, pourcentage de gradient et pourcentage de M-value pour chaque segment d'une plongée =
# =============================================================================================================================
analyse <- function(ex) {

    Mval <- apply(ex$M.values,2,max)
    Tensions <- apply(ex$tensions,2,max)
    Pc.Grad <- apply(ex$PC.Gradient,2,max)
    Pc.Val <- ex$tensions * 100 / ex$M.values
    Pc.Mval <- apply(Pc.Val,2,max)

    out <- cbind(Mval,Tensions,Pc.Grad,Pc.Mval)
    rownames(out) <- names(ex$Duree.segments)
    out

}

# ==============================================================================
# = Calcul de la décompression quand un gaz est perdu au moment de la remontée =
# ==============================================================================
# Cette fonction est trés similaire é BIGONE. Seuls les gaz disponibles au départ du fond modifiés
Perte.gaz <- function (prof.max=90,temps=20,vit.desc=30,vit.rem=10,gradient=c(0.20,0.75),gaz=rbind(c(15,53),c(40,0)),durete=3,PO2.min=0.18,PO2.max=1.61,gaz.perdu=1,pas=0.5,manu.down=F,gaz.change.down=c(10,45),manu.up=F,gaz.change.up=c(60,30,6)){

    # ================================================================================================
    # = ===================================== CREATION DES VARIABLES UTILES ======================== =
    # ================================================================================================

    # Data frame des M-Values pour l'azote et l'hélium pour les 16 compartiments (on utilise le 1b et non le 1)
    ZHL16.C <- data.frame(
        comp = 1:16,
        N2.periode = c(5,8,12.5,18.5,27,38.3,54.3,77,109,146,187,239,305,390,498,635),
        a.N2 = c(11.696,10.000,8.618,7.562,6.667,5.600,4.947,4.500,4.187,3.798,3.497,3.223,2.850,2.737,2.523,2.327),
        b.N2 = c(0.5577867,0.6513809,0.7221781,0.7824726,0.8126117,0.8433837,0.8692629,0.8910274,0.9091736,0.9221689,0.9318796,0.9402915,0.9476876,0.9543806,0.9602458,0.9653441),
        He.periode = c(1.88,3.02,4.72,6.99,10.21,14.48,20.53,29.11,41.20,55.19,70.69,90.34,115.29,147.42,188.24,240.03),
        a.He = c(16.189,13.830,11.919,10.458,9.22,8.205,7.305,6.502,5.95,5.545,5.333,5.189,5.181,5.176,5.172,5.119),
        b.He = c(0.4770082,0.5747126,0.6526989,0.7222824,0.7582076,0.7956715,0.8278831,0.8552857,0.8757334,0.8903134,0.8996851,0.9072764,0.9121591,0.9170946,0.9216590,0.9266982)
    )

    # Recupération des données de periode et des coef a et b pour les 2 gaz
    periode <- cbind(ZHL16.C$He.periode,ZHL16.C$N2.periode)
    a <- cbind(ZHL16.C$a.He,ZHL16.C$a.N2)
    b <- cbind(ZHL16.C$b.He,ZHL16.C$b.N2)

    # Initialisation des paramétres de gradient et de profondeur
    GF.Lo <- gradient[1]
    GF.Hi <- gradient[2]

    # Création de la matrice des mélanges gazeux
    if (is.vector(gaz)) gaz <- t(as.matrix(gaz))
    gaz			<- gaz[order(gaz[,1]),]
    mel		 	<- prof.melange(prof.max=prof.max,melange=gaz,heliair=F,PO2.min=PO2.min,PO2.max=PO2.max)

    if (manu.down)	{
        mel.desc 	<- melange.matrix.desc.manuel(mel,p.max=prof.max,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.down)
    } else {	mel.desc 	<- melange.matrix.desc(mel,p.max=prof.max,PO2.min=PO2.min,PO2.max=PO2.max) }

    nb.mel.desc	<- nrow(mel.desc)
    mel.fond	<- mel.desc[nb.mel.desc,]
    if (is.vector(mel.desc)) mel.desc <- t(as.matrix(mel.desc))
    if (is.vector(mel.fond)) mel.fond <- t(as.matrix(mel.fond))

    # Création et initialisation des objets oé seront stockées les valeurs de chaque segment.
    ppart.tot 	<- NULL
    Mval.tot  	<- NULL
    pc.grad.tot	<- NULL
    tension.tot	<- NULL
    temps.tot	<- NULL
    nb.segment <- 0


    # ==========================================================================================
    # = =================================== DEBUT DE LA PLONGEE ============================== =
    # ==========================================================================================
    # Saturation é la descente
    ppart.desc <- NULL
    for (i in 1:(nb.mel.desc-1)) {
        prof.init 		<- mel.desc[i,4]
        prof.fin 		<- mel.desc[i,5]
        gaz 			<- mel.desc[i,1:3]
        tps.desc		<- (prof.fin-prof.init) / vit.desc

        ppart.desc		<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.desc,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.desc,Pprec=ppart.desc)
        Mval.desc		<- M.val(a,b,ppart.desc,prof.fin)
        pc.grad.desc	<- pc.grad(ppart.desc,prof.fin,Mval.desc)
        tension.desc	<- apply(ppart.desc,1,sum)

        ppart.tot 	<- cbind(ppart.tot,ppart.desc)
        Mval.tot	<- cbind(Mval.tot,Mval.desc)
        pc.grad.tot	<- cbind(pc.grad.tot,pc.grad.desc)
        tension.tot	<- cbind(tension.tot,tension.desc)
        temps.tot	<- c(temps.tot,tps.desc)
        nb.segment <- nb.segment + 1
    }

    #  Saturation au fond
    prof.init		<- mel.desc[nb.mel.desc,4]
    prof.fin		<- mel.desc[nb.mel.desc,5]
    gaz				<- as.matrix(mel.fond[1:3])
    tps.fond		<- temps - sum(temps.tot)

    ppart.fond		<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.desc,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.fond,Pprec=ppart.desc)
    Mval.fond		<- M.val(a,b,ppart.desc,prof.max)
    pc.grad.fond	<- pc.grad(ppart.fond,prof.max,Mval.fond)
    tension.fond	<- apply(ppart.fond,1,sum)

    ppart.tot		<- cbind(ppart.tot,ppart.fond)
    Mval.tot 		<- cbind(Mval.tot,Mval.fond)
    pc.grad.tot		<- cbind(pc.grad.tot,pc.grad.fond)
    tension.tot		<- cbind(tension.tot,tension.fond)
    temps.tot		<- c(temps.tot,tps.fond)
    nb.segment		<- nb.segment + 1


    # ========================================================================================
    # = ====== DETERMINATION DE LA ZONE DE DECO ET DE LA PROFONDEUR DU PREMIER PALIER ====== =
    # ========================================================================================
    # Au départ du fond, détermination de la zone oé commence la décompression. é partir de cette zone, ralentir la remontée é au moins 10 m/min
    dec.zone <- deco.zone(ppart.fond,prof.max,vit.rem,gaz,periode)

    #  Au départ du fond, détermination de la profondeur du premier palier é effectuer.
    #  Attention, pendant la remontée jusqu'au premier palier, les tissus les plus rapides désaturent. Le premier palier réel sera donc probablement moins profond et sera recalculé.
    first.stop <- GF.paliers(a,b,ppart.fond,GF.Lo,prof.max)[[1]]

    #  Perte d'un gaz de déco au moment oé on quitte le fond
    mel <- mel[-(1+gaz.perdu),]

    k=0
    repeat {
        k = k+1
        #  Détermination des mélanges é utiliser jusqu'é la profondeur du SECOND palier
        prof.stops <- seq(first.stop-3,0,-3)
        if (manu.up) {
            mel.rem		<- melange.matrix.paliers.2(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.up)
        } else {
            mel.rem		<- melange.matrix.paliers(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max)
        }
        mel.pal		<- mel.rem[[2]]
        mel.rem		<- mel.rem[[1]]
        nb.mel.rem 	<- nrow(mel.rem)
        ppart.rem	<- ppart.fond
        # IL EST LA LE PB !!! QUAND IL Y A PLUSIEURS MELANGES, POUR LE SECOND, LE TROISIEME ETC, IL NE FAUT PAS REPARTIR DE PPART.FOND !!!

        #  Calcul de la décompression jusqu'au SECOND palier
        for (i in 1:nb.mel.rem) {
            prof.init	<- mel.rem[i,4]
            prof.fin 	<- mel.rem[i,5]
            gaz			<- mel.rem[i,1:3]
            tps.rem		<- (prof.init-prof.fin) / vit.rem

            ppart.rem	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.rem,Pprec=ppart.rem)
            Mval.rem	<- M.val(a,b,ppart.rem,prof.fin)
            pc.grad.rem	<- pc.grad(ppart.rem,prof.fin,Mval.rem)
            tension.rem	<- apply(ppart.rem,1,sum)
        }

        #  Si la tension d'un des tissus dépasse GF.Lo, alors on est allé trop loin et il faut faire un palier é la profondeur du premier palier.
        #  Sinon, on recommence toute la procédure en montant de 3 m
        if ((max(pc.grad.rem)/100) < GF.Lo) {
            first.stop <- first.stop - 3
        } else {
            break
        }

    }

    prof.stops 	<- seq(first.stop,0,-3)
    if (manu.up) {
        mel.rem		<- melange.matrix.paliers.2(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max,prof.chang=gaz.change.up)
    } else {
        mel.rem		<- melange.matrix.paliers(mel,p.max=prof.max,prof.paliers=prof.stops,PO2.min=PO2.min,PO2.max=PO2.max)
    }
    mel.pal		<- mel.rem[[2]]
    mel.rem		<- mel.rem[[1]]
    nb.mel.rem 	<- nrow(mel.rem)

    # =========================================================================================
    # = ================ CALCUL DE LA DECOMPRESSION JUSQU'AU PREMIER PALIER ================= =
    # =========================================================================================

    for (i in 1:nb.mel.rem) {

        prof.init	<- mel.rem[i,4]
        prof.fin 	<- mel.rem[i,5]
        gaz			<- mel.rem[i,1:3]
        tps.rem		<- (prof.init-prof.fin) / vit.rem

        ppart.rem	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=tps.rem,Pprec=ppart.fond)
        Mval.rem	<- M.val(a,b,ppart.rem,prof.fin)
        pc.grad.rem	<- pc.grad(ppart.rem,prof.fin,Mval.rem)
        tension.rem	<- apply(ppart.rem,1,sum)

        ppart.tot 	<- cbind(ppart.tot,ppart.rem)
        Mval.tot	<- cbind(Mval.tot,Mval.rem)
        pc.grad.tot	<- cbind(pc.grad.tot,pc.grad.rem)
        tension.tot	<- cbind(tension.tot,tension.rem)
        temps.tot	<- c(temps.tot,tps.rem)
        nb.segment <- nb.segment + 1
    }

    # ==========================================================================================
    # = =================================== DEBUT DES PALIERS ================================ =
    # ==========================================================================================

    #  Determination de la profondeur de tous les paliers.
    #  La durée de chaque palier est initialement fixée é 0 minute.
    #  Les incréments de temps seront d'une minute
    minute 		<- rep(0,length(prof.stops))
    names(minute) <- prof.stops

    #  Calcul de tous les GF. Ces valeurs sont les pourcentages de gradient é ne pas dépasser pour chaque palier.
    all.GF <- seq(GF.Lo,GF.Hi,length.out=length(prof.stops))
    # if (length(all.GF)>1) all.GF <- c(all.GF[-2],GF.Hi)

    #  Détermination de la durée de chaque palier
    if (first.stop != 0) {
        ppart.up <- ppart.rem
        tps.rem	<- 3 / vit.rem
        for (i in 1:(length(prof.stops)-1)) {

            prof.init 	<- mel.pal[(2*i - 1),4]
            prof.fin	<- mel.pal[(2*i - 1),5]
            gaz 		<- mel.pal[(2*i - 1),1:3]

            ppart.palier 	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=minute[i],Pprec=ppart.up)
            Mval.palier		<- M.val(a,b,ppart.palier,prof.fin)
            pc.grad.palier	<- pc.grad(ppart.palier,prof.fin,Mval.palier)
            tension.palier	<- apply(ppart.palier,1,sum)

            prof.init 	<- mel.pal[(2*i),4]
            prof.fin	<- mel.pal[(2*i),5]
            gaz 		<- mel.pal[(2*i),1:3]

            ppart.up.new	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=(3/vit.rem),Pprec=ppart.palier)
            Mval.up			<- M.val(a,b,ppart.up.new,prof.fin)
            pc.grad.up		<- pc.grad(ppart.up.new,prof.fin,Mval.up)
            tension.up		<- apply(ppart.up.new,1,sum)

            while ((max(pc.grad.up)/100) > all.GF[i]) {
                minute[i] 	<- minute[i] + pas
                prof.init 	<- mel.pal[(2*i - 1),4]
                prof.fin	<- mel.pal[(2*i - 1),5]
                gaz 		<- mel.pal[(2*i - 1),1:3]

                ppart.palier 	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=minute[i],Pprec=ppart.up)
                Mval.palier		<- M.val(a,b,ppart.palier,prof.fin)
                pc.grad.palier	<- pc.grad(ppart.palier,prof.fin,Mval.palier)
                tension.palier	<- apply(ppart.palier,1,sum)

                prof.init 	<- mel.pal[(2*i),4]
                prof.fin	<- mel.pal[(2*i),5]
                gaz 		<- mel.pal[(2*i),1:3]

                ppart.up.new	<- both(half.t=periode,prof.init=prof.init,prof.fin=prof.fin,vitesse=vit.rem,gaz=gaz,molec=c("azote","helium"),dur.qr=durete,t=(3/vit.rem),Pprec=ppart.palier)
                Mval.up			<- M.val(a,b,ppart.up.new,prof.fin)
                pc.grad.up		<- pc.grad(ppart.up.new,prof.fin,Mval.up)
                tension.up		<- apply(ppart.up.new,1,sum)
            }

            ppart.tot 		<- cbind(ppart.tot,ppart.palier,ppart.up.new)
            Mval.tot		<- cbind(Mval.tot,Mval.palier,Mval.up)
            pc.grad.tot		<- cbind(pc.grad.tot,pc.grad.palier,pc.grad.up)
            tension.tot		<- cbind(tension.tot,tension.palier,tension.up)
            nb.segment		<- nb.segment + 2
            temps.tot		<- c(temps.tot,minute[i],tps.rem)
            mess			<- "Profile de decompression"
            ppart.up		<- ppart.up.new
        }
        # recover()
        # ===========================================================================================
        # = =================================== RESULTATS DES CALCULS ============================= =
        # ===========================================================================================



    } # else {
    # 		Mval.up		<- M.val(a,b,ppart.up,0)
    # 		Mval 		<- cbind(Mval,Mval.up)
    # 		mess		<- "Plongee sans palier"
    # 		}

    out <- list(pressions.parielles=ppart.tot,Valeurs.Gradient=all.GF,Premier.palier=first.stop,Profondeur.paliers=prof.stops,Zone.decompression=dec.zone,Duree.paliers=minute,M.values=Mval.tot,PC.Gradient=pc.grad.tot,tensions=tension.tot,Duree.segments=temps.tot,Nombre.segments=nb.segment,Mel.descente=mel.desc,Mel.fond=mel.fond,Mel.rem=mel.rem,Mel.pal=mel.pal)
    out

}

# ==================================================================================================
# = Fonction qui produit un tableau récapitulatif quand plusieurs profils de plongée sont calculés =
# ==================================================================================================
recap <- function(ex.tot,pro.tot,nb.gaz,prof,depassement.prof,temps,depassement.temps) {

    plong.temps.tot	<- numeric()
    plong.prem.pal	<- numeric()
    plong.tps.pal	<- numeric()
    plong.PO2.fond	<- numeric()
    plong.PO2.max	<- numeric()
    plong.narc.fond	<- numeric()
    plong.narc.max	<- numeric()

    noms.mel <- paste("deco",1:(nb.gaz-1))

    plong.nb	<- 1:(nb.gaz+6)
    plong.prof	<- c(rep(prof,nb.gaz),prof+depassement.prof,prof,prof+depassement.prof,prof-depassement.prof,prof,prof-depassement.prof)
    plong.temps	<- c(rep(temps,nb.gaz+1),temps+depassement.temps,temps+depassement.temps,temps,temps-depassement.temps,temps-depassement.temps)
    plong.perdu <- c("aucun",noms.mel,rep("aucun",6))

    for (i in 1:(nb.gaz+6)) {
        plong.temps.tot <- c(plong.temps.tot,max(pro.tot[[i]][,8]))
        plong.prem.pal  <- c(plong.prem.pal,ex.tot[[i]][[3]])
        plong.tps.pal	<- c(plong.tps.pal,sum(ex.tot[[i]]$Duree.pal))
        plong.PO2.fond	<- c(plong.PO2.fond,pro.tot[[i]][which.max(pro.tot[[i]][,4]),9])
        plong.PO2.max	<- c(plong.PO2.max,max(pro.tot[[i]][,9:10]))
        plong.narc.fond	<- c(plong.narc.fond,pro.tot[[i]][which.max(pro.tot[[i]][,4]),11])
        plong.narc.max	<- c(plong.narc.max,max(pro.tot[[i]][,11:12]))
    }

    recap.plongees <- data.frame(plong.prof,plong.temps,plong.perdu,plong.temps.tot,plong.tps.pal,plong.prem.pal,plong.PO2.fond,plong.PO2.max,plong.narc.fond,plong.narc.max)
    colnames(recap.plongees) <- c("Pmax","Tps","Gaz perdu","Tps tot","Tps pal","Pr pal 1","PO2 fd","PO2 max","Narc fd","Narc max")
    rownames(recap.plongees) <- paste("Pl.",1:(nb.gaz+6))
    recap.plongees
}

# ===========================================================================
# = Fonction qui compte le nombre de gaz disponibles au cours de la plongée =
# ===========================================================================
nbgaz <- function(gaz) {
    if (is.vector(gaz)) gaz <- t(as.matrix(gaz))
    nb.gaz=nrow(gaz)
}
# gaz	<- as.matrix(gaz[order(gaz[,1]),])


# ================================
# = Création de toutes les décos =
# ================================
toutes.decos <- function(prof,temps,vitesse.descente,vitesse.remontee,gradient,gaz,nb.gaz,Durete,PO2.min,PO2.max,pas,depassement.temps,depassement.prof,manu.down,gaz.change.down,manu.up,gaz.change.up) {

    ex.tot <-  list()

    ex.tot[[1]] <- BIGONE(prof.max=prof,temps=temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    if (nb.gaz >1) {
        for (i in 1:(nb.gaz-1)) {
            ex.tot[[i+1]] <- Perte.gaz(prof.max=prof,temps=temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,gaz.perdu=i,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up[-i])
        }
    }

    ex.tot[[nb.gaz+1]] <- BIGONE(prof.max=prof+depassement.prof,temps=temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot[[nb.gaz+2]] <- BIGONE(prof.max=prof,temps=temps+depassement.temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot[[nb.gaz+3]] <- BIGONE(prof.max=prof+depassement.prof,temps=temps+depassement.temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot[[nb.gaz+4]] <- BIGONE(prof.max=prof-depassement.prof,temps=temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot[[nb.gaz+5]] <- BIGONE(prof.max=prof,temps=temps-depassement.temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot[[nb.gaz+6]] <- BIGONE(prof.max=prof-depassement.prof,temps=temps-depassement.temps,vit.desc=vitesse.descente,vit.rem=vitesse.remontee,gradient=gradient,gaz=gaz,durete=Durete,PO2.min=PO2.min,PO2.max=PO2.max,pas=pas,manu.down=manu.down,gaz.change.down=gaz.change.down,manu.up=manu.up,gaz.change.up=gaz.change.up)

    ex.tot

}



# =======================
# = Calculs de gonflage =
# =======================
calc.pression <- function(gaz,press.service) {

    if (is.vector(gaz)) gaz <- t(as.matrix(gaz))
    gaz	<- as.matrix(gaz[order(gaz[,1]),])
    gaz <- cbind(gaz,100-apply(gaz,1,sum)) # é ce stade, un ligne = un mélange et les mélanges sont classés du moins oxygéné au plus oxygéné

    nb.gaz <- nrow(gaz)

    for (i in 1:nb.gaz) {
        gg <- gaz[i]
    }

}

ggg <- function(gg=c(100,0,0),press.service=220,stick=T,compo.reste=c(0,0,0),pression.reste=0) {

    vol.reste <- compo.reste * pression.reste
    vol.souhaite <- gg * press.service
    vol.ness <- vol.souhaite - vol.reste

    gg <- vol.ness *100 / sum(vol.ness)
    press.service <- press.service - pression.reste

    press.He		<- 0

    if (gg[2] != 0) {
        press.He		<- press.service * gg[2] / 100
        press.service	<- press.service - press.He
        gg[2]			<- 0
        gg[1]			<- 100 * gg[1] / (gg[1] + gg[3])
        gg[3]			<- 100 - gg[1]
    }

    if (gg[1] == 100) {
        press.O2		<- press.service
        press.air		<- 0
        press.nitrox	<- 0
    }

    if (gg[1] < 40 & gg[2] == 0 & stick == T) {
        press.O2		<- 0
        press.air		<- 0
        press.nitrox	<- press.service
    }

    if ((gg[1] >= 40 & gg[2] == 0) | (gg[1] < 40 & gg[2] == 0 & stick == F)) {
        N2 <- press.service * gg[3]/100
        press.air		<- N2 / 0.79
        press.O2		<- press.service - press.air
        press.nitrox	<- 0
    }

    Compo.Nitrox <- round(c(gg[1],gg[3]),1)
    names(Compo.Nitrox) <- c("O2","N2")
    press.gaz <- c(press.air,press.O2,press.nitrox,press.He,Compo.Nitrox)

    if (press.nitrox != 0) out <- list(O2=press.O2,Air=press.air,He=press.He,Nitrox=press.nitrox,Compo.nitrox=Compo.Nitrox)
    if (press.nitrox == 0) out <- list(O2=press.O2,Air=press.air,He=press.He,Nitrox=press.nitrox)
    if (any(press.gaz < 0) | Compo.Nitrox[2]>79) out <- paste("Erreur : melange impossible... Ajouter du N2 pur, ou commencer par vider le bloc")

    out
}

#
# BIGONE(temps=15,gaz=c(15,53))->ex
#
# matplot(xa,t(ex$PC),type="l",col=rainbow(16),lty=1,ylim=c(0,105),xlim=c(0,54),axes=F,xlab="Profondeur",ylab="Pourcentage de gradient")
# lines(c(0,54),c(75,20),lty=2)
# abline(h=100,lwd=3)
# abline(h=00,lwd=3)
# abline(v=00,lwd=3)
# abline(v=54,lwd=3)
# axis(side=1,at=seq(0,54,3))
# axis(side=2,at=seq(0,100,10))
# text(9,70,labels="Gradient",srt=-20)
#
#
# as.matrix(ex$PC)->PC
# PC[PC<0]<-0
# barplot(PC,beside=T,col=rainbow(16),names=c(90,90,rep(seq(51,3,-3),each=2),0))


# BIGONE(gaz=c(21,0),prof=20,temps=30,dur=3,grad=c(.85,.85))->ex
# ex
# nb.seg <- ncol(ex$pre)
# tens.He <- ex$pre[,seq(1,nb.seg,2)]
# tens.N2 <- ex$pre[,seq(2,nb.seg,2)]
# tension <- tens.He + tens.N2
# ex$M
# pc.M.val <- round(pre*100/ex$M,1)
#
# Pam <- matrix(c(40,40,30),ncol=3,nrow=16,byr=T)
# Pam
# (pre-Pam)*100/(ex$M-Pam)
# (Pam-pre)/(ex$M-Pam)
# (pre-Pam)/(ex$M-Pam)
# Pam <- matrix(c(30,30,10),ncol=3,nrow=16,byr=T)
# pc.grad <- (pre-Pam)*100/(ex$M-Pam)


# # ==============
# # = Graphiques =
# # ==============
# par(mfrow=c(1,2))
# plot(c(0,50),c(0,50),type="n",xlab="Pression ambiante, absolue",ylab="Pression absolue du gaz inerte dans le compartiment",main="M-Values Béhlman pour l'azote")
# for (i in 1:16) abline(as.matrix(ZHL16.C[i,c(3,4)]),col=rainbow(17)[i])
# abline(0,1,lwd=2)
# abline(v=1,lty=2)
# legend("bottomright",legend=ZHL16.C[,2],lty=1,col=rainbow(17),ncol=2,title="Période des 17 compartiments",cex=.7)
#
# plot(c(0,50),c(0,50),type="n",xlab="Pression ambiante, absolue",ylab="Pression absolue du gaz inerte dans le compartiment",main="M-Values Béhlman pour l'hélium")
# for (i in 1:16) abline(as.matrix(ZHL16.C[i,c(6,7)]),col=rainbow(17)[i])
# abline(0,1,lwd=2)
# abline(v=1,lty=2)
# legend("bottomright",legend=ZHL16.C[,5],lty=1,col=rainbow(17),ncol=2,title="Période des 17 compartiments",cex=.7)
# par(mfrow=c(1,1))
#
# plot(c(0,50),c(0,50),type="n",xlab="Pression ambiante, absolue",ylab="Pression absolue du gaz inerte dans le compartiment",main="M-Values Béhlman pour l'azote")
# axis(side=2,at=seq(0,50,5))
# for (i in 1:17) abline(as.matrix(coef[i,]),col=rainbow(17)[i])
# abline(0,1,lwd=2)
# abline(v=1,lty=2)
# legend("bottomright",legend=row.names(coef),lty=1,col=rainbow(17),ncol=2,title="Période des 17 compartiments",cex=.7)
#
# prof.max=98
# BIGONE(prof=prof.max,temps=18)->ex
# rep(ex[[2]],each=2)->pro
# pro<-c(prof.max,pro[1:(length(pro)-1)]) +10
# ex[[4]][,seq(1,dim(ex[[4]])[2],2)]->the1
# ex[[4]][,seq(2,dim(ex[[4]])[2],2)]->the2
# the <- the1+the2
# the <- the[,-1]
#
# matplot(pro,t(the),type="l",lty=1,col=rainbow(16),ylim=range(0,the),xlim=range(0,pro),xlab="Pression absolue",ylab="Tension de gaz inerte dans les tissus")
# abline(0,1)
# abline(v=0)
# abline(v=10,lty=2)
# coef <- ZHL16.C[,c(3,4)]
# # coef[,2] <- 1/coef[,2]
# for (i in 1:16) abline(as.matrix(coef[i,]),col=rainbow(16)[i])
# legend("bottomright",legend=c("1b",2:16),lty=1,col=rainbow(17),title="compartiments",cex=.7,ncol=8)


# 0.16 < PPO2 < 1.6
# sous le Siroco du séchoir é cheveux, la petite garce laisse choir 'je veux'
# Gérer le bug de la plongee sans palier
# Revoir le probléme des tensions tolerées de la fonction GF é remplacer par une fonction beaucoup plus simple "tension" qui fait simplement la somme des ppart des différents gaz inertes des chaque tissus et pour chaque segment de la plongée
# Voir comment prndre en compte les gradients de facteurs dans la condition while !





