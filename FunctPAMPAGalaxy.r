#Rscript


##################################################################################################################################
############## ?????????????????????????????????????????????????????????????????????????????????    ##############################
##################################################################################################################################

#### Based on Romain Lorrillière R script
#### Modified by Alan Amosse and Benjamin Yguel for integrating within Galaxy-E

##### workes with the R version 3.5.1 (2018-07-02)
##### Package used with the version:
#suppressMessages(library(lme4))  version 1.1.18.1
#suppressMessages(library(ggplot2))  version 3.0.0
#suppressMessages(library(speedglm))  version 0.3.2
#suppressMessages(library(arm))  version 1.10.1
#suppressMessages(library(reshape))  version 0.8.8
#suppressMessages(library(data.table))  version 1.12.0
#suppressMessages(library(reshape2))   version 1.4.3

######################################### start of the function fact.def.f called by FunctExeCalcCommIndexesGalaxy.r and FunctExeCalcPresAbsGalaxy.r
####### Define the finest aggregation with the observation table

fact.det.f <- function (Obs,
                        size.class="size.class",
                        code.especes="species.code",
                        unitobs="observation.unit")
{
    if (any(is.element(c(size.class), colnames(obs))) && all(! is.na(obs[, size.class])))
        {
            factors <- c(unitobs, code.especes, size.class)
        }else{
            factors <- c(unitobs, code.especes)
        }
    return(factors)
}

######################################### end of the function fact.def.f 

######################################### start of the function def.typeobs.f called by FunctExeCalcCommIndexesGalaxy.r and FunctExeCalcPresAbsGalaxy.r
####### Define observation type from colnames

def.typeobs.f <- function(Obs)
{
    if (any(is.element(c("rotation","rot","rotate"),colnames(obs))))
    {
        ObsType <- "SVR"
    }else{
        ObsType <- "other"
    }
    return(ObsType)
}
######################################### end of the function fact.def.f 

######################################### start of the function create.unitobs called by FunctExeCalcCommIndexesGalaxy.r and FunctExeCalcPresAbsGalaxy.r
####### Create unitobs column when inexistant
create.unitobs <- function(data,year="year",point="point", unitobs="observation.unit")
{
    if (!is.element("observation.unit",colnames(data)))
    {
        unitab <- unite(data,col="observation.unit",c(year,point))
    }else{ 
        unitab <- data
    }
    return(unitab)
}
######################################### start of the function create.unitobs

######################################### start of the function create.year.point called by FunctExeCalcCommIndexesGalaxy.r and FunctExeCalcPresAbsGalaxy.r
####### separate unitobs column when existant
create.year.point <- function(data,year="year",point="point", unitobs="observation.unit")
{
    if (all(grepl("[1-2][0|8|9][0-9]{2}_.*",data[,unitobs]))==TRUE)
    {
        tab <- separate(data,col=unitobs,into=c(year,point),sep="_")
        tab <- cbind(tab,data[,unitobs])
    }else{
        tab <- separate(data,col=unitobs,into=c("site1", year,"obs"),sep=c(2,4))
        tab <- unite(tab, col=point, c("site1","obs"))
        tab <- cbind(tab,data[,unitobs])
    }
    return(tab)
}
######################################### start of the function create.unitobs

######################################### start of the function check_file called by every Galaxy Rscripts

# Fonction pour verifier les données d'entrée / General function to check integrity of input file. Will check numbers and contents of variables(colnames). 
#return an error message and exit if mismatch detected
#Faut rentrer le nom du jeu de données, le nbre et le nom des variables / Enter dataset name,  expected number and names of variables. + an exit error message to guide user.

check_file<-function(dataset,err_msg,vars,nb_vars){
    if(ncol(dataset) < nb_vars){ #Verifiction de la présence du bon nb de colonnes, si c'est pas le cas= message d'erreur / checking for right number of columns in the file if not = error message
        cat("\nerr nb var\n") 
        stop(err_msg, call.=FALSE)
    }

    for(i in vars){
        if(!(i %in% names(dataset))){
            stop(err_msg,call.=FALSE)
        }
    }
}

######################################### end of the function check_file


######################################### start of the function statRotationsNumber.f called by calc.numbers.f

statRotationsNumber.f <- function(factors, obs)
{
    ## Purpose: Calcul des statistiques des abondances (max, sd) par rotation
    ##          en se basant sur des données déjà interpolées.
    ## ----------------------------------------------------------------------
    ## Arguments: factors : vecteur des noms de factors d'agrégation
    ##                       (résolution à laquelle on travaille).
    ##            obs : données d'observation.
    ##            dataEnv : environnement des données (pour sauvegarde de
    ##                      résultats intermédiaires).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 29 oct. 2012, 16:01

    ## Identification des rotations valides :
    if (is.element("observation.unit", factors))
    {
        ## Rotations valides (les vides doivent tout de même être renseignés) :
        rotations <- tapply(obs$rotation,
                            as.list(obs[ , c("observation.unit", "rotation"), drop=FALSE]),
                            function(x)length(x) > 0)

        ## Les rotations non renseignés apparaissent en NA et on veut FALSE :
        rotations[is.na(rotations)] <- FALSE
    }else{
        #stop(mltext("statRotations.err.1"))
    }

    ## ###########################################################
    ## Nombres par rotation avec le niveau d'agrégation souhaité :
    nombresR <- tapply(obs$number,
                       as.list(obs[ , c(factors, "rotation"), drop=FALSE]),
                       function(x,...){ifelse(all(is.na(x)), NA, sum(x,...))},
                       na.rm = TRUE)

    ## Les NAs sont considérés comme des vrais zéros lorsque la rotation est valide :
    nombresR <- sweep(nombresR,
                      match(names(dimnames(rotations)), names(dimnames(nombresR)), nomatch=NULL),
                      rotations,        # Tableau des secteurs valides (booléens).
                      function(x, y)
                  {
                      x[is.na(x) & y] <- 0 # Lorsque NA et secteur valide => 0.
                      return(x)
                  })

    ## ##################################################
    ## Statistiques :

    ## Moyennes :
    nombresMean <- apply(nombresR, which(is.element(names(dimnames(nombresR)), factors)),
                         function(x,...){ifelse(all(is.na(x)), NA, mean(x,...))}, na.rm=TRUE)

    ## Maxima :
    nombresMax <- apply(nombresR, which(is.element(names(dimnames(nombresR)), factors)),
                        function(x,...){ifelse(all(is.na(x)), NA, max(x,...))}, na.rm=TRUE)

    ## Déviation standard :
    nombresSD <- apply(nombresR, which(is.element(names(dimnames(nombresR)), factors)),
                       function(x,...){ifelse(all(is.na(x)), NA, sd(x,...))}, na.rm=TRUE)

    ## Nombre de rotations valides :
    nombresRotations <- apply(rotations, 1, sum, na.rm=TRUE)

    ## Retour des résultats sous forme de liste :
    return(list(nombresMean=nombresMean, nombresMax=nombresMax, nombresSD=nombresSD,
                nombresRotations=nombresRotations, nombresTot=nombresR))
}

######################################### end of the function statRotationsNumber.f 

######################################### start of the function calcNumber.default.f called by calc.numbers.f
## Calcul des nombres au niveau d'agrégation le plus fin.

calcNumber.default.f <- function(obs,
                                 factors=c("observation.unit", "species.code", "size.class"),
                                 nbName="number")
{
    ## Arguments: obs : table des observations (data.frame).
    ##            factors : les factors d'agrégation.
    ##            nbName : nom de la colonne nombre.
    ##
    ## Output: un array avec autant de dimensions que de factors
    ##         d'agrégation.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 19 déc. 2011, 13:38

    ## Somme des nombres d'individus :
    nbr <- tapply(obs[ , nbName],
                  as.list(obs[ , factors]),
                  sum, na.rm = TRUE)

    ## Absences considérée comme "vrais zéros" :
    nbr[is.na(nbr)] <- 0

    return(nbr)
}

######################################### end of the function calcNumber.default.f

######################################### start of the function calc.numbers.f

calc.numbers.f <- function(obs, ObsType="", factors=c("observation.unit", "species.code", "size.class"), nbName="number")
{
    ## Purpose: Produit la data.frame qui va servir de table, à partir du
    ##          tableau de nombres produit par calcNumber.default.f().
    ## ----------------------------------------------------------------------
    ## Arguments: nbr : array de nombres avec autant de dimensions que de
    ##                  factors d'agrégations.
    ##            nbName : nom de la colonne nombre.
    ##
    ## Output: une data.frame avec (nombre de factors d'agrégation + 1)
    ##         colonnes.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 19 déc. 2011, 13:46

    if (ObsType == "SVR")
    {
         ## Calcul des statistiques de nombres :
         statRotations <- statRotationsNumber.f(factors=factors,
                                                  obs=obs)

         ## Moyenne pour les vidéos rotatives (habituellement 3 rotation) :
         nbr <- statRotations[["nombresMean"]]

    }else{

         nbr <- calcNumber.default.f(obs, factors, nbName)
    }

    res <- as.data.frame(as.table(nbr), responseName=nbName)
    ## res$unitobs <- res$observation.unit # Pour compatibilité uniquement !!!

    if (is.element("size.class", colnames(res)))
    {
        res$size.class[res$size.class == ""] <- NA
    }else{}

    ## Si les nombres sont des entiers, leur redonner la bonne classe :
    if (isTRUE(all.equal(res[ , nbName], as.integer(res[ , nbName]))))
    {
        res[ , nbName] <- as.integer(res[ , nbName])
    }else{}

    if (ObsType == "SVR")
    {
        ## Statistiques sur les nombres :
        res$number.max <- as.vector(statRotations[["nombresMax"]])
        res$number.sd <- as.vector(statRotations[["nombresSD"]])
              
    }else{}

    return(res)
}

######################################### end of the function calc.numbers.f

######################################### start of the function presAbs.f called by calcBiodiv.f

presAbs.f <- function(nombres, logical=FALSE)
{
    ## Purpose: Renvoie les présences/absences d'après les nombres.
    ## ----------------------------------------------------------------------
    ## Arguments: nombres : vecteur de nombre d'individus.
    ##            logical : faut-il renvoyer les résultats sous forme de
    ##                      booléens, ou 0/1 (booléen).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 29 oct. 2010, 10:20

    if (any(nombres < 0, na.rm=TRUE))
    {
        stop("Negative abundances!")
    }else{}

    if (logical)
    {
        return(nombres > 0)
    }else{
        nombres[nombres > 0] <- 1
        return(nombres)
    }
}

######################################### end of the function presAbs.f

######################################### start of the function betterCbind called by agregations.generic.f

betterCbind <- function(..., dfList=NULL, deparse.level = 1)
{
    ## Purpose: Appliquer un cbind à des data.frames qui ont des colonnes
    ##          communes en supprimant les redondances (comme un merge mais
    ##          les lignes doivent être en mêmes nombres et
    ##          dans le même ordre)
    ## ----------------------------------------------------------------------
    ## Arguments: ceux de cbind...
    ##            dfList : une liste de data.frames (evite un do.call
    ##                     supplémentaire).
    ##                     ... est utilisé à la place si NULL.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 17 janv. 2012, 21:10

    if (is.null(dfList))
    {
        dfList <- list(...)
    }else{}

    return(do.call(cbind,
                   c(list(dfList[[1]][ , c(tail(colnames(dfList[[1]]), -1),
                                           head(colnames(dfList[[1]]), 1))]),
                     lapply(dfList[-1],
                            function(x, colDel)
                        {
                            return(x[ , !is.element(colnames(x),
                                                    colDel),
                                     drop=FALSE])
                        },
                            colDel=colnames(dfList[[1]])),
                     deparse.level=deparse.level)))
}

######################################### end of the function betterCbind

######################################### start of the function agregation.f called by agregations.generic.f

agregation.f <- function(metric, Data, factors, casMetrique, dataEnv,
                         nbName="number")
{
    ## Purpose: Agrégation d'une métrique.
    ## ----------------------------------------------------------------------
    ## Arguments: metric: nom de la colonne de la métrique.
    ##            Data: table de données non-agrégées.
    ##            factors: vecteur des noms de factors d'agrégation.
    ##            casMetrique: vecteur nommé des types d'observation en
    ##                         fonction de la métrique choisie.
    ##            dataEnv: environnement des données.
    ##            nbName : nom de la colonne nombre.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 20 déc. 2011, 14:29

    switch(casMetrique[metric],
           "sum"={
               res <- tapply(Data[ , metric],
                             as.list(Data[ , factors, drop=FALSE]),
                             function(x)
                         {
                             ifelse(all(is.na(x)),
                                    NA,
                                    sum(x, na.rm=TRUE))
                         })
           },
           "w.mean"={
               res <- tapply(1:nrow(Data),
                             as.list(Data[ , factors, drop=FALSE]),
                             function(ii)
                         {
                             ifelse(all(is.na(Data[ii, metric])),
                                    NA,
                                    weighted.mean(Data[ii, metric],
                                                  Data[ii, nbName],
                                                  na.rm=TRUE))
                         })
           },
           "w.mean.colonies"={
               res <- tapply(1:nrow(Data),
                             as.list(Data[ , factors, drop=FALSE]),
                             function(ii)
                         {
                             ifelse(all(is.na(Data[ii, metric])),
                                    NA,
                                    weighted.mean(Data[ii, metric],
                                                  Data[ii, "colonies"],
                                                  na.rm=TRUE))
                         })
           },
           "w.mean.prop"={
               res <- tapply(1:nrow(Data),
                             as.list(Data[ , factors, drop=FALSE]),
                             function(ii)
                         {
                             ifelse(all(is.na(Data[ii, metric])) || sum(Data[ii, "nombre.tot"], na.rm=TRUE) == 0,
                                    NA,
                                    ifelse(all(na.omit(Data[ii, metric]) == 0), # Pour ne pas avoir NaN.
                                           0,
                                           (sum(Data[ii, nbName][ !is.na(Data[ii, metric])], na.rm=TRUE) /
                                            sum(Data[ii, "nombre.tot"], na.rm=TRUE)) *
                                           ## Correction si la classe de taille n'est pas un facteur d'agrégation
                                           ## (sinon valeur divisée par le nombre de classes présentes) :
                                           ifelse(is.element("size.class", factors),
                                                  100,
                                                  100 * length(unique(Data$size.class)))))
                         })

           },
           "w.mean.prop.bio"={
               res <- tapply(1:nrow(Data),
                             as.list(Data[ , factors, drop=FALSE]),
                             function(ii)
                         {
                             ifelse(all(is.na(Data[ii, metric])) || sum(Data[ii, "tot.biomass"], na.rm=TRUE) == 0,
                                    NA,
                                    ifelse(all(na.omit(Data[ii, metric]) == 0), # Pour ne pas avoir NaN.
                                           0,
                                           (sum(Data[ii, "biomass"][ !is.na(Data[ii, metric])], na.rm=TRUE) /
                                            sum(Data[ii, "tot.biomass"], na.rm=TRUE)) *
                                           ## Correction si la classe de taille n'est pas un facteur d'agrégation
                                           ## (sinon valeur divisée par le nombre de classes présentes) :
                                           ifelse(is.element("size.class", factors),
                                                  100,
                                                  100 * length(unique(Data$size.class)))))
                         })

           },
           "pres"={
               res <- tapply(Data[ , metric],
                             as.list(Data[ , factors, drop=FALSE]),
                             function(x)
                         {
                             ifelse(all(is.na(x)), # Cas où il n'y a que des NAs.
                                    NA,
                                    ifelse(any(x > 0, na.rm=TRUE), # Sinon...
                                           1, # ...présence si au moins une observation dans le groupe.
                                           0))
                         })
           },
           "nbMax"={
               ## Récupération des nombres brutes avec sélections :
               nbTmp <- getReducedSVRdata.f(dataName=".NombresSVR", data=Data, dataEnv=dataEnv)

               ## Somme par croisement de facteur / rotation :
               nbTmp2 <- apply(nbTmp,
                             which(is.element(names(dimnames(nbTmp)), c(factors, "rotation"))),
                             function(x)
                         {
                             ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))
                         })

               ## Somme par croisement de facteur :
               res <- as.array(apply(nbTmp2,
                                     which(is.element(names(dimnames(nbTmp)), factors)),
                                     function(x)
                                 {
                                     ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE))
                                 }))
           },
           "nbSD"={
               ## Récupération des nombres brutes avec sélections :
               nbTmp <- getReducedSVRdata.f(dataName=".NombresSVR", data=Data, dataEnv=dataEnv)

               ## Somme par croisement de facteur / rotation :
               nbTmp2 <- apply(nbTmp,
                             which(is.element(names(dimnames(nbTmp)), c(factors, "rotation"))),
                             function(x)
                         {
                             ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))
                         })

               ## Somme par croisement de facteur :
               res <- as.array(apply(nbTmp2,
                                     which(is.element(names(dimnames(nbTmp)), factors)),
                                     function(x)
                                 {
                                     ifelse(all(is.na(x)), NA, sd(x, na.rm=TRUE))
                                 }))
           },
           "densMax"={
               ## Récupération des nombres brutes avec sélections :
               densTmp <- getReducedSVRdata.f(dataName=".DensitesSVR", data=Data, dataEnv=dataEnv)

               ## Somme par croisement de facteur / rotation :
               densTmp2 <- apply(densTmp,
                                 which(is.element(names(dimnames(densTmp)), c(factors, "rotation"))),
                                 function(x)
                             {
                                 ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))
                             })

               ## Somme par croisement de facteur :
               res <- as.array(apply(densTmp2,
                                     which(is.element(names(dimnames(densTmp)), factors)),
                                     function(x)
                                 {
                                     ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE))
                                 }))
           },
           "densSD"={
               ## Récupération des nombres brutes avec sélections :
               densTmp <- getReducedSVRdata.f(dataName=".DensitesSVR", data=Data, dataEnv=dataEnv)

               ## Somme par croisement de facteur / rotation :
               densTmp2 <- apply(densTmp,
                                 which(is.element(names(dimnames(densTmp)), c(factors, "rotation"))),
                                 function(x)
                             {
                                 ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))
                             })

               ## Somme par croisement de facteur :
               res <- as.array(apply(densTmp2,
                                     which(is.element(names(dimnames(densTmp)), factors)),
                                     function(x)
                                 {
                                     ifelse(all(is.na(x)), NA, sd(x, na.rm=TRUE))
                                 }))
           },
           "%.nesting"={
               res <- tapply(1:nrow(Data),
                             as.list(Data[ , factors, drop=FALSE]),
                             function(ii)
                         {
                             ifelse(all(is.na(Data[ii, metric])),
                                    NA,
                                    weighted.mean(Data[ii, metric],
                                                  Data[ii, "readable.tracks"],
                                                  na.rm=TRUE))
                         })
           },
           stop("Not implemented!")
           )

    ## Nom des dimensions
    names(dimnames(res)) <- c(factors)

    ## Transformation vers format long :
    reslong <- as.data.frame(as.table(res), responseName=metric)
    reslong <- reslong[ , c(tail(colnames(reslong), 1), head(colnames(reslong), -1))] # métrique en première.

    return(reslong)
}

######################################### end of the function agregation.f

######################################### start of the function agregations.generic.f called y calcBiodiv.f in FucntExeCalcCommIndexesGalaxy.r

agregations.generic.f <- function(Data, metrics, factors, listFact=NULL, unitSpSz=NULL, unitSp=NULL, info=FALSE,
                                  dataEnv=.GlobalEnv, nbName="number")
{
    ## Purpose: Agréger les données selon un ou plusieurs factors.
    ## ----------------------------------------------------------------------
    ## Arguments: Data : Le jeu de données à agréger.
    ##            metrics : la métrique agrégée.
    ##            factors : les factors
    ##            listFact : noms des factors supplémentaires (agrégés et
    ##                       ajoutés à la table de sortie).
    ##            unitSpSz : Table de métriques par unitobs/esp/CT.
    ##            unitSp : Table de métriques par unitobs/esp
    ##            info : affichage des infos ?
    ##            nbName : nom de la colonne nombre.
    ##
    ## Output : une data.frame agrégée.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 18 oct. 2010, 15:47

    ## Informations (l'étape peut être longue) :
    #if (info) WinInfo <- agregation.info.f()

    ## traitements selon le type de métrique :
    casMetrique <- c("number"="sum",
                     "mean.length"="w.mean",
                     "taille_moy"="w.mean",
                     "biomass"="sum",
                     "Biomass"="sum",
                     "weight"="sum",
                     "mean.weight"="w.mean",
                     "density"="sum",
                     "Density"="sum",
                     "CPUE"="sum",
                     "CPUE.biomass"="sum",
                     "pres.abs"="pres",
                     "abundance.prop.SC"="w.mean.prop", # Pas bon [!!!] ?
                     "biomass.prop.SC"="w.mean.prop.bio",  # Pas bon [!!!] ?
                     ## Benthos :
                     "colonies"="sum",
                     "coverage"="sum",
                     "mean.size.colonies"="w.mean.colonies",
                     ## SVR (expérimental) :
                     "number.max"="nbMax",
                     "number.sd"="nbSD",
                     "density.max"="densMax",
                     "density.sd"="densSD",
                     "biomass.max"="sum",
                     "spawning.success"="%.nesting",
                     "spawnings"="sum",
                     "readable.tracks"="sum",
                     "tracks.number"="sum")

    ## Ajout de "readable.tracks" pour le pourcentage de ponte :
    if (any(casMetrique[metrics] == "%.nesting"))
    {
        if (is.element("size.class", colnames(Data)))
        {
            if (is.null(unitSpSz)) stop("unitSpSz doit être défini")

            Data <- merge(Data,
                          unitSpSz[ , c("species.code", "observation.unit", "size.class", "readable.tracks")],
                          by=c("species.code", "observation.unit", "size.class"),
                          suffixes=c("", ".y"))
        }else{
            if (is.null(unitSp)) stop("unitSp must be defined")

            Data <- merge(Data,
                          unitSp[ , c("species.code", "observation.unit", "readable.tracks")],
                          by=c("species.code", "observation.unit"),
                          suffixes=c("", ".y"))
        }
    }else{}

    ## Ajout du champ nombre pour le calcul des moyennes pondérées s'il est absent :
    if (any(casMetrique[metrics] == "w.mean" | casMetrique[metrics] == "w.mean.prop"))
    {
        if (is.element("size.class", colnames(Data)))
        {
            if (is.null(unitSpSz)) stop("unitSpSz must be defined")

            Data <- merge(Data,
                          unitSpSz[ , c("species.code", "observation.unit", "size.class", nbName)],
                          by=c("species.code", "observation.unit", "size.class"),
                          suffixes=c("", ".y"))

            ## Ajout de l'abondance totale /espèce/unité d'observation :
            nbTot <- tapply(unitSpSz[ , nbName],
                            as.list(unitSpSz[ , c("species.code", "observation.unit")]),
                            sum, na.rm=TRUE)

            Data <- merge(Data,
                          as.data.frame(as.table(nbTot), responseName="nombre.tot"))
        }else{
            if (is.null(unitSp)) stop("unitSp must be defined")

            Data <- merge(Data,
                          unitSp[ , c("species.code", "observation.unit", nbName)], # [!!!] unitSpSz ?
                          by=c("species.code", "observation.unit"),
                          suffixes=c("", ".y"))
        }
    }else{}

    ## Ajout du champ biomasse pour les proportions de biomasses par classe de taille :
    if (any(casMetrique[metrics] == "w.mean.prop.bio"))
    {
        if (is.null(unitSpSz)) stop("unitSpSz doit être défini")

        Data <- merge(Data,
                      unitSpSz[ , c("species.code", "observation.unit", "size.class", "biomass")],
                      by=c("species.code", "observation.unit", "size.class"),
                      suffixes=c("", ".y"))

        ## Ajout de la biomasse totale /espèce/unité d'observation :
        biomTot <- tapply(unitSpSz$biomass,
                          as.list(unitSpSz[ , c("species.code", "observation.unit")]),
                          function(x)
                      {
                          ifelse(all(is.na(x)),
                                 NA,
                                 sum(x, na.rm=TRUE))
                      })

        Data <- merge(Data,
                      as.data.frame(as.table(biomTot), responseName="tot.biomass"))
    }

    ## Ajout du champ colonie pour le calcul des moyennes pondérées s'il est absent :
    if (any(casMetrique[metrics] == "w.mean.colonies" & ! is.element("colonies", colnames(Data))))
    {
        Data$colonies <- unitSp[match(apply(Data[ , c("species.code", "observation.unit")],
                                           1, paste, collapse="*"),
                                     apply(unitSp[ , c("species.code", "observation.unit")],
                                           1, paste, collapse="*")), "colonies"]
    }else{}


    ## Agrégation de la métrique selon les factors :
    reslong <- betterCbind(dfList=lapply(metrics,   # sapply utilisé pour avoir les noms.
                                         agregation.f,
                                         Data=Data, factors=factors, casMetrique=casMetrique, dataEnv=dataEnv,
                                         nbName=nbName))

    ## Agrégation et ajout des factors supplémentaires :
    if ( ! (is.null(listFact) || length(listFact) == 0))
    {
        reslong <- cbind(reslong,
                         sapply(Data[ , listFact, drop=FALSE],
                                function(fact)
                            {
                                tapply(fact,
                                       as.list(Data[ , factors, drop=FALSE]),
                                       function(x)
                                   {
                                       if (length(x) > 1 && length(unique(x)) > 1) # On doit n'avoir qu'une seule
                                        # modalité...
                                       {
                                           return(NULL)                  # ...sinon on retourne NULL
                                       }else{
                                           unique(as.character(x))
                                       }
                                   })
                            }))
    }else{}

    ## Si certains factors ne sont pas de classe facteur, il faut les remettre dans leur classe d'origine :
    if (any(tmp <- sapply(reslong[ , listFact, drop=FALSE], class) != sapply(Data[ , listFact, drop=FALSE], class)))
    {
        for (i in which(tmp))
        {
            switch(sapply(Data[ , listFact, drop=FALSE], class)[i],
                   "integer"={
                       reslong[ , listFact[i]] <- as.integer(as.character(reslong[ , listFact[i]]))
                   },
                   "numeric"={
                       reslong[ , listFact[i]] <- as.numeric(as.character(reslong[ , listFact[i]]))
                   },
                   reslong[ , listFact[i]] <- eval(call(paste("as", sapply(Data[ , listFact, drop=FALSE], class)[i], sep="."),
                                                        reslong[ , listFact[i]]))
                   )
        }
    }else{}

    ## Rétablir l'ordre initial des nivaux de factors :
    reslong <- as.data.frame(sapply(colnames(reslong),
                                    function(x)
                                {
                                    if (is.factor(reslong[ , x]))
                                    {
                                        return(factor(reslong[ , x], levels=levels(Data[ , x])))
                                    }else{
                                        return(reslong[ , x])
                                    }
                                }, simplify=FALSE))


    ## Fermeture de la fenêtre d'information
    if (info) close.info.f(WinInfo)

    ## Vérification des factors supplémentaires agrégés. Il ne doit pas y avoir d'élément nul (la fonction précédente
    ## renvoie NULL si plusieurs niveaux de factors, i.e. le facteur est un sous ensemble d'un des factors
    ## d'agrégation des observations) :
    if (any(sapply(reslong[ , listFact], function(x){any(is.null(unlist(x)))})))
    {
        warning(paste("One of the suppl. factors is probably a subset",
                      " of the observations grouping factor(s).", sep=""))
        return(NULL)
    }else{
        return(reslong)
    }
}

######################################### end of the function agregations.generic.f

######################################### start of the function dropLevels.f called y calcBiodiv.f in FucntExeCalcCommIndexesGalaxy.r and modeleLineaireWP2.unitobs.f in FunctExeCalcGLMGalaxy.r
dropLevels.f <- function(df, which=NULL)
{
    ## Purpose: Supprimer les 'levels' non utilisés des factors d'une
    ##          data.frame.
    ## ----------------------------------------------------------------------
    ## Arguments: df : une data.frame
    ##            which : indice des colonnes à inclure (toutes par défaut).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 10 août 2010, 13:29

    if (class(df) != "data.frame")
    {
        stop("'df' must be a data.frame")
    }else{
        if (is.null(which))
        {
            x <- as.data.frame(sapply(df, function(x)
                                  {
                                      return(x[ ,drop=TRUE])
                                  }, simplify=FALSE),
                               stringsAsFactors=FALSE)
        }else{                          # Cas où seulement certaines colonnes sont traitées.
            x <- df

            x[ , which] <- as.data.frame(sapply(df[ , which, drop=FALSE],
                                                function(x)
                                            {
                                                return(x[ , drop=TRUE])
                                            }, simplify=FALSE),
                                         stringsAsFactors=FALSE)
        }

        return(x)
    }
}
######################################### end of the function dropLevels.f

######################################### start of the function subsetToutesTables.f called by modeleLineaireWP2.unitobs.f in FunctExeCalcGLMGalaxy.r

subsetToutesTables.f <- function(metrique, tabMetrics, facteurs, selections,
                                 tabUnitobs, refesp, tableMetrique="", nbName="number", ObsType = "",
                                 exclude=NULL, add=c("species.code", "observation.unit"))
{
    ## Purpose: Extraire les données utiles uniquement, d'après les métrique
    ##          et facteur(s) séléctionnés, ainsi que leur(s) sélection(s) de
    ##          modalité(s).
    ## ----------------------------------------------------------------------
    ## Arguments: metrique : la métrique choisie.
    ##            facteurs : les facteurs sélectionnés (tous)
    ##            selections : les sélections de modalités correspondantes
    ##                         (liste).
    ##            tableMetrique : le nom de la table des métriques.
    ##            exclude : niveau de facteur à ne pas prendre en compte pour
    ##                      le subset.
    ##            add : champ(s) (de la table de métrique) à ajouter aux
    ##                  données.
    ##            dataEnv : l'environnement des données.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  6 août 2010, 16:46

    ## Si pas de table de métrique disponible ou déjà calculée
    ## ("TableOcurrences" est calculée à partir de la sélection) :
    if (is.element(tableMetrique, c("", "TableOccurrences", "TablePresAbs")))
    {
        tableMetrique <- "unitSp"
    }else{}

    casTables <- c("unitSp"="unitSp",
                   "TablePresAbs"="unitSp",
                   "unitSpSz"="unitSpSz")

    ## Récupération de la table de métriques :
    dataMetrique <- tabMetrics
    unitobs <- tabUnitobs
    refesp <- refesp

    ## Si pas de métrique disponible ou déjà calculée ("freq.occurrence" est calculée à partir de la sélection) :
    if (is.element(metrique, c("", "occurrence.frequency")))
    {
        metrique <- "tmp"
        dataMetrique$tmp <- 0
        dataMetrique$tmp[dataMetrique[ , nbName] > 0] <- 1
    }else{}

    if (!is.null(add))
    {
        metriques <- c(metrique, add[is.element(add, colnames(dataMetrique))])
    }else{
        metriques <- metrique
    }

    ## Subset en fonction de la table de métrique
    switch(casTables[tableMetrique],
           ## Cas de la table d'observation ou des tables de présence :
           unitSp={
                restmp <- cbind(dataMetrique[!is.na(dataMetrique[ , metrique]) , metriques, drop=FALSE],
                                unitobs[match(dataMetrique$observation.unit[!is.na(dataMetrique[ , metrique])],
                                              unitobs$observation.unit), # ajout des colonnes sélectionnées d'unitobs
                                        facteurs[is.element(facteurs, colnames(unitobs))], drop=FALSE],
                                refesp[match(dataMetrique$species.code[!is.na(dataMetrique[ , metrique])],
                                             refesp$species.code),        # ajout des colonnes sélectionnées d'especes
                                       facteurs[is.element(facteurs, colnames(refesp))], drop=FALSE])
            },
           ## Cas de la table d'observations par classes de taille :
           unitSpSz={
               restmp <- cbind(dataMetrique[!is.na(dataMetrique[ , metrique]) ,
                                            c(metriques, "size.class"), drop=FALSE],
                               unitobs[match(dataMetrique$observation.unit[!is.na(dataMetrique[ , metrique])],
                                             unitobs$observation.unit), # ajout des colonnes sélectionnées d'unitobs
                                       facteurs[is.element(facteurs, colnames(unitobs))], drop=FALSE],
                               refesp[match(dataMetrique$species.code[!is.na(dataMetrique[ , metrique])],
                                            refesp$species.code),        # ajout des colonnes sélectionnées d'especes
                                      facteurs[is.element(facteurs, colnames(refesp))], drop=FALSE])
           },
           ## Autres cas :
           restmp <- cbind(dataMetrique[!is.na(dataMetrique[ , metrique]) , metriques, drop=FALSE],
                           unitobs[match(dataMetrique$observation.unit[!is.na(dataMetrique[ , metrique])],
                                         unitobs$observation.unit), # ajout des colonnes sélectionnées d'unitobs.
                                   facteurs[is.element(facteurs, colnames(unitobs))], drop=FALSE])
           )

    selCol <- which(!is.na(selections))
    if (!is.null(exclude))
    {
        selCol <- selCol[selCol != exclude]
    }

    ####Subset des modalités de facteurs a dégager mais fonctionne pas à corriger !!!!
    #for (i in selCol)
    #{
    #    restmp <- subset(restmp, is.element(restmp[ , facteurs[i]], selections[[i]]))
    #}

    ## Traitement particulier des classes de taille (mise en facteur avec ordre défini selon le context) :
    if (is.element("size.class", colnames(restmp)))
    {
        if (length(grep("^[[:digit:]]*[-_][[:digit:]]*$", unique(as.character(restmp$size.class)), perl=TRUE)) ==
            length(unique(as.character(restmp$size.class))))
        {
            restmp$size.class <-
                factor(as.character(restmp$size.class),
                       levels=unique(as.character(restmp$size.class))[
                               order(as.numeric(sub("^([[:digit:]]*)[-_][[:digit:]]*$",
                                                    "\\1",
                                                    unique(as.character(restmp$size.class)),
                                                    perl=TRUE)),
                                     na.last=FALSE)])
        }else{
            restmp$size.class <- factor(restmp$size.class)
        }
    }else{}

    ## Conversion des biomasses et densités -> /100m² :
    if (any(is.element(colnames(restmp), c("biomass", "density",
                                           "biomass.max", "density.max",
                                           "biomass.sd", "density.sd"))) && ObsType != "fishing")
    {
        restmp[ , is.element(colnames(restmp),
                             c("biomass", "density",
                               "biomass.max", "density.max",
                               "biomass.sd", "density.sd"))] <- 100 *
                                   restmp[, is.element(colnames(restmp),
                                                       c("biomass", "density",
                                                         "biomass.max", "density.max",
                                                         "biomass.sd", "density.sd"))]
    }else{}

    return(restmp)
}

######################################### end of the function subsetToutesTables.f

######################################### start of the function calcLM.f called by modeleLineaireWP2.unitobs.f in FunctExeCalcGLMGalaxy.r
calcLM.f <- function(loiChoisie, formule, metrique, Data)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 17 sept. 2010, 14:49

     switch(loiChoisie,
            ## Modèle linéaire :
            NO={
                res <- lm(formule, data=Data)
            },
            ## Modèle linéaire, données log-transformées :
            LOGNO={
                ## Ajout d'une constante à la métrique si contient des zéros :
                if (sum(Data[ , metrique] == 0, na.rm=TRUE))
                {
                    Data[ , metrique] <- Data[ , metrique] +
                        ((min(Data[ , metrique], na.rm=TRUE) + 1) / 1000)
                }else{}

                res <- lm(formule, data=Data)
            },
            ## GLM, distribution Gamma :
            GA={
                ## Ajout d'une constante à la métrique si contient des zéros :
                if (sum(Data[ , metrique] == 0, na.rm=TRUE))
                {
                    Data[ , metrique] <- Data[ , metrique] +
                        ((min(Data[ , metrique], na.rm=TRUE) + 1) / 1000)
                }else{}

                res <- glm(formule, data=Data, family="Gamma")
            },
            ## GLM, distribution de Poisson :
            PO={
                res <- glm(formule, data=Data, family="poisson")
            },
            ## GLM, distribution binomiale négative :
            NBI={
                res <- glm.nb(formule, data=Data)
            },
            ## GLM, distribution binomiale (présences/absences) :
            BI={
                res <- glm(formule, data=Data, family="binomial")
            },
            )
     return(res)
}



######################################### end of the function calcLM.f

######################################### start of the function sortiesLM.f called by modeleLineaireWP2.unitobs.f in FunctExeCalcGLMGalaxy.r
sortiesLM.f <- function(objLM, formule, metrique, factAna, modSel, listFact, listFactSel, Data, dataEnv,
                        Log=FALSE, sufixe=NULL, type="espece", baseEnv=.GlobalEnv)
{
    ## Purpose: Formater les résultats de lm et les écrire dans un fichier
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : un objet de classe lm
    ##            formule : la formule utilisée (pas lisible dans le call).
    ##            metrique : la métrique choisie.
    ##            factAna : le facteur de séparation des analyses.
    ##            modSel : la modalité courante.
    ##            listFact : liste du (des) facteur(s) de regroupement.
    ##            Data : les données utilisées.
    ##            Log : données log-transformées ou non (booléen).
    ##            sufixe : un sufixe pour le nom de fichier.
    ##            type : type d'analyse, pour traitement conditionnel des
    ##                   titres et noms de fichiers.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 25 août 2010, 16:19


    ## Ajout d'une constante si des zéros dans la métrique + transformation 'log' :
    if (sum(Data[ , metrique] == 0, na.rm=TRUE) & Log)
    {
        Data[ , metrique] <- Data[ , metrique] +
            ((min(Data[ , metrique], na.rm=TRUE) + 1) / 1000)
    }else{}

    ## Formule de modèle lisible:
    objLM$call$formula <- formule
    formule <<- formule
    resLM <<- objLM

    ## Chemin et nom de fichier :
    #resFile <- resFileLM.f(objLM=objLM, metrique=metrique, factAna=factAna, modSel=modSel, listFact=listFact,
     #                      dataEnv=dataEnv, Log=Log, sufixe=sufixe, type=type)
    #on.exit(tryCatch(close(resFile), error=function(e){}), add=TRUE)

    resFile <- "GLMSummary.txt"
    ## Informations et statistiques globales sur le modèle :
    infoStatLM.f(objLM=objLM, resFile=resFile)


    ## Anova globale du modèle + significativité des coefficients :
    signifParamLM.f(objLM=objLM, metrique=metrique, listFact=listFact, resFile=resFile)


    ## ##################################################
    ## Valeurs prédites par le modèle :
    valPreditesLM.f(objLM=objLM, Data=Data, listFact=listFact, resFile=resFile)

    ## ##################################################
    ## Comparaisons multiples :

    ## if (all(is.element(c("year", "protection.status"), listFact)))
    if (length(listFact) == 2)
    {
 

        ## compMultiplesLM.f(objLM=objLM, Data=Data, factSpatial="protection.status", factTemp="year", resFile=resFile)
        compMultiplesLM.f(objLM=objLM, Data=Data, fact1=listFact[1], fact2=listFact[2],
                          resFile=resFile,Log=Log)

        ## Représentation des interactions :suppr
        
    }else{
        if (length(listFact) == 1)
        {
            compSimplesLM.f(objLM=objLM, Data=Data, fact=listFact,
                            resFile=resFile, Log=Log)
        }else{}
    }

    # suppr

    ## ##################################################
    ## Sauvegarde des données :
    filename <- "GLMSummaryFull.txt"

 #   close(resFile)                      # Maintenant seulement on peut fermer ce fichier.

#    if ( ! isTRUE(sufixe == "(red)"))
 #   {
  #      writeData.f(filename=filename, Data=Data,
   #                 cols=NULL)
    #}else{}

    ## Sauvegarde des infos sur les données et statistiques :
    if ( ! isTRUE(sufixe == "(red)"))
    {
        infoStats.f(filename=filename, Data=Data, agregLevel=type, type="stat",
                    metrique=metrique, #factGraph=factAna, factGraphSel=modSel,
                    listFact=listFact)#, listFactSel=listFactSel,
                    #dataEnv=dataEnv, baseEnv=baseEnv)
    }else{}

    ## flush.console()
}


######################################### end of the function sortiesLM.f

######################################### start of the function infoStatLM.f called by sortiesLM.f

infoStatLM.f <- function(objLM, resFile)
{
    ## Purpose: Écrit les informations sur le modèle insi que les
    ##          statistiques globale dans un fichier résultat
    ## ----------------------------------------------------------------------
    ## Arguments: objLM un objet de classe 'lm' ou 'glm'.
    ##            resFile : une connection pour les sorties.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  8 sept. 2010, 16:57

    ## [!!!] Attention, il arrive que les calculs bloquent ici lors du premier lancement (origine inconnue)
    sumLM <- switch(class(objLM)[1],
                    lm = summary.lm(objLM),
                    glm = summary.glm(objLM),
                    negbin = MASS:::summary.negbin(objLM),
                    summary(objLM))

    ## Informations sur le modèle :
    cat("Fitted model:", file=resFile, fill=1,append=TRUE)
    cat("\t", deparse(objLM$call), "\n\n\n", file=resFile, sep="",append=TRUE)

    ## Stats globales :
    if (length(grep("^glm", objLM$call)) == 0)
    {
        cat("Global Fisher's statistics and R^2:", "\n\n", file=resFile,append=TRUE)
        cat("\t", "Multiple R^2: ", format(sumLM$r.squared, digits=3),
            " ",
            "\t", "Adjusted R^2 ", format(sumLM$adj.r.squared, digits=3), "\n", file=resFile, sep="",append=TRUE)

        cat("\t", "F-statistics:",
            paste(sapply(sumLM$fstatistic, format, digits=4, nsmall=0),
                  " over and DF,"),
            "\t", "P-value: ",
            format.pval(pf(sumLM$fstatistic[1L], sumLM$fstatistic[2L], sumLM$fstatistic[3L], lower.tail = FALSE),
                        digits=4),
            "\n\n\n", file=resFile, sep="",append=TRUE)
    }else{}
}

######################################### end of the function infoStatLM.f

######################################### start of the function signifParamLM.f called by sortiesLM.f

signifParamLM.f <- function(objLM, metrique, listFact, resFile)
{
    ## Purpose: Écrire les résultats de l'anova globale du modèle et
    ##          l'estimation de significativités des coefficients du modèle.
    ## ----------------------------------------------------------------------
    ## Arguments: objLM un objet de classe 'lm' ou 'glm'.compM
    ##            resFile : une connection pour les sorties.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  8 sept. 2010, 17:07

    ## Anovas et résumés :
    if (length(grep("^glmmTMB", objLM$call)) > 0) #GLMMTMB
    {

        anovaLM <- NULL
        sumLM <- summary(objLM)
        var <- as.numeric(sumLM$varcor$cond)
        varstd <- as.data.frame(capture.output(sumLM$varcor)[-c(1:2)])
        varstd$var <- c("Variance",var)
        colnames(varstd) <- rep("",length(colnames(varstd)))

        cat("---------------------------------------------------------------------------", 
            "\nSummary table :",
            "\n\nFamily : ", sumLM$family, ", link : ", sumLM$link,
            "\nResponse : ", metrique,
            "\n\nAIC table : \n\n AIC\tBIC\tlogLik\tdeviance\tdf.resid\n", sumLM$AICtab,
            file=resFile,append=TRUE)
        cat("\n\nAnalysis of Variance table (random effects) :\n",file=resFile,append=TRUE)
        capture.output(varstd,file=resFile,append=TRUE)     
        cat("\nNumber of obs : ", sumLM$nobs,", groups : ",sumLM$ngrps$cond,
            file=resFile,append=TRUE)
        ## Significativités des paramètres :
        cat("\n\n", "Parameter significances (fixed effects) :", #"\n(only the significant factors/interactions are shown):",
            "\n\n",
            file=resFile,append=TRUE)

        capture.output(printCoefmat.red(sumLM$coef$cond, anovaLM=anovaLM, objLM=objLM), file=resFile,append=TRUE)  
    }else{
        if (length(grep("^glm", objLM$call)) > 0) # Pour les GLMs.
        {
            anovaLM <- anova(objLM, test="Chisq") 
       
        }else{
            anovaLM <- anova(objLM) # Pour les LMs.
        }

    ## Anova globale du modèle :
    capture.output(print.anova.ml(anovaLM), file=resFile,append=TRUE)
    

    sumLM <- summary(objLM)
    ## Significativités des paramètres :
    cat("\n\n", "Parameter significances :", #"\n(only the significant factors/interactions are shown):",
        "\n\n",
        file=resFile,append=TRUE)

    capture.output(printCoefmat.red(sumLM$coef, anovaLM=anovaLM, objLM=objLM), file=resFile,append=TRUE)
    }  
}

######################################### end of the function signifParamLM.f 

######################################### start of the function print.anova.ml called by signifParamLM.f 

print.anova.ml <- function(x, digits = max(getOption("digits") - 2, 3), signif.stars = getOption("show.signif.stars"),
                           ...)
{
    ## Purpose: Hack de la méthode print.anova pour (franciser les sorties et)
    ##          supprimer les infos inutiles.
    ## ----------------------------------------------------------------------
    ## Arguments: ceux de print.anova
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 26 août 2010, 11:36

    #attr(x, "row.names")[attr(x, "row.names") == "Residuals"] <- mltext("print.anova.ml.KW.resid")

    ## Françisation des en-têtes (gsub itératif) :
    #attr(x, "heading") <- iter.gsub(pattern=c("Analysis of Deviance Table",
     #                                         "Analysis of Variance Table",
      #                                        "Model:",
       #                                       "Negative Binomial",
        #                                      "binomial",
         #                                     "Terms added sequentially \\(first to last\\)",
          #                                    "Response:",
           #                                   "link:"),
            #                        replacement=c(paste0("\n--------------------------------------",
             #                                            "-------------------------------------\n",
              #                                           c(mltext("print.anova.ml.KW.devTab"),
               #                                            mltext("print.anova.ml.KW.varTab"))),
                #                                  mltext("print.anova.ml.KW.family"),
                 #                                 mltext("print.anova.ml.KW.NB"),
                  #                                mltext("print.anova.ml.KW.B"),
                   #                               mltext("print.anova.ml.KW.termSeq"),
                    #                              mltext("print.anova.ml.KW.response"),
                     #                             mltext("print.anova.ml.KW.link")),
                      #              x=attr(x, "heading"), fixed=TRUE)

    ## Définitions issues de la fonction originale :
    if (!is.null(heading <- attr(x, "heading")))
    {
        cat("\n---------------------------------------------------------------------------\n", heading, sep = "\n",append=TRUE)
    }else{}

    nc <- dim(x)[2L]
    if (is.null(cn <- colnames(x)))
    {
        stop("'anova' object must have colnames")
    }else{}
    has.P <- grepl("^(P|Pr)\\(", cn[nc])
    zap.i <- 1L:(if (has.P)
             {
                 nc - 1
             }else{
                 nc
             })
    i <- which(substr(cn, 2, 7) == " value")
    i <- c(i, which(!is.na(match(cn, c("F", "Cp", "Chisq")))))
    if (length(i))
    {
        zap.i <- zap.i[!(zap.i %in% i)]
    }else{}

    tst.i <- i
    if (length(i <- grep("Df$", cn)))
    {
        zap.i <- zap.i[!(zap.i %in% i)]
    }else{}

    printCoefmat(x, digits = digits, signif.stars = signif.stars,
                 signif.legend=FALSE,
                 has.Pvalue = has.P, P.values = has.P, cs.ind = NULL,
                 zap.ind = zap.i, tst.ind = tst.i, na.print = "", ...)
    invisible(x)
}

######################################### end of the function print.anova.ml

######################################### start of the function printCoefmat.red called by signifParamLM.f

printCoefmat.red <- function(x, digits = max(3, getOption("digits") - 2),
                             signif.stars = getOption("show.signif.stars"),
                             signif.legend = signif.stars, dig.tst = max(1, min(5, digits - 1)),
                             cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(0),
                             P.values = NULL,
                             has.Pvalue = nc >= 4 &&
                                          substr(colnames(x)[nc], 1, 3) == "Pr(", eps.Pvalue = .Machine$double.eps,
                             na.print = "NA",
                             anovaLM=NULL,
                             objLM=NULL,
                             ...)
{
    ## Purpose: Modification de printCoefmat pour n'afficher que les z-values
    ##          et p-values, et pour les facteurs significatife uniquement.
    ## ----------------------------------------------------------------------
    ## Arguments: ceux de printCoefmat
    ##            + anovaLM : résultat d'anova globale du modèle (pour les
    ##                        facteurs et intéractions significatifs).
    ##            objLM : objet de classe 'lm' ou 'glm'
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 31 août 2010, 10:46

    ## Sélection des coefficients à montrer (pour effets/interactions significatifs) :
    #x <- x[selRowCoefmat(x, anovaLM, objLM), , drop=FALSE]

    ## Définitions issues de la fonction originale :
    if (is.null(d <- dim(x)) || length(d) != 2L)
        stop("'x' must be coefficient matrix/data frame")
    nc <- d[2L]
    if (is.null(P.values)) {
        scp <- getOption("show.coef.Pvalues")
        if (!is.logical(scp) || is.na(scp)) {
            warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
            scp <- TRUE
        }
        P.values <- has.Pvalue && scp
    }
    else if (P.values && !has.Pvalue)
        stop("'P.values' is TRUE, but 'has.Pvalue' is not")
    if (has.Pvalue && !P.values) {
        d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
        nc <- nc - 1
        has.Pvalue <- FALSE
    }
    else xm <- data.matrix(x)
    k <- nc - has.Pvalue - (if (missing(tst.ind))
        1
    else length(tst.ind))
    if (!missing(cs.ind) && length(cs.ind) > k)
        stop("wrong k / cs.ind")
    Cf <- array("", dim = d, dimnames = dimnames(xm))
    ok <- !(ina <- is.na(xm))
    for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
    if (length(cs.ind)) {
        acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
        if (any(ia <- is.finite(acs))) {
            digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
                floor(log10(range(acs[acs != 0], finite = TRUE)))
            else 0
            Cf[, cs.ind] <- format(round(coef.se, max(1, digits -
                digmin)), digits = digits)
        }
    }
    if (length(tst.ind))
        Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
            digits = digits)
    if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc))))
        for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
    okP <- if (has.Pvalue)
        ok[, -nc]
    else ok
    x1 <- Cf[okP]
    dec <- getOption("OutDec")
    if (dec != ".")
        x1 <- chartr(dec, ".", x1)
    x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
    if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
        Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1,
            digits - 1))
    }
    if (any(ina))
        Cf[ina] <- na.print
    if (P.values) {
        if (!is.logical(signif.stars) || is.na(signif.stars)) {
            warning("option \"show.signif.stars\" is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
        if (any(okP <- ok[, nc])) {
            pv <- as.vector(xm[, nc])
            Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst,
                eps = eps.Pvalue)
            signif.stars <- signif.stars && any(pv[okP] < 0.1)
            if (signif.stars) {
                Signif <- symnum(pv, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
                Cf <- cbind(Cf, format(Signif))
            }
        }
        else signif.stars <- FALSE
    }
    else signif.stars <- FALSE

    ## Sélection de colonnes :
    Cf <- Cf[ , ncol(Cf) - c(2:0)]

    print.default(Cf, quote = FALSE, right = TRUE, na.print = na.print,
        ...)
    if (signif.stars && signif.legend)
        cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    invisible(x)
}

######################################### end of the function printCoefmat.red

######################################### start of the function valPreditesLM.f called by sortiesLM.f

valPreditesLM.f <- function(objLM, Data, listFact, resFile)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : objet de classe 'lm' ou 'glm'.
    ##            Data : les données utilisées pour ajuster le modèle.
    ##            listFact : un vecteur donnant la liste des noms de
    ##                       facteurs.
    ##            resFile : la connection au fichier résultat
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  8 sept. 2010, 16:12


    ## ##################################################
    ## Valeurs prédites :
    OrdreNivFact <- sapply(unique(Data[ , listFact]), as.numeric)

    if (!is.matrix(OrdreNivFact))       # Si un seul facteur, on transforme le vecteur d'ordre des niveaux en matrice.
    {
        OrdreNivFact <- matrix(OrdreNivFact, ncol=1, dimnames=list(NULL, listFact))
    }else{}

    ## Valeurs prédites pour chaque combinaison réalisée des facteurs :
    if (length(grep("^glm", objLM$call)) > 0)
    {
        valPredites <- predict(objLM, newdata=unique(Data[ , listFact, drop=FALSE]), type="response")
    }else{
        valPredites <- predict(objLM, newdata=unique(Data[ , listFact, drop=FALSE]))
    }

    ## Noms des valeurs prédites (combinaisons des différents niveaux de facteurs) :
    nomCoefs <- unique(apply(Data[ , listFact, drop=FALSE], 1, paste, collapse=":"))
    names(valPredites) <- nomCoefs

    ## On remet les modalités en ordre :
    valPredites <- valPredites[eval(parse(text=paste("order(",
                                          paste("OrdreNivFact[ , ", 1:ncol(OrdreNivFact), "]", sep="", collapse=", "),
                                          ")", sep="")))]

    ## Écriture de l'en-tête :
    cat("\n\n\n---------------------------------------------------------------------------",
        "\n", "Values predicted by the model:", "\n\n",
        file=resFile,append=TRUE)

    ## Écriture du résultat :
    capture.output(print(valPredites), file=resFile,append=TRUE)

}

######################################### end of the function valPreditesLM.f

######################################### start of the function compMultiplesLM.f called by sortiesLM.f

compMultiplesLM.f <- function(objLM, Data, fact1, fact2, resFile, exclude="", Log=FALSE)
{
    ## Purpose: Calculer et écrire les résultats des comparaisons multiples.
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : objet de classe 'lm' ou 'glm'.
    ##            Data : les données utilisées pour ajuster le modèle.
    ##            fact1 : le nom du premier facteur utilisé pour les
    ##                    comparaisons multiples.
    ##            fact2 : le nom du second facteur utilisé pour les
    ##                    comparaisons multiples.
    ##            resFile : la connection pour les sorties textes.
    ##            exclude : facteur non analysé.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 oct. 2010, 09:54

    facts <- c(fact1, fact2)

    ## écriture des en-têtes :
    cat("\n\n\n---------------------------------------------------------------------------",
        "\nMultiple comparisons:",
        file=resFile,append=TRUE)

    ## Avertissement concernant les estimations de différences :
    compMultiplesAvertissement.f(objLM=objLM, Log=Log, resFile=resFile)

    ## Test si un facteur est temporel :
    tempFact <- is.temporal.f(facts, unitobs)

    ## Calculs des matrices de différences :
    for (i in seq(along=facts))
    {
        ## fact <- get(paste("fact", i, sep=""))
        if (tempFact[i])                # Si le facteur inclus est temporel :
        {
            assign(paste("diff", i, sep=""),
                   diffTemporelles.f(objLM=objLM,
                                     factSpatial=facts[-i],
                                     factTemp=facts[i],
                                     Data=Data,
                                     exclude=exclude))
        }else{                          # ... sinon :
            difftmp <- diffSpatiales.f(objLM=objLM,
                                       factSpatial=facts[i],
                                       factTemp=facts[-i],
                                       Data=Data,
                                       exclude=exclude)
            ## On réordonne d'après le second facteur (plus lisible) :
            assign(paste("diff", i, sep=""),
                   difftmp[order(sub("^([^:]+) :.+$", "\\1", row.names(difftmp))), ])
            rm(difftmp)
        }
    }

    ## Si des coefs n'ont pu être calculés, glht plante... à moins que :
    if (any(is.na(coef(objLM))))
    {
        ## Avertissement :
        cat("\n\n\t",
            "Warning: difference matrices reduced to account for",
            "\n\tnot calculable coefficients (missing data for some levels",
            "\n\tof factors/interactions).", "\n",
            file=resFile,append=TRUE)

        ## Réduction des matrices de différences :
        diff1 <- diff1[ , !is.na(coef(objLM))]
        diff2 <- diff2[ , !is.na(coef(objLM))]

        objLM$coefficients <- objLM$coefficients[!is.na(coef(objLM))]
    }

    for (i in seq(along=facts))
    {
        ## Résultats des comparaisons spatiales/de statut :
        cat(paste("\n\n", "Comparisons for differences in '", facts[i], "' ",
                  ifelse(tempFact[i],
                         paste0("(",
                                "temporal",
                                ") "),
                         ""),
                  "par '", facts[-i], "' ",
                  ifelse(tempFact[-i],
                         paste0("(",
                                "temporal",
                                ") "),
                         ""),
                  ":\n", sep=""),
            file=resFile,append=TRUE)
        #if (all(get(paste("diff", i, sep="") > 0)))
        #{
         #   capture.output(print.summary.glht.red(summary(glht(objLM,
          #                                                     linfct=get(paste("diff", i, sep="")),
           #                                                    alternative="two.sided"))),
            #               file=resFile,append=TRUE)
        #}else{}
    }
}

######################################### end of the function compMultiplesLM.f

######################################### start of the function is.temporal.f called by compMultiplesLM.f and compSimplesLM.f
is.temporal.f <- function(facteur, table)
{
    ## Purpose: test si un facteur est temporel ou non
    ## ----------------------------------------------------------------------
    ## Arguments: facteur : le nom (chaîne de caractères) du facteur.
    ##            table : la table dans laquelle se trouve le champ.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 oct. 2010, 10:01

    res <- sapply(facteur,
                  function(x)
              {
                  switch(x,
                         year={           # An est toujours censé être temporel.
                             TRUE
                         },
                         annee.campagne={           # Vérifié en amont.
                             TRUE
                         },
                         geogr.descriptor2={ # Dépend du format.
                             ifelse(all(grepl("^[cC]?[[:digit:]]{4}$",
                                              as.character(table[ , "geogr.descriptor2"])), na.rm=TRUE),
                                    TRUE,
                                    FALSE)
                         },
                         FALSE)
              })
    return(res)
}

######################################### end of the function is.temporal.f

######################################### start of the function compMultiplesAvertissement.f called by compMultiplesLM.f and compSimplesLM.f

compMultiplesAvertissement.f <- function(objLM, Log, resFile)
{
    ## Purpose: Afficher un avertissement concernant les différences (dans la
    ##          fonction de lien).
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : objet de classe (G)LM.
    ##            Log : booléen indiquant la log-transfomation des données.
    ##            resFile : fichier de sortie.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 31 janv. 2011, 14:11

    cat("\n",
        switch(modelType.f(objLM=objLM, Log=Log),
               "LM"={
                   ""
               },
               "LM-log"={
                   paste("\tWarning: differences are estimated on the logarithms:",
                         "\n\t(log(A) - log(B))", sep="")
               },
               "GLM-NB"={
                   paste("\tWarning: differences estimated in the link space (log):",
                         "\n\tlog(A) - log(B)", sep="")
               },
               "GLM-P"={
                   paste("\tWarning: differences estimated in the link space (log):",
                         "\n\tlog(A) - log(B)", sep="")
               },
               "GLM-B"={
                   paste("\tWarning: differences estimated in the link space (logit):",
                         "", sep="")
               },
               "GLM-Ga"={
                   paste("\tWarning: differences estimated in the link space (inverse):",
                         "\n\t(1/A) - (1/B)\t=>\t*", "inverse the sign of differences*", sep="")
               },
               ""),
        file=resFile,append=TRUE)
}

######################################### end of the function compMultiplesAvertissement.f 

######################################### start of the function modelType.f called by compMultiplesAvertissement.f 

modelType.f <- function(objLM, Log)
{
    ## Purpose: Fournir un prefix décrivant le modèle utilisé.
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : un objet de classe LM ou GLM.
    ##            Log : log-transformation des données (boolean).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 14 oct. 2010, 16:29

    return(ifelse(length(grep("^lm\\(", deparse(objLM$call), perl=TRUE)) > 0,
                  paste("LM", ifelse(Log, "-log", ""), sep=""),
                  ifelse(length(grep("^glm\\.nb", deparse(objLM$call), perl=TRUE)) > 0,
                         "GLM-NB",
                         ifelse(length(grep("^glm.*poisson", deparse(objLM$call), perl=TRUE)) > 0,
                                "GLM-P",
                                ifelse(length(grep("^glm.*\"binomial\"", deparse(objLM$call), perl=TRUE)) > 0,
                                       "GLM-B",
                                       ifelse(length(grep("family[[:blank:]]*=[[:blank:]]*\"Gamma\"", deparse(objLM$call), perl=TRUE)) > 0,
                                              "GLM-Ga",
                                              "Unknown-model"))))))
}

######################################### end of the function modelType.f 

######################################### start of the function diffTemporelles.f called by compMultiplesLM.f

diffTemporelles.f <- function(objLM, factSpatial, factTemp, Data, exclude)
{
    ## Purpose: Calcule et retourne la matrice de différences temporelles par
    ##          statut(pour une utilisation avec la fonction 'glht').
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : un objet de classe 'lm' ou 'glm'.
    ##            factSpatial : nom du facteur spatial.
    ##            factTemp : nom du facteur temporel.
    ##            Data : données utilisées pour ajuster le modèle.
    ##            exclude : facteur non analysé.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  8 sept. 2010, 11:11

    ## Coefficients :
    if (length(grep("^glmmTMB", objLM$call)) > 0)
    {
        theta <- c(levels(objLM$frame[,factTemp]),levels(objLM$frame[,factSpatial]))
    }else{
        theta <- names(coef(objLM))
    }

    tDiff <- paste(c(head(rev(levels(Data[ , factTemp])), 1), head(rev(levels(Data[ , factTemp])),  - 1)),
                   c(tail(rev(levels(Data[ , factTemp])), 1), tail(rev(levels(Data[ , factTemp])),  - 1)),
                   sep=" - ")

    ## Matrice pour construire les différences entre coefficients :
    Dtemp <- matrix(0,
                    ## Il y a autant de différences temporelles que de niveaux pour la variable temporelle (en raison de
                    ## la différence supplémentaire final - initial) :
                    nrow=nlevels(Data[ , factTemp]) * nlevels(Data[ , factSpatial]),
                    ncol=length(theta))

    ## Noms des colonnes (pas obligatoire mais utile pour vérification) :

    row.names(Dtemp) <- paste(rep(levels(Data[ , factSpatial]), each=nlevels(Data[ , factTemp])),
                              tDiff,
                              sep=" : ")


    colnames(Dtemp) <- theta

    ## Noms de colonnes ordonnés :
    if (length(grep("^glmmTMB", objLM$call)) > 0)
    {
        namesCol <- c(colnames(Data)[1],
                      colnames(objLM$frame[ , -1]))
    }else{
        namesCol <- c(colnames(Data)[1],
                      colnames(objLM$model[ , -1]))
    }

    namesCol <- namesCol[! is.element(namesCol, exclude)] # - colonne exclue

    ## Calculs des nombres de colonnes des facteurs et intéraction :
    nlev <- combn(sapply(Data[ , namesCol],
                         function(x)
                     {
                         ifelse(is.factor(x),
                                nlevels(x) - 1,
                                1)
                     }),
                  2)

    ## Nombre de colonnes par type de facteur/interaction :
    nCol <- apply(nlev, 2, prod)

    ## Position de la première colonne
    premiereCol <- cumsum(c(1, nCol[- length(nCol)])) + 1

    ## Position des facteurs d'intérêt et leur interaction,
    ## dans l'ordre de l'ensemble des facteurs et interactions :
    facts <- c(factSpatial, factTemp)
    posTemp <- 1
    posSpatial <- 2
    #posInteraction <- which(is.element(attr(objLM$terms, "term.labels"),
     #                                  paste(facts, rev(facts), sep=":")))

    ## Différences sur l'effet temporel seul :
    d1 <- rbind(c(-1, rep(0, nCol[posTemp] - 1), 1),
                cbind(0, diag(1, nCol[posTemp])[ , seq(nCol[posTemp], 1)]) +
                cbind(diag(-1, nCol[posTemp])[ , seq(nCol[posTemp], 1)], 0))[ , -1]


    Dtemp[ , seq(from=premiereCol[posTemp],
                 length.out=nCol[posTemp])] <- sapply(as.data.frame(d1), rep, nlevels(Data[ , factSpatial]))


    ## Différences sur les interactions :
#    d2 <- Dtemp[ , seq(from=premiereCol[posInteraction],
 #                       length.out=nCol[posInteraction]), drop=FALSE]

  #  l <- nlevels(Data[ , factTemp]) + 1
   # for (i in seq(from=0, length.out=nCol[posSpatial]))
    #{
     #   if (posSpatial > posTemp)       # traitement différent selon l'imbrication des facteurs :
      #  {                               # Cas où le facteur temporel est en premier :
       #     d2[seq(from=l, length.out=nlevels(Data[ , factTemp])) ,
        #       seq(from=1, length.out=nCol[posTemp]) + i * nCol[posTemp]] <- d1
#        }else{                          #... cas où il est en second :
 #           d2[seq(from=l, length.out=nlevels(Data[ , factTemp])) ,
  #             seq(from=1 + i, by=nCol[posSpatial], length.out=nCol[posTemp])] <- d1
   #     }
#
 #       l <- l + nlevels(Data[ , factTemp])
  #  }
#
 #   Dtemp[ , seq(from=premiereCol[posInteraction],
  #               length.out=nCol[posInteraction])] <- d2

    return(Dtemp)

}

######################################### end of the function diffTemporelles.f

######################################### start of the function diffSpatiales.f called by compMultiplesLM.f

diffSpatiales.f <- function(objLM, factSpatial, factTemp, Data, exclude)
{
    ## Purpose: Calcule et retourne la matrice de différences spatiales par
    ##          année (pour une utilisation avec la fonction 'glht').
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : un objet de classe 'lm' ou 'glm'.
    ##            factSpatial : nom du facteur spatial.
    ##            factTemp : nom du facteur temporel.
    ##            Data : données utilisées pour ajuster le modèle.
    ##            exclude : facteur non analysé.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  7 sept. 2010, 16:15

    ## Coefficients :
    if (length(grep("^glmmTMB", objLM$call)) > 0)
    {
        theta <- c(levels(objLM$frame[,factTemp]),levels(objLM$frame[,factSpatial]))
    }else{
        theta <- names(coef(objLM))
    }
    ## Nom des différences spatiales (statut de protection) :

    sDiff <- apply(combn(unique(Data[ , factSpatial]), 2),
                   2,
                   function(x){paste(rev(x), collapse = " - ")})

    ## Matrice pour construire les différences entre coefficients :
    Dspat <- matrix(0,
                    nrow=nlevels(Data[ , factTemp]) * choose(nlevels(Data[ , factSpatial]), 2),
                    ncol=length(theta))

    ## Noms des colonnes (pas obligatoire mais utile pour vérification) :
    #row.names(Dspat) <- paste(levels(Data[ , factTemp]),
     #                         rep(sDiff, each=nlevels(Data[ , factTemp])), sep=" : ")
    colnames(Dspat) <- theta

    ## Noms de colonnes ordonnés :
    if (length(grep("^glmmTMB", objLM$call)) > 0)
    {
        namesCol <- c(colnames(Data)[1],
                      colnames(objLM$frame[ , -1]))
    }else{
        namesCol <- c(colnames(Data)[1],
                      colnames(objLM$model[ , -1]))
    }

    namesCol <- namesCol[! is.element(namesCol, exclude)] # - colonne exclue

    ## Calculs des nombres de colonnes des facteurs et intéraction :
    nlev <- combn(sapply(Data[ , namesCol],
                         function(x)
                     {
                         ifelse(is.factor(x),
                                nlevels(x) - 1,
                                1)
                     }),
                  2)

    ## Nombre de colonnes par type de facteur/interaction :
    nCol <- apply(nlev, 2, prod)

    ## Position de la première colonne
    premiereCol <- cumsum(c(1, nCol[- length(nCol)])) + 1

    ## Position des facteurs d'intérêt et leur interaction,
    ## dans l'ordre de l'ensemble des facteurs et interactions :
    facts <- c(factSpatial, factTemp)
    posTemp <- 1
    posSpatial <- 2
    #posInteraction <- which(is.element(attr(objLM$terms, "term.labels"),
     #                                  paste(facts, rev(facts), sep=":")))

    ## Différences entres les effets statuts (sans intéraction temporelles) :

    tmp <- sapply(as.data.frame(combn(1:nlevels(Data[ , factSpatial]), 2)),
                  function(x)
              {
                  m <- matrix(0,
                              ncol=nlevels(Data[ , factSpatial]),
                              nrow=nlevels(Data[ , factTemp]))
                  m[ , x] <- matrix(c(-1, 1),
                                    nrow=nlevels(Data[ , factTemp]),
                                    ncol=2,
                                    byrow=TRUE)
                  return(m)
              }, simplify=FALSE)

    m <- tmp[[1]][NULL, ]
    for(i in 1:length(tmp))
    {
        m <- rbind(m, tmp[[i]])
    }

    Dspat[ , premiereCol[posSpatial] - 1 + 1:nCol[posSpatial]] <- m[ , -1]

    ## Ajout des intéractions :
  #  tmp2 <- Dspat[ , seq(from=premiereCol[posInteraction], length.out=nCol[posInteraction]), drop=FALSE]
#
  #  l <- 1
   # for (i in as.data.frame(combn(0:nCol[posSpatial], 2))) # pour chaque combinaison de statut :
    #{
     #   if(i[1] != 0)
      #  {
       #     d1 <- rbind(0, diag(-1, nrow=nCol[posTemp]))
        #    if (posSpatial > posTemp)   # facteur spatial après le facteur temporel...
         #   {
          #      tmp[seq(from=l, length.out=nlevels(Data[ , factTemp])),
           #          seq(from=(i[1] - 1) * nCol[posTemp] + 1, length.out=nCol[posTemp])] <- d1
            #}else{                      # ... avant le facteur temporel.
             #   tmp[seq(from=l, length.out=nlevels(Data[ , factTemp])),
              #       seq(from=i[1], by=nCol[posSpatial] , length.out=nCol[posTemp])] <- d1
            #}
        #}else{}

      #  d2 <- rbind(0, diag(1, nrow=nCol[posTemp]))

       # if (posSpatial > posTemp)       # facteur spatial après le facteur temporel...
        #{
         #   stop(tmp)
      #      tmp[seq(from=l, length.out=nlevels(Data[ , factTemp])),
       #          seq(from=(i[2] - 1) * nCol[posTemp] + 1, length.out=nCol[posTemp])] <- d2
        #}else{                          # ... avant le facteur temporel.
         #   tmp[seq(from=l, length.out=nlevels(Data[ , factTemp])),
          #       seq(from=i[2], by=nCol[posSpatial], length.out=nCol[posTemp])] <- d2
        #}

     #   l <- l + nlevels(Data[ , factTemp])
    #}
#
 #   ## Stockage des différences d'interactions :
  #  Dspat[ , seq(from=premiereCol[posInteraction], length.out=nCol[posInteraction])] <- tmp2

    return(Dspat)
}

######################################### end of the function diffSpatiales.f

######################################### start of the function print.summary.glht.red called by compMultiplesLM.f and compSimplesLM.f

print.summary.glht.red <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    ## cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(x$observation.type))
        cat("Multiple comparisons of the means:", x$observation.type,
            "Contrasts", "\n\n\n")
    call <- if (isS4(x$model))
            {
                x$model@call
            }else{
                x$model$call
            }
    ## if (!is.null(call)) {
    ##     cat("Fit: ")
    ##     print(call)
    ##     cat("\n")
    ## }
    cat("\n")
    pq <- x$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    error <- attr(pq$pvalues, "error")
    pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df ==
        0, "z", "t"), ")", sep = ""), greater = paste("Pr(>",
        ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|",
        ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
    colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==
        0, "z value", "t value"), pname)
    type <- ifelse(is.null(pq$observation.type) & ! is.null(pq$type), pq$type, pq$observation.type)
    if (!is.null(error) && error > .Machine$double.eps)
    {
        sig <- which.min(abs(1/error - (10^(1:10))))
        sig <- 1/(10^sig)
    }else{
        sig <- .Machine$double.eps
    }
    cat("Linear hypotheses:", "\n")
    alt <- switch(x$alternative, two.sided = "==", less = ">=",
        greater = "<=")
    rownames(mtests) <- paste(rownames(mtests), alt, x$rhs)
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE,
                 P.values = TRUE, eps.Pvalue = sig)
    switch(type, univariate = cat("(Univariate P-values)"),
        `single-step` = cat("(Adjusted P-values -- 'single step' method)"),
        Shaffer = cat("(Adjusted P-values -- Shaffer method)"),
        Westfall = cat("(Adjusted P-values -- Westfall method)"),
        cat("(Adjusted P-values --", type, "method)"))
    cat("\n\n")
    invisible(x)
}

######################################### end of the function print.summary.glht.red

######################################### start of the function graphTitle.f called by sortiesLM.f

graphTitle.f <- function(metrique, modGraphSel, factGraph, listFact, model=NULL, type="espece",
                         lang = getOption("P.lang"))
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 14 oct. 2010, 15:44
    return(paste(ifelse(is.null(model),
                        "Values of ",
                        paste(model,
                              " for",
                              sep="")),
                 metrique,
                 ifelse(is.element(type, c("espece", "unitobs", "CL_espece", "unitobs(CL)")),
                        paste("aggregated"),
                        ""),
                 switch(type,
                        "espece"=" per species and station",
                        "CL_espece"=" per size class, species and station",
                        "unitobs"=" per station",
                        "unitobs(CL)"=" per station",
                        "CL_unitobs"=" per size class and station",
                        "biodiv"=" per station",
                        ""),
                 switch(type,
                        "espece"={
                            ifelse(modGraphSel == "", # Facteur de séparation uniquement si défini.
                                   "",
                                   paste("\nfor the field",
                                         " '", factGraph, "' = ", modGraphSel, sep=""))
                        },
                        "CL_espece"={
                            ifelse(modGraphSel == "", # Facteur de séparation uniquement si défini.
                                   "",
                                   paste("\nfor the field",
                                         " '", factGraph, "' = ", modGraphSel, sep=""))
                        },
                        "unitobs"={
                            ifelse(modGraphSel[1] == "", # Facteur de séparation uniquement si défini.
                                   "\nfor all species",
                                   paste("\nfor all species matching",
                                         " '", factGraph, "' = (",
                                         paste(modGraphSel, collapse=", "), ")", sep=""))
                        },
                        "unitobs(CL)"={
                            ifelse(modGraphSel[1] == "", # Facteur de séparation uniquement si défini.
                                   "\nfor all size classes",
                                   paste("\nfor size classes matching",
                                         " '", factGraph, "' = (",
                                         paste(modGraphSel, collapse=", "), ")", sep=""))
                        },
                        "CL_unitobs"={
                            ifelse(modGraphSel[1] == "", # Facteur de séparation uniquement si défini.
                                   "\nfor all species",
                                   paste("\nfor all species matching",
                                         " '", factGraph, "' = (",
                                         paste(modGraphSel, collapse=", "), ")", sep=""))
                        },
                        "biodiv"={
                            ifelse(modGraphSel[1] == "", # Facteur de séparation uniquement si défini.
                                   "",
                                   paste("\nfor stations matching",
                                         " '", factGraph, "' = (",
                                         paste(modGraphSel, collapse=", "), ")", sep=""))
                        },
                        ""),
                 "\n by ",
                 paste(sapply(listFact[length(listFact):1],
                              function(x)paste(c(## varNames.f(x, "article"),
                                                 "",
                                                 varNames.f(x, "nom")), collapse="")),
                       collapse=" and"),
                 "\n", sep=""))
}

######################################### end of the function graphTitle.f

######################################### start of the function Capitalize.f called by sortiesLM.f

Capitalize.f <- function(x, words=FALSE)
{
    ## Purpose: Mettre en majuscule la première lettre de chaque mot
    ## ----------------------------------------------------------------------
    ## Arguments: x : une chaîne de caractères
    ##            words : tous les mots (TRUE), ou juste le premier.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  9 août 2010, 21:08

    if (words)
    {
        s <- strsplit(x, " ")[[1]]
    }else{
        s <- x
    }

    return(paste(toupper(substring(s, 1,1)), substring(s, 2),
                 sep="", collapse=" "))
}

######################################### end of the function Capitalize.f

######################################### start of the function compSimplesLM.f called by sortiesLM.f

compSimplesLM.f <- function(objLM, Data, fact, resFile, Log=FALSE)
{
    ## Purpose: Calculer et écrire les résultats des comparaisons simples.
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : objet de classe 'lm' ou 'glm'.
    ##            Data : les données utilisées pour ajuster le modèle.
    ##            fact : le nom du facteur utilisé.
    ##            resFile : la connection pour les sorties textes.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 oct. 2010, 15:51

    ## écriture des en-têtes :
    cat("\n\n\n---------------------------------------------------------------------------",
        "\n", "Comparisons of levels:",
        file=resFile,append=TRUE)

    ## Avertissement concernant les estimations de différences :
    compMultiplesAvertissement.f(objLM=objLM, Log=Log, resFile=resFile)

    if (is.temporal.f(fact, unitobs))
    {
        ## Suite en-tête :
        cat(paste("\n\n\t", "Factor",
                  " '", fact, "' (",
                  "temporal",
                  ") :\n", sep=""),
            file=resFile,append=TRUE)

        ## Comparaisons temporelles :
        compSimple <- glht(objLM,
                           linfct=diffTempSimples.f(objLM=objLM, fact=fact, Data=Data),
                           alternative="two.sided")

        ## Écriture des résultats :
        capture.output(print.summary.glht.red(summary(compSimple)),
                   file=resFile,append=TRUE)

    }else{
        ## Suite en-tête :
        cat(paste("\n\n", "Factor",
                  " '", fact, "' :\n", sep=""),
            file=resFile,append=TRUE)

        ## Comparaisons de toutes les paires ("Tukey") :
        compSimple <- glht(objLM,
                           linfct=eval(parse(text=paste("mcp(", fact, "=\"Tukey\")"))),
                           alternative="two.sided")

        ## Écriture des résultats :
        capture.output(print.summary.glht.red(summary(compSimple)),
                   file=resFile,append=TRUE)
    }
}

######################################### end of the function compSimplesLM.f

######################################### start of the function diffTempSimples.f called by compSimplesLM.f
diffTempSimples.f <- function(objLM, fact, Data)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments: objLM : un objet de classe 'lm' ou 'glm'.
    ##            factSpatial : nom du facteur spatial.
    ##            factTemp : nom du facteur temporel.
    ##            Data : données utilisées pour ajuster le modèle.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 oct. 2010, 16:28

    tDiff <- paste(c(head(rev(levels(Data[ , fact])), 1), head(rev(levels(Data[ , fact])),  - 1)),
                   c(tail(rev(levels(Data[ , fact])), 1), tail(rev(levels(Data[ , fact])),  - 1)),
                   sep=" - ")

    ## Coefficients :
    if (length(grep("^glmmTMB", objLM$call)) > 0)
    {
        theta <- levels(objLM$frame[,fact])
    }else{
        theta <- names(coef(objLM))
    }
    diffDim <- length(theta)

    diffMat <- matrix(0, ncol=diffDim, nrow=diffDim,
                      dimnames=list(tDiff, theta))


    ## Différences sur l'effet temporel seul :
    diffMat[ , -1] <- rbind(c(-1, rep(0, diffDim - 2), 1),
                            cbind(0, diag(1, diffDim - 1)[ , seq(diffDim - 1, 1)]) +
                            cbind(diag(-1, diffDim - 1)[ , seq(diffDim - 1, 1)], 0))[ , -1]

    return(diffMat)
}
######################################### end of the function diffTempSimples.f 

######################################### start of the function infoStats.f called by sortiesLM.f

infoStats.f <- function(filename, Data, agregLevel=c("species", "unitobs"), type=c("graph", "stat"),
                        metrique, factGraph, factGraphSel, listFact, listFactSel,
                        dataEnv, baseEnv=.GlobalEnv)
{
    ## Purpose: Écrire les infos et statistique sur les données associées à
    ##          un graphique ou analyse.
    ## ----------------------------------------------------------------------
    ## Arguments: filename : chemin du fichier de résultats.
    ##            Data : données du graphique/de l'analyse.
    ##            agregLevel : niveau d'agrégation de la fonction appelante.
    ##            type : type de fonction appelante (grapique ou analyse).
    ##            metrique : la métrique choisie.
    ##            factGraph : le facteur sélection des espèces.
    ##            factGraphSel : la sélection de modalités pour ce dernier
    ##            listFact : liste du (des) facteur(s) de regroupement
    ##            listFactSel : liste des modalités sélectionnées pour ce(s)
    ##                          dernier(s)
    ##            dataEnv : environnement de stockage des données.
    ##            baseEnv : environnement de l'interface.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 10 sept. 2012, 15:26

    ## Ajout de l'extension si besoin :
    #if ( ! grepl("\\.stats$", filename[1], ignore.case=TRUE))
    #{
     #   filename <- paste(filename, ".stats", sep="")
    #}else{}

    ## Ouverture du fichier :
    File <- file(description=filename,
                 open="w", encoding="latin1")

    ## Si erreur, on referme le fichier à la sortie de fonction :
    on.exit(if (exists("filename") &&
                tryCatch(isOpen(File),
                         error=function(e)return(FALSE))) close(File))

    ## Informations générales sur les données :
    #printGeneralDataInfo.f(dataEnv=dataEnv, baseEnv=baseEnv, File=File) ##Utilise les données de l'environnement et le tcltk donc à voir comment corriger

    ## Informations sur les métriques et facteurs du graphique :
    printSelectionInfo.f(metrique=metrique, #factGraph=factGraph, factGraphSel=factGraphSel,
                         listFact=listFact, #listFactSel=listFactSel, 
                         File=File,
                         agregLevel=agregLevel, type=type)

    ## Statistiques :
    if (class(Data) == "list")
    {
        cat("\n###################################################",
            "\nStatistics per level of splitting factor:\n",
            sep="", file=File,append=TRUE)

        invisible(sapply(1:length(Data),
                         function(i)
                     {
                         printStats.f(Data=Data[[i]], metrique=metrique, listFact=listFact, File=File,
                                      headline=factGraphSel[i])
                     }))
    }else{
        printStats.f(Data=Data, metrique=metrique, listFact=listFact, File=File,
                     headline=NULL)
    }

    ## Fermeture du fichier :
    close(File)

}

######################################### end of the function infoStats.f

######################################### start of the function printGeneralDataInfo.f called by infoStats.f

#printGeneralDataInfo.f <- function(dataEnv, baseEnv, File)
#{
    ## Purpose: Écrire dans un fichier les informations générales sur le jeu
    ##          de données (inclue les sélections au niveau de la
    ##          plateforme).
    ## ----------------------------------------------------------------------
    ## Arguments: dataEnv : environnement des données.
    ##            baseEnv : environnement de l'interface principale.
    ##            File : connection du fichier où écrire les informations.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 11 sept. 2012, 10:36

    ## Informations sur les fichiers de données :
    #cat(paste("####################\nDatasets:\n",
     #         "\n  * Identification of the case study: ", paste(getOption("P.MPA"), collapse=", "),
              #"\n  * Data directory: ", dataEnv$fileNames["ws"], "/Data/",
      #        mltext("printGeneralDataInfo.f.4"), dataEnv$fileNames["obs"],
       #       mltext("printGeneralDataInfo.f.5"), dataEnv$fileNames["unitobs"],
        #      mltext("printGeneralDataInfo.f.6"), dataEnv$fileNames["refesp"],
         #     mltext("printGeneralDataInfo.f.7"), dataEnv$fileNames["refspa"], "\n",
          #    sep=""),
        #file=File,append=TRUE)

    ## Sélections au niveau de la plateforme :
    #cat(ifelse((tmp <- evalq(tclvalue(tkcget(MonCritere, "-text")), envir=.baseEnv)) == "Tout",
     #          "\nNo general selection on data.\n",
      #         paste("\nGeneral selection(s):\n\n", tmp, "\n", sep="")),
       # file=File,append=TRUE)
#}

######################################### end of the function printGeneralDataInfo.f

######################################### start of the function printSelectionInfo.f called by infoStats.f

printSelectionInfo.f <- function(metrique, factGraph, factGraphSel, listFact, listFactSel,
                                 File,
                                 agregLevel=c("species", "unitobs"), type=c("graph", "stat"))
{
    ## Purpose: Écrire dans un fichier les informations sur la sélection de
    ##          données (obtenue par l'interface standard de sélection).
    ## ----------------------------------------------------------------------
    ## Arguments: metrique : la métrique choisie.
    ##            factGraph : le facteur sélection des espèces.
    ##            factGraphSel : la sélection de modalités pour ce dernier
    ##            listFact : liste du (des) facteur(s) de regroupement
    ##            listFactSel : liste des modalités sélectionnées pour ce(s)
    ##                          dernier(s)
    ##            File : connection du fichier où écrire les informations.
    ##            agregLevel : niveau d'agrégation de la fonction appelante.
    ##            type : type de fonction appelante (grapique ou analyse).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 11 sept. 2012, 10:41

    cat("\n##################################################\n",
        "Metrics and factors (and possible units/selections):\n",
        sep="", file=File,append=TRUE)

    ## Informations sur la métrique :
    cat("\n Metrics:", metrique,
        "\n", file=File,append=TRUE)

    ## Niveau d'agrégation :
    cat("            aggregated per ",
        switch(agregLevel,
               "CL_espece"=,"CL_unitobs"=,"spCL_unitobs"=,"spCL_espece"={
                   "size class / "
               }),
        switch(agregLevel,
               "CL_espece"=,"spCL_espece"=,"species"=,"spSpecies"=,"spEspece"={
                   "species / "
               }),
        switch(agregLevel,
               "spUnitobs"=,"spCL_unitobs"=,"spCL_espece"=,"spUnitobs(CL)"=,"spSpecies"=,"spEspece"={
                   paste(listFact, " (mean over ", sep="")
              }),
        "observation units",
        switch(agregLevel,
               "spUnitobs"=,"spCL_unitobs"=,"spCL_espece"=,"spUnitobs(CL)"=,"spSpecies"=,"spEspece"={
                   ")"
              }),
        ".\n",
        sep="", file=File,append=TRUE)

    ## Facteurs de séparation de graphiques/analyses ou sélection d'observations :
#    switch(agregLevel,
 #          "species"=,"CL_espece"=,"espece"={ # Adapté également pour les LMs.
  #             cat("\n",
   #                switch(type,
    #                      "graph"="Graphics separation factor",
     #                     "stat"="Analyses separation factor"),
      #             " : ",
       #            ifelse(factGraph == "", "printSelectionInfo.f.11",
        #                  ifelse(is.na(factGraphSel[1]),
         #                        paste(varNames.f(factGraph, "nom"), "none!"),
          #                       paste(varNames.f(factGraph, "nom"), " (",
           #                            paste(factGraphSel, collapse=", "), ")", sep=""))), "\n",
            #       sep="", file=File,append=TRUE)
#           },
 #          "unitobs"=,"CL_unitobs"=,"unitobs(CL)"=,"spUnitobs"={
  #             cat("(warning: no selection!!!)",
   #                ifelse(factGraph == "", "\nSelection factor for aggregation of observations: ",
    #                      ifelse(is.na(factGraphSel[1]),
     #                            paste(varNames.f(factGraph, "nom"), "none (all species/size classes)!"),
      #                           paste(varNames.f(factGraph, "nom"), " (",
       #                                paste(factGraphSel, collapse=", "), ")", sep=""))), "\n",
        #           sep="", file=File,append=TRUE)
         #  })

    ## Facteurs de regroupements :
    if (is.element(agregLevel, c("spCL_unitobs", "spCL_espece", "spSpecies", "spEspece",
                                 "spUnitobs", "spUnitobs(CL)"))) {type <- "spatialGraph"}

    cat(switch(type,
               "graph"="\nGrouping factor(s): \n * ",
               "stat"="\nAnalyses factor(s): \n * ",
               "spatialGraph"="\nSpatial aggregation factor(s): \n * "),
        paste(listFact,collaspe="\n * "),"\n",file=File,append=TRUE)

#    invisible(sapply(1:length(listFact),
 #                    function(i)
  #               {
   #                  cat("\n  * ",
    #                     ifelse(is.na(listFactSel[[i]][1]),
     #                                  paste(varNames.f(listFact[i], "nom"), "(no selection)"),
      #                                 paste(varNames.f(listFact[i], "nom"), " (",
       #                                      paste(listFactSel[[i]], collapse=", "), ")", sep="")), "\n",
        #                 sep="", file=File,append=TRUE)
         #        }))
}

######################################### end of the function printSelectionInfo.f

######################################### start of the function varNames.f called by printSelectionInfo.f

varNames.f <- function(fields, info="name", quote=TRUE)
{
    ## Purpose: revoyer les informations (en particulier nom) sur le nom
    ##          "d'usage" d'un ou plusieurs champ(s).
    ## ----------------------------------------------------------------------
    ## Arguments: fields : champ(s) recherché(s).
    ##            info : type d'info ("name", "article", "gender", "unit")
    ##            quote : faut-il mettre des guillemets pour les noms de
    ##                    champs tels-quels (pas de nom d'usage défini).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 21 févr. 2013, 18:21

    info <- info[1]

    if (is.element(info, c("nom", "name")))
    {
        ## S'il n'est pas définit, le nom d'usage est remplacé par le nom de champ plutôt que par NA :
        res <- ifelse(is.na(tmp <- varNames[fields, "nom"]),
                      paste(ifelse(quote, "\"", ""),
                            fields,
                            ifelse(quote, "\"", ""), sep=""),
                      tmp)
    }else{
        ## Possibilité de nommer les infos en français et anglais:
        res <- ifelse(is.na(varNames[fields, info]),
                      "",
                      varNames[fields,
                               switch(info,
                                      "article"="article",
                                      "genre"=,
                                      "gender"="genre",
                                      "unite"=,
                                      "unit"="unite",
                                      "nom")])
    }

    return(res)
}

######################################### end of the function varNames.f

######################################### start of the function printStats.f called by infoStats.f

printStats.f <- function(Data, metrique, listFact, File, headline=NULL)
{
    ## Purpose: Écrire les tableaux de statistiques générales et par
    ##          croisement de facteur dans un fichier.
    ## ----------------------------------------------------------------------
    ## Arguments: Data : les données du graphique/de l'analyse.
    ##            metrique : nom de la métrique.
    ##            listFact : liste des facteurs de regroupement/de l'analyse.
    ##            File : la connection du fichier où écrire.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 11 sept. 2012, 10:09

    ## Ligne d'en-tête (si besoin : traitement par espèces uniquement) :
    if ( ! is.null(headline))
    {
        cat("\n", rep("#", nchar(headline) + 3), "\n",
            "## ", headline, "\n",
            sep="", file=File,append=TRUE)
    }else{}

    cat("\n########################\nBase statistics:\n\n", file=File,append=TRUE)

    capture.output(print(summary.fr(Data[ , metrique])), file=File, append=TRUE)

    if ( ! is.null(listFact))
    {
        cat("\n#########################################",
            "\nStatistics per combination of factor levels:\n\n", file=File, sep="",append=TRUE)

        ## Calcul du summary pour chaque croisement (existant) de facteur :
        res <- with(Data,
                    tapply(eval(parse(text=metrique)),
                           INDEX=do.call(paste,
                                         c(lapply(listFact,
                                                  function(y)eval(parse(text=y))),
                                           sep=".")),
                           FUN=summary.fr))

        ## Assemblage du résultat dans un tableau
        capture.output(print(do.call(rbind, res)),
                       file=File, append=TRUE)
    }else{}

    ## Ligne vide (pour l'esthétique) :
    cat("\n", file=File,append=TRUE)
}

######################################### end of the function printStats.f

######################################### start of the function errorLog.f called by modeleLineaireWP2.unitobs.f in FunctExeCalcGLMCalaxy.r

errorLog.f <- function(error, niv=-3)
{
    ## Purpose: Écrire les erreurs dans un fichier log + avertissement de
    ##          l'utilisateur
    ## ----------------------------------------------------------------------
    ## Arguments: error : erreur (récupérée par la fonction tryCatch).
    ##            niv : niveau de l'appel pour retrouver la fonction
    ##                  appelante (-3 par défaut pour tryCatch).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 22 déc. 2010, 11:54

    on.exit(if (exists("logFile") &&
                tryCatch(isOpen(logFile),
                         error=function(e)return(FALSE))) close(logFile))

    ## Test d'existance et éventuelle création du dossier de logs :
    #if (!isTRUE(file.info("./logs")$isdir))
     # {
      #    dir.create("./logs")
      #}

    ## Test d'existance et éventuelle création du fichier de log du jour :
    logFileName <- "Errors.txt"

#    if (!file.exists(paste("./logs/", logFileName, sep="")) ||
 #       isTRUE(file.info(paste("./logs/", logFileName, sep=""))$isdir))
  #    {
   #       file.create(paste("./logs/", logFileName, sep=""))
    #  }

    #logFile <- file(description=paste("./logs/", logFileName, sep=""),
     #               open="a", encoding="latin1")


    callingFct <- sys.call(niv)

    cat(paste("\n", format(Sys.time(), "[%H:%M:%S]"), "\n",
              paste(deparse(callingFct), collapse="\n\t"), " :\n", sep=""),
        file=logFile,append=TRUE)
    capture.output(print(error), file=logFile,append=TRUE)
    cat("\n", file=logFile,append=TRUE)

    close(logFile)

    message("\n\tThere was an error.", "\n\tPlease see the log file: ", logFileName, "\n")
}

######################################### end of the function errorLog.f

######################################### start of the function summary.fr called by printStats.f
summary.fr <- function(object, digits = max(3, getOption("digits") - 3),...)
{
    ## Purpose: Adding SD and N to summary
    ## ----------------------------------------------------------------------
    ## Arguments: object : objet à résumer.
    ##            ... : argument supplémentaires passés à summary().
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 13 sept. 2012, 15:47

    if ( ! is.numeric(object)) stop("Programming error")

    ## Calcul du résumé :
    res <- c(summary(object=object, digits, ...), "sd"=signif(sd(x=object), digits=digits), "N"=length(object))

    return(res)
}

######################################### start of the function summary.fr


