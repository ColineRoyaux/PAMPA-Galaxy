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

######################################### start of the function fact.def.f called by FunctExeCalcCommIndexesGalaxy.r
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
        warning(paste(mltext("agregations.generic.war.1"),
                      mltext("agregations.generic.war.2"), sep=""))
        return(NULL)
    }else{
        return(reslong)
    }
}

######################################### end of the function agregations.generic.f

######################################### start of the function dropLevels.f called y calcBiodiv.f in FucntExeCalcCommIndexesGalaxy.r
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
    for (i in selCol)
    {
        restmp <- subset(restmp, is.element(restmp[ , facteurs[i]], selections[[i]]))
    }

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
    resFile <- resFileLM.f(objLM=objLM, metrique=metrique, factAna=factAna, modSel=modSel, listFact=listFact,
                           dataEnv=dataEnv, Log=Log, sufixe=sufixe, type=type)
    on.exit(tryCatch(close(resFile), error=function(e){}), add=TRUE)


    ## Informations et statistiques globales sur le modèle :
    infoStatLM.f(objLM=objLM, resFile=resFile)


    ## Anova globale du modèle + significativité des coefficients :
    signifParamLM.f(objLM=objLM, resFile=resFile)


    ## ##################################################
    ## Valeurs prédites par le modèle :
    valPreditesLM.f(objLM=objLM, Data=Data, listFact=listFact, resFile=resFile)

    ## ##################################################
    ## Comparaisons multiples :

    ## if (all(is.element(c("year", "protection.status"), listFact)))
    if (length(listFact) == 2)
    {
        WinInfo <- tktoplevel()
        on.exit(tkdestroy(WinInfo))
        tkwm.title(WinInfo, mltext("sortiesLM.W.Title"))

        tkgrid(tklabel(WinInfo, text="\t "),
               tklabel(WinInfo, text=paste0("\n", mltext("sortiesLM.Wminfo.Info.1"), "\n")),
               tklabel(WinInfo, text="\t "),
               sticky="w")

        tkgrid(tklabel(WinInfo, text="\t "),
               tklabel(WinInfo,
                       text=paste(mltext("sortiesLM.Wminfo.Info.2"),
                                  mltext("sortiesLM.Wminfo.Info.3"),
                                  "\n", sep="")),
               sticky="w")

        tkfocus(WinInfo)
        winSmartPlace.f(WinInfo)

        ## compMultiplesLM.f(objLM=objLM, Data=Data, factSpatial="protection.status", factTemp="year", resFile=resFile)
        compMultiplesLM.f(objLM=objLM, Data=Data, fact1=listFact[1], fact2=listFact[2],
                          resFile=resFile, exclude=factAna, Log=Log)

        ## Représentation des interactions :
        mainTitle <- graphTitle.f(metrique=metrique,
                                  modGraphSel=modSel, factGraph=factAna,
                                  listFact=listFact,
                                  model=mltext("sortiesLM.Graph.Title",
                                               language = getOption("P.lang")),
                                  type=type)

        eval(call(winFUN, pointsize=ifelse(isTRUE(getOption("P.graphPaper")), 14, 12)))
        par(mar=c(5, 4,
                  ifelse(isTRUE(getOption("P.graphPaper")), 2, 5),
                  2) + 0.1)
        with(Data,
             if (Log)                   # Les sens de variations peuvent changer en log (sur les moyennes) =>
                                        # besoin d'un graphique adapté :
         {
             eval(parse(text=paste("interaction.plot(", listFact[1], ", ", listFact[2],
                        ", log(", metrique, "), ylab=\"",
                        paste(mltext("sortiesLM.Graph.ylab.pfx",
                                     language = getOption("P.lang")),
                              "log(", Capitalize.f(varNames[metrique, "nom"]), ")",
                              mltext("sortiesLM.Graph.ylab.sfx",
                                     language = getOption("P.lang")), sep=""),
                        "\", xlab=\"", Capitalize.f(varNames[listFact[1], "nom"]),
                        "\", main=\"",
                        ifelse((! isTRUE(getOption("P.graphPaper"))) && isTRUE(getOption("P.title")), mainTitle, ""),
                        "\", trace.label=\"", Capitalize.f(varNames[listFact[2], "nom"]),
                        "\", cex.main=0.9)", sep="")))

         }else{
             eval(parse(text=paste("interaction.plot(", listFact[1], ", ", listFact[2],
                        ", ", metrique, ", ylab=\"",
                        paste(mltext("sortiesLM.Graph.ylab.pfx",
                                     language = getOption("P.lang")),
                              Capitalize.f(varNames[metrique, "nom"]),
                              mltext("sortiesLM.Graph.ylab.sfx",
                                     language = getOption("P.lang")),
                              switch(varNames[metrique, "genre"],
                                     "f"=, # Double the consonnant in French!
                                     "fp"=mltext("sortiesLM.Graph.ylab.sfxDblCons",
                                                       language = getOption("P.lang")),
                                     ""),
                              switch(varNames[metrique, "genre"],
                                     "f"=mltext("graphTitle.f",
                                                language = getOption("P.lang")),
                                     "fp"=mltext("graphTitle.fp",
                                                 language = getOption("P.lang")),
                                     "mp"=mltext("graphTitle.mp",
                                                 language = getOption("P.lang")),
                                     ""), # "moyen", moyens,
                                        # "moyenne" ou "moyennes" selon le genre.
                              sep=""),
                        "\", xlab=\"", Capitalize.f(varNames[listFact[1], "nom"]),
                        "\", main=\"", ifelse((! isTRUE(getOption("P.graphPaper"))) && isTRUE(getOption("P.title")),
                                              mainTitle, ""),
                        "\", trace.label=\"", Capitalize.f(varNames[listFact[2], "nom"]),
                        "\", cex.main=0.9)", sep="")))
         })

        tkdestroy(WinInfo)
    }else{
        if (length(listFact) == 1)
        {
            compSimplesLM.f(objLM=objLM, Data=Data, fact=listFact,
                            resFile=resFile, Log=Log)
        }else{}
    }

    subTitle <- graphTitle.f(metrique=metrique,
                             modGraphSel=modSel, factGraph=factAna,
                             listFact=listFact, model=modelType.f(objLM=objLM, Log=Log),
                             type=type)

    eval(call(winFUN, width=45, height=35))
    par(mfrow=c(2, 2), oma=c(0, 0, 4.7, 0))
    hist(objLM$residuals,
         xlab=mltext("sortiesLM.Graph.hist.xlab",
                     language = getOption("P.lang")),
         ylab= mltext("sortiesLM.Graph.hist.ylab",
                      language = getOption("P.lang")),
         main=NULL)
    mtext(mltext("sortiesLM.Graph.hist.title",
                 language = getOption("P.lang")), side=3, cex=0.8)

    ## Titre général :
    mtext(mltext("sortiesLM.Graph.diag.title",
                 language = getOption("P.lang")),
          side=3, outer=TRUE, line=3.4, cex=1.2)
    mtext(subTitle, side=3, outer=TRUE, line=-2.4, cex=1.1)

    ## Essayer glm.diag.plots('glm')...
    plot.lm.ml(objLM, which=2, cex.caption=0.8)
    plot.lm.ml(objLM, which=c(1, 4), cex.caption=0.8)

    ## ##################################################
    ## Sauvegarde des données :
    filename <- summary(resFile)$description

    close(resFile)                      # Maintenant seulement on peut fermer ce fichier.

    if (getOption("P.saveData") &&  ! isTRUE(sufixe == "(red)"))
    {
        writeData.f(filename=filename, Data=Data,
                    cols=NULL)
    }else{}

    ## Sauvegarde des infos sur les données et statistiques :
    if (getOption("P.saveStats") &&  ! isTRUE(sufixe == "(red)"))
    {
        infoStats.f(filename=filename, Data=Data, agregLevel=type, type="stat",
                    metrique=metrique, factGraph=factAna, factGraphSel=modSel,
                    listFact=listFact, listFactSel=listFactSel,
                    dataEnv=dataEnv, baseEnv=baseEnv)
    }else{}

    ## flush.console()
}


######################################### end of the function sortiesLM.f
