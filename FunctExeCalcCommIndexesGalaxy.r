#Rscript 

#####################################################################################################################
#####################################################################################################################
################################# Calculate community indexes from observation data #################################
#####################################################################################################################
#####################################################################################################################

###################### Packages R base

###################### Load arguments and declaring variables

args = commandArgs(trailingOnly=TRUE)
#options(encoding = "UTF-8")

if (length(args) < 4) {
    stop("At least one argument must be supplied, an input dataset file (.tabular).", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1

} else {
    Importdata<-args[1] ###### Nom du fichier importé avec son extension / file name imported with the file type ".filetype"  
    index <- args[2] ###### List of selected metrics to calculate
    source(args[3]) ###### Import functions

}
#### Data must be a dataframe with at least 3 variables : unitobs representing location and year ("observation.unit"), species code ("species.code") and abundance ("number")


#Import des données / Import data 
obs<- read.table(Importdata,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
obs[obs == -999] <- NA 
factors <- fact.det.f(Obs=obs)

ObsType <- def.typeobs.f(Obs=obs)

vars_data<-c("observation.unit","species.code","number")
err_msg_data<-"The input dataset doesn't have the right format. It need to have at least the following 3 variables :\n- observation.unit\n- species.code\n- number\n"
check_file(obs,err_msg_data,vars_data,3)



####################################################################################################
################## create community metrics table ## Function : calcBiodiv.f #######################
####################################################################################################

########################################################################################################################
calcBiodiv.f <- function(Data, 
                         #refesp, 
                         MPA, unitobs="observation.unit", code.especes="species.code", nombres="number",
                         indices=index, global=FALSE, printInfo=FALSE
                         #, dataEnv=.GlobalEnv
                         )
{
    ## Purpose: calcul des indices de biodiversité
    ## ----------------------------------------------------------------------
    ## Arguments: Data : les données à partir desquelles calculer les
    ##                   indices. Doivent comporter au minimum (colones) :
    ##                     * unités d'observations/sites
    ##                     * espèces présentes
    ##                     * nombre d'individus /espèce/unitobs.
    ##            refesp : le référentiel espèces.
    ##            MPA : l'AMP (chaîne de charactères).
    ##            unitobs : nom de la colone d'unités d'observation.
    ##            code.especes : nom de la colone d'espèces.
    ##            nombres : nom de la colone de nombres.
    ##            indices : liste des indices à calculer
    ##                      (vecteur de caractères)
    ##            global : est-ce que les résultats doivent être exportés
    ##                     globalement (booléen).
    ##            printInfo : affichage des infos (chargement) ? (booléen).
    ##            dataEnv : environnement des données
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 29 oct. 2010, 08:58

    ## Unitobs appartenant a l'AMP courante:
    #unitobsData <- get("unitobs", envir=dataEnv)

    #Data <- subset(Data,
                   #is.element(Data[ , unitobs],
                              #unitobsData[unitobsData[ , getOption("P.MPAfield")] == MPA ,
                                          #unitobs]))

    #DataTmp <- Data

    ## Supression de tout ce qui n'a pas d'espèce précisee (peut être du non biotique ou identification >= genre) :

    notspline <- grep("(sp\\.)$|([1-9])$|^(Absencemacrofaune)$|^(NoID)$|^(Acrobranc)$|^(Acrodigit)$|^(Acroencr)$|^(Acrosubm)$|^(Acrotabu)$|^(Adredure)$|^(Adremoll)$|^(Algaturf)$|^(Balimona)$|^(Corablan)$|^(CoradurV)$|^(Coraenal)$|^(Coramor1)$|^(Coramor2)$|^(Coramou)$|^( Dallcora)$|^(Debrcora)$|^(Debris)$|^(Hare)$|^(HexaChar)$|^(MuraCong)$|^(Nacrbran)$|^(Nacrcham)$|^(Nacrencr)$|^(Nacrfoli)$|^(Nacrmass)$|^(Nacrsubm)$|^(Recrcora)$|^(Roche)$|^(Sable)$|^(Vase)$",Data[, code.especes], value=FALSE)
    Data <- Data[-notspline, ]
    
    #NotSpecies <- c("Absencemacrofaune","NoID","Acrobranc","Acrodigit","Acroencr","Acrosubm","Acrotabu","Adredure","Adremoll","Algaturf","Balimona","Corablan","CoradurV","Coraenal","Coramou"," Dallcora","Debrcora","Debris","Hare","HexaChar","MuraCong","Nacrbran","Nacrcham","Nacrencr","Nacrfoli","Nacrmass","Nacrsubm","Recrcora","Roche","Sable","Vase")
    

    #if (! nrow(Data <- Data[(spTmp <- refesp$species[match(Data[ , code.especes], refesp$species.code)]) != "sp." &
     #                       !is.na(spTmp), ]))
    #{

     #   return(Data)
    #}else{}

    ## Suppression des niveaux de facteur inutilisés :
    Data <- dropLevels.f(df=Data)


    ## Si les données ne sont pas encore agrégées /espèce/unitobs on le fait ici :
    if (nrow(Data) > nrow(expand.grid(unique(Data[ , unitobs]), unique(Data[ , code.especes]))))
    {
        Data <- agregations.generic.f(Data=Data, metrics=nombres,
                                      factors=c(unitobs, code.especes),
                                      listFact=NULL)#, dataEnv=dataEnv)
    }else{}

    df.biodiv <- as.data.frame(as.table(tapply(Data[ , nombres],
                                               Data[ , unitobs],
                                               sum, na.rm=TRUE)))

    colnames(df.biodiv) <- c(unitobs, nombres)

## ##################################################
    ## Richesse spécifique :
    Data$pres.abs <- presAbs.f(nombres=Data[ , nombres], logical = FALSE)

    df.biodiv$species.richness <- as.vector(tapply(Data$pres.abs,
                                                   Data[ , unitobs], sum, na.rm=TRUE),
                                            "integer")
    ## ... as.vector to avoid the class "array".

 ## ##################################################
    ## Indices de Simpson et Shannon et dérivés :

    matNombres <- tapply(Data[ , nombres], # Matrice de nombres d'individus /espèce/unitobs.
                         list(Data[ , unitobs], Data[ , code.especes]),
                         sum, na.rm=TRUE)

    matNombres[is.na(matNombres)] <- 0  # Vrais zéros

    ## Proportion d'individus de chaque espèce dans l'unitobs :
    propIndiv <- sweep(matNombres, 1,                           #
                       apply(matNombres, 1, sum, na.rm = TRUE), # Nombre d'individus / unitobs ; équiv df.biodiv$nombre.
                       FUN="/")

    ## Indices de Simpson :
    df.biodiv$simpson <- apply(propIndiv^2, 1, sum, na.rm=TRUE)

    if (any(is.element(c("all", "simpson.l"), indices)))
    {
        df.biodiv$simpson.l <- 1 - df.biodiv$simpson
    }

    ## calcul de l'indice de Shannon :
    df.biodiv$shannon <- -1 * apply(propIndiv * log(propIndiv), 1, sum, na.rm=TRUE)

    ## calcul de l'indice de Pielou :
    if (any(is.element(c("all", "pielou"), indices)))
    {
        df.biodiv$pielou <- df.biodiv$shannon / log(df.biodiv$species.richness)
    }

    ## calcul de l'indice de Hill :
    if (any(is.element(c("all", "hill"), indices)))
    {
        df.biodiv$hill <- (1 - df.biodiv$simpson) / exp(df.biodiv$shannon)
                                        # équiv df.biodiv$l.simpson / exp(df.biodiv$shannon)
    }

    ## suppression de l'indice de shannon (non pertinent)
    #df.biodiv$shannon <- NULL

    ## ## On retablit les niveaux de facteurs:
    ## colFact <- colnames(df.biodiv)[is.element(sapply(df.biodiv, class), "factor")]

    ## for (col in colFact)
    ## {
    ##     levels(df.biodiv[ , colFact]) <- levels(DataTmp[ , colFact])
    ## }


    return(df.biodiv)
}

################# Analysis

res <- calc.numbers.f(obs, ObsType=ObsType , factors=factors, nbName="number")

tableCommunityIndexes <- calcBiodiv.f(res, #refesp, 
                         MPA, unitobs="observation.unit", code.especes="species.code", nombres="number",
                         indices=index, global=FALSE, printInfo=FALSE #, dataEnv=.GlobalEnv
                         )
#Save dataframe in a tabular format

filenameComm <- "TabCommunityIndexes.tabular"
write.table(tableCommunityIndexes, filenameComm, row.names=FALSE, sep="\t", dec=".",fileEncoding="UTF-8")
cat(paste("\nWrite table with Community indexes. \n--> \"",filenameComm,"\"\n",sep=""))

