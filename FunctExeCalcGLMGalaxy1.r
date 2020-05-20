#Rscript 

#####################################################################################################################
#####################################################################################################################
################################# Compute a Generalized Linear Model from your data #################################
#####################################################################################################################
#####################################################################################################################

###################### Packages
#suppressMessages(library(MASS))
suppressMessages(library(multcomp))
suppressMessages(library(glmmTMB)) ###Version: 0.2.3

###################### Load arguments and declaring variables

args = commandArgs(trailingOnly=TRUE)
#options(encoding = "UTF-8")

if (length(args) < 11) {
    stop("At least one argument must be supplied, an input dataset file (.tabular).", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1

} else {
    Importdata <- args[1] ###### file name : metrics table
    ImportUnitobs <- args[2] ###### file name : unitobs informations
    colmetric <- as.numeric(args[3]) ###### Selected interest metric for GLM
    listFact <- strsplit(args [4],",")[[1]] ###### Selected response factors for GLM
    listRand <- strsplit(args [5],",")[[1]] ###### Selected randomized response factors for GLM
    FactAna <- args[6] ####### (optional) Selected splitting factors for GLMs
    Distrib <- args[7] ###### (optional) Selected distribution for GLM 
    log <- args[8] ###### (Optional) Log on interest metric ?
    aggreg <- args[9] ###### Aggregation level of the data table
    SupprOutlay <- args[10] ####### TRUE/FALSE : suppress oulayers ?
    source(args[11]) ###### Import functions

}
#### Data must be a dataframe with at least 3 variables : unitobs representing location and year ("observation.unit"), species code ("species.code") and abundance ("number")


#Import des données / Import data 
obs<- read.table(Importdata,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
obs[obs == -999] <- NA 
metric <- colnames(obs)[colmetric]
tabUnitobs <- read.table(ImportUnitobs,sep="\t",dec=".",header=TRUE,encoding="UTF-8")
tabUnitobs[tabUnitobs == -999] <- NA 
#factors <- fact.det.f(Obs=obs)

vars_data<-c("observation.unit","species.code","number")
err_msg_data<-"The input dataset doesn't have the right format. It need to have at least the following 3 variables :\n- observation.unit\n- species.code\n- number\n"
#check_file(obs,err_msg_data,vars_data,3)


####################################################################################################
########## Computing Generalized Linear Model ## Function : modeleLineaireWP2.unitobs.f ############
####################################################################################################

modeleLineaireWP2.unitobs.f <- function(metrique, listFact, listRand, FactAna, Distrib, log=FALSE, tabMetrics, tableMetrique, tabUnitobs, unitobs="observation.unit", outresiduals = FALSE, nbName="number")
{
    ## Purpose: Gestions des différentes étapes des modèles linéaires.
    ## ----------------------------------------------------------------------
    ## Arguments: metrique : la métrique choisie.
    ##            factAna : le facteur de séparation des graphiques.
    ##            factAnaSel : la sélection de modalités pour ce dernier
    ##            listFact : liste du (des) facteur(s) de regroupement
    ##            listFactSel : liste des modalités sélectionnées pour ce(s)
    ##                          dernier(s)
    ##            tabMetrics : table de métriques.
    ##            tableMetrique : nom de la table de métriques.
    ##            dataEnv : environnement de stockage des données.
    ##            baseEnv : environnement de l'interface.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 18 août 2010, 15:59

    tmpData <- tabMetrics

    if (listRand[1] != "None")
    {
        if (all(is.element(listFact,listRand)) || listFact[1] == "None")
        {
            RespFact <- paste("(1|",paste(listRand,collapse=") + (1|"),")")
            listFact <- listRand
        }else{
            listF <- listFact[!is.element(listFact,listRand)]
            RespFact <- paste(paste(listF, collapse=" + ")," + (1|",paste(listRand,collapse=") + (1|"),")")
            listFact <- c(listF,listRand)
        }   
    }else{

        RespFact <- paste(listFact, collapse=" + ")
    }

    ##Creating model's expression :

    #if (log == FALSE) {
        exprML <- eval(parse(text=paste(metrique, "~", RespFact)))
    #}else{
     #   exprML <- eval(parse(text=paste("ln(",metrique,")", "~", RespFact)))
    #}

    ##Creating analysis table :
    listFactTab <- c(listFact, FactAna)
    if(tableMetrique == "unit")
    {
        col <- c(unitobs,metrique)
        tmpData <- cbind(tmpData[,col], tabUnitobs[match(tmpData[,unitobs],tabUnitobs[,unitobs]),listFactTab])
        colnames(tmpData) <- c(col,listFactTab)

        for (i in listFactTab) {
            tmpData[,i] <- as.factor(tmpData[,i])
            #switch(i,
             #     "year" = {tmpData[,"year"] <- as.factor(tmpData[,"year"])},
              #    "site" = {tmpData[,"site"] <- as.factor(tmpData[,"site"])},
               #   "habitat1" = {tmpData[,"habitat1"] <- as.factor(tmpData[,"habitat1"])},
                #  "statut_protection" = {tmpData[,"statut_protection"] <- as.factor(tmpData[,"statut_protection"])},
                 # )
         }
    }else{
        stop("Warning")
    }

    ## Suppression des 'levels' non utilisés :
    tmpData <- dropLevels.f(tmpData)

    ## Aide au choix du type d'analyse :
    if (Distrib == "None") 
    {
        switch(class(tmpData[,metrique]),
              "integer"={
                         loiChoisie <- "poisson"
              },
              "numeric"={
                         loiChoisie <- "gaussian"
              },
              )
    }else{
        loiChoisie <- Distrib
    }

    Allfact <- c(metrique,listFact)
    ## Compute Model :    

    if (listRand[1] != "None")
    {
        res <- glmmTMB(exprML,family=loiChoisie, data=tmpData)
    }else{
        res <- glm(exprML,data=tmpData,family=loiChoisie)
    }
    

    ## Écriture des résultats formatés dans un fichier :
 
    sortiesLM.f(objLM=res, formule=exprML, metrique=metrique,
                #factAna=factAna, modSel=iFactGraphSel, listFactSel=listFactSel,
                listFact=listFact,
                Data=tmpData, #Log=Log,
                type=ifelse(tableMetrique == "unitSpSz" && factAna != "size.class",
                            "CL_unitobs",
                            "unitobs"))
    

}

################# Analysis

modeleLineaireWP2.unitobs.f(metrique=metric, listFact=listFact, listRand=listRand, FactAna=FactAna, Distrib=Distrib, log=log, tabMetrics=obs, tableMetrique=aggreg, tabUnitobs=tabUnitobs, outresiduals=SupprOutlay, nbName="number")
