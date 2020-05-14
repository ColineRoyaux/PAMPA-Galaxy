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
    ImportRefEsp <- args[2] ###### file name : species referential
    ImportUnitobs <- args[3] ###### file name : unitobs informations
    colmetric <- as.numeric(args[4]) ###### Selected interest metric for GLM
    listFact <- strsplit(args [5],",")[[1]] ###### Selected response factors for GLM
    listRand <- strsplit(args [6],",")[[1]] ###### Selected randomized response factors for GLM
    listFactSel <- args[7] ####### Selected modalities used in response factors for GLM
    factAna <- args[8] ###### Selected splitting factors for GLMs
    factAnaSel <- args[9] ###### Selected modalities used in splitting factors
    aggreg <- args[10] ###### Aggregation level of the data table
    SupprOutlay <- args[11] ####### TRUE/FALSE : suppress oulayers ?
    source(args[12]) ###### Import functions

}
#### Data must be a dataframe with at least 3 variables : unitobs representing location and year ("observation.unit"), species code ("species.code") and abundance ("number")


#Import des données / Import data 
obs<- read.table(Importdata,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
obs[obs == -999] <- NA 
metric <- colnames(obs)[colmetric]
tabUnitobs <- read.table(ImportUnitobs,sep="\t",dec=".",header=TRUE,encoding="UTF-8")
tabUnitobs[tabUnitobs == -999] <- NA 
refesp <- read.table(ImportRefEsp,sep="\t",dec=".",header=TRUE,encoding="UTF-8")
refesp[refesp == -999] <- NA 
#factors <- fact.det.f(Obs=obs)

vars_data<-c("observation.unit","species.code","number")
err_msg_data<-"The input dataset doesn't have the right format. It need to have at least the following 3 variables :\n- observation.unit\n- species.code\n- number\n"
#check_file(obs,err_msg_data,vars_data,3)


####################################################################################################
########## Computing Generalized Linear Model ## Function : modeleLineaireWP2.unitobs.f ############
####################################################################################################

modeleLineaireWP2.unitobs.f <- function(metrique, listFact, listRand,tabMetrics, tableMetrique, tabUnitobs, unitobs="observation.unit", refesp, outresiduals = FALSE, nbName="number")
                                        #,dataEnv, baseEnv=.GlobalEnv)
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

    if (length(listRand) != 0 && listRand[1] != "")
    {
        if (all(is.element(listFact,listRand)))
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
    exprML <- eval(parse(text=paste(metrique, "~", RespFact)))

    ##Creating analysis table :
    if(is.element("factor(year)",listFact))
    {
        listFact[[grep("factor",listFact)]]<-"year"
    }else{}

    if(tableMetrique == "unit")
    {
        col <- c(unitobs,metrique)
        tmpData <- cbind(tmpData[,col], tabUnitobs[match(tmpData[,unitobs],tabUnitobs[,unitobs]),listFact])
        colnames(tmpData) <- c(col,listFact)
    }else{
        stop("Warning")
    }

    ## Suppression des 'levels' non utilisés :
    tmpData <- dropLevels.f(tmpData)

    ## Aide au choix du type d'analyse :
    if (metrique == "pres.abs")
    {
        loiChoisie <- "binomial"
    }else{
        switch(class(tmpData[,metrique]),
               "integer"={
                          loiChoisie <- "poisson"
               },
               "numeric"={
                          loiChoisie <- "gaussian"
               },
               )
    }

    Allfact <- c(metrique,listFact)
    ## Compute Model :

    if (length(listRand) != 0 && listRand[1] != "")
    {
        res <- glmmTMB(exprML,family=loiChoisie, data=tmpData)
    }else{
        res <- glm(exprML,data=tmpData,family=loiChoisie)
    }


        ## Écriture des résultats formatés dans un fichier :
        #tryCatch(
                 sortiesLM.f(objLM=res, formule=exprML, metrique=metrique,
                             #factAna=factAna, modSel=iFactGraphSel, listFactSel=listFactSel,
                             listFact=listFact,
                             Data=tmpData, #Log=Log,
                             type=ifelse(tableMetrique == "unitSpSz" && factAna != "size.class",
                                         "CL_unitobs",
                                         "unitobs"))
         #        ,error=errorLog.f)
    #cat(summary(res),file="Errors.txt")

    #return(sum)

}

################# Analysis

modeleLineaireWP2.unitobs.f(metrique=metric, listFact=listFact, listRand=listRand, tabMetrics=obs, tableMetrique=aggreg, tabUnitobs=tabUnitobs, refesp=refesp, outresiduals=SupprOutlay, nbName="number")
#filename <- "Errors.txt"
#write.table(res, filename, row.names=FALSE, sep="\t", dec=".",fileEncoding="UTF-8")
