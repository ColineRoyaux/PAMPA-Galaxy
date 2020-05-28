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
    stop("At least 4 arguments must be supplied : \n- two input dataset files (.tabular) : metrics table and unitobs table \n- Interest variable field from metrics table \n- Response variable from unitobs table.", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1

} else {
    Importdata <- args[1] ###### file name : metrics table
    ImportUnitobs <- args[2] ###### file name : unitobs informations
    colmetric <- as.numeric(args[3]) ###### Selected interest metric for GLM
    listFact <- strsplit(args [4],",")[[1]] ###### Selected response factors for GLM
    listRand <- strsplit(args [5],",")[[1]] ###### Selected randomized response factors for GLM
    colFactAna <- args[6] ####### (optional) Selected splitting factors for GLMs
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

if (colFactAna != "None")
{
    FactAna <- colnames(tabUnitobs)[as.numeric(colFactAna)]
    if (class(tabUnitobs[FactAna]) == "numeric" || FactAna == "observation.unit"){stop("Wrong chosen separation factor : Analysis can't be separated by observation unit or numeric factor")}
}else{
    FactAna <- colFactAna
}


#factors <- fact.det.f(Obs=obs)

vars_data1<- NULL
err_msg_data1<-"The input metrics dataset doesn't have the right format. It needs to have at least the following 2 variables :\n- observation.unit (or year and site)\n- numeric or integer metric\n"
check_file(obs,err_msg_data1,vars_data1,2)

vars_data2 <- c(listFact,listRand)
err_msg_data2<-"The input unitobs dataset doesn't have the right format. It needs to have at least the following 2 variables :\n- observation.unit (or year and site)\n- factors used in GLM (habitat1, year and/or site)\n"
check_file(tabUnitobs,err_msg_data2,vars_data2,2)

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
     #   exprML <- eval(parse(text=paste("log(",metrique,")", "~", RespFact)))
    #}

    ##Creating analysis table :
    listFactTab <- c(listFact, FactAna)
    listFactTab <- listFactTab[listFactTab != "None"]

    if (all(is.na(match(tmpData[,unitobs],tabUnitobs[,unitobs])))) {stop("Observation units doesn't match in the two input tables")}

    if(! is.element("species.code",colnames(tmpData)))
    {
        col <- c(unitobs,metrique)
        tmpData <- cbind(tmpData[,col], tabUnitobs[match(tmpData[,unitobs],tabUnitobs[,unitobs]),listFactTab])
        colnames(tmpData) <- c(col,listFactTab)

        for (i in listFactTab) {
            tmpData[,i] <- as.factor(tmpData[,i])
         }
    }else{
        stop("Warning : wrong data frame, data frame should be aggregated by observation unit (year and site)")
    }

    ## Suppression des 'levels' non utilisés :
    tmpData <- dropLevels.f(tmpData)

    ## Aide au choix du type d'analyse :
    if (Distrib == "None") 
    {
        switch(class(tmpData[,metrique]),
              "integer"={loiChoisie <- "poisson"},
              "numeric"={loiChoisie <- "gaussian"},
              stop("Selected metric class doesn't fit, you should select an integer or a numeric variable"))
    }else{
        loiChoisie <- Distrib
    }

    ## Compute Model(s) :
    if (FactAna != "None" && nlevels(tmpData[,FactAna]) > 1)
    {
        Anacut <- levels(tmpData[,FactAna])
    }else{
        Anacut <- NULL
    }

    resFile <- "GLMSummary.txt"
    
    for (cut in Anacut) 
    {
        cutData <- tmpData[grep(cut,tmpData[,FactAna]),]
        cutData <- dropLevels.f(cutData)

        cat("--------------------------------------------------------------------------------\n",
            "--------------------------------------------------------------------------------\n Analysis for level ",cut,
            " :\n--------------------------------------------------------------------------------\n--------------------------------------------------------------------------------\n",
            sep="",file=resFile,append=TRUE)

        res <-""

        if (listRand[1] != "None")
        {
            res <- tryCatch(glmmTMB(exprML,family=loiChoisie, data=cutData), error=function(e){})
        }else{
            res <- tryCatch(glm(exprML,data=cutData,family=loiChoisie), error=function(e){})
        }

          ## Écriture des résultats formatés dans un fichier :
         if (! is.null(res))
         {
            sortiesLM.f(objLM=res, formule=exprML, metrique=metrique,
                        factAna=factAna, #modSel=iFactGraphSel, listFactSel=listFactSel,
                        listFact=listFact,
                        Data=cutData, #Log=Log,
                        type=ifelse(tableMetrique == "unitSpSz" && factAna != "size.class",
                                    "CL_unitobs",
                                    "unitobs"))

        }else{
            cat("\nOne or more factor(s) have only one level \n\n",file=resFile,append=TRUE)
        }

    }

    ## Global analysis : 


    cat("--------------------------------------------------------------------------------\n",
        "--------------------------------------------------------------------------------\n Global analysis",
        " :\n--------------------------------------------------------------------------------\n--------------------------------------------------------------------------------\n",
        sep="",file=resFile,append=TRUE)

    if (listRand[1] != "None")
    {
        resG <- glmmTMB(exprML,family=loiChoisie, data=tmpData)
    }else{
        resG <- glm(exprML,data=tmpData,family=loiChoisie)
    }

    ## Écriture des résultats formatés dans un fichier :
    sortiesLM.f(objLM=resG, formule=exprML, metrique=metrique,
                #factAna=factAna, modSel=iFactGraphSel, listFactSel=listFactSel,
                listFact=listFact,
                Data=tmpData, #Log=Log,
                type=ifelse(tableMetrique == "unitSpSz" && factAna != "size.class",
                            "CL_unitobs",
                            "unitobs"))

}

################# Analysis

modeleLineaireWP2.unitobs.f(metrique=metric, listFact=listFact, listRand=listRand, FactAna=FactAna, Distrib=Distrib, log=log, tabMetrics=obs, tableMetrique=aggreg, tabUnitobs=tabUnitobs, outresiduals=SupprOutlay, nbName="number")
