#Rscript 

#####################################################################################################################
#####################################################################################################################
###################################### Create a plot from your community data #######################################
#####################################################################################################################
#####################################################################################################################

###################### Packages
suppressMessages(library(ggplot2))


###################### Load arguments and declaring variables

args = commandArgs(trailingOnly=TRUE)
#options(encoding = "UTF-8")

if (length(args) < 10) {
    stop("At least 4 arguments must be supplied : \n- two input dataset files (.tabular) : metrics table and unitobs table \n- Interest variable field from metrics table \n- Response variable from unitobs table.", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1

} else {
    Importdata <- args[1] ###### file name : metrics table
    ImportUnitobs <- args[2] ###### file name : unitobs informations
    colmetric <- as.numeric(args[3]) ###### Selected interest metric for plot (y)
    colFactP <- as.numeric(args[4]) ###### Selected primary response factor for plot (x)
    colFactS <- args[5] ###### (optional) Selected secondary response factor for plot (col)
    plotType <- args[6] ####### Selected plot type
    colFactAna <- args[7] ####### (optional) Selected splitting factors for GLMs
    log <- args[8] ###### (Optional) Log on interest metric ?
    SupprOutlay <- args[9] ####### TRUE/FALSE : suppress oulayers ?
    source(args[10]) ###### Import functions

}
#### Data must be a dataframe with at least 3 variables : unitobs representing location and year ("observation.unit"), species code ("species.code") and abundance ("number")


#Import des données / Import data 
obs<- read.table(Importdata,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
obs[obs == -999] <- NA 
metric <- colnames(obs)[colmetric]
tabUnitobs <- read.table(ImportUnitobs,sep="\t",dec=".",header=TRUE,encoding="UTF-8")
tabUnitobs[tabUnitobs == -999] <- NA 
FactP <- colnames(tabUnitobs)[colFactP]

if (colFactS != "None")
{
    FactS <- colnames(tabUnitobs)[as.numeric(colFactS)]
    if (class(tabUnitobs[FactS]) == "numeric" || FactS == "observation.unit"){stop("Wrong chosen separation factor : Analysis can't be separated by observation unit or numeric factor")}
}else{
    FactS <- colFactS
}
FactS <- colnames(tabUnitobs)[colFactS]

if (colFactAna != "None")
{
    FactAna <- colnames(tabUnitobs)[as.numeric(colFactAna)]
    if (class(tabUnitobs[FactAna]) == "numeric" || FactAna == "observation.unit"){stop("Wrong chosen separation factor : Analysis can't be separated by observation unit or numeric factor")}
}else{
    FactAna <- colFactAna
}

vars_data1<- NULL
err_msg_data1<-"The input metrics dataset doesn't have the right format. It needs to have at least the following 2 variables :\n- observation.unit (or year and site)\n- numeric or integer metric\n"
check_file(obs,err_msg_data1,vars_data1,2)

vars_data2 <- NULL
err_msg_data2<-"The input unitobs dataset doesn't have the right format. It needs to have at least the following 2 variables :\n- observation.unit (or year and site)\n- factors used in plot \n"
check_file(tabUnitobs,err_msg_data2,vars_data2,2)

####################################################################################################
######################### Creating plot community ## Function : CommPlot.f #########################
####################################################################################################

CommPlot.f <- function(metriqueY, FactX, FactCol, FactAna, tabMetrics, tabUnitobs, plotType, unitobs="observation.unit", nbName="number")
{
    ## Purpose: Creating plot with community metrics
    ## ----------------------------------------------------------------------
    ## Arguments: metriqueY : chosen metric (plot Y)
    ##            FactX : Chosen X factor (plot X)
    ##            FactCol : Chosen colours factor (plot col)
    ##            FactAna : Chosen splitting factor for plural plots
    ##            tabMetrics : metrics table
    ##            tabUnitobs : unitobs table
    ##            plotType : type of plot wanted
    ##            unitobs : name of field unitobs (primary key linking input tables)
    ##            nbName : name of abundance field
    ## ----------------------------------------------------------------------
    ## Author: Coline ROYAUX, Date: 02 june 2020, 15:59

    tmpData <- tabMetrics

    ##Creating analysis table :
    listFactTab <- c(FactX,FactCol,FactAna)
    listFactTab <- na.omit(listFactTab[listFactTab != "None"])

    if (all(is.na(match(tmpData[,unitobs],tabUnitobs[,unitobs])))) {stop("Observation units doesn't match in the two input tables")}

    if(! is.element("species.code",colnames(tmpData)))
    {
        col <- c(unitobs,metriqueY)
        tmpData <- cbind(tmpData[,col], tabUnitobs[match(tmpData[,unitobs],tabUnitobs[,unitobs]),listFactTab])
        colnames(tmpData) <- c(col,listFactTab)

        if (FactX == "year"){tmpData[,"year"] <- as.integer(tmpData[,"year"])}
    }else{
        stop("Warning : wrong data frame, data frame should be aggregated by observation unit (year and site)")
    }

    ## Suppression des 'levels' non utilisés :
    tmpData <- dropLevels.f(tmpData)

    ##Creating main ggplot expression :
    gg <- ggplot(data=tmpData)
    #stop(colnames(tmpData))
    if (!is.na(FactS))
    {
        switch(plotType,
               "plot"={gg <- gg + geom_point(aes_string(x=FactX,y=metriqueY,colour=FactS)) + geom_line(aes_string(x=FactX,y=metriqueY,colour=FactS))},
               "boxplot"={gg <- gg + geom_boxplot(aes_string(x=FactX,y=metriqueY,colour=FactS))},
               "barplot"={gg <- gg + geom_bar(stat='identity', aes_string(colour=FactS,fill=FactS))},
              )
    }else{
        switch(plotType,
               "plot"={gg <- gg + geom_point(colour='darkblue') + geom_line(colour='darkblue')},
               "boxplot"={gg <- gg + geom_boxplot(colour='darkblue')},
               "barplot"={gg <- gg + geom_bar(stat='identity', colour='darkblue',fill='darkblue')},
              )
    }

    if (FactX == "year" || plotType == "plot") {gg <- gg + scale_x_continuous("time",breaks=unique(tmpData[,FactX]))}

    ggsave("plot.png",gg,width=15,height=9,units="cm")
    print("plot.png")

    ## Create plot(s) :
    
#    for (sp in levels(tmpData[,FactAna])) 
 #   {
  #      cutData <- tmpData[grep(sp,tmpData[,FactAna]),]
   #     cutData <- dropLevels.f(cutData)
#
 #       cat("--------------------------------------------------------------------------------\n",
  #          "--------------------------------------------------------------------------------\n Analysis for species ",sp,
   #         " :\n--------------------------------------------------------------------------------\n--------------------------------------------------------------------------------\n",
    #        sep="",file=resFile,append=TRUE)
#
 #       res <-""
#
 #       if (listRand[1] != "None")
  #      {
   #         res <- tryCatch(glmmTMB(exprML,family=loiChoisie, data=cutData), error=function(e){})
    #    }else{
     #       res <- tryCatch(glm(exprML,data=cutData,family=loiChoisie), error=function(e){})
      #  }
#
          ## Écriture des résultats formatés dans un fichier :
 #        if (! is.null(res))
  #       {
   #         sortiesLM.f(objLM=res, formule=exprML, metrique=metrique,
    #                    factAna=factAna, #modSel=iFactGraphSel, listFactSel=listFactSel,
     #                   listFact=listFact,
      #                  Data=cutData, #Log=Log,
       #                 type=ifelse(tableMetrique == "unitSpSz" && factAna != "size.class",
        #                            "CL_unitobs",
         #                           "unitobs"))
#
 #       }else{
  #          cat("\nCannot compute GLM. Check if one or more factor(s) have only one level, or try with another distribution for the model in advanced settings \n\n",file=resFile,append=TRUE)
   #     }

    #}
#
}

################# Analysis


CommPlot.f(metriqueY=metric, FactX=FactP, FactCol=FactS, FactAna=FactAna, tabMetrics=obs, tabUnitobs=tabUnitobs, plotType=plotType)
