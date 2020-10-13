#Rscript 

#####################################################################################################################
#####################################################################################################################
###################################### Create a plot from your community data #######################################
#####################################################################################################################
#####################################################################################################################

###################### Packages
suppressMessages(library(ggplot2))
suppressMessages(library(boot))

###################### Load arguments and declaring variables

args = commandArgs(trailingOnly=TRUE)
#options(encoding = "UTF-8")

if (length(args) < 2) {
    stop("At least 3 arguments must be supplied input dataset file with GLM results (.tabular)", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1

} else {
    Importdata <- args[1] ###### file name : glm results table
    DataTab <- args[2]
    UnitobsTab <- args[3]
    source(args[4]) ###### Import functions

}
#### Data must be a dataframe with at least 3 variables : unitobs representing location and year ("observation.unit"), species code ("species.code") and abundance ("number")


#Import des données / Import data 
glmtable <- read.table(Importdata,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
datatable <- read.table(DataTab,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #
unitobs <- read.table(UnitobsTab,sep="\t",dec=".",header=TRUE,encoding="UTF-8") #

vars_data1 <- c("analysis","Interest.var","distribution")
err_msg_data1<-"The input GLM results dataset doesn't have the right format. It needs to have at least the following 3 variables :\n- analysis\n- Interest.var\n- distribution\n"
check_file(glmtable,err_msg_data1,vars_data1,4)

if (length(grep("[0-2][0|9][0-9][0-9].Estimate",colnames(glmtable))) == 0){stop("The input GLM results dataset doesn't have the right format or informations. This tool is to represent temporal trends, if your GLM doesn't take the year variable as a fixed effect this tool is not proper to make any representation of it. It needs to have at least estimates of every year from your time series GLM as columns with name formated as : yyyy Estimate (example : 2020 Estimate).")}



####################################################################################################
######################### Creating plot community ## Function : CommPlot.f #########################
####################################################################################################


       

############################################################################################################ fonction graphique appelée par main.glm / function called by main.glm for graphical output
ggplot.glm <- function(glmtable, datatable,unitobs,serie=NULL,sp,description=TRUE,
                          tendanceSurFigure=TRUE,seuilOccu=14,assessIC=TRUE) 
{

    seuilSignif <- 0.05
    metric <- as.character(glmtable[1,"Interest.var"])
    distrib <- as.character(glmtable[1,"distribution"])
    col <- c("observation.unit","location",metric)

    if (colnames(glmtable)[length(glmtable)]=="separation")
    {
        cut <- as.character(glmtable[1,"separation"])
        if(cut != "None")
        {
            datatable <- cbind(datatable[,col], unitobs[match(datatable[,"observation.unit"],unitobs[,"observation.unit"]),c("year",cut)])
            colnames(datatable) <- c(col,"year",cut)
        }else{
            datatable <- cbind(datatable[,col], unitobs[match(datatable[,"observation.unit"],unitobs[,"observation.unit"]),"year"])
            colnames(datatable) <- c(col,"year")
        }

    }else{
        cut <- "species.code"
        col <- c(col,cut)
        datatable <- cbind(datatable[,col], unitobs[match(datatable[,"observation.unit"],unitobs[,"observation.unit"]),"year"])
        colnames(datatable) <- c(col,"year")
    }

    ##vpan vecteur des panels de la figure  ###### POUR FAIRE LES GRAPHIQUES
    switch(as.character(metric),
           "number" = vpan <- c("Abundance variation","Raw abundance"),
           "pres.abs" = vpan <- c("Presence-absence variation","% presence in location"),
           vpan <- c(paste(metric," variation"),paste("Mean ", metric)))

    ##Cut table for 1 species
    glmtab <- glmtable[glmtable[,1]==sp,]

    glmtab <- glmtab[,grep("FALSE",is.na(glmtab[1,]))]
    ## specifications des variables temporelles necesaires pour les analyses / specification of temporal variable necessary for the analyses
    an <- as.numeric(unlist(strsplit(gsub("X","",paste(colnames(glmtab)[grep("[0-2][0|9][0-9][0-9].Estimate",colnames(glmtab))],collapse=" ")),split=".Estimate")))

    year <- sort(c(min(an)-1,an))
    nbans <- length(year)
    pasdetemps <- nbans-1
    firstY <- min(year)
    lastY <- max(year)
	
    nbSp <- length(glmtable[,1])

    coefan <- glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].Estimate",colnames(glmtab))] ## tendences sur la periode = coefficient regression de variable year  / tendency of
    coefan <- unlist(coefan[grep("FALSE",is.na(coefan))])
    erreuran <- glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].Std.Err",colnames(glmtab))]
    erreuran <- unlist(erreuran[grep("FALSE",is.na(erreuran))]) #### erreur standard sur le coefficient de regression de la variable year  / standard error on the regression coefficient of the year 

    switch(distrib,
           "poisson"={coefyear <- c(1,exp(as.numeric(coefan)))
                      erreuryear1 <- c(0,as.numeric(erreuran)*exp(as.numeric(coefan)))
                      if(assessIC) 
                      {
                          ic_inf_sim <- c(1,exp(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])))
                          ic_sup_sim <- c(1,exp(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])))
                      } else {
                          ic_inf_sim <- NA
                          ic_sup_sim <- NA
                      }},
           "quasipoisson"={coefyear <- c(1,exp(as.numeric(coefan)))
                           erreuryear1 <- c(0,as.numeric(erreuran)*exp(as.numeric(coefan)))
                           if(assessIC) 
                           {
                               ic_inf_sim <- c(1,exp(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])))
                               ic_sup_sim <- c(1,exp(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])))
                           } else {
                               ic_inf_sim <- NA
                               ic_sup_sim <- NA
                           }},
           "inverse.gaussian"={coefyear <- c(1,as.numeric(coefan)^(-1/2))
                               erreuryear1 <- c(0,as.numeric(erreuran)*(as.numeric(coefan)^(-1/2)))
                               if(assessIC) 
                               {
                                   ic_inf_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])^(-1/2))
                                   ic_sup_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])^(-1/2))
                               } else {
                                   ic_inf_sim <- NA
                                   ic_sup_sim <- NA
                               }},
           "binomial"={coefyear <- c(1,inv.logit(as.numeric(coefan)))
                       erreuryear1 <- c(0,as.numeric(erreuran)*inv.logit(as.numeric(coefan)))
                       if(assessIC) 
                       {
                           ic_inf_sim <- c(1,inv.logit(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])))
                           ic_sup_sim <- c(1,inv.logit(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])))
                       } else {
                           ic_inf_sim <- NA
                           ic_sup_sim <- NA
                       }},
           "quasibinomial"={coefyear <- c(1,inv.logit(as.numeric(coefan)))
                            erreuryear1 <- c(0,as.numeric(erreuran)*inv.logit(as.numeric(coefan)))
                            if(assessIC) 
                            {
                                ic_inf_sim <- c(1,inv.logit(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])))
                                ic_sup_sim <- c(1,inv.logit(as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])))
                            } else {
                                ic_inf_sim <- NA
                                ic_sup_sim <- NA
                            }},
           "Gamma"={coefyear <- c(1,as.numeric(coefan)^(-1))
                    erreuryear1 <- c(0,as.numeric(erreuran)*(as.numeric(coefan)^(-1)))
                    if(assessIC) 
                    {
                        ic_inf_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))])^(-1))
                        ic_sup_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))])^(-1))
                    } else {
                        ic_inf_sim <- NA
                        ic_sup_sim <- NA
                    }},
           {coefyear <- c(1,as.numeric(coefan))
            erreuryear1 <-  c(0,as.numeric(erreuran)*as.numeric(coefan))## erreur standard par année / the standard error per year
            if(assessIC) 
            {
                ic_inf_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_inf",colnames(glmtab))]))
                ic_sup_sim <- c(1,as.numeric(glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].IC_up",colnames(glmtab))]))
            } else {
                ic_inf_sim <- NA
                ic_sup_sim <- NA
            }})

    pval <- glmtab[glmtab[,1]==sp,grep("[0-2][0|9][0-9][0-9].Pvalue",colnames(glmtab))]
    pval <- c(1,unlist(pval[grep("FALSE",is.na(pval))])) ###### p value

    tab1 <- data.frame(year,val=coefyear,  ## tab1 table pour la realisation des figures / table for the graphical outputs
                       LL=unlist(ic_inf_sim),UL=unlist(ic_sup_sim),
                       catPoint=ifelse(pval<seuilSignif,"significatif",NA),pval,
                       courbe=vpan[1],
                       panel=vpan[1])
    ## netoyage des intervalle de confiance mal estimés et qd donnees pas suffisantes pour calcul d'IC /cleaning of wrong or biaised measures of the confidence interval
    if(assessIC) 
    {
        tab1$UL <-  ifelse(tab1$UL == Inf, NA,tab1$UL)
        tab1$UL <-  ifelse(tab1$UL > 1.000000e+20, NA,tab1$UL)
        tab1$UL[1] <- 1
        tab1$val <-  ifelse(tab1$val > 1.000000e+20,1.000000e+20,tab1$val)
    }

    coefancontinu <- as.numeric(as.character(glmtab[glmtab[,1]==sp,grep("year.Estimate",colnames(glmtab))])) ## tendences sur la periode = coefficient regression de variable year  / tendency of population evolution on the studied period = regression coefficient of the variable year 
    switch(distrib,
           "poisson"={trend <- round(exp(as.numeric(coefancontinu)),3) 
                      pourcentage <- round((exp(as.numeric(coefancontinu)*as.numeric(pasdetemps))-1)*100,2)},
           "quasipoisson"={trend <- round(exp(as.numeric(coefancontinu)),3)
                           pourcentage <- round((exp(as.numeric(coefancontinu)*as.numeric(pasdetemps))-1)*100,2)},
           "inverse.gaussian"={trend <- round(as.numeric(coefancontinu)^(-1/2),3)
                               pourcentage <- round((((as.numeric(coefancontinu)*as.numeric(pasdetemps))^(-1/2))-1)*100,2)},
           "binomial"={trend <- round(inv.logit(as.numeric(coefancontinu)),3)
                       pourcentage <- round((inv.logit(as.numeric(coefancontinu)*as.numeric(pasdetemps))-1)*100,2)},
           "quasibinomial"={trend <- round(inv.logit(as.numeric(coefancontinu)),3)
                            pourcentage <- round((inv.logit(as.numeric(coefancontinu)*as.numeric(pasdetemps))-1)*100,2)},
           "Gamma"={trend <- round(as.numeric(coefancontinu)^(-1),3)
                    pourcentage <- round((((as.numeric(coefancontinu)*as.numeric(pasdetemps))^(-1))-1)*100,2)},
           {trend <- round(as.numeric(coefancontinu),3)
            pourcentage <-  round((((as.numeric(coefancontinu)*as.numeric(pasdetemps)))-1)*100,2)})
        
    pval <- as.numeric(as.character(glmtab[glmtab[,1]==sp,grep("year.Pvalue",colnames(glmtab))]))
        
    erreuran <- as.numeric(as.character(glmtab[glmtab[,1]==sp,grep("year.Std.Err",colnames(glmtab))])) #### récuperer l'erreur standard / retrieve the error 

    ## tab1t table utile pour la realisation des figures  / table used for the figures
    tab1t <- NULL
    if (length(pval) > 0)
    {
        tab1t <- data.frame(Est=trend,
                            #LL , UL,
                            pourcent=pourcentage,signif=pval<seuilSignif,pval)
    }
    ## Tableau 2

    if(sp == "global")
    {
        datatablecut <- datatable[grep("FALSE",is.na(datatable[,as.character(metric)])),]

    }else{

        datatablecut <- datatable[datatable[,as.character(cut)]==sp,]
        datatablecut <- datatablecut[grep("FALSE",is.na(datatablecut[,as.character(metric)])),]
    }

    switch(as.character(metric),
           "number" = {valplot <- lapply(sort(year), FUN=function(x){sum(na.omit(as.numeric(subset(datatablecut, year==x)[,as.character(metric)])))})},
           "pres.abs" = {nb_loc <- lapply(sort(year), FUN=function(x){length(unique(subset(datatablecut,year == x)[,"location"]))}) ## nb_loc nombre de loc suivie par year / number of plots per year
                         nb_loc_presence <- lapply(sort(year), FUN=function(x){length(unique(subset(datatablecut[datatablecut[,metric] > 0,], year == x)[,"location"]))}) ## nb_loc_presence nombre de location de presence par year / number the plots where the species were observed
                         valplot <- (na.omit(as.numeric(nb_loc_presence)) / na.omit(as.numeric(nb_loc)))*100},
           {valplot <- lapply(sort(year), FUN=function(x){mean(na.omit(as.numeric(subset(datatablecut, year==x)[,as.character(metric)])))})}
          )

    tab2 <- data.frame(year,val = round(as.numeric(valplot),2),
                       LL=NA,UL=NA,catPoint=NA,pval=NA,
                       courbe=vpan[2],
                       panel=vpan[2])

    ## les figures     

    dgg <- tab1
  
    figname<- paste(sp,".png",sep = "")

  ## coordonnee des ligne horizontal de seuil pour les abondances et les occurences
    hline.data1 <- data.frame(z = c(1), panel = c(vpan[1]),couleur = "variation abondance",type="variation abondance")
    hline.data3 <- data.frame(z = 0, panel = vpan[2] ,couleur = "seuil",type="seuil")  
    hline.data <- rbind(hline.data1,hline.data3)
    titre <- paste(sp)#,"\n",min(year)," - ",max(year),sep="")

  ## texte de la tendance / text for the population evolution trend

    pasdetemps <- max(dgg$year) - min(dgg$year) + 1
    if (! is.null(tab1t))
    {
        if(assessIC){
            txtPente1 <- paste("Global trend : ", tab1t$Est,
                               ifelse(tab1t$signif," *",""),
                               ifelse(tab1t$signif,paste("\n",ifelse(tab1t$pourcent>0,"+ ","- "), 
                                                         abs(tab1t$pourcent)," % in ",pasdetemps," years",sep=""),""),sep="")
        }else{  
            txtPente1 <- ifelse(tab1t$signif,paste("\n",ifelse(tab1t$pourcent>0,"+ ","- "),
                                                   abs(tab1t$pourcent)," % in ",pasdetemps," years",sep=""),"")
        }
    }else{
        tendanceSurFigure <- FALSE
    }

  ## table du texte de la tendance / table of the text for the population evolution trend
    tabTextPent <- data.frame(y=c(max(c(dgg$val,dgg$UL),na.rm=TRUE)*.9),
                              x=median(dgg$year),
                              txt=ifelse(tendanceSurFigure,c(txtPente1),""),
                              courbe=c(vpan[1]),panel=c(vpan[1]))
    dgg <- rbind(tab1,tab2)
  ## les couleurs / the colors
    vecColPoint <- c("#ffffff","#eeb40f","#ee0f59")
    names(vecColPoint) <- c("significatif","infSeuil","0")
    vecColCourbe <- c("#3c47e0","#5b754d","#55bb1d","#973ce0")
    names(vecColCourbe) <- c(vpan[1],"loc","presence",vpan[2])
    vecColHline <- c("#ffffff","#e76060")
    names(vecColHline) <- c("variation abondance","seuil")
  
    col <- c(vecColPoint,vecColCourbe,vecColHline)
    names(col) <- c(names(vecColPoint),names(vecColCourbe),names(vecColHline))
  
  ## si description graphique en 3 panels

    if(description) {
        p <- ggplot(data = dgg, mapping = aes(x = year, y = val))
    ## Titre, axes ...
        p <- p + facet_grid(panel ~ ., scale = "free") +
        theme(legend.position="none",
              panel.grid.minor=element_blank(),
              panel.grid.major.y=element_blank())  +
              ylab("") + xlab("year")+ ggtitle(titre) +
              scale_colour_manual(values=col, name = "" ,
                                  breaks = names(col))+
              scale_x_continuous(breaks=min(dgg$year):max(dgg$year))
        p <- p + geom_hline(data =hline.data,mapping = aes(yintercept=z, colour = couleur,linetype=type ),
                        alpha=1,size=1.2)
    if(assessIC){ ############# ONLY FOR THE CONFIDENCE INTERVAL
    p <- p + geom_ribbon(mapping=aes(ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2) 
    p <- p + geom_pointrange(mapping= aes(y=val,ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2)
	}
    p <- p + geom_line(mapping=aes(colour=courbe),size = 1.5)
    p <- p + geom_point(mapping=aes(colour=courbe),size = 3)
    p <- p + geom_point(mapping=aes(colour=catPoint,alpha=ifelse(!is.na(catPoint),1,0)),size = 2)
    p <-  p + geom_text(data=tabTextPent, mapping=aes(x,y,label=txt) ,parse=FALSE,color=col[vpan[1]],fontface=2, size=4)
    ggsave(figname, p,width=16,height=15, units="cm")
	#print (figname)  ##### CAN BE REMOVED IF YOU DO NOT WANT THE GRAPH TO BE PLOTTED

  } else {

    p <- ggplot(data = subset(dgg,panel=="Variation abondance"), mapping = aes(x = year, y = val))
    ## Titre, axes ...

    p <- p + facet_grid(panel ~ ., scale = "free") +
      theme(legend.position="none",
            panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank())  +
      ylab("") + xlab("year")+ ggtitle(titre) +
      scale_colour_manual(values=col, name = "" ,
                          breaks = names(col))+
      scale_x_continuous(breaks=min(dgg$year):max(dgg$year))
    p <- p + geom_hline(data =subset(hline.data,panel==vpan[1]),mapping = aes(yintercept=z, colour = couleur,linetype=type ),
                        alpha=1,size=1.2)
    
   if(assessIC){ ############# ONLY FOR THE CONFIDENCE INTERVAL
    p <- p + geom_ribbon(mapping=aes(ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2) 
    p <- p + geom_pointrange(mapping= aes(y=val,ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2)
	}
    p <- p + geom_line(mapping=aes(colour=courbe),size = 1.5)
    p <- p + geom_point(mapping=aes(colour=courbe),size = 3)
    p <- p + geom_point(mapping=aes(colour=catPoint,alpha=ifelse(!is.na(catPoint),1,0)),size = 2)
    p <-  p + geom_text(data=tabTextPent, mapping=aes(x,y,label=txt),parse=FALSE,color=col[vpan[1]],fontface=2, size=4)
    ggsave(figname, p,width=15,height=9,units="cm")
  #print (figname) ##### CAN BE REMOVED IF YOU DO NOT WANT THE GRAPH TO BE PLOTTED
    
   #return(p)

  }
}
############################################################################################################ fin fonction graphique / end of function for graphical output

################# Analysis
#plots <- list()

for (sp in glmtable[,1]) 
{

    if (!all(is.na(glmtable[glmtable[,1]==sp,4:(length(glmtable)-1)]))) ##ignore lines with only NA
    { 
        #p <- 
        ggplot.glm(glmtable = glmtable, datatable=datatable,unitobs=unitobs, serie=NULL, sp=sp, description=TRUE, tendanceSurFigure=TRUE, seuilOccu=14, assessIC=TRUE)
        #plots <- list(plots,p)
    }
}

#plot <- rbind_ggplot_timeseries(ggplot_list=plots, limits = c(dmy("2008","2017"))) ##Bind several timeseries plots : to investigate
#ggsave("tot.png", plot ,width=16,height=15, units="cm")

