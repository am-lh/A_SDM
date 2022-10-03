# ______________________________________________________________________________
# CREATION DE GRAPHIQUES DESCRIPTIFS DES DONNEES PROGRAMME REPONSES PREDICTEURS
# ______________________________________________________________________________
# NB animation graphique en time series https://www.datanovia.com/en/blog/
# gganimate-how-to-create-plots-with-beautiful-animation-in-r/

rm(list=ls())
#________________________________________________________________
# PACKAGES USED ----
library(tidyverse); 
library(ggExtra); library(gridExtra) # ;ggMarginal; grid.arrange
library(cowplot); library(ggpubr); library(grid) ; # plot_grid; ggscatterhist
library(GGally); library(RColorBrewer)

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "E:/" #"C:/Users/lehuen201/Nextcloud/" "E:/" 
tsk <- "A_SDM_NEO/"
wdpath <- paste(pc,"Melting Pot/BDD/",tsk,sep=""); 
wdwork <- paste(wdpath,"Matrices/",sep="")
wdgraph <- paste(wdpath,"Graphiques/",sep="")
wdres <- paste(wdpath,"Resultats/",sep="")
setwd(wdpath)

#________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
prgm <- "CSLN" # CSLN Mabes Geco Beaug
etude <- paste(prgm,"_Mars",sep="")
load(paste(wdwork,etude,"_BDD.RData",sep=""))
# load(paste(wdwork,"CSLN_Mars_GIPSAZone_BDD",".RData", sep=""))

# Choix espece etudiee, reponse...
spe <- 1#:nrow(species) # 1:CERED 2:CORVO 3:HEDDI 4:LIMBA 5:PERUL 6:SCRPL
reponse<-reponse[1:3,] # 1:Biomass_gAFDWm2 2:Density_indm2 3:MSRtot 4:Itot  
answ <- 1:nrow(reponse)
sai <- 1#:nrow(saison) # 1:Year 2:Winter 3:Summer
pred<- 1:nrow(predict)

# __________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colDarj <- function(x) {wes_palette("Darjeeling2",x, type = "continuous")}
colZiss <- function(x) {wes_palette("Zissou1",x, type = "continuous")}
colSpec <- colorRampPalette(brewer.pal(8, "Spectral")); 
colDark <- colorRampPalette(brewer.pal(8, "Dark2"));
Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
Scale_fill <- function(x) {scale_fill_manual(values=colDarj(x))}

#________________________________________________________________
###### FONCTION GGSCATTERHIST MAISON POUR NE PAS AVOIR PB AVEC GGSAVE ET GGARRANGE ######
# Scatter plot colored by groups ("Species")
ggscatterhist<-function (df, x = x, y = y,color = z,
                         title=title, ylab = ylab, xlab = xlab,
                         palette = palette,legend = legend){
    pmain <- ggplot(df, aes_string(x = x, y = y, color = z))+
                    geom_point()+
                    labs(title=title,x=xlab, y=ylab) +
                    scale_fill_manual(values=palette)+
                    theme(legend.position=legend,legend.title = element_blank())+
                    theme_bw()
    pmain<-ggMarginal(pmain, groupColour = TRUE, groupFill = TRUE)
      # xbox <- axis_canvas(pmain, axis = "x", coord_flip = TRUE) +
      #   geom_density(df, aes_string(x = x, color = z)) 
      # ybox <- axis_canvas(pmain, axis = "y") +
      #   geom_density(df, aes_string(x =  y, color = z)) +
      #   coord_flip()
      # p1 <- insert_xaxis_grob(pmain, xbox, grid::unit(2, "in"), position = "top")
      # p2 <- insert_yaxis_grob(p1, ybox, grid::unit(2, "in"), position = "right")
      # ggdraw(p2)
}

#________________________________________________________________
# GLOBAL VISION OF DATA ----
# Visualisation d'ensemble matriciel ----
for (sp in spe){
  for (sa in sai){
    df <- CSLN_Mars %>% filter(SPCourt == species[sp,1]) %>%
      select(paste(predict[,1],saison[sa,1],sep=""),reponse[,1])
    colnames(df)<-c(predict[,2],reponse[,2])
    titreG <- paste("Mars3D Predictors for ",species[sp,2]," in ",saison[sa,2]," - ", prgm, sep="")
    tp <- ggpairs(df, title=titreG) +
      theme(plot.title = element_text(size=18,face="bold"),
            strip.text.x = element_text(size=12,face="bold"),
            strip.text.y = element_text(size=10,face="bold"))
    ggsave(paste(wdgraph,species[sp,1],"/",etude,"_MatRP_",species[sp,1],saison[sa,1],".png",sep=""), plot = tp, width = 16, height = 9)
    # Matrice des correlations simple ----
    cp <- ggcorr(df, method = c("pairwise", "pearson"),
                 low="#FC4E07", mid="white", high="#00AFBB", nbreaks = 5,# palette=RdBu,
                 geom="tile", hjust = 0.85, angle = 0, layout.exp=3, size=4,
                 label=TRUE,label_alpha = TRUE,label_round=2,label_size=4)+
      labs(title=paste("Mars3D Predictors for ",species[sp,2]," in ",saison[sa,2]," - ", prgm, sep=""))+
      theme(plot.title = element_text(size=18,face="bold"))+
      theme_bw()
    ggsave(paste(wdgraph,species[sp,1],"/",etude,"_MatRP_corr_",species[sp,1],saison[sa,1],".png",sep=""), plot = cp, width = 9, height = 9)
}}
# ----

# Diagramme de densite des predicteurs Mars3D ----
df<-CSLN_Mars
mlist<-vector(mode = "list", length = length(pred))
for (pre in pred){ # pre=1
  x <- predict[pre,1]
  z <- "Zone"
  dp <- ggplot(df)+
    # geom_density(aes(flow_m,fill=Zone), alpha=0.5) +
    geom_density(aes_string(x,fill=z,colour=z), alpha=0.5) + #, na.rm = TRUE
    xlab(paste(predict[pre,2]," (",predict[pre,3],")",sep="")) +
    theme(legend.position="bottom",legend.title = element_blank()) +
    Scale_col(CSLN_unique$Zone) + Scale_fill(nbcol)
  assign(paste('dp', pre, sep=''), dp + theme(legend.position="none"))
  dp <- dp + ggtitle(paste("Mars3D Data Density", " - ", prgm, sep=""))
  mlist[[pre]] <- eval(parse(text = paste('dp', pre, sep='')))
  ggsave(paste(wdgraph,"Predictors/",etude,"_GphR_","Mars3D","_Density","_",predict[pre,1],".png",sep=""), plot = dp, width = 12, height = 8)
}
if (length(pred)>1 ){
    titreG <- paste("Mars3D Data Density", " - ", prgm, sep="")
    nncol = 3; nnrow=ceiling(nrow(predict)/nncol)
    dp <- ggarrange(plotlist=mlist, ncol=nncol, nrow=nnrow, labels="AUTO",legend="bottom",common.legend = TRUE)#
    dp <- annotate_figure(dp, top = text_grob(titreG, face = "bold", size = 14))+bgcolor("white"); dp
    ggsave(paste(wdgraph,etude,"_GphR_Mars3D","_Density",".png",sep=""), plot = dp, width = 12, height = 12)
  # dev.off()
}
# ----

#________________________________________________________________
# Reponse vs Predicteurs avec densite en marge ----
  for (sp in spe){ # sp=1
    df <- CSLN_Mars[which(CSLN_Mars$SPCourt == species[sp,1]),]
      for (rep in answ){  # rep=3
        mlist<-vector(mode = "list", length = length(pred))
        for (sa in sai){ # sa=1
          for (pre in pred){  # pre=1
            x <- paste(predict[pre,1],saison[sa,1],sep="")
            y <- reponse[rep,1]
            z <- "Zone" ;
            titreG <- paste(reponse[rep,2]," for " ,species[sp,2]," vs ",predict[pre,2],
                         " (",predict[pre,3],")"," in ",saison[sa,2], " - ", prgm, sep="")
            dp <- ggscatterhist(df, x = x, y = y,color = z, 
                                title=titreG, 
                                ylab = paste(reponse[rep,2], " (",reponse[rep,3],")",sep=""), 
                                xlab = paste(predict[pre,2], " (",predict[pre,3],")"," in ",saison[sa,2],sep=""),
                                palette = colSpec(CSLN_unique$Zone),legend = "bottom")
            ggsave(paste(wdgraph,species[sp,1],"/Detail/",etude,"_GphRP_",species[sp,1],"_",predict$Couche[pre],"_",y,"_",x,".png",sep=""),plot = dp, width = 8, height = 6) #plot = dp,
            dp <- ggscatterhist(df, x = x, y = y,color = z, 
                                title="", ylab = "", 
                                xlab = paste(predict[pre,2], " (",predict[pre,3],")"," in ",saison[sa,2],sep=""),
                                palette = colSpec(CSLN_unique$Zone),legend = "none")
            assign(paste('dp', pre, sep=''), dp)
            mlist[[pre]] <- eval(parse(text = paste('dp', pre, sep='')))
        }
        if (length(pred)>1 ){
            titreG <- paste(predict$Couche[pre]," Data for ",species[sp,2]," in ",saison[sa,2], " - ", prgm, sep="")
            nncol = 3; nnrow=ceiling(nrow(predict)/CSLN_unique$Zone)
            dp <- grid.arrange(grobs=mlist, ncol=CSLN_unique$Zone)
            # dp <- ggarrange(plotlist=mlist, ncol=nncol, nrow=nnrow, labels="AUTO",legend="bottom",common.legend = TRUE) #
            # METTRE LA GESTION DE LA LEGENDE DES COULEURS A TROUVER
            dp <- annotate_figure(dp,top = text_grob(titreG, face = "bold", size = 14), left= paste(reponse[rep,2], " (",reponse[rep,3],")",sep=""))+bgcolor("white")
            # ggsave(paste(wdgraph,species[sp,1],"/",etude,"_GphRP_",species[sp,1],"_",predict$Couche[pre],"_",y,saison[sa,1],".png",sep=""), plot = dp, width = 12, height = 8)
            # save_plot(paste(wdgraph,species[sp,1],"/",etude,"_GphRP_",species[sp,1],"_",predict$Couche[pre],"_",y,saison[sa,1],".png",sep=""),dp)
            png(file=paste(wdgraph,species[sp,1],"/",etude,"_GphRP_",species[sp,1],"_",y,saison[sa,1],".png",sep=""),width=1000, height=1200,bg = "white")
            print(dp)
            # dev.off()
          }
}}}
# ----


#________________________________________________________________
# GRAPHES MNT ----
#________________________________________________________________
# Diagramme de densite des predicteurs ----
df<-CSLN_Mars
mlist<-vector(mode = "list", length = nrow(predictMNT))
for (pre in 1:nrow(predictMNT)){ #
  x <- predictMNT[pre,1]
  y <- "Annee"
  z <- "Zone"
  dp <- ggplot(df)+
    geom_density(aes_string(x,fill=z,colour=z), alpha=0.5, na.rm = TRUE) +
    xlab(paste(predictMNT[pre,2]," (",predictMNT[pre,3],")",sep="")) +
    theme(legend.position="bottom",legend.title = element_blank()) +
    Scale_col(CSLN_unique$Zone) + Scale_fill(nbcol)
  assign(paste('dp', pre, sep=''), dp + theme(legend.position="none"))
  dp <- dp + ggtitle(paste("MNT Data Density", " - ", prgm, sep=""))
  mlist[[pre]] <- eval(parse(text = paste('dp', pre, sep='')))
  ggsave(paste(wdgraph,"Predictors/",etude,"_GphR_","MNT","_Density","_",predictMNT[pre,1],".png",sep=""), plot = dp, width = 12, height = 8)
}
if (length(pred)>1 ){
  titreG <- paste("MNT Data Density", " - ", prgm, sep="")
  nncol = 3; nnrow=ceiling(nrow(predictMNT)/nncol)
  dp <- ggarrange(plotlist=mlist, ncol=nncol, nrow=nnrow, labels="AUTO",legend="bottom",common.legend = TRUE)#
  dp <- annotate_figure(dp, top = text_grob(titreG, face = "bold", size = 14))+bgcolor("white") ; dp
  ggsave(paste(wdgraph,etude,"_GphR_","MNT","_Density",".png",sep=""), plot = dp, width = 12, height = 6)
  # dev.off()
}
# ----

# Visualisation d'ensemble matriciel ----
for (sp in spe){
  df <- CSLN_Mars %>% filter(SPCourt == species[sp,1]) %>%
    select(predictMNT[,1],reponse[,1])
  colnames(df)<-c(predictMNT[,2],reponse[,2])
  titreG <- paste("MNT Predictors for ",species[sp,2]," - ", prgm, sep="")
  # tp <- ggpairs(df, title=titreG) + 
  #   theme(plot.title = element_text(size=18,face="bold"),
  #         strip.text.x = element_text(size=12,face="bold"),
  #         strip.text.y = element_text(size=10,face="bold"))
  # ggsave(paste(wdgraph,species[sp,1],"/",etude,"_MatRP_",species[sp,1],saison[sa,1],".png",sep=""), plot = tp, width = 16, height = 9)
  # Matrice des correlations simple ----
  cp <- ggcorr(df, method = c("pairwise", "pearson"),
               low="#FC4E07", mid="white", high="#00AFBB", nbreaks = 5,# palette=RdBu,
               geom="tile", hjust = 0.85, angle = 0, layout.exp=3, size=4,
               label=TRUE,label_alpha = TRUE,label_round=2,label_size=4)+
    labs(title=paste("HMS Predictors for ",species[sp,2]," in ",saison[sa,2]," - ", prgm, sep=""))+
    theme(plot.title = element_text(size=18,face="bold"))
  ggsave(paste(wdgraph,species[sp,1],"/",etude,"_Mat_",species[sp,1],"_MNT RepPredict corr",".png",sep=""), plot = cp, width = 9, height = 9)
}
# ----

#________________________________________________________________
# Reponse vs Predicteurs avec densite en marge ----
for (sp in spe){
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == species[sp,1]) & !is.na(CSLN_Mars$Zone),]
  for (rep in answ){
    mlist<-vector(mode = "list", length = nrow(predictMNT))
    for (pre in 1:nrow(predictMNT)){
      x <- predictMNT[pre,1]
      y <- reponse[rep,1]
      z <- "Zone"
      titreG <- paste(reponse[rep,2]," for " ,species[sp,2]," vs ",predictMNT[pre,2],
                      " (",predictMNT[pre,3],")", " - ", prgm, sep="")
      dp <- ggscatterhist(df, x = x, y = y,color = z, 
                          title=titreG, 
                          ylab = paste(reponse[rep,2], " (",reponse[rep,3],")",sep=""), 
                          xlab = paste(predictMNT[pre,2], " (",predictMNT[pre,3],")"," in ",saison[sa,2],sep=""),
                          palette = colSpec(CSLN_unique$Zone),legend = "bottom")
      ggsave(paste(wdgraph,species[sp,1],"/Detail/",etude,"_GphRP_",species[sp,1],"_",y,"_",x,".png",sep=""),plot = dp, width = 8, height = 6) #plot = dp,
      dp <- ggscatterhist(df, x = x, y = y,color = z, 
                          title="", ylab = "", 
                          xlab = paste(predictMNT[pre,2], " (",predictMNT[pre,3],")"," in ",saison[sa,2],sep=""),
                          palette = colSpec(CSLN_unique$Zone),legend = "none")
      assign(paste('dp', pre, sep=''), dp)
      mlist[[pre]] <- eval(parse(text = paste('dp', pre, sep='')))
    }
    titreG <- paste("MNT Data for ",species[sp,2], " - ", prgm, sep="")
    nncol = 3; nnrow=ceiling(nrow(predictMNT)/CSLN_unique$Zone)
    dp <- grid.arrange(grobs=mlist, ncol=CSLN_unique$Zone)
    dp <- annotate_figure(dp,top = text_grob(titreG, face = "bold", size = 14), left= paste(reponse[rep,2], " (",reponse[rep,3],")",sep=""))+bgcolor("white")
    png(file=paste(wdgraph,species[sp,1],"/",etude,"_GphRP_",species[sp,1],"_MNT_",y,saison[sa,1],".png",sep=""),width=1000, height=600,bg = "white")
    print(dp)
    dev.off()
  }}
# ----

# # ________________________________________________________________
# # COMPARAISON MUD CONTENT DE CSLN_Mars ET SILT ET ARGILES DE MNT ----
# # ________________________________________________________________
# # mud <-CSLN_Mars[,c("idStationUnique","siltsArgiles","mud_content_med","Zone")]
# mud <-na.omit(mud)
# mud <-unique(mud)
# nbcol <- length(levels(factor(mud$Zone)))
# mud<-ggplot(mud, aes_string("siltsArgiles","mud_content_med",color="Zone"))+
#   geom_point() +  Scale_col(length(nbcol)) +
#   labs(title="MNT mud and silt content versus MARS3D mud content calculated via Csed",
#        caption="n=937",
#        y="Mud content from Mars3D", x="Mud content from MNT (%)")
# ggsave(paste(wdgraph,etude,"_GphR_","_Mud_MNT_MARS",".png",sep=""), plot = mud, width = 8, height = 6)
# # ----

