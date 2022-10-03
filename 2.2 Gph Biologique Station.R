# ______________________________________________________________________________
# CREATION DE GRAPHIQUES DESCRIPTIFS DES DONNEES PROGRAMME REPONSES PREDICTEURS
# ______________________________________________________________________________
# NB animation graphique en time series https://www.datanovia.com/en/blog/
# gganimate-how-to-create-plots-with-beautiful-animation-in-r/

rm(list=ls())
#________________________________________________________________
# PACKAGES USED ----
library(tidyverse)
library(ggpubr) ; library(grid) ; library(gridExtra)
library(RColorBrewer) ; library(wesanderson) 
library(treemap)

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "E:/" #"C:/Users/lehuen201/Nextcloud/" "E:/" 
tsk <- "A_SDM_NEO/"
wdpath <- paste(pc,"Melting Pot/BDD/",tsk,sep=""); 
wdwork <- paste(wdpath,"Matrices/",sep="")
wdgraph <- paste(wdpath,"Graphiques/",sep="")
wdres <- paste(wdpath,"Resultats/",sep="")
setwd(paste(wdpath,"Scripts/",sep=""))

#________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
prgm <- "CSLN" # CSLN Mabes Geco Beaug
etude <- paste(prgm,"_Stat",sep="")
load(paste(wdwork,"CSLN_Mars_GIPSAStation_BDD",".RData", sep=""))


spe <- 1:6#nrow(species) # 1:CERED 2:CORVO 3:HEDDI 4:LIMBA 5:PERUL 6:SCRPL
reponse<-reponse[1:3,] # 1:Biomass_gAFDWm2 2:Density_indm2 3:MSRtot 4:Itot  
answ <- 1:nrow(reponse)

#________________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colDarj <- function(x) {wes_palette("Darjeeling2",x, type = "continuous")}
colZiss <- function(x) {wes_palette("Zissou1",x, type = "continuous")}
colSpec <- colorRampPalette(brewer.pal(8, "Spectral")); 
colDark <- colorRampPalette(brewer.pal(8, "Dark2"));
Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
Scale_fill <- function(x) {scale_fill_manual(values=colDarj(x))}

#________________________________________________________________
# GLOBAL VISION OF DATA ----
# Nb record per species
dp<-CSLN %>% count(SP,SPCourt,Taxon_SNa)%>%group_by(SP)%>%
              arrange(desc(n),.by_group = TRUE)%>%slice(1:14) %>%
              ggplot(aes(x=reorder(Taxon_SNa,n),y=n,fill = SP))+geom_col() +
              geom_text(aes(y=60, label = n)) + #
              scale_fill_manual(values=colDarj(CSLN_unique$Annee))+
              labs(title="List of top 20 occurences species, including model species", x="Species",y="Occurence")+
              theme(axis.text = element_text(size=15,face="bold"))+
              coord_flip()+theme_bw();dp
ggsave(sprintf("%s%s_Species_Occr_List20.png",wdgraph,etude), plot = dp, width = 10, height = 8)


# Treemap 2 levels : Tidal levels distribution ----
titre = paste("Total Density (ind.m-2) of Species"," - ", etude,sep="")
Sumry<-CSLN %>% group_by(Tidal_level,SPCourt) %>%
  summarise(Smry=sum(Density_indm2,na.rm=TRUE),Nb=n()) #n_distinct(idStationUnique)
Sumry$label <- paste(Sumry$SPCourt, Sumry$Nb, sep = "\n")
png(filename=paste(wdgraph,etude,"_TreeMap2 Tidal Level Density",".png",sep=""),width=1000, height=600)
treemap(Sumry,index=c("Tidal_level","label"),vSize="Smry",vColor="Tidal_level",type="index",
        title=titre, palette=colDarj(CSLN_unique$Annee), fontface.labels=c(2,2),
        border.col=c("black","white"), border.lwds=c(3,1),bg.labels=0,
        fontsize.labels=c(12,9), fontcolor.labels=c("black","white"),
        align.labels=list( c("left","bottom"), c("right","center")))
dev.off()

titre = paste("Total Biomass (gAFDW.m-2) of Species"," - ", etude,sep="")
Sumry<-CSLN %>% group_by(Tidal_level,SPCourt) %>% 
  summarise(Smry=sum(Biomass_gAFDWm2,na.rm=TRUE),Nb=n()) #n_distinct(idStationUnique)
Sumry$label <- paste(Sumry$SPCourt, Sumry$Nb, sep = "\n")
png(filename=paste(wdgraph,etude,"_TreeMap2 Tidal Level Biomass",".png",sep=""),width=1000, height=600)
treemap(Sumry,index=c("Tidal_level","label"),vSize="Smry",vColor="Tidal_level",type="index",
        title=titre, palette=colDarj(CSLN_unique$Annee), fontface.labels=c(2,2),
        border.col=c("black","white"), border.lwds=c(3,1),bg.labels=0,
        fontsize.labels=c(12,9), fontcolor.labels=c("black","white"),
        align.labels=list( c("left","bottom"), c("right","center")))
dev.off()

# Treemap 3 levels: Distribution spatial ----
titre = paste("Total Density (ind.m-2) of Species"," - ", etude,sep="")
Sumry<-CSLN %>% group_by(Zone,Tidal_level,SPCourt) %>%
  summarise(Smry=sum(Density_indm2,na.rm=TRUE),Nb=n())
Sumry$label <- paste(Sumry$SPCourt, Sumry$Nb, sep = "\n")
png(filename=paste(wdgraph,etude,"_TreeMap3 Zone Density",".png",sep=""),width=1000, height=600)
treemap(Sumry,index=c("Zone","Tidal_level","label"),vSize="Smry",vColor="Zone",type="index",
        title=titre, palette=colDarj(CSLN_unique$Annee), fontface.labels=c(2,2,1),
        border.col=c("black","white","gray81"), border.lwds=c(4,3,1),bg.labels=0,
        fontsize.labels=c(14,14,10), fontcolor.labels=c("black","white","gray81"),
        align.labels=list( c("left","bottom"), c("left","top"), c("right","center")))
dev.off()
# ----

#________________________________________________________________
# FONCTION boxplot with facets
boxplot_amlh <- function(df,x,y,fil,face,xlab,ylab,title,palette){
  bp<-ggplot(df) + geom_boxplot(aes(x=x,y=y,fill=fil)) + 
    facet_wrap(face) + # scales="free_y", margins=TRUE 
    scale_fill_manual(values=palette) +
    labs(title=title, x=xlab, y=ylab, fill="") +
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))
}
#FONCTION Bargraph with error bars and facets
bargraph_sd <- function(df,x,y,sd,fil,face,xlab,ylab,title,palette){
  empilement = .9 # .9 pour des groupes par zone cote cote, .6 pour empilement
  bp <-ggplot(df, aes(x=x, y=y, fill=fil))+
    geom_bar(stat="identity", position=position_dodge(width=empilement))+
    geom_errorbar(aes(ymin=y-sd, ymax=y+sd),
                  width=.2,position=position_dodge(empilement))+
    labs(title=title,x=xlab, y=ylab, fill="")+
    facet_wrap(face)+
    theme(legend.position="bottom",legend.title = element_blank())+
    scale_fill_manual(values=palette)
}

#________________________________________________________________
# BOXPLOT BY PERIOD AND FACET ZONE
titre <- paste("Periodic Evolution of Specific Richness",sep="")
bp<-boxplot_amlh(CSLN_sm,CSLN_sm$Period,CSLN_sm$SR,CSLN_sm$Period,CSLN_sm$Groupe_simprof,
                 "Years","Specific Richness",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) ;bp 
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp1<-bp+ labs(title="Specific Richness",x="",y="") + theme(legend.position='none')

titre <- paste("Periodic Evolution of Pielou",sep="")
bp<-boxplot_amlh(CSLN_sm,CSLN_sm$Period,CSLN_sm$Pielou,CSLN_sm$Period,CSLN_sm$Groupe_simprof,
                 "Years","Index of Equitability",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) ;bp 
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp2<-bp+ labs(title="Index of Equitability",x="",y="") + theme(legend.position='none')

titre <- paste("Periodic Evolution of Total Biomass",sep="")
bp<-boxplot_amlh(CSLN_sm,CSLN_sm$Period,CSLN_sm$Biomass_t,CSLN_sm$Period,CSLN_sm$Groupe_simprof,
                 "Years","Total Biomass (gAFDW/m²)",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) +
  ylim(0,750);bp 
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp3<-bp+ labs(title="Total Biomass (gAFDW/m²)",x="",y="") + theme(legend.position='none')

titre <- paste("Periodic Evolution of Metabolic Rate",sep="")
bp<-boxplot_amlh(CSLN_sm,CSLN_sm$Period,CSLN_sm$MSRtot_t,CSLN_sm$Period,CSLN_sm$Groupe_simprof,
                 "Years","Total MSRtot (mW/m²)",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) ;bp 
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp4<-bp+ labs(title="Total MSRtot (mW/m²)",x="",y="") + theme(legend.position='none')

titre <- paste("Periodic Evolution of Biodiversity",sep="")
bp <- ggarrange(bp1,bp2,bp3,bp4, ncol=2, nrow=2, labels="AUTO", legend="bottom", common.legend = TRUE)
bp <- annotate_figure(bp,top = textGrob(titre, gp=gpar(fontsize=15,font=1)),
                      bottom=xlab)+bgcolor("white"); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 15, height = 15)

#________________________________________________________________
## BARPLOT BY PERIOD AND FACET ZONE
# df  <-  CSLN_sm %>% unite(AS,Annee,Season,sep = "_", remove = FALSE)
df <- CSLN_sm %>% group_by(Groupe_simprof,Period) %>% #
  summarise(SR_m=median(SR,na.rm=TRUE),SR_sd=sd(SR,na.rm=TRUE),
            Biomass_t_m=median(Biomass_t,na.rm=TRUE),Biomass_t_sd=sd(Biomass_t,na.rm=TRUE),
            Pielou_m=median(Pielou,na.rm=TRUE),Pielou_sd=sd(Pielou,na.rm=TRUE),
            MSRtot_t_m=median(MSRtot_t,na.rm=TRUE),MSRtot_t_sd=sd(MSRtot_t,na.rm=TRUE))
df$Groupe_simprof<-as.factor(df$Groupe_simprof)

xlab <- ""
titre <- paste("Temporal Evolution of Specific Richness",sep="")
bp <-bargraph_sd(df,df$Period,df$SR_m,df$SR_sd,df$Period,df$Groupe_simprof,
                 xlab,"Specific Richness",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp1<-bp+ labs(title="Specific Richness",y="") + theme(legend.position='none')

titre <- paste("Temporal Evolution of Density",sep="")
bp <-bargraph_sd(df,df$Period,df$Biomass_t_m,df$Biomass_t_sd,df$Period,df$Groupe_simprof,
                 xlab,"Total Biomass (gAFDW/m²)",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp2<-bp+ labs(title="Total Biomass (gAFDW/m²)",y="") + theme(legend.position='none')

titre <- paste("Temporal Evolution of Pielou index",sep="")
bp <-bargraph_sd(df,df$Period,df$Pielou_m,df$Pielou_sd,df$Period,df$Groupe_simprof,
                 xlab,"Index of Equitability",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp3<-bp+ labs(title="Index of Equitability",y="") + theme(legend.position='none')

titre <- paste("Temporal Evolution of metabolic rate",sep="")
bp <-bargraph_sd(df,df$Period,df$MSRtot_t_m,df$MSRtot_t_sd,df$Period,df$Groupe_simprof,
                 xlab,"Total MSRtot (mW/m²)",titre,colDarj(length(unique(CSLN_mds_dat$Groupe_simprof)))) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 12, height = 6)
bp4<-bp+ labs(title="Total MSRtot (mW/m²)",y="") + theme(legend.position='none')

# Un plot pour les 3 differents graphiques
titre=paste("Temporal Evolution of Benthic fauna",sep="")
bp <- ggarrange(bp1,bp2,bp3,bp4, ncol=2, nrow=2, labels="AUTO", legend="bottom", common.legend = TRUE)
bp <- annotate_figure(bp,top = textGrob(titre, gp=gpar(fontsize=15,font=1)),
                      bottom=xlab)+bgcolor("white"); bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 15, height = 15)

#________________________________________________________________
# CHOSEN SPECIES FOCUS ----
CSLN_sp <- CSLN_sp #%>% filter(SPCourt %in% species[1:6,1])
titre <- paste("Periodic Evolution of Biomass for chosen species",sep="")
bp<-ggplot(CSLN_sp)+geom_boxplot(aes(x=Period,y=Biomass_gAFDWm2,fill=Period))+ 
  facet_grid(SPCourt~Groupe_simprof)+ # Period, scales="free_y", margins=TRUE
  scale_y_continuous(trans = "log10") + 
  theme(legend.position="bottom",
        axis.text.x=element_blank())+ #element_text(angle=-90, vjust=0.4,hjust=1)) +
  labs(title=titre,y="Biomass (gAFDW/m²)") +
  Scale_fill(CSLN_unique$Annee) ;bp
ggsave(paste(wdgraph,etude," ",titre,".png",sep=""), plot = bp, width = 9, height = 9)

# Boxplot des reponses par annee pour 1 esp avec facet ----
for (sp in spe){ #sp=1
  df <- CSLN[which(CSLN$SPCourt == species[sp,1]),]
  mlist<-vector(mode = "list", length = length(answ))
  for (rep in answ){ #rep=1
    x <- "Annee"; y <- reponse[rep,1]
    z <- "Period"
    titreG <- paste(y," time evolution for ",species[sp,2], " - ", etude,sep="")
    bp <- ggplot(df, aes_string(x,y,fill=z)) + geom_boxplot() + 
      facet_wrap(~Zone) + #,scales="free_y"
      theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
      labs(title=titreG, fill="")+
      theme(legend.position="bottom",legend.title = element_blank())+
      scale_y_continuous(trans = "log10") + annotation_logticks(sides="lr") +
      Scale_fill(CSLN_unique$Annee) #bp
    ggsave(paste(wdgraph,species[sp,1],"/",etude,"_Bxp_",species[sp,1],"_",y,".png",sep=""), plot = bp, width = 12, height = 6)
    bp <- bp + labs(title=sprintf("%s (%s)",reponse[rep,2],reponse[rep,3]) ,y="") + theme(legend.position='none')
    assign(paste('bp', rep, sep=''), bp + theme(legend.position="none"))
    mlist[[rep]] <- eval(parse(text = paste('bp', rep, sep='')));
  }
  if (length(answ)>1 ){
    titreG <- paste("Responses time evolution for ",species[sp,2], " - ", etude,sep="")
    nncol = length(answ); nnrow=ceiling(length(answ)/nncol)
    bp <- ggarrange(plotlist=mlist,ncol=nncol, nrow=nnrow, labels="AUTO",legend="bottom",common.legend = TRUE)
    bp <- annotate_figure(bp, top = text_grob(titreG, face = "bold", size = 14))+bgcolor("white"); print(bp)
    ggsave(paste(wdgraph,species[sp,1],"/",etude,"_Bxp_",species[sp,1],"_Answ",".png",sep=""), plot = bp, width = 16, height = 6)
  }}

# # Boxplot des Densites moy par annee pour 1 esp avec facet ----
# for (sp in spe){ 
#   df <- CSLN_sp[which(CSLN_sp$SPCourt == species[sp,1]),]
#   titreG <- paste("Time evolution for ",species[sp,2], " - ", etude,sep="")
#   bp <- ggplot(df, aes(Annee,Density_indm2,fill=Period)) + geom_boxplot() + 
#     facet_wrap(~Zone) +
#     labs(title=titreG) + 
#     scale_y_continuous(trans = "log10") + annotation_logticks(sides="lr") + 
#     theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
#     Scale_fill(CSLN_unique$Annee); bp
#   ggsave(paste(wdgraph,species[sp,1],"/",etude,"_Bxp1Zone_",species[sp,1],"_","Density_indm2",".png",sep=""), plot = bp, width = 12, height = 6)
# }