#_______________________________________________________________________________
# TREATMENT OF BIOLOGIC DATA FROM CSLN BEFORE SDM TREATMENT
# Amelie LEHUEN Janvier 2022
#_______________________________________________________________________________
# library(diffr)
# setwd("C:/Users/lehuen201/Nextcloud/Melting Pot/BDD/A_SDM_NEO/Scripts")
# diffr("1.2 Community Analysis GIPSA Zone.R","1.2 Community Analysis GIPSA Station.R")

rm(list=ls())
#________________________________________________________________
# PACKAGES USED ----
# install.packages("devtools"); library(devtools);Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(readxl) ; library(openxlsx) # Edition d'un fichier Excel
library(tidyverse) #The toolbox indispensable
library(vegan) # traitement des donnees envir (diversity, adonis2,metaMDS,vegdist)
library(labdsv) # Analyse des communautes (invdal)
library(pairwiseAdonis) #(pairwise.adonis)
library(pastecs) ; library(clustsig) ; library(dendextend) # # Analyse TPA (abund); Analyse SIMPROF (simprof); library(ggdendro) ggdendrogramm
# Graphics packages
library(RColorBrewer) ; library(wesanderson) # Palettes de couleurs
library(ggpubr);library(gridExtra) ; library(grid) # Mozaic of graphs tools
library(treemap) # pour les graphiques treemap
# GIS Packages
library(sf); library(sfheaders); # st_as_sf ; sf_to_df
library(tmap) # tmap_mode; for static and interactive maps
library(htmlwidgets) # library(leaflet) # saveWidget ; for interactive maps
library(rnaturalearth) # ne_states install_github("ropensci/rnaturalearthhires")
# library(spData) ; library(spDataLarge) #remotes::install_github("Nowosad/spDataLarge")
# library(raster) 
options(scipen=999) # Empeche affichage scientifique des nombres

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "C:/Users/lehuen201/Nextcloud/" # "E:/" # 
tsk <- "A_SDM_NEO/"
wdtask <- paste(pc,"Melting Pot/BDD/",tsk,sep="")
wdsource <- paste(wdtask,"Sources/Faune/CSLN/",sep="")
wdwork <- paste(wdtask,"Matrices/",sep="")
wdgraph <- paste(wdtask,"Graphiques/",sep="")
wdres <- paste(wdtask,"Resultats/",sep="")
wdGIS <- paste(pc,"Melting Pot/SIG/",sep="");
wdlogos <- paste(pc,"Melting Pot/Rapport - Presentations/Images/Logos/",sep="") # To add logos on map
wdscript <- (paste(pc,"Melting Pot/BDD/Scripts/",sep=""))
# wdmsr <- (paste(wdscript,"MSR/MSR.R",sep=""))
# setwd(paste(wdtask,"Scripts/",sep=""))
setwd(wdtask)
load(file = paste(wdwork,"CSLN_Mars_BDD",".RData", sep=""))
# If exists
# load(file = paste(wdwork,"CSLN_Mars_GIPSAZone_BDD",".RData", sep=""))

#________________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colDarj <- function(x) {wes_palette("Darjeeling2",x, type = "continuous")}
colZiss <- function(x) {wes_palette("Zissou1",x, type = "continuous")}
colSpec <- colorRampPalette(brewer.pal(11, "Spectral")); 
colDark <- colorRampPalette(brewer.pal(8, "Dark2"));
Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
Scale_fill <- function(x) {scale_fill_manual(values=colDarj(x))}

#________________________________________________________________
# SELECTION OF DATA ONLY MUDFLAT FOR GIPSA ----
CSLN_Mars<-CSLN_Mars %>% filter(grepl("Mudflat",Zone))
CSLN<-CSLN_Mars %>%
    select(idStationUnique,Station_originelle,Zone,Period,Tidal_level,Annee,SP,SPCourt,Taxon_SNa,Biomass_gAFDWm2,Density_indm2,MSRtot)

#________________________________________________________________
# COMPARATIVE ANALYSIS ON ZONE AND PERIODS ----
# FAUNA SUMMARY BY ZONE : SPECIES CHOSEN
CSLN_sp <- CSLN %>% select(Zone,Period,Annee,SP,SPCourt,Density_indm2,Biomass_gAFDWm2,MSRtot) %>% # ,Tidal_level
                    filter(SP=="SpCh") %>% select(-SP) %>%
                    complete(nesting(Zone,Period,Annee),nesting(SPCourt), # ,Tidal_level
                               fill=list(Biomass_gAFDWm2=0,Density_indm2=0))
CSLN_spm <- CSLN_sp %>% 
                      group_by(Zone,Period,Annee,SPCourt) %>% # ,Tidal_level
                      summarise(Biomass_m=median(Biomass_gAFDWm2,na.rm=TRUE),Biomass_sd=sd(Biomass_gAFDWm2,na.rm=TRUE),
                                Density_m=median(Density_indm2,na.rm=TRUE),Density_sd=sd(Density_indm2,na.rm=TRUE),
                                MSRtot_m=median(MSRtot,na.rm=TRUE),MSRtot_sd=sd(MSRtot,na.rm=TRUE))
# DIVERSITY TABLE BY ZONE
CSLN_sm <- CSLN %>%
                      select(Zone,Period,Annee,SPCourt,Density_indm2,Biomass_gAFDWm2,MSRtot) %>% # ,Tidal_level
                      group_by(Zone,Period,Annee) %>% # ,Tidal_level
                      summarise(SR=n_distinct(SPCourt), n_records=n(),
                                Biomass_t=sum(Biomass_gAFDWm2,na.rm=TRUE),
                                MSRtot_t=sum(MSRtot,na.rm=TRUE)) %>%
                      unite("code",Zone,Period,Annee,sep = "_",remove = FALSE) # ,Tidal_level
# CONTINGENCY TABLE
CSLN_cont <- CSLN %>%
                      select(Zone,Period,Annee,SPCourt,Density_indm2) %>% # ,Tidal_level
                      complete(nesting(Zone,Period,Annee),nesting(SPCourt), # ,Tidal_level
                               fill=list(Density_indm2=0)) %>%
                      group_by(Zone,Period,Annee,SPCourt) %>% # ,Tidal_level
                      summarise(Density_m_=median(Density_indm2,na.rm=TRUE)) %>%
                      pivot_wider(names_from = SPCourt, values_from = Density_m_, values_fill = 0) %>%
                      unite("code",Zone,Period,Annee,sep = "_",remove = TRUE) # ,Tidal_level
CSLN_cont_mat<-as.matrix(CSLN_cont[,-1])
row.names(CSLN_cont_mat)<-CSLN_cont$code

# SHANNON & PIELOU INDEXES ----
# Pielou's index of equitability (J ): normalization of the Shannon-Wiener index (H'), 
# a value of taxonomic diversity as a function of the number of taxa per area and the 
# abundance of individuals within each taxon; 0 means that one taxon dominates the others,
# 1 means that there is an equitable distribution of individuals between taxa
CSLN_sm$Shannon<-diversity(CSLN_cont_mat,index="shannon")
CSLN_sm$Pielou<-CSLN_sm$Shannon/log2(CSLN_sm$SR)

CSLN_cont_name<-data.frame(code=row.names(CSLN_cont_mat))  %>% 
  separate(code, c("Zone", "Period", "Annee"), sep = "_",remove = TRUE) %>% # ,"Tidal_level"
  unite("ZP",Zone,Period,sep = "_",remove = FALSE) %>%  # ,Tidal_level
  unite("ZA",Zone,Annee,sep = "_",remove = FALSE) %>%  # ,Tidal_level
  unite("ZPA",Zone,Period,Annee,sep = "_",remove = FALSE)  # ,Tidal_level

#________________________________________________________________
# RARE SPECIES - TPA ----
# David, V., 2017. Traitement de donn?es en sciences environnementales, ISTE Editions. ed, Ecologie. Page 44
# The Abundance Sorting method  (TPA :Tri Par Abondance) was adapted from 
# Ibanez et al (1993) in the pastecs package (abund). It takes into account both the 
# number of 0 present in the database for a given species, but also the rare 
# species (absent from many stations) but abundant for some stations.
# A coefficient f, between 0 and 1, allows to adjust the weight given to the 
# frequency of null values for a species (f=0 if only this criterion is taken .
# into account) and to the abundance of species expressed in log

CSLN_abd<-abund(CSLN_cont_mat,f=0.2) #f=0,2 recommande
# Find where our species are in the TPA order to keep them
# intersect(colnames(extract(CSLN_abd,length(CSLN_abd$vr))),species[,1])
# Find the last rank of our species
minSpecies<-max(which(colnames(extract(CSLN_abd,length(CSLN_abd$vr))) %in% species[,1]))
# dev.print(device = png, file = "CSLN_TPA.png", width = 1000, height=600) # If save of graph
plot(CSLN_abd,xlab="Taxon",ylab="Abond rel", dpos=c(40,100),
     main=paste("Methode TPA f=",CSLN_abd$f,sep=""))
# dev.off()
CSLN_abd$n<-max(72,minSpecies) # identify(CSLN_abd) # Fin du plateau
CSLN_tpa_SPCourt<-colnames(extract(CSLN_abd,CSLN_abd$n))
CSLN_cont_tpa<-as.matrix(extract(CSLN_abd,CSLN_abd$n))
CSLN_cont_tpa<-subset(CSLN_cont_tpa,rowSums(CSLN_cont_mat)!=0)

# DELETION OF RARE SPECIES IN THE COMPLETE DATABASE -----
# CSLN <- CSLN[CSLN$SPCourt %in% CSLN_tpa_SPCourt,]

#________________________________________________________________
# PERMANOVA ----
  # Permanova, a permutation multivariate anova, is used to determine whether
  # there is a difference between the station-year pairs. This statistical test 
  # compares more than two variables in the same way as an Anova, but using a 
  # permutation of the data. The post-hoc pairwise test is used to find out which 
  # stations and/or years have a significant difference
  # The average contingency table per area reduced by TPA allows the establishment 
  # of a dissimilarity matrix by the Bray-Curtis method, which measures the differences
  # between factor1-factor2 pairs according to the taxa found there
# If Pr (p-value) is : Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# the factor variations explain R^2 % of the dissimilarities
#________________________________________________________________
# # ADONIS and PAIRWISE ADONIS - ZONE ----
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # # WARNING FOLLOWING OPERATION CAN LAST VERY VERY LONG
# ad<-adonis2(CSLN_cont_mat~CSLN_cont_name$Zone*CSLN_cont_name$Period,method="bray") #*CSLN_cont_name$Annee
# CSLN_adonis<-as.data.frame(cbind(test=row.names(ad),ad))
# padZ<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$Zone)
# padP<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$Period)
# # padA<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$Annee)
# padZP<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$ZP)
# # padZA<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$ZA)
# # padZPA<-pairwise.adonis(CSLN_cont_mat,CSLN_cont_name$ZPA)
# CSLN_pairado<-as.data.frame(rbind(padZ,padP,padZP)) #,padZA,padZPA
# CSLN_pairado<-CSLN_pairado %>% filter(sig != "") # Keep only significant pairs
# # # SAVE BECAUSE SO LONG TO CALCULATE
# save(ad,CSLN_adonis,CSLN_pairado,file = paste(wdwork,"CSLN_GIPSAZone_Ado.RData", sep=""))
# # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load(paste(wdwork,"CSLN_GIPSAZone_Ado.RData", sep=""))

#________________________________________________________________
# ANALYSE DES COMMUNAUTES ----
#________________________________________________________________
# Hierarchical Ascending Classification (HAC or dendrogram) ----
  # constructed with the UPGMA agglomeration method to define 'clustchoix' groups 
  # at an equal cut-off level

CSLN_cont_ZP <- CSLN %>%
  select(Zone,Period,SPCourt,Density_indm2) %>% # ,Tidal_level
  complete(nesting(Zone,Period),nesting(SPCourt), # ,Tidal_level
           fill=list(Density_indm2=0)) %>%
  group_by(Zone,Period,SPCourt) %>% # ,Tidal_level
  summarise(Density_m_=median(Density_indm2,na.rm=TRUE)) %>%
  pivot_wider(names_from = SPCourt, values_from = Density_m_, values_fill = 0) %>%
  unite("code",Zone,Period,sep = "_",remove = TRUE) # ,Tidal_level 
CSLN_cont_ZP_mat<-as.matrix(CSLN_cont_ZP[,-1])
row.names(CSLN_cont_ZP_mat)<-CSLN_cont_ZP$code
CSLN_cont_ZP_mat<-subset(CSLN_cont_ZP_mat,rowSums(CSLN_cont_ZP_mat)!=0)
CSLN_cont_ZP_name<-data.frame(code=row.names(CSLN_cont_ZP_mat))  %>% 
  separate(code, c("Zone", "Period"), sep = "_",remove = TRUE) %>% # ,"Tidal_level"
  unite("ZP",Zone,Period,sep = "_",remove = FALSE) # ,Tidal_level

CSLN_cont_t <- CSLN_cont_ZP_mat # CSLN_cont_tpa # CSLN_cont_mat # 
# any(is.nan(CSLN_cont_t))
# SIMPROF Determination of significantly different groups
CSLN_simprof<-simprof(CSLN_cont_t, num.expected=1000, num.simulated=999,
                      method.cluster="average",method.distance="braycurtis") # summary(CSLN_simprof)
# groupe.clust<-cutree(CSLN_simprof$hclust, CSLN_simprof$numgroups)
# groupesdend<-as.data.frame(cbind(code=row.names(as.data.frame(groupe.clust)),Groupe_simprof=groupe.clust))
groupesdend <- CSLN_simprof$significantcluster %>% 
  set_names(seq_along(.)) %>% enframe %>% unnest(cols=c(value)) %>%
  rename(Groupe_Simprof=name, code=value) %>% 
  mutate(Groupe_Simprof = parse_integer(Groupe_Simprof))
groupe.clust <- groupesdend$Groupe_Simprof; names(groupe.clust) <- groupesdend$code
simprof.plot(CSLN_simprof)
CSLN_sm <-CSLN_sm %>% unite(ZP,Zone,Period,sep="_") %>%
  left_join(groupesdend, by=c("ZP"="code")) %>% relocate(Groupe_simprof, .before="SR")
CSLN_sm$Groupe_simprof<-as.factor(CSLN_sm$Groupe_simprof)

col_Period<-colZiss(CSLN_unique$Period)[factor(CSLN_cont_ZP_name$Period)]
col_Zone<-colZiss(length(unique(CSLN_cont_ZP_name$Zone)))[factor(CSLN_cont_ZP_name$Zone)]
ggd1 <- CSLN_simprof$hclust %>% as.dendrogram %>% ladderize
ggd1 <- ggd1 %>% set("branches_k_color", k=CSLN_simprof$numgroups,value=colDarj(CSLN_simprof$numgroups)) %>% #
          set("branches_lwd", .8) %>%
          set("labels_col", col_Zone[order.dendrogram(ggd1)]) %>%
          set("labels_cex", .6) %>%
          # set("hang_leaves", 1) %>%
          set("leaves_pch", 19) %>%
          set("leaves_col", col_Period[order.dendrogram(ggd1)])
ggd1<- as.ggdend(ggd1)
ggd1 <- ggplot(ggd1, horiz = TRUE, labels = TRUE)
ggd1 <- ggd1 + ggtitle("Simprof on mudflats") ; ggd1
ggsave(paste(wdgraph,"CSLN_Simprof_Zone",".png",sep=""), plot = ggd1, width = 12, height = 8)

# Significance of the groups
ad_CAH<-adonis2(CSLN_cont_t~groupe.clust,method="bray"); ad_CAH
pad_CAH<-pairwise.adonis(CSLN_cont_t,groupe.clust); pad_CAH

#________________________________________________________________
# INDVAL ----
  # For the characteristic taxa of each of these groups, the Indval has been 
  # calculated (Species Indicator Values). This index combines the relative 
  # abundance of each species and its occurrence in the samples to define an index 
  # between 0 and 1. The higher the IndVal, the more characteristic the species 
  # is of the community structure of the group.
#________________________________________________________________
CSLN_cont_indval_mat<- CSLN_cont_ZP_mat # CSLN_cont_mat # CSLN_cont_tpa # 
CSLN_cont_indval_mat<-subset(CSLN_cont_indval_mat,select=colSums(CSLN_cont_indval_mat)!=0)
indval_mat <- indval(CSLN_cont_indval_mat,groupe.clust)
indvalresults <- data.frame(cbind(group=indval_mat$maxcls,
                                  indval=round(indval_mat$indcls,2),
                                  pvalue=round(indval_mat$pval,3)))
indvalresults <- cbind(sp=row.names(indvalresults),indvalresults)[order(indvalresults$group, indvalresults$pvalue, -indvalresults$indval),]
CSLN_indval_10<-indvalresults %>% group_by(group) %>% arrange(group) %>% slice(1:10) %>%
  pivot_wider(names_from = group, values_from = indval)

# IndVal list reduced to significant species
alpha=0.05 #0.05
gr <- indval_mat$maxcls[indval_mat$pval<=alpha]
iv <- indval_mat$indcls[indval_mat$pval<=alpha]
pv <- indval_mat$pval[indval_mat$pval<=alpha]
# fr <- apply(CSLN_cont_CAH[,-1]>0, 2, sum)[indval_mat$pval<=alpha]
indvalsummary <- data.frame(group=gr, indval=round(iv,2), pvalue=round(pv,3)) #, freq=fr
indvalsummary <- cbind(sp=row.names(indvalsummary),indvalsummary)[order(indvalsummary$group, -indvalsummary$indval),]
indres<-as.data.frame(indval_mat$indcls)
indres<-cbind(sp=row.names(indres),indres)

# save.image(file = paste(wdwork,"CSLN_BDD_GIPSAZone",".RData", sep=""))

#________________________________________________________________
# nMDS-----
  # An non-metric MultiDimensional Scaling represents the distances between
  # elements in a two-dimensional space. It is an iterative process of positioning
  # the data on a plane to best approximate the dissimilarity matrix. The calculated
  # stress is the indicator of the difference between the representation on the 
  # plane and the matrix. It should be at most 0.2, ideally less than 0.1, 
  # in order to reflect ecological and/or hydroclimatic factors.
#________________________________________________________________
#  NMDS by ZONE
CSLN_cont_nMDS<-CSLN_cont
CSLN_cont_nMDS_mat<-CSLN_cont_mat

CSLN_mds <- metaMDS(CSLN_cont_nMDS_mat,distance = "bray", k = 2,try = 1000, autotransform =FALSE)

CSLN_mds_dat <- data.frame(scores(CSLN_mds)) %>% # Using the scores function from vegan to extract the code scores
                    mutate(code=row.names(scores(CSLN_mds))) %>% relocate(code) %>%
                    separate(code, c("Zone", "Period","Annee"), sep = "_",remove = FALSE) %>%
                    unite(ZP,Zone,Period, sep = "_",remove = FALSE) %>%
                    group_by(Zone,Period) %>% #,Tidal_level
                    mutate(NMDS1.m = mean(NMDS1, na.rm = TRUE),
                           NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_mds_dat <- CSLN_mds_dat %>% left_join(groupesdend, by=c("ZP"="code"))

CSLN_mds_dat$Groupe_simprof<-as.factor(CSLN_mds_dat$Groupe_simprof)
rownames(CSLN_mds_dat) <- c()
CSLN_mds_mean <- CSLN_mds_dat %>% 
  unite("ZP",Zone,Period, sep = "_",remove = FALSE) %>%
  group_by(ZP) %>% 
  summarise(NMDS1.m = mean(NMDS1, na.rm = TRUE),
            NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_mds_species<-as.data.frame(CSLN_mds$species)
CSLN_mds_species<-cbind(sp=rownames(CSLN_mds_species),CSLN_mds_species)  # create a column of species, from the rownames of CSLN_mds_species

nbA = length(unique(CSLN_mds_dat$Zone)) # dim(CSLN_mds_mean)[1]
nbB = length(unique(CSLN_mds_dat$ZP))
palette=colZiss # colDarj colZiss colSpec
bp<-ggplot(data=CSLN_mds_dat) +
# Add species name
    # geom_point(data=CSLN_mds_species,aes(x=MDS1,y=MDS2),shape=8,size=1,colour="gray50") + # add the point markers
    geom_text(data=CSLN_mds_species,aes(x=MDS1,y=MDS2,label=sp),
              size=2.5,colour="gray50",check_overlap = TRUE, fontface = "italic") +
# Base scatter plot 
    geom_point(aes(x=NMDS1,y=NMDS2,shape=Period,colour=ZP),size=3) +
    geom_text(aes(x=NMDS1,y=NMDS2,colour=ZP,label=Annee),size=2.5,hjust=-0.3, fontface = "plain") +  # add the site labels
    scale_colour_manual(values=palette(nbB)) +
    stat_ellipse(aes(x=NMDS1,y=NMDS2,group=Groupe_simprof), #,linetype = groupe.clust,color=groupe.clust
                 type = "norm",level = 0.9,linetype = 3,show.legend = FALSE) + #
# Lines to connect the same Zone and label it
    geom_segment(aes(x=NMDS1.m, y=NMDS2.m, xend=NMDS1, yend=NMDS2, linetype = Period),
                 size=.1,color="grey") +
    annotate(geom="label",x=CSLN_mds_mean$NMDS1.m,y=CSLN_mds_mean$NMDS2.m,label=CSLN_mds_mean$ZP,
             fill=palette(nbB),alpha=0.4) +
# Graphic visuals
    annotate("text",x=min(CSLN_mds_dat$NMDS1),y=min(CSLN_mds_dat$NMDS2),hjust=-0.2,
             label=paste("Stress =",round(CSLN_mds$stress,4)), size = 4)+
    ggtitle(label="nMDS CSLN data by Station_originelle")+
    theme(axis.text.x = element_blank(),  # remove x-axis text
          axis.text.y = element_blank(), # remove y-axis text
          axis.ticks = element_blank());bp #+theme_bw()
ggsave(paste(wdgraph,"nMDS Macrofaune by Zone",".png",sep=""), plot = bp, width = 12, height = 8)

# #________________________________________________________________
# # ANALYSES STATISTIQUES ----
# CSLN_stats <- data.frame(matrix(NA,ncol = 10, nrow = length(anneesData))) #NULL
# colnames(CSLN_stats) <- c("Categorie","Variable","Annee","Mois","Shapiro","A","B","C","D","Kruskal") #NULL
# k=1
# for (j in 1:length(anneesBiv)){
  # df<-with(CSLN_summary,CSLN_summary[Annee==anneesBiv[j],])
  # shap<-by(df$Densite.m2.mea,df$Zone,shapiro.test)#Si p>alpha normalite
#   krus<-kruskal.test(Densite.m2.mea~as.factor(Zone),data=df)#Si p<alpha diff entre groupes
#   CSLN.zone.stat.s[k,]<-c("CSLN","Densite",as.character(anneesData[j]),"Mois","Shapiro",
#                           round(shap$A$p.value,4),round(shap$B$p.value,4),
#                           round(shap$C$p.value,4),round(shap$D$p.value,4),round(krus$p.value,4))
#   k<-k+1
# }
# #glob<-pairwise.wilcox.test(df$Densite.m2,df$Zone,p.adj="bonf")

# _______________________________________________________________
# OUTPUT SAVE ----
  wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep=""))
  writeData(wb, sheet = "Ado_Zone", x = CSLN_adonis, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "PairAdo_Zone", x = CSLN_pairado, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "indval_Zone", x = CSLN_indval_10, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Smry_Zone", x = CSLN_sm, startCol = 1, startRow = 1,withFilter = FALSE)
  saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)
  
  save.image(file = paste(wdwork,"CSLN_Mars_GIPSAZone_BDD",".RData", sep=""))
  
# GIS SAVE
  # write.csv(CSLN,file=paste(wdres,"CSLN_Mars", ".csv",sep=""), na = "",row.names = FALSE)
  # write.csv(CSLN_pur,file=paste(wdres,"CSLN_data", ".csv",sep=""), na = "",row.names = FALSE)
  write.csv(CSLN_Stations,file=paste(wdres,"CSLN_Stations_Zone", ".csv",sep=""), na = "",row.names = FALSE)
  
