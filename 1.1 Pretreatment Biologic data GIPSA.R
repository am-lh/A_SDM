#_______________________________________________________________________________
# TREATMENT OF BIOLOGIC DATA FROM CSLN BEFORE SDM TREATMENT
# Amelie LEHUEN Janvier 2022
#_______________________________________________________________________________
diffr("1.0 Pretreatment Biologic data.R","1.1 Pretreatment Biologic data GIPSA.R")

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
library(pastecs) ; library(clustsig)  # Analyse TPA (abund); Analyse SIMPROF (simprof)
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
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colDarj <- function(x) {wes_palette("Darjeeling2",x, type = "continuous")}
colZiss <- function(x) {wes_palette("Zissou1",x, type = "continuous")}
colSpec <- colorRampPalette(brewer.pal(8, "Spectral")); 
colDark <- colorRampPalette(brewer.pal(8, "Dark2"));
Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
Scale_fill <- function(x) {scale_fill_manual(values=colDarj(x))}

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "E:/" #"C:/Users/lehuen201/Nextcloud/" # 
tsk <- "A_SDM_NEO/"
wdtask <- paste(pc,"Melting Pot/BDD/",tsk,sep="")
wdsource <- paste(wdtask,"Sources/Faune/CSLN/",sep="")
wdwork <- paste(wdtask,"Matrices/",sep="")
wdgraph <- paste(wdtask,"Graphiques/",sep="")
wdres <- paste(wdtask,"Resultats/",sep="")
wdGIS <- paste(pc,"Melting Pot/SIG/",sep="");
wdlogos <- paste(pc,"Melting Pot/Rapport - Presentations/Images/Logos/",sep="") # To add logos on map
wdscript <- (paste(pc,"Melting Pot/BDD/Scripts/",sep=""))
wdmsr <- (paste(wdscript,"MSR/MSR.R",sep=""))
# setwd(paste(wdtask,"Scripts/",sep=""))
setwd(wdtask)
# If exists
# load(file = paste(wdwork,"CSLN_Mars_GIPSA_BDD",".RData", sep=""))

# ________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
# to be upgraded
species <-data.frame(SPCourt=c("CERED","CORVO","HEDDI","MACBA","PERUL","SCRPL"),
                         Taxon_SNa=c("Cerastoderma edule","Corophium volutator","Hediste diversicolor","Macoma balthica","Peringia ulvae","Scrobicularia plana"))
reponse <-data.frame(rvar=c("Biomass_gAFDWm2","Density_indm2","MSRtot","Itot"),
                          rdescr=c("Biomass","Density","Specific Respiration Rate","Metabolic Rate"),
                          runit=c("gAFDW/m2","ind/m2","mW/m2","mW/m2"))
predictMNT <- data.frame(pvar=c("moyenneMu","modeMu","medianeMu","siltsArgiles","sablesFins","sablesMoyens","sablesGrossiers","graviers"),
                          pdescr=c("Granulometric Mean","Granulometric Mode","Granulometric Median","Mud and Silts","Light Sands","Medium Sands","Coarse Sands","Gravels"),
                          punit=c("µm","µm","µm","%","%","%","%","%"))

#________________________________________________________________
# LOADING DATA ----
# FAUNA
fauna_file <- paste(wdsource,"CSLN_Biology_source.xlsx",sep="")
  CSLN <- as.data.frame(read_excel(fauna_file,sheet = "Biology_station", na = ""))
  granulo <- as.data.frame(read_excel(fauna_file,sheet = "Granulo", na = ""))
  # summary(CSLN)
# MARS
mars_file <- paste(wdres,"ES_Ncf_BDD.xlsx",sep="")
  Varnames <- as.data.frame(read_excel(mars_file,sheet = "Varnames", na = "",col_names = c("Var", "Desc", "Unit", "NAval", "Dim", "Varid")))
  predict <- as.data.frame(read_excel(mars_file,sheet = "Predicteurs", na = "",col_names = c("Var", "Desc", "Unit","Couche")))
  # MEVE <- as.data.frame(read_excel(mars_file,sheet = "MEVE", na = "",col_names = c("Var", "Desc", "Unit")))
  saison <- as.data.frame(read_excel(mars_file,sheet = "Saison",col_names = c("Suff", "M_Def", "M_Per")))
saison[1,1] <- "" # Remplacement du NaN pour l'annee
# GIS
ES_Areas<-paste(wdGIS,"Layers made/ES_Areas_WGS.shp",sep="")
  ES_Areas<-st_read(ES_Areas,quiet=TRUE,crs=4326) %>% st_transform(2154)
  ES_Areas<-ES_Areas %>% filter(Zone!="NA") # | Zone !="Bay"
Mars_csv<-paste(wdres,"ES_Mars_Maps_sh.csv",sep="")
  Mars_csv<-read.csv(Mars_csv); Mars_csv[Mars_csv == "NaN"] <- NA
  Mars_csv<-Mars_csv %>% filter(!is.na(Lat) & !is.na(Lon) & !is.na(flow_m))# %>% filter(Lat!="NaN" | Lon!="NaN")
  # TIDAL_LEVEL CALCULATION : should be in mars matlab script ? ----
  Mars_csv$Tidal_level<-NA
  Mars_csv$Tidal_level[Mars_csv$inunt>=0 & Mars_csv$inunt<0.25]<-"Supratidal"
  Mars_csv$Tidal_level[Mars_csv$inunt>=0.25 & Mars_csv$inunt<0.75]<-"Intertidal"
  Mars_csv$Tidal_level[Mars_csv$inunt>=0.75 & Mars_csv$inunt<=1]<-"Infratidal"
  Mars_csv$Tidal_level<-as.factor(Mars_csv$Tidal_level)
  # Definition of temporal periods
  Mars_csv$Period<-"2011-2018"
  Mars_csv$Period[Mars_csv$Annee %in% 1990:1999]<-"1990-1999"
  Mars_csv$Period[Mars_csv$Annee %in% 2000:2010]<-"2000-2010"
  Mars_csv$Period<-as.factor(Mars_csv$Period)
  Mars_csv_sf<-st_as_sf(Mars_csv, coords=c("Lon","Lat"),crs=4326,na.fail=FALSE,remove = FALSE) %>% st_transform(2154)
Mars_dat_sf<-paste(wdGIS,"Layers made/ES_Maille_nc.shp",sep="")
  Mars_dat_sf<-st_read(Mars_dat_sf,quiet=TRUE,crs=4326) %>% st_transform(2154) %>%
               rename(NINJ=NINJ_v, Lon=Lon_c, Lat=Lat_c) %>% select(-id_maille)
  # TIDAL_LEVEL CALCULATION : should be in mars matlab script ? ----
  Mars_dat_sf$Tidal_level<-NA
  Mars_dat_sf$Tidal_level[Mars_dat_sf$inunt>=0 & Mars_dat_sf$inunt<0.25]<-"Supratidal"
  Mars_dat_sf$Tidal_level[Mars_dat_sf$inunt>=0.25 & Mars_dat_sf$inunt<0.75]<-"Intertidal"
  Mars_dat_sf$Tidal_level[Mars_dat_sf$inunt>=0.75 & Mars_dat_sf$inunt<=1]<-"Infratidal"
  Mars_dat_sf$Tidal_level<-as.factor(Mars_dat_sf$Tidal_level)
  # Definition of temporal periods
  Mars_dat_sf$Period<-"2011-2018"
  Mars_dat_sf$Period[Mars_dat_sf$Annee %in% 1990:1999]<-"1990-1999"
  Mars_dat_sf$Period[Mars_dat_sf$Annee %in% 2000:2010]<-"2000-2010"
  Mars_dat_sf$Period<-as.factor(Mars_dat_sf$Period)
  anMars<-unique(Mars_csv$Annee)
# Mars_mesh<-paste(wdGIS,"Maillage Mars3D/Maillage_mars3D_WGS.shp",sep="")
  # Mars_mesh<-st_read(Mars_mesh,quiet=TRUE,crs=4326) %>% st_transform(2154)
  # cadreES_Mars<-Mars_mesh[Mars_mesh$NINJ=="Cadre_ES",] # Keep aside ES frame in shp
  # Mars_mesh<-Mars_mesh[Mars_mesh$NINJ!="Cadre_ES",]# Remove of ES frame in shp
  
#________________________________________________________________
# PRELIMINARY TREATMENT ----
colnames(CSLN)[which(colnames(CSLN) =="ScientificName_accepted")]<- "Taxon_SNa"
CSLN <-CSLN %>% left_join(granulo[,-c(21:24)]) %>% 
                arrange(idStationUnique, Taxon_SNa) %>% 
                mutate(IndBodySize_gAFDW=Biomass_gAFDWm2/Density_indm2)
CSLN$IndBodySize_gAFDW[CSLN$Density_indm2==0]<-0
CSLN$SP<-"OTHER"; CSLN$SP[CSLN$SPCourt %in% species$SPCourt]<-"SpCh";
# CSLN$SP<-"OTHER"; CSLN$SP[CSLN$SPCourt %in% species$SPCourt]<-"SpCh"; CSLN$SP[CSLN$SPCourt %in% speciesB$SPCourt]<-"SpMSR";
# Definition of seasons and temporal periods
CSLN$Season[CSLN$Mois %in% 1:3]   <- "Q1"
CSLN$Season[CSLN$Mois %in% 4:6]   <- "Q2"
CSLN$Season[CSLN$Mois %in% 7:10]  <- "Q3"
CSLN$Season[CSLN$Mois %in% 11:12] <- "Q4"
CSLN$Period<-"2011-2018"
CSLN$Period[CSLN$Annee %in% 1990:1999]<-"1990-1999"
CSLN$Period[CSLN$Annee %in% 2000:2010]<-"2000-2010"

# Deleting unused columns
suppr_col<-c("Source","Site","Campagne","Jour","id_station","Taxon", #"Station_originelle",
             "Engin","maille_tamis","Forme_maille","nb_replicat_protocole","nb_replicat",
             "surface_engin","surface_station","Abondance","PSLCUnitaire",
             "identite_tri","identite_det","Operateur","commentaire","Methode")
CSLN<-CSLN %>% select(-all_of(suppr_col))

# Patch pour sauver les INTERMUD (avant 2000) ----
tmp<-CSLN %>% filter(Filtre=="INTERMUD") %>% select(Taxon_SNa) %>% 
                              left_join(CSLN %>% group_by(Taxon_SNa) %>% 
                                       summarise(medInd=median(IndBodySize_gAFDW, na.rm = TRUE))) %>%
                              select(medInd)
CSLN$Biomass_gAFDWm2[CSLN$Filtre=="INTERMUD"]<-as.vector(tmp$medInd)

# Deleting unusable records
CSLN<-CSLN %>% filter(!is.na(Density_indm2), !is.na(Biomass_gAFDWm2), !is.na(AphiaID_accepted),
                      longitude!="NaN",latitude!="NaN",
                      ((Biomass_gAFDWm2!=0 & Density_indm2!=0) | (Biomass_gAFDWm2==0 & Density_indm2==0)))
# Removal of duplicates on the idStationUnique and SPCourt pairs
CSLN<-CSLN %>% distinct(idStationUnique,SPCourt, .keep_all = TRUE) %>% arrange(Annee)
# summary(CSLN)

#________________________________________________________________
# INTEGRATION OF MARS IN FAUNA DATA -----
  # # CSLN_sf <- CSLN_sf %>% st_join(Mars_csv_sf,join=st_nearest_feature,by=c("NINJ","Annee")) %>%
  # #                        setNames(gsub('\\.x$', '', names(.)))# %>% # Rename those with .x suffixes in original name
  # #                        # select(-c("Annee.y","NINJ.y")) # Remove extra fields
# CONVERSION OF FAUNA DATA INTO GEOGRAPHIC TABLE
CSLN_sf <- st_as_sf(CSLN, coords=c("longitude","latitude"),crs=4326,remove = FALSE) %>% st_transform(2154) # transform to planar as required by st_intersection()
# st_crs(CSLN_sf) # TO CHECK GEOG DATA

# ASSIGNMENT OF AREAS TO MARS3D
Mars_csv_sf <- st_intersection(Mars_csv_sf,ES_Areas)
Mars_csv_sf$Zone<-as.factor(Mars_csv_sf$Zone)
Mars_dat_sf <- st_intersection(Mars_dat_sf,ES_Areas)
Mars_dat_sf$Zone<-as.factor(Mars_dat_sf$Zone)

# ADD NINJ OF MARS_CSV FOR EACH SAMPLING POINT IN FAUNA
CSLN_sf$NINJ <- Mars_csv_sf$NINJ[st_nearest_feature(CSLN_sf,Mars_csv_sf)]

# ASSIGNMENT OF AREAS TO SAMPLING POINTS
CSLN_sf <- CSLN_sf %>%  mutate(intersection = as.integer(st_intersects(geometry, ES_Areas)),
                              Zone = if_else(is.na(intersection),"NA", ES_Areas$Zone[intersection]),
                              Type = if_else(is.na(intersection),"NA", ES_Areas$Type[intersection])) # Type : french areas name
CSLN_sf$Zone[CSLN_sf$Zone=="NA"]<-NA; CSLN_sf$Type[CSLN_sf$Type=="NA"]<-NA;
# # SELECTION OF DATA POINT ONLY IN ES_Areas
# CSLN_sf <- st_intersection(CSLN_sf, st_union(ES_Areas)) %>% select(-intersection)

# FAUNA TABLE WITH MARS DATA
CSLN_Mars<-sf_to_df(CSLN_sf,fill=TRUE)
CSLN_Mars<-CSLN_Mars %>% left_join(Mars_csv,by=c("NINJ"="NINJ","Annee"="Annee","Period"))
# CSLN_Mars<-CSLN_Mars %>% relocate(c("Lon","Lat","x","y"),.after=NINJ)

# # # TIDAL_LEVEL CALCULATION : should be in mars matlab script ? ----
# # Add missing tidal levels on other years : TO BE DELETED WHEN ALL RUN ON MARS ?
# CSLN_Mars[CSLN_Mars$NINJ %in% CSLN_Mars$NINJ[which(CSLN_Mars$Tidal_level=="Intertidal")],
#           "Tidal_level"]<-"Intertidal"
# CSLN_Mars[CSLN_Mars$NINJ %in% CSLN_Mars$NINJ[which(CSLN_Mars$Tidal_level=="Infratidal")],
#           "Tidal_level"]<-"Infratidal"
# CSLN_Mars[CSLN_Mars$NINJ %in% CSLN_Mars$NINJ[which(CSLN_Mars$Tidal_level=="Supratidal")],
#           "Tidal_level"]<-"Supratidal"

# Suppression of missing information for levels
CSLN_Mars<-CSLN_Mars %>% filter(!is.na(Zone)) #filter(!is.na(Tidal_level) & !is.na(Zone))

#________________________________________________________________
# METABOLIC RATE CALCULATION (with yearly med temp for the moment) ----
CSLN_Mars$MSR_mW<-NA
setwd(wdscript) # when the package is done, with repertory in it, no more needed
source(wdmsr)
for (i in which(!is.na(CSLN_Mars$temp_m) & CSLN_Mars$SP=="SpCh")){
# for (i in which(!is.na(CSLN_Mars$temp_m) & CSLN_Mars$SP!="OTHER")){
  CSLN_Mars$MSR_mW[i]<-msr(CSLN_Mars$Density_indm2[i],CSLN_Mars$Biomass_gAFDWm2[i],
                           CSLN_Mars$AphiaID_accepted[i],CSLN_Mars$temp_m[i],1,"w")*1000 #Mean MSR (mW/ind)
}
setwd(wdtask) # when the package is done, with repertory in it, no more needed
CSLN_Mars$MSRtot<-CSLN_Mars$MSR_mW*CSLN_Mars$Density_indm2
# CSLN_Mars$Itot<-CSLN_Mars$MSRtot

#________________________________________________________________
# FINAL SET FOR BASIS TABLE ----
CSLN_Mars<-CSLN_Mars %>% select(-all_of(c("AphiaID_accepted","sfg_id","point_id","Lon","Lat","x","y","intersection")))
facto_col <-c("Filtre","Zone","Type","Tidal_level",
              "Period","Season","Annee","Mois",
              "idStationUnique","Station_originelle","NINJ",
              "SP","Taxon_SNa","SPCourt") #
CSLN_Mars<-CSLN_Mars %>% group_by(list(facto_col)) %>% arrange(list(facto_col), .by_group = TRUE) %>% ungroup
CSLN_Mars<-CSLN_Mars %>% select(-`list(facto_col)`)
# MOVE DESCRIPTIVE FIELDS AT BEGINNING OF TABLE
CSLN_Mars<-CSLN_Mars %>% relocate(c(all_of(facto_col)))
CSLN_Mars<-CSLN_Mars %>% relocate(c("IndBodySize_gAFDW","MSR_mW","MSRtot"),.after=Biomass_gAFDWm2) #,"Itot"
# FACTORISATION FOR BASIS TABLE ----
CSLN_Mars[,facto_col] <- lapply(CSLN_Mars[,facto_col], as.factor)

# Creation of vector with number of unique values for all field in databasis
CSLN_unique<-CSLN_Mars %>% summarise(across(everything(),n_distinct))
annees <- CSLN_Mars %>% distinct(Annee) %>% arrange(Annee) %>% pull(Annee)
mois <- CSLN_Mars %>% distinct(Mois) %>% arrange(Mois) %>% pull(Mois)

CSLN_pur<-CSLN_Mars # CSLN_Mars<-CSLN_pur #summary(CSLN_pur) # Save before more modifications
# summary(CSLN_Mars)

#________________________________________________________________
# SELECTION OF DATA ONLY MUDFLAT FOR GIPSA ----
CSLN_Mars<-CSLN_Mars %>% filter(grepl("Mudflat",Zone)) %>% 
    select(idStationUnique,Station_originelle,Zone,Period,Tidal_level,Annee,SP,SPCourt,Taxon_SNa,Biomass_gAFDWm2,Density_indm2)

#________________________________________________________________
# COMPARATIVE ANALYSIS ON ZONE AND PERIODS ----
# FAUNA SUMMARY BY ZONE : SPECIES CHOSEN
CSLN_Mars_spZ <- CSLN_Mars %>% 
                      select(Zone,Period,Annee,SP,SPCourt,Density_indm2,Biomass_gAFDWm2) %>% # ,Tidal_level
                      filter(SP=="SpCh") %>% select(-SP) %>%
                      complete(nesting(Zone,Period,Annee),nesting(SPCourt), # ,Tidal_level
                               fill=list(Biomass_gAFDWm2=0,Density_indm2=0))
CSLN_Mars_spmZ <- CSLN_Mars_spZ %>% 
                      group_by(Zone,Period,Annee,SPCourt) %>% # ,Tidal_level
                      summarise(Biomass_m_Z=median(Biomass_gAFDWm2,na.rm=TRUE),Biomass_sd_Z=sd(Biomass_gAFDWm2,na.rm=TRUE),
                               Density_m_Z=median(Density_indm2,na.rm=TRUE),Density_sd_Z=sd(Density_indm2,na.rm=TRUE))
# DIVERSITY TABLE BY ZONE
CSLN_Mars_smZ <- CSLN_Mars %>%
                      select(Zone,Period,Annee,SPCourt,Density_indm2,Biomass_gAFDWm2) %>% # ,Tidal_level
                      group_by(Zone,Period,Annee) %>% # ,Tidal_level
                      summarise(SRZ=n_distinct(SPCourt), n_records=n(),
                                Biomass_t_Z=sum(Biomass_gAFDWm2,na.rm=TRUE))
# CONTINGENCY TABLE
CSLN_contZ <- CSLN_Mars %>%
                      select(Zone,Period,Annee,SPCourt,Density_indm2) %>% # ,Tidal_level
                      complete(nesting(Zone,Period,Annee),nesting(SPCourt), # ,Tidal_level
                               fill=list(Density_indm2=0)) %>%
                      group_by(Zone,Period,Annee,SPCourt) %>% # ,Tidal_level
                      summarise(Density_m_Z=median(Density_indm2,na.rm=TRUE)) %>%
                      pivot_wider(names_from = SPCourt, values_from = Density_m_Z, values_fill = 0) %>%
                      unite("code",Zone,Period,Annee,sep = "_",remove = TRUE) # ,Tidal_level
CSLN_contZ_mat<-as.matrix(CSLN_contZ[,-1])
row.names(CSLN_contZ_mat)<-CSLN_contZ$code

# SHANNON & PIELOU INDEXES ----
# Pielou's index of equitability (J ): normalization of the Shannon-Wiener index (H'), 
# a value of taxonomic diversity as a function of the number of taxa per area and the 
# abundance of individuals within each taxon; 0 means that one taxon dominates the others,
# 1 means that there is an equitable distribution of individuals between taxa
CSLN_Mars_smZ$Shannon<-diversity(CSLN_contZ_mat,index="shannon")
CSLN_Mars_smZ$Pielou<-CSLN_Mars_smZ$Shannon/log2(CSLN_Mars_smZ$SRZ)

CSLN_contZ_mat<-subset(CSLN_contZ_mat,rowSums(CSLN_contZ_mat)!=0)
CSLN_contZ_name<-data.frame(code=row.names(CSLN_contZ_mat))  %>% 
  separate(code, c("Zone", "Period", "Annee"), sep = "_",remove = TRUE) %>% # ,"Tidal_level"
  unite("ZP",Zone,Period,sep = "_",remove = FALSE) %>%  # ,Tidal_level
  unite("ZA",Zone,Annee,sep = "_",remove = FALSE) %>%  # ,Tidal_level
  unite("ZPA",Zone,Period,Annee,sep = "_",remove = FALSE)  # ,Tidal_level

#________________________________________________________________
# GRAPHICS
titre <- paste("Periodic Evolution of Specific Richness",sep="")
bp<-ggplot(CSLN_Mars_smZ)+geom_col(aes(x=Annee,y=SRZ,fill=Period))+ 
  facet_grid(.~Zone)+ # Period, scales="free_y", margins=TRUE 
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
  labs(title=titre,y="Specific Richness") +
  Scale_fill(length(annees)) ;bp 
# ggsave(paste(wdgraph,"Graph ",titre,".png",sep=""), plot = bp, width = 18, height = 6)
# bp1<-bp+ labs(title="",x="") + theme(legend.position='none')
titre <- paste("Periodic Evolution of Pielou",sep="")
bp<-ggplot(CSLN_Mars_smZ)+geom_col(aes(x=Annee,y=Pielou,fill=Period))+ 
  facet_grid(.~Zone)+ # Period, scales="free_y", margins=TRUE
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
  labs(title=titre,y="Index of Equitability") +
  Scale_fill(length(annees)) ;bp 
# ggsave(paste(wdgraph,"Graph ",titre,".png",sep=""), plot = bp, width = 18, height = 6)
# bp1<-bp+ labs(title="",x="") + theme(legend.position='none')
titre <- paste("Periodic Evolution of Total Biomass",sep="")
bp<-ggplot(CSLN_Mars_smZ)+geom_col(aes(x=Annee,y=Biomass_t_Z,fill=Period))+ 
  facet_grid(.~Zone)+ # Period, scales="free_y", margins=TRUE
  ylim(0,750)+ # Attention gros cluster de MYTED avec gros score gomm? par ylim
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
  labs(title=titre,y="Total Biomass (gAFDW/m?)") +
  Scale_fill(length(annees)) ;bp #
# ggsave(paste(wdgraph,"Graph ",titre,".png",sep=""), plot = bp, width = 18, height = 6)
# bp2<-bp+ labs(title="",x="") + theme(legend.position='none')
# titre <- paste("Periodic Evolution of Biodiversity",sep="")
# bp <- ggarrange(bp1,bp2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
# bp <- annotate_figure(bp,top = textGrob(titre, gp=gpar(fontsize=15,font=1)),
#                       bottom=xlab)+bgcolor("white")
# ggsave(paste(wdgraph,"Graph ",titre,".png",sep=""), plot = bp, width = 15, height = 5)

titre <- paste("Periodic Evolution of Biomass for chosen species",sep="")
bp<-ggplot(CSLN_Mars_spZ)+geom_boxplot(aes(x=Annee,y=Biomass_gAFDWm2,fill=Period))+ 
  facet_grid(SPCourt~Zone, margins=TRUE)+ # Period, scales="free_y"
  scale_y_continuous(trans = "log10") + 
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) +
  labs(title=titre,y="Biomass (gAFDW/m?)") +
  Scale_fill(length(levels(annees))) ;bp # length(levels(CSLN_Mars_spZ$Period))
# ggsave(paste(wdgraph,"Graph ",titre,".png",sep=""), plot = bp, width = 18, height = 6)

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

CSLN_abd<-abund(CSLN_contZ_mat,f=0.4) #f=0,2 recommande
# Find where our species are in the TPA order to keep them
# intersect(colnames(extract(CSLN_abd,length(CSLN_abd$vr))),species[,1])
# Find the last rank of our species
minSpecies<-max(which(colnames(extract(CSLN_abd,length(CSLN_abd$vr))) %in% species[,1]))
# dev.print(device = png, file = "CSLN_TPA.png", width = 1000, height=600) # If save of graph
plot(CSLN_abd,xlab="Taxon",ylab="Abond rel", dpos=c(40,100),
     main=paste("Methode TPA f=",CSLN_abd$f,sep=""))
# dev.off()
CSLN_abd$n<-max(92,minSpecies) # identify(CSLN_abd) # Fin du plateau
CSLN_Ztpa_SPCourt<-colnames(extract(CSLN_abd,CSLN_abd$n))
CSLN_contZ_tpa<-as.matrix(extract(CSLN_abd,CSLN_abd$n))

# DELETION OF RARE SPECIES IN THE COMPLETE DATABASE -----
# CSLN_Mars <- CSLN_Mars[CSLN_Mars$SPCourt %in% CSLN_Ztpa_SPCourt,]

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
# # suppress station where nothing found in TPA species, required for adonis:
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # WARNING FOLLOWING OPERATION CAN LAST VERY VERY LONG
ad<-adonis2(CSLN_contZ_mat~CSLN_contZ_name$Zone*CSLN_contZ_name$Period*CSLN_contZ_name$Annee,method="bray")
CSLN_Zadonis<-as.data.frame(cbind(test=row.names(ad),ad))
padZ<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$Zone)
padP<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$Period)
padA<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$Annee)
padZP<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$ZP)
padZA<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$ZA)
padZPA<-pairwise.adonis(CSLN_contZ_mat,CSLN_contZ_name$ZPA)
CSLN_Zpairado<-as.data.frame(rbind(padZ,padP,padZP,padZA,padZPA))
CSLN_Zpairado<-CSLN_Zpairado %>% filter(sig != "") # Keep only significant pairs
# # SAVE BECAUSE SO LONG TO CALCULATE
save(ad,CSLN_Zadonis,CSLN_Zpairado,file = paste(wdwork,"CSLN_Mars_GIPSA_Ado.RData", sep=""))
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# load(paste(wdwork,"CSLN_Mars_GIPSA_Ado.RData", sep=""))

#________________________________________________________________
# ANALYSE DES COMMUNAUTES ----
#________________________________________________________________
# Hierarchical Ascending Classification (HAC or dendrogram) ----
  # constructed with the UPGMA agglomeration method to define 'clustchoix' groups 
  # at an equal cut-off level

# SIMPROF Determination of significantly different groups
CSLN_Zsimprof<-simprof(CSLN_contZ_mat, num.expected=1000, num.simulated=999,
                      method.cluster="average",method.distance="braycurtis") # summary(CSLN_simprof)
simprof.plot(CSLN_Zsimprof, leafcolors=NA, plot=TRUE, fill=TRUE,
             leaflab="perpendicular", siglinetype=1)
# dev.print(device = png, file = paste(wdgraph,"CSLN_Mars_Simprof.png",sep=""), width = 1000, height=600)
# dev.off()
# tmp<-
# CSLN_Mars_smZ$Groupe_Simprof[CSLN_Mars_smZ$] CSLN_Zsimprof[["significantclusters"]][[9]]

# CAH manual version 
matdist <-vegdist(CSLN_contZ_mat, method="bray")
dendro <-hclust(matdist,method="average")
plot(dendro, ylab="Bray-Curtis dissimilarity", xlab="Stations",
     main="CAH on Bray-Curtis dissimilarities / Mean link",
     hang = -1)
clustchoix<-3#CSLN_simprof$numgroups
rect.hclust(dendro, k=clustchoix, border="red")
# dev.print(device = png, file = paste(wdgraph,"CSLN_Mars_Dendro.png",sep=""), width = 1000, height=600)
# dev.off()
# Creation of similarity groups
groupe.clust <- cutree(dendro, k=clustchoix)
groupesdend<-cbind(Groupe=row.names(as.data.frame(groupe.clust)),as.data.frame(groupe.clust))[order(groupe.clust),]
# Significance of the groups
adZ_CAH<-adonis2(CSLN_contZ_mat~groupe.clust,method="bray"); adZ_CAH
padZ_CAH<-pairwise.adonis(CSLN_contZ_mat,groupe.clust); padZ_CAH

#________________________________________________________________
# INDVAL ----
  # For the characteristic taxa of each of these groups, the Indval has been 
  # calculated (Species Indicator Values). This index combines the relative 
  # abundance of each species and its occurrence in the samples to define an index 
  # between 0 and 1. The higher the IndVal, the more characteristic the species 
  # is of the community structure of the group.
#________________________________________________________________
CSLN_cont_indval_mat<-CSLN_contZ_mat
CSLN_cont_indval_mat<-subset(CSLN_cont_indval_mat,rowSums(CSLN_cont_indval_mat)!=0)
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
CSLN_contZ_nMDS<-CSLN_contZ
CSLN_contZ_nMDS_mat<-CSLN_contZ_mat

CSLN_Zmds <- metaMDS(CSLN_contZ_nMDS_mat,distance = "bray", k = 2,try = 1000, autotransform =FALSE)
CSLN_Zmds_points <- data.frame(cbind(code=row.names(CSLN_Zmds$points),CSLN_Zmds$points)) %>% 
                    separate(code, c("Zone", "Period","Annee"), sep = "_",remove = FALSE)
rownames(CSLN_Zmds_points) <- c()

CSLN_Zmds_dat <- data.frame(scores(CSLN_Zmds)) %>% # Using the scores function from vegan to extract the code scores
                    mutate(code=row.names(scores(CSLN_Zmds))) %>% relocate(code) %>%
                    separate(code, c("Zone", "Period","Annee"), sep = "_",remove = FALSE) %>%
                    unite(ZP,Zone,Period, sep = "_",remove = FALSE) %>%
                    group_by(Zone,Period) %>% #,Tidal_level
                    mutate(NMDS1.m = mean(NMDS1, na.rm = TRUE),
                           NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_Zmds_dat <- CSLN_Zmds_dat %>% left_join(groupesdend, by = c("code" = "Groupe"))
CSLN_Zmds_dat$groupe.clust<-as.factor(CSLN_Zmds_dat$groupe.clust)
rownames(CSLN_Zmds_dat) <- c()
CSLN_Zmds_mean <- CSLN_Zmds_dat %>% 
  unite("ZP",Zone,Period, sep = "_",remove = FALSE) %>%
  group_by(ZP) %>% 
  summarise(NMDS1.m = mean(NMDS1, na.rm = TRUE),
            NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_Zmds_species<-as.data.frame(CSLN_Zmds$species)
CSLN_Zmds_species<-cbind(sp=rownames(CSLN_Zmds_species),CSLN_Zmds_species)  # create a column of species, from the rownames of CSLN_Zmds_species

nbA = dim(CSLN_Zmds_mean)[1]
# nbB = length(unique(CSLN_Zmds_dat$Period))
palette=colZiss # colDarj colZiss colSpec
bp<-ggplot(data=CSLN_Zmds_dat) +
# Add species name
    geom_point(data=CSLN_Zmds_species,aes(x=MDS1,y=MDS2),shape=8,size=1,colour="gray50") + # add the point markers
    geom_text(data=CSLN_Zmds_species,aes(x=MDS1,y=MDS2,label=sp),
              size=2.5,colour="gray50",check_overlap = TRUE, fontface = "italic") +
# Base scatter plot 
    geom_point(aes(x=NMDS1,y=NMDS2,shape=Period,colour=ZP),size=3) +
    geom_text(aes(x=NMDS1,y=NMDS2,colour=ZP,label=Annee),size=2.5,hjust=-0.3, fontface = "plain") +  # add the site labels
    scale_colour_manual(values=palette(nbA)) +
    stat_ellipse(aes(x=NMDS1,y=NMDS2,group=groupe.clust), #,linetype = groupe.clust,color=groupe.clust
                 type = "norm",level = 0.9,linetype = 3,show.legend = FALSE) + #
# Lines to connect the same Zone and label it
    geom_segment(aes(x=NMDS1.m, y=NMDS2.m, xend=NMDS1, yend=NMDS2, linetype = Period),
                 size=.1,color="grey") +
    annotate(geom="label",x=CSLN_Zmds_mean$NMDS1.m,y=CSLN_Zmds_mean$NMDS2.m,label=CSLN_Zmds_mean$ZP,
             fill=palette(nbA),alpha=0.5) +
    # geom_label(aes(x=CSLN_Zmds_mean$NMDS1.m,y=CSLN_Zmds_mean$NMDS2.m,label=CSLN_Zmds_mean$Zone),
    #            fill = palette(nbA),size = 3,label.size =0,alpha=0.5)+
# Graphic visuals
    annotate("text",x=min(CSLN_Zmds_dat$NMDS1),y=min(CSLN_Zmds_dat$NMDS2),hjust=-0.2,
             label=paste("Stress =",round(CSLN_Zmds$stress,4)), size = 4)+
    ggtitle(label="nMDS CSLN data by Zone")+
    theme(axis.text.x = element_blank(),  # remove x-axis text
          axis.text.y = element_blank(), # remove y-axis text
          axis.ticks = element_blank(),  # remove axis ticks
          plot.background = element_blank())+theme_bw();bp
# ggsave(paste(wdgraph,"nMDS Macrofaune by Zone",".png",sep=""), plot = bp, width = 12, height = 8)


# #________________________________________________________________
# # ANALYSES STATISTIQUES ----
# # COMPARAISON DES Moyennes SPATIALES
# CSLN.zone.stat.s <- data.frame(matrix(NA,ncol = 10, nrow = length(anneesData))) #NULL
# colnames(CSLN.zone.stat.s) <- c("Categorie","Variable","Annee","Mois","Shapiro","A","B","C","D","Kruskal") #NULL
# k=1
# for (j in 1:length(anneesBiv)){
#   df<-with(CSLN_summary,CSLN_summary[Annee==anneesBiv[j],])
#   shap<-by(df$Densite.m2.mea,df$Zone,shapiro.test)#Si p>alpha normalite
#   krus<-kruskal.test(Densite.m2.mea~as.factor(Zone),data=df)#Si p<alpha diff entre groupes
#   CSLN.zone.stat.s[k,]<-c("CSLN","Densite",as.character(anneesData[j]),"Mois","Shapiro",
#                           round(shap$A$p.value,4),round(shap$B$p.value,4),
#                           round(shap$C$p.value,4),round(shap$D$p.value,4),round(krus$p.value,4))
#   k<-k+1
# }
# #glob<-pairwise.wilcox.test(df$Densite.m2,df$Zone,p.adj="bonf")

#________________________________________________________________
# SHP CREATION OF FAUNA & MARS DATA ----
CSLN_sf <- st_as_sf(CSLN_pur,coords= c("longitude","latitude"),crs=4326,remove = FALSE) #c("x","y"),crs=2154,remove = FALSE) %>% st_transform(4326) #"Lon","Lat" in WGS=4326
Sel_col<-c(facto_col,"Density_indm2","Biomass_gAFDWm2",
           "IndBodySize_gAFDW","MSR_mW","MSRtot") #,"Itot"
CSLN_sf <- CSLN_sf %>% select(all_of(Sel_col)) # lighten the data to help build the map
# GIS HTML MAP CREATION ----
# https://thinkr.fr/sil-te-plait-dessine-moi-carte-r/
# https://thinkr.fr/cartographie-interactive-comment-visualiser-mes-donnees-spatiales-de-maniere-dynamique-avec-leaflet/
# https://thinkr.fr/cartographie-interactive-avec-r-la-suite/
# https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html
# https://bookdown.org/nicohahn/making_maps_with_r5/docs/tmap.html
#________________________________________________________________
# boxbds=c(xmin=-0.1, ymin=49.3, xmax=0.45, ymax=49.7)
# bds<-st_crop(ES_Areas,boxbds, crs = st_crs(4326)) %>% # decoupe zone interet 
#   st_transform(2154) # transform to planar as required by st_intersection() 
# bds_points <- cbind(bds, st_coordinates(st_centroid(bds$geometry))) # def centroides pour placer les noms
# # pol1 = st_polygon(list(rbind(c(0,0),c(1,0),c(1,1),c(0,1),c(0,0))))

tmap_mode("view") #tmap_mode("plot") # ttm() toggle #
tm_Bio<-tm_basemap("OpenStreetMap.HOT") + # OpenStreetMap.HOT .Mapnik .France - Stamen.Watercolor #https://leaflet-extras.github.io/leaflet-providers/preview/
  # tm_logo(paste(wdlogos,"logo-MIE.png",sep=""), height = 2) +
  # tm_logo(c(paste(wdlogos,"logoMP.png",sep=""),
  #           paste(wdlogos,"Logo Borea.png",sep="")), height = 2) +
  tm_scale_bar(position = c("left", "bottom"), width = 0.15)+ #SCALE
  tm_compass(position = c("left", "top"), size = 2)+          #NORTH COMPASS
  tm_shape(ES_Areas) +
  tm_fill(col = "Zone", palette = "Spectral", alpha = 0.6) +
  tm_borders("white", lwd = 1) +
  tm_shape(CSLN_sf) +
  # tm_bubbles(size="Density_indm2") +
  tm_dots(size=0.001) +
  tm_layout(legend.outside = TRUE)
# tm_Bio # WARNING THIS OPERATION CAN LAST LONG
tmlf_Bio<-tmap_leaflet(tm_Bio) # conversion to leaflet object, quicker??
# tmlf_Bio

# # Simple and basic map of France
# france <- ne_states(country = "France", returnclass = "sf") %>%
#   filter(region %in% c("Haute-Normandie","Basse-Normandie"))
# map_france<-tm_shape(france) + tm_polygons()#+
# # tm_shape(mybb)+tm_borders("white", lwd = 1) # ajout du cadre rouge sur la zone interet
# print(map_france, vp = grid::viewport(0.9, 0.7, width = 0.2, height = 0.2))
# # map_france<-qtm("france") # creation tres rapide d'une carte

# _______________________________________________________________
# OUTPUT SAVE ----
  wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep=""))
  writeData(wb, sheet = "CSLN", x = CSLN, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "CSLN_Mars", x = CSLN_Mars, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Sumry_ZTYS", x = Sumry_ZTYS, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Varnames", x = Varnames, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Predicteurs", x = predict, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "MEVE", x = MEVE, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Saison", x = saison, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Reponse", x = reponse, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Species", x = species, startCol = 1, startRow = 1,withFilter = FALSE)
  saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)
  
  save.image(file = paste(wdwork,"CSLN_Mars_BDD_GIPSA",".RData", sep=""))
  
# GIS SAVE
  write.csv(CSLN_Mars,file=paste(wdres,"CSLN_Mars", ".csv",sep=""), na = "",row.names = FALSE)
  write.csv(CSLN_pur,file=paste(wdres,"CSLN_data", ".csv",sep=""), na = "",row.names = FALSE)
  # WARNING FOLLOWING OPERATION CAN LAST LONG
  tmap_save(tm_Bio, filename = paste(wdres,'CSLN_Mars Map.html',sep=""))
  # saveWidget(tmlf_Bio, paste(wdres,'CSLN_Mars Map.html',sep=""), selfcontained = TRUE)
  
#________________________________________________________________
# # OPTION LEAFLET----
# pal <- colorFactor(
#   palette = "viridis",na.color = NA,
#   levels = factor(bds$Type))
# map_leaflet <- leaflet() %>%
#   addProviderTiles("OpenStreetMap.HOT") %>% #addTiles()
#   # setView(lng = 2.80, lat = 46.80, zoom = 5) %>%
#   # addMarkers(data = CSLN_sf) %>%  # ATTENTION HYPER LONG !!! 
#   addPolygons(data = bds,
#               label = ~Type, # En passant la souris
#               popup = ~Zone, # En cliquant sur l'icone
#               fill = TRUE, 
#               fillColor = ~pal(Type),
#               fillOpacity = 0.8,
#               highlightOptions = highlightOptions(color = "white", weight = 2)) %>% 
#    addRectangles(
#     lng1 = boxbds[1], lat1 = boxbds[2],
#     lng2 = boxbds[3], lat2 = boxbds[4],
#     color = "green",
#     fill = FALSE) %>%
#   addLegend(
#     title = "Zones",
#     pal = pal, values = bds$Type)
# map_leaflet
# saveWidget(map_leaflet, 'test_leaflet.html', selfcontained = TRUE)
