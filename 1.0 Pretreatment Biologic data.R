#_______________________________________________________________________________
# TREATMENT OF BIOLOGIC DATA FROM CSLN BEFORE SDM TREATMENT
# Amelie LEHUEN Janvier 2022
#_______________________________________________________________________________

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
library(pastecs) ; library(clustsig) ; library(dendextend) # # Analyse TPA (abund); Analyse SIMPROF (simprof); library(ggdendro) ggdendrogramm# Graphics packages
library(RColorBrewer) ; library(wesanderson) # Palettes de couleurs
library(ggpubr);library(gridExtra) ; library(grid) # Mozaic of graphs tools
library(treemap) # pour les graphiques treemap
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
wdmsr <- (paste(wdscript,"MSR/MSR.R",sep=""))
# setwd(paste(wdtask,"Scripts/",sep=""))
setwd(wdtask)
# If exists
# load(file = paste(wdwork,"CSLN_Mars_BDD",".RData", sep=""))

#________________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colDarj <- function(x) {wes_palette("Darjeeling2",x, type = "continuous")}
colZiss <- function(x) {wes_palette("Zissou1",x, type = "continuous")}
colSpec <- colorRampPalette(brewer.pal(8, "Spectral")); 
colDark <- colorRampPalette(brewer.pal(8, "Dark2"));
Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
Scale_fill <- function(x) {scale_fill_manual(values=colDarj(x))}

# ________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
# to be upgraded
species <-data.frame(SPCourt=c("CERED","CORVO","HEDDI","MACBA","PERUL","SCRPL"),
                     Taxon_SNa=c("Cerastoderma edule","Corophium volutator","Hediste diversicolor","Macoma balthica","Peringia ulvae","Scrobicularia plana"))
speciesB <-data.frame(SPCourt=c("AREMA","BATPI","BATSA","CYACA","ETELO","NEPCI","NEPHO","PYGEL"),
                     Taxon_SNa=c("Arenicola marinea","Bathyporea pilosa","Bathyporea sarsi","Cyathura carinata","Eteone longa","Nephtys cirrosa","Nephtys hombergii","Pygospio elegans"))
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
CSLN$SP<-"OTHER"; CSLN$SP[CSLN$SPCourt %in% species$SPCourt]<-"SpCh"; CSLN$SP[CSLN$SPCourt %in% speciesB$SPCourt]<-"SpMSR";
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
for (i in which(!is.na(CSLN_Mars$temp_m) & CSLN_Mars$SP!="OTHER")){
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
CSLN_Stations <- CSLN_Mars %>% select(Station_originelle,Filtre,Zone,Period,longitude,latitude) %>% 
  group_by(Station_originelle,Filtre,Zone,Period) %>%
  summarise(longitude=mean(as.numeric(longitude),na.rm=TRUE),latitude=mean(as.numeric(latitude),na.rm=TRUE)) %>%
  unite(SP, Station_originelle,Period,remove=FALSE)

CSLN_pur<-CSLN_Mars # CSLN_Mars<-CSLN_pur #summary(CSLN_pur) # Save before more modifications
# summary(CSLN_Mars)

#________________________________________________________________
# FAUNA SUMMARY BY ZONE : SPECIES CHOSEN
CSLN_Mars_sp <- CSLN_Mars %>% 
  select(Zone,Tidal_level,Period,Annee,Season,SP,SPCourt,Density_indm2,Biomass_gAFDWm2,MSRtot) %>% # 
  filter(SP=="SpCh") %>% select(-SP) %>%
  complete(nesting(Zone,Tidal_level,Period,Annee,Season),nesting(SPCourt), # 
           fill=list(Biomass_gAFDWm2=0,Density_indm2=0))
CSLN_Mars_spm <- CSLN_Mars_sp %>% 
  group_by(Zone,Tidal_level,Period,Annee,Season,SPCourt) %>% # 
  summarise(Biomass_m=median(Biomass_gAFDWm2,na.rm=TRUE),Biomass_sd=sd(Biomass_gAFDWm2,na.rm=TRUE),
            Density_m=median(Density_indm2,na.rm=TRUE),Density_sd=sd(Density_indm2,na.rm=TRUE),
            MSRtot_m=median(MSRtot,na.rm=TRUE),MSRtot_sd=sd(MSRtot,na.rm=TRUE))
# DIVERSITY TABLE BY ZONE
CSLN_Mars_sm <- CSLN_Mars %>%
  select(Zone,Tidal_level,Period,Annee,Season,SPCourt,Density_indm2,Biomass_gAFDWm2,MSRtot) %>% # 
  group_by(Zone,Tidal_level,Period,Annee,Season) %>% # 
  summarise(SR=n_distinct(SPCourt), n_records=n(),
            Biomass_t=sum(Biomass_gAFDWm2,na.rm=TRUE),
            MSRtot_t=sum(MSRtot,na.rm=TRUE)) %>%
  unite("code",Zone,Tidal_level,Period,Annee,Season,sep = "_",remove = FALSE) # 
# CONTINGENCY TABLE
CSLN_cont <- CSLN_Mars %>%
  select(Zone,Tidal_level,Period,Annee,Season,SPCourt,Density_indm2) %>% # 
  complete(nesting(Zone,Tidal_level,Period,Annee,Season),nesting(SPCourt), # 
           fill=list(Density_indm2=0)) %>%
  group_by(Zone,Tidal_level,Period,Annee,Season,SPCourt) %>% # 
  summarise(Density_m=median(Density_indm2,na.rm=TRUE)) %>%
  pivot_wider(names_from = SPCourt, values_from = Density_m, values_fill = 0) %>%
  unite("code",Zone,Tidal_level,Period,Annee,Season,sep = "_",remove = TRUE) # 
CSLN_cont_mat<-as.matrix(CSLN_cont[,-1])
row.names(CSLN_cont_mat)<-CSLN_cont$code

# SHANNON & PIELOU INDEXES ----
# Pielou's index of equitability (J ): normalization of the Shannon-Wiener index (H'), 
# a value of taxonomic diversity as a function of the number of taxa per area and the 
# abundance of individuals within each taxon; 0 means that one taxon dominates the others,
# 1 means that there is an equitable distribution of individuals between taxa
CSLN_Mars_sm$Shannon<-diversity(CSLN_cont_mat,index="shannon")
CSLN_Mars_sm$Pielou<-CSLN_Mars_sm$Shannon/log2(CSLN_Mars_sm$SR)

CSLN_cont_mat<-subset(CSLN_cont_mat,rowSums(CSLN_cont_mat)!=0)
CSLN_cont_name<-data.frame(code=row.names(CSLN_cont_mat))  %>% 
  separate(code, c("Zone", "Tidal_level", "Period", "Annee","Season"), sep = "_",remove = TRUE) %>% # 
  unite("ZP",Zone,Period,sep = "_",remove = FALSE) %>% 
  unite("ZT",Zone,Tidal_level,sep = "_",remove = FALSE) %>% 
  unite("ZA",Zone,Annee,sep = "_",remove = FALSE) %>%
  unite("ZS",Zone,Season,sep = "_",remove = FALSE) %>%
  unite("PT",Period,Tidal_level,sep = "_",remove = FALSE) %>%
  unite("PS",Period,Season,sep = "_",remove = FALSE) %>%
  unite("ZPA",Zone,Period,Annee,sep = "_",remove = FALSE) %>%  
  unite("ZTP",Zone,Tidal_level,Period,sep = "_",remove = FALSE) %>%
  unite("ZPS",Zone,Period,Season,sep = "_",remove = FALSE) %>%  
  unite("ZTP",Zone,Tidal_level,Season,sep = "_",remove = FALSE) %>%
  unite("ZTPA",Zone,Tidal_level,Period,Annee,sep = "_",remove = FALSE)
  unite("ZTPS",Zone,Tidal_level,Period,Season,sep = "_",remove = FALSE)

# save.image(file = paste(wdwork,"CSLN_Mars_BDD",".RData", sep=""))
  
# _______________________________________________________________
# OUTPUT SAVE ----
  wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep=""))
  writeData(wb, sheet = "CSLN", x = CSLN, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "CSLN_Mars", x = CSLN_Mars, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "CSLN_Mars_sm", x = CSLN_Mars_sm, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Varnames", x = Varnames, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Predicteurs", x = predict, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Saison", x = saison, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Reponse", x = reponse, startCol = 1, startRow = 1,withFilter = FALSE)
  writeData(wb, sheet = "Species", x = species, startCol = 1, startRow = 1,withFilter = FALSE)
  saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)

  save.image(file = paste(wdwork,"CSLN_Mars_BDD",".RData", sep=""))
  
  write.csv(CSLN_Stations,file=paste(wdres,"CSLN_Stations", ".csv",sep=""), na = "",row.names = FALSE)
  