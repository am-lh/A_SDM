#######SCRIPT MAXIME ######## ----

# Environnement de travail----

pc <- "C:/Users/maxcb/OneDrive/"  
wdtask <- paste(pc,"Bureau/StageM1/BDD/",sep="")
wdsource <- paste(wdtask,"Sources/Faune/CSLN/",sep="")
wdwork <- paste(wdtask,"Matrices/",sep="")
wdgraph <- paste(wdtask,"Graphiques/",sep="")
wdres <- paste(wdtask,"Resultats/",sep="")
wdGIS <- paste(pc,"Bureau/StageM1/SIG/",sep="");
wdscript <- (paste(pc,"Bureau/StageM1/Scripts/",sep=""))
wdmsr <- (paste(wdscript,"MSR/MSR.R",sep=""))
# setwd(paste(wdtask,"Scripts/",sep=""))
setwd(wdtask) 
#load(file = paste(wdwork,"BBDD_CSLN_Annee_periode_zone_nt",".RData", sep = ""))

# load(file = paste(wdwork,"Script Mudflat",".RData", sep = ""))

# Library ----
library(readxl); library(openxlsx)
library(tidyverse)
library(dplyr)
#install.packages("tidyverse")
library(vegan)
library(pastecs) #CAH
library(ggdendro)
library(clustsig) #SIMPROF
library(labdsv)

#GIS 
library(sf)
library(sfheaders)
#install.packages("devtools"); library(devtools);Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#Couleur
library(viridis)
library(RColorBrewer)  #display.brewer.all()
# Utilisation dans ggplot : 
# scale_fill_brewer() pour box plot, bar plot, violin plot, dot plot, etc
# scale_color_brewer() pour les lignes et les points
#  brewer.pal() est utilisée pour générer un vecteur de couleurs BARPLOT
# https://www.datanovia.com/en/fr/blog/top-palettes-de-couleurs-r-a-connaitre-pour-une-meilleur-visualisation-des-donnees/



# Chargement des données----

fauna_file <- paste(wdsource,"CSLN_Biology_source.xlsx",sep="")
CSLN.data <- as.data.frame(read_excel(fauna_file,sheet = "Biology_station", na = ""))
granulo <- as.data.frame(read_excel(fauna_file,sheet = "Granulo", na = ""))

# Mise en forme des données ----
#Ajout granulo
CSLN<-CSLN.data %>% left_join(granulo[,-c(21:24)]) %>% 
  dplyr::arrange(idStationUnique, Taxon) %>%dplyr::mutate(IndBodySize_gAFDW= 
                                                     Biomass_gAFDWm2/Density_indm2)#calcul de la densite INTERMUD en fonction de l'abondance
CSLN$IndBodySize_gAFDW[CSLN$Density_indm2==0]<-0

# ajout  Biomasse sur Intermud
tmp<-CSLN%>% filter(Filtre=="INTERMUD") %>% 
  select(ScientificName_accepted) %>% left_join(CSLN %>% group_by(ScientificName_accepted) %>% 
                                                  summarise(Mediane.IndBodySize = median
                                                            (IndBodySize_gAFDW,na.rm = T))) %>%
  select(Mediane.IndBodySize)

CSLN$Biomass_gAFDWm2[CSLN$Filtre=="INTERMUD"]<-as.vector(tmp$Mediane.IndBodySize)

# Mise en forme longitude latitude
CSLN<-CSLN %>% 
  filter(longitude!="NaN",latitude!="NaN", 
         !is.na(Density_indm2),!is.na(Biomass_gAFDWm2),!is.na(AphiaID_accepted),
         (Biomass_gAFDWm2!=0 & Density_indm2!=0 )| 
           (Biomass_gAFDWm2==0 & Density_indm2==0 ))

# Supprimer les doublons de IdstationUnique, SPCourt
CSLN<-CSLN %>% distinct(idStationUnique, SPCourt, .keep_all = TRUE)

# Suppression de colonnes
suppr_col<-c("Source","Site","Campagne","Jour","id_station","Engin",
             "maille_tamis","Forme_maille","nb_replicat_protocole","nb_replicat",
             "surface_engin","surface_station","Abondance","PSLCUnitaire",
             "identite_tri","identite_det","Operateur","commentaire","Methode")
CSLN<-CSLN %>% select(-all_of(suppr_col))

# Definition des especes interessantes
CSLN$specix<-"Autre especes"
CSLN$specix[CSLN$ScientificName_accepted %in% c("Eteone longa", "Crangon crangon", "Oligochaeta", "Carcinus maenas", "Mesopodopsis slabberi",
                                                "Nephtys cirrosa", "Macoma balthica"  ,"Pholoe baltica","Nephtys hombergii", 
                                                "Kirkegaardia", "Glycera tridactyla", "Nemertea" ,"Fabulina fabula",
                                                "Mediomastus fragilis" , "Corophium volutator", "Cerastoderma edule", "Abra alba",
                                                "Petricolaria pholadiformis", "Donax vittatus", "Pygospio elegans" ,
                                                "Gastrosaccus spinifer", "Scrobicularia plana","Bathyporeia pilosa",
                                                "Spio martinensis",  "Nucula nitidosa",  "Hediste diversicolor", "Kurtiella bidentata",    
                                                "Ensis leei"             ,    "Streblospio benedicti",      "Spisula subtruncata",       
                                                "Polydora ciliata"        ,   "Diastylis bradyi"      ,     "Cyathura carinata"   ,      
                                                "Peringia ulvae"           ,  "Phaxas pellucidus"      ,    "Alitta succinea"      ,     
                                                "Nephtys caeca"             , "Haustorius arenarius"    ,   "Magelona johnstoni"    ,    
                                                "Acrocnida brachiata"      ,  "Mya arenaria"           ,   
                                                "Varicorbula gibba"          ,"Tritia reticulata"         , "Schistomysis kervillei"  ,  
                                                "Tharyx"            ,         "Chaetozone"                 ,"Nephtys kersivalensis"    , 
                                                "Idotea linearis"    ,        "Copepoda"               ,    "Mactra stultorum"          ,
                                                "Gammarus salinus"    ,          "Magelona filiformis"       ,
                                                "Asterias rubens"      ,      "Spio decorata"            ,  "Spiophanes bombyx"         ,
                                                "Glycinde nordmanni"    ,     "Phyllodoce mucosa"         , "Sphaeroma serratum"        ,
                                                "Ophiura ophiura"        ,    "Achelia hispida"     )]<- "Les especes +++"

# suprression d'especes

CSLN<- CSLN %>% filter(SPCourt!="SEMBA",SPCourt!="AMPIM", SPCourt!="AUSMO",SPCourt!="BALCR",SPCourt!="MYTED", SPCourt!="BIVAL", SPCourt!="ANNEL") 
CSLN<- CSLN %>% filter(!grepl("débris",ScientificName_accepted))

# Definition par periode 1996-1999/2000-2010/2010-2018

CSLN$periode<-NA 
CSLN$periode[CSLN$Annee <2000 ]<-"1996-1999" 
CSLN$periode[CSLN$Annee >=2000 & CSLN$Annee <= 2005]<-"2000-2005"
CSLN$periode[CSLN$Annee > 2005 & CSLN$Annee <= 2010]<-"2006-2010"
CSLN$periode[CSLN$Annee > 2010 & CSLN$Annee <=2015]<-"2011-2015"
CSLN$periode[CSLN$Annee > 2015 & CSLN$Annee <= 2019]<-"2015-2019"
CSLN <- CSLN %>% filter(periode!= "1996-1999")

CSLN$periode<-as.factor(CSLN$periode)

#Definition des niveaux tidaux
# CSLN$niv.tidal<-NA 
# CSLN$niv.tidal[CSLN$inunt >= 0 & CSLN$inunt<= 0.25]<-"Supratidal"
# CSLN$niv.tidal[CSLN$inunt >0.25 & CSLN$inunt<= 0.75]<-"Intertidal"
# CSLN$niv.tidal[CSLN$inunt > 0.75 & CSLN$inunt<= 1]<-"Infratidal"
# CSLN$niv.tidal<-as.factor(CSLN$niv.tidal)

#Completer les niveaux tidaux manquant
# CSLN[CSLN$NINJ %in% CSLN$NINJ[which(CSLN$niv.tidal=="Intertidal")],
#        "niv.tidal"]<-"Intertidal"
# CSLN[CSLN$NINJ %in% CSLN$NINJ[which(CSLN$niv.tidal=="Infratidal")],
#        "niv.tidal"]<-"Infratidal"
# CSLN[CSLN$NINJ %in% CSLN$NINJ[which(CSLN$niv.tidal=="Supratidal")],
#        "niv.tidal"]<-"Supratidal"

#Suppression des NA par niveau tidal et zone
CSLN<-CSLN %>% filter(!is.na(Zone)) #& !is.na(niv.tidal))



#Données GEOgraphique----

CSLN.geo<-CSLN  %>% st_as_sf(coords = c("longitude","latitude"), #%>% filter(Zone=="North Downstream Mudflat", "North Median Mudflat", "North Upstream Mudflat", "South Mudflat")
                            crs=4326) %>% st_transform(2154) 
Coord.uniq_CSLN<-CSLN %>% distinct(longitude,latitude) %>% 
  st_as_sf(coords= c("longitude","latitude"),crs=4326) %>% 
  st_transform(2154) # transformer WGS en lambert
station_uniq<- CSLN %>% filter(grepl("Mudlfat",Zone)) %>% distinct(Station_originelle) 


#GIS Integration Mars dans les données biologiques

Mars.csv<-read.csv("C:/Users/maxcb/OneDrive/Bureau/StageM1/BDD/Resultats/ES_Mars_Maps_sh.csv")
Mars.csv[Mars.csv== "NaN"]<-NA 
Mars.shp<-Mars.csv %>% filter(!is.na(Lon),!is.na(Lat),!is.na(flow_m)) %>% 
  st_as_sf(coords = c("Lon","Lat"),crs=4326) %>% st_transform(2154)
CSLN.geo$NINJ<-Mars.shp$NINJ[st_nearest_feature(CSLN.geo,Mars.shp)] #ajoute le NINJ au CSV

#Integration de la zone dans les données biologiques
Gis.shp<-read_sf("C:/Users/maxcb/OneDrive/Bureau/StageM1/SIG/Layers_made/ES_Areas_WGS.shp")
Gis.shp<-Gis.shp%>%  st_transform(2154) %>% filter(Zone!="NA" & Zone!="Bay"  & 
                                                     Zone!="OffShore" & Zone!="Octeville" & 
                                                     Zone!="Cote Fleurie" &
                                                     Zone!="Ilot Oiseaux") #grepl("Mudlaf", Zone)

# Intersection entre zone et coordonnée des points de prélèvements pour identifier dans quelle zone sont les points
CSLN.geo<- CSLN.geo %>% mutate(intersection  =as.integer(st_intersects(geometry,Gis.shp)),
                               Zone = if_else(is.na(intersection),"NA", Gis.shp$Zone[intersection]),
                               Type = if_else(is.na(intersection),"NA", Gis.shp$Type[intersection]))
CSLN.geo$Zone[CSLN.geo$Zone=="NA"]<-NA; CSLN.geo$Type[CSLN.geo$Type=="NA"]<-NA #trnasformation du champs NA en donnée manquante
CSLN.geo<- CSLN.geo %>% filter(Zone!="NA" & Zone!="Bay"  & 
                                     Zone!="OffShore" & Zone!="Octeville" & 
                                     Zone!="Cote Fleurie" &
                                     Zone!="Ilot Oiseaux") #grepl("Mudlaf", Zone)
CSLN<-sf_to_df(CSLN.geo, fill=TRUE)
CSLN$Zone<-as.factor(CSLN$Zone)

#réaliser une jointure entre le Mars.CSV et forme 2 qui ont un champs en commun le NINJ et l'année
CSLN<-CSLN %>% left_join(Mars.csv, by = c("Annee","NINJ"))

#Representation en carte 
# plot(Coord.uniq_CSLN)
# plot(Gis.shp$geometry)
# carte<-ggplot(Gis.shp)+geom_sf(aes(fill=Zone))+ 
#   geom_sf(data=Coord.uniq_CSLN) ; carte
#  
# fr <- c(left = -3, bottom = 49, right =0 , top = 50)
# get_stamenmap(fr, zoom = 5,"toner-lite") %>% ggmap()

# library(leaflet)
# 
# m <- leaflet() %>%
#   addTiles() %>%  # Add default OpenStreetMap map tiles
#   addMarkers(lng=0.12641, lat=49.42174, popup="The birthplace of R") ; m 
# 
# devtools::install_github("ropensci/rnaturalearth")
# devtools::install_github("ropensci/rnaturalearthdata")
# devtools::install_github("ropensci/rnaturalearthhires")

# boxbds=c(xmin=-0.1, ymin=49.3, xmax=0.45, ymax=49.7)
# bds<-st_crop(ES_Areas,boxbds, crs = st_crs(4326)) %>% # decoupe zone interet
#   st_transform(2154) # transform to planar as required by st_intersection()
# bds_points <- cbind(bds, st_coordinates(st_centroid(bds$geometry))) # def centroides pour placer les noms
# # pol1 = st_polygon(list(rbind(c(0,0),c(1,0),c(1,1),c(0,1),c(0,0))))
library(tmap)
tmap_mode("view") 
#tmap_mode("plot") # ttm() toggle #
tm_Bio<-tm_basemap("OpenStreetMap.HOT") + # OpenStreetMap.HOT .Mapnik .France - Stamen.Watercolor #https://leaflet-extras.github.io/leaflet-providers/preview/
  # tm_logo(paste(wdlogos,"logo-MIE.png",sep=""), height = 2) +
  # tm_logo(c(paste(wdlogos,"logoMP.png",sep=""),
  #           paste(wdlogos,"Logo Borea.png",sep="")), height = 2) +
  tm_scale_bar(position = c("left", "bottom"), width = 0.15)+ #SCALE
  tm_compass(position = c("left", "top"), size = 2)+          #NORTH COMPASS
  tm_shape(Gis.shp) +
  tm_fill(col = "Zone", palette = "Spectral", alpha = 0.6) +
  tm_borders("white", lwd = 1) +
  tm_shape(CSLN.geo) +
  # tm_bubbles(size="Density_indm2") +
  tm_dots(size=0.001) +
  tm_layout(legend.outside = TRUE)
tm_Bio # WARNING THIS OPERATION CAN LAST LONG
# tmlf_Bio<-tmap_leaflet(tm_Bio) # conversion to leaflet object, quicker??
#  tmlf_Bio
# library(rnaturalearth)
# # Simple and basic map of France
# france <- ne_states(country = "France", returnclass = "sf") %>%
#   filter(region %in% c("Haute-Normandie","Basse-Normandie"))
# map_france<-tm_shape(france) + tm_polygons()+
#  tm_shape(mybb)+tm_borders("white", lwd = 1) # ajout du cadre rouge sur la zone interet
# 
# print(map_france, vp = grid::viewport(0.9, 0.7, width = 0.2, height = 0.2))
# map_france<-qtm("france") # creation tres rapide d'une carte

#ggsave(file = paste(wdgraph,"Carte.png", sep=""),carte)




# Calcul des indices BTA ----
  
  #taux METABOLIQUE #BREY (with yearly med temp for the moment)

source(wdmsr)
CSLN$MSR_mW<-NA
setwd(wdscript) # when the package is done, with repertory in it, no more needed

for (i in which(!is.na(CSLN$temp_m) & CSLN$specix=="Les especes +++")){
  CSLN$MSR_mW[i]<-msr(CSLN$Density_indm2[i],CSLN$Biomass_gAFDWm2[i],
                      CSLN$AphiaID_accepted[i],CSLN$temp_m[i],1,"w")*1000 #Mean MSR (mW/ind)
}
setwd(wdtask) # when the package is done, with repertory in it, no more needed
CSLN$MSRtot<-CSLN$MSR_mW*CSLN$Density_indm2

# CSLN$Itot<-CSLN$MSRtot

Tb.MSR<- CSLN %>%dplyr::select(Annee,periode,Zone, idStationUnique, Station_originelle, ScientificName_accepted,
                         SPCourt, Density_indm2,Biomass_gAFDWm2, specix, MSR_mW, MSRtot) %>% 
  filter((specix== "Les especes +++"),grepl("Mudflat",Zone)) #,!is.na(MSRtot))

Tb.MSR2<- Tb.MSR %>% group_by(idStationUnique, Zone, Annee, periode)  %>% summarise(MSR.an.station= sum(MSRtot,na.rm=T))

Tb.MSR3<- Tb.MSR2 %>% group_by(Annee,periode,Zone)%>% summarise(MSR =mean(MSR.an.station, na.rm=T)) %>% filter(periode!="2000-2005")

library(ggplot2)
graph.msr<-ggplot(Tb.MSR3)+ geom_boxplot(aes(y= MSR, fill= periode)) +coord_cartesian( ylim = c(0, 100)) +
   facet_grid( .~Zone,scales = "free_y") +labs(title="(A)", y="MSR (mW/m²)") +
theme(axis.text.x = element_text(size=0))  ;graph.msr 

  # Potentiel de bioturbation (Solan & quieros)
potentiel.de.bioturbation <- read.csv("C:/Users/maxcb/OneDrive/Bureau/potentiel de bioturbation.csv", sep=";")
CSLN<- CSLN %>% left_join(potentiel.de.bioturbation)
CSLN$solan <-sqrt(CSLN$Biomass_gAFDWm2/CSLN$Density_indm2)*CSLN$Density_indm2*CSLN$Mi*CSLN$Ri

pbiot<- CSLN %>% dplyr::select(idStationUnique, Annee, periode, Zone, ScientificName_accepted, Density_indm2, 
                               Biomass_gAFDWm2, Ri,Mi, specix) %>% filter(specix=="Les especes +++",Zone!="Channel")
pbiot$solan<- sqrt(pbiot$Biomass_gAFDWm2/pbiot$Density_indm2)*pbiot$Density_indm2*pbiot$Mi*pbiot$Ri
Tb.pBiot<- pbiot %>% group_by(idStationUnique, Annee, periode, Zone) %>% summarise(pbiot.an.station= sum(solan,na.rm=T)) 
Tb.pBiot2<- Tb.pBiot %>% group_by(Annee, periode,Zone)%>% summarise(Potentiel_de_bioturbation=mean(pbiot.an.station, na.rm=T))


graph.pbiot<-ggplot(Tb.pBiot2)+ geom_boxplot(aes(y= Potentiel_de_bioturbation, fill= periode)) + #coord_cartesian( ylim = c(0, 100)) +
  facet_grid( .~Zone,scales = "free_y") + labs(title="(B)", y="Potentiel de bioturbation") + theme(axis.text.x = element_text(size=0))  ;graph.pbiot

library(ggpubr)
ggarrange(graph.msr,graph.pbiot, common.legend = TRUE, legend = "bottom")

# statistique BTA----

CSLN.BTA <- CSLN %>% filter(periode!="2000-2005", Zone!="Channel", ScientificName_accepted=="Hediste diversicolor"|ScientificName_accepted =="Scrobicularia plana" |
                              ScientificName_accepted =="Macoma balthica"|ScientificName_accepted =="Cerastoderma edule" |ScientificName_accepted =="Peringia ulvae" |ScientificName_accepted =="Corophium volutator",!is.na(MSRtot), !is.na(solan),solan<8000, MSRtot<75) %>%
  dplyr::select(idStationUnique, Annee, periode, Zone,ScientificName_accepted, Density_indm2, Biomass_gAFDWm2, specix, MSRtot, solan)

    #MSR
a.MSR <- aov(MSRtot~periode*Zone, data = CSLN.BTA) #Tb.MSR3
shapiro.test(a.MSR$residuals)
ggqqplot(residuals(a.MSR))
bartlett.test(MSRtot~periode, data = CSLN.BTA)
bartlett.test(MSR~Zone, data = Tb.MSR3)
msr_con<-as.matrix(Tb.MSR3[,-c(1:3)])
row.names(msr_con)<-paste(Tb.MSR3$periode, Tb.MSR3$Annee, Tb.MSR3$Zone, sep = "_")

# a.BP <- aov(solan~periode*Zone, data = CSLN.BTA)
# shapiro.test(a.BP$residuals)
# ggqqplot(residuals(a.MSR))
# bartlett.test(solan~periode, data = CSLN.BTA)
# bartlett.test(solan~Zone, data = CSLN.BTA)
# BP_con<-as.matrix(Tb.pBiot2[,-c(1:3)])
# row.names(BP_con)<-paste(Tb.pBiot2$periode, Tb.pBiot2$Annee, Tb.pBiot2$Zone, sep = "_")


msr_con<-subset(msr_con,rowSums(msr_con)!=0) #
msr_con_names<-data.frame(Names=row.names(msr_con)) %>%
  separate(Names,c("periode", "Annee","Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F) %>%
  unite("ZA",Zone,Annee,sep = "_", remove = F) %>%
  unite("ZPA",Zone,periode,Annee,sep = "_", remove = F) #  %>%
# unite("ZNT",Zone,niv.tidal,sep = "_", remove = F) %>%
# unite("ZNTP",Zone,niv.tidal, periode,sep = "_", remove = F) %>%
# unite("ZNTA",Zone,niv.tidal,Annee,sep = "_", remove = F) %>%
# unite("ZNTPA",Zone,niv.tidal, periode, Annee,sep = "_", remove = F)
# CSLN_cont.mat<-CSLN_cont.mat %>% filter(grepl(Code,"North Median Mudflat"))

perm<- adonis2(Tb.MSR3~Zone)
pairwise.adonis(Tb.MSR3,Tb.MSR3$periode)

djk<-adonis2(msr_con ~ msr_con_names$periode*msr_con_names$Zone,  #*CSLN_cont.mat_names$Annee, 
             method = "bray", permutations=999)
MSR_padZ<-pairwise.adonis(msr_con,msr_con_names$Zone)
MSR_padP<-pairwise.adonis(msr_con,msr_con_names$periode)
# CSLN_padA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$Annee)
#CSLN_padNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$niv.tidal)
MSR_padZP<-pairwise.adonis(msr_con,msr_con_names$ZP)
MSR_tabAD<-as.data.frame(rbind(MSR_padZ,MSR_padP,MSR_padZP)) %>% filter(sig!="")

    #potentiel de bioturbation

a.BP <- aov(solan~periode*Zone, data = CSLN.BTA)
shapiro.test(a.BP$residuals)
ggqqplot(residuals(a.MSR))
bartlett.test(solan~periode, data = CSLN.BTA)
bartlett.test(solan~Zone, data = CSLN.BTA)
BP_con<-as.matrix(Tb.pBiot2[,-c(1:3)])
row.names(BP_con)<-paste(Tb.pBiot2$periode, Tb.pBiot2$Annee, Tb.pBiot2$Zone, sep = "_")


BP_con<-subset(BP_con,rowSums(BP_con)!=0) #
BP_con_names<-data.frame(Names=row.names(BP_con)) %>%
  separate(Names,c("periode", "Annee","Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F) %>%
  unite("ZA",Zone,Annee,sep = "_", remove = F) %>%
  unite("ZPA",Zone,periode,Annee,sep = "_", remove = F) #  %>%
# unite("ZNT",Zone,niv.tidal,sep = "_", remove = F) %>%
# unite("ZNTP",Zone,niv.tidal, periode,sep = "_", remove = F) %>%
# unite("ZNTA",Zone,niv.tidal,Annee,sep = "_", remove = F) %>%
# unite("ZNTPA",Zone,niv.tidal, periode, Annee,sep = "_", remove = F)
# CSLN_cont.mat<-CSLN_cont.mat %>% filter(grepl(Code,"North Median Mudflat"))


bp.adon<-adonis2(BP_con ~ BP_con_names$periode*BP_con_names$Zone,  #*CSLN_cont.mat_names$Annee, 
                 method = "bray", permutations=999)
BP_padZ<-pairwise.adonis(BP_con,BP_con_names$Zone)
BP_padP<-pairwise.adonis(BP_con,BP_con_names$periode)
# CSLN_padA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$Annee)
#CSLN_padNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$niv.tidal)
BP_padZP<-pairwise.adonis(BP_con,BP_con_names$ZP)
BP_tabAD<-as.data.frame(rbind(MSR_padZ,MSR_padP,MSR_padZP)) %>% filter(sig!="")

#ANCOVA BTA----
regr<- ggplot(CSLN.BTA) + geom_point(aes(x=MSRtot, y=solan, colour = ScientificName_accepted));regr #, na.rm = T
Table.MSR.sp<- CSLN.BTA %>% filter(ScientificName_accepted=="Peringia ulvae")
  cor.test(Table.MSR.sp$MSRtot,Table.MSR.sp$solan)


  hgv<-lm(CSLN.BTA$solan~CSLN.BTA$MSRtot*CSLN.BTA$ScientificName_accepted)
graph<- ggplot(CSLN.BTA, aes(y=solan, x=MSRtot, colour=ScientificName_accepted)) + geom_point()+
  geom_smooth(method="lm") + scale_color_brewer(palette="Set1");graph
graph1<- graph + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = ScientificName_accepted)) +
  labs( y="Potentiel de bioturbation", x="MSR (mW/m²)", color="Espece");graph1 #R² de 0,90

# # Calculer le modèle, la covariable passe en premier
# model<-lm(CSLN.BTA$solan ~ CSLN.BTA$MSRtot + CSLN.BTA$ScientificName_accepted)
# # Inspecter les paramètres de diagnostic du modèle
# model.metrics <- augment(model)  # Supprimer les détails
# head(model.metrics, 3)
# # Évaluer la normalité des résidus à l'aide du test de Shapiro-Wilk
# shapiro_test(model.metrics$.resid)
# model.metrics %>% levene_test(.resid ~ CSLN.BTA$ScientificName_accepted)
# 
# kruskal.test(MSRtot ~ ScientificName_accepted, data = CSLN.BTA)
# kruskal.test(solan ~ ScientificName_accepted, data = CSLN.BTA)
# 
# 
# pairwise.wilcox.test(CSLN.BTA$solan~CSLN.BTA$MSRtot, CSLN.BTA$ScientificName_accepted)
# 
# sm.ancova(CSLN.BTA$MSRtot, CSLN.BTA$solan, CSLN.BTA$ScientificName_accepted)
# CSLN.BTA<-as.vector(CSLN.BTA$MSRtot)
# CSLN.BTA.NUM <- CSLN.BTA %>% filter (Zone=="North Downstream Mudflat")
# regr2<- ggplot(CSLN.BTA.NUM) + geom_point(aes(x=MSRtot, y=solan)); regr2
# regrrr<- ggplot() + geom_point(aes(x=Tb.MSR3$MSR.moy, y =Tb.pBiot2$pbiot.moy));regrrr
# 
# 
# regrelin<- CSLN %>% dplyr::select(Annee,periode, Zone, ScientificName_accepted, Density_indm2,specix, MSRtot, solan) %>% filter(Zone!="Channel",specix == "Les especes +++")
# reg_lin <- lm(MSRtot~solan, data=regrelin)
# summary(regrelin)

Fauna_Brey.spe <- read.csv("C:/Users/maxcb/OneDrive/Bureau/StageM1/Scripts/MSR/Fauna_Brey spe.csv", sep=";")


CSLN<-CSLN %>% left_join(Fauna_Brey.spe)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!PEut etre enlever les NA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



######Difference entre Zones (Mudflat)/Tidal/periode/annees#############----------

    #Indice par Année/Periode ####AJOUTER niv.tidal a chaque select,group, filter

CSLN_Gip<-CSLN %>% dplyr::select( Annee,periode,Zone, Station_originelle, ScientificName_accepted,
                          SPCourt, Density_indm2,Biomass_gAFDWm2,specix) %>% 
 filter(grepl("Mudflat",Zone)) ###Changer pour autre zone
sttaion<- CSLN_Gip %>% distinct(SPCourt, ScientificName_accepted)

CSLN_Gip$periode<-as.factor(CSLN_Gip$periode)



CSLN_Gip_sp<- CSLN_Gip %>% select(Annee,periode,Zone,ScientificName_accepted, Density_indm2, Biomass_gAFDWm2, specix) %>% # filter(specix=="Les 6 especes") %>%
  complete(nesting(Zone,Annee,periode),nesting(ScientificName_accepted),
           fill = list(Biomass_gAFDWm2 = 0,Density_indm2=0)) %>%
  group_by(Zone,Annee, periode,ScientificName_accepted) %>%
  summarise(Biomasse.moy = mean(Biomass_gAFDWm2, na.rm = TRUE),
            Biomasse.sd = sd(Biomass_gAFDWm2, na.rm = TRUE),
            Densite.moy= mean(Density_indm2, na.rm = T),
            Densite.sd.moy = sd(Density_indm2, na.rm=T),
            nb.stat=n())

    #Occurence

OCC.mudflat <- CSLN_Gip %>% dplyr::select(Annee,periode,Zone, Station_originelle, ScientificName_accepted,
                                   SPCourt, Density_indm2,Biomass_gAFDWm2,specix) %>% count(specix, SPCourt,ScientificName_accepted) %>%
  arrange(desc(n)) %>% # summarise(nn=sum(n))
  slice(1:20)

OCC.mudflat$Pourcentage<- (OCC.mudflat$n /9384)*100
OCC.mudflat$Pourcentage<- round(OCC.mudflat$Pourcentage, digits = 1)

OCC.mudflat.grapho <- OCC.mudflat %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(aes(y=0.5, label = Pourcentage)) + #
  labs(title="(A)", x="espèce",y="Occurrence (%)")+ #
  theme(axis.text.y = element_text(size=1,face="bold"))+
  coord_flip()+theme_bw();OCC.mudflat.grapho

    #Densité

 OCC.mudflatDensity <- CSLN_Gip %>% select(Annee,periode,Zone, ScientificName_accepted,
                                     SPCourt, Density_indm2,Biomass_gAFDWm2,specix) %>% group_by(ScientificName_accepted) %>%
 summarise(Densite_moyenne = mean(Density_indm2)) %>%  arrange(desc(Densite_moyenne))%>% slice(1:20) 
 
 OCC.mudflatDensity$Densite_moyenne<- round( OCC.mudflatDensity$Densite_moyenne, digits = 1)
 
OCC.mudflat.graphD <- OCC.mudflatDensity %>% #pour avoir la table lancer avant pipe
  ggplot2::ggplot(aes(x=reorder(ScientificName_accepted,Densite_moyenne),y=Densite_moyenne))+geom_col(fill  ="grey",colour = "grey") +
  geom_text(aes(y=50, label = Densite_moyenne)) + # Xlab("Densité (ind/m²)") #
    labs(title="(B)", x="espèce",y="Densité moyenne (ind/m²)")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.mudflat.graphD

#Biomasse
OCC.mudflatbiom <- CSLN_Gip %>% select(Annee,periode,Zone, ScientificName_accepted,
                                          SPCourt, Density_indm2,Biomass_gAFDWm2,specix) %>% group_by( ScientificName_accepted) %>%
  summarise(biom_moyenne = mean(Biomass_gAFDWm2)) %>%  arrange(desc(biom_moyenne))%>% slice(1:20) 
OCC.mudflatbiom$biom_moyenne<- round( OCC.mudflatbiom$biom_moyenne, digits = 1)

OCC.mudflat.graphB <- OCC.mudflatbiom %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,biom_moyenne),y=biom_moyenne))+geom_col(fill  ="grey",colour = "grey") +
  geom_text(aes(y=0.25, label = biom_moyenne)) + #
  labs(title="(C)", x="espèce",y="Biomasse moyenne (g.matière sèche sans cendre/m²)")+ #
  theme(axis.text = element_text(size=1,face="bold"))+
  coord_flip()+theme_bw();OCC.mudflat.graphB
library(ggpubr)
ggarrange(OCC.mudflat.grapho, OCC.mudflat.graphD, OCC.mudflat.graphB, ncol = 1, nrow = 3)

           
write.xlsx(occ.mudflat.tb, wdgraph, asTable = FALSE, overwrite = TRUE)

#graphique Biomasse/Densite PAR ESPECES -----

# BIOMASSE
# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=periode, y=Biomasse.moy, fill= periode))+ ylim=200 +
#   facet_grid( .~Zone,scales = "free_y") +  theme(axis.text.x = element_text(angle = -90,size=4)) ;graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Biomasse par zone et annee, Specix.png",
#                                     sep = ""),width=1000,height = 600)
# # 
# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=periode, y=Biomasse.moy, fill= periode)) +
#   facet_grid(SPCourt~Zone,scales = "free_y") +  theme(axis.text.x = element_text(angle = -90,size=4));graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Biomasse par zone et periode, specix.png",
#                                     sep = ""),width=1000, height = 600)
#     #niv.tidal
# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=Annee, y=Biomasse.moy, fill= Annee)) +
#   facet_grid(SPCourt~Zone+niv.tidal,scales = "free_y") ;graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Biomasse par zone,niveau tidal et annee, specix.png",
#                                     sep = ""),width=1000, height = 600)

# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=periode, y=Biomasse.moy, fill= periode)) +
#   facet_grid(SPCourt~Zone+niv.tidal,scales = "free_y")+ 
#   theme(axis.text.x = element_text(angle = -90,size=4));graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Biomasse par zone,niveau tidal et periode, specix.png",
#                                     sep = ""),width=1000, height = 600)

#DENSITE

# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=Annee, y=Densite.moy, fill= Annee)) +
#   facet_grid(SPCourt~Zone,scales = "free_y") +  theme(axis.text.x = element_text(angle = -90,size=4)) ;graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Densite par zone et annee, specix.png",
#                                     sep = ""),width=1000,height = 600)

graph<-ggplot(CSLN_Gip_sp)+ geom_col(aes(x=periode, y=Densite.moy, fill= periode)) + 
facet_grid(.~Zone,scales = "free_y")  +  theme(axis.text.x = element_text(angle = -90,size=4)) ;graph  #Attention à l'echelle
#dev.print(device = png, file= paste(wdgraph,"CSLN_Densite moy  Mudflat et periode.png",
                                    sep = ""),width=1000, height = 600)

graph<-ggplot(CSLN_Gip)+ geom_col(aes(x=ScientificName_accepted, y= Density_indm2, fill= periode)) +
  facet_grid( .~Zone,scales = "free_y") +  theme(axis.text.x = element_text(angle = -90,size=6))  ;graph 
#niv.tidal
# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=Annee, y=Densite.moy, fill= Annee)) +
#   facet_grid(SPCourt~Zone+niv.tidal,scales = "free_y") ;graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Densite par zone,niveau tidal et annee,Specix.png",
#                                     sep = ""),width=1000, height = 600)

# graph<-ggplot(CSLN_Gip_sp)+ geom_boxplot(aes(x=periode, y=Densite.moy, fill= periode)) +
#   facet_grid(SPCourt~Zone+niv.tidal,scales = "free_y")+ 
#   theme(axis.text.x = element_text(angle = -90,size=4));graph  #Attention à l'echelle
# dev.print(device = png, file= paste(wdgraph,"CSLN_Densite par zone,niveau tidal et periode, Specix.png",
#                                     sep = ""),width=1000, height = 600)


#Indices Ecologiques (feuille d'entrée)----
    
      #Par Zone
Indice.ecoZ<- CSLN_Gip %>% group_by( Annee,periode, Zone) %>%  #niv.tidal
  summarise(Richesse.spe = n_distinct(ScientificName_accepted),
            n.indiv= n(),
            Biomasse.sum.ZNT = sum(Biomass_gAFDWm2),
            Densite_moyenne=mean(Density_indm2))


      #niv.tidal
# Indice.ecoZNT<-CSLN_Gip %>% group_by(Annee,periode, Zone, niv.tidal) %>%  #niv.tidal
#   filter(Zone !="Channel",Zone!="Cote Fleurie", Zone!="Ilot Oiseaux",
#          Zone!="Octeville",Zone!="OffShore") %>%
#   summarise(Richesse.spe = n_distinct(SPCourt),
#             n.indiv= n(),
#             Biomasse.sum.ZNT = sum(Biomass_gAFDWm2))

    #PAR ANNEE et Zone #
graph<-ggplot(Indice.ecoZ) + geom_col(aes (x=Annee, y= Richesse.spe,fill=Annee))+ 
  facet_grid(.~Zone)+ #~niv.tidal , margins = T
  theme(axis.text.x = element_text(angle = -90,size=4)) labs(title="Densités des espèces sur les zones vasières");graph #,scales="free_y"
# dev.print(device = png, file= paste(wdgraph,"CSLN_RS par Mediane Mudflat et annee.png",
#                                    sep = ""),width=1000, height = 600)
Indice.ecoZ$Annee<- as.factor(Indice.ecoZ$Annee)

graph<-ggplot(Indice.ecoZ) + geom_col(aes(x=Annee, y= Biomasse.sum.ZNT,fill=Annee))+
  facet_grid(.~Zone) + #~niv.tidal
  theme(axis.text.x = element_text(angle = -90,size=4));graph #, scales="free"
# dev.print(device = png, file= paste(wdgraph,"CSLN_biomasse tot Mediane Mudflat et annee.png",
#                                     sep = ""),width=1000, height = 600)

      #PAR PERIODE et zone
graph<-ggplot(Indice.ecoZ) + geom_boxplot(aes (x=periode, y= Richesse.spe,fill=periode))+
  facet_grid(.~Zone);graph #,scales="free_y"
# dev.print(device = png, file= paste(wdgraph,"CSLN_RS par Mediane Mudflat et periode.png",
#                                     sep = ""),width=1000, height = 600)


 graph<-ggplot(Indice.ecoZ) + geom_boxplot(aes(x=Annee, y= Biomasse.sum.ZNT,fill=Annee))+
  facet_grid(.~Zone);graph #, scales="free"
# dev.print(device = png, file= paste(wdgraph,"CSLN_biomasse tot par Mediane Mudflat et periode.png",
#                                     sep = ""),width=1000, height = 600)

      #PAR annee et niv.tidal

# graph<-ggplot(Indice.ecoZNT) + geom_col(aes (x=Annee, y= Richesse.spe,fill=Annee))+ #Indice.ecoZNT
#   facet_grid(niv.tidal~Zone,margins = T) 
#   theme(axis.text.x = element_text(angle = -90,size=4));graph #,scales="free_y"
# dev.print(device = png, file= paste(wdgraph,"CSLN_RS par zone, NT et annee.png",
#                                     sep = ""),width=1000, height = 600)
# 
# 
# graph<-ggplot(Indice.ecoZNT) + geom_col(aes(x=Annee, y= Biomasse.sum.ZNT,fill=Annee))+
#   facet_grid(niv.tidal~Zone,margins = T) + 
#   theme(axis.text.x = element_text(angle = -90,size=4));graph #, scales="free"
# dev.print(device = png, file= paste(wdgraph,"CSLN_biomasse tot par zone, NT et annee.png",
#                                     sep = ""),width=1000, height = 600)
      #PAR PERIODE
# graph<-ggplot(Indice.ecoZNT) + geom_col(aes (x=periode, y= Richesse.spe,fill=periode))+
#   facet_grid(Zone~niv.tidal,margins = T);graph #,scales="free_y"
# dev.print(device = png, file= paste(wdgraph,"CSLN_RS par zone,niveau tidal et periode.png",
#                                     sep = ""),width=1000, height = 600)
# 
# graph<-ggplot(Indice.ecoZNT) + geom_col(aes(x=periode, y= Biomasse.sum.ZNT,fill=periode))+
#   facet_grid(Zone~niv.tidal,margins = T);graph #, scales="free"
# dev.print(device = png, file= paste(wdgraph,"CSLN_biomasse tot par zone,niveau tidal et periode.png",
#                                     sep = ""),width=1000, height = 600)

#Table de contingence----
CSLN_cont<- CSLN_Gip_sp %>% dplyr::select(Zone,Annee, periode,ScientificName_accepted ,Densite.moy) %>% #niv.tidal
  group_by(Zone,Annee, periode, ScientificName_accepted) %>%
  pivot_wider(names_from = ScientificName_accepted, values_from = 
                Densite.moy, values_fill = 0)

CSLN_cont.mat<-as.matrix(CSLN_cont[,-c(1:3)])
row.names(CSLN_cont.mat)<-paste(CSLN_cont$periode, CSLN_cont$Annee, CSLN_cont$Zone, sep = "_")


CSLN_cont.mat<-subset(CSLN_cont.mat,rowSums(CSLN_cont.mat)!=0) #
CSLN_cont.mat_names<-data.frame(Names=row.names(CSLN_cont.mat)) %>%
  separate(Names,c("periode", "Annee","Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F) %>%
  unite("ZA",Zone,Annee,sep = "_", remove = F) %>%
  unite("ZPA",Zone,periode,Annee,sep = "_", remove = F) #  %>%
# unite("ZNT",Zone,niv.tidal,sep = "_", remove = F) %>%
# unite("ZNTP",Zone,niv.tidal, periode,sep = "_", remove = F) %>%
# unite("ZNTA",Zone,niv.tidal,Annee,sep = "_", remove = F) %>%
# unite("ZNTPA",Zone,niv.tidal, periode, Annee,sep = "_", remove = F)
# CSLN_cont.mat<-CSLN_cont.mat %>% filter(grepl(Code,"North Median Mudflat"))

#INDICE ECOLOGIQUE DE SHANNON/PIELOU----

Indice.ecoZm<- Indice.ecoZ %>% group_by(Zone) %>% summarise(moy.Rs = mean(Richesse.spe))
class(Indice.ecoZ$Richesse.spe)
Indice.ecoZ$Richesse.spe <- as.numeric(Indice.ecoZ$Richesse.spe)
Indice.ecoZ$Shannon<- diversity(CSLN_cont.mat, index = "shannon")

Indice.ecoZ$Pielou<-Indice.ecoZ$Shannon/log2(Indice.ecoZ$Richesse.spe)
Indice.ecoZ$Diversite.simpson<- diversity(CSLN_cont.mat, index = "simpson")
Indice.ecoZ$E_Simpson<-(1 - Indice.ecoZ$Diversite.simpson - min
                        (1 - Indice.ecoZ$Diversite.simpson )) /
  ( max(1 - Indice.ecoZ$Diversite.simpson) - min( 1 - Indice.ecoZ$Diversite.simpson)) #Équitabililité de Simpson noté E 
Indice.ecoZ$D_Hill<-diversity(CSLN_cont.mat,index = "invsimpson")/exp(Indice.ecoZ$Shannon)

#write.xlsx(Indice.ecoZ, wdgraph, asTable = FALSE, overwrite = TRUE)

      #graphique des indices pour le niveau tidal #Indice.ecoZNT + facet_grid(Zone~niv.tial)
#ligne184 <- RS
devtools::install_github('thomasp85/gganimate')
library(gganimate)
means <- aggregate(Densite_moyenne ~  periode, Indice.ecoZ, mean)
graphDens<-ggplot(Indice.ecoZ) + geom_boxplot(aes (x=periode,y= Densite_moyenne,fill=periode))+ coord_cartesian( ylim = c(0, 1035)) + 
  facet_grid(.~Zone)+ xlab(" ") + ylab("Densité moyenne (ind/m²)") + #labs(title="Richesses spécifiques par zone et période") + #~niv.tidal , margins = T 
  theme(axis.text.x = element_text(angle = -90,size= 0))+ transition_time(periode) +
  labs(title = "Year: {frame_time}");graphDens

graphrs<-ggplot(Indice.ecoZ) + geom_boxplot(aes (y= Richesse.spe,fill=periode))+ 
  facet_grid(.~Zone)+ ylab("Richesse spécifique")+ #labs(title="Richesses spécifiques par zone et période") + #~niv.tidal , margins = T 
  theme(axis.text.x = element_text(angle = -90,size= 0));graphrs #,scales="free_y"
# dev.print(device = png, file= paste(wdgraph,"CSLN_RS par Mediane Mudflat et annee.png",
#                                    sep = ""),width=1000, height = 600)

GraphS<-ggplot(Indice.ecoZ)+geom_boxplot(aes(y=Shannon,fill=periode))+ #periode
  facet_grid(.~Zone) + ylab("Indice de Shannon") + #labs(title="Shannon par zone et période") +
  theme(axis.text.x=element_text(size = 0)); GraphS
# dev.print(device = png, file= paste(wdgraph,"CSLN_Shannon  par période et zone.png",
                                    # sep = ""),width=1000, height = 600)

GraphP<-ggplot(Indice.ecoZ)+geom_boxplot(aes(y=Pielou,fill=periode))+ #periode
  facet_grid(.~Zone)+ ylab("Indice d'équitabilité Pielou") + #labs(title="Pielou par zone et période") +
  theme(axis.text.x=element_text(size = 0)); GraphP
# dev.print(device = png, file= paste(wdgraph,"CSLN_Pielou par zone et période.png",
                                    # sep = ""),width=1000, height = 600)

ggarrange(graphDens,graphrs, GraphS, GraphP,common.legend = TRUE, legend = "bottom")


# Graph<-ggplot(Indice.ecoZ)+geom_col(aes(x=Annee,y=Diversite.simpson,fill=Annee))+ #periode
#   facet_grid(.~Zone) + theme(axis.text.x=element_text(angle=-90, vjust=0.4)); Graph
# dev.print(device = png, file= paste(wdgraph,"CSLN_div_Simpson par zone et période.png",
#                                     sep = ""),width=1000, height = 600)
# 
# Graph<-ggplot(Indice.ecoZ)+geom_boxplot(aes(x=periode,y=E_Simpson,fill=periode))+ #periode
#   facet_grid(.~Zone) + theme(axis.text.x=element_text(angle=-90, vjust=0.4)); Graph
# dev.print(device = png, file= paste(wdgraph,"CSLN_E_Simpson par zone et période.png",
#                                     sep = ""),width=1000, height = 600)




#sp rare----
CSLN.abd <- abund(CSLN_cont.mat, f=0.4) #Application de la méthode TPA 
summary(CSLN.abd)
plot(CSLN.abd, dpos=c(40,100),xlab="Genres",ylab= "Abondances relatives", 
     main="Méthode TPA f=0.4") #Représentation graphique
CSLN.abd$n <- 30 # identify(CSLN.abd) #Un petit plateau est visible, cliquer sur la fin du plateau Number of variables extracted: 129 on a total of 616 
CSLN.TPA <- pastecs::extract(CSLN.abd, CSLN.abd$n) 
colnames(extract(CSLN.abd, CSLN.abd$n))

row.names(CSLN.TPA)<-paste(CSLN_cont$periode, CSLN_cont$Annee, CSLN_cont$Zone, sep = "_")


CSLN.TPA<-subset(CSLN.TPA,rowSums(CSLN.TPA)!=0)
CSLN.TPA_names<-data.frame(Names=row.names(CSLN.TPA)) %>%
  separate(Names,c("periode","Annee", "Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F) %>%
  unite("ZA",Zone,Annee,sep = "_", remove = F) %>%
  unite("ZPA",Zone,periode,Annee,sep = "_", remove = F)

# CSLN.TPA<-CSLN.TPA %>% filter(grepl(Code,"North Median Mudflat"))

#PERMANOVA : ADONIS2 et PAIRWISE----
CSLN_adonis<-adonis2(CSLN.TPA ~ CSLN.TPA_names$periode*CSLN.TPA_names$Zone,  #*CSLN_cont.mat_names$Annee, 
                     method = "bray", permutations=999) #CSLN_cont.mat_names$niv.tidal


# 
# 
# perm<-PerMANOVA.Simple( CSLN.TPA_names$periode, CSLN.TPA,nperm = 999, seed = NULL, C = NULL)
# 
# adonis(CSLN.TPA_names$periode,distance.CAH, permutations = 999, method = "bray", strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", parallel = getOption("mc.cores"))
# 
# CSLN.TPA_names$periode<-as.factor(CSLN.TPA_names$periode)
#!!!!!!!!!!!! #CSLN_pairwise_Z<-pairwise(CSLN_cont.mat ,CSLN_cont.mat_names$Zone) ###MARCHE PAS probleme sur ANNEE

CSLN_padZ<-pairwise.adonis(CSLN.TPA,CSLN.TPA_names$Zone)
CSLN_padP<-pairwise.adonis(CSLN.TPA,CSLN.TPA_names$periode)
# CSLN_padA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$Annee)
#CSLN_padNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$niv.tidal)
CSLN_padZP<-pairwise.adonis(CSLN.TPA,CSLN.TPA_names$ZP)
# CSLN_padZA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZA)
# CSLN_padZPA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZPA)
#CSLN_padZNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNT)
# CSLN_padZNTP<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTP)
# CSLN_padZNTA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTA)
# CSLN_padZNTPA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTPA)
CSLN_tabAD<-as.data.frame(rbind(CSLN_padZ,CSLN_padP,CSLN_padZP)) %>% filter(sig!="")

#write.xlsx(CSLN_adonis, wdgraph, asTable = F, overwrite = TRUE)


#Analyse des communautes (CAH)----

  ## SImprof

# Table de contingence  pour le simprof
CSLN_Gip.per<-CSLN %>% dplyr::select(periode,Zone, ScientificName_accepted, Density_indm2,Biomass_gAFDWm2,specix) %>% #, MSR_mW, MSRtot
  filter(grepl("Mudflat",Zone)) ###Changer pour autre zone

CSLN_Gip$periode<-as.factor(CSLN_Gip$periode)



CSLN_Gip_sp.per<- CSLN_Gip.per %>% dplyr::select(periode,Zone,ScientificName_accepted, Density_indm2, Biomass_gAFDWm2, specix) %>% # , MSR_mW, MSRtot
  complete(nesting(Zone,periode),nesting(ScientificName_accepted),
           fill = list(Biomass_gAFDWm2 = 0,Density_indm2=0)) %>%
  group_by(Zone, periode,ScientificName_accepted) %>%
  summarise(Biomasse.moy = mean(Biomass_gAFDWm2, na.rm = TRUE),
            Biomasse.sd = sd(Biomass_gAFDWm2, na.rm = TRUE),
            Densite.moy= mean(Density_indm2, na.rm = T),
            Densite.sd.moy = sd(Density_indm2, na.rm=T)) #, MSRtot_m=median(MSRtot,na.rm=TRUE),MSRtot_sd=sd(MSRtot,na.rm=TRUE))

CSLN_cont.per<- CSLN_Gip_sp.per %>% dplyr::select(Zone,periode,ScientificName_accepted,Densite.moy) %>% #niv.tidal
  group_by(Zone,periode, ScientificName_accepted) %>%
  pivot_wider(names_from = ScientificName_accepted, values_from = 
                Densite.moy, values_fill = 0)

CSLN_cont.mat.per<-as.matrix(CSLN_cont.per[,-c(1:2)])
row.names(CSLN_cont.mat.per)<-paste(CSLN_cont.per$periode,CSLN_cont.per$Zone, sep = "_")


CSLN_cont.mat.per<-subset(CSLN_cont.mat.per,rowSums(CSLN_cont.mat.per)!=0) #
CSLN_cont.mat_names.per<-data.frame(Names=row.names(CSLN_cont.mat.per)) %>%
  separate(Names,c("periode","Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F)


#sp rare par periode
CSLN.abd.per <- abund(CSLN_cont.mat.per, f=0.4) #Application de la méthode TPA 
summary(CSLN.abd.per)
plot(CSLN.abd.per, dpos=c(40,100),xlab="Genres",ylab= "Abondances relatives", 
     main="Méthode TPA f=0.4") #Représentation graphique
CSLN.abd.per$n <- 60 # identify(CSLN.abd) #Un petit plateau est visible, cliquer sur la fin du plateau Number of variables extracted: 129 on a total of 616 
colnames(extract(CSLN.abd.per, CSLN.abd.per$n))

CSLN.TPA.per <- pastecs::extract(CSLN.abd.per, CSLN.abd.per$n) 


# medianglob<-apply(CSLN_cont.mat.per,2,median) #Vecteur avec les médianes par genre
# boxplot (medianglob, xlab="Mediane globale", ylab="Abondances (ind/L)")
# subset(t(CSLN_cont.mat.per),medianglob>median(medianglob))->phyto2 #Sélection des genres dont les médians sont supérieurs à la médiane des médianes
# dim(phyto2)
# phytosimp<-t(phyto2)



row.names(CSLN.TPA.per)<-paste(CSLN_cont.per$periode, CSLN_cont.per$Zone, sep = "_")


CSLN.TPA.per<-subset(CSLN.TPA.per,rowSums(CSLN.TPA.per)!=0)
CSLN.TPA_names.per<-data.frame(Names=row.names(CSLN.TPA.per)) %>%
  separate(Names,c("periode","Zone"),sep = "_", remove = T) %>%
  unite("ZP",Zone, periode,sep = "_", remove = F) #%>%
 # unite("ZA",Zone,Annee,sep = "_", remove = F) #%>%
  #unite("ZPA",Zone,periode,Annee,sep = "_", remove = F)



Indice.ecoZP<- CSLN_Gip.per %>% group_by(periode, Zone) %>%  #niv.tidal
  summarise(Richesse.spe = n_distinct(ScientificName_accepted),
            n.indiv= n())

Indice.ecoZP$Shannon<- diversity(CSLN.TPA.per, index = "shannon")

Indice.ecoZP$Equi.pielou<-Indice.ecoZP$Shannon/log2(Indice.ecoZP$Richesse.spe)
# Indice.ecoZ$Diversite.simpson<- diversity(CSLN_cont.mat, index = "simpson")
# Indice.ecoZ$E_Simpson<-(1 - Indice.ecoZ$Diversite.simpson - min
#                         (1 - Indice.ecoZ$Diversite.simpson )) /
#   ( max(1 - Indice.ecoZ$Diversite.simpson) - min( 1 - Indice.ecoZ$Diversite.simpson)) #Équitabililité de Simpson noté E 
# Indice.ecoZ$D_Hill<-diversity(CSLN_cont.mat,index = "invsimpson")/exp(Indice.ecoZ$Shannon)

#PERMANOVA : ADONIS2 et PAIRWISE
CSLN_adonis<-adonis2(CSLN_cont.mat ~ CSLN_cont.mat_names$periode*CSLN_cont.mat_names$Zone,  #*CSLN_cont.mat_names$Annee, 
                     method = "bray", permutations=999) #CSLN_cont.mat_names$niv.tidal


# 
# 
# perm<-PerMANOVA.Simple( CSLN.TPA_names$periode, CSLN.TPA,nperm = 999, seed = NULL, C = NULL)
# 
# adonis(CSLN.TPA_names$periode,distance.CAH, permutations = 999, method = "bray", strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", parallel = getOption("mc.cores"))
# 
# CSLN.TPA_names$periode<-as.factor(CSLN.TPA_names$periode)
#!!!!!!!!!!!! #CSLN_pairwise_Z<-pairwise(CSLN_cont.mat ,CSLN_cont.mat_names$Zone) ###MARCHE PAS probleme sur ANNEE

CSLN_padZ<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$Zone)
CSLN_padP<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$periode)
# CSLN_padA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$Annee)
#CSLN_padNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$niv.tidal)
CSLN_padZP<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZP)
CSLN_padZP$sig[CSLN_padZP$p.value <= 0.05 ]<-"*"
# CSLN_padZA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZA)
# CSLN_padZPA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZPA)
#CSLN_padZNT<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNT)
# CSLN_padZNTP<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTP)
# CSLN_padZNTA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTA)
# CSLN_padZNTPA<-pairwise.adonis(CSLN_cont.mat,CSLN_cont.mat_names$ZNTPA)
CSLN_tabAD<-as.data.frame(rbind(CSLN_padZ,CSLN_padP)) %>% filter(sig!="")


      #CAH station / zone
simprof.Bray <- simprof(data=CSLN_cont.mat.per, num.expected=1000,
                        num.simulated=999,method.cluster="average",
                        method.distance="braycurtis",
                        sample.orientation ="row", undef.zero = T)

SS<-simprof.plot(simprof.Bray, leaflab="perpendicular", siglinetype = 1, ) 
# groupesdend<-as.data.frame(cbind(code=unlist(simprof.Bray[["significantclusters"]]),Groupe_Simprof=simprof.Bray[["significantclusters"]]))
groupesdend<- simprof.Bray$significantclusters %>% set_names(seq_along(.)) %>%
  enframe %>% unnest(cols=c(value)) %>% rename(Groupe_simprof=name,code=value) %>% 
  mutate(Groupe_simprof = parse_integer(Groupe_simprof))
grp.CAH<-groupesdend$Groupe_simprof
names(grp.CAH)<-groupesdend$code

CSLN.ckd<- CSLN %>% dplyr::select(Zone, periode,ScientificName_accepted,Density_indm2, specix) %>%
  group_by(Zone,periode, ScientificName_accepted,specix) %>% summarise(d.moy= mean(Density_indm2)) %>% filter(Zone!="Channel", specix=="Les especes +++")

graphndm<-ggplot(Indice.ecoZ) + geom_boxplot(aes (y= Densite_moyenne,fill=periode))+ coord_cartesian( ylim = c(0, 1035)) + 
  facet_grid(.~Zone)+ #labs(title="Richesses spécifiques par zone et période") + #~niv.tidal , margins = T 
  theme(axis.text.x = element_text(angle = -90,size= 0));graphDens

# dev.print(device = png, file = paste(wdgraph,"CSLN_Simprof.png",sep=""), width = 1000, height=600)
# dev.off()
#  library(dendextend)
# ggd1 <- simprof.Bray$hclust %>% as.dendrogram %>% ladderize 
# ggd1<- color_branches(ggd1, grp.CAH, groupLabels = TRUE)
# ggd1 <- ggd1 %>% #set("branches_k_color", k=simprof.Bray$numgroups,value=colorRampPalette(brewer.pal(8, "Dark2"))(grp.CAH)) %>%
#   set("branches_lwd", .8) %>%
#   set("labels_col", k=simprof.Bray$numgroups, value=colorRampPalette(brewer.pal(8, "Dark2"))(simprof.Bray$numgroups)) %>%
#   set("labels_cex", .6) %>%
#     set("hang_leaves", -1) %>% 
#    set("leaves_col")
# ggd1<- as.ggdend(ggd1)
# 
# ggd1 <- ggplot(ggd1, horiz = TRUE, labels = TRUE); ggd1 #+ geom_text(size = 2)




dev.print(device = png, file = paste(wdgraph,"CSLN_Simprof3.png",sep=""), width = 1000, height=600)
dev.off()

 library(cluster); library(vegan); library(RColorBrewer)

# finalclust.CSLN<- reorder(simprof.Bray$significantclusters,CSLN.TPA) #Re-ordination des stations selon les contraintes de l'analyse de groupement
# dend<-as.dendrogram(simprof.Bray$hclust)
# or<-vegemite(CSLN.TPA, simprof.Bray$hclust, scale="Hill")
# heatmap(t(log1p(phytoS)[rev(or$species)]), Rowv=NA, Colv=dend, col=c("white",brewer.pal (5, "Greens")),
#         margin=c(4,4), ylab="genres (moyennes ponderées par sites)", xlab="Stations")


# ggg<- ggdendrogram(simprof.Bray$hclust, theme_dendro = FALSE); ggg

# simprof.plot(CSLN_simprof, leafcolors=NA, plot=TRUE, fill=TRUE,
#              leaflab="perpendicular", siglinetype=1)
# dev.print(device = png, file = paste(wdgraph,"CSLN_Simprof.png",sep=""), width = 1000, height=600)
# dev.off()

# plot(simprof.Bray$hclust)
# dev.print(device = png, file= paste(wdgraph,"CSLN_Med_mudflat_Simprof.png",
# sep = ""),width=1000,  height = 600)
##########   6 GROUPES

#Indice.ecoZ$groupe.CAH[<-simprof.Bray[["significantclusters"]]

#dev.off() # Pour nettoyer les graphiques
#OU
# distance.CAH<-vegdist(CSLN.TPA.per, method="bray")
# CAH.bray<-hclust(distance.CAH, method = "average")
# plot(CAH.bray,hang=-1, main="Liens complets", xlab="Zone & niveau tidal", ylab="Dissimilarités de Bray-Curtis")
# Choix.groupe<-simprof.Bray$hclust
# rect.hclust(CAH.bray,k=Choix.groupe)
# grp.CAH<-cutree(CAH.bray,k=Choix.groupe)
#
# Indice.ecoZ$grp.CAH<- grp.CAH # %>% relocate(Choix.groupe, .before = Richesse.spe) !!!!!!!!!!#NE SUIT PAS LES DONNEES SIGNIF DE SIMPROF

#ou
# ggdendrogram(simprof.Bray$hclust,colour=simprof.Bray$hclust)
# rect.hclust(simprof.Bray$hclust)
#En utilisant CSLN.TPA

# Significativite des groupes
# adonis_CAH<-adonis2(CSLN_cont.mat.per~grp.CAH,method="bray"); adonis_CAH
# pad_CAH<-pairwise.adonis(CSLN_cont.mat.per,grp.CAH); pad_CAH

# Contribution specifique (par groupe CAH)----

OCC.com1 <-CSLN %>% 
  filter(Zone=="North Downstream Mudflat", specix =="Les especes +++") %>% #, periode=="2006-2010"| periode=="2000-2005"
  group_by(ScientificName_accepted) %>% summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc)) 
OCC.com1$Pourcentage<- (OCC.com1$Densitemoc /sum(OCC.com1$Densitemoc))*100
OCC.com1$Pourcentage<- round(OCC.com1$Pourcentage, digits = 1)
OCC.com12 <- OCC.com1 %>% dplyr::slice(1:15)
# SP<- OCC.com12$ScientificName_accepted
OCC.com1.graph <- OCC.com12 %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(size = 2.8,aes(y=0.5, label = Pourcentage)) + #
  labs(title="Groupe 1", x="especes",y="Occurrence")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.com1.graph
# OCC.com11 <-CSLN %>% 
#   filter(Zone=="North Downstream Mudflat", ScientificName_accepted=SP) %>% #, periode=="2006-2010"| periode=="2000-2005"
#   group_by(ScientificName_accepted) %>% summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc)) 


OCC.com2 <-CSLN %>%  filter(Zone=="North Downstream Mudflat", periode=="2006-2010", specix=="Les especes +++") %>% group_by(ScientificName_accepted) %>%
  summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc))
OCC.com2$Pourcentage<- (OCC.com2$Densitemoc /sum(OCC.com2$Densitemoc))*100
OCC.com2$Pourcentage<- round(OCC.com2$Pourcentage, digits = 1)
OCC.com2 <- OCC.com2 %>% dplyr::slice(1:15)
OCC.com2.graph <- OCC.com2 %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(size = 2.8,aes(y=0.5, label = Pourcentage)) + #
  labs(title="Groupe 2", x="especes",y="Occurrence")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.com2.graph

OCC.com3 <-CSLN %>%  filter(Zone=="North Upstream Mudflat",specix=="Les especes +++") %>% group_by(ScientificName_accepted) %>%
  summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc))
OCC.com3$Pourcentage<- (OCC.com3$Densitemoc /sum(OCC.com3$Densitemoc))*100
OCC.com3$Pourcentage<- round(OCC.com3$Pourcentage, digits = 1)
OCC.com3 <- OCC.com3 %>% dplyr::slice(1:15)
OCC.com3.graph <- OCC.com3 %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(size = 2.8,aes(y=1, label = Pourcentage)) + #
  labs(title="Groupe 3", x="especes",y="Occurrence")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.com3.graph

OCC.com4 <-CSLN %>%  filter(Zone=="North Median Mudflat",specix=="Les especes +++") %>% group_by(ScientificName_accepted) %>%
  summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc))
OCC.com4$Pourcentage<- (OCC.com4$Densitemoc /sum(OCC.com4$Densitemoc))*100
OCC.com4$Pourcentage<- round(OCC.com4$Pourcentage, digits = 1)
OCC.com4 <- OCC.com4 %>% dplyr::slice(1:15)
OCC.com4.graph <- OCC.com4 %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(size = 2.8,aes(y=0.5, label = Pourcentage)) + #
  labs(title="Groupe 4", x="especes",y="Occurrence")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.com4.graph

OCC.com5 <-CSLN %>%  filter(Zone=="South Mudflat") %>% group_by(ScientificName_accepted) %>%
  summarise(Densitemoc = mean(Density_indm2)) %>% arrange(desc(Densitemoc)) 
OCC.com5$Pourcentage<- (OCC.com5$Densitemoc /sum(OCC.com5$Densitemoc))*100
OCC.com5$Pourcentage<- round(OCC.com5$Pourcentage, digits = 1)
OCC.com5 <- OCC.com5 %>% dplyr::slice(1:15)
OCC.com5.graph <- OCC.com5 %>% #pour avoir la table lancer avant pipe
  ggplot(aes(x=reorder(ScientificName_accepted,Pourcentage),y=Pourcentage))+geom_col(fill = "grey", colour = "grey") +
  geom_text(size = 2.8, aes(y=0.5, label = Pourcentage)) + #
  labs(title="Groupe 5", x="especes",y="Occurrence")+ #
  theme(axis.text = element_text(size=15,face="bold"))+
  coord_flip()+theme_bw();OCC.com5.graph

OCC.com<- ggarrange(OCC.com1.graph, OCC.com2.graph,OCC.com3.graph, OCC.com4.graph,OCC.com5.graph); OCC.com


# Indval ----
#### la considération d'une espèce indicatrice comme une espèce spécifique à un groupe donné et fidèle à toutes les stations/dates de ce groupe (IndVal
CSLN_cont.mat.per2 <- CSLN_cont.mat.per[, (!apply(CSLN_cont.mat.per==0,2,all))] 
CSLN.Indval<-indval(CSLN.TPA.per,grp.CAH) #CSLN_cont.mat.per2
GR<-CSLN.Indval$maxcls #Groupe de station où chaque genre a un indval max *
gr<-CSLN.Indval$indcls #Groupe de station avec indval correspondant aux indval max 

pv.indval<-CSLN.Indval$pval #p-value correspondantes 
tab.CSLN.Indval<-data.frame(group=GR,indval=gr,pvalue=pv.indval) #Tableau avec pour les espèces indicatrices, leur appartenance aux groupes et la p-value correspondante
tab.CSLN.Indval<-cbind(sp=row.names(tab.CSLN.Indval),tab.CSLN.Indval)[order(tab.CSLN.Indval$group,tab.CSLN.Indval$pvalue) ,] #Trier les indVal selon la p-value obtenue par les tests de permutation

tab.CSLN.Indval$specix<-"Autre especes"
tab.CSLN.Indval$specix[tab.CSLN.Indval[,1] %in% c("CERED","CORVO","HEDDI","MACBA","PERUL","SCRPL","AREMA","BATPI","BATSA","CYACA","NEPCI",
                                                  "NEPHO","PYGEL","ETELO")]<-"Les especes++"
write.xlsx(tab.CSLN.Indval, wdgraph, asTable = F, overwrite = TRUE)

#SIMPER----
#déterminer les contributions de chaque espèce à la dissimilarité entre stations

Simp<-simper(CSLN_cont.mat.per,CSLN_mds_dat$Groupe_Simprof , permutations = 0, trace = FALSE,  parallel = getOption("mc.cores"))
CC<-summary(Simp, ordered = TRUE,
        digits = max(3,getOption("digits") - 3))



#NMDS----
CSLN.mds<-metaMDS( CSLN_cont.mat, k=2, trymax = 100,autotransform=F) 

# CSLN.mds.point<-data.frame(cbind(Names=row.names(CSLN.mds$points),CSLN.mds$points))%>%
#                              separate(Names,c("Zone","period","Annee"),sep = "-", remove = F)
# rownames(CSLN.mds)<-c()
# 
# CSLN_mds_dat<-data.frame(scores(CSLN.mds)) %>%
#   mutate(Names=row.names(scores(CSLN.mds$points),CSLN.mds$points)) %>%
#   separate(Names,c("Zone","period","Annee"), sep ="-", remove = F)

yy<-plot(CSLN.mds, type = "n");yy
ordiplot(CSLN.mds,type="n")
orditorp(CSLN.mds,display="species",col="red",cex = 0.25,air=0.01,)
orditorp(CSLN.mds,display="sites",cex=0.5,air=0.01,)

#__________________________

#  NMDS par ZONE
CSLN_cont_nMDS<-CSLN.TPA # CSLN_cont
CSLN_cont_nMDS_mat<- CSLN_cont.mat

CSLN.mds <- metaMDS(CSLN.TPA.per,distance = "bray", k = 2,try = 10000, autotransform =FALSE)
# CSLN_mds_points <- data.frame(cbind(code=row.names(CSLN_mds$points),CSLN_mds$points)) %>% 
#   separate(code, c("Zone", "periode","Annee"), sep = "_",remove = FALSE)
# rownames(CSLN_mds_points) <- c()

CSLN_mds_dat <- data.frame(scores(CSLN.mds)) %>% # Using the scores function from vegan to extract the code scores
  mutate(code=row.names(scores(CSLN.mds))) %>% relocate(code) %>%
  separate(code, c("periode","Zone"), sep = "_",remove = FALSE) %>%
  unite(ZP,Zone,periode, sep = "_",remove = FALSE) %>%
  group_by(Zone,periode) %>% #,Tidal_level
  mutate(NMDS1.m = mean(NMDS1, na.rm = TRUE),
         NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_mds_dat <- CSLN_mds_dat %>% left_join(groupesdend)
CSLN_mds_dat$Groupe_simprof<-as.factor(CSLN_mds_dat$Groupe_Simprof)
rownames(CSLN_mds_dat) <- c()
CSLN_mds_mean <- CSLN_mds_dat %>% 
  unite("ZP",Zone,periode, sep = "_",remove = FALSE) %>%
  group_by(ZP) %>% 
  summarise(NMDS1.m = mean(NMDS1, na.rm = TRUE),
            NMDS2.m = mean(NMDS2, na.rm = TRUE))
CSLN_mds_species<-as.data.frame(CSLN.mds$species)
CSLN_mds_species<-cbind(sp=rownames(CSLN_mds_species),CSLN_mds_species)  # create a column of species, from the rownames of CSLN_mds_species

nbA = dim(CSLN_mds_mean)[1]
nbB = length(unique(CSLN_mds_dat$periode))
palette= colorRampPalette(brewer.pal(8, "Dark2")) # colDarj colZiss colSpec
bp<-ggplot(data=CSLN_mds_dat) +
  # Add species name
  # geom_point(data=CSLN_mds_species,aes(x=MDS1,y=MDS2),shape=8,size=1,colour="gray50")+ # add the point markers
  geom_text(data=CSLN_mds_species,aes(x=MDS1,y=MDS2,label=sp),
            size=2.5,colour="gray50",check_overlap = TRUE, fontface = "bold")  +
  # Base scatter plot 
  geom_point(aes(x=NMDS1,y=NMDS2,shape=periode,colour=Zone),size=3) +
  #geom_text(aes(x=NMDS1,y=NMDS2,colour=Zone,label=Zone),size=10.5,hjust=-0.1, fontface = "plain")+  # add the site labels
  scale_colour_manual(values=palette(nbA)) +
 # stat_ellipse(aes(x=NMDS1,y=NMDS2,group=Groupe_simprof), #linetype = Groupe_simprof,color=Groupe_simprof,
             #  type = "norm",level = 0.9,linetype = 3,show.legend = FALSE) + 
  # geom_label(aes(x=CSLN_mds_mean$NMDS1.m,y=CSLN_mds_mean$NMDS2.m,label=CSLN_mds_mean$Zone),
  #            fill = palette(nbA),size = 3,label.size =0,alpha=0.5)+
  # Graphic visuals
  annotate("text",x=min(CSLN_mds_dat$NMDS1),y=min(CSLN_mds_dat$NMDS2),hjust=-2,
           label=paste("Stress =",round(CSLN.mds$stress,4)), size = 4)  +
 # ggtitle(label="nMDS par zone Mudflat et periode")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        plot.background = element_blank());bp
ggsave(paste(wdgraph,"nMDS Macrofaune by Zone3",".png",sep=""), plot = bp, width = 12, height = 8)


#ou
# CSLN.mds<-data.frame(CSLN.mds)
# CSLN_mds_dat<-scores(CSLN.mds) %>% #score permet d'extraire la colonne Names
#   mutate(Names=row.names(scores(CSLN.mds))) %>%
#   ggplot(aes(x = NMDS1, y = NMDS2)) +
#   geom_point(aes(size = Richness, color = Group)) +
#   stat_ellipse(geom = "polygon", aes(group = Group, color = Group, fill = Group), alpha = 0.3) +
#   annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(CSLN.mds$stress, digits = 4)), hjust = 0) +
#   theme_bw()
#ou

# Treemap a 3 niveaux: Distribution spatial ----
library(treemap)

Sumry<-CSLN %>%filter(grepl("Mudflat", Zone))  %>% group_by(Zone,SPCourt, periode) %>%
  summarise(Smry=sum(Density_indm2,na.rm=TRUE),Nb=n())

Sumry$label <- paste(Sumry$SPCourt, Sumry$Nb, sep = "\n")
png(filename=paste(wdgraph,"_TreeMap3 Zone Density",".png",sep=""),width=1000, height=600)
treemp <-treemap(Sumry,index=c("Zone","periode", "label"),vSize="Smry",vColor="Zone",type="index",
  fontface.labels=c(2,2,1),
        border.col=c("black","white","gray81"), border.lwds=c(4,3,1),bg.labels=0,
        fontsize.labels=c(14,14,10), fontcolor.labels=c("black","white","gray81"),
        align.labels=list( c("left","bottom"), c("left","top"), c("right","center")));treemp
dev.off()

install.packages('influencer')
devtools::install_github("timelyportfolio/d3treeR")
library(d3treeR)
inter <- d3tree2( treemp ,  rootname = "General" )
library(htmlwidgets)
saveWidget(inter, file=paste(wdgraph, "Treemap_anim",".html",sep=""))  
##Analyse statisitique----


#----------------------------------------------------------------
  
#GRAPHIQUE-----
 #https://ericmarcon.github.io/MesuresBioDiv2/entropie.html#d%C3%A9finition-de-lentropie simpson shannon

  
  
  ####SAVE-----
save.image(file = paste(wdwork,"BBDD_CSLN_Annee_periode_zone_nt",".RData", sep = ""))

save.image(file = paste(wdwork,"Script Mudflat",".RData", sep = ""))
