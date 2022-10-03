# ____________________________________________________________________________________________________
# CUBIC B-SPLINES QUANTILE REGRESSION ANALYSIS WITH KOENKER'S QUANTREG PACKAGE
# ____________________________________________________________________________________________________
  
rm(list=ls())
#________________________________________________________________
# PACKAGES USED ----
library(readxl) ; library(openxlsx) # Edition d'un fichier Excel
library(tidyverse)
library(GGally); library(ggpubr)
library(RColorBrewer);
library(quantreg);# library(visreg)
library(splines); library(ggforce)# bs BSplines function
library(plotly); library(plot3D);  # graphiques 3D plot 3D for mesh library(pracma) 
# GIS Packages
library(sf); library(sfheaders); # st_as_sf ; sf_to_df
library(tmap) # tmap_mode; for static and interactive maps
library(htmlwidgets) # library(leaflet) # saveWidget ; for interactive maps

#________________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
colorsS <- colorRampPalette(brewer.pal(8, "Spectral")); 
coldiscS <- function(x) {scale_colour_manual(values=colorsS(x))}
coldiscFS <- function(x) {scale_fill_manual(values=colorsS(x))}
colorsD <- colorRampPalette(brewer.pal(8, "Dark2")); coldiscD <- function(x) {scale_colour_manual(values=colorsD(x))}
collist<-c("#8b6508","#458b74","#9a32cd","#458b00","#cdcd00","#005a7c")
blank<-alpha("#dae8ed",0.1); colIC<-"#800040" ; prettyred ="#CD2626"
# ----

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "E:/" #"C:/Users/lehuen201/Nextcloud/" "E:/" 
tsk <- "A_SDM_NEO/"
wdpath <- paste(pc,"Melting Pot/BDD/",tsk,sep=""); 
wdwork <- paste(wdpath,"Matrices/",sep="")
wdgraph <- paste(wdpath,"Graphiques/",sep="")
wdgraphEx<-paste(pc,"Copie-HD/Melting Potes/",tsk,"Graphiques/",sep="")
wdres <- paste(wdpath,"Resultats/",sep="")
wdGIS <- paste(pc,"Melting Pot/SIG/",sep="");
setwd(wdpath)

#________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
prgm <- "CSLN" # 1:CSLN 2:Mabes 3:Geco 4:Beaug
etude <- paste(prgm,"_Mars",sep="")
load(paste(wdwork,etude,"_BDD.RData",sep=""))
  
choixttt <- 123 # treatment choices 1: RQ 1 to 4 factors coeff and AICm; 2: rq 2d graphs or combinaisons

# Choix des predicteurs utilises selon etude autocorrelation, taus...
# Var_choosen<-c("flow_mxd","inunt","flow_m","sal_m","tenfon_m","mudrate_m") #,"bathy"
Var_choosen<-c("flow_mxd","inunt","sal_dtd","temp_m","tenfon_m","mudrate_m")
pred_red <- predict  %>% subset(Var %in% Var_choosen) %>% arrange(match(Var,Var_choosen))
taus <- c(0.5,0.85,0.9,0.95,0.975,0.99) #c(1:19/20,0.975,0.99)
spe <- 1#:nrow(species) # 1:CERED 2:CORVO 3:HEDDI 4:LIMBA 5:PERUL 6:SCRPL
# reponse<-reponse[1:3,] # 1:Biomass_gAFDWm2 2:Density_indm2 3:MSRtot 4:Itot  
answ <- 1#nrow(reponse)
sai <- 1#:nrow(saison) # 1:Year 2:Winter 3:Summer
explo <- 1:nrow(pred_red) # Choix du domaine d'exploration des predicteurs reduits
anm <- 1:length(anMars) # balayage des annees dispo Mars

# Definition des SDM choisis
name<-c('SDM1','SDM2','SDM3','SDM4','SDM5','SDM6')
SDM_desc<-c('Cubic BSpline two factors','Cubic BSpline two factors','Cubic BSpline two factors','Cubic BSpline two factors','Cubic BSpline two factors','Cubic BSpline two factors')
SDM_tau<-c(0.99,0.975,0.975,0.975,0.975,0.975)
yt<-c('Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2')
x1t<-c('flow_mxd','flow_mxd','flow_mxd','inunt','sal_dtd','tenfon_m')
x2t<-c('inunt','mudrate_m','tenfon_m','mudrate_m','mudrate_m','mudrate_m')
sdmlist<-data.frame(name,SDM_desc,SDM_tau,yt,x1t,x2t)
sdm<-1:nrow(sdmlist)
# ----
# if exists
# load(paste(wdwork,etude,"_RQ_BDD",".RData", sep=""))

#________________________________________________________________
# 1 : ONE to TWO FACTOR Cubic BSpline RQ with AICm and Rone calculation ----
if (choixttt==1 | choixttt==12 | choixttt==13 | choixttt==123){
tmp<-NA
indrq_l <- as.data.frame(matrix(nrow=0,ncol=7)); 
smrq_l <- as.data.frame(matrix(nrow=0,ncol=8))
for (sp in spe) { # sp=1
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == species[sp,1]),]
  indrq_s <- as.data.frame(matrix(nrow=0,ncol=7)); 
  smrql_s <- as.data.frame(matrix(nrow=0,ncol=8))
  for (rep in answ){ # rep=1
    indrq <- as.data.frame(matrix(nrow=0,ncol=7)); 
    smrql <- as.data.frame(matrix(nrow=0,ncol=8))
    for (sa in sai) {# sa=1
      indrq_tmp <- as.data.frame(matrix(nrow=0,ncol=7)); 
      smrq_tmp <- as.data.frame(matrix(nrow=0,ncol=9))
      for (k in explo) { # k=1 
        
        # ONE FACTOR RQ ----
        yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
        zt = "Zone"
        x1t = sprintf("%s%s",pred_red[k,1],saison[sa,1]); x1l = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
        xt<-x1t; xl<-x1l
        xxx <-c("Beta 0",x1l)
        dfrq <-df[,c(yt,x1t,zt)]; dfrq <-as.data.frame(na.omit(dfrq))
        y <- jitter(dfrq[,1]); x1 <- jitter(dfrq[,2])
        z <- dfrq[,3]; nbcol <- length(levels(factor(z)))
        tryCatch({
          modelq<-rq(y~bs(x1,degree=3,knots=median(x1)), tau=taus)
          modelq0<-rq(y ~ 1, tau=taus)
          modx<-seq(from = min(x1), to = max(x1),length.out=length(x1))
          modelqpred<-as.data.frame(cbind(modx,predict(modelq,newdata=data.fram
                                                       e(x = modx))))
          tmp <- c(species[sp,1],yl,x1l,"Simple",
                          round(median(1 - modelq$rho/modelq0$rho),5), #Rone
                          round(median(AIC(modelq)),1))
          indrq_tmp <- rbind(indrq_tmp,data.frame(t(tmp))); tmp<-NA
          smrq <- summary(modelq, se="boot") #, se="boot" or "boot" # avec boot, p value et std err
          for (t in 1:length(taus)){
            tmp<-data.frame(Sp=species[sp,1],reponse=yl,predict=xl,mode="Simple",tau=taus[t],
                            Var=row.names(smrq[[t]][["coefficients"]]),as.data.frame(smrq[[t]][["coefficients"]])) #modelq[["coefficients"]][,t]
            smrq_tmp <- rbind(smrq_tmp,tmp); tmp<-NA
          }
          # ONE FACTOR LINEAR SUMMARY GRAPHS ----
          titreG <- sprintf("%s for %s in %s - %s",yl,species[sp,2],saison[sa,2],prgm)
          png(file=sprintf("%s%s/RQ BSpline/Summaries/%s_CbsRq_sm_%s_%s_%s.png",
                            wdgraph,species[sp,1],prgm,species[sp,1],yt,xt),width=600, height=600)
          sm<-plot(smrq, cex=.7,pch=19,lcol=colorsS(3)[3],col=c(colorsS(3)[1],colorsS(3)[2]),xlab = "tau", ylab = yt, )
          title(main=paste(titreG,"\n\n\n",sep=""))
          dev.off()
          modelq<-NULL; # smrq<-NULL;     
        },error = function(e) {print(e)})#,finally = {})
        
        #tryCatch({ ONE FACTOR LINEAR GRAPHS ----
        titreG <- sprintf("%s - in %s - %s\n%s vs %s",species[sp,2],saison[sa,2],prgm,yl,xl)
        dp <-  ggplot(dfrq, aes_string(x=x1t, y=yt, color=zt)) + geom_point() + coldiscS(nbcol) +
                    labs(title=titreG,x=paste(x1l," (",pred_red[k,3],")",sep=""),y=yl)
          for (t in 1:length(taus)){
          dp <- dp +
                    geom_line(aes(x=modelqpred$modx,y=modelqpred[,t]))+
                    theme(legend.position="bottom")+
                    annotate("text",x=max(x1,na.rm=TRUE)*0.9,y=max(y)*((t/length(taus))*0.9),
                             colour=colorsD(length(taus))[t],size=4,fontface=2,hjust=0, label = paste(taus[t]))
        }
        ggsave(sprintf("%s%s/RQ BSpline/RQL Simple/%s_CbsRq_%s_%s_%s_%s.png",
                       wdgraph,species[sp,1],prgm,species[sp,1],yt,x1t,saison[sa,2]),plot = dp, width = 8, height = 8)
        # },error = function(e) {print(e)})#,finally = {})
        # ----
               
        # Two factors loop ----
        for (k2 in explo[-c(1:k)]) { #k2=2
          yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
          x1t = sprintf("%s%s",pred_red[k,1],saison[sa,1]); x1l = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
          x2t = sprintf("%s%s",pred_red[k2,1],saison[sa,1]); x2l = sprintf("%s%s (%s)",pred_red[k2,2],saison[sa,1],pred_red[k2,3])
          
          # Two factors loop Addition----
          xt<-paste(x1t,x2t,sep='+'); xl<-paste(x1l,x2l,sep='+')
          xxx <-c("Beta 0",x1l,x2l)
          dfrq <-df[,c(yt,x1t,x2t)]; dfrq <-as.data.frame(na.omit(dfrq))
          y <- jitter(dfrq[,1]); x1 <- jitter(dfrq[,2]); x2 <- jitter(dfrq[,3]);
          tryCatch({
            modelq<-rq(y~bs(x1,degree=3,knots=median(x1))+bs(x2,degree=3,knots=median(x2)), tau=taus)
            modelq0<-rq(y ~ 1, tau=taus)
            tmp <- c(species[sp,1],yl,xl,"Addition2",
                            round(median(1 - modelq$rho/modelq0$rho),5),
                            round(median(AIC(modelq)),1))
            indrq_tmp <- rbind(indrq_tmp,data.frame(t(tmp))); tmp<-NA
            smrq <- summary(modelq, se="boot") 
            for (t in 1:length(taus)){
              tmp<-data.frame(Sp=species[sp,1],reponse=yl,predict=xl,mode="Addition2",tau=taus[t],
                              Var=xxx,as.data.frame(smrq[[t]][["coefficients"]])) #modelq[["coefficients"]][,t]
              smrq_tmp <- rbind(smrq_tmp,tmp); tmp<-NA
            }
            smrq<-NULL; modelq<-NULL        
          },error = function(e) {print(e)})#,finally = {})
        }
        # ----
      } # explo
    indrq <- rbind(indrq,indrq_tmp)
    smrql <- rbind(smrql,as.data.frame(smrq_tmp))
    } # saison
    indrq_s <- rbind(indrq_s,indrq)
    smrql_s <- rbind(smrql_s,smrql)
  } # reponse
  indrq_l <- rbind(indrq_l,indrq_s)
  smrq_l <- rbind(smrq_l,smrql_s)
} # species
colnames(indrq_l)<-c("Sp","reponse","predict","mode","Rone","AICm");
# indrq_l$AICm<-as.numeric(levels(indrq_l$AICm))[indrq_l$AICm] # AICm se retrouve en factor, conversion
# indrq_l$Rone<-as.numeric(levels(indrq_l$Rone))[indrq_l$Rone]
indrq_l$AICm<-as.numeric(indrq_l$AICm)
indrq_l$Rone<-as.numeric(indrq_l$Rone)
indrq_l <- indrq_l %>% arrange(reponse,AICm)
# ----

# Graphic representation of model performance
dp<-indrq_l %>% group_by(Sp) %>% filter(!str_detect(mode, "4") & !str_detect(mode, "3")) %>% filter(!str_detect(mode, "Add")) %>% 
  arrange(AICm,.by_group = TRUE) %>% slice(1:20) %>% 
  group_by(Sp,reponse,mode) %>% arrange(AICm,.by_group = TRUE) %>% #slice(1:3) %>%
  ggplot(aes(x=reorder(predict,-AICm),y=AICm,fill = mode))+geom_col(position = "dodge") +
  geom_text(aes(y=AICm, label = AICm),position = position_dodge(width = .9), hjust = 1) + #
  scale_fill_manual(values=colorsS(2)) +
  labs(title="Cubic BSpline Quantile Regression AICm scores", x="Model",y="AICm") +
  theme(axis.text = element_text(size=15,face="bold")) +
  coord_flip()+theme_bw();dp
ggsave(sprintf("%s%s/RQ BSpline/%s_CbsRq_%s_AIC_scores.png",wdgraph,species[sp,1],prgm,species[sp,1]), plot = dp, width = 8, height = 8)
# dev.off()

# Sauvegarde des outputs ----
wb <- loadWorkbook(paste(wdres,prgm,"_BDD",".xlsx", sep="")) # addWorksheet(wb, sheetName = "rql")
writeData(wb, sheet = "CbsRq", x = indrq_l, startCol = 1, startRow = 1,withFilter = FALSE)
writeData(wb, sheet = "CbsRq_coeff", x = smrq_l, startCol = 1, startRow = 1,withFilter = FALSE)
saveWorkbook(wb,file=paste(wdres,prgm,"_BDD",".xlsx", sep=""), overwrite = TRUE)

save.image(file = paste(wdwork,etude,"_RQ_BDD",".RData", sep=""))
} #end if choixttt 1
# ----

#________________________________________________________________
# 2 : Graphiques des SDM 1 à 3 : Deux facteurs ----
if (choixttt==2 | choixttt==12 | choixttt==23 | choixttt==123){
  Mars_SDM<-Mars_dat_sf %>% dplyr::select(c(NINJ,Lon,Lat,Annee,pred_red$Var))
  
for (sp in spe) { #sp=1
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == species[sp,1]),]
  for (sdi in sdm) { #sdi=1
    for (sa in sai) {# sa=1
      k=which(pred_red$Var==sdmlist$x1t[sdi]);k2=which(pred_red$Var==sdmlist$x2t[sdi]);
      rep<-which(reponse$rvar==sdmlist$yt[sdi])
      yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
      zt = "Zone"
      x1t = sprintf("%s%s",pred_red[k,1],saison[sa,1]); x1l = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
      x2t = sprintf("%s%s",pred_red[k2,1],saison[sa,1]); x2l = sprintf("%s%s (%s)",pred_red[k2,2],saison[sa,1],pred_red[k2,3])
      xt<-paste(x1t,x2t,sep=sdmlist$type[sdi]); xl<-paste(x1l,x2l,sep=sdmlist$type[sdi])
      xxx <-c("Beta 0",x1l,x2l,paste(x1l,x2l,sep=sdmlist$type[sdi]))
      dfrq <-df[,c(yt,x1t,x2t,zt)]; dfrq <-as.data.frame(na.omit(dfrq))
      y <- jitter(dfrq[,1]); x1 <- jitter(dfrq[,2]); x2 <- jitter(dfrq[,3]); 
      z <- dfrq[,4]; nbcol <- length(levels(factor(z)))
      
      tryCatch({ # GRAPH SUMMARY ----
      modelq<-eval(parse(text=paste("rq(y~x1",sdmlist$type[sdi],"x2, tau=taus)",sep="")))
      modelq0<-rq(y ~ 1, tau=taus)
      smrq <- summary(modelq, se="boot") #, se="boot"
      psmrq<- plot(smrq,main=xxx,cex=.7,pch=19,lcol=colorsS(3)[3],col=c(colorsS(3)[1],colorsS(3)[2]))
      mlista<-vector(mode = "list", length = dim(psmrq[[1]])[1])
      for (i in 1:dim(psmrq[[1]])[1]){ #i=2
        psmrq_ <- as.data.frame(cbind(taus,t(psmrq[[1]][i,,-(length(taus)+1)])))
        sma <- ggplot(psmrq_,aes_string(x="taus"))+
          geom_ribbon(aes_string(ymin="`lower bd`", ymax="`upper bd`"), fill=colorsS(3)[2], alpha=0.8)+
          geom_point(aes_string(y="coefficients"),col=colorsS(3)[1],size=1)+geom_line(aes_string(y="coefficients"),col=colorsS(3)[1],size=.5, alpha=0.8)+
          geom_hline(yintercept = psmrq[[1]][i,1,(length(taus)+1)],col=colorsS(3)[3])+
          geom_hline(yintercept = psmrq[[1]][i,2:3,(length(taus)+1)],col=colorsS(3)[3],linetype = "dashed",size=.5)+
          labs(title=xxx[i], y=NULL, x=NULL) + 
          theme(plot.title = element_text(hjust = 0.5,size=9,face="bold"),
                axis.text.y = element_text(size=6,face="bold"),axis.text.x = element_text(size=8,face="bold"),
                panel.grid.minor = element_blank())
        assign(paste('sma', i, sep=''), sma)
        mlista[[i]] <- eval(parse(text = paste('sma', i, sep='')))
      }
      nncol = 2; nnrow=ceiling(dim(psmrq[[1]])[1]/nncol)
      sm <- ggarrange(plotlist=mlista, ncol=nncol, nrow=nnrow,legend="bottom",common.legend = TRUE) + bgcolor("white")
      titreG <- sprintf("%s of %s for %s in %s - %s", sdmlist$SDM_desc[sdi],yl,species[sp,2],saison[sa,2],prgm)
      sm <- annotate_figure(sm,top = text_grob(titreG, face = "bold", size = 10))+bgcolor("white");# print(sm)
      ggsave(sprintf("%s%s/RQ BSpline/Summaries/%s_CbsRq_sm_%s_%s_%s_%s.png",
                              wdgraph,species[sp,1],prgm,species[sp,1],yt,sdmlist$name[sdi],saison[sa,2]),plot = sm, width = 10, height = 5)
      },error = function(e) {print(e)})#,finally = {})

      # SDM-NEO GRAPHS ----
      # Experimental points on 3D graphic
      titreG <- sprintf("%s - %s - %s\n%s vs %s",species[sp,2],saison[sa,2],sdmlist$SDM_desc[sdi],yl,xl)
      dp1 <- plot_ly(showlegend=F) %>% 
                  add_trace(x = x1, y = x2, z = y,mode = "markers", type = "scatter3d",
                            marker = list(size = 2, color = "blue", symbol = 104)) %>%
                  layout(title = titreG, scene = list(xaxis = list(title = x1l), yaxis = list(title = x2l), zaxis = list(title = yl)))
      # Surface model of initial conditions
      x1mod <- seq(min(x1),max(x1),length.out=length(x1)); x2mod <- seq(min(x2),max(x2),length.out=length(x1))
      grid<-mesh(x1mod,x2mod)
      dp4<-dp1
      CbsRqMod <- list()
      for (t in 1:length(taus)){ # t=6 length(taus)
        tryCatch({
          modelq<-eval(parse(text=paste("rq(y~x1",sdmlist$type[sdi],"x2, tau=taus[t])",sep="")))
          # Definition of point over model
          rqlim<-coef(modelq)[1] + coef(modelq)[2]*x1 + coef(modelq)[3]*x2 + coef(modelq)[4]*x1*x2
          x1sup<-x1[y>rqlim]; x2sup<-x2[y>rqlim]; ysup<-y[y>rqlim]
          x1inf<-x1[y<=rqlim]; x2inf<-x2[y<=rqlim]; yinf<-y[y<=rqlim]; rm(rqlim)
          # Definition of model surface
          CbsRqMod[[t]]<-coef(modelq)[1] + coef(modelq)[2]*grid[["x"]] + coef(modelq)[3]*grid[["y"]] +
                        coef(modelq)[4]*grid[["x"]]*grid[["y"]]
          CbsRqMod[[t]][which(CbsRqMod[[t]]<0)]<-NA
          
          # 2D static graphic : RASTER
          dp2d <- ggplot() +
            geom_raster(aes(x = grid[["x"]], y = grid[["y"]], fill = CbsRqMod[[t]])) +
            geom_point(aes(x = x1inf, y = x2inf), shape = 8, size = 2, color = colIC, alpha=0.5) +
            geom_point(aes(x = x1sup, y = x2sup), shape = 1, size = 2, color = prettyred, alpha=0.5) +
            labs(title=titreG,x=x1l, y=x2l,fill = paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep="")) +
            scale_fill_gradient(low=blank, high=collist[t]); dp2d # +scale_fill_distiller(palette = "Spectral")
          ggsave(sprintf("%s%s/RQ BSpline/%s_CbsRq2d_%s_%s_%s_%s_%s_%s.png",
                         wdgraph,species[sp,1],prgm,species[sp,1],yt,sdmlist$name[sdi],x1t,x2t,taus[t]),plot = dp2d, width = 10, height = 9)
          
          # 3D plot with model surface and initial condition
          dp3 <- dp1 %>% add_surface(x = grid[["x"]], y = grid[["y"]], z = CbsRqMod[[t]],
                                     opacity = 0.7,
                                     colorscale = list(c(0,1),c(blank,collist[t])),
                                     colorbar=list(title=list(text=paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep="")))) %>%
                        add_trace(x = x1sup, y = x2sup, z = ysup, mode = "markers", type = "scatter3d",
                                  marker = list(size = 2, color = prettyred, symbol = 104))%>%
                        layout(title = titreG)
          htmlwidgets::saveWidget(dp3,sprintf("%s%s/RQ BSpline/Graphes 3D/%s_CbsRq2d_%s_%s_%s_%s_%s_%s.html",
                                  wdgraphEx,species[sp,1],prgm,species[sp,1],yt,sdmlist$name[sdi],x1t,x2t,taus[t]), selfcontained = F, libdir = "lib")
          # 3D plot with all taus models surface added
          dp4 <- dp4 %>% add_surface(x = grid[["x"]], y = grid[["y"]], z = CbsRqMod[[t]],
                                     opacity = 0.5, colorscale = list(c(0,1),c(blank,collist[t])),
                                     colorbar=list(title=list(text=paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep=""))))# %>%
          
          # 3D static graphic
          # data conversion
          # zmod<-CbsRqMod[[t]]; rownames(zmod) = x1; colnames(zmod) = x2
          # zmod<-as.data.frame(zmod) %>% rownames_to_column(var="x1") %>% 
          #   pivot_longer(-x1) %>% mutate(x1=as.numeric(x1),name=as.numeric(name)) %>%
          #   rename(x1z=x1,x2z=name,yz=value)
          # jet.colors <- colorRampPalette( c("blue", "green") ); color <- jet.colors(100)
          # persp(x = grid[["x"]], y = grid[["y"]], z = CbsRqMod[[t]], col = color, phi = 30, theta = -30)
          
          # ----
          
        },error = function(e) {print(e)})#,finally = {})
      } # taus
      # 3D plot with all taus models surface and initial condition
      dp4 <- dp4 %>% add_trace(x = x1sup, y = x2sup, z = ysup, mode = "markers", type = "scatter3d",
                               marker = list(size = 2, color = prettyred, symbol = 104)) 
      htmlwidgets::saveWidget(dp4, sprintf("%s%s/RQ BSpline/Graphes 3D/%s_CbsRq2dStack_%s_%s_%s_%s_%s.html",
                                           wdgraphEx,species[sp,1],prgm,species[sp,1],yt,sdmlist$name[sdi],x1t,x2t), selfcontained = F, libdir = "lib")
      # MARS3D SDM CALCULATION ----
      tau <- sdmlist$SDM_tau[sdi]
      coeff_rql <- smrq_l[smrq_l$Sp==species[sp,1] & smrq_l$reponse==yl & smrq_l$predict==xl & smrq_l$tau==tau, c(7,8)]
      x1M <- pull(Mars_dat_sf, x1t); x2M <- pull(Mars_dat_sf, x2t); 
      Mars_SDM$SDM<-inter(x1M,x2M,coeff_rql$Value[1],coeff_rql$Value[2],
                          coeff_rql$Value[3],coeff_rql$Value[4])
      Mars_SDM$SDM[which(Mars_SDM$SDM<0)]<-NA
      Mars_SDM<-Mars_SDM%>%rename(!!paste(species[sp,1],sdi,saison[sa,1],sep=""):=SDM) # !!'string': to interpret text as variable
      # ----
    }#sai
  }#sdm
  st_write(Mars_SDM, sprintf("%sLayers made/SDM_NEO_CbsRq_%s.shp",wdGIS,species[sp,1]),append=FALSE)
}#spe
  # Sauvegarde des outputs ----
  wb <- loadWorkbook(paste(wdres,prgm,"_BDD",".xlsx", sep="")) # addWorksheet(wb, sheetName = "RQ_line")
  writeData(wb, sheet = "CbsRq_sdmlist", x = sdmlist, startCol = 1, startRow = 1,withFilter = FALSE)
  saveWorkbook(wb,file=paste(wdres,prgm,"_BDD",".xlsx", sep=""), overwrite = TRUE)
  save.image(file = paste(wdwork,etude,"_CbsRq_BDD",".RData", sep=""))
} #end if choixttt 2
# ----

