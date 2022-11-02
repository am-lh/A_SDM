# ____________________________________________________________________________________________________
# GAUSSIAN QUANTILE REGRESSION ANALYSIS WITH KOENKER'S QUANTREG PACKAGE
# ____________________________________________________________________________________________________

rm(list=ls())
#________________________________________________________________
# PACKAGES USED ----
library(readxl) ; library(openxlsx) # Edition d'un fichier Excel
library(tidyverse); library(reshape2); library(rlist) # the one; melt; list.append
library(ggpubr); #library(GGally); # stat_compare_means ;
library(scales); library(RColorBrewer); library(wesanderson); library(grafify); library(colorspace); library(ggsci)# show_col and colors colors colors!
library(quantreg);# library(visreg)
library(splines); library(ggforce)# bs BSplines function
library(plotly); library(plot3D);  # graphiques 3D plot 3D for mesh library(pracma) 
# # GIS Packages
library(sf); library(sfheaders); # st_as_sf ; sf_to_df
# library(tmap) # tmap_mode; for static and interactive maps
library(htmlwidgets) # library(leaflet) # saveWidget ; for interactive maps

#________________________________________________________________
# DEFINITION OF THE GRAPHIC CHARTER ----
theme_set(theme_bw()) # theme_gray() theme_bw() theme_light()
pal_brew <- colorRampPalette(brewer.pal(8, "Spectral")); 
Scale_brew <- function() {scale_colour_brewer(palette="Spectral",aesthetics=c("colour","fill"))}
pal_ggscc <- pal_material("teal"); # show_col(pal_ggscc(9)) # ggsci
pal_cspx <- function(x) {divergingx_hcl(x,palette = "Zissou 1")}; show_col(pal_cspx(6)) # colorspace

blank<-alpha("#dae8ed",0.1); colRQ<-pal_csp(6); 
colInliers<-pal_material("teal")(10)[5] ; colOutliers = pal_material("deep-orange")(10)[5]  # show_col(colOutliers)
# ----

#________________________________________________________________
# WORKING ENVIRONMENT AND LOADING OF BASIC DATA ----
pc <- "C:/Users/lehuen201/Nextcloud/" # "E:/" # 
tsk <- "A_SDM_NEO/"
wdpath <- paste(pc,"Melting Pot/BDD/",tsk,sep=""); 
wdwork <- paste(wdpath,"Matrices/",sep="")
wdgraph <- paste(wdpath,"Graphiques/",sep="")
wdgraphEx<-wdgraph #paste(pc,"Copie-HD/Melting Potes/",tsk,"Graphiques/",sep="")
wdres <- paste(wdpath,"Resultats/",sep="")
wdGIS <- paste(pc,"Melting Pot/SIG/",sep="");
setwd(wdpath) 

#________________________________________________________________
# DEFINITION OF BASIC VARIABLES ----
prgm <- "CSLN" # 1:CSLN 2:Mabes 3:Geco 4:Beaug
etude <- "CSLN_Mars"
load(paste(wdwork,etude,"_BDD.RData",sep=""))
analysis <- "RQ Nonlineaire"
# if exists
# load(paste(wdwork,etude,"_nlRQ_BDD",".RData", sep=""))
# load(paste(wdwork,etude,"_nlRQ2d_BDD",".RData", sep=""))

choixttt<-123 # treatment choices 1 : simple nlrq; 2 : nlrq 2d coeff; 3 : nlrq 2d graphs or combinaisons

# Choice of predictors used according to autocorrelation study, taus...
# Var_choosen<-c("flow_mxd","inunt","flow_m","sal_m","tenfon_m","mudrate_m") #,"bathy"
Var_choosen<-c("flow_mxd","inunt","sal_dtd","temp_m","tenfon_m","mudrate_m")
pred_red <- predict  %>% subset(Var %in% Var_choosen) %>% arrange(match(Var,Var_choosen))
taus <- c(0.5,0.85,0.9,0.95,0.975,0.99) #c(1:19/20,0.975,0.99)
spe <- 1#:nrow(species) # 1:CERED 2:CORVO 3:HEDDI 4:LIMBA 5:PERUL 6:SCRPL
reponse<-reponse[1:2,] # 1:Biomass_gAFDWm2 2:Density_indm2
answ <- 1:nrow(reponse)
sai <- 1#:nrow(saison) # 1:Year 2:Winter 3:Summer
explo <- 1:nrow(pred_red) # 1 for flow_maxd

# Model equation tested ----
# Gaussian equation one factor with Initial conditions (CI) vector
gaussf <- function(x,A,mu,sigma) {A*exp(-((x-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi)))}
gaussfnom = "1D Gaussian"; 
FactA_CI <- c(0.99,0.99,0.99,0.99,0.99,0.99) #(0.99, bathy) modif facteur A regarding results of RQ, same length as pred_red
SDM_tau1<-c(0.99,0.975,0.975,0.975,0.95,0.975);

# Gaussian equation two factors with Initial conditions (CI) vector
# epsilon<-1*10^-6;
# gauss2d <- function(x1,x2,A,mu1,sigma1,mu2,sigma2,epsilon) {A*exp(-(((x1-mu1)^2/(2*sigma1^2))+((x2-mu2)^2/(2*sigma2^2))))+epsilon}
gauss2d <- function(x1,x2,A,mu1,sigma1,mu2,sigma2) {A*exp(-(((x1-mu1)^2/(2*sigma1^2))+((x2-mu2)^2/(2*sigma2^2))))}
gauss2dnom = "2D Gaussian"; 
FactA_CI2d <- 0.99
#----

# Definition of the selected SDMs
name<-c('SDM1','SDM2','SDM3','SDM4','SDM5','SDM6')
SDM_desc<-c('Gaussian two factors','Gaussian two factors','Gaussian two factors','Gaussian two factors','Gaussian two factors','Gaussian two factors')
SDM_tau<-c(0.99,0.975,0.975,0.975,0.975,0.975)
yt<-c('Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2','Biomass_gAFDWm2')
x1t<-c('flow_mxd','flow_mxd','flow_mxd','inunt','sal_dtd','tenfon_m')
x2t<-c('inunt','mudrate_m','tenfon_m','mudrate_m','mudrate_m','mudrate_m')
sdm_choice<-data.frame(name,SDM_desc,SDM_tau,yt,x1t,x2t)
sdm<-1:nrow(sdm_choice)
SDM_tau_all<-c(0.99,0.975,0.975,0.95,0.975,0.975,0.975,0.975,0.99,0.975,0.95,0.95,0.95,0.95,0.9);
# ----



#________________________________________________________________
# 1 : Non-linear One factor Gaussian Quantile Regression ----
if (choixttt==1 | choixttt==12 | choixttt==13 | choixttt==123){
# Indic, coeff et graphes d'un coup ----
funcname <- gaussf  ; funcnom = gaussfnom;
indnlrq_1 <- as.data.frame(matrix(nrow=0,ncol=1));
smnlrq_1 <- as.data.frame(matrix(nrow=0,ncol=1))
indic=1 # SDM increment number
Mars_SDM<-Mars_dat_sf %>% dplyr::select(c(NINJ,Lon,Lat,Zone,Tidal_level,Annee,Period,pred_red$Var))
for (sp in spe) { #sp=1
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == speciesMP$SPCourt[sp]),]
  indnlrq_2 <- as.data.frame(matrix(nrow=0,ncol=1));
  smnlrq_2 <- as.data.frame(matrix(nrow=0,ncol=1))
  smnlrq_tmp <- array(0, dim=c(3, 4, length(taus)))
  for (sa in sai) {# sa=1
    for (rep in answ){ #rep=1
      for (k in explo) { #k=1 explo
        yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
        zt = "Zone"
        xt = sprintf("%s%s",pred_red[k,1],saison[sa,1]); xl = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
        dfrq <-df[,c(yt,xt,zt)]; dfrq <-as.data.frame(na.omit(dfrq))
        y <- jitter(dfrq[,1]); x <- jitter(dfrq[,2]); 
        z <- dfrq[,3]; nbcol <- length(levels(factor(z)))
        lci <- list(A=quantile(y,FactA_CI[k]), mu=median(x), sigma=sd(x)); 
        conditions<- data.frame(name = paste("SDM",indic,sep=""),Sp = speciesMP$SPCourt[sp],
                                reponse=yl, predict = xl, Season=saison[sa,2],
                                yt=yt,x1t=xt,sai=saison[sa,1],
                                SDM_desc = funcnom, type = paste(deparse(body(funcname)[[2]]), collapse = ''))
        
        # tryCatch({ # INDICATORS and COEFFICIENTS EXTRACTION ----
        indnlrq_3<- vector()
        for (t in 1:length(taus)){ #t=5
          tryCatch({
            # nlRQ PARAMETER CONTROL TO AVOID R ABORT WHEN SUMMARY(MODEL) : InitialStepSize=0 (ex=1) ----
            cc<-nlrq.control(maxiter=100, k=2, InitialStepSize = 0, big=1e+20, eps=1e-06, beta=0.97)
            modnlq<-NULL
            modnlq<-nlrq(y~funcname(x,A,mu,sigma), start = lci, tau=taus[t],control=cc, method="BFGS")
            indnlrq_3 <- c(AIC(modnlq),indnlrq_3)
            
            smnlrq_tmp[,,t]<-summary(modnlq)[["coefficients"]]
            Var<-rownames(summary(modnlq)[["coefficients"]])
            colvar<-colnames(summary(modnlq)[["coefficients"]])
            smnlrq <- summary(modnlq) #, se="boot"
            tmp <- as.data.frame(smnlrq[["coefficients"]]);
            tmp$coeff <- row.names(tmp); tmp$tau <- taus[t];
            tmp<-cbind(tmp,conditions)
            smnlrq_2 <- rbind(smnlrq_2,tmp)

            },error = function(e) {print(e)})#,finally = {})
        }
        tmp <- data.frame(AICm=round(median(indnlrq_3,na.rm=TRUE),1));
        tmp<-cbind(tmp,conditions)
        indnlrq_2 <- rbind(indnlrq_2,tmp);
        # },error = function(e) {print(e)})#,finally = {}) 
        #----
      
        tryCatch({ # ONE FACTOR GAUSSIAN SUMMARY GRAPHS ----
          # # ONE FACTOR LINEAR SUMMARY GRAPHS ----
          # titreG <- sprintf("%s for %s in %s - %s",yl,speciesMP$Taxon_SNa[sp],saison[sa,2],prgm)
          # png(file=sprintf("%s%s/RQ Lineaire/Summaries/%s_Rq_sm_%s_%s_%s.png",
          #                  wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],yt,xt),width=600, height=600)
          # sm<-plot(smrq, main=xxx,cex=.7,pch=19,lcol=pal_brew(3)[3],col=c(pal_brew(3)[1],pal_brew(3)[2]),xlab = "tau", ylab = yt, )
          # title(main=paste(titreG,"\n\n\n",sep=""))
          # dev.off()
          
          mlista<-vector(mode = "list", length = dim(smnlrq_tmp)[1])
          for (i in 1:dim(smnlrq_tmp)[1]){ #i=1
            smnlrq_tmp_df <- as.data.frame(cbind(taus,t(smnlrq_tmp[i,,])))
            colnames(smnlrq_tmp_df)<-c("taus",colvar)
            smnlrq_tmp_df$lowerbd<- smnlrq_tmp_df$Value-smnlrq_tmp_df$`Std. Error`
            smnlrq_tmp_df$upperbd<- smnlrq_tmp_df$Value+smnlrq_tmp_df$`Std. Error`
            sma <- ggplot(smnlrq_tmp_df,aes(x=taus))+
              geom_ribbon(aes(ymin=lowerbd, ymax=upperbd), fill=pal_brew(3)[2], alpha=0.8)+
              geom_point(aes_string(y="Value"),col=pal_brew(3)[1],size=1)+geom_line(aes_string(y="Value"),col=pal_brew(3)[1],size=.5, alpha=0.8)+
              geom_hline(yintercept = 0,col=pal_brew(3)[3])+
              labs(title=Var[i], y=NULL, x=NULL) + theme_bw()+
              theme(plot.title = element_text(hjust = 0.5,size=9,face="bold"),
                    axis.text.y = element_text(size=6,face="bold"),axis.text.x = element_text(size=8,face="bold"),
                    panel.grid.minor = element_blank()) + bgcolor("white")
            assign(paste('sma', i, sep=''), sma)
            mlista[[i]] <- eval(parse(text = paste('sma', i, sep='')))
          }
        nncol = 1; nnrow=ceiling(dim(smnlrq_tmp)[1]/nncol)
        sm <- ggarrange(plotlist=mlista, ncol=nncol, nrow=nnrow,legend="bottom",common.legend = TRUE)
        titreG <- sprintf("%s - in %s - %s\n%s of %s vs %s",speciesMP$Taxon_SNa[sp],saison[sa,2],prgm,funcnom,yl,xl)
        sm <- annotate_figure(sm,top = text_grob(titreG, face = "bold", size = 10))#; print(sm)
        ggsave(sprintf("%s%s/RQ Nonlineaire/Summaries/%s_nlRQ_sm_%s_%s_%s.png",
                                wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],yt,xt),plot = sm, width = 5, height = 5)
        },error = function(e) {print(e)})#,finally = {})

        # tryCatch({ # ONE FACTOR GAUSSIAN GRAPHS ----
        dp <-  ggplot(data=dfrq) +
                      geom_point(aes(x=x, y=y), shape = 16, size = 1.5, color = colOutliers, alpha=0.5) +
                      labs(x=xl, y=yl) + Scale_brew(nbcol)#+ ylim(0, 100) # Limit y
        dp<-dp + stat_function(fun=funcname,color=pal_brew(length(taus))[6],linetype = "dotted",size=1,args=lci); #print(dp)
        dp1<-dp;
        for (t in 1:length(taus)){ # t=1
          tryCatch({
          modnlq<-nlrq(y~funcname(x,A,mu,sigma), start = lci, tau=taus[t])
          dp1<-dp1 + stat_function(fun=funcname,color=pal_ggscc(length(taus))[t],size=.7,
                                   args=list(A=coef(modnlq)[1],mu=coef(modnlq)[2],sigma=coef(modnlq)[3]))+
          annotate("text",x=max(x,na.rm=TRUE)*0.9,y=max(y)*((t/length(taus))*0.9),colour=pal_ggscc(length(taus))[t],size=4,fontface=2,hjust=0,
                     label = paste(taus[t]))
          },error = function(e) {print(e)})#,finally = {})
        }
        titreG <- sprintf("%s - in %s - %s",speciesMP$Taxon_SNa[sp],saison[sa,2],prgm)
        dp1<-dp1 + labs(title = titreG, caption = "");print(dp1) #subtitle = sstitreG,
        ggsave(sprintf("%s%s/RQ Nonlineaire/nlRQ Simple/%s_nlRQ_%s_%s_%s.png",
                       wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],yt,xt),plot = dp1, width = 8, height = 8)
        # },error = function(e) {print(e)})#,finally = {})

        # MARS3D SDM CALCULATION ----
        tau <- SDM_tau1[rep]
        coeff_nlrqt <- smnlrq_2[smnlrq_2$name==paste("SDM",indic,sep="") & smnlrq_2$tau==tau,]
        x1M <- pull(Mars_dat_sf, xt);
        Mars_SDM$SDM<-funcname(x=x1M,A=coeff_nlrqt$Value[1],
                               mu=coeff_nlrqt$Value[2],sigma=coeff_nlrqt$Value[3])
        Mars_SDM$SDM[which(Mars_SDM$SDM<0)]<-NA
        Mars_SDM<-Mars_SDM%>%rename(!!paste(speciesMP$SPCourt[sp],indic,saison[sa,1],sep=""):=SDM) # !!'string': to interpret text as variable
        # ----
        
        indic<-indic+1
      } # k explo
    } #reponse
  } # saison
  indnlrq_1 <- rbind(indnlrq_1,indnlrq_2)
  # indnlrq_1$AICm<-as.numeric(levels(indnlrq_1$AICm))[indnlrq_1$AICm] # AICm se retrouve en factor, conversion
  indnlrq_1$AICm<-as.numeric(indnlrq_1$AICm)
  indnlrq_1 <- indnlrq_1 %>% arrange(reponse,AICm)
  smnlrq_1 <- rbind(smnlrq_1,smnlrq_2)
  st_write(Mars_SDM, sprintf("%sLayers made/SDM_NEO_nlrq_%s.shp",wdGIS,speciesMP$SPCourt[sp]),append=FALSE)
} #species
sdmlist<-indnlrq_1[,-sai];sdmlist$SDM_tau<-SDM_tau1;

# Graphic representation of model performance
dp<-indnlrq_1 %>% group_by(Sp) %>% arrange(AICm,.by_group = TRUE) %>% slice(1:20) %>% #filter(str_detect(mode, "Inter"))%>%
  group_by(Sp,reponse) %>% arrange(AICm,.by_group = TRUE) %>% #slice(1:3) %>%
  ggplot(aes(x=reorder(predict,-AICm),y=AICm,fill = reponse))+geom_col(position = "dodge") +
  geom_text(aes(y=AICm, label = AICm),position = position_dodge(width = .9), hjust = 1) + #
  scale_fill_manual(values=pal_brew(7)) +
  labs(title="One Factor Gaussian Quantile Regression AICm scores", x="Model",y="AICm") +
  theme(axis.text = element_text(size=15,face="bold")) +
  coord_flip()+theme_bw();dp
ggsave(sprintf("%s%s/RQ Nonlineaire/%s_%s_nlRQ_AIC_scores.png",wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp]), plot = dp, width = 8, height = 8)

# Sauvegarde des outputs ----
wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep="")) # addWorksheet(wb, sheetName = "RQ_line")
writeData(wb, sheet = "nlRQ", x = indnlrq_1, startCol = 1, startRow = 1,withFilter = FALSE)
writeData(wb, sheet = "nlRQ_coeff", x = smnlrq_1, startCol = 1, startRow = 1,withFilter = FALSE)
saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)
save.image(file = paste(wdwork,etude,"_nlRQ_BDD",".RData", sep=""))
} #end if choixttt 1 ----

#________________________________________________________________
# 2 : Non-linear Two factor Gaussian Quantile Regression ----
if (choixttt==2 | choixttt==12 | choixttt==23 | choixttt==123){
funcname <- gauss2d  ; funcnom = gauss2d;  #
indnlrq_1 <- as.data.frame(matrix(nrow=0,ncol=7));
smnlrq_1 <- as.data.frame(matrix(nrow=0,ncol=7))
indic=1 # SDM increment number
for (sp in spe) { #sp=1
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == speciesMP$SPCourt[sp]),]
  indnlrq_2 <- as.data.frame(matrix(nrow=0,ncol=7));
  smnlrq_2 <- as.data.frame(matrix(nrow=0,ncol=7))
  smnlrq_tmp <- array(0, dim=c(5, 4, length(taus)))
  for (rep in answ){ #rep=1
    for (sa in sai) {# sa=1
      for (k in explo) { # k=1 explo
      for (k2 in explo[-c(1:k)]) { # k2=2 explo[-c(1:k)]
        yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
        zt = "Zone"
        x1t = sprintf("%s%s",pred_red[k,1],saison[sa,1]); x1l = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
        x2t = sprintf("%s%s",pred_red[k2,1],saison[sa,1]); x2l = sprintf("%s%s (%s)",pred_red[k2,2],saison[sa,1],pred_red[k2,3])
        xt<-paste(x1t,x2t,sep="&"); xl<-paste(x1l,x2l,sep=" & ")
        xxx <-c("Beta 0",x1l,x2l,paste(x1l,x2l,sep=" & "))
        dfrq <-df[,c(yt,x1t,x2t,zt)]; dfrq <-as.data.frame(na.omit(dfrq))
        y <- jitter(dfrq[,1]); x1 <- jitter(dfrq[,2]); x2 <- jitter(dfrq[,3]);
        z <- dfrq[,4]; nbcol <- length(levels(factor(z)))
        lci <- list(A=quantile(y,FactA_CI2d),mu1=median(x1),sigma1=sd(x1),mu2=median(x2),sigma2=sd(x2)); #,epsilon=0
        conditions<- data.frame(name = paste("SDM",indic,sep=""),Sp = speciesMP$SPCourt[sp],
                                reponse=yl, predict = xl, Season=saison[sa,2],
                                yt=yt,x1t=x1t,x2t=x2t,sai=saison[sa,1],
                                SDM_desc = funcnom, type = paste(deparse(body(funcname)[[2]]), collapse = ''))
        # tryCatch({ # INDICATORS and COEFFICIENTS EXTRACTION ----
        indnlrq_3<- vector()
        for (t in 1:length(taus)){ #t=3
          tryCatch({
            modnlq<-NULL
            # nlRQ PARAMETER CONTROL TO AVOID R ABORT WHEN SUMMARY(MODEL) : InitialStepSize=0 (ex=1) ----
            cc<-nlrq.control(maxiter=100, k=2, InitialStepSize = 0, big=1e+20, eps=1e-06, beta=0.97)
            modnlq<-nlrq(y~funcname(x1,x2,A,mu1,sigma1,mu2,sigma2), start=lci, tau=taus[t],control=cc, method="BFGS") #,epsilon

            indnlrq_3 <- c(AIC(modnlq),indnlrq_3)
            smnlrq <- summary(modnlq) #, se="boot"
            smnlrq_tmp[,,t]<-smnlrq[["coefficients"]]
            Var<-rownames(smnlrq[["coefficients"]])
            colvar<-colnames(smnlrq[["coefficients"]])
            tmp <- as.data.frame(smnlrq[["coefficients"]]);
            tmp$coeff <- row.names(tmp); tmp$SDM_tau <- taus[t];
            tmp<-cbind(tmp,conditions)
            smnlrq_2 <- rbind(smnlrq_2,tmp)
          },error = function(e) {print(e)})#,finally = {})
        }
        tmp <- data.frame(AICm=round(median(indnlrq_3,na.rm=TRUE),1));
        tmp <- cbind(tmp,conditions)
        indnlrq_2 <- rbind(indnlrq_2,tmp);
        # },error = function(e) {print(e)})#,finally = {})
        #----

        tryCatch({ # GRAPH SUMMARY ----
          mlista<-vector(mode = "list", length = dim(smnlrq_tmp)[1])
          for (i in 1:dim(smnlrq_tmp)[1]){ #i=1
            smnlrq_tmp_df <- as.data.frame(cbind(taus,t(smnlrq_tmp[i,,])))
            colnames(smnlrq_tmp_df)<-c("taus",colvar)
            smnlrq_tmp_df$lowerbd<- smnlrq_tmp_df$Value-smnlrq_tmp_df$`Std. Error`
            smnlrq_tmp_df$upperbd<- smnlrq_tmp_df$Value+smnlrq_tmp_df$`Std. Error`
            sma <- ggplot(smnlrq_tmp_df,aes(x=taus))+
              geom_ribbon(aes(ymin=lowerbd, ymax=upperbd), fill=pal_brew(3)[2], alpha=0.8)+
              geom_point(aes_string(y="Value"),col=pal_brew(3)[1],size=1)+geom_line(aes_string(y="Value"),col=pal_brew(3)[1],size=.5, alpha=0.8)+
              geom_hline(yintercept = 0,col=pal_brew(3)[3])+
              labs(title=Var[i], y=NULL, x=NULL) + theme_bw()+
              theme(plot.title = element_text(hjust = 0.5,size=9,face="bold"),
                    axis.text.y = element_text(size=6,face="bold"),axis.text.x = element_text(size=8,face="bold"),
                    panel.grid.minor = element_blank())
            assign(paste('sma', i, sep=''), sma)
            mlista[[i]] <- eval(parse(text = paste('sma', i, sep='')))
          }
          nncol = 2; nnrow=ceiling(dim(smnlrq_tmp)[1]/nncol)
          sm <- ggarrange(plotlist=mlista, ncol=nncol, nrow=nnrow,legend="bottom",common.legend = TRUE)
          titreG <- sprintf("%s - in %s - %s\n%s of %s vs %s",speciesMP$Taxon_SNa[sp],saison[sa,2],prgm,funcnom,yl,xl)
          sm <- annotate_figure(sm,top = text_grob(titreG, face = "bold", size = 10))+bgcolor("white")#; print(sm)
          ggsave(sprintf("%s%s/RQ Nonlineaire/Summaries/%s_nlRQ2d_sm_%s_%s_%s.png",
                         wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],yt,xt),plot = sm, width = 5, height = 5)
        },error = function(e) {print(e)})#,finally = {})
        indic<-indic+1
      } #k2
      } #k explo
    } # saison
  } #reponse
  indnlrq_1 <- rbind(indnlrq_1,indnlrq_2)
  indnlrq_1$AICm<-as.numeric(indnlrq_1$AICm)
  indnlrq_1 <- indnlrq_1 %>% arrange(reponse,AICm)
  smnlrq_1 <- rbind(smnlrq_1,smnlrq_2)
} #species

# Graphic representation of model performance
dp<-indnlrq_1 %>% group_by(Sp) %>% arrange(AICm,.by_group = TRUE) %>% slice(1:20) %>% #filter(str_detect(mode, "Inter"))%>%
  group_by(Sp,reponse) %>% arrange(AICm,.by_group = TRUE) %>% #slice(1:3) %>%
  ggplot(aes(x=reorder(predict,-AICm),y=AICm,fill = reponse))+geom_col(position = "dodge") +
  geom_text(aes(y=AICm, label = AICm),position = position_dodge(width = .9), hjust = 1) + #
  scale_fill_manual(values=pal_brew(7)) +
  labs(title="Two Factor Gaussian Quantile Regression AICm scores", x="Model",y="AICm") +
  theme(axis.text = element_text(size=15,face="bold")) +
  coord_flip()+theme_bw();dp
ggsave(sprintf("%s%s/RQ Nonlineaire/%s_%s_nlRQ2d_AIC_scores.png",wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp]), plot = dp, width = 8, height = 8)

# Sauvegarde des outputs ----
wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep="")) # addWorksheet(wb, sheetName = "RQ_line")
writeData(wb, sheet = "nlRQ2d", x = indnlrq_1, startCol = 1, startRow = 1,withFilter = FALSE)
writeData(wb, sheet = "nlRQ2d_coeff", x = smnlrq_1, startCol = 1, startRow = 1,withFilter = FALSE)
saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)
save.image(file = paste(wdwork,etude,"_nlRQ2d_BDD",".RData", sep=""))
} #end if choixttt 2 ----

#________________________________________________________________
# 3 : TRACE DES GRAPHES DOUBLE GAUSSIENS ----
if (choixttt==3 | choixttt==13 | choixttt==23 | choixttt==123){
sdmlist<-sdm_choice
  # sdmlist<-indnlrq_1[,-sai]; sdmlist$SDM_tau<-SDM_tau_all; sdm<-1:length(indnlrq_1$name) # to scan all possibilities
funcname <- gauss2d 
Mars_SDM2d<-Mars_dat_sf %>% dplyr::select(c(NINJ,Lon,Lat,Zone,Tidal_level,Annee,pred_red$Var))
for (sp in spe) { #sp=1
  df <- CSLN_Mars[which(CSLN_Mars$SPCourt == speciesMP$SPCourt[sp]),]
  for (sdi in sdm){ #sdi=1
    for (sa in sai) {# sa=1
      rep<-which(reponse$rvar==sdmlist$yt[sdi])
      k=which(pred_red$Var==sdmlist$x1t[sdi]); k2=which(pred_red$Var==sdmlist$x2t[sdi]);
      funcnom = sdmlist$SDM_desc[sdi]; 
      
      yt = reponse[rep,1] ; yl = sprintf("%s (%s)",reponse[rep,2],reponse[rep,3])
      zt = "Zone"
      x1t = sprintf("%s%s",pred_red[k,1],saison[sa,1]); x1l = sprintf("%s%s (%s)",pred_red[k,2],saison[sa,1],pred_red[k,3])
      x2t = sprintf("%s%s",pred_red[k2,1],saison[sa,1]); x2l = sprintf("%s%s (%s)",pred_red[k2,2],saison[sa,1],pred_red[k2,3])
      xt<-paste(x1t,x2t,sep="&"); xl<-paste(x1l,x2l,sep=" & ")
      xxx <-c("Beta 0",x1l,x2l,paste(x1l,x2l,sep=" & "))
      dfrq <-df[,c(yt,x1t,x2t,zt)]; dfrq <-as.data.frame(na.omit(dfrq))
      y <- jitter(dfrq[,1]); x1 <- jitter(dfrq[,2]); x2 <- jitter(dfrq[,3]); 
      z <- dfrq[,4]; nbcol <- length(levels(factor(z)))
      coeff_nlrq <- smnlrq_1[smnlrq_1$name==sdmlist$name[sdi],]
      lci <- list(A=quantile(y,FactA_CI2d),mu1=median(x1),sigma1=sd(x1),mu2=median(x2),sigma2=sd(x2)); #,epsilon=0
    
      # Experimental points on 3D graphic
      titreG <- sprintf("%s - %s - %s\n%s vs %s",speciesMP$Taxon_SNa[sp],saison[sa,2],funcnom,yl,xl)
      dp1 <- plot_ly(showlegend=F) %>% add_trace(x = x1, y = x2, z = y, 
                                    mode = "markers", type = "scatter3d",
                                    marker = list(size = 2, color = "blue", symbol = 104))%>%
                                    layout(title = titreG, scene = list(xaxis = list(title = x1l), 
                                           yaxis = list(title = x2l), zaxis = list(title = yl)))
      # Surface model of initial conditions
      x1mod <- seq(min(x1),max(x1),length.out=length(x1)); x2mod <- seq(min(x2),max(x2),length.out=length(x1))
      grid<-mesh(x1mod,x2mod)
      Zinit<-funcname(x1=grid[["x"]],x2=grid[["y"]],A=lci[[1]],mu1=lci[[2]],sigma1=lci[[3]],mu2=lci[[4]],sigma2=lci[[5]]) #,epsilon=lci[[6]]
      # dp1bis <- dp1 %>% add_surface(x = ~grid[["x"]], y = ~grid[["y"]],z = ~Zinit,
      #               opacity = 0.6,colorscale = list(c(0,1),c(blank,colInliers)),
      #               colorbar=list(title=list(text="Initial conditions"),limits=c(0,max(Zinit))))
      # htmlwidgets::saveWidget(dp1bis,sprintf("%s%s/RQ Nonlineaire/Graphes 3D/%s_nlRQ2d_Initial_%s_%s_%s_%s_%s.html",
      #                                     wdgraphEx,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],sdmlist$name[sdi],yt,xt,taus[t]), selfcontained = F, libdir = "lib")

      dp4<-dp1
      nlrqMod <- list()
      for (t in 1:length(taus)){ # t=6 length(taus)
        coeff_nlrqt <- coeff_nlrq[coeff_nlrq$SDM_tau==taus[t],]
        tryCatch({
        # Definition of point over model
        nlrqlim<-funcname(x1=x1,x2=x2,A=coeff_nlrqt$Value[1],
                          mu1=coeff_nlrqt$Value[2],sigma1=coeff_nlrqt$Value[3],
                          mu2=coeff_nlrqt$Value[4],sigma2=coeff_nlrqt$Value[5]) #,epsilon=coeff_nlrqt$Value[6]
        x1sup<-x1[y>nlrqlim]; x2sup<-x2[y>nlrqlim]; ysup<-y[y>nlrqlim]
        x1inf<-x1[y<=nlrqlim]; x2inf<-x2[y<=nlrqlim]; yinf<-y[y<=nlrqlim]; rm(nlrqlim)
        # Definition of model surface
        nlrqMod[[t]]<-funcname(x1=grid[["x"]],x2=grid[["y"]],A=coeff_nlrqt$Value[1],
                          mu1=coeff_nlrqt$Value[2],sigma1=coeff_nlrqt$Value[3],
                          mu2=coeff_nlrqt$Value[4],sigma2=coeff_nlrqt$Value[5]) #,epsilon=coeff_nlrqt$Value[6]
        nlrqMod[[t]][which(nlrqMod[[t]]<0)]<-NA
        
        # 2D static graphic : RASTER
        dp2d <- ggplot() +
          geom_raster(aes(x = grid[["x"]], y = grid[["y"]], fill = nlrqMod[[t]])) +
          # geom_point(data=dfrq, aes(x = x1, y = x2), shape = 1, size = 2, color = colInliers) +
          geom_point(aes(x = x1inf, y = x2inf), shape = 8, size = 2, color = colInliers, alpha=0.5) +
          geom_point(aes(x = x1sup, y = x2sup), shape = 1, size = 2, color = colOutliers, alpha=0.5) +
          labs(title=titreG,x=x1l, y=x2l,fill = paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep="")) +
          scale_fill_gradient(low=blank, high=colRQ[t],limits = c(0,max(nlrqMod[[t]]))); dp2d # +scale_fill_distiller(palette = "Spectral")
        ggsave(sprintf("%s%s/RQ Nonlineaire/%s_nlRQ2d_%s_%s_%s_%s_%s.png",
                       wdgraph,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],sdmlist$name[sdi],yt,xt,taus[t]),plot = dp2d, width = 10, height = 9)
    
        # 3D plot with model surface and initial condition
        dp3 <- dp1 %>% add_surface(x = grid[["x"]], y = grid[["y"]], z = nlrqMod[[t]],
                                   opacity = 0.7, colorscale = list(c(0,1),c(blank,colRQ[t])),
                                   colorbar=list(title=list(text=paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep=""),
                                                 limits=c(0,max(nlrqMod[[t]]))))) %>%
                       add_trace(x = x1sup, y = x2sup, z = ysup,
                                   mode = "markers", type = "scatter3d",
                                   marker = list(size = 2, color = colOutliers, symbol = 104))%>%
                       add_surface(x = ~grid[["x"]], y = ~grid[["y"]],z = ~Zinit,
                                   opacity = 0.4,colorscale = list(c(0,1),c(blank,colInliers)),
                                   colorbar=list(title=list(text="Initial conditions"),limits=c(0,max(Zinit)))) %>%
                       layout(title = titreG)
        # htmlwidgets::saveWidget(dp3,sprintf("%s%s/RQ Nonlineaire/Graphes 3D/%s_nlRQ2d_%s_%s_%s_%s_%s.html",
        #                                     wdgraphEx,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],sdmlist$name[sdi],yt,xt,taus[t]), selfcontained = F, libdir = "lib")
        # 3D plot with all taus models surface added
        dp4 <- dp4 %>% add_surface(x = grid[["x"]], y = grid[["y"]], z = nlrqMod[[t]],
                                   opacity = 0.6,
                                   colorscale = list(c(0,1),c(blank,colRQ[t])),
                                   colorbar=list(title=list(text=paste("SDM-NEO\n",reponse[rep,3],"\nTau=",taus[t],sep=""),
                                                            limits=c(0,max(nlrqMod[[t]])))))# %>%
        
        # # 3D static graphic
        # data conversion
        # zmod<-nlrqMod[[t]]; rownames(zmod) = x1; colnames(zmod) = x2
        # zmod<-as.data.frame(zmod) %>% rownames_to_column(var="x1") %>% 
        #   pivot_longer(-x1) %>% mutate(x1=as.numeric(x1),name=as.numeric(name)) %>%
        #   rename(x1z=x1,x2z=name,yz=value)
        # jet.colors <- colorRampPalette( c("blue", "green") ); color <- jet.colors(100)
        # persp(x = grid[["x"]], y = grid[["y"]], z = nlrqMod[[t]], col = color, phi = 30, theta = -30)
        
        },error = function(e) {print(e)})#,finally = {})
      } #taus
      # 3D plot with all taus models surface and initial condition
      dp4 <- dp4 %>% add_trace(x = x1sup, y = x2sup, z = ysup,
                      mode = "markers", type = "scatter3d",
                      marker = list(size = 2, color = colOutliers, symbol = 104))%>%
            add_surface(x = ~grid[["x"]], y = ~grid[["y"]],z = ~Zinit,
                        opacity = 0.4,colorscale = list(c(0,1),c(blank,colInliers)),
                        colorbar=list(title=list(text="Initial conditions"),limits=c(0,max(Zinit))))
      # htmlwidgets::saveWidget(dp4, sprintf("%s%s/RQ Nonlineaire/Graphes 3D/%s_nlRQ2dStack_%s_%s_%s_%s.html",
      #                                      wdgraphEx,speciesMP$SPCourt[sp],prgm,speciesMP$SPCourt[sp],sdmlist$name[sdi],yt,xt), selfcontained = F, libdir = "lib")
      
      # MARS3D SDM CALCULATION ----
      tau <- sdmlist$SDM_tau[sdi]
      coeff_nlrqt <- coeff_nlrq[coeff_nlrq$SDM_tau==tau,]
      x1M <- pull(Mars_dat_sf, x1t); x2M <- pull(Mars_dat_sf, x2t); 
      Mars_SDM2d$SDM<-funcname(x1=x1M,x2=x2M,A=coeff_nlrqt$Value[1],
               mu1=coeff_nlrqt$Value[2],sigma1=coeff_nlrqt$Value[3],
               mu2=coeff_nlrqt$Value[4],sigma2=coeff_nlrqt$Value[5]) #,epsilon=coeff_nlrqt$Value[6]
      Mars_SDM2d$SDM[which(Mars_SDM2d$SDM<0)]<-NA
      Mars_SDM2d<-Mars_SDM2d%>%rename(!!paste(speciesMP$SPCourt[sp],sdi,saison[sa,1],sep=""):=SDM) # !!'string': to interpret text as variable
      # ----
    } #sai
  } #sdm
  st_write(Mars_SDM2d, sprintf("%sLayers made/SDM_NEO_nlrq2d_%s.shp",wdGIS,speciesMP$SPCourt[sp]),append=FALSE)
} #spe
# Sauvegarde des outputs ----
wb <- loadWorkbook(paste(wdres,"CSLN_BDD",".xlsx", sep=""))
writeData(wb, sheet = "nlRQ2d_sdmlist", x = sdmlist, startCol = 1, startRow = 1,withFilter = FALSE)
saveWorkbook(wb,file=paste(wdres,"CSLN_BDD",".xlsx", sep=""), overwrite = TRUE)
save.image(file = paste(wdwork,etude,"_nlRQ2d_BDD",".RData", sep=""))
} #end if choixttt 3

#________________________________________________________________
# 4 : TEMPORAL EVOLUTION OF SDM FOR GIP SEINE AVAL ----
# # Definition of temporal periods
# Mars_SDM2d$Period<-"2011-2018"
# Mars_SDM2d$Period[Mars_SDM2d$Annee %in% 1990:1999]<-"1990-1999"
# Mars_SDM2d$Period[Mars_SDM2d$Annee %in% 2000:2010]<-"2000-2010"
# Mars_SDM2d$Period<-as.factor(Mars_SDM2d$Period)
# 
# # TIDAL_LEVEL CALCULATION : should be in mars matlab script ? ----
# Mars_SDM2d$Tidal_level<-NA
# Mars_SDM2d$Tidal_level[Mars_SDM2d$inunt>=0 & Mars_SDM2d$inunt<0.25]<-"Supratidal"
# Mars_SDM2d$Tidal_level[Mars_SDM2d$inunt>=0.25 & Mars_SDM2d$inunt<0.75]<-"Intertidal"
# Mars_SDM2d$Tidal_level[Mars_SDM2d$inunt>=0.75 & Mars_SDM2d$inunt<=1]<-"Infratidal"
# Mars_SDM2d$Tidal_level<-as.factor(Mars_SDM2d$Tidal_level)

# ASSIGNMENT OF AREAS TO SAMPLING POINTS
# ADD NINJ OF MARS_CSV FOR EACH SAMPLING POINT IN FAUNA
# Mars_SDM2d <- st_intersection(Mars_SDM2d,ES_Areas)
# Mars_SDM2d$Zone<-as.factor(Mars_SDM2d$Zone)

titreG <- paste(y," Estuary evolution for ",speciesMP$Taxon_SNa[sp], " - ", prgm,sep="")
bp <- ggplot(Mars_SDM2d, aes(CERED1,fill=Period)) + geom_boxplot() + 
  facet_grid(Tidal_level~Zone) +
  Scale_brew(length(levels(Mars_SDM2d$Period))) + 
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)); bp
bp <- bp + ggtitle(titreG)
ggsave(paste(wdgraph,speciesMP$SPCourt[sp],"/",etude,"_Bxp_",speciesMP$SPCourt[sp],"_",y,".png",sep=""), plot = bp, width = 12, height = 6)

