#| label: load-packages ----

library(conflicted)
library(readxl) ; library(openxlsx); library(beepr) # Edition d'un fichier Excel
library(tidyverse); library(reshape2); library(rlist) # the one; melt; list.append
library(data.table)
library(rstatix); library(Hmisc)  # corr and pvalue calculation
library(boot); library(splines)

library(ggpubr); #library(GGally); # stat_compare_means ;
library(scales); library(RColorBrewer); library(wesanderson); library(grafify); library(colorspace); library(ggsci)# show_col and colors colors colors!
library(quantreg);# library(visreg)
library(plotly); library(plot3D);  # graphiques 3D plot 3D for mesh library(pracma)
library(patchwork)
library(ggridges)
# GIS Packages
library(sf); library(sfheaders); # st_as_sf ; sf_to_df
library(tmap); library(tmaptools)
library(htmlwidgets) # library(leaflet) # saveWidget ; for interactive maps
library(introdataviz) # geom_split_violin # devtools::install_github("psyteachr/introdataviz")

conflict_prefer("transpose", "purrr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("layout", "plotly")

options(scipen=999) # Prevents scientific display of numbers
set.seed(123)

# GENERAL FUNCTIONS ----

#| label: graphchart ----

theme_set(theme_bw(base_size = 12)) # theme_gray() theme_bw() theme_light()
theme_update(plot.title = element_text(hjust = 0.5), # by default titles are centered
             strip.background = element_blank(),
             plot.tag = element_text(size = 12,face = "bold"))

pal_cspx <- function(x) {divergingx_hcl(x,palette = "Zissou 1")}; # show_col(pal_cspx(6)) # colorspace
colRQ<-pal_cspx(4); colBin<-c(pal_cspx(6)[3],pal_cspx(6)[5]); blank<-alpha("#dae8ed",0.1); # show_col(colRQ)
f.Scalc_rq<- function() {scale_colour_manual(values=colRQ)}
f.Scalf_rq2d <- function() {scale_fill_gradientn(colours=colRQ)} # function(x) {scale_fill_material("teal")} # 

f.Scale_brew <- function()
{scale_colour_brewer(palette="Spectral",aesthetics=c("colour","fill"))}

pal_ggscc <- pal_material("teal"); # show_col(pal_ggscc(9)) # ggsci
colSum <- c(pal_ggscc(9)[3],pal_ggscc(9)[7],pal_ggscc(9)[1]) # show_col(colSum)

colInliers<-pal_material("teal")(10)[5] ; colOutliers = pal_material("deep-orange")(10)[5] # show_col(colOutliers)
colSpec <- colorRampPalette(brewer.pal(8, "Spectral")); # show_col(colSpec(2))

pal_csps <- function(x) {sequential_hcl(x,palette = "Batlow")} # show_col(pal_csps(8))
f.Scalc_csps <- function() {scale_color_continuous_sequential(palette = "Batlow")}
f.Scalf_csps <- function() {scale_fill_continuous_sequential(palette = "Batlow")}
# f.Scalc_csps() + f.Scalf_csps()

# colDarj <- function(x) {wes_palette("Cavalcanti1",x, type = "continuous")}; colDarj(8)
# Scalc_misc <- function() {scale_color_discrete_diverging(palette = "Hawaii")};
# Scale_col <- function(x) {scale_colour_manual(values=colDarj(x))}
# Scale_fill <- function() 
#   {scale_fill_manual(palette=colorRampPalette(brewer.pal(11, "Spectral")))}


#| label: functmade ----

f.loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(here::here(fileName))
  mget(ls()[ls() != "fileName"])
}

f.AICc<-function(modelq){ 
    df<-dim(modelq$coefficients)[1]-1 
    AICc<-AIC(modelq)+2*df*(df+1)/(length(modelq$y)-df-1) 
    return(AICc) 
} 
f.AICc_nl<-function(modelq,smrq){
    df<-dim(smrq$coefficients)[1]-1
    aicc<-AIC(modelq)+2*df*(df+1)/(length(modelq$y)-df-1)
    return(aicc)
}

f.replace_outliers <- function(x, coefv=1.5) {
  replace(x, x %in% boxplot.stats(x, coef = coefv)$out, NA) #mean(x)
}

f.gaussf <- function(x,A,mu,sigma) {
  A*exp(-((x-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi)))}
f.gauss2d <- function(x1,x2,A,mu1,sigma1,mu2,sigma2) {
  A*exp(-(((x1-mu1)^2/(2*sigma1^2))+((x2-mu2)^2/(2*sigma2^2))))}

# Function for corr ad pvalue table x was a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" was supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
f.corstars <-function(x, method=c("pearson", "spearman"), 
                    removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing was important.
  mystars <- ifelse(p < .0001, "****", 
                    ifelse(p < .001, "***", 
                           ifelse(p < .01, "**", 
                                  ifelse(p < .05, "*", " "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper/lower triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  }
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- as.data.frame(Rnew)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
}

# ggpairs custom colors on corr
f.corr_col <- function(data, mapping, method="p", use="pairwise", ...){
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  # calculate colour based on correlation value
  colFn <- colorRampPalette(c(colSpec(2)[1], "white", colSpec(2)[2]), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  ggally_text(
    label = as.character(round(corr, 2)),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    ...) + #  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
  
} #wrap(cor_func,method = 'spearman', symbol = "Corr:\n")


#| label: data_set_description ----

f.data_descr <- function (dataset,vars) {
  # dataset=Mars_SDM %>% st_drop_geometry %>% filter(!Zone %in% c("Octeville"))
  # vars ="inunt"
  # dataset= CSLN_df  ; vars ="Biomass_gAFDWm2"
  
  descr<-list(geog=list(stat.test=NULL,plot=NULL),
              time=list(stat.test=NULL,plot=NULL))
  
  # dataset %>% 
  #   levene_test(formula(paste0(vars," ~ Period")))
  
  descr$geog$stat.test <- dataset %>% 
    # group_by(Zone) %>% 
    wilcox_test(formula(paste0(vars," ~ Period"))) %>%
    # levene_test(formula(paste0(vars," ~ Period"))) %>% 
    add_xy_position(x = "Period",group = "Zone",
                    dodge = 0,step.increase = 0.05) %>%
    filter(p.adj.signif!="ns")
  descr$geog$plot <- ggplot(dataset) + 
    # geom_boxplot(aes(x=Period, y = .data[[vars]], fill = Period)) +
    # geom_density_ridges(aes(y = Period , fill=Period,
    #                       x=.data[[vars]]), alpha=0.3) +
    geom_flat_violin(aes(x=Period, fill = Period,
                         y = .data[[vars]] ),
                     alpha = 0.5,size = .1,show.legend = FALSE) +
    stat_summary(aes(x=Period, fill = Period,
                     y = .data[[vars]]),  
                 size=.4, fun.data = mean_cl_boot, 
                 shape=21, stroke =.5, linewidth = .5,
                 inherit.aes = FALSE, show.legend = FALSE) +
    # stat_pvalue_manual(descr$geog$stat.test,label = "p.adj.signif",
    #                    size = 3,tip.length = 0, coord.flip = TRUE) +
    coord_flip() +
    facet_wrap( ~ Zone, 
                nrow = 1, #scales = "free_y", 
                labeller = label_wrap_gen(width = 20)) + 
    labs(y = vars, x="") +
    f.Scale_brew() +
    theme(legend.position="none",
          strip.background = element_blank())
  
  descr$time$stat.test <- dataset %>% 
    group_by(Period) %>%
    wilcox_test(formula(paste0(vars," ~ Zone"))) %>% 
    add_xy_position(x = "Zone",group = "Period",
                    dodge = 0,step.increase = 0.05) %>%
    filter(p.adj.signif!="ns")
  
  descr$time$plot <- ggplot(dataset) + 
    # geom_boxplot(aes(x=Zone, y = .data[[vars]], fill = Zone)) +
    geom_flat_violin(aes(x=Zone, fill = Zone,
                         y = .data[[vars]] ),
                     alpha = 0.5,size = .1,show.legend = FALSE) +
    stat_summary(aes(x=Zone, fill = Zone,
                     y = .data[[vars]]),  
                 size=.4, fun.data = mean_cl_boot, 
                 shape=21, stroke = .5, linewidth = .5,
                 inherit.aes = FALSE, show.legend = FALSE) +
    # stat_pvalue_manual(descr$time$stat.test,label = "p.adj.signif", 
    #                    size = 3,tip.length = 0, coord.flip = TRUE) +
    coord_flip() +
    facet_wrap( ~ Period, 
                nrow = 1, #scales = "free_y", 
                labeller = label_wrap_gen(width = 20)) + 
    labs(y = vars, x="") +
    f.Scale_brew() +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.position="none",
          strip.background = element_blank())
  
  return(descr)
}


#| label: hms_distr ----

f.pl_var_st_board<-function(data,vars){
  # data=Mars_tmp
  # vars=pred_red
  
  descr<-list(pred=NULL, all=NULL)
  descr$pred <-imap(setNames(as.list(vars$Var),vars$whole), ~{
    data %>% 
      ggplot() + 
      geom_density_ridges(aes(y = Period , fill=Period,
                              height = after_stat(count), # to have a comparative height between periods
                              x=.data[[.x]]), 
                          stat="density", alpha=0.3, size = .2) +
      facet_grid(. ~ Zone, 
                 # scales = "free_y",
                 labeller = label_wrap_gen(width = 20)) + 
      labs(x=.y, y="") +
      f.Scale_brew() } )
  
  descr$all<-wrap_plots(descr$pred, ncol=1) & 
    theme(legend.position="none",
          strip.background = element_blank(),
          text=element_text(size=10)) &
    plot_annotation(tag_levels = c('A')) 
  return(descr)
}

#| label: hms_map ----

f.pl_var_map<-function(data, vars){
  # data=Mars_tmp
  # vars=pred_red %>% transpose %>% setNames(.,pred_red$Var); vars=vars[[1]]
  
  pl_map<-list(all=NULL, one=NULL)
  
  pl_map$all <- tm_shape(osm_df) +  
    tm_rgb(saturation = .6, alpha= 0.8) +
    # tm_scale_bar(position = c("left", "bottom"), width = 0.15) + #SCALE
    # tm_compass(position = c("left", "top"), size = 2) +          #NORTH COMPASS
    tm_shape(data) +
    tm_fill(col = vars$Var, midpoint = NA,
            style="cont", palette = "-Spectral", alpha = 0.6) +
    # style="cont", palette = areas_pal, alpha = 0.6) +
    tm_facets(by="Period", nrow = 1) +
    tm_shape(tm_areas) +  
    tm_borders(alpha = 0.3, lwd = 1) +
    tm_layout(
      main.title=vars$whole,
      main.title.size= .8,
      panel.label.bg.color = 'white',
      panel.label.size=.8,
      legend.outside = TRUE, 
      legend.outside.size = .08,
      legend.title.size= .8,
      legend.text.size= .6,
      legend.format=list(text.separator="-"),
      outer.margins=c(0,0,0,0), inner.margins=c(0,0,0,0))
  
  peri<-setNames(as.list(levels(data$Period)), levels(data$Period))
  # Define map for each Period
  pl_map$one <- map(peri, ~{
    tm_shape(osm_df) +  
      tm_rgb(saturation = .6, alpha= 0.8) +
      tm_shape(data %>% 
                 filter(Period %in% .x)) +
      tm_fill(col = vars$Var, 
              title=.x, midpoint = NA,
              style="cont", palette = "-Spectral", alpha = 0.6) +
      # style="cont", palette = areas_pal, alpha = 0.6) +
      tm_shape(tm_areas) +  
      tm_borders(alpha = 0.3, lwd = 1) +
      tm_layout(
        main.title=vars$whole,
        main.title.size= .8,
        panel.label.size=.8,
        legend.outside = TRUE,
        legend.outside.size = .08,
        legend.title.size= .8,
        legend.text.size= .6,
        legend.format=list(text.separator="-"),
        outer.margins=c(0,0,0,0), inner.margins=c(0,0,0,0))
  } ) # map peri
  return(pl_map)
}

# QUANTILE REGRESSION FUNCTIONS ----

#| label: rq_Mod_calc ----

f.rq_Lin <- function(df,biolo,sdmod,taus,typ) {
      # biolo=reponse_l[[1]]; typ=type[[2]]
      # sdmod=pred_red_comb$x1[[1]]
    typetxt<-paste0("RQ",length(sdmod$id),typ[2])
    sdmname<-sprintf("%s_%g%g%g%s%s",typetxt,sp,sai,
                     biolo$repi,
                     paste0(sdmod$id,collapse = ""),
                     ifelse(length(sdmod$id)==1,"0","") )
    yt<- biolo$rvar
    yl<- biolo$compl
    xt<-paste(sdmod$Var, collapse=typ[1])
    xl<-paste(sdmod$whole, collapse=paste0(" ",typ[1]," "))
    dfrq <-df[,c(yt,sdmod$Var)] %>% 
      drop_na() %>%
      setnames(old=c(yt,sdmod$Var),
               new=c("y",paste0("x",1:length(sdmod$Var))) ) %>% 
      mutate(across(everything(),jitter))
    dftmp<-df %>% 
      drop_na(all_of(c(yt,sdmod$Var))) %>%
      select(all_of(c("Zone","Tidal_level","Period","Season","Annee","Mois","NINJ")))
    formumod<-as.formula(paste("y",paste0("x",1:length(sdmod$Var),collapse=typ[1]),sep=" ~ "))
  
    # Model grid definition
    xmod<-map( #x1x2mod
              map(as.list(names(dfrq %>% select(-y))),
                  ~{ dfrq %>% select(-y) %>% 
                    summarise(across(everything(),
                                     list(~min(.x, na.rm = TRUE), ~max(.x, na.rm = TRUE)),
                                     .names = "{.fn}.{.col}")) %>% 
                    select(contains(.x)) %>% 
                    unlist(., use.names=FALSE) } ) %>% 
                setNames(names(dfrq %>% select(-y))) ,
              ~seq(.[1],.[2],length.out=graphfine) ) %>% 
          bind_cols()
    gridx <- expand.grid(xmod)
    if(length(sdmod$Var)==2){
      gridxmat <- mesh(xmod$x1,xmod$x2)
      } else{ gridxmat<-NULL }
    mod_grid<-list(xmod=xmod,gridx=gridx,gridxmat=gridxmat)
    mod_t<-NULL
    
    tryCatch({
      modelq0<-rq(y ~ 1, data=dfrq, tau=taus)
      
      # # --- BOOTSTRAP OF RQ MODEL TO GET GOOD LCI THAT WILL GO INTO THE VALIDATED MODEL ---- 
      # calc_rq<-function(data,indices, formumod, taut){
      #   d <- data[indices,] # d=dfrq
      #   tryCatch({
      #     fit<-rq(formumod, data=d, tau=taut)
      #     return(coef(fit)) },
      #   error = function(e) NULL) }
      # boot_res <- boot(dfrq, calc_rq, R = 1000, 
      #                  formumod=formumod, taut=taus)
      # lci <- boot_res$t %>%
      #   as.data.frame() %>% 
      #   setNames(expand.grid(rownames(boot_res[["t0"]]), taus) %>% 
      #               unite(Var, Var1, Var2) %>% pull) %>% 
      #   summarise(across(.cols=everything(),
      #                    list(mean = ~mean(.,na.rm=TRUE), sd = ~sd(.,na.rm=TRUE)))) %>% 
      #   pivot_longer(everything()) %>% 
      #   separate(name,c("Var","tau","stat"), sep="_") %>% 
      #   mutate(tau=sprintf("tau= %.3f",as.numeric(tau))) %>% 
      #   pivot_wider(names_from = tau, values_from = value) %>% 
      #   filter(stat=="mean") %>% select(-stat) %>% 
      #   column_to_rownames(., var = "Var") %>% 
      #   as.matrix
      # 
      # modelq<-rq(formumod, data=dfrq, tau=taus) # , coef=lci # ca marche pas
      # modelq$coefficients<-lci #what about the fitted values?
      # # --- END OF BOOTSTRAP ---
      
      modelq<-rq(formumod, data=dfrq, tau=taus) # should be commented if bootstrap
      smrq <- summary(modelq, se="boot") %>% 
        setNames(sprintf("%.3f",taus))
      meta <- list(type=typ[1],typetxt=typetxt,
                   Sp=espece,Season=saison[sai,2],
                   reponse=yt,reponset=yl,
                   unit=biolo$runit,
                   predict=xt,predictt=xl,
                   predictl=sdmod$whole,
                   predfile=paste(sdmod$Var, collapse="_"))
      valid <- map(as.list(seq(length(smrq))),
                  ~list(AICc=round(f.AICc(modelq)[.x],1),
                       Rone=round((1 - modelq$rho[.x]/modelq0$rho[.x]),5)) ) %>% 
        setNames(sprintf("%.3f",taus))    
      
      smrq_t <- map2(smrq, taus_l,
                ~{as.data.frame(.x$coefficients) %>%
                  mutate(Var=rownames(.)) %>%
                  relocate(Var) %>%
                  remove_rownames(.) %>%
                  mutate(tau=.y, taust=sprintf("tau= %s",.y),
                         sdmname=sdmname,
                         formula=format(modelq$formula),
                         type=typ[1],typetxt=typetxt,
                         Sp=espece, Season=saison[sai,2],
                         reponse=yt, reponset=yl, unit=biolo$runit,
                         # Predictor1=pred_red[k,1], Predictor2=pred_red[k2,1],
                         predict=xt, predictt=xl,
                         predfile=paste(sdmod$Var, collapse="_")) } )
      smrq_t<-smrq_t %>% map2(valid,., append)
  
      # Out of limits points calculation
      rqlim_t <- modelq %>% 
        predict %>% 
        as.data.frame %>%
        rename_all(~sprintf("%.3f",taus)) %>% 
        mutate(across(everything(), function(x){replace(x, which(x<0), NA)})) %>% 
        bind_cols(tibble(dfrq)) %>% 
        bind_cols(tibble(dftmp)) %>% 
       melt(id.vars=c(names(dfrq),names(dftmp)),
                      variable.name = "tau",value.name = "RqLim") %>% 
        mutate(sdmname=sdmname,
               taust=sprintf("tau= %s",tau),
               tau=as.numeric(as.character(tau))) %>% 
        mutate(status=case_when(y>RqLim ~"over",y<=RqLim ~"under" ))
  
      # Surface model calculation
      mod_pred_t <- as.data.frame(predict(modelq,newdata=gridx)) %>%
        rename_all(~sprintf("%.3f",taus)) %>% 
        mutate(across(everything(), 
                      function(x){replace(x, which(x<0), NA)})) %>% 
        bind_cols(gridx) 
      mod_pred_t <- melt(mod_pred_t,id.vars=names(dfrq %>% select(-y)),
                      variable.name = "tau",value.name = "RqMod") %>% 
        mutate(sdmname=sdmname,
               taust=sprintf("tau= %s",tau),
               tau=as.numeric(as.character(tau)))

    mod_t<-list(sdmname=sdmname, meta=meta,
                modelq=modelq, smrq=smrq, smrq_t=smrq_t,
                rqlim_t=rqlim_t, mod_grid=mod_grid, mod_pred_t=mod_pred_t)
    },error = function(e) {print(e)})#,finally = {}) # end try catch
    return(mod_t)
  } # end func

f.rq_nLin <- function(df,biolo,sdmod,typ,taus) {
  # biolo=reponse_l[[1]]; typ=type[[1]]
  # sdmod=pred_red_comb[[1]][[2]]
  typetxt<-paste0("RQ",length(sdmod$id),typ[2])
  sdmname<-sprintf("%s_%g%g%g%s%s",typetxt,sp,sai,
                   biolo$repi,
                   paste0(sdmod$id,collapse = ""),
                   ifelse(length(sdmod$id)==1,"0","") )
  
  yt<- biolo$rvar
  yl<- biolo$compl
  xt<-paste(sdmod$Var, collapse=typ[1])
  xl<-paste(sdmod$whole, collapse=paste0(" ",typ[1]," "))
  dfrq <-df[,c(yt,sdmod$Var)] %>% 
    drop_na() %>%
    setnames(old=c(yt,sdmod$Var),
             new=c("y",paste0("x",1:length(sdmod$Var))) ) %>% 
    mutate(across(everything(),jitter))
  dftmp<-df %>% 
    drop_na(all_of(c(yt,sdmod$Var))) %>%
    select(all_of(c("Zone","Tidal_level","Period","Season","Annee","Mois","NINJ")))
  formumod<-as.formula(ifelse(length(sdmod$id)==1,"y~f.gaussf(x1,A,mu,sigma)",
                              "y~f.gauss2d(x1,x2,A,mu1,sigma1,mu2,sigma2)"))
  
  lci <- if(length(sdmod$id)==1){
    list(A=quantile(dfrq$y,FactA_CI), mu=median(dfrq$x1), sigma=sd(dfrq$x1))
  } else {
    list(A=quantile(dfrq$y,FactA_CI), mu1=median(dfrq$x1), sigma1=sd(dfrq$x1),
         mu2=median(dfrq$x2), sigma2=sd(dfrq$x2))}
  
  # Model grid definition
  xmod<-map( #x1x2mod
    map(as.list(names(dfrq %>% select(-y))),
        ~{ dfrq %>% select(-y) %>% 
            summarise(across(everything(),
                             list(~min(.x, na.rm = TRUE), ~max(.x, na.rm = TRUE)),
                             .names = "{.fn}.{.col}")) %>% 
            select(contains(.x)) %>% 
            unlist(., use.names=FALSE) } ) %>% 
      setNames(names(dfrq %>% select(-y))) ,
    ~seq(.[1],.[2],length.out=graphfine) ) %>% 
    bind_cols()
  gridx <- expand.grid(xmod)
  if(length(sdmod$Var)==2){
    gridxmat <- mesh(xmod$x1,xmod$x2)
  } else{ gridxmat<-NULL }
  mod_grid<-list(xmod=xmod,gridx=gridx,gridxmat=gridxmat)
  
  modt<-map(as.list(taus), \(taut) { #taut=taus[[4]]
    tryCatch({
      # nlRQ PARAMETER CONTROL TO AVOID R ABORT WHEN SUMMARY(MODEL) : InitialStepSize=0 (ex=1)
      cc<-nlrq.control(maxiter=100, k=2, 
                       InitialStepSize = 0, 
                       big=1e+20, eps=1e-06, beta=0.97)
      modelq0<-rq(y ~ 1, data=dfrq, tau=taut)
      
      # # --- BOOTSTRAP OF RQ MODEL TO GET GOOD LCI THAT WILL GO INTO THE VALIDATED MODEL ---- 
      # calc_rq<-function(data,indices, formumod, lci, taut, cc){
      #   d <- data[indices,] # d=dfrq
      #   tryCatch({
      #     fit<-nlrq(formumod, data=d, 
      #                start = lci, tau=taut,
      #                control=cc, method="L-BFGS-B")
      #     return(coef(fit)) },
      #   error = function(e) NA ) }
      # boot_res <- boot(dfrq, calc_rq, R = 1000, 
      #                  formumod=formumod, lci=lci, taut=taut, cc=cc)
      # lci <- boot_res$t %>%
      #   as.data.frame() %>% 
      #   setNames(names(lci)) %>% 
      #   summarise(across(.cols=everything(),list(mean = ~mean(.,na.rm=TRUE), sd = ~sd(.,na.rm=TRUE)))) %>% 
      #   rename_with(.,~gsub("_mean","",.x)) %>% 
      #   select(names(lci)) %>% 
      #   as.list
      # # --- END OF BOOTSTRAP ---
      
      modelq<-nlrq(formumod, data=dfrq, 
                   start = lci, tau=taut,
                   control=cc, method="L-BFGS-B")
      smrq <- summary(modelq)  #, se="boot"
      meta <- list(type=typ[1],typetxt=typetxt,
                   Sp=espece,Season=saison[sai,2],
                   reponse=yt,reponset=yl,
                   unit=biolo$runit,
                   predict=xt,predictt=xl,
                   predictl=sdmod$whole,
                   predfile=paste(sdmod$Var, collapse="_"),
                   lci=lci)
      valid <- list(AICc=round(f.AICc_nl(modelq,smrq),1),
                    Rone=round((1 - modelq$m$rho/modelq0$rho),5))    
      smrq_t <- smrq$coefficients %>% 
        as.data.frame %>%
        mutate(Var=rownames(.)) %>%
        relocate(Var) %>%
        remove_rownames(.) %>%
        mutate(tau=taut, taust=sprintf("tau= %s",tau),
               sdmname=sdmname,
               formula=paste0(format(formumod),collapse=""), #format(modelq$formula),
               type=typ[1],typetxt=typetxt,
               Sp=espece, Season=saison[sai,2],
               reponse=yt, reponset=yl, unit=biolo$runit,
               # Predictor1=pred_red[k,1], Predictor2=pred_red[k2,1],
               predict=xt, predictt=xl,
               predfile=paste(sdmod$Var, collapse="_"))
      smrq_t<-smrq_t %>% append(valid)
      # Out of limits points calculation
      rqlim_t <- modelq %>% 
        predict %>% 
        as.data.frame %>%
        rename_all(~sprintf("%s",taut)) %>% 
        mutate(across(everything(), function(x){replace(x, which(x<0), NA)})) %>% 
        bind_cols(tibble(dfrq)) %>% 
        bind_cols(tibble(dftmp)) %>% 
        melt(id.vars=c(names(dfrq),names(dftmp)),
             variable.name = "tau",value.name = "RqLim") %>% 
        mutate(sdmname=sdmname,
               taust=sprintf("tau= %s",tau),
               tau=as.numeric(as.character(tau))) %>% 
        mutate(status=case_when(y>RqLim ~"over",y<=RqLim ~"under" ))
      # Surface model calculation
      mod_pred_t <- as.data.frame(predict(modelq,newdata=gridx)) %>%
        rename_all(~sprintf("%s",taut)) %>% 
        mutate(across(everything(), 
                      function(x){replace(x, which(x<0), NA)})) %>% 
        bind_cols(gridx) 
      mod_pred_t <- mod_pred_t %>% 
        rename_all(~c("RqMod",names(dfrq %>% select(-y))) ) %>% 
        mutate(sdmname=sdmname,
               tau=taut,
               taust=sprintf("tau= %s",tau))
      
      modl<-list(sdmname=sdmname, meta=meta,
                 modelq=modelq, smrq=smrq, smrq_t=smrq_t,
                 rqlim_t=rqlim_t, mod_grid=mod_grid, mod_pred_t=mod_pred_t)
    },error = function(e) {print(e)})#,finally = {}) # end try catch
  } ) # end map taus
  
  modt <- modt %>%
    setNames(sprintf("%.3f",taus)) %>% 
    transpose %>% 
    modify_at(.,c("sdmname","meta","mod_grid"),~.x[[1]]) %>% 
    modify_at(.,c("smrq_t","rqlim_t","mod_pred_t"),~.x %>% bind_rows())
  return(modt)
} # end func

f.rq_Bsp <- function(df,biolo,sdmod,taus,typ) {
  # biolo=reponse_l[[1]]; typ=type[[1]]
  # sdmod=pred_red_comb$x2[[1]]
  typetxt<-paste0("RQ",length(sdmod$id),typ[2])
  sdmname<-sprintf("%s_%g%g%g%s%s",typetxt,sp,sai,
                   biolo$repi,
                   paste0(sdmod$id,collapse = ""),
                   ifelse(length(sdmod$id)==1,"0","") )
  
  yt<- biolo$rvar
  yl<- biolo$compl
  xt<-paste(sdmod$Var, collapse=typ[1])
  xl<-paste(sdmod$whole, collapse=paste0(" ",typ[1]," "))
  dfrq <-df[,c(yt,sdmod$Var)] %>% 
    drop_na() %>%
    setnames(old=c(yt,sdmod$Var),
             new=c("y",paste0("x",1:length(sdmod$Var))) ) %>% 
    mutate(across(everything(),jitter))
  dftmp<-df %>% 
    drop_na(all_of(c(yt,sdmod$Var))) %>%
    select(all_of(c("Zone","Tidal_level","Period","Season","Annee","Mois","NINJ")))
  formumod<-as.formula(paste("y",paste0("bs(x",1:length(sdmod$Var),
                                        ",degree=3,knots=median(x",1:length(sdmod$Var),"))",
                                        collapse=typ[1]),sep=" ~ "))
  
  # Model grid definition
  xmod<-map( #x1x2mod
    map(as.list(names(dfrq %>% select(-y))),
        ~{ dfrq %>% select(-y) %>% 
            summarise(across(everything(),
                             list(~min(.x, na.rm = TRUE), ~max(.x, na.rm = TRUE)),
                             .names = "{.fn}.{.col}")) %>% 
            select(contains(.x)) %>% 
            unlist(., use.names=FALSE) } ) %>% 
      setNames(names(dfrq %>% select(-y))) ,
    ~seq(.[1],.[2],length.out=graphfine) ) %>% 
    bind_cols()
  gridx <- expand.grid(xmod)
  if(length(sdmod$Var)==2){
    gridxmat <- mesh(xmod$x1,xmod$x2)
  } else{ gridxmat<-NULL }
  mod_grid<-list(xmod=xmod,gridx=gridx,gridxmat=gridxmat)
  mod_t<-NULL
  
  tryCatch({
    modelq0<-rq(y ~ 1, data=dfrq, tau=taus)
    
    # # --- BOOTSTRAP OF RQ MODEL TO GET GOOD LCI THAT WILL GO INTO THE VALIDATED MODEL ---- 
    # calc_rq<-function(data,indices, formumod, taut){
    #   d <- data[indices,] # d=dfrq
    #   tryCatch({
    #     fit<-rq(formumod, data=d, tau=taut)
    #     return(coef(fit)) },
    #   error = function(e) NULL) }
    # boot_res <- boot(dfrq, calc_rq, R = 1000, 
    #                  formumod=formumod, taut=taus)
    # lci <- boot_res$t %>%
    #   as.data.frame() %>% 
    #   setNames(expand.grid(rownames(boot_res[["t0"]]), taus) %>% 
    #               unite(Var, Var1, Var2) %>% pull) %>% 
    #   summarise(across(.cols=everything(),list(mean = ~mean(.,na.rm=TRUE), sd = ~sd(.,na.rm=TRUE)))) %>% 
    #   pivot_longer(everything()) %>% 
    #   separate(name,c("Var","tau","stat"), sep="_") %>% 
    #   mutate(tau=sprintf("tau= %.3f",as.numeric(tau))) %>% 
    #   pivot_wider(names_from = tau, values_from = value) %>% 
    #   filter(stat=="mean") %>% select(-stat) %>% 
    #   column_to_rownames(., var = "Var") %>% 
    #   as.matrix
    # 
    # modelq<-rq(formumod, data=dfrq, tau=taus) # , coef=lci # ca marche pas
    # modelq$coefficients<-lci #what about the fitted values?
    # # --- END OF BOOTSTRAP ---
    
    modelq<-rq(formumod, data=dfrq, tau=taus) # should be commented if bootstrap
    smrq <- summary(modelq, se="boot") %>% 
      setNames(sprintf("%.3f",taus))
    meta <- list(type=typ[1],typetxt=typetxt,
                 Sp=espece,Season=saison[sai,2],
                 reponse=yt,reponset=yl,
                 unit=biolo$runit,
                 predict=xt,predictt=xl,
                 predictl=sdmod$whole,
                 predfile=paste(sdmod$Var, collapse="_"))
    valid <- map(as.list(seq(length(smrq))),
                 ~list(AICc=round(f.AICc(modelq)[.x],1),
                       Rone=round((1 - modelq$rho[.x]/modelq0$rho[.x]),5)) ) %>% 
      setNames(sprintf("%.3f",taus))    
    
    smrq_t <- map2(smrq, as.list(taus),
                   ~{as.data.frame(.x$coefficients) %>%
                       mutate(Var=rownames(.)) %>%
                       relocate(Var) %>%
                       remove_rownames(.) %>%
                       mutate(tau=.y, taust=sprintf("tau= %s",.y),
                              sdmname=sdmname,
                              formula=paste0(format(modelq$formula),collapse=""),
                              type=typ[1],typetxt=typetxt,
                              Sp=espece, Season=saison[sai,2],
                              reponse=yt, reponset=yl, unit=biolo$runit,
                              # Predictor1=pred_red[k,1], Predictor2=pred_red[k2,1],
                              predict=xt, predictt=xl,
                              predfile=paste(sdmod$Var, collapse="_")) } )
    smrq_t<-smrq_t %>% map2(valid,., append)
    
    # Out of limits points calculation
    rqlim_t <- modelq %>% 
      predict %>% 
      as.data.frame %>%
      rename_all(~sprintf("%.3f",taus)) %>% 
      mutate(across(everything(), function(x){replace(x, which(x<0), NA)})) %>% 
      bind_cols(tibble(dfrq)) %>% 
      bind_cols(tibble(dftmp)) %>% 
      melt(id.vars=c(names(dfrq),names(dftmp)),
           variable.name = "tau",value.name = "RqLim") %>% 
      mutate(sdmname=sdmname,
             taust=sprintf("tau= %s",tau),
             tau=as.numeric(as.character(tau))) %>% 
      mutate(status=case_when(y>RqLim ~"over",y<=RqLim ~"under" ))
    
    # Surface model calculation
    mod_pred_t <- as.data.frame(predict(modelq,newdata=gridx)) %>%
      rename_all(~sprintf("%.3f",taus)) %>% 
      mutate(across(everything(), 
                    function(x){replace(x, which(x<0), NA)})) %>% 
      bind_cols(gridx) 
    mod_pred_t <- melt(mod_pred_t,id.vars=names(dfrq %>% select(-y)),
                       variable.name = "tau",value.name = "RqMod") %>% 
      mutate(sdmname=sdmname,
             taust=sprintf("tau= %s",tau),
             tau=as.numeric(as.character(tau)))
    mod_t<-list(sdmname=sdmname, meta=meta,
                modelq=modelq, smrq=smrq, smrq_t=smrq_t,
                rqlim_t=rqlim_t, mod_grid=mod_grid, mod_pred_t=mod_pred_t)
  },error = function(e) {print(e)})#,finally = {}) # end try catch
  return(mod_t)
} # end func

#| label: rq_Mod_mars ----
f.rq_Mod_mars <- function(rqMod,MarsSDM,taus) {
      # MarsSDM=Mars_SDM_per
      # rqMod=rq_Mod_sel[[1]][[1]]
    typetxt<-rqMod$meta$typetxt
    vars<-unlist(strsplit(rqMod$meta$predict,rqMod$meta$type,fixed = TRUE))
    # sdmname<-rqMod$sdmname

    xmars <- MarsSDM %>%
      st_drop_geometry() %>% # as.data.frame() %>%
      ungroup %>%
      select(all_of(vars)) %>%
      rename_all(~paste0("x",1:length(vars)) )

    MarsSDM<-MarsSDM %>%
      select(all_of(c("NINJ","Lon","Lat","Zone","Tidal_level","Period","Annee",vars))) #%>%
    RqModMars <- as.data.frame(predict(rqMod$modelq,newdata=xmars)) %>%
      rename_all(~sprintf("t%.3f",taus)) %>%
      # rename_all(~sprintf("%s_t%.3f",sdmname,taus)) %>%
      mutate(across(everything(), function(x){replace(x, which(x<0), NA)}))
    MarsSDM<-MarsSDM %>% cbind(RqModMars)
    return(MarsSDM)
}

f.rq_Mod_mars_nl <- function(rqMod,MarsSDM,taus) {
  # MarsSDM=Mars_SDM
  # rqMod=rq_Mod_sel[[1]][[2]]
  typetxt<-rqMod$meta$typetxt
  vars<-unlist(strsplit(rqMod$meta$predict,rqMod$meta$type,fixed = TRUE))
  taurel<-names(rqMod$modelq)
  
  xmars <- MarsSDM %>%
    st_drop_geometry() %>% # as.data.frame() %>%
    ungroup %>%
    select(all_of(vars)) %>%
    rename_all(~paste0("x",1:length(vars)))
  
  MarsSDM<-MarsSDM %>%
    select(all_of(c("NINJ","Lon","Lat","Zone","Tidal_level","Period","Annee",vars))) #%>%
  
  RqModMars <- as.data.frame(map(rqMod$modelq,
                                 ~predict(.x,,newdata=xmars) )) %>%
    rename_all(~sprintf("t%s",taurel)) %>%
    mutate(across(everything(), function(x){replace(x, which(x<0), NA)}))
  
  MarsSDM<-MarsSDM %>% cbind(RqModMars)
  return(MarsSDM)
}

#| label: rq_Mod_mars_si ----

f.suit_index<-function(rqmodmars,rqMod,taus){
   # rqmodmars=nlrqdata$rq_Mod_mars_sel[[1]][[1]]
   # rqMod=nlrqdata$rq_Mod_sel[[1]][[1]]

  taul<-as.list(
    rqMod$mod_pred_t %>% 
      group_by(tau) %>% 
      summarise(max=max(RqMod,na.rm=TRUE)) %>% 
      pull) %>%  setNames(sprintf("t%.3f",taus))

  imap_dfc(taul, 
          ~ rqmodmars %>%
            st_drop_geometry %>% 
            transmute(!! str_c('si_', .y) := .data[[.y]]/.x)) %>%
    bind_cols(rqmodmars, .)
}

# GRAPHIC FUNCTIONS ----

#| label: rq_Mod_plot_sum ----

f.pl_rq_Mod_sum_rq<-function(rqMod,titleG){
    # biolo=reponse_l[[1]]
    # rqMod=rq_Mod_sel[[biolo$rdescr]][[1]]
  typetxt<-substr(rqMod$sdmname, 1, 6) 
  png(file=sprintf("%s/%s/%s_%s_sm.tiff",
                   graph_path, typetxt, espece, rqMod$sdmname),
      width=3600, height=3600, res=400)
  plot(rqMod$smrq, main=c("Beta 0",rqMod$meta$predictt),cex=.7,pch=19,lcol=colSum[1],col=colSum[2:3])
  title(sub = titleG)
  dev.off()
  
  # rqMod$smrq_t$lowerbd<- rqMod$smrq_t$Value-rqMod$smrq_t$`Std. Error`
  # rqMod$smrq_t$upperbd<- rqMod$smrq_t$Value+rqMod$smrq_t$`Std. Error`
  # ggplot(rqMod$smrq_t,aes(x=tau,y=Value, group = 1)) +
  #   geom_ribbon(aes(ymin=lowerbd, ymax=upperbd), 
  #               fill=colSum[3], alpha=0.8) +
  #   geom_point(col=colSum[2],size=1) +
  #   geom_line(col=colSum[2],linewidth=.5, alpha=0.8) +
  #   facet_wrap(~Var,scales="free",labeller = label_wrap_gen(width = 30)) +
  #   labs(title = titleG) +
  #   theme(text=element_text(size=8),
  #         strip.text = element_text(size=6))

}

f.pl_rq_Mod_sum<-function(rqMod,titleG){
  # rqMod=rq_Mod_sel[[1]][[1]]
  rqMod$smrq_t$lowerbd<- rqMod$smrq_t$Value-rqMod$smrq_t$`Std. Error`
  rqMod$smrq_t$upperbd<- rqMod$smrq_t$Value+rqMod$smrq_t$`Std. Error`
  ggplot(rqMod$smrq_t,aes(x=tau,y=Value, group = 1)) +
    geom_ribbon(aes(ymin=lowerbd, ymax=upperbd), 
                fill=colSum[3], alpha=0.8) +
    geom_point(col=colSum[2],size=1) +
    geom_line(col=colSum[2],linewidth=.5, alpha=0.8) +
    facet_wrap(~Var,scales="free",labeller = label_wrap_gen(width = 30)) +
    labs(title = titleG) +
    theme(text=element_text(size=8),
          strip.text = element_text(size=6))
}

#| label: rq_Mod_plot_aic ----

f.pl_rq_Mod_aic<-function(smrql,titleG){
# smrql <- rq_Mod_x1 %>%
#   map_depth(.,2, ~map_dfr(.$smrq_t, ~.) ) %>%
#   map(., ~map_dfr(., ~.) ) %>%
#   map_dfr(., ~.)
  smrql %>%
    group_by(tau) %>%
    ggplot(aes(x=reorder(predict,-Delta_AICc),
               y=Delta_AICc,
               color = factor(tau),shape=reponse))+
    geom_point(alpha=0.8, size= 5) +
    scale_color_manual(values=colRQ) +
    labs(title=paste("Delta AICc scores for ",titleG,sep=""),
         x="Model",y="Delta AICc",
         color="Quantile",shape="Biologic") +
    theme(axis.text = element_text(size=10,face="bold")) +
    coord_flip()
}


#| label: prefig-poplot_boot ----

f.rq_po_plot <- function(rqMod) { # without bootstrap
    # rqMod=rq_Mod_sel[[1]][[1]]
  
   rqMod$rqlim_t %>% 
    mutate(taustx=as.character(tau)) %>% 
    nest(data=-taustx) %>% 
    mutate(data=map(data, ~{.x %>% 
                    arrange(.,RqLim) %>% 
                    mutate(.,points_bin = ntile(RqLim, n=10)) %>% 
                    group_by(points_bin) %>% 
                    summarise(Observed = quantile(y, unique(tau), na.rm=TRUE),
                              Predict = mean(RqLim, na.rm=TRUE),
                              tau = unique(tau)) %>%
                    mutate(tau=as.factor(tau)) 
      })) %>% 
    unnest(data) %>% 
    ggplot(aes(x=Observed,y=Predict,col=tau,fill=tau,shape=tau)) +
      geom_point(size=2) +
      geom_abline(slope=1,intercept=0) +
      geom_smooth(method = "lm", formula= y~x, aes(fill = tau), alpha=0.05) +
      stat_cor(aes(label=after_stat(rr.label)), 
               label.x.npc = "left", label.y.npc = "top", #label.x.npc = .6,label.y.npc=.3,
               show.legend = FALSE) +
      # stat_regline_equation(label.x.npc = .7,label.y.npc=.4,show.legend = FALSE) +
      xlim(c(0, max(max(data$Observed,na.rm=TRUE),max(data$Predicted,na.rm=TRUE)) )) +
      ylim(c(0, max(max(data$Observed,na.rm=TRUE),max(data$Predicted,na.rm=TRUE)) )) +
      labs(title=sprintf("%s (%s)\n%s", 
                            rqMod$meta$typetxt, 
                            rqMod$meta$unit,
                            rqMod$meta$predictt)) +
      theme(legend.position="bottom",
            aspect.ratio=1,
            text=element_text(size=10)#, plot.title=element_text(size=10)
            ) +
      f.Scale_brew()
}

f.rq_po_plot <- function(data, metas) { # for the bootstrap
    # data=rq_mod_boot_tmp
    # metas= rq_Mod_sel[[1]][[1]]$meta
   data %>% 
    ggplot(aes(x=Observed, y=Predicted, col=tau, fill=tau, shape=tau)) +
      # geom_pointrange(data=data,aes(#xmin=Observed-Observed_se,
      #                   #xmax=Observed+Observed_se,
      #                   ymin=Predicted-Predicted_se,
      #                   ymax=Predicted+Predicted_se),
      #               colour="gray50") +
      geom_point(size=2) +
      geom_abline(slope=1,intercept=0) +
      geom_smooth(method = "lm", formula= y~x, aes(fill = tau), alpha=0.05) +
      stat_cor(aes(label=after_stat(rr.label)), 
               label.x.npc = "left", label.y.npc = "top",
               show.legend = FALSE) +
      # stat_regline_equation(label.x.npc = .7,label.y.npc=.4,show.legend = FALSE) +
      xlim(c(0, max(max(data$Observed,na.rm=TRUE),max(data$Predicted,na.rm=TRUE)) )) +
      ylim(c(0, max(max(data$Observed,na.rm=TRUE),max(data$Predicted,na.rm=TRUE)) )) +
      labs(title=sprintf("%s (%s)\n%s",
                            metas$typetxt,
                            metas$unit,
                            metas$predictt)) +
      theme(legend.position="bottom",
            aspect.ratio=1,
            text=element_text(size=10) #, plot.title=element_text(size=10)
            ) +
      f.Scale_brew()
}

f.calc_boot <- function(data, indices, taut) {
  d <- data[indices,]
  Observed <- quantile(d$y, taut, na.rm=TRUE) 
  Predicted <- mean(d$RqLim, na.rm=TRUE)
  return(c(Observed, Predicted))
}

f.bootstrap_table <- function(data, taut, bin, n_bootstraps) {
    # data=rq_mod_tmp; taut=0.95; bin=9; n_bootstraps=10
  subset_data <- data[data$tau == taut & 
                      data$points_bin == bin,]

  bootstrap_results <- boot(subset_data, f.calc_boot, R = n_bootstraps, taut=taut)
  Observed <- bootstrap_results$t[,1]
  Predicted <- bootstrap_results$t[,2]
  Observed_ci <- boot.ci(bootstrap_results, type = "basic")$basic[2:3]
  Predicted_ci <- boot.ci(bootstrap_results, type = "basic")$basic[4:5]
  
  results_table <- data.frame(
    Observed = mean(Observed), Observed_se = sd(Observed), Observed_ci,
    Predicted = mean(Predicted), Predicted_se = sd(Predicted), Predicted_ci,
    tau=taut, points_bin=bin)

  return(results_table)
}

f.rq_po_plot_boot<- function(rqMod) {
     # rqMod=rq_Mod_sel[[1]][[1]]

    rq_mod_tmp <- rqMod$rqlim_t %>%
    mutate(taustx=as.character(tau)) %>%
    drop_na(RqLim) %>% 
    nest(data=-taustx) %>%
    mutate(data=map(data, ~{.x %>%
                    arrange(.,RqLim) %>%
                    mutate(.,points_bin = ntile(RqLim, n=10)) })) %>% 
    unnest(data) %>% 
    select(y, RqLim, tau, points_bin) %>% 
    group_by(tau, points_bin)
  
  as.list(expand_grid(
    taut=unique(rq_mod_tmp$tau),
    bins=unique(rq_mod_tmp$points_bin))) %>% 
    transpose %>%
    map(~f.bootstrap_table(data = rq_mod_tmp, 
                         taut = .x$taut, bin = .x$bins, 
                         n_bootstraps = 1000)) %>% 
    bind_rows() %>%
    mutate(tau=as.factor(tau)) %>% 
    f.rq_po_plot(.,metas=rqMod$meta)
}


#| label: rq_Mod_plot_res ----

f.pl_rq_Mod_res<-function(rqMod,titleG){
    # rqMod=rq_Mod_sel[[1]][[1]]
  rqResid <- melt(residuals(rqMod$modelq)) %>% 
    rename(tau=Var2, residuals=value) 
  subtitleG <- sprintf("%s vs %s",rqMod$meta$reponset,rqMod$meta$predictt)
  ggplot(rqResid, aes(sample=residuals, color=tau)) + 
    stat_qq() + stat_qq_line() + 
    facet_wrap(~tau) +
    f.Scalc_rq() +
    theme(legend.position="none") +
    labs(title=paste("Residuals for ",titleG,sep=""),
         subtitle=subtitleG)
}

#| label: rq_Mod_plot_1d ----

f.pl_rq_Mod_1d<-function(rqMod,titleG){
    # rqMod=rq_Mod_all[[1]][[1]]

  ggplot() +
    geom_line(data=rqMod$mod_pred_t, 
              aes(x=x1, y=RqMod, col=as.factor(tau))) +
    geom_point(data=rqMod$rqlim_t %>% 
                 filter(tau==max(.$tau,na.rm=TRUE)), 
               aes(x=x1, y=y,shape=status),
               alpha=0.5) +
    f.Scalc_rq() +
    scale_shape_manual(values=c(8, 21)) +
    guides(shape = "none") +
    theme(legend.position="bottom") +
    labs(title=titleG,
       color="Quantile", 
         x=rqMod$meta$predictt, y=rqMod$meta$reponset)
}


#| label: rq_Mod_plot_2d ----
# that could be nice to have the density plots in margin to see the max etc ?

f.pl_rq_Mod_2d<-function(rqMod,titleG){
    # biolo=reponse_l[[1]]
    # rqMod=rq_Mod_sel[[biolo$rdescr]][[1]]
  
  pl_2d<-list(all=NULL, one=NULL)

  pl_2d$all<-ggplot(rqMod$mod_pred_t) +
    geom_raster(aes(x = x1, y = x2, fill = RqMod),alpha=0.8) + 
    # geom_point(data=rqMod$rqlim_t %>% filter(status=="under"),aes(x=x1, y=x2), 
    #            shape=21, size=1, color="gray50", #colInliers, #fill=colInliers, 
    #            alpha=.4, na.rm = TRUE) +
    geom_density_2d(data=rqMod$rqlim_t %>% filter(status=="under"),
                    aes(x=x1, y=x2),
                    alpha=0.4,color="gray50",bins = 30) +
    geom_point(data=rqMod$rqlim_t %>% filter(status=="over"),aes(x=x1, y=x2), 
               shape=8, size=1, color=colOutliers, #fill=colOutliers,
               alpha=.2, na.rm = TRUE) +
    labs(title=titleG,
         x=rqMod$meta$predictl[1],
         y=rqMod$meta$predictl[2],
         fill = paste("SDM-NEO\n",rqMod$meta$unit,sep="")) +
    guides(alpha = "none") +
    f.Scalf_rq2d() +
    facet_wrap(taust ~ .) +
    theme()
  
    # Define surfaces for each tau
    pl_2d$one<- map(taus_l, ~{ #taut=taus_l[[1]]
      rqMod$mod_pred_t %>% 
        filter(tau==.x) %>% 
      ggplot() +
        geom_raster(aes(x = x1, y = x2, fill = RqMod),alpha=0.8) + 
        # geom_point(data=rqMod$rqlim_t %>% filter(status=="under"),
        #            aes(x=x1, y=x2), 
        #            shape=21, size=1, color="gray50", #colInliers, #fill=colInliers, 
        #            alpha=.4, na.rm = TRUE) +
        geom_density_2d(data=rqMod$rqlim_t %>% filter(status=="under"),
                      aes(x=x1, y=x2),
                      alpha=0.4,color="gray50",bins = 30) +
        geom_point(data=rqMod$rqlim_t %>% filter(status=="over"),
                   aes(x=x1, y=x2), 
                   shape=8, size=1, color=colOutliers, #fill=colOutliers,
                   alpha=.2, na.rm = TRUE) +
        labs(title=sprintf("tau = %0.3f",.x),
             x=rqMod$meta$predictl[1],
             y=rqMod$meta$predictl[2],
             fill = paste("SDM-NEO\n",rqMod$meta$unit,sep="")) +
        guides(alpha = "none") +
        f.Scalf_rq2d()
      
    } ) # map taus
  return(pl_2d)
}


#| label: rq_Mod_plot_3d ----

f.pl_rq_Mod_3d<-function(rqMod,taus_l,titleG){
    # biolo=reponse_l[[1]]
    # rqMod=rq_Mod_sel[[1]][[1]]
  
  yl<-rqMod$meta$reponset
  xl<-rqMod$meta$predictl
  myscene<-list(camera = list(eye = list(x = -1.5, y = 1.5, z = 0.3)),
                aspectmode='cube', # define standard layout scene
                xaxis = list(title = xl[1]), 
                yaxis = list(title = xl[2]), 
                zaxis = list(title = yl))
  pl_3d<-list(all=NULL, one=NULL)
  
  rqlim <- rqMod$rqlim_t %>% 
    filter(tau==taus_l[[length(taus_l)]])
  rqsup<-rqlim %>% filter(status=="over")
  rqinf<-rqlim %>% filter(status=="under")
  
  # Experimental points on 3D graphic and all surfaces
  pl_3d$all<-plot_ly(showlegend=FALSE) %>% 
      add_trace(x = rqsup$x1, y = rqsup$x2, z = rqsup$y, 
                mode = "markers", type = "scatter3d",
                marker = list(size = 3, color = colOutliers, 
                              opacity = 0.7, symbol = "diamond"))%>%
      add_trace(x = rqinf$x1, y = rqinf$x2, z = rqinf$y, 
                mode = "markers", type = "scatter3d",
                marker = list(size = 3, color = colInliers, 
                              opacity = 0.7, symbol = "circle")) %>%
    layout(title = titleG, scene = myscene)
  for (t in 1:length(taus_l)){ # t=4
      gridxmat<-rqMod$mod_grid$gridxmat
      RqMod_mat <- rqMod$mod_pred_t %>% 
        filter(tau==taus_l[[t]])
      RqMod_mat <- array(RqMod_mat$RqMod,dim=c(graphfine,graphfine))
      pl_3d$all<-pl_3d$all %>% 
          add_surface(x = gridxmat$x, y = gridxmat$y, z = RqMod_mat,
                      opacity = 0.8, colorscale = list(c(0,1),c(blank,colRQ[t])))
  } # taus
  unname(taus_l)
  # Define surfaces for each tau
  pl_3d$one<-map2(taus_l,seq(1:length(taus_l)), ~{ #.x=taus_l[[1]]
    gridxmat<-rqMod$mod_grid$gridxmat
    RqMod_mat <- rqMod$mod_pred_t %>% 
      filter(tau==.x)
    RqMod_mat <- array(RqMod_mat$RqMod,dim=c(graphfine,graphfine))
    rqlim <- rqMod$rqlim_t %>% 
      filter(tau==.x)
    rqsup<-rqlim %>% filter(status=="over")
    rqinf<-rqlim %>% filter(status=="under")
      
    plot_ly(showlegend=F) %>% #, scene=paste("scene",.y,sep="")
      add_surface(x = gridxmat$x, 
                  y = gridxmat$y, 
                  z = RqMod_mat,
                  opacity = 0.9, colorscale = list(c(0,1),c(blank,colRQ[.y])),
                  colorbar=list(title=list(
                      text=paste("SDM-NEO\n",
                                 rqMod$meta$unit,"\nTau=",.x,sep="")))) %>%
      add_trace(x = rqsup$x1, y = rqsup$x2, z = rqsup$y, 
                mode = "markers", type = "scatter3d",
                marker = list(size = 3, color = colOutliers, 
                              opacity = 0.7, symbol = "diamond"))%>%
      add_trace(x = rqinf$x1, y = rqinf$x2, z = rqinf$y, 
                mode = "markers", type = "scatter3d",
                marker = list(size = 3, color = colInliers, 
                              opacity = 0.7, symbol = "circle")) %>%
      layout(title = titleG, scene = myscene)
    } ) # map taus
  return(pl_3d)
}


#| label: rq_Mod_mars_stat ----

f.st_rq_Mod_map<-function(rqmodmars,rqMod,tauchoice){
    # rqmodmars=rq_Mod_mars_sel[[1]][[1]]
    # rqMod=rq_Mod_sel[[1]][[1]]
  vars<-unlist(strsplit(rqMod$meta$predict,rqMod$meta$type,fixed = TRUE))
  varsl<-rqMod$meta$predictl

  rqmodmars <- rqmodmars %>% 
    st_drop_geometry() %>% 
    filter(grepl("Mudflat",.$Zone)) %>%
    # filter(!Zone %in% c("Bay","Cote Fleurie","Ilot Oiseaux","Octeville")) %>% 
    mutate(across(any_of(c("bathy","tenfon_mxd")), ~f.replace_outliers(.,1.5))) %>%
    mutate(Zone=fct_drop(.$Zone)) %>% 
    mutate(across(c(Zone,Period),factor))
  
  biorange<-max(
    max(f.replace_outliers(rqmodmars %>% select(sprintf("t%.3f",tauchoice))), na.rm=TRUE),
    max(f.replace_outliers(rqMod$rqlim_t$y), na.rm=TRUE)) %>% 
    ceiling(.)
  
  descr<-list(pred=NULL, mod=NULL, res=NULL, all=NULL)

  # pred and mod could be done with a facet_grid, but how to put the dots only on the model?
  
  descr$pred <-map2(as.list(vars), as.list(varsl), ~{
    rqmodmars %>% 
    ggplot() + 
    geom_density_ridges(aes(y = Period , fill=Period,
                            height = after_stat(count), # to have a comparative height between periods
                            x=.data[[.x]]), 
                        stat="density", alpha=0.3, size = .2) +
    facet_grid(. ~ Zone, 
               # scales = "free_y",
               labeller = label_wrap_gen(width = 20)) + 
    labs(x=.y, y="") +
    f.Scale_brew() } )

  tmp<-rqMod$rqlim_t %>%
    filter(grepl("Mudflat",.$Zone)) %>%
    mutate(Zone=fct_drop(.$Zone)) %>%
    mutate(across(c(Zone,Period),factor)) %>% 
    group_by(Zone, Period) %>% 
    summarise(Observed = quantile(y, unique(tauchoice), na.rm=TRUE),
              n=n())

  descr$mod <-rqmodmars %>% 
    ggplot() + 
    geom_density_ridges(aes(y = Period , fill=Period, 
                            height = after_stat(count), # to have a comparative height between periods
                            x=.data[[sprintf("t%.3f",tauchoice)]]), 
                        stat="density", alpha=0.25, size = .2) +
    geom_point(data=tmp,aes(y=Period,x=Observed, fill=Period, size=n),
               shape=21, alpha=0.7) +
    xlim(c(0,biorange)) +
    facet_grid(. ~ Zone, 
               # scales = "free_y",
               labeller = label_wrap_gen(width = 20)) +
    guides(size="none") + 
    scale_size_continuous(range = c(1, 3)) +
    labs(y="",x=sprintf("%s (%s) t%.3f - %s", 
                          rqMod$meta$typetxt, 
                          rqMod$meta$unit,tauchoice,
                          rqMod$meta$predictt)) +
    f.Scale_brew()
  
  # descr$res <- rqMod$rqlim_t %>%
  #   filter(grepl("Mudflat",.$Zone)) %>%
  #   mutate(Zone=fct_drop(.$Zone)) %>% 
  #   mutate(across(c(Zone,Period),factor)) %>% 
  #   ggplot() +
  #   geom_density_ridges(aes(y = Period , fill=Period,
  #                         x=y), alpha=0.3, size = .2) +
  #   xlim(c(0,biorange)) +
  #   facet_grid(. ~ Zone,
  #              # scales = "free_y",
  #              labeller = label_wrap_gen(width = 20)) +
  #   labs(y="",x=rqMod$meta$reponset) +
  #   f.Scale_brew()
 
 
  descr$all<-wrap_plots(descr$pred, ncol=1) / descr$mod & #/ descr$res &
    # plot_layout(guides = "collect") & 
    theme(legend.position="none") &
    plot_annotation(tag_levels = c('A')) 
  return(descr)
}


#| label: rq_Mod_plot_map ----

f.pl_rq_Mod_map<-function(rqmodmars,rqMod,tauchoice){
   # rqmodmars=rq_Mod_mars_sel[[1]][[1]]
   # rqMod=rq_Mod_sel[[1]][[1]]

  pl_map<-list(all=NULL, one=NULL)

  vars<-unlist(strsplit(rqMod$meta$predict,rqMod$meta$type,fixed = TRUE))
  varsl<-rqMod$meta$predictl

  # # option ggplot
  # sdmname<-rqMod$sdmname
  # tmp_sf <-
    # ggplot(data = rqBmars) +
  #   geom_sf(aes(fill = SDM), color = NA) +
  #   facet_wrap(~Annee) +
  #   labs(x="Latitude",y="Longitude",
  #        fill = sprintf("SDM-NEO\n%s",sdmname)) +
  #   theme(plot.margin = margin(0.05,0.05,0.05,0.05, "cm")) +
  #   scale_fill_distiller(palette = "Spectral")
  
  pl_map$all <- tm_shape(osm_df) +  
    tm_rgb(saturation = .6, alpha= 0.8) +
  # tm_scale_bar(position = c("left", "bottom"), width = 0.15) + #SCALE
  # tm_compass(position = c("left", "top"), size = 2) +          #NORTH COMPASS
    tm_shape(rqmodmars) +
      tm_fill(col = sprintf("t%.3f",tauchoice),
              style="cont", palette = "-Spectral", alpha = 0.7) + # pal_csps(9)
    tm_facets(by="Period", nrow = 1) +
    tm_shape(tm_areas) +  
      tm_borders(alpha = 0.3, lwd = 1) +
    tm_layout(
      main.title=sprintf("%s (%s) t%.3f - %s", 
                      rqMod$meta$typetxt, 
                      rqMod$meta$unit,tauchoice,
                      rqMod$meta$predictt),
      main.title.size= .8,
      panel.label.bg.color = 'white',
      panel.label.size=.8,
      legend.outside = TRUE, 
      legend.outside.size = .08,
      legend.title.size= .8,
      legend.text.size= .6,
      legend.format=list(text.separator="-"),
      outer.margins=c(0,0,0,0), inner.margins=c(0,0,0,0))
  
  peri<-setNames(as.list(levels(rqmodmars$Period)), levels(rqmodmars$Period))
  # Define map for each Period
  pl_map$one <- map(peri, ~{
    tm_shape(osm_df) +  
      tm_rgb(saturation = .6, alpha= 0.8) +
    tm_shape(rqmodmars %>% 
               filter(Period %in% .x)) +
      tm_fill(col = sprintf("t%.3f",tauchoice), 
              title=.x, 
              style="cont", palette = "-Spectral", alpha = 0.7) +
    tm_shape(tm_areas) +  
      tm_borders(alpha = 0.3, lwd = 1) +
    tm_layout(
              main.title=sprintf("%s (%s) t%.3f - %s", 
                              rqMod$meta$typetxt, 
                              rqMod$meta$unit,tauchoice,
                              rqMod$meta$predictt),
              main.title.size= .8,
              panel.label.size=.8,
              legend.outside = TRUE,
              legend.outside.size = .08,
              legend.title.size= .8,
              legend.text.size= .6,
              legend.format=list(text.separator="-"),
              outer.margins=c(0,0,0,0), inner.margins=c(0,0,0,0))
  } ) # map peri
  return(pl_map)
}


#| label: rq_Mod_plot_si ----

f.pl_suit_idx<-function(rqmodmars,rqMod){
#  NB: https://cran.r-project.org/web/packages/Hmisc/Hmisc.pdf for the smean.cl.boot calculation mode
    # rqmodmars=nlrqdata$rq_Mod_mars_sel[[1]][[1]]
    # rqMod=nlrqdata$rq_Mod_sel[[1]][[1]]
  # rqmodmars=rq_Mod_mars_sel[[1]][[1]]
  # rqMod=rq_Mod_sel[[1]][[1]]
  
  pl_tmp<- rqmodmars %>% 
    st_drop_geometry() %>% 
    group_by(Zone, Period, #Annee, 
             NINJ) %>% 
    # mutate(Annee=as.factor(Annee)) %>% 
    summarise(across(starts_with("si"),~mean(.x,na.rm=TRUE))) %>% 
    pivot_longer(cols = starts_with("si"),
                 names_to = "tau",
                 names_prefix = "si_",
                 values_to = "suit_index") %>% 
    filter(!Zone %in% c("Ilot Oiseaux")) %>% #,"Cote Fleurie","Octeville"
    filter(tau==sprintf("t%.3f",tauchoice)) %>% 
    ggplot() +
    # stat_summary(aes(x=Period, group= Zone, fill=Zone,
    #                  y=suit_index),
    #              fun.data = mean_se, geom = "errorbar",
    #              position=position_dodge(width = .3),
    #              width=.3,
    #              linewidth = .4) +
    stat_summary(aes(x=Period, group= Zone, col=Zone, 
                     y=suit_index),  
                 fun.data = mean_cl_boot, geom = "line",
                 position=position_dodge(width = .3),
                 linewidth = 1) +
    stat_summary(aes(x=Period, group= Zone, fill=Zone,
                     y=suit_index),
                 fun.data = mean_cl_boot, geom = "pointrange",
                 position=position_dodge(width = .3),
                 shape=21, size=.6, stroke = .3,
                 linewidth = .4) +
    labs(title=sprintf("%s", #%s (%s) t%.3f
                       # rqMod$meta$typetxt, 
                       # rqMod$meta$unit,tauchoice,
                       rqMod$meta$predictt),
         color="Areas", fill="Areas") +
    ylab("Suitability Index") +
    theme(legend.position="bottom",
          text=element_text(size=10)) +
    f.Scale_brew() ; pl_tmp

  # ggsave(pl_tmp,
  #        width = 8, height = 8, dpi=400,
  #        filename = sprintf("%sfig-suit_idx_A.tiff",wdgraph) )
}

#| label: rq_Mod_plot_map_gif ----

f.pl_rq_Mod_map_gif<-function(rqmodmars,sdmname,tauchoice){
    # rqmodmars=rq_Mod_mars_sel[[1]][[1]]
    # sdmname=names(rq_Mod_mars_sel[[1]])[1]
  
  tmap_mode("view") # "plot" "view"
  tmp_sf <-
  tm_basemap(leaflet::providers$OpenStreetMap.HOT) +
  tm_shape(rqmodmars) +
    tm_fill(col = sprintf("t%.3f",tauchoice),
            palette = "-Spectral", alpha = 0.7) +
    # tm_borders("white", lwd = 0) +
    tm_facets(along = "Period") +
  tm_layout(legend.outside = TRUE,
            title.size=10)

  tmap_animation(
    tmp_sf, filename = sprintf("%s/%s/%s_%s_map.gif",
                    graph_path, substr(sdmname,1,6), espece, sdmname),
    fps = 1, width = 2000, height = 2000, dpi=600 )
  # tmap_animation(
  #   tmp_sf, filename = sprintf("%s/%s_map.gif",
  #                              wdgraph, vars$Var),
  #   fps = 1, width = 2000, height = 2000, dpi=600 )
  
}


#| label: auto_tabset

# ::: panel-tabset
# ```{r}
# #| include: true
# #| results: asis
# #| fig-cap: "Scenarios computed with linear with interaction (top, numbered 1) and nonlinear with gaussian equation (bottom numbered 2), the isometric curve represents the observed HMS data"
# #| fig-width: 16
# #| fig-height: 8
# 
# iwalk(pl_board_rq_nlrq_real, 
#       ~ {
#         cat('#### ', .y, ' {#sec-sdm_real}\n\n')
#         print(.x)
#         cat('\n\n')
#       })
# ```
# :::

# iwalk(nlrqdata$pl_rq_sel_3d,
#       ~ {
#         cat('#### ', .y, '\n\n')
#         map(.x, ~print(.x$all))
#         cat('\n\n')
#         })

# imap(nlrqdata$pl_rq_sel_3d,
#       ~ {
#         knitr::knit_child(text = c(
#               "#### ", .y, "\n\n",
#               "```{r}", "\n\n",
#               "#| layout-ncol: 2",
#               "#| label: fig-stat descr ", .y, "\n\n",
#               "#| fig-cap: try about this ", .y, "\n\n",
#               "map(.x, ~print(.x$all))",
#               "\n\n", "```", "\n\n"
#             ), envir = globalenv(), quiet = TRUE)
#         })
