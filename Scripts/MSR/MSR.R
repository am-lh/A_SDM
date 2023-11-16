# ____________________________________________________________________________________________________
# Mass specific respiration rate (MSR) of aquatic invertebrates
# Compiled by Amelie Lehuen, 2020
# Based on excel RespirationANN01.xlsx (Brey, 2010) :
#   http://www.thomas-brey.de/science/virtualhandbook/navlog/index.html
# # Input variables : CANNOT BE VECTORS MUST BE SINGLE VALUES
  # dens : density in ind/m2
  # biom : Biomass (gAFDW/m2)
  # aphia : AphiaID of species in Repertory
  # temp : temperature in degC, converted in K in function
  # depth : water height in m, by default 1m for intertidal
  # starved : is the species has been fed or not, but defaut yes, can be omitted
# # Inner data
  # Energy density <- 21.4469 (J/mgAFDW) (Brey et al., 2010 Body composition in aquatic organisms - A global data bank of relationships between mass, elemental composition and energy content)
  # Mj <- biom/dens : individual mass energy (J)
# # Output variables :
  # MSR.moy  MSR.lower  MSR.upper 
  # Unit option out : "w" in W/ind or "j" in J/J/d, "j" by default

# Folder MSR MUST be in work directory defined in called script
# Result is a list of result of mean and upper and lower confidence interval
# INPUT MUST BE A SCALAR NOT A VECTOR

# ____________________________________________________________________________________________________

msr <- function(dens,biom,aphia,tempe,depth=1,unitout,starved) {
  EnerDens<-21.4469*10^3; # Energy density in J/gAFDW
  Mj=0
  if (missing(starved)) {starved=0}
  if (missing(unitout)) {unitout="j"}
  if (length(aphia)>1){ # test if input is only one species (not yet input as a vector)
    error('Input argument aphia must be scalar')
    MSR.moy <- NA; MSR.upper <- NA; MSR.lower <- NA
  } else if (dens==0 || biom==0) {  
    MSR.moy<-0; MSR.lower<-0; MSR.upper<-0 
  } else if (aphia==0) {  
    MSR.moy<-0; MSR.lower<-0; MSR.upper<-0 
  } else {
    load ("BDD_Brey.RData")
    squish <- function(x) {1/(1+exp(-x))}
    Mj=biom/dens*EnerDens
    
    # Creation of the input vector with parameters corresponding to specie
    esp <- Repertory[Repertory$AphiaIdAccepted == aphia,]
    if (nrow(esp)==0) {
      warning('AphiaID not founded in Repertory')
      MSR.moy <- NA; MSR.upper <- NA; MSR.lower <- NA
    } else {   
      if (nrow(esp)>1) {
        warning('AphiaID founded more than once in Repertory, first one got selected')
        esp <- esp[1,]}
        esp_input <- rep(1,29) # The first factor for intercept
        esp_input[2:4] <- c(1/(273+tempe),log10(depth),log10(Mj))
        esp_input[ann_factors==esp$TaxaBrey] <- -1
        esp_input[ann_factors==esp$Mobility] <- -1
        esp_input[ann_factors==esp$Trophic] <- -1
        if (esp$Vision!=0 & !is.na(esp$Vision)) {esp_input[ann_factors=="Vision"] <- -1}
        if (starved==1) {esp_input[ann_factors=="Starved"] <- -1}
      
      # Calculation of the 4 nods for each 5 ANN in the third dimension
      # H1 = Squish(b0 + b1 * 1/T +  b2 * log(D) + b3 * log(M) + b4 *Porifera ... + b28 * Starved)
      h = squish(apply(esp_input*ann, MARGIN=c(2,3), sum));
    
      # log(MSR) = (a0 + a1 * H1 + a2 * H2 + a3 * H3 + a4 * H4) * a5 + a6	
      logMSR <- (factOut[1,,]+
                   factOut[2,,]*h[1,]+
                   factOut[3,,]*h[2,]+
                   factOut[4,,]*h[3,]+
                   factOut[5,,]*h[4,])*
                factOut[6,,]+factOut[7,,]
    
      # Respiration Rate : Mean MSR	J/J/d  95% Confidence Limits	(lower)	(upper)
      MSR.moy<-10^(mean(logMSR))
      MSR.lower<-10^(mean(logMSR)-sd(logMSR)*0.953)
      MSR.upper<-10^(mean(logMSR)+sd(logMSR)*0.953)
    } 
    }
  # Conversion in W/ind
  if (unitout=="w"){
    MSR.moy <- MSR.moy*Mj/(24*60*60) # J/J/d * BodyMass (J/ind) / (days->seconds) => J/s/ind => W/ind
    MSR.lower <- MSR.lower*Mj/(24*60*60)
    MSR.upper <- MSR.upper*Mj/(24*60*60)
  }
  return (MSR.moy) #(MSR=list(MSR.moy,MSR.lower,MSR.upper))
}# end of function

# # Update Brey table ----
# library(readxl)
# source <- "Fauna_Brey.xlsx" # Main file with coefficients and fauna data from Worms completed
# Repertory <- as.data.frame(read_excel(source,sheet = "Repertory"))
# factbrey <- as.data.frame(read_excel(source,sheet = "Brey"))
# # Factors names extraction for each node
# ann_factors = factbrey$Factor[factbrey[,1]=="H1"]
# # Factor extraction to compile the 4 nodes for each 5 ANN in 3rd dim
#     factOut <- cbind(factbrey$ANN4_1[factbrey[,1]=="Output"],
#                       factbrey$ANN4_2[factbrey[,1]=="Output"],
#                       factbrey$ANN4_3[factbrey[,1]=="Output"],
#                       factbrey$ANN4_4[factbrey[,1]=="Output"],
#                       factbrey$ANN4_5[factbrey[,1]=="Output"]);
#     dim(factOut) <- c(length(factbrey$ANN4_1[factbrey[,1]=="Output"]), 1, 5)
# # Factor extraction of the 4 nodes for each 5 ANN in 3rd dim
# ann = cbind(cbind(factbrey$ANN4_1[factbrey[,1]=="H1"],factbrey$ANN4_1[factbrey[,1]=="H2"],
#                   factbrey$ANN4_1[factbrey[,1]=="H3"],factbrey$ANN4_1[factbrey[,1]=="H4"]),
#             cbind(factbrey$ANN4_2[factbrey[,1]=="H1"],factbrey$ANN4_2[factbrey[,1]=="H2"],
#                   factbrey$ANN4_2[factbrey[,1]=="H3"],factbrey$ANN4_2[factbrey[,1]=="H4"]),
#             cbind(factbrey$ANN4_3[factbrey[,1]=="H1"],factbrey$ANN4_3[factbrey[,1]=="H2"],
#                   factbrey$ANN4_3[factbrey[,1]=="H3"],factbrey$ANN4_3[factbrey[,1]=="H4"]),
#             cbind(factbrey$ANN4_4[factbrey[,1]=="H1"],factbrey$ANN4_4[factbrey[,1]=="H2"],
#                   factbrey$ANN4_4[factbrey[,1]=="H3"],factbrey$ANN4_4[factbrey[,1]=="H4"]),
#             cbind(factbrey$ANN4_5[factbrey[,1]=="H1"],factbrey$ANN4_5[factbrey[,1]=="H2"],
#                   factbrey$ANN4_5[factbrey[,1]=="H3"],factbrey$ANN4_5[factbrey[,1]=="H4"]))
# dim(ann) <- c(length(factbrey$ANN4_1[factbrey[,1]=="H1"]),4,5)
# save(Repertory,factbrey,ann_factors,ann,factOut,file = "BDD_Brey.RData")

# # To test and use
# source("MSR.R")
# dens <- 63.694267515923560
# biom <- 39.238681844261215
# aphia <- 138998
# tempe <- 18
# depth <- 1
# myMSR = msr(dens,biom,aphia,tempe,depth,"j")
# # RespirationANN01.xlsx results MSR.moy = 0,001323782  MSR.lower = 0,001199948  MSR.upper =  0,001460396

# dens <- 0
# biom <- 0
# aphia <- 138998
# tempe <- 18
# depth <- 1
# myMSR = msr(dens,biom,aphia,tempe,depth,"j")