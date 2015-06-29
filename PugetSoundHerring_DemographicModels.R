
## ============================================================================
# MARSS models for testing population structure within Puget Sound Pacific herring
# Megsie Siple
# Last edit: 6/29/15
# This file contains all the models we compared to find out which population structure
# best describes the biomass dynamics for Puget Sound Pacific herring. 
# We used AICbp values to compare model fits. Model descriptions are in Table 1 of the
# manuscript, and detailed desciptions of all the matrices are in Electronic Supplementary
# Material, Appendix 3.
## ============================================================================

wd <- "/Users/mcsiple/Dropbox/Herring Data_3/Demographic models"
setwd(wd)

require(MARSS)
require(reshape2)

#Read in data
load("Eggs.RData")      # Eggs.RData and tEggs.RData are log-transformed biomass 
load("tEggs.RData")     # Egg biomass data. Sites are in rows, years are columns

#Remove sites outside Puget Sound
eggs <- eggs[,-c(19,21)]
teggs <- teggs[-c(18,20),]

# -------------------------------------------------------------------------
# Basic plot of lob(biomass) over time ------------------------------------
# -------------------------------------------------------------------------

# Make a plot of the 21 time series
stockID <- colnames(eggs)[2:ncol(eggs)]
NN<-1:length(stockID)
par(mfrow=c(1,1),mai=c(rep(0.5,4)))

plot(eggs[,1],eggs[,3],
     xlim=c(1972,2024),ylim=c(0,max(eggs$CP)),  
     xlab="Year",ylab="log(Biomass)",main="",type="n",xaxt="n")

      for (i in 2:length(stockID)){
        points(eggs[,1],eggs[,i],col=rainbow(23)[i],pch=16)
        lines(eggs[,1],eggs[,i],col=rainbow(23)[i],pch=16)
      }
      axis(1, at=seq(1973,2012,by=6),labels=seq(1973,2012,by=6),tck=-.01)
      legend('topright',inset=0.01, legend=stockID, col=rainbow(23), pch=16)


###########################################################################
#SECTION 1:  Demographic pop structure
###########################################################################

# -------------------------------------------------------------------------
# Model 1: Pops have same trend and the same process variance -------------
# -------------------------------------------------------------------------

model.h1$Q="diagonal and equal"   
model.h1$R="diagonal and equal"
model.h1$Z= "identity"      #This means Z is a n x n identity matrix and m = n.
model.h1$A="scaling"        
model.h1$U="equal"

kem.h1 = MARSS(teggs, model=model.h1,control=list(maxit=2000,allow.degen=TRUE))

# kemh1.with.AICb = MARSSaic(kem.h3, output = "AICbp",
#                            Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 2: Stocks covary with same growth rate (m=n=19) -------------------
# -------------------------------------------------------------------------

model.h2=list() #h9 is based on model.h2
model.h2$R="diagonal and equal"
model.h2$A="scaling"     
model.h2$Q="equalvarcov"      
model.h2$Z="identity" 
model.h2$U="equal"

kem.h2 = MARSS(teggs, model=model.h2,control=list(maxit=2000,allow.degen=TRUE))

#kemh9.with.AICb = MARSSaic(kem.h9, output = "AICbp",
#                           Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 3: Pops separated by NOAA-defined basins (m=5)---------------------
# -------------------------------------------------------------------------

sites <- rownames(teggs)
NOAA_basins <- 0

for (i in 1:length(sites)){
  if (sites[i]=="FIDBAY"
      |sites[i]=="SAMPORT" | sites[i]=="INTSJ"
      |sites[i]=="NWSJ"|sites[i]=="SEMI"
      |sites[i]=="CP" |sites[i]=="DISCO" 
      |sites[i]=="DUNG") NOAA_basins[i]="NPS" 
  
  if (sites[i]=="SKAGIT" | sites[i]=="PS"
      |sites[i]=="HH") NOAA_basins[i]="WHIDBEY"
  
  if (sites[i]=="QM" | sites[i]=="POPM"
      |sites[i]=="KH") NOAA_basins[i]="MAIN"
  
  if (sites[i]=="SHOOD" | sites[i]=="QUIL"
      |sites[i]=="PG") NOAA_basins[i]="HOODCANAL"
  
  if (sites[i]=="SQUAXIN" | sites[i]=="WOLLO") NOAA_basins[i]="SPS"
  }

model.h3=list()
model.h3$R="diagonal and equal"
model.h3$A="scaling"
model.h3$Z = factor(NOAA_basins)
model.h3$U = matrix(c("uNorth","uWhidbey","uMain","uHood","uSouth"),5,1)       #Data do not have to be in this order
model.h3$Q="diagonal and equal"

kem.h3 = MARSS(teggs, model=model.h3,control=list(allow.degen=TRUE,maxit=10000))

#kemh3.with.AICb = MARSSaic(kem.h3, output = "AICbp",
#                            Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 4: Herring are one panmictic population (m=1)----------------------
# -------------------------------------------------------------------------

model.h4=list()
model.h4$Q="diagonal and equal"
model.h4$R="diagonal and equal"
model.h4$Z=matrix(1,19,1)
model.h4$A="scaling"
model.h4$U="unequal"

kem.h4 = MARSS(teggs, model=model.h4, control=list(minit=100))

#kemh4.with.AICb = MARSSaic(kem.h4, output = "AICbp",
#                           Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 5: Pops separated by genetics -------------------------------------
# -------------------------------------------------------------------------

#Herring pops distinguished by genetic differences (Small et. al 2005) 
# Groups: Cherry Point, Squaxin Pass, and all other

model.h5=list()
model.h5$Q="diagonal and equal"
model.h5$R="diagonal and equal"
model.h5$Z=cbind(c(1,rep(0,times=18)),c(rep(0,times=17),1,0),c(0,rep(1,16),0,1))
model.h5$A="scaling"
model.h5$U=matrix(c("CherryPoint","Squaxin Pass","Others"),3,1)

kem.h5 = MARSS(teggs, model=model.h5)

kemh5.with.AICb = MARSSaic(kem.h5, output = "AICbp",
                           Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 6: Pops separated by regions (m = 3) ------------------------------
# -------------------------------------------------------------------------

sites <- rownames(teggs)
regions <- 0
for (i in 1:length(sites)){
  if (sites[i]=="FIDBAY" | sites[i]=="SAMISH"
      |sites[i]=="SAMPORT" | sites[i]=="INTSJ"
      |sites[i]=="NWSJ"|sites[i]=="SEMI"
      |sites[i]=="CP") regions[i]="NPS"
  if (sites[i]=="SQUAXIN" | sites[i]=="WOLLO"
      |sites[i]=="QM" | sites[i]=="POPM"
      |sites[i]=="SHOOD"|sites[i]=="QUIL"
      |sites[i]=="PG" |sites[i]=="KH"
      |sites[i]=="PS" |sites[i]=="HH"
      |sites[i]=="SKAGIT") regions[i]="SPS"
  if (sites[i]=="DISCO" | sites[i]=="DUNG"
      |sites[i]=="SEQUIM") regions[i]="STRAIT"
}


model.h6$R="diagonal and equal"
model.h6$A="scaling"
model.h6$Z = factor(regions)
model.h6$U = matrix(c("uNorth","uStrait","uSouth"),3,1)       #Data do not have to be in this order 
model.h6$Q="diagonal and equal"

kem.h6 = MARSS(teggs, model=model.h6,control=list(allow.degen=TRUE,maxit=1000))

kemh6.with.AICb = MARSSaic(kem.h6, output = "AICbp",
                             Options = list(nboot = 100, silent=FALSE))

# -------------------------------------------------------------------------
# Model 7: Pops separated by contaminants (West et al 2008) ---------------
# -------------------------------------------------------------------------

# Groups are: Cherry Point, Squaxin Pass, Others

model.h7=list()
model.h7$Q="diagonal and equal"
model.h7$R="diagonal and equal"
model.h7$Z=cbind(c(1,rep(0,times=18)),c(rep(0,times=14),1,rep(0,times=4)),c(0,rep(1,13),0,rep(1,times=4)))
model.h7$A="scaling"
model.h7$U=matrix(c("uCherryPoint","uSemiahmoo","uOthers"),3,1)

kem.h7 = MARSS(teggs, model=model.h7,control=list(allow.degen=TRUE,maxit=1000))

# kemh7.with.AICb = MARSSaic(kem.h7, output = "AICbp",
#                           Options = list(nboot = 100, silent=FALSE))


# -------------------------------------------------------------------------
# Model 8: Pops are separate, with different trends and the same process variance --------
# -------------------------------------------------------------------------


model.h8=list()
model.h8$Q="diagonal and equal"   
model.h8$R="diagonal and equal"
model.h8$Z= "identity"
model.h8$A="scaling"  
model.h8$U="unequal"    

kem.h8 = MARSS(teggs, model=model.h8,control=list(maxit=500,allow.degen=FALSE)) 

#kemh8.with.AICb = MARSSaic(kem.h8, output = "AICbp",
#                          Options = list(nboot = 100, silent=FALSE))



# -------------------------------------------------------------------------
# New section - after reviews from Oecologia, 1/26/15 ---------------------
# Reviewer #1 recommended trying to estimate all the elements of Q, the 
# var-cov matrix. 
# -------------------------------------------------------------------------

# Model 9: each pop has a unique variance, each pair of pops has unique covariance ----
# Add more iterations to estimate Q (maxit=2000)

model.h9=list()
model.h9$Q="unconstrained"   
model.h9$R="diagonal and equal"
model.h9$Z= "identity"      #Z is a n x n identity matrix and m = n.
model.h9$A="scaling"        
model.h9$U="unequal"        

kem.h9 = MARSS(teggs, model=model.h9,control=list(maxit=2000,allow.degen=FALSE)) 

# kemh9.with.AICb = MARSSaic(kem.h9, output = "AICbp",
#                           Options = list(nboot = 100, silent=FALSE))

# Print full Q matrix
(fullQ <- coef(kem.h9,type="matrix")$Q)

# ------------------------------------------------------------------------


