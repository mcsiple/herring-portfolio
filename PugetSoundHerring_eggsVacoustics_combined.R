###################################################
#SECTION 2:  Egg vs. acoustic surveys
###################################################
# To compare different models that describe observation error for two survey types
# by WDFW, we set up 4 different models. We compared the fits of each to the 
# survey data, and then compared the estimated observation errors for each type of survey 
# in the best fit model. Detailed model structure is the Electronic Supplementary Materials.


# Hypotheses:
# 1- Eggs and acoustics measure the same process, obs error has same magnitude at each beach
# 2- Eggs and acoustics measure different processes, obs error has same magnitude at each beach
# 3- Eggs and acoustics measure the same process, obs error varies by beach
# 4- Eggs and acoustics measure different processes, obs error varies by beach

setwd(wd)
library(MARSS)

load("tTotal.RData")
ttotal <- ttotal[-c(18,20,39,41),]  #Take out stocks outside Puget Sound
op <- par(oma=c(5,7,1,1))

# -------------------------------------------------------------------------
# 1 - Two survey types measure the same process; mmt error the same 
# across all beaches (e.g., observation error from egg surveys is the same across all beaches)
# ------------------------------------------------------------------------

model.h1=list() 
model.h1$Q = "equalvarcov"

R <- matrix(list(0),38,38)
for (i in 1:19) R[i,i]=paste("eggsurvey")
for (i in 20:38) R[i,i]=paste("acousticsurvey")
model.h1$R <- R

Z <- factor(c(seq(1,19), seq(1,19)))
model.h1$Z=Z  
model.h1$U="equal"
#B=1 is the default

kem.h1 = MARSS(ttotal, model=model.h1,control=list(safe=TRUE,maxit=3000),silent=FALSE)
kem.h1.with.AICbp = MARSSaic(kem.h1, output = "AICbp",
                               Options = list(nboot = 100, silent=FALSE))


# -------------------------------------------------------------------------
# 2- Eggs and acoustics measure different processes, obs error magnitude same at each beach
# -------------------------------------------------------------------------
model.h2 <- model.h1
R <- matrix(list(0),38,38)
for (i in 1:19) R[i,i]=paste("eggsurvey")
for (i in 20:38) R[i,i]=paste("acousticsurvey")
model.h2$R <- R

#Need to hand-make a Q matrix that has two different diagonals
Qll = Qur <- matrix(list(0),19,19)

Qul <- matrix(list("covar.eggs"),19,19)
diag(Qul) <- rep("var.eggs",times=19)

Qlr <- matrix(list("covar.acoustics"),19,19)
diag(Qlr) <- rep("var.acoustics",times=19)

Q <- cbind(rbind(Qul,Qll),rbind(Qur,Qlr))
model.h2$Q = Q

Z1 = factor(seq(1:38))
model.h2$Z = Z1
model.h2$U = "equal"

kem.h2 = MARSS(ttotal, model=model.h2,control=list(safe=TRUE,maxit=1000))
kem.h2.AICbp = MARSSaic(kem.h2, output = "AICbp",
                          Options = list(nboot = 100, silent=FALSE))



#  ------------------------------------------------------------------------
# 3 - Two survey types measure the same process; mmt error varies by beach
# ------------------------------------------------------------------------

model.h3=list()
model.h3$Q = "equalvarcov"

R <- matrix(list(0),38,38)
for (i in 1:19) R[i,i]=paste("eggsurvey",rownames(ttotal)[i],sep="_")
for (i in 20:38) R[i,i]=paste("acousticsurvey",rownames(ttotal)[i],sep="_")
model.h3$R <- R
Z <- factor(c(seq(1,19), seq(1,19)))
model.h3$Z=Z  
model.h3$U="equal"
#B=1 is the default
#Let MARSS decide what A will be; surveys may be covering a different portion of the state

kem.h3 = MARSS(ttotal, model=model.h3,control=list(safe=TRUE,maxit=1000))


kem.h3.with.AICbp = MARSSaic(kem.h3, output = "AICbp",
                             Options = list(nboot = 100, silent=FALSE))


# -------------------------------------------------------------------------
# 4 - Eggs and acoustics measure different processes, obs error varies by beach
# ------------------------------------------------------------------------
model.h4=model.h1

Z1 = factor(seq(1:38))
model.h4$Z = Z1
model.h4$U = "equal"   

#Use special Q matrix from H2:
model.h4$Q = model.h2$Q

R <- matrix(list(0),38,38)
for (i in 1:19) R[i,i]=paste("eggsurvey",rownames(ttotal)[i],sep="_")
for (i in 20:38) R[i,i]=paste("acousticsurvey",rownames(ttotal)[i],sep="_")
model.h4$R <- R


kem.h4 = MARSS(ttotal, model=model.h4,control=list(safe=TRUE,maxit=500))
kem.h4.with.AICbp = MARSSaic(kem.h4, output = "AICbp",
                             Options = list(nboot = 100, silent=FALSE))


