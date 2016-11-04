source('./MSmixture/R/simulate_MSmixture.R')

#call simulate_MSmixture to simulate multistate mixture capture-recapture data
set.seed(12345)
t_steps=10
n_strata=3
n_mixtures=2
N_releases = array(100,dim=c(n_mixtures,n_strata,t_steps))
Psi=vector("list",n_mixtures)
Psi[[1]]=matrix(0.1,3,3)
diag(Psi[[1]])=0.8
Psi[[2]]=matrix(1/3,3,3)
Phi=matrix(0.8,t_steps,n_mixtures)
Phi[,2]=0.9
P = matrix(0.6,t_steps,n_strata)
P[,2]=0.7
P[,3]=0.8
enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

Hists = simulate_MSmixture(t_steps=t_steps,N_releases=N_releases,Psi=Psi,Phi=Phi,P=P,enc_hist_format=enc_hist_format,fname=fname)
Hists$record = 1  #try to decrease dimensionality of dp call while retaining data frame
library(marked)
dp=process.data(Hists,model="mvmscjs",strata.labels=list(area=c("A","B","C"),pop=c("1","2","u")))
ddl=make.design.data(dp)
# fix delta = 1 for unknown mixture (true mixture never known)
ddl$delta$fix = ifelse(ddl$delta$obs.pop=="u", 1, 0)
# Set Psi to 0 for cases which are not possible - moving between population subgroups
ddl$Psi$fix[as.character(ddl$Psi$pop)=="1"&as.character(ddl$Psi$topop)=="2"]=0
ddl$Psi$fix[as.character(ddl$Psi$pop)=="2"&as.character(ddl$Psi$topop)=="1"]=0
# Create indicator variables for transitioning between states
ddl$Psi$AtoB=ifelse(ddl$Psi$area=="A"&ddl$Psi$toarea=="B",1,0)  # A to B movement
ddl$Psi$BtoA=ifelse(ddl$Psi$area=="B"&ddl$Psi$toarea=="A",1,0)  # B to A movement
ddl$Psi$AtoC=ifelse(ddl$Psi$area=="A"&ddl$Psi$toarea=="C",1,0)  # A to C movement
ddl$Psi$CtoA=ifelse(ddl$Psi$area=="C"&ddl$Psi$toarea=="A",1,0)  # C to A movement
ddl$Psi$BtoC=ifelse(ddl$Psi$area=="B"&ddl$Psi$toarea=="C",1,0)  # B to C movement
ddl$Psi$CtoB=ifelse(ddl$Psi$area=="C"&ddl$Psi$toarea=="B",1,0)  # C to B movement
# formulas
#Psi.1=list(formula=~-1+AtoB+AtoC+BtoC+BtoA+CtoA+CtoB)
#Psi.1=list(formula=~-1+AtoB+ AtoB:pop + AtoC + AtoC:pop + BtoC + BtoC:pop + BtoA + BtoA:pop + CtoA + CtoA:pop + CtoB + CtoB:pop)
Psi.1=list(formula=~-1+pop:stratum:tostratum)
p.1=list(formula=~area)
delta.1=list(formula= ~0)
Phi.1=list(formula=~pop)
#Phi.1=list(formula=~1)
pi.1 = list(formula = ~0+stratum*pop)

#initial= list(Phi=1.35,p=c(-1,0,1),Psi=c(-1,-1,1,3,3,1))
initial=NULL
mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)

#mod_admb=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)

#pars on real scale
cat(paste("pA = ",expit(mod_optim$results$beta$p["(Intercept)"]),"\n"))
cat(paste("pB = ",expit(mod_optim$results$beta$p["(Intercept)"]+mod_optim$results$beta$p["areaB"]),"\n"))
cat(paste("pC = ",expit(mod_optim$results$beta$p["(Intercept)"]+mod_optim$results$beta$p["areaC"]),"\n"))
cat(paste("S1 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]),"\n"))
cat(paste("S2 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]+mod_optim$results$beta$Phi["pop2"]),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop1:stratumA1:tostratumB1"])+exp(mod_optim$results$beta$Psi["pop1:stratumA1:tostratumC1"])
cat(paste("PsiAA_1 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiAB_1 = ",exp(mod_optim$results$beta$Psi["pop1:stratumA1:tostratumB1"])/sum(1+sum(Tmp)),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop1:stratumB1:tostratumA1"])+exp(mod_optim$results$beta$Psi["pop1:stratumB1:tostratumC1"])
cat(paste("PsiBB_1 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiBA_1 = ",exp(mod_optim$results$beta$Psi["pop1:stratumB1:tostratumA1"])/sum(1+sum(Tmp)),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop1:stratumC1:tostratumA1"])+exp(mod_optim$results$beta$Psi["pop1:stratumC1:tostratumB1"])
cat(paste("PsiCC_1 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiCA_1 = ",exp(mod_optim$results$beta$Psi["pop1:stratumC1:tostratumA1"])/sum(1+sum(Tmp)),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop2:stratumA2:tostratumB2"])+exp(mod_optim$results$beta$Psi["pop2:stratumA2:tostratumC2"])
cat(paste("PsiAA_2 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiAB_2 = ",exp(mod_optim$results$beta$Psi["pop2:stratumA2:tostratumB2"])/sum(1+sum(Tmp)),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop2:stratumB2:tostratumA2"])+exp(mod_optim$results$beta$Psi["pop2:stratumB2:tostratumC2"])
cat(paste("PsiBB_2 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiBA_2 = ",exp(mod_optim$results$beta$Psi["pop2:stratumB2:tostratumA2"])/sum(1+sum(Tmp)),"\n"))
Tmp=exp(mod_optim$results$beta$Psi["pop2:stratumC2:tostratumA2"])+exp(mod_optim$results$beta$Psi["pop2:stratumC2:tostratumB2"])
cat(paste("PsiCC_2 = ",1/sum(1+sum(Tmp)),"\n"))
cat(paste("PsiCA_2 = ",exp(mod_optim$results$beta$Psi["pop2:stratumC2:tostratumA2"])/sum(1+sum(Tmp)),"\n"))


