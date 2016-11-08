source('./MSmixture/R/simulate_MSmixture.R')
library(marked)

#call simulate_MSmixture to simulate multistate mixture capture-recapture data
n_sims=10000
t_steps=10
n_strata=2
n_mixtures=2
N_releases = array(100,dim=c(n_mixtures,n_strata,t_steps))
N_releases[2,1,1]=50
N_releases[,2,]=0
Psi=vector("list",n_mixtures)
Psi[[1]]=Psi[[2]]=matrix(0,2,2)
Psi[[1]][1,1]=1.0
Psi[[2]][2,1]=Psi[[2]][1,2]=1
Phi=matrix(0.9,t_steps,n_mixtures)
Phi[,1]=0.95
P = matrix(0.4,t_steps,n_strata)
P[,2]=0
enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

Result=matrix(0,n_sims,5)

for(isim in 1:n_sims){
  set.seed(10000+isim)
  
  Hists = simulate_MSmixture(t_steps=t_steps,N_releases=N_releases,Psi=Psi,Phi=Phi,P=P,enc_hist_format=enc_hist_format,fname=fname)
  Hists$record = 1  #try to decrease dimensionality of dp call while retaining data frame
  dp=process.data(Hists,model="mvmscjs",strata.labels=list(area=c("A","B"),pop=c("1","2","u")))
  ddl=make.design.data(dp)
  # fix delta = 1 for unknown mixture (true mixture never known)
  ddl$delta$fix = ifelse(ddl$delta$obs.pop=="u", 1, 0)
  # Set Psi to 0 for cases which are not possible - moving between population subgroups
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$topop=="2"]=0
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$topop=="1"]=0
  # set Psi = 0 for residents moving and migrants not moving
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$toarea=="B"]=0
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$area=="B"&ddl$Psi$toarea=="A"&ddl$Psi$topop=="1"]=0
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$area=="B"&ddl$Psi$toarea=="B"&ddl$Psi$topop=="1"]=1 #note mixture 1 can't move to B, but transitions from B still need to sum to 1.0
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$area=="A"&ddl$Psi$toarea=="A"&ddl$Psi$topop=="1"]=1
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="A" & ddl$Psi$toarea=="A"] =0
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="B" & ddl$Psi$toarea=="B"] =0 #note these animals can't move to B, but sum from B must sum to 1.0
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="A" & ddl$Psi$toarea=="B" & ddl$Psi$topop=="2"] =1
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="B" & ddl$Psi$toarea=="A" & ddl$Psi$topop=="2"] =1
  ddl$p$fix[ddl$p$area=="B"]=0
  # Create indicator variables for transitioning between states
  #ddl$Psi$AtoB=ifelse(ddl$Psi$area=="A"&ddl$Psi$toarea=="B",1,0)  # A to B movement
  #ddl$Psi$BtoA=ifelse(ddl$Psi$area=="B"&ddl$Psi$toarea=="A",1,0)  # B to A movement
  #ddl$Psi$AtoC=ifelse(ddl$Psi$area=="A"&ddl$Psi$toarea=="C",1,0)  # A to C movement
  #ddl$Psi$CtoA=ifelse(ddl$Psi$area=="C"&ddl$Psi$toarea=="A",1,0)  # C to A movement
  #ddl$Psi$BtoC=ifelse(ddl$Psi$area=="B"&ddl$Psi$toarea=="C",1,0)  # B to C movement
  #ddl$Psi$CtoB=ifelse(ddl$Psi$area=="C"&ddl$Psi$toarea=="B",1,0)  # C to B movement
  # formulas
  #Psi.1=list(formula=~-1+AtoB+AtoC+BtoC+BtoA+CtoA+CtoB)
  #Psi.1=list(formula=~-1+AtoB+ AtoB:pop + AtoC + AtoC:pop + BtoC + BtoC:pop + BtoA + BtoA:pop + CtoA + CtoA:pop + CtoB + CtoB:pop)
  Psi.1=list(formula=~-1)
  p.1=list(formula=~pop)
  delta.1=list(formula= ~0)
  Phi.1=list(formula=~pop)
  #Phi.1=list(formula=~1)
  pi.1 = list(formula = ~0+pop)
  
  #initial= list(Phi=1.35,p=c(-1,0,1),Psi=c(-1,-1,1,3,3,1))
  #initial=list
  initial=NULL
  
  mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)
  
  
}
#mod_admb=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)

#pars on real scale
#cat(paste("p_mix1 = ",expit(mod_optim$results$beta$p["(Intercept)"]),"\n"))
#cat(paste("p_mix2 = ",expit(mod_optim$results$beta$p["(Intercept)"]+mod_optim$results$beta$p["pop2"]),"\n"))
#cat(paste("S1 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]),"\n"))
#cat(paste("S2 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]+mod_optim$results$beta$Phi["pop2"]),"\n"))


