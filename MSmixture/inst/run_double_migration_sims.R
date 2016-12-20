#run this one in output directory
setwd('c:/users/paul.conn/git/MSmixture/output')
source('../MSmixture/R/simulate_MSmixture.R')
library(marked)
#call simulate_MSmixture to simulate multistate mixture capture-recapture data
n_sims=1000
t_steps=10
n_strata=3
n_mixtures=2
N_releases = array(0,dim=c(n_mixtures,n_strata,t_steps))
P = matrix(0,t_steps,n_strata)
for(it in 1:5){
  N_releases[1:2,1,it*2-1] = 50
  #N_releases[1,2,it*2]= 50
  #N_releases[2,3,it*2]=50
  P[it*2-1,1]=0.4
  P[it*2,2]=P[it*2,3]=0.2
}
Psi=vector("list",n_mixtures)
Psi[[1]]=Psi[[2]]=matrix(0,3,3)
Psi[[1]][2,1]=Psi[[1]][1,2]=1
Psi[[2]][3,1]=Psi[[2]][1,3]=1
Phi=matrix(0.9,t_steps,n_mixtures)
Phi[,1]=0.95
enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

Result=matrix(NA,n_sims,19)
colnames(Result) = c("S1","S2","p1","p2","p3","pi","S1.se","S2.se","p1.se","p2.se","p3.se","pi.se","Converge","S1.cover","S2.cover","p1.cover","p2.cover","p3.cover","pi.cover")

#stuff for delta method calcs
Par_DM = diag(6)   
Par_DM[2,1]=1
Par_DM[4,3]=1
Par_DM[5,3]=1

expit.deriv <- function(x){
  exp(x)/(1+exp(x))^2
}

mlog.deriv <- function(x){  #note that in the bivariate case that mlogit is equal to logit
  #-exp(x)/(1+exp(x))^2
  exp(x)/(1+exp(x))^2
}

link.deriv <- function(x,Link) {
  Derivs = rep(0,length(x))
  for(i in 1:length(x)){
    if(Link[i]=="logit")Derivs[i]=expit.deriv(x[i])
    if(Link[i]=="mlog")Derivs[i]=mlog.deriv(x[i])
  }
  Derivs
}

real.pars <-function(x,Link){
  Reals = rep(0,length(x))
  for(i in 1:length(x)){
    if(Link[i]=="logit")Reals[i]=expit(x[i])
    if(Link[i]=="mlog")Reals[i]=expit(x[i])   #1/(1+exp(x[i]))
  }
  Reals
}

Link = c("logit","logit","logit","logit","logit","mlog")

Link_pars =  c(logit(.95),logit(.9),logit(.4),logit(.2),logit(.2),0)   #link scale pars for coverage calculations


for(isim in 1:n_sims){
  set.seed(10000+isim)
  
  Hists = simulate_MSmixture(t_steps=t_steps,N_releases=N_releases,Psi=Psi,Phi=Phi,P=P,enc_hist_format=enc_hist_format,fname=fname)
  Hists$record = 1  #try to decrease dimensionality of dp call while retaining data frame
  dp=process.data(Hists,model="mvmscjs",strata.labels=list(area=c("A","B","C"),pop=c("1","2","u")))
  ddl=make.design.data(dp)
  # fix delta = 1 for unknown mixture (true mixture never known)
  ddl$delta$fix = ifelse(ddl$delta$obs.pop=="u", 1, 0)
  # Set Psi to 0 for cases which are not possible - moving between population subgroups
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$topop=="2"]=0
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$topop=="1"]=0
  # set Psi = 0 for mixture 1 moving to 3, mixture 2 moving to 2, and staying in the same place
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$toarea=="C"]=0
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$toarea=="B"]=0
  ddl$Psi$fix[ddl$Psi$area=="A"&ddl$Psi$toarea=="A"]=0
  ddl$Psi$fix[ddl$Psi$area=="B"&ddl$Psi$toarea=="B"]=0
  ddl$Psi$fix[ddl$Psi$area=="C"&ddl$Psi$toarea=="C"]=0
  ddl$Psi$fix = 0
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$area=="A" & ddl$Psi$toarea=="B" & ddl$Psi$topop=="1"] =1
  ddl$Psi$fix[ddl$Psi$pop=="1"&ddl$Psi$area=="B" & ddl$Psi$toarea=="A" & ddl$Psi$topop=="1"] =1
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="A" & ddl$Psi$toarea=="C" & ddl$Psi$topop=="2"] =1
  ddl$Psi$fix[ddl$Psi$pop=="2"&ddl$Psi$area=="C" & ddl$Psi$toarea=="A" & ddl$Psi$topop=="2"] =1
  ddl$Psi$fix[ddl$Psi$stratum=="B2"&ddl$Psi$tostratum=="A2" & (as.numeric(as.character(ddl$Psi$time))%%2)==0]=1
  ddl$Psi$fix[ddl$Psi$stratum=="B1"&ddl$Psi$tostratum=="A1" & (as.numeric(as.character(ddl$Psi$time))%%2)==0]=1
  ddl$Psi$fix[ddl$Psi$stratum=="C2"&ddl$Psi$tostratum=="A2" & (as.numeric(as.character(ddl$Psi$time))%%2)==0]=1
  ddl$Psi$fix[ddl$Psi$stratum=="C1"&ddl$Psi$tostratum=="A1" & (as.numeric(as.character(ddl$Psi$time))%%2)==0]=1
  #need to make sure transitions from states all sum to one even when there are no animals there
  ddl$Psi$fix[(ddl$Psi$area=="B" | ddl$Psi$area=="C") & (as.numeric(as.character(ddl$Psi$time))%%2)==1]=0  
  ddl$Psi$fix[(ddl$Psi$area=="B" | ddl$Psi$area=="C") & (as.numeric(as.character(ddl$Psi$time))%%2)==1 & ddl$Psi$tostratum=="B1"]=1  
  ddl$p$fix = 0
  ddl$p$fix[ddl$p$area=="A" & (ddl$p$Time%%2 == 1)]=NA
  ddl$p$fix[ddl$p$area=="B" & (ddl$p$Time%%2 == 0)]=NA
  ddl$p$fix[ddl$p$area=="C" & (ddl$p$Time%%2 == 0)]=NA
  ddl$pi$fix[ddl$pi$area=="B" & ddl$pi$pop=="2" & ddl$pi$Time%%2==1 & is.na(ddl$pi$fix)==1]=0
  ddl$pi$fix[ddl$pi$area=="C" & ddl$pi$pop=="1" & ddl$pi$Time%%2==1 & ddl$pi$fix==1]=0
  ddl$pi$fix[ddl$pi$area=="C" & ddl$pi$pop=="2" & ddl$pi$Time%%2==1 & is.na(ddl$pi$fix)==1]=1
  

  Psi.1=list(formula=~-1)
  p.1=list(formula=~area)
  delta.1=list(formula= ~0)
  Phi.1=list(formula=~pop)
  pi.1 = list(formula = ~0+pop)
  
  #initial= list(Phi=1.35,p=c(-1,0,1),Psi=c(-1,-1,1,3,3,1))
  #initial=list
  initial=NULL
  
  mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)
  if(mod_optim$results$convergence==1)mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)

    
  Result[isim,"Converge"]=mod_optim$results$convergence
  if(Result[isim,"Converge"]==0){

    Pars_link = c(mod_optim$results$beta$Phi[1],sum(mod_optim$results$beta$Phi),mod_optim$results$beta$p[1],sum(mod_optim$results$beta$p[1:2]),sum(mod_optim$results$beta$p[c(1,3)]),mod_optim$results$beta$pi)
    Result[isim,1:6]=real.pars(Pars_link,Link=Link)
    V_link = diag(Par_DM %*% mod_optim$results$beta.vcv %*% t(Par_DM))  #variance on link scale
    Interv = cbind(Pars_link - 1.96*sqrt(V_link),Pars_link + 1.96*sqrt(V_link))
    Result[isim,7:12]= sqrt(link.deriv(Pars_link[1:6],Link=Link)^2*V_link)
    Result[isim,14:19]= (Link_pars>Interv[,1] & Link_pars<Interv[,2])
  }
}

save(Result,file="double_migration_results.Rdata")
#mod_admb=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)

#pars on real scale
#cat(paste("p_mix1 = ",expit(mod_optim$results$beta$p["(Intercept)"]),"\n"))
#cat(paste("p_mix2 = ",expit(mod_optim$results$beta$p["(Intercept)"]+mod_optim$results$beta$p["pop2"]),"\n"))
#cat(paste("S1 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]),"\n"))
#cat(paste("S2 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]+mod_optim$results$beta$Phi["pop2"]),"\n"))


