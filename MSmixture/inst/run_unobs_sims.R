source('./MSmixture/R/simulate_MSmixture.R')
library(marked)
#call simulate_MSmixture to simulate multistate mixture capture-recapture data
n_sims=100
t_steps=10
n_strata=3
n_mixtures=2
N_releases = array(0,dim=c(n_mixtures,n_strata,t_steps))
P = matrix(0,t_steps,n_strata)
N_releases[,1:2,]=50
P[,1] = 0.6
P[,2] = 0.4

Psi=vector("list",n_mixtures)
Psi[[1]]=matrix(c(0.7,0.5,0.6,0.1,0.3,0.1,0,0,0),3,3)
Psi[[1]][,3]=1-rowSums(Psi[[1]])
Psi[[2]]=matrix(c(0.4,0,0.2,0,0.3,0,0,0,0),3,3)
Psi[[2]][,3]=1-rowSums(Psi[[2]])

Phi=matrix(0.9,t_steps,n_mixtures)
Phi[,1]=0.95
enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

Result=matrix(NA,n_sims,16)
colnames(Result) = c("S1","S2","p1","p2","pi","S1.se","S2.se","p1.se","p2.se","pi.se","Converge","S1.cover","S2.cover","p1.cover","p2.cover","pi.cover")

#stuff for delta method calcs
Par_DM = diag(5)   
Par_DM[2,1]=1
Par_DM[4,3]=1


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

Link = c("logit","logit","logit","logit","mlog","mlog","mlog","mlog","mlog","mlog","mlog","mlog","mlog")

Link_pars = c(logit(.95),logit(.9),logit(.4),logit(.4),0)   #link scale pars for coverage calculations


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
  # set Psi = 0 for moving to B for pop 2
  ddl$Psi$fix[ddl$Psi$toarea=="B"& ddl$Psi$pop=="2"]=0
  # set values for moving from B for pop 2 so that they're fixed and sum to one
  ddl$Psi$fix[ddl$Psi$area=="B" & ddl$Psi$toarea=="A"& ddl$Psi$pop=="2" & ddl$Psi$topop=="2"]=1
  ddl$Psi$fix[ddl$Psi$area=="B" & ddl$Psi$toarea=="C"& ddl$Psi$pop=="2" & ddl$Psi$topop=="2"]=0
  #create a factor for each possible estimated transition type (6 for pop 1, 2 for pop 2)
  ddl$Psi$estpar = paste0(ddl$Psi$area,ddl$Psi$toarea,ddl$Psi$pop)
  ddl$Psi$estpar[ddl$Psi$estpar %in% c("AA1","BB1","CC1","AA2","BB2","CC2","AB2","BA2","BC2","CB2")]="na"
  ddl$Psi$estpar = factor(ddl$Psi$estpar)

  ddl$p$fix[ddl$p$pop=="2" & ddl$p$area=="B"]=0  #2nd population never goes to winter area
  ddl$p$fix[ddl$p$area=="C"]=0  #unobservable state
  
  ddl$pi$fix[ddl$pi$area=="B" & ddl$pi$pop=="2"]=0

  Psi.1=list(formula=~0+estpar)
  p.1=list(formula=~area)
  delta.1=list(formula= ~0)
  Phi.1=list(formula=~pop)
  pi.1 = list(formula = ~0+pop*area)
  
  #initial= list(Phi=1.35,p=c(-1,0,1),Psi=c(-1,-1,1,3,3,1))
  #initial=list
  initial=NULL
  
  mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)
  
  Result[isim,"Converge"]=mod_optim$results$convergence
  if(Result[isim,"Converge"]==0){

    Pars_link = c(mod_optim$results$beta$Phi[1],sum(mod_optim$results$beta$Phi),mod_optim$results$beta$p[1],sum(mod_optim$results$beta$p[1:2]),sum(mod_optim$results$beta$p[c(1,3)]),mod_optim$results$beta$pi)
    Result[isim,1:6]=real.pars(Pars_link,Link=Link)
    V_link = diag(Par_DM %*% mod_optim$results$beta.vcv %*% t(Par_DM))  #variance on link scale
    Interv = cbind(Pars_link - 1.96*sqrt(V_link),Pars_link + 1.96*sqrt(V_link))
    Result[isim,7:12]= sqrt(link.deriv(Par_DM%*%Result[isim,1:6],Link=Link)^2*V_link)
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


