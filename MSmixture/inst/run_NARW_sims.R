#### run NARW sims

source('./MSmixture/R/simulate_MSmixture_NARW.R')
library(marked)
#call simulate_MSmixture to simulate multistate mixture capture-recapture data
n_sims=1
t_steps=22  #6 mo intervals
t_steps_yr = t_steps/2
n_strata=5   # includes unobservable state
n_mixtures=6  #females, males * fundy, non-fundy, offshore
N_releases = array(0,dim=c(n_mixtures,n_strata,t_steps))
n_unique=320 #approx number of unique animals observed in sigthing efforts (368 for 1980-1987)
n_unique=2000 # temporarily increase to look at identifiability

N_obs = matrix(0,18,9)    # number of observations reported by Wade and Clapham 1980-1997
#areas are SE U.S., Mass Bay, Grt S. Channel, Bay of Fundy, Scotian Shelf, Mid-Atlantic, Jeffrey's Ledge, Gulf of Maine, Other
N_obs[1,] = c(0,5,19,23,13,2,0,17,0)
N_obs[2,] = c(2,8,62,56,0,0,2,2,1)
N_obs[3,] = c(3,13,20,52,30,2,3,4,0)
N_obs[4,] = c(3,30,7,26,26,2,1,3,0)
N_obs[5,] = c(15,28,22,54,29,5,6,5,0)
N_obs[6,] = c(6,52,31,31,11,2,6,8,0)
N_obs[7,] = c(23,54,27,42,80,5,5,1,0)
N_obs[8,] = c(14,39,40,31,83,2,0,1,1)
N_obs[9,] = c(9,29,94,44,118,0,4,6,0)
N_obs[10,] = c(41,21,39,74,115,4,6,3,0)
N_obs[11,] = c(32,29,0,68,47,1,6,2,1)
N_obs[12,] = c(20,33,8,51,90,1,1,0,0)
N_obs[13,] = c(26,45,0,72,17,0,4,4,0)
N_obs[14,] = c(53,41,1,148,1,5,1,6,0)
N_obs[15,] = c(24,26,2,185,0,3,11,1,0)
N_obs[16,] = c(24,63,2,183,5,1,0,1,0)
N_obs[17,] = c(93,78,4,161,0,5,0,45,0)
N_obs[18,] = c(24,50,8,169,2,0,2,5,0)


#strata are SE US, Gulf of Maine + Mass Bay, Offshore, Bay of Fundy, Unobservable
N_obs_model = cbind(N_obs[c(2:12),1],N_obs[c(2:12),2]+N_obs[c(2:12),8],N_obs[c(2:12),3]+N_obs[c(2:12),5],N_obs[c(2:12),4],0)

P = matrix(0,t_steps,n_strata)
P[c(1:t_steps_yr)*2,1] = N_obs_model[,1]/80
P[c(1:t_steps_yr)*2-1,2] = N_obs_model[,2]/75
P[c(1:t_steps_yr)*2-1,3] = N_obs_model[,3]/230
P[c(1:t_steps_yr)*2-1,4] = N_obs_model[,4]/120

Obs_by_year = apply(N_obs,1,'sum')
New_IDs = c(46,53,24,16,14,17,22,17,18,20,16,9,12,12,10,3,13,0)
Prob_new_year = New_IDs[2:12]/sum(New_IDs[2:12])
Releases_by_year = round(Prob_new_year * n_unique)

Psi=vector("list",n_mixtures)  #1st list element gives mixture, 2nd gives to/from breeding areas
#fundy females
Psi[[1]][[1]] = t(matrix(c(1,0,0,0,0,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7),5,5))
Psi[[1]][[2]] = t(matrix(c(0,0.3,0.1,0.6,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.3,0.1,0.6,0),5,5))
#non-fundy females 
Psi[[2]][[1]] = t(matrix(c(1,0,0,0,0,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7),5,5))
Psi[[2]][[2]] = t(matrix(c(0,0.8,0.2,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.8,0.2,0,0),5,5))
#offshore females
Psi[[3]][[1]] = t(matrix(c(1,0,0,0,0,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7,0.3,0,0,0,0.7),5,5))
Psi[[3]][[2]] = t(matrix(c(0,0.2,0.6,0.2,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.2,0.6,0.2,0),5,5))
#fundy males
Psi[[4]][[1]] = t(matrix(c(1,0,0,0,0,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9),5,5))
Psi[[4]][[2]] = t(matrix(c(0,0.3,0.1,0.6,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.3,0.1,0.6,0),5,5))
#non-fundy males 
Psi[[5]][[1]] = t(matrix(c(1,0,0,0,0,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9),5,5))
Psi[[5]][[2]] = t(matrix(c(0,0.8,0.2,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.8,0.2,0,0),5,5))
#offshore males
Psi[[6]][[1]] = t(matrix(c(1,0,0,0,0,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9,0.1,0,0,0,0.9),5,5))
Psi[[6]][[2]] = t(matrix(c(0,0.2,0.6,0.2,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0.2,0.6,0.2,0),5,5))

Phi=matrix(0,t_steps,n_mixtures)
for(it in 1:(t_steps/2)){
  for(im in 1:n_mixtures){
    Phi[it*2-1,im]=Phi[it*2,im] = expit(5 + (it-1)*-0.1 + (im==3 | im==6)*.4 + (im>3)*0.4)
  }
}
Phi[,1:2]=0.85
Phi[,3:5]=0.9
Phi[,6]=0.95

Release_prob_by_mixture = matrix(0,n_mixtures,n_strata)
Release_prob_by_mixture[1,] = c(0.3,0.3,0.1,0.6,0)
Release_prob_by_mixture[2,] = c(0.3,0.8,0.2,0,0)
Release_prob_by_mixture[3,] = c(0.3,0.2,0.6,0.2,0)
Release_prob_by_mixture[4,] = c(0.1,0.3,0.1,0.6,0)
Release_prob_by_mixture[5,] = c(0.1,0.8,0.2,0,0)
Release_prob_by_mixture[6,] = c(0.1,0.2,0.6,0.2,0)
Release_prob_by_mixture = Release_prob_by_mixture/rowSums(Release_prob_by_mixture)


enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

#Result=matrix(NA,n_sims,19)
#colnames(Result) = c("S1","S2","p1","p2","p3","pi","S1.se","S2.se","p1.se","p2.se","p3.se","pi.se","Converge","S1.cover","S2.cover","p1.cover","p2.cover","p3.cover","pi.cover")

#stuff for delta method calcs
#Par_DM = diag(6)   
#Par_DM[2,1]=1
#Par_DM[4,3]=1
#Par_DM[5,3]=1

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

Mixture_proportions = c(0.2,0.2,0.1,0.15,0.2,0.15)
Releases_per_mixture = matrix(0,n_mixtures,t_steps_yr)
#Link = c("logit","logit","logit","logit","logit","mlog")

#Link_pars =  c(logit(.95),logit(.9),logit(.4),logit(.2),logit(.2),0)   #link scale pars for coverage calculations

for(isim in 1:n_sims){
  set.seed(10000+isim)

  #base N_releases by mixture and strata on total # of releases
  #number of releases per mixture
  for(it in 1:t_steps_yr)Releases_per_mixture[,it]=tabulate(sample(x=c(1:6),size=Releases_by_year[it],prob=Mixture_proportions,replace=TRUE))
  N_releases = array(0,dim=c(n_mixtures,n_strata,t_steps))
  for(imix in 1:n_mixtures){
    for(it in 1:t_steps_yr){
      if(Releases_per_mixture[imix,it]>0){
        Cur_release = tabulate(sample(c(1:n_strata),Releases_per_mixture[imix,it],prob=Release_prob_by_mixture[imix,],replace=TRUE),nbins=n_strata)
        N_releases[imix,,it*2-1]=c(0,Cur_release[2:n_strata])
        N_releases[imix,1,it*2]=Cur_release[1]
      }
    }
  }
    
  t_steps = 8   #temporarily decrease to 8 periods to get it to run faster
  P[2,1] = P[4,1] = P[6,1] = P[8,1] = 0.5
  P[1,2] = P[3,2] = P[5,2] = P[7,2] = 0.6
  P[1,3] = P[3,3] = P[5,3] = P[7,3] = 0.4
  P[1,4] = P[3,4] = P[5,4] = P[7,4] = 0.8

  Hists = simulate_MSmixture_NARW(t_steps=t_steps,N_releases=N_releases,Psi=Psi,Phi=Phi,P=P,enc_hist_format=enc_hist_format,fname=fname)
  Hists=Hists[-which(Hists[,1]=="0,0,0,0,0,0,0,0"),]
  
  Hists$record = 1  #try to decrease dimensionality of dp call while retaining data frame
  dp=process.data(Hists,model="mvmscjs",strata.labels=list(area=c("A","B","C","D","E"),pop=c("1","2","3","4","5","6","u")))
  cpu_time1 =  proc.time()
  ddl=make.design.data(dp)
  proc.time() - cpu_time1
  # fix delta = 1 for unknown mixture (true mixture never known)
  ddl$delta$fix = ifelse(ddl$delta$obs.pop=="u", 1, 0)
  
  ddl$Psi$fix = 0
  #first, set the state transitions to get by subtraction.  This is "B", "C", or "D" to "E" and  "A" to "B" and "E" to "B"
  ddl$Psi$fix[ddl$Psi$area %in% c("B","C","D") & ddl$Psi$toarea=="E"] = 1
  ddl$Psi$fix[ddl$Psi$area == "A" & ddl$Psi$toarea=="B"] = 1
  ddl$Psi$fix[ddl$Psi$area == "E" & ddl$Psi$toarea=="B"] = 1
  
  #set the state transitions to estimate.  B,C,D --> A in odd years, A,E --> C,D in even years
  ddl$Psi$fix[ddl$Psi$area %in% c("B","C","D") & ddl$Psi$toarea=="A" & ddl$Psi$Time%%2==0] = NA
  ddl$Psi$fix[ddl$Psi$area %in% c("A","E") & ddl$Psi$toarea %in% c("C","D") & ddl$Psi$Time%%2==1] = NA

  #non-fundy whales can't move to fundy
  ddl$Psi$fix[ddl$Psi$pop %in% c("2","5") & ddl$Psi$toarea=="D"] = 0
    
  #set all transitions between subpopulations to zero
  ddl$Psi$fix[ddl$Psi$pop=="1" & ddl$Psi$topop!="1"]=0
  ddl$Psi$fix[ddl$Psi$pop=="2" & ddl$Psi$topop!="2"]=0
  ddl$Psi$fix[ddl$Psi$pop=="3" & ddl$Psi$topop!="3"]=0
  ddl$Psi$fix[ddl$Psi$pop=="4" & ddl$Psi$topop!="4"]=0
  ddl$Psi$fix[ddl$Psi$pop=="5" & ddl$Psi$topop!="5"]=0
  ddl$Psi$fix[ddl$Psi$pop=="6" & ddl$Psi$topop!="6"]=0
  

  ddl$p$fix = 0
  ddl$p$fix[ddl$p$area == "E"] = 0
  ddl$p$fix[ddl$p$area=="A" & (ddl$p$Time%%2 == 0)]=NA
  ddl$p$fix[ddl$p$area %in% c("B","C","D") & (ddl$p$Time%%1 == 0)]=NA
  
  #create another column to specify psi parameters to be estimated 
  ddl$Psi$estpar = paste0(ddl$Psi$area,ddl$Psi$toarea,ddl$Psi$pop)
  ddl$Psi$estpar[ddl$Psi$estpar %in% c("BA1","CA1","DA1","BA2","CA2","DA2","BA3","CA3","DA3")]="FemToA"    #set all female movement rates to SE US the same
  #set return rates to summer strata equal between SEUS and 'unkown'
  ddl$Psi$estpar[ddl$Psi$estpar == "EC1"] = "AC1"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED1"] = "AD1"
  ddl$Psi$estpar[ddl$Psi$estpar == "EC2"] = "AC2"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED2"] = "AD2"
  ddl$Psi$estpar[ddl$Psi$estpar == "EC3"] = "AC3"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED3"] = "AD3"
  ddl$Psi$estpar[ddl$Psi$estpar == "EC4"] = "AC4"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED4"] = "AD4"
  ddl$Psi$estpar[ddl$Psi$estpar == "EC5"] = "AC5"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED5"] = "AD5"
  ddl$Psi$estpar[ddl$Psi$estpar == "EC6"] = "AC6"
  ddl$Psi$estpar[ddl$Psi$estpar == "ED6"] = "AD6"
  ddl$Psi$estpar = factor(ddl$Psi$estpar)
  
  #create another column to specify pi parameters to be estimated
  ddl$pi$estpar = paste0(ddl$pi$area,ddl$pi$pop)
  #set pi for fundy whales to zero for strata D ("fundy")
  ddl$pi$fix[ddl$pi$area=="D" & ddl$pi$pop %in% c("2","5")] = 0
  
  
  Psi.1=list(formula=~0+estpar)
  p.1=list(formula=~area)
  delta.1=list(formula= ~0)
  Phi.1=list(formula=~pop)
  pi.1 = list(formula = ~0+estpar)
  
  initial= list(Phi=c(2.2,0,0,0,0,0),p=c(0,1,0,2),Psi=c(-1,-1,1,-1,-1,1,1,0,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1),pi=c(0,0,-1,-1,-1,-1,-1,0,-1,-1,-1,1,-1,-1,1,-1,1,-1))
  #initial=list
  #initial=NULL
  
  mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)
  #if(mod_optim$results$convergence==1)mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)
  
  
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


