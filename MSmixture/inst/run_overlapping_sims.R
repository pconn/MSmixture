source('./MSmixture/R/simulate_MSmixture.R')
library(marked)
#call simulate_MSmixture to simulate multistate mixture capture-recapture data
n_sims=1000
t_steps=10
n_strata=3
n_mixtures=2
N_releases = array(0,dim=c(n_mixtures,n_strata,t_steps))
P = matrix(0,t_steps,n_strata)
N_releases[,1:3,]=20
P[,1] = 0.7
P[,2] = 0.5
P[,3] = 0.4

Psi=vector("list",n_mixtures)
Psi[[1]]=matrix(c(0.7,0.5,0.6,0.1,0.3,0.1,0,0,0),3,3)
Psi[[1]][,3]=1-rowSums(Psi[[1]])
Psi[[2]]=matrix(c(0.3,0.2,0.2,0.5,0.6,0.4,0,0,0),3,3)
Psi[[2]][,3]=1-rowSums(Psi[[2]])

Phi=matrix(0.9,t_steps,n_mixtures)
Phi[,1]=0.95
enc_hist_format='marked'
fname=NULL
#enc_hist_format='esurge'
#fname='multistate.inp'

Result=matrix(NA,n_sims,61)
colnames(Result) = c("S1","S2","pA","pB","pC","psiAB1","psiAB2","psiAC1","psiAC2","psiBA1","psiBA2","psiBC1","psiBC2","psiCA1","psiCA2","psiCB1","psiCB2","pi.2A","pi.2B","pi.2C",
  "S1.se","S2.se","pA.se","pB.se","pC.se","psiAB1.se","psiAB2.se","psiAC1.se","psiAC2.se","psiBA1.se","psiBA2.se","psiBC1.se","psiBC2.se","psiCA1.se","psiCA2.se","psiCB1.se","psiCB2.se","pi.2A.se","pi.2B.se","pi.2C.se",
  "S1.cover","S2.cover","pA.cover","pB.cover","pC.cover","psiAB1.cover","psiAB2.cover","psiAC1.cover","psiAC2.cover","psiBA1.cover","psiBA2.cover","psiBC1.cover","psiBC2.cover","psiCA1.cover","psiCA2.cover","psiCB1.cover","psiCB2.cover","pi.2A.cover","pi.2B.cover","pi.2C.cover","Converge")

Opposite_order = c(2,1,3,4,5,7,6,9,8,11,10,13,12,15,14,17,16)  #for reaarranging S, p, psi if label switching happens

#stuff for logit delta method calcs
n_logit=5
Par_DM = diag(5)   
Par_DM[2,1]=1
Par_DM[4,3]=1
Par_DM[5,3]=1
Par_DM_all = diag(20)
Par_DM_all[1:5,1:5]=Par_DM

#stuff for multinom logit delta method
Denom = vector("list",15)
Denom[[1]]=c(1,3)
Denom[[2]]=c(2,4)
Denom[[3]]=c(1,3)
Denom[[4]]=c(2,4)
Denom[[5]]=c(5,7)
Denom[[7]]=c(5,7)
Denom[[6]]=c(6,8)
Denom[[8]]=c(6,8)
Denom[[9]]=c(9,11)
Denom[[11]]=c(9,11)
Denom[[10]]=c(10,12)
Denom[[12]]=c(10,12)
Denom[[13]]=13
Denom[[14]]=14
Denom[[15]]=15



expit.deriv <- function(x){
  exp(x)/(1+exp(x))^2
}

mlog.deriv <- function(X,num){  #X gives vector of m-logit pars, num gives the value that's in the numerator
  #-exp(x)/(1+exp(x))^2
  (1+sum(exp(X))-exp(num))*exp(num)/(1+sum(exp(X)))^2
}

link.deriv <- function(x,Link,Denom,n_logit) {
  Derivs = rep(0,length(x))
  for(i in 1:length(x)){
    if(Link[i]=="logit")Derivs[i]=expit.deriv(x[i])
    if(Link[i]=="mlog")Derivs[i]=mlog.deriv(X=x[Denom[[i-n_logit]]+n_logit],num=x[i])
  }
  Derivs
}

real.pars <-function(x,Link,Denom,n_logit){
  Reals = rep(0,length(x))
  for(i in 1:length(x)){
    if(Link[i]=="logit")Reals[i]=expit(x[i])
    if(Link[i]=="mlog")Reals[i]=exp(x[i])/(1+sum(exp(x[Denom[[i-n_logit]]+n_logit])))   
  }
  Reals
}

Link = c(rep("logit",5),rep("mlog",15))

#link scale pars for coverage calculations
Link_pars = c(logit(.95),logit(.9),logit(.7),logit(.5),logit(.4),
      log(Psi[[1]][1,2]/Psi[[1]][1,1]),log(Psi[[2]][1,2]/Psi[[2]][1,1]),log(Psi[[1]][1,3]/Psi[[1]][1,1]),log(Psi[[2]][1,3]/Psi[[2]][1,1]),   
      log(Psi[[1]][2,1]/Psi[[1]][2,2]),log(Psi[[2]][2,1]/Psi[[2]][2,2]),log(Psi[[1]][2,3]/Psi[[1]][2,2]),log(Psi[[2]][2,3]/Psi[[2]][2,2]),   #link scale pars for coverage calculations
      log(Psi[[1]][3,1]/Psi[[1]][3,3]),log(Psi[[2]][3,1]/Psi[[2]][3,3]),log(Psi[[1]][3,2]/Psi[[1]][3,3]),log(Psi[[2]][3,2]/Psi[[2]][3,3]),   #link scale pars for coverage calculations
      0, 0, 0)

for(isim in 1:n_sims){
  set.seed(10000+isim)
  
  Hists = simulate_MSmixture(t_steps=t_steps,N_releases=N_releases,Psi=Psi,Phi=Phi,P=P,enc_hist_format=enc_hist_format,fname=fname)
  Hists$record = 1  #try to decrease dimensionality of dp call while retaining data frame
  dp=process.data(Hists,model="mvmscjs",strata.labels=list(area=c("A","B","C"),pop=c("1","2","u")))
  ddl=make.design.data(dp)
  # fix delta = 1 for unknown mixture (true mixture never known)
  ddl$delta$fix = ifelse(ddl$delta$obs.pop=="u", 1, 0)
  
  #create a factor for each possible estimated transition type (6 for pop 1, 2 for pop 2)
  ddl$Psi$fix[ddl$Psi$pop=="1" & ddl$Psi$topop=="2"]=0
  ddl$Psi$fix[ddl$Psi$pop=="2" & ddl$Psi$topop=="1"]=0
  ddl$Psi$estpar = paste0(ddl$Psi$area,ddl$Psi$toarea,ddl$Psi$pop)
  ddl$Psi$estpar[ddl$Psi$estpar %in% c("AA1","BB1","CC1","AA2","BB2","CC2")]="na"
  ddl$Psi$estpar = factor(ddl$Psi$estpar)


  Psi.1=list(formula=~0+estpar)
  p.1=list(formula=~area)
  delta.1=list(formula= ~0)
  Phi.1=list(formula=~pop)
  pi.1 = list(formula = ~0+area)
  
  #initial= list(Phi=1.35,p=c(-1,0,1),Psi=c(-1,-1,1,3,3,1))
  #initial=list
  initial=NULL
  
  mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)
  #if(mod_optim$results$convergence==1)mod_optim=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=F,run=T)
  
  Result[isim,"Converge"]=mod_optim$results$convergence
  if(Result[isim,"Converge"]==0){
    Pars_link = c(Par_DM%*%c(mod_optim$results$beta$Phi,mod_optim$results$beta$p),mod_optim$results$beta$Psi,mod_optim$results$beta$pi)
    Reals = real.pars(x=Pars_link,Link=Link,Denom=Denom,n_logit=n_logit)
    V_link = diag(Par_DM_all %*% mod_optim$results$beta.vcv %*% t(Par_DM_all))  #variance on link scale
    Interv = cbind(Pars_link - 1.96*sqrt(V_link),Pars_link + 1.96*sqrt(V_link))
    SE =  sqrt(link.deriv(x=Pars_link,Link=Link,Denom=Denom,n_logit=n_logit)^2*V_link)
    Coverage = (Link_pars>Interv[,1] & Link_pars<Interv[,2])
    if((Reals[7]+Reals[9])<(Reals[6]+Reals[8])){  #define label switching if psi_1^AA < psi_2^AA
      Reals[1:17] = Reals[Opposite_order]
      Reals[18:20] = 1-Reals[18:20]
      SE[1:17] = SE[Opposite_order]
      Coverage[1:17] = (Link_pars[Opposite_order]>Interv[1:17,1] & Link_pars[Opposite_order]<Interv[1:17,2])
      Coverage[1:17] = Coverage[Opposite_order]
    }
    Result[isim,1:20]=Reals
    Result[isim,21:40]= SE
    Result[isim,41:60]= Coverage
  }
}

save(Result,file="overlapping_results.Rdata")
#mod_admb=crm(dp,ddl,model="mvmscjs",initial=initial,model.parameters=list(Psi=Psi.1,pi=pi.1,p=p.1,delta=delta.1,Phi=Phi.1),hessian=TRUE,use.admb=T,run=T)

#pars on real scale
#cat(paste("p_mix1 = ",expit(mod_optim$results$beta$p["(Intercept)"]),"\n"))
#cat(paste("p_mix2 = ",expit(mod_optim$results$beta$p["(Intercept)"]+mod_optim$results$beta$p["pop2"]),"\n"))
#cat(paste("S1 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]),"\n"))
#cat(paste("S2 = ",expit(mod_optim$results$beta$Phi["(Intercept)"]+mod_optim$results$beta$Phi["pop2"]),"\n"))

#summary taking out NA estimates
is.any.na <-function(x){
  cur.na = sum(is.na(x))
  if(cur.na>0)cur.na=1
  cur.na
}
NA_ind = apply(Result,1,is.any.na)
colMeans(Result[-which(NA_ind==1),])
