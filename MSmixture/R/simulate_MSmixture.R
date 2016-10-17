#' function to simulate multistate mixture capture-recapture data
#' @param t_steps Number of time steps
#' @param N_releases A 3-d array specifying the number of releases by (suppopulation,stratum,time) 
#' @param Psi A vectored list, where each list element gives a transition matrix for each subpopulation (mixture)
#' @param Phi A matrix specifying survival of each population subgroup (column) by time (rows).
#' @param P A matrix specifying detection probability by stratum (column) and time (rows)
#' @param enc_hist_format Outputs encounter histories in format for multivariate state analysis in the marked package (enc_hist_format = 'marked'; default) or E-SURGE (enc_hist_format='esurge')
#' @param fname  A character string giving the filename for encounter history output; if null, just returns encounter histories in a character vector
#' @return A vector of character strings for each encounter history
#' @export
#' @keywords encounter history, simulation
#' @author Paul B. Conn
simulate_MSmixture <- function(t_steps,N_releases,Psi,Phi,P,enc_hist_format="marked",fname=NULL){
  n_strata = nrow(Psi[[1]])
  if(n_strata>10)cat("ERROR: currently limited to 10 strata")
  n_mix = length(Psi)
  E_hists = matrix(0,sum(N_releases),t_steps)
  State_list=c('A','B','C','D','E','F','G','H','I','J') #max 10 for now
  counter=0
  for(imix in 1:n_mix){
    for(istr in 1:n_strata){
      for(it in 1:t_steps){
        if(N_releases[imix,istr,it]>0){
          for(iind in 1:N_releases[imix,istr,it]){
            E_hists[counter+iind,it]=istr
            i_alive=1
            if(it<t_steps){
              cur_st=istr
              for(irem in it:(t_steps-1)){
                if(i_alive==1){
                  i_alive = rbinom(1,1,Phi[it,imix])
                  cur_st = sample(c(1:n_strata),1,prob=Psi[[imix]][cur_st,])
                  if(i_alive==1)E_hists[counter+iind,irem+1]=cur_st*rbinom(1,1,P[irem+1,cur_st])
                }
              }
            }
          }
          counter=counter+N_releases[imix,istr,it]
        }
      }
    }
  }
  # output results in desired format
  E_hists=matrix(as.character(E_hists),sum(N_releases),t_steps)
  if(enc_hist_format=='esurge'){
    #remove spaces 
    my_fun<-function(str)paste(str,collapse='')
    E_hists=apply(E_hists,1,'my_fun')
    E_hists=cbind(E_hists,rep("1;",sum(N_releases)))
    if(is.null(fname)==FALSE)write(t(E_hists),file=fname,ncolumns=ncol(E_hists))
  }
  else{
    for(i in 1:n_strata){
      Cur_which = which(E_hists == as.character(i))
      if(length(Cur_which)>0)E_hists[Cur_which]=State_list[i]
    }
    Which_obs = which(E_hists!='0')
    E_hists[Which_obs]=paste(E_hists[Which_obs],'u',sep='')
    my_fun<-function(str)paste0(str,collapse=',')
    E_hists=apply(E_hists,1,'my_fun')
    E_hists=data.frame(ch=E_hists,record=c(1:length(E_hists)))
    E_hists$ch=as.character(E_hists$ch)
    if(!is.null(fname))save(E_hists,file=fname)
  }
  E_hists
}