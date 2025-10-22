## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hmcdm)

## -----------------------------------------------------------------------------
N = dim(Design_array)[1]
J = nrow(Q_matrix)
K = ncol(Q_matrix)
L = dim(Design_array)[3]

## -----------------------------------------------------------------------------
tau <- numeric(K)
for(k in 1:K){
  tau[k] <- runif(1,.2,.6)
}
R = matrix(0,K,K)
# Initial alphas
p_mastery <- c(.5,.5,.4,.4)
Alphas_0 <- matrix(0,N,K)
for(i in 1:N){
  for(k in 1:K){
    prereqs <- which(R[k,]==1)
    if(length(prereqs)==0){
      Alphas_0[i,k] <- rbinom(1,1,p_mastery[k])
    }
    if(length(prereqs)>0){
      Alphas_0[i,k] <- prod(Alphas_0[i,prereqs])*rbinom(1,1,p_mastery)
    }
  }
}
Alphas <- sim_alphas(model="indept",taus=tau,N=N,L=L,R=R,alpha0=Alphas_0)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
Svec <- runif(K,.1,.3)
Gvec <- runif(K,.1,.3)

Y_sim <- sim_hmcdm(model="NIDA",Alphas,Q_matrix,Design_array,
                   Svec=Svec,Gvec=Gvec)

## -----------------------------------------------------------------------------
output_NIDA_indept = hmcdm(Y_sim, Q_matrix, "NIDA_indept", Design_array,
                           100, 30, R = R)
output_NIDA_indept
summary(output_NIDA_indept)
a <- summary(output_NIDA_indept)
head(a$ss_EAP)

## -----------------------------------------------------------------------------
AAR_vec <- numeric(L)
for(t in 1:L){
  AAR_vec[t] <- mean(Alphas[,,t]==a$Alphas_est[,,t])
}
AAR_vec

PAR_vec <- numeric(L)
for(t in 1:L){
  PAR_vec[t] <- mean(rowSums((Alphas[,,t]-a$Alphas_est[,,t])^2)==0)
}
PAR_vec


## -----------------------------------------------------------------------------
a$DIC
head(a$PPP_total_scores)
head(a$PPP_item_means)
head(a$PPP_item_ORs)

