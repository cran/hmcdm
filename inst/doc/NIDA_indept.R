## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hmcdm)

## -----------------------------------------------------------------------------
N = length(Test_versions)
J = nrow(Q_matrix)
K = ncol(Q_matrix)
L = nrow(Test_order)
Jt = J/L

## -----------------------------------------------------------------------------
Test_versions_sim <- sample(1:5,N,replace = L)
nClass = 2^K
R <- matrix(0,K,K)
tau <- numeric(K)
for(k in 1:K){
  tau[k] <- runif(1,.2,.6)
}
# Initial Alphas
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
# Subsequent Alphas
Alphas <- simulate_alphas_indept(tau,Alphas_0,L,R)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
Gmats <- array(NA,c(Jt,K,L))
Smats <- array(NA,c(Jt,K,L))
for(k in 1:K){
  Smats[,k,] <- runif(1,.1,.3)
  Gmats[,k,] <- runif(1,.1,.3)
}

Y_sim = simNIDA(Alphas,Smats[1,,1],Gmats[1,,1],Q_matrix,Test_order,Test_versions_sim)

## -----------------------------------------------------------------------------
output_NIDA_indept = hmcdm(Y_sim,Q_matrix,"NIDA_indept",Test_order,Test_versions_sim,100,30,
                                    R = R)
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

