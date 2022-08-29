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
ETAs <- ETAmat(K, J, Q_matrix)
TP <- TPmat(K)
Omega_true <- rOmega(TP)
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
for(i in 1:N){
  Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
Alphas <- simulate_alphas_FOHM(Omega_true, Alphas_0,L)
itempars_true <- array(runif(Jt*2*L,.1,.2), dim = c(Jt,2,L))

Y_sim <- simDINA(Alphas,itempars_true,ETAs,Test_order,Test_versions)

## -----------------------------------------------------------------------------
output_FOHM = hmcdm(Y_sim,Q_matrix,"DINA_FOHM",Test_order,Test_versions,100,30)
output_FOHM
summary(output_FOHM)
a <- summary(output_FOHM)
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

