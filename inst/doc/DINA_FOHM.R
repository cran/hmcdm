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
TP <- TPmat(K)
Omega_true <- rOmega(TP)
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
for(i in 1:N){
 Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
Alphas <- sim_alphas(model="FOHM", Omega = Omega_true, N=N, L=L)
itempars_true <- matrix(runif(J*2,.1,.2), ncol=2)

Y_sim <- sim_hmcdm(model="DINA",Alphas,Q_matrix,Design_array,
                   itempars=itempars_true)

## -----------------------------------------------------------------------------
output_FOHM = hmcdm(Y_sim,Q_matrix,"DINA_FOHM",Design_array,100,30)
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


