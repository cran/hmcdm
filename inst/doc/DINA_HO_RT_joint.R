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
Q_examinee <- Q_list(Q_matrix, Test_order, Test_versions)
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
mu_thetatau = c(0,0)
Sig_thetatau = rbind(c(1.8^2,.4*.5*1.8),c(.4*.5*1.8,.25))
Z = matrix(rnorm(N*2),N,2)
thetatau_true = Z%*%chol(Sig_thetatau)
thetas_true = thetatau_true[,1]
taus_true = thetatau_true[,2]
G_version = 3
phi_true = 0.8
for(i in 1:N){
  Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
lambdas_true <- c(-2, .4, .055)       # empirical from Wang 2017
Alphas <- simulate_alphas_HO_joint(lambdas_true,thetas_true,Alphas_0,Q_examinee,L,Jt)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
itempars_true <- array(runif(J*2,0.1,0.3), dim = c(Jt,2,L))
RT_itempars_true <- array(NA, dim = c(Jt,2,L))
RT_itempars_true[,2,] <- rnorm(Jt*L,3.45,.5)
RT_itempars_true[,1,] <- runif(Jt*L,1.5,2)

Y_sim <- simDINA(Alphas,itempars_true,ETAs,Test_order,Test_versions)
L_sim <- sim_RT(Alphas,RT_itempars_true,Q_matrix,taus_true,phi_true,ETAs,G_version,Test_order,Test_versions)

## -----------------------------------------------------------------------------
output_HMDCM_RT_joint = hmcdm(Y_sim,Q_matrix,"DINA_HO_RT_joint",Test_order,Test_versions,100,30,
                                 Latency_array = L_sim, G_version = G_version,
                                 theta_propose = 2,deltas_propose = c(.45,.25,.06))
output_HMDCM_RT_joint
summary(output_HMDCM_RT_joint)
a <- summary(output_HMDCM_RT_joint)
a
a$ss_EAP
head(a$ss_EAP)

## -----------------------------------------------------------------------------
(cor_thetas <- cor(thetas_true,a$thetas_EAP))
(cor_taus <- cor(taus_true,a$response_times_coefficients$taus_EAP))

(cor_ss <- cor(as.vector(itempars_true[,1,]),a$ss_EAP))
(cor_gs <- cor(as.vector(itempars_true[,2,]),a$gs_EAP))

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
head(a$PPP_total_RTs)
head(a$PPP_item_means)
head(a$PPP_item_mean_RTs)
head(a$PPP_item_ORs)

