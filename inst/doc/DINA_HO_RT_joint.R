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

## -----------------------------------------------------------------------------
ETAs <- ETAmat(K, J, Q_matrix)
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
Alphas <- sim_alphas(model="HO_joint", 
                    lambdas=lambdas_true, 
                    thetas=thetas_true, 
                    Q_matrix=Q_matrix, 
                    Design_array=Design_array)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
itempars_true <- matrix(runif(J*2,.1,.2), ncol=2)
RT_itempars_true <- matrix(NA, nrow=J, ncol=2)
RT_itempars_true[,2] <- rnorm(J,3.45,.5)
RT_itempars_true[,1] <- runif(J,1.5,2)

Y_sim <- sim_hmcdm(model="DINA",Alphas,Q_matrix,Design_array,
                   itempars=itempars_true)
L_sim <- sim_RT(Alphas,Q_matrix,Design_array,
                  RT_itempars_true,taus_true,phi_true,G_version)

## -----------------------------------------------------------------------------
output_HMDCM_RT_joint = hmcdm(Y_sim,Q_matrix,"DINA_HO_RT_joint",Design_array,100,30,
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

(cor_ss <- cor(as.vector(itempars_true[,1]),a$ss_EAP))
(cor_gs <- cor(as.vector(itempars_true[,2]),a$gs_EAP))

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

