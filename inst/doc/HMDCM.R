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
Q_examinee <- Q_list(Q_matrix, Test_order, Test_versions)
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
thetas_true = rnorm(N)
for(i in 1:N){
  Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
lambdas_true = c(-1, 1.8, .277, .055)
Alphas <- simulate_alphas_HO_sep(lambdas_true,thetas_true,Alphas_0,Q_examinee,L,Jt)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
itempars_true <- array(runif(Jt*2*L,.1,.2), dim = c(Jt,2,L))
ETAs <- ETAmat(K, J, Q_matrix)
Y_sim <- simDINA(Alphas,itempars_true,ETAs,Test_order,Test_versions)

## -----------------------------------------------------------------------------
output_HMDCM = hmcdm(Y_sim,Q_matrix,"DINA_HO",Test_order,Test_versions,
                     chain_length=100,burn_in=30,
                     theta_propose = 2,deltas_propose = c(.45,.35,.25,.06))

output_HMDCM

summary(output_HMDCM)
a <- summary(output_HMDCM)
a$ss_EAP
a$lambdas_EAP
a$PPP_total_scores
a$PPP_item_means
a$PPP_item_ORs
mean(a$PPP_total_scores)
mean(upper.tri(a$PPP_item_ORs))
mean(a$PPP_item_means)

## -----------------------------------------------------------------------------
AAR_vec <- numeric(L)
for(t in 1:L){
  AAR_vec[t] <- mean(Alphas[,,t]==a$Alphas_est[,,t])
}
AAR_vec

## -----------------------------------------------------------------------------
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
library(bayesplot)
pp_check(output_HMDCM)
pp_check(output_HMDCM, plotfun="dens_overlay", type="item_mean")
pp_check(output_HMDCM, plotfun="hist", type="item_OR")
pp_check(output_HMDCM, plotfun="stat_2d", type="item_mean")
pp_check(output_HMDCM, plotfun="scatter_avg", type="total_score")
pp_check(output_HMDCM, plotfun="error_scatter_avg", type="total_score")

