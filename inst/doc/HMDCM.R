## ---- include = FALSE---------------------------------------------------------
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
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
for(i in 1:N){
 Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
thetas_true = rnorm(N)
lambdas_true = c(-1, 1.8, .277, .055)
Alphas <- sim_alphas(model="HO_sep", 
                    lambdas=lambdas_true, 
                    thetas=thetas_true, 
                    Q_matrix=Q_matrix, 
                    Design_array=Design_array)
table(rowSums(Alphas[,,5]) - rowSums(Alphas[,,1])) # used to see how much transition has taken place
itempars_true <- matrix(runif(J*2,.1,.2), ncol=2)

Y_sim <- sim_hmcdm(model="DINA",Alphas,Q_matrix,Design_array,
                   itempars=itempars_true)

## -----------------------------------------------------------------------------
output_HMDCM = hmcdm(Y_sim,Q_matrix,"DINA_HO",Test_order = Test_order, Test_versions = Test_versions,
                     chain_length=100,burn_in=30,
                     theta_propose = 2,deltas_propose = c(.45,.35,.25,.06))

output_HMDCM = hmcdm(Y_sim,Q_matrix,"DINA_HO",Design_array,
                     chain_length=100,burn_in=30,
                     theta_propose = 2,deltas_propose = c(.45,.35,.25,.06))

output_HMDCM

summary(output_HMDCM)
a <- summary(output_HMDCM)
a$ss_EAP
a$lambdas_EAP
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

## -----------------------------------------------------------------------------
# output_HMDCM1 = hmcdm(Y_sim, Q_matrix, "DINA_HO", Design_array,
#                      chain_length=100, burn_in=30,
#                      theta_propose = 2, deltas_propose = c(.45,.35,.25,.06))
# output_HMDCM2 = hmcdm(Y_sim, Q_matrix, "DINA_HO", Design_array,
#                      chain_length=100, burn_in=30,
#                      theta_propose = 2, deltas_propose = c(.45,.35,.25,.06))
# 
# library(coda)
# 
# x <- mcmc.list(mcmc(t(rbind(output_HMDCM1$ss, output_HMDCM1$gs, output_HMDCM1$lambdas))),
#                mcmc(t(rbind(output_HMDCM2$ss, output_HMDCM2$gs, output_HMDCM2$lambdas))))
# 
# gelman.diag(x, autoburnin=F)


