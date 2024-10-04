setwd("~/Trans-Dynam-Whit-Co")
source("Whit_Sim_Mech_Model.R")

############# pMCMC Setup ############# 

### Pick variances for the proposal

standard.devs <-c(1, 1, 0.1, 0.5)
rho <- 0.76
params.est <- c("betai","betaj","betaij","k")
rw.sd <- setNames(standard.devs, params.est)
# proposal for the first iteration
# proposal <- mvn_diag_rw(rw.sd)

# Coviance matrix calculated after the initial run
rw.var <- matrix(c(0.279665318,  0.002480186, -0.01192454,  0.010768938,
0.002480186,  0.727658376, -0.09595934,  0.006475304,
-0.011924538, -0.095959344,  0.03121689, -0.018194206,
0.010768938,  0.006475304, -0.01819421,  0.769467891), nrow = 4, dimnames = list(c("betai","betaj","betaij","k"), c("betai","betaj","betaij","k")))

proposal <- mvn_rw_adaptive(
  rw.var=rw.var,
  scale.start = 1,
  scale.cooling = 0.999,
  shape.start = 500,
  target = 0.234,
  max.scaling = 50
)

### Set starting values for unknown parameters for n chains
n=5
### Latin Hypercube Sampling
names=c("betai","betaj","betaij","k")

X <- randomLHS(n, length(names))
X[,1] <- qunif(X[,1], 10,15)
X[,2] <- qunif(X[,2], 5,10)
X[,3] <- qunif(X[,3], 0.01,1)
X[,4] <- qunif(X[,4], 1,4)

params <- data.frame(X)
colnames(params) <- names


### Known parameter starting values

sigma=rep(7/3.59,n)
gamma=rep(7/3.56,n)

Ni=rep(14254,n)
Ei0=rep(floor(79/rho),n)
Ii0=rep(floor(79/rho),n)
Ri0=rep(floor(0.041*14254),n)
Nj=rep(20785,n)
Ej0=rep(floor(2/rho),n)
Ij0=rep(floor(2/rho),n)
Rj0=rep(floor(0.041*20785),n)


### Create starting values data frame

theta.start <- data.frame(cbind(params, sigma, gamma, Ni, Ei0, Ii0, Ri0, Nj, Ej0, Ij0, Rj0))
### Run pMCMC

M=10000 # the number of mcmc iterations to run
start_time <- Sys.time()
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE) %dopar% {
  library(pomp)
  library(magrittr)
  set.seed(320)
  SEIR %>% pmcmc(Nmcmc=M,
                 proposal=proposal,
                 Np=2000,
                 params=theta.start
  ) -> pmcmc
  results <- as.data.frame(traces(pmcmc))
} -> results_pmcmc
end_time <- Sys.time()
end_time - start_time    # capture computational time

list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
chain.no <- seq(1:n)
iter.no <- seq(1, M+1, by=1)

for(i in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[i]]$chain <- rep(chain.no[i],nrow(list_results_pmcmc[[i]]))
  list_results_pmcmc[[i]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)
### Export pMCMC Results
#write.csv(posterior, file="posterior_end_rhofixed.csv")

##### OR calculate the covariance matrix
post <- posterior %>% 
  dplyr::select(betai,betaj,betaij, k, chain, iter)
post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())
### Fills the list with MCMC objects
for(i in seq_along(post.list)){
  mcmc.list[[i]] <- mcmc(post.list[[i]])
}
processed <- window(mcmc.list, start=50, end=M+1, thin=50) 
processed <- data.frame(do.call(rbind, processed))
# calculate covariance matrix after initial run
# df <- processed[,-c(5,6)]
# rw.var=cov(df)
