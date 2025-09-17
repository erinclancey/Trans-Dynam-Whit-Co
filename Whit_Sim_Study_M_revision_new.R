setwd("~/Trans-Dynam-Whit-Co-main")
source("Whit_Sim_Mech_Model_2_m.R")


############# Make an empty pomp object ############# 
Week <- seq(34, 52, 1)
Data <- as.data.frame(Week)
Data$reports_i <- rep(NA, length(Week))
Data$reports_j <- rep(NA, length(Week))

empty_data <- pomp(data=Data,
                  times="Week",t0=33,
                  rprocess=euler(SEIR.compart,delta.t=1/7),
                  rinit=SEIR_rinit,
                  rmeasure=rmeas,
                  accumvars=accumvars,
                  statenames=statenames,
                  paramnames=paramnames)
###############set up true  parameter grid#########
n=25
names=c("betai","betaj","betaij","k","m")
Y <- randomLHS(n, length(names))
Y[,1] <- qunif(Y[,1], 9.81, 12.24)
Y[,2] <- qunif(Y[,2], 7.68, 10.49)
Y[,3] <- qunif(Y[,3], 0.00, 0.35)
Y[,4] <- qunif(Y[,4], 1.45, 5.85)
Y[,5] <- qunif(Y[,5], 2.00, 18.00)
pars <- data.frame(Y)
colnames(pars) <- names
sigma=rep(7/3.59,n)
gamma=rep(7/3.56,n)
sigma=rep(7/3.59,n)
gamma=rep(7/3.56,n)
rho <- 0.76
Ni=rep(14254,n)
Ei0=rep(floor(79/rho),n)
Ii0=rep(floor(79/rho),n)
Ri0=rep(floor(0.041*14254),n)
Nj=rep(20785,n)
Ej0=rep(floor(2/rho),n)
Ij0=rep(floor(2/rho),n)
Rj0=rep(floor(0.041*20785),n)
pars <- cbind(pars, sigma, gamma,Ni, Ei0, Ii0, Ri0, Nj, Ej0, Ij0, Rj0)

############# pMCMC Setup ############# 

### Pick variances for the adaptive proposal
rw.var <- matrix(c(1,0,0,0,0,
                   0,1,0,0,0,
                   0,0,0.1,0,0,
                   0,0,0,0.5,0,
                   0,0,0,0,1), 
                 nrow = 5, dimnames = list(c("betai","betaj","betaij","k","m"), 
                                           c("betai","betaj","betaij","k","m")))


proposal <- mvn_rw_adaptive(
  rw.var=rw.var,
  scale.start = 1,
  scale.cooling = 0.999,
  shape.start = 500,
  target = 0.234,
  max.scaling = 50
)

### Set starting values for unknown parameters for n chains
chain=5
### Latin Hypercube Sampling
names=c("betai","betaj","betaij","k","m")

X <- randomLHS(chain, length(names))
X[,1] <- qunif(X[,1], 10,15)
X[,2] <- qunif(X[,2], 5,10)
X[,3] <- qunif(X[,3], 0.01,1)
X[,4] <- qunif(X[,4], 1,4)
X[,5] <- qunif(X[,5], 10,15)
params <- data.frame(X)
colnames(params) <- names
### Known parameter starting values
sigma=rep(7/3.59,chain)
gamma=rep(7/3.56,chain)
rho <- 0.76
Ni=rep(14254,chain)
Ei0=rep(floor(79/rho),chain)
Ii0=rep(floor(79/rho),chain)
Ri0=rep(floor(0.041*14254),chain)
Nj=rep(20785,chain)
Ej0=rep(floor(2/rho),chain)
Ij0=rep(floor(2/rho),chain)
Rj0=rep(floor(0.041*20785),chain)
### Create starting values data frame
theta.start <- data.frame(cbind(params, sigma, gamma, Ni, Ei0, Ii0, Ri0, Nj, Ej0, Ij0, Rj0))

sim_study <- data.frame(matrix(ncol = 15, nrow = 0))
x <- c("betai_mode","betaj_mode","betaij_mode","k_mode","m_mode",
       "betai_low","betaj_low","betaij_low","k_low","m_low",
       "betai_hi","betaj_hi","betaij_hi","k_hi","m_hi")
colnames(sim_study) <- x

##########################Run Sim Study
for (i in 1:nrow(pars)){
##########Simulate data
empty_data %>%
  simulate(
    params=pars[i,],
    nsim=1,format="data.frame",include.data=FALSE) %>% 
  select(Week,reports_i, reports_j) -> sim_cases

sim_cases %>%
  pomp(
    times="Week",t0=33,
    rprocess=euler(SEIR.compart,delta.t=1/7),
    rinit=SEIR_rinit,
    rmeasure=rmeas,             
    dmeasure=dmeas,              
    accumvars=accumvars,
    statenames=statenames,
    paramnames=paramnames
  ) -> SEIR

##########Run pMCMC#######
M=10000 # the number of mcmc iterations to run
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE) %dopar% {
  library(pomp)
  library(magrittr)
  SEIR %>% pmcmc(Nmcmc=M,
                 proposal=proposal,
                 Np=1000,
                 params=theta.start
  ) -> pmcmc
  results <- as.data.frame(traces(pmcmc))
} -> results_pmcmc

list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
chain.no <- seq(1:chain)
iter.no <- seq(1, M+1, by=1)

for(k in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[k]]$chain <- rep(chain.no[k],nrow(list_results_pmcmc[[k]]))
  list_results_pmcmc[[k]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)
write.csv(posterior, file= paste('result_',i,".csv",sep = ""))

post <- posterior %>% 
  dplyr::select(betai,betaj,betaij,k,m, chain, iter)
post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())
for(j in seq_along(post.list)){
  mcmc.list[[j]] <- mcmc(post.list[[j]])
}
processed <- window(mcmc.list, start=100, end=M+1, thin=50) 
processed <- data.frame(do.call(rbind, processed))

betai_mode=as.vector(posterior.mode(mcmc(processed$betai), adjust=1))
betaj_mode=as.vector(posterior.mode(mcmc(processed$betaj), adjust=1))
betaij_mode=as.vector(posterior.mode(mcmc(processed$betaij), adjust=1))
k_mode=as.vector(posterior.mode(mcmc(processed$k), adjust=1))
m_mode=as.vector(posterior.mode(mcmc(processed$m), adjust=1))

betai_low=ci(processed$betai, ci=0.95, method = "HDI")$CI_low
betai_hi=ci(processed$betai, ci=0.95, method = "HDI")$CI_high
betaj_low=ci(processed$betaj, ci=0.95, method = "HDI")$CI_low
betaj_hi=ci(processed$betaj, ci=0.95, method = "HDI")$CI_high
betaij_low=ci(processed$betaij, ci=0.95, method = "HDI")$CI_low
betaij_hi=ci(processed$betaij, ci=0.95, method = "HDI")$CI_high
k_low=ci(processed$k, ci=0.95, method = "HDI")$CI_low
k_hi=ci(processed$k, ci=0.95, method = "HDI")$CI_high
m_low=ci(processed$m, ci=0.95, method = "HDI")$CI_low
m_hi=ci(processed$m, ci=0.95, method = "HDI")$CI_high

result <- data.frame(cbind(betai_mode,betaj_mode,betaij_mode,k_mode,m_mode,
                          betai_low,betaj_low,betaij_low,k_low,m_low,
                          betai_hi,betaj_hi,betaij_hi,k_hi,m_hi))
sim_study[i,] <- result
print(i)

}
sim_study_df <- cbind(sim_study, pars)

sim_study_df <- read.csv(file="sim_study_df.csv", header=TRUE)


umod <- summary(lm(betai_mode ~ betai , data = sim_study_df))
cmod <- summary(lm(betaj_mode ~ betaj , data = sim_study_df))
mixmod <- summary(lm(betaij_mode ~ betaij , data = sim_study_df))
kmod <- summary(lm(k_mode ~ k , data = sim_study_df))
mmod <- summary(lm(m_mode ~ m , data = sim_study_df))
umod
cmod
mixmod
kmod
mmod

round(c(umod$coefficients[2,1]-umod$coefficients[2,2]*1.96,umod$coefficients[2,1]+umod$coefficients[2,2]*1.96),2)
round(c(cmod$coefficients[2,1]-cmod$coefficients[2,2]*1.96,cmod$coefficients[2,1]+cmod$coefficients[2,2]*1.96),2)
round(c(mixmod$coefficients[2,1]-mixmod$coefficients[2,2]*1.96,mixmod$coefficients[2,1]+mixmod$coefficients[2,2]*1.96),2)
round(c(kmod$coefficients[2,1]-kmod$coefficients[2,2]*1.96,kmod$coefficients[2,1]+kmod$coefficients[2,2]*1.96),2)
round(c(mmod$coefficients[2,1]-mmod$coefficients[2,2]*1.96,mmod$coefficients[2,1]+mmod$coefficients[2,2]*1.96),2)

round(c(umod$coefficients[1,1]-umod$coefficients[2,1]*1.96,umod$coefficients[1,1]+umod$coefficients[2,1]*1.96),2)/100000
round(c(cmod$coefficients[1,1]-cmod$coefficients[2,1]*1.96,cmod$coefficients[1,1]+cmod$coefficients[2,1]*1.96),2)/100000
round(c(mixmod$coefficients[1,1]-mixmod$coefficients[2,1]*1.96,mixmod$coefficients[2,1]+mixmod$coefficients[2,2]*1.96),2)/100000
round(c(kmod$coefficients[1,1]-kmod$coefficients[2,1]*1.96,kmod$coefficients[1,1]+kmod$coefficients[2,1]*1.96),2)
round(c(mmod$coefficients[1,1]-mmod$coefficients[2,1]*1.96,mmod$coefficients[1,1]+mmod$coefficients[2,1]*1.96),2)

p1 <- ggplot(sim_study_df, aes(x=betai/100000, y=betai_mode/100000)) +
  geom_point(color="#CC79A7", shape=20, size=3) +
  xlab(TeX("$\\beta_u$")) + ylab(TeX("$\\hat{\\beta}_u$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = umod$coefficients[1,1]/100000, slope = umod$coefficients[2,1], 
              linetype=1, color="#CC79A7", size=1)+
  xlim(c(0.000075,0.00014))+ylim(c(0.000075,0.00014))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )


p2 <- ggplot(sim_study_df, aes(x=betaj/100000, y=betaj_mode/100000)) +
  geom_point( color="#0072B2",shape=20, size=3) +
  xlab(TeX("$\\beta_c$")) + ylab(TeX("$\\hat{\\beta}_c$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = cmod$coefficients[1,1]/100000, slope = cmod$coefficients[2,1], 
              linetype=1, color="#0072B2", size=1)+
  xlim(c(0.00001,0.00013))+ylim(c(0.00001,0.00013))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p3 <- ggplot(sim_study_df, aes(x=betaij/100000, y=betaij_mode/100000)) +
  geom_point( color="#009E73",shape=20, size=3) +
  xlab(TeX("$\\beta_m$")) + ylab(TeX("$\\hat{\\beta}_m$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = mixmod$coefficients[1,1]/100000, slope = mixmod$coefficients[2,1], 
              linetype=1, color="#009E73", size=1)+
  xlim(c(0,0.000008))+ylim(c(0,0.000008))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p4 <- ggplot(sim_study_df, aes(x=k, y=k_mode)) +
  geom_point(color="#999999", shape=20, size=3) +
  xlab(TeX("$k$")) + ylab(TeX("$\\hat{k}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = kmod$coefficients[1,1], slope = kmod$coefficients[2,1], 
              linetype=1, color="#999999", size=1)+
  xlim(c(0,15))+ylim(c(0,15))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p5 <- ggplot(sim_study_df, aes(x=m, y=m_mode)) +
  geom_point(color="#D55E00", shape=20, size=3) +
  xlab(TeX("$m$")) + ylab(TeX("$\\hat{m}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = mmod$coefficients[1,1], slope = mmod$coefficients[2,1], 
              linetype=1, color="#D55E00", size=1)+
  xlim(c(0,30))+ylim(c(0,30))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

ggarrange(p1,p2,p3,p4,p5,
          ncol = 3, nrow = 2)



betaiCI <- vector()
betajCI <- vector()
betaijCI <- vector()
kCI <- vector()
mCI <- vector()
for (i in 1:nrow(sim_study_df)){
  if(sim_study_df$betai[i] >=sim_study_df$betai_low[i] && sim_study_df$betai[i]<=sim_study_df$betai_hi[i]){
    betaiCI[i]=1} else{
      betaiCI[i]=0
    }
  if(sim_study_df$betaj[i]>=sim_study_df$betaj_low[i] && sim_study_df$betaj[i]<=sim_study_df$betaj_hi[i]){
    betajCI[i]=1} else{
      betajCI[i]=0
    }
  if(sim_study_df$betaij[i]>=sim_study_df$betaij_low[i] && sim_study_df$betaij[i]<=sim_study_df$betaij_hi[i]){
    betaijCI[i]=1} else{
      betaijCI[i]=0
    }
  if(sim_study_df$k[i]>=sim_study_df$k_low[i] && sim_study_df$k[i]<=sim_study_df$k_hi[i]){
    kCI[i]=1} else{
      kCI[i]=0
    }
  if(sim_study_df$m[i]>=sim_study_df$m_low[i] && sim_study_df$m[i]<=sim_study_df$m_hi[i]){
    mCI[i]=1} else{
      mCI[i]=0
    }
}

sum(betaiCI)
sum(betajCI)
sum(betaijCI)
sum(kCI)
sum(mCI)


