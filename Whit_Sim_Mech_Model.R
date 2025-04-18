###########Load Packages
library(pomp)
library(tidyverse)
library(cowplot)
library(EpiEstim)
library(reshape2)
library(latex2exp)
library(moments)
library(mousetrap)


#####
library(foreach)
library(iterators)
library(parallel)
library(rngtools)
library(doParallel)
library(doRNG)
registerDoParallel()
registerDoRNG(2488820)
library(ggpubr)
library(bayestestR) 
library(coda)
library(MCMCglmm) 
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(tidyr)
library(tidyselect)
library(tibble)
library(readr)
library(stringr)
library(forcats)
library(lhs)
set.seed(320)

# dir.create("tmp")
# options(pomp_cdir="./tmp")



#############################Code Model in Pomp
SEIR.compart <- Csnippet("
  double dN_SiEi = rbinom(Si,1-exp((-betai/100000*Ii-betaij/100000*Ij)*dt));
  double dN_SjEj = rbinom(Sj,1-exp((-betaj/100000*Ij-betaij/100000*Ii)*dt));

  double dN_EiIi = rbinom(Ei,1-exp(-sigma*dt));
  double dN_EjIj = rbinom(Ej,1-exp(-sigma*dt));
  
  double dN_IiRi = rbinom(Ii,1-exp(-gamma*dt));
  double dN_IjRj = rbinom(Ij,1-exp(-gamma*dt));
  
  Si -= dN_SiEi;
  Sj -= dN_SjEj;
  
  Ei += dN_SiEi - dN_EiIi;
  Ej += dN_SjEj - dN_EjIj;
  
  Ii += dN_EiIi - dN_IiRi;
  Ij += dN_EjIj - dN_IjRj;
  
  Ri += dN_IiRi;
  Rj += dN_IjRj;
  
  Hi += dN_EiIi;
  Hj += dN_EjIj;
")

SEIR_rinit <- Csnippet("
  Si = Ni-Ei0-Ii0-Ri0;
  Sj = Nj-Ej0-Ij0-Rj0;
  Ei = Ei0;
  Ej = Ej0;
  Ii = Ii0;
  Ij = Ij0;
  Ri = Ri0;
  Rj = Rj0;
  Hi = Ii;
  Hj = Ij;
")

rmeas <- Csnippet("
  reports_i = rnbinom_mu(k,0.76*Hi);
  reports_j = rnbinom_mu(k,0.76*Hj);
  ")

dmeas <- Csnippet("if (betai<0||betaj<0||betaij<0||k<0) {
  lik = (give_log) ? R_NegInf : 0.0;
} else {
  lik=dnbinom_mu(reports_i,k,(0.76*Hi),give_log)+dnbinom_mu(reports_j,k,(0.76*Hj),give_log);
}")


accumvars=c("Hi","Hj")
statenames=c("Si","Sj","Ei","Ej","Ii","Ij","Ri","Rj","Hi","Hj")
paramnames=c("betai","betaj","betaij","sigma","gamma","k",
             "Ni","Ei0","Ii0","Ri0","Nj","Ej0","Ij0","Rj0")


############# Weekly Case Data - THE REAL DATA ############# 
# NOTE: population i = WSU students; population j = community members
read.csv("Whit_CO_week.csv", header=TRUE) %>%
  select(Week,reports_i=WSU_cases, reports_j=PUL_cases) -> WSU_PUL
WSU_PUL %>% as.data.frame() %>% head()
cases.plot <- ggplot(WSU_PUL, aes(x = Week)) + theme_minimal()+ ylab("Case Reports") + xlab("Week")+
  geom_line(aes(y = reports_i, color = 'reports_i')) + 
  geom_line(aes(y = reports_j, color = 'reports_j')) + guides(colour="none") + ggtitle("The Data")

plot(cases.plot)


############# Input Real Data into a pomp object ############# 

WSU_PUL %>%
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

############ Function to calculate R0
R_t_func <- function(Xi,Xj,betai,betaj,betaij,gamma){
  (Xi*betai + Xj*betaj + 
     sqrt(Xi^2*betai^2 + 4*Xi*Xj*betaij^2 - 2*Xi*Xj*betai*betaj + Xj^2*betaj^2))/
    (2*gamma)
}
