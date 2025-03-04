setwd("~/Covid_Mat/EC MS/Trans-Dynam-Whit-Co-main_revision")
source("Whit_Sim_Mech_Model.R")


posterior <- read.csv(file="posterior_revision.csv", header=TRUE)
post <- posterior %>% 
  dplyr::select(betai,betaj,betaij,k)
head(post)
nrow(post)
n=1000
sample_post <- post[sample(nrow(post), n), ]
head(sample_post)

sigma=rep(7/3.59,n)
gamma=rep(7/3.56,n)
rho <-rep(0.76,n) 
Ni=rep(14254,n)
Ei0=rep(floor(79/rho),n)
Ii0=rep(floor(79/rho),n)
Ri0=rep(floor(0.041*14254),n)
Nj=rep(20785,n)
Ej0=rep(floor(2/rho),n)
Ij0=rep(floor(2/rho),n)
Rj0=rep(floor(0.041*20785),n)

pars <- data.frame(cbind(sample_post, sigma, gamma, rho, Ni, Ei0, Ii0, Ri0, Nj, Ej0, Ij0, Rj0))

simdf <- data.frame(matrix(ncol = 15, nrow = 0))
x <- c("Week", ".id","Si","Sj","Ei","Ej","Ii","Ij","Ri","Rj","Hi","Hj","reports_i","reports_j","simID")
colnames(simdf) <- x
for(i in 1:n){
SEIR %>%
  simulate(
    params=pars[i,],
    nsim=1,format="data.frame",include.data=FALSE) -> sims
  simID <- rep(i,nrow(sims))
  sims <- cbind(sims, simID)
  simdf <- rbind(simdf, sims)
}

A <- ggplot() + theme_minimal() + 
  geom_line(data=simdf, aes(x=Week,y=reports_i, group=simID), color="grey20",alpha=0.1) + 
  geom_line(data=WSU_PUL, aes(x=Week, y=reports_i), color="tomato",size=0.75)+
  ylab("Case Reports")+xlab("Week of Year")+
  ggtitle("University Student Subpopulation")+
  scale_x_continuous(n.breaks=12)+
  scale_y_continuous(n.breaks=7, limits=c(0,500))


B <- ggplot() + theme_minimal() +
  geom_line(data=simdf, aes(x=Week,y=reports_j, group=simID), color="grey20",alpha=0.1) + 
  geom_line(data=WSU_PUL, aes(x=Week, y=reports_j), color="turquoise3", size=0.75)+
  ylab("Case Reports")+xlab("Week of Year")+
  ggtitle("Community Subpopulation")+
  scale_x_continuous(n.breaks=12)+
  scale_y_continuous(n.breaks=5, limits=c(0,300))

plot <- plot_grid(A,B, ncol = 1, nrow =2 , rel_heights=c(1,1), labels = c('A', 'B'))
plot(plot)
