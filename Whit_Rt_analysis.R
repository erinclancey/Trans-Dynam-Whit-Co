library(tidyverse)
library(cowplot)
library(EpiEstim)
library(extrafont)
library(latex2exp)

setwd("~/Trans-Dynam-Whit-Co")

###parameter values for latency and recovery and generation interval
alpha=1/3.59
gamma=1/3.56

GI_Whit=1/gamma + 1/alpha 
sd_GI_Whit=sqrt(1/gamma^2 + 1/alpha^2)

J_date <- seq(237,362,1)

Whit_Co_weekly <- read.csv("Whit_CO_week.csv", header=TRUE)
Whit_Co_weekly <- subset(Whit_Co_weekly, Week>34)
C <- ggplot(Whit_Co_weekly, aes(x=Week,y=Cases)) + theme_minimal()+ 
  geom_bar(stat="identity", color="darkolivegreen4", fill="darkolivegreen4", alpha=0.25)+
  ylab("Case Reports")+xlab("Week of Year")+
  ggtitle("Whitman Co. Weekly Case Reports")+
  scale_x_continuous(n.break=8)+
  scale_y_continuous(n.breaks=10)
plot(C)

agg_R_Whit <- estimate_R(incid = Whit_Co_weekly$Cases,
                         dt = 7L,
                         dt_out = 7L,
                         recon_opt = "naive",
                         iter = 10L,
                         tol = 1e-6,
                         grid = list(precision = 0.001, min = -1, max = 1),
                         method="parametric_si",
                         config = make_config(list(
                           mean_si = GI_Whit, 
                           std_si = sd_GI_Whit)))
agg_R_Whit_df <- agg_R_Whit$R
#agg_R_Whit_df <- agg_R_Whit$R[-c(1:7),]
agg_R_Whit_df$J_date <- J_date[-c(1:13)]
filter(agg_R_Whit_df, `Mean(R)` <1)[1,1]

D <- ggplot(agg_R_Whit_df , aes(x=J_date, y=`Mean(R)`)) + theme_minimal()+
  xlab("Day of Year")+geom_line(color="darkolivegreen4")+
  ylab(TeX("\\textit{R}$_t$"))+
  ggtitle("Whitman Co. Total Population Reproductive Numbers")+
  scale_x_continuous(limits=c(250,365), n.breaks=9)+
  scale_y_continuous(limits=c(0,10), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")
plot(D)

##############Separate Week

E <- ggplot(Whit_Co_weekly, aes(x=Week,y=PUL_cases)) + theme_minimal()+ 
  geom_bar(stat="identity", color="steelblue",fill="steelblue", alpha=0.6)+
  geom_bar(aes(x=Week,y=WSU_cases),stat="identity",color="darkred", fill="darkred", alpha=0.3)+
  ylab("Case Reports")+xlab("Week of Year")+
  ggtitle("Whitman Co. Subpopulation Weekly Case Reports")+
  scale_x_continuous(n.breaks=8)+
  scale_y_continuous(n.breaks=10)
plot(E)


WSU_R_t <- estimate_R(incid = Whit_Co_weekly$WSU_cases,
                         dt = 7L,
                         dt_out = 7L,
                         recon_opt = "naive",
                         iter = 10L,
                         tol = 1e-6,
                         grid = list(precision = 0.001, min = -1, max = 1),
                         method="parametric_si",
                         config = make_config(list(
                           mean_si = GI_Whit, 
                           std_si = sd_GI_Whit)))
WSU_R_t_df <- WSU_R_t$R
WSU_R_t_df$J_date <- J_date[-c(1:13)]
filter(WSU_R_t_df, `Mean(R)` <1)[1,1]

PUL_R_t <- estimate_R(incid = Whit_Co_weekly$PUL_cases,
                      dt = 7L,
                      dt_out = 7L,
                      recon_opt = "naive",
                      iter = 10L,
                      tol = 1e-6,
                      grid = list(precision = 0.001, min = -1, max = 1),
                      method="parametric_si",
                      config = make_config(list(
                        mean_si = GI_Whit, 
                        std_si = sd_GI_Whit)))
PUL_R_t_df <- PUL_R_t$R
PUL_R_t_df$J_date <- J_date[-c(1:13)]
filter(PUL_R_t_df, `Mean(R)` <1)[1,1]

G <- ggplot()+
  geom_line(data = WSU_R_t_df, aes(x=J_date, y=`Mean(R)`), color = "darkred") +
  geom_line(data = PUL_R_t_df, aes(x=J_date, y=`Mean(R)`), color = "steelblue") +
  xlab("Day of Year")+ theme_minimal()+
  ylab(TeX("\\textit{R}$_t$"))+
  ggtitle("Whitman Co. Subpopuation Reproductive Numbers")+
  scale_x_continuous(limits=c(250,365), n.breaks=9)+
  scale_y_continuous(limits=c(0,10), n.breaks=5)+
  geom_ribbon(data = WSU_R_t_df,aes(x=J_date, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_ribbon(data = PUL_R_t_df,aes(x=J_date, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")
plot(G)

plot3 <- plot_grid(C,D,E,G, ncol = 2, nrow =2 , rel_heights=c(1,1,1,1,1,1), labels = c('A','', 'B'))
plot(plot3)
