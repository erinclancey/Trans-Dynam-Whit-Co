setwd("~/Covid_Mat/EC MS/Trans-Dynam-Whit-Co-main_revision/Final code and figures for revision")
source("Whit_Sim_Mech_Model_2_m.R")

# Read-in Joint Posterior Distribution
#posterior <- read.csv(file="posterior_revision.csv", header=TRUE)
# Split the data frame into two parts
# posterior_first_half <- posterior[1:(floor(nrow(posterior)/2)+1), ]
# posterior_second_half <- posterior[floor((nrow(posterior)/2) + 2):nrow(posterior), ]

posterior1 <- read.csv(file="posterior_migrate_first_half.csv", header=TRUE)
posterior2 <- read.csv(file="posterior_migrate_second_half.csv", header=TRUE)
posterior <- rbind(posterior1,posterior2)

posterior$m <- round(posterior$m)
M=50000

### Post-process the chains
post <- posterior %>% 
  dplyr::select(betai,betaj,betaij,k,m, chain, iter)
post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())
### Fills the list with MCMC objects
for(i in seq_along(post.list)){
  mcmc.list[[i]] <- mcmc(post.list[[i]])
}

### Perform Diagnistic Tests for convergence and burn-in
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
geweke.diag(mcmc.list, frac1=0.5, frac2=0.25)
gelman.diag(mcmc.list, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)

# Review the autocorrelations for each parameter
acf(post.list[[1]]$betai, lag.max = 100)
acf(post.list[[1]]$betaj, lag.max = 100) 
acf(post.list[[1]]$betaij, lag.max = 100)
acf(post.list[[1]]$k, lag.max = 100)
acf(post.list[[1]]$m, lag.max = 100)

#Post-process the chains
processed <- window(mcmc.list, start=30, end=M+1, thin=80) 
processed <- data.frame(do.call(rbind, processed))
processed.long <- melt(processed, id=c("chain", "iter"))
nrow(processed)

# Re-scale transmission parameters for plotting
post$betai <- post$betai/100000
post$betaj <- post$betaj/100000
post$betaij <- post$betaij/100000
chains <- post %>% mutate_at("chain", as.character)
chains.long <- melt(chains, id=c("chain", "iter"))
### Plot the five unprocessed pMCMC chains
# Quick plot
### Plot the five unprocessed pMCMC chains
plot_names <- as_labeller(c('betai' = "paste(beta)[u]",
                            'betaj' = "paste(beta)[c]",
                            'betaij' = "paste(beta)[m]",
                            'k' = "paste(theta)",
                            'm' = "paste(m(t))"),label_parsed)
upc <- ggplot(chains.long, aes(x = iter, y = value, group=chain )) + 
  geom_line(aes(color=chain)) + 
  theme_minimal() +
  scale_color_manual(values = c("#F0E442","#000000", "#009E73","#999999","#0072B2")) + 
  facet_wrap(vars(variable),labeller = plot_names, scales='free', ncol = 3)+
  labs(x="Iteration",y="Parameter Value")+ 
  theme(strip.text = element_text(size=20), axis.title = element_text(size = 18),panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
plot(upc)

### Calculate mode values
modes <- processed.long %>% 
  group_by(variable) %>% 
  summarise(mode = posterior.mode(mcmc(value), adjust=1))
modes <- data.frame(modes)
modes

### Calculate Credible Intervals for Parameters
ci(processed$betai, ci=0.95, method = "HDI") 
ci(processed$betaj, ci=0.95, method = "HDI")
ci(processed$betaij, ci=0.95, method = "HDI")
ci(processed$k, ci=0.95, method = "HDI")
ci(processed$m, ci=0.95, method = "HDI")


# Calculate posterior for R0 to find modes and CI's
R0 <- R_t_func(Xi=14254, Xj=20785,
               betai=processed$betai/100000,betaj=processed$betaj/100000,betaij=processed$betaij/100000,
               gamma=7/3.56)
R0_i <- (14254*processed$betai/100000)/(7/3.56)
R0_j <- (20785*processed$betaj/100000)/(7/3.56)

R0 <- data.frame(R0)
R0_df <- data.frame(R0_i,R0_j,R0)
R0_df.long <- melt(R0_df)
R0_mode <- R0_df.long %>% 
  group_by(variable) %>% 
  summarise(mode = posterior.mode(mcmc(value), adjust=1))
R0_mode <- data.frame(R0_mode)
R0_mode

### Calculate Credible Intervals for R0's
ci(R0_i, ci=0.95, method = "HDI")
ci(R0_j, ci=0.95, method = "HDI")
ci(R0, ci=0.95, method = "HDI")


processed.long.2 <- subset(processed.long, variable !="k" & variable !="m")
mode <- modes[-c(4,5),]
mode$mode <- mode$mode/100000
processed.long.2$value <- processed.long.2$value/100000

P <- ggplot(processed.long.2, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=60, colour="#0072B2",fill="#0072B2")+
  geom_density(alpha=.2, fill="#0072B2", color="#0072B2", adjust=1)+
  geom_vline(data=mode, aes(xintercept=mode), color="black", linewidth=0.75, linetype=2)+
  facet_wrap(vars(variable),labeller = as_labeller(plot_names), scales='free_y', ncol = 1)+labs(x="Transmission Rate",y="Density")+
  scale_x_continuous(n.breaks = 6)+scale_y_continuous(n.breaks = 5)+
  theme(strip.text = element_text(size=20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 18),  
        panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))

R0_names <- as_labeller(c('R0' = "paste(R)[0[total]]",
                            'R0_i' = "paste(R)[0[u]]",
                            'R0_j' = "paste(R)[0[c]]"),label_parsed)
Q <- ggplot(R0_df.long, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=40, colour="black",fill="black")+
  geom_density(alpha=.2, colour="black",fill="black", adjust=1)+
  geom_vline(data=R0_mode, aes(xintercept=mode), color="black", linewidth=0.75, linetype=2)+
  geom_vline(aes(xintercept=1), color="darkred", linewidth=0.75, linetype=3,alpha=.75)+
  facet_wrap(vars(variable),labeller = as_labeller(R0_names), scales='free_y', ncol = 1)+labs(x="Reproductive Number",y="")+
  scale_x_continuous(n.breaks = 6)+scale_y_continuous(n.breaks = 5)+
  theme(strip.text = element_text(size=20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 18), 
        panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))

plot <- plot_grid(P,Q, ncol = 2, rel_widths=c(1,1))
plot(plot)


plot_names.2 <- as_labeller(c('k' = "paste(theta)",'m' = "paste(m(t))"),label_parsed)
mode.2 <- modes[-c(1,2,3),]
processed.long.3 <- subset(processed.long, variable =="k" | variable =="m")
P.2 <- ggplot(processed.long.3, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=35, colour="#0072B2",fill="#0072B2")+
  geom_density(alpha=.2, fill="#0072B2", color="#0072B2", adjust=1)+
  geom_vline(data=mode.2, aes(xintercept=mode), color="black", linewidth=0.75, linetype=2)+
  facet_wrap(vars(variable),labeller = as_labeller(plot_names.2), scales='free_y', ncol = 2)+labs(x="Parameter Value",y="Density")+
  scale_x_continuous(n.breaks = 8)+scale_y_continuous(n.breaks = 8)+
  theme(strip.text = element_text(size=20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 18), 
        panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
plot(P.2)

