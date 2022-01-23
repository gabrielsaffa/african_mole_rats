#########################
### African mole-rats ###


setwd ("C:/R_folder/african mole rats")

library (phytools)
library (rethinking)

# setting recommended by Stan

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
memory.limit(size=2e6)


###############################################
### varying intercept, varying slope models ###

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# log-transform the data

log_IFA <- log (data$IFA)
log_BM <- log (data$BM)

d <- list(
   log_IFA=log_IFA,
   log_BM=log_BM,
   species=species
)

# IFA

m_IFA_ml <- ulam(
   alist(
      # Gaussian model for the outcome
      log_IFA ~ normal (mu_y, sigma_y), 
      
      # linear model for the mean
      mu_y <- a[species] + bBM[species]*log_BM,
      
      # adaptive priors 
      c(a,bBM)[species] ~ multi_normal (c(mu_a,mu_b), Rho, sigma_ab),
      
      # (hyper)priors
      sigma_y ~ exponential (1),
      mu_a ~ normal (0,0.5),
      mu_b ~ normal (0,0.5),
      sigma_ab ~ exponential (1),
      Rho ~ lkj_corr (2)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_IFA_ml_prec <- precis (m_IFA_ml,2)
m_IFA_ml_prec <- data.frame (round (m_IFA_ml_prec,2))
write.csv (data.frame (m_IFA_ml_prec), file="m_IFA_ml_prec.csv")

m_IFA_ml_post <- extract.samples (m_IFA_ml) # get posterior samples
write.csv (data.frame (m_IFA_ml_post), file="m_IFA_ml_post.csv")


### load the data again and drop rows with missing values and naked mole rats

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# log-transform the data

log_RDP <- log (data$RDP)
log_TJI <- log (data$TJI)
log_BM <- log (data$BM)

# create a data list

d <- list(
   log_RDP=log_RDP,
   log_TJI=log_TJI,
   log_BM=log_BM,
   species=species
)

# RDP

m_RDP_ml <- ulam(
   alist(
      # Gaussian model for the outcome
      log_RDP ~ normal (mu_y, sigma_y), 
      
      # linear model for the mean
      mu_y <- a[species] + bBM[species]*log_BM,
      
      # adaptive priors 
      c(a,bBM)[species] ~ multi_normal (c(mu_a,mu_b), Rho, sigma_ab),

      # (hyper)priors
      sigma_y ~ exponential (1),
      mu_a ~ normal (0,0.5),
      mu_b ~ normal (0,0.5),
      sigma_ab ~ exponential (1),
      Rho ~ lkj_corr (2)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_RDP_ml_prec <- precis (m_RDP_ml,2)
m_RDP_ml_prec <- data.frame (round (m_RDP_ml_prec,2))
write.csv (data.frame (m_RDP_ml_prec), file="m_RDP_ml_prec.csv")

m_RDP_ml_post <- extract.samples (m_RDP_ml) # get posterior samples
write.csv (data.frame (m_RDP_ml_post), file="m_RDP_ml_post.csv")


# TJI

m_TJI_ml <- ulam(
   alist(
      # Gaussian model for the outcome
      log_TJI ~ normal (mu_y, sigma_y), 
      
      # linear model for the mean
      mu_y <- a[species] + bBM[species]*log_BM,
      
      # adaptive priors 
      c(a,bBM)[species] ~ multi_normal (c(mu_a,mu_b), Rho, sigma_ab),
      
      # (hyper)priors
      sigma_y ~ exponential (1),
      mu_a ~ normal (0,0.5),
      mu_b ~ normal (0,0.5),
      sigma_ab ~ exponential (1),
      Rho ~ lkj_corr (2)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_TJI_ml_prec <- precis (m_TJI_ml,2)
m_TJI_ml_prec <- data.frame (round (m_TJI_ml_prec,2))
write.csv (data.frame (m_TJI_ml_prec), file="m_TJI_ml_prec.csv")

m_TJI_ml_post <- extract.samples (m_TJI_ml) # get posterior samples
write.csv (data.frame (m_TJI_ml_post), file="m_TJI_ml_post.csv")


#####################################################################################
### species-level regressions with measurement error on the outcome and predictor ###

tree <- read.nexus ("mole-rat_tree.nex") # read tree without out-groups
tree <- ls.consensus (tree) # create a consensus tree 
tree <- drop.tip (tree, tip=c("Bathyergus_janetta")) # drop B. janetta
spp_i <- c("Heterocephalus_glaber","Heliophobius_argenteocinereus", "Fukomys_mechowi", "Fukomys_damarensis", "Cryptomys_hottentotus", "Georychus_capensis", "Bathyergus_suillus")
tree_trimmed <- keep.tip (tree, spp_i) # combine tips of the phylogeny with societies in the sample
Dmat <- cophenetic (tree_trimmed) # compute distance matrix
Dmat <- Dmat[spp_i,spp_i]/max(Dmat) # normalize the distances to scale from 0 to 1 and include it to the data list

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# log-transform the data

log_IFA <- log (data$IFA)
log_BM <- log (data$BM)

# compute species means and standard errors for each variable
# define standard error function

se_fun <- function (x) {sd(x)/sqrt(length(x))}

IFA_mu <- tapply (log_IFA, data$species, mean)
IFA_se <- tapply (log_IFA, data$species, se_fun)

BM_mu <- tapply (log_BM, data$species, mean)
BM_se <- tapply (log_BM, data$species, se_fun)


# create a data list

d <- list(
   IFA_mu=IFA_mu,
   BM_mu=BM_mu,
   IFA_se=IFA_se,
   BM_se=BM_se,
   N=nrow(IFA_mu),
   Dmat=Dmat
)


### IFA
### measurement error model without phylogeny

m_IFA_me <- ulam(
   alist(
      # model for the observed IFA
      IFA_mu ~ normal (IFA_true, IFA_se),
      
      # model for the "true" RDP
      vector[N]:IFA_true ~ normal (mu,sigma),
      mu <- a + bBM*BM_true[i],
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),
      
      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      sigma ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_IFA_me_prec <- precis (m_IFA_me)
m_IFA_me_prec <- data.frame (round (m_IFA_me_prec,2))
write.csv (data.frame (m_IFA_me_prec), file="m_IFA_me_prec.csv")

m_IFA_me_post <- extract.samples (m_IFA_me) # get posterior samples
write.csv (data.frame (m_IFA_me_post), file="m_IFA_me_post.csv")


### measurement error model with phylogeny

m_IFA_me_phy <- ulam(
   alist(
      # model for the observed RDP
      IFA_mu ~ normal (IFA_true, IFA_se),
      
      # model for the "true" RDP with phylogenetic covariance matrix
      vector[N]:IFA_true ~ multi_normal (mu,SIGMA),
      mu <- a + bBM*BM_true[i],
      
      # Gaussian process prior for the covariance in RDP
      matrix[N,N]:SIGMA <- cov_GPL1 (Dmat, etasq, rhosq, 0.01),
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),
      
      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      etasq ~ exponential (1),
      rhosq ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_IFA_me_phy_prec <- precis (m_IFA_me_phy)
m_IFA_me_phy_prec <- data.frame (round (m_IFA_me_phy_prec,2))
write.csv (data.frame (m_IFA_me_phy_prec), file="m_IFA_me_phy_prec.csv")

m_IFA_me_phy_post <- extract.samples (m_IFA_me_phy) # get posterior samples
write.csv (data.frame (m_IFA_me_phy_post), file="m_IFA_me_phy_post.csv")


### repeat but without naked mole rats

tree <- read.nexus ("mole-rat_tree.nex") # read tree without out-groups
tree <- ls.consensus (tree) # create a consensus tree 
tree <- drop.tip (tree, tip=c("Heterocephalus_glaber","Bathyergus_janetta")) # drop H. glaber and B. janetta
spp_i <- c("Heliophobius_argenteocinereus", "Fukomys_mechowi", "Fukomys_damarensis", "Cryptomys_hottentotus", "Georychus_capensis", "Bathyergus_suillus")
tree_trimmed <- keep.tip (tree, spp_i) # combine tips of the phylogeny with societies in the sample
Dmat <- cophenetic (tree_trimmed) # compute distance matrix
Dmat <- Dmat[spp_i,spp_i]/max(Dmat) # normalize the distances to scale from 0 to 1 and include it to the data list

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# log-transform the data

log_RDP <- log (data$RDP)
log_TJI <- log (data$TJI)
log_BM <- log (data$BM)

# compute species means and standard errors for each variable
# define standard error function

se_fun <- function (x) {sd(x)/sqrt(length(x))}

RDP_mu <- tapply (log_RDP, data$species, mean)
RDP_se <- tapply (log_RDP, data$species, se_fun)

TJI_mu <- tapply (log_TJI, data$species, mean)
TJI_se <- tapply (log_TJI, data$species, se_fun)

BM_mu <- tapply (log_BM, data$species, mean)
BM_se <- tapply (log_BM, data$species, se_fun)


# create a data list

d <- list(
   RDP_mu=RDP_mu,
   TJI_mu=TJI_mu,
   TJI_se=TJI_se,
   BM_mu=BM_mu,
   RDP_se=RDP_se,
   TJI_se=TJI_se,
   BM_se=BM_se,
   N=nrow(RDP_mu),
   Dmat=Dmat
)


### RDP
### measurement error model without phylogeny

m_RDP_me <- ulam(
   alist(
      # model for the observed RDP
      RDP_mu ~ normal (RDP_true, RDP_se),
      
      # model for the "true" RDP
      vector[N]:RDP_true ~ normal (mu,sigma),
      mu <- a + bBM*BM_true[i],
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),
      
      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      sigma ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_RDP_me_prec <- precis (m_RDP_me)
m_RDP_me_prec <- data.frame (round (m_RDP_me_prec,2))
write.csv (data.frame (m_RDP_me_prec), file="m_RDP_me_prec.csv")

m_RDP_me_post <- extract.samples (m_RDP_me) # get posterior samples
write.csv (data.frame (m_RDP_me_post), file="m_RDP_me_post.csv")


### measurement error model with phylogeny

m_RDP_me_phy <- ulam(
   alist(
      # model for the observed RDP
      RDP_mu ~ normal (RDP_true, RDP_se),
      
      # model for the "true" RDP with phylogenetic covariance matrix
      vector[N]:RDP_true ~ multi_normal (mu,SIGMA),
      mu <- a + bBM*BM_true[i],
      
      # Gaussian process prior for the covariance in RDP
      matrix[N,N]:SIGMA <- cov_GPL1 (Dmat, etasq, rhosq, 0.01),
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),
      
      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      etasq ~ exponential (1),
      rhosq ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_RDP_me_phy_prec <- precis (m_RDP_me_phy)
m_RDP_me_phy_prec <- data.frame (round (m_RDP_me_phy_prec,2))
write.csv (data.frame (m_RDP_me_phy_prec), file="m_RDP_me_phy_prec.csv")

m_RDP_me_phy_post <- extract.samples (m_RDP_me_phy) # get posterior samples
write.csv (data.frame (m_RDP_me_phy_post), file="m_RDP_me_phy_post.csv")


### TJI
### measurement error model without phylogeny

m_TJI_me <- ulam(
   alist(
      # model for the observed RDP
      TJI_mu ~ normal (TJI_true, TJI_se),
      
      # model for the "true" RDP
      vector[N]:TJI_true ~ normal (mu,sigma),
      mu <- a + bBM*BM_true[i],
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),
      
      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      sigma ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_TJI_me_prec <- precis (m_TJI_me)
m_TJI_me_prec <- data.frame (round (m_TJI_me_prec,2))
write.csv (data.frame (m_TJI_me_prec), file="m_TJI_me_prec.csv")

m_TJI_me_post <- extract.samples (m_TJI_me) # get posterior samples
write.csv (data.frame (m_TJI_me_post), file="m_TJI_me_post.csv")


### measurement error model with phylogeny

m_TJI_me_phy <- ulam(
   alist(
      # model for the observed RDP
      TJI_mu ~ normal (TJI_true, TJI_se),
      
      # model for the "true" RDP with phylogenetic covariance matrix
      vector[N]:TJI_true ~ multi_normal (mu,SIGMA),
      mu <- a + bBM*BM_true[i],
      
      # Gaussian process prior for the covariance in RDP
      matrix[N,N]:SIGMA <- cov_GPL1 (Dmat, etasq, rhosq, 0.01),
      
      # model for the observed BM
      BM_mu ~ normal (BM_true,BM_se),
      
      # model for the "true" BM
      vector[N]:BM_true ~ normal (0,0.5),

      # priors
      a ~ normal (0,0.5),
      bBM ~ normal (0,0.5),
      etasq ~ exponential (1),
      rhosq ~ exponential (1)
      
   ), data=d, chains=4, cores=4, iter=5000,
   control=list(adapt_delta=0.99,
                max_treedepth=15))

m_TJI_me_phy_prec <- precis (m_TJI_me_phy)
m_TJI_me_phy_prec <- data.frame (round (m_TJI_me_phy_prec,2))
write.csv (data.frame (m_TJI_me_phy_prec), file="m_TJI_me_phy_prec.csv")

m_TJI_me_phy_post <- extract.samples (m_TJI_me_phy) # get posterior samples
write.csv (data.frame (m_TJI_me_phy_post), file="m_TJI_me_phy_post.csv")


##################################
### posterior predictive plots ###

### RDT ml
# load posteriors

RDP_ml_post <- read.csv ("m_RDP_ml_post.csv", header=TRUE)
TJI_ml_post <- read.csv ("m_TJI_ml_post.csv", header=TRUE)
IFA_ml_post <- read.csv ("m_IFA_ml_post.csv", header=TRUE)
RDP_me_post <- read.csv ("m_RDP_me_post.csv", header=TRUE)
TJI_me_post <- read.csv ("m_TJI_me_post.csv", header=TRUE)
IFA_me_post <- read.csv ("m_IFA_me_post.csv", header=TRUE)
RDP_me_phy_post <- read.csv ("m_RDP_me_phy_post.csv", header=TRUE)
TJI_me_phy_post <- read.csv ("m_TJI_me_phy_post.csv", header=TRUE)
IFA_me_phy_post <- read.csv ("m_IFA_me_phy_post.csv", header=TRUE)

{

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id
   
par (mfrow=c(2,2), mai=c(0.7,0.7,0.2,0.05))

### RDT
# B. suillus

plot (NULL, xlab="log (body mass)", ylab="log (RDT) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.75,-0.35), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.75,-0.35, length.out=5), labels=c("-0.75","-0.65","-0.55","-0.45","-0.35"), cex.axis=1.5)
mtext ("A", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (5.5,7.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.1 + RDP_ml_post$bBM.1*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points (log(data$BM[species==1]), log(data$RDP[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq (3.5,5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.2 + RDP_ml_post$bBM.2*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points (log(data$BM[species==2]), log(data$RDP[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq (4,5.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.3 + RDP_ml_post$bBM.3*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points (log(data$BM[species==3]), log(data$RDP[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq (5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.4 + RDP_ml_post$bBM.4*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points (log(data$BM[species==4]), log(data$RDP[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq (4,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.5 + RDP_ml_post$bBM.5*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points (log(data$BM[species==5]), log(data$RDP[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq (4.5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) RDP_ml_post$a.6 + RDP_ml_post$bBM.6*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points (log(data$BM[species==6]), log(data$RDP[species==6]), pch=5, col="turquoise2")


### TJI
# B. suillus

plot (NULL, xlab="log (body mass)", ylab="log (TJI) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.85,-0.45), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.85,-0.45, length.out=5), labels=c("-0.85","-0.75","-0.65","-0.55","-0.45"), cex.axis=1.5)
mtext ("B", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (5.5,7.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.1 + TJI_ml_post$bBM.1*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points (log(data$BM[species==1]), log(data$TJI[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq (3.5,5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.2 + TJI_ml_post$bBM.2*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points (log(data$BM[species==2]), log(data$TJI[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq (4,5.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.3 + TJI_ml_post$bBM.3*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points (log(data$BM[species==3]), log(data$TJI[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq (5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.4 + TJI_ml_post$bBM.4*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points (log(data$BM[species==4]), log(data$TJI[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq (4,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.5 + TJI_ml_post$bBM.5*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points (log(data$BM[species==5]), log(data$TJI[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq (4.5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) TJI_ml_post$a.6 + TJI_ml_post$bBM.6*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points (log(data$BM[species==6]), log(data$TJI[species==6]), pch=5, col="turquoise2")


### IFA

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# B. suillus

plot (NULL, xlab="log (body mass)", ylab="log (IFA) +/- 89% CI", xlim=c(2.5,7.5), ylim=c(-1.6,-1), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(2.5,7.5, length.out=6), labels=c("2.5","3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-1.6,-1, length.out=4), labels=c("-1.60","-1.40","-1.20","-1.00"), cex.axis=1.5)
mtext ("C", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (5.5,7.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.1 + IFA_ml_post$bBM.1*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points (log(data$BM[species==1]), log(data$IFA[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq (3.5,5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.2 + IFA_ml_post$bBM.2*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points (log(data$BM[species==2]), log(data$IFA[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq (4,5.5, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.3 + IFA_ml_post$bBM.3*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points (log(data$BM[species==3]), log(data$IFA[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq (5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.4 + IFA_ml_post$bBM.4*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points (log(data$BM[species==4]), log(data$IFA[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq (4,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.5 + IFA_ml_post$bBM.5*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points (log(data$BM[species==5]), log(data$IFA[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq (4.5,6, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.6 + IFA_ml_post$bBM.6*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points (log(data$BM[species==6]), log(data$IFA[species==6]), pch=5, col="turquoise2")

# Heterocephalus glaber

x_seq <- seq (2.5,4, length.out=50)
RDT_ml <- sapply (x_seq, function(x) IFA_ml_post$a.7 + IFA_ml_post$bBM.7*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("magenta2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="magenta2", lwd=2)
points (log(data$BM[species==7]), log(data$IFA[species==7]), pch=0, col="magenta2")

}

par (xpd=NA)
sp_names <- c("B. suillus","C. hottentotus","F. damarensis","F. mechowii","G. capensis","H. argenteocinereus","H. glaber")
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2","magenta2")
pts <- c(4,6,2,2,1,5,0)
legend (9.5,-1.15, lwd=2, col=cols, pch=pts, legend=sp_names, text.font=3, cex=1.5, box.col=NA)

dev.off()


### me

{

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

log_RDP <- log (data$RDP)
log_TJI <- log (data$TJI)
log_BM <- log (data$BM)

RDP_mu <- tapply (log_RDP, data$species, mean)
TJI_mu <- tapply (log_TJI, data$species, mean)
BM_mu <- tapply (log_BM, data$species, mean)

   
par (mfrow=c(2,2), mai=c(0.7,0.7,0.2,0.05))

### RDT

plot (NULL, xlab="mean (log (body mass))", ylab="mean (log (RDT))", xlim=c(3.5,7.5), ylim=c(-0.65,-0.25), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(4,7, length.out=7), labels=c("4.0","4.5","5.0","5.5","6.0","6.5","7.0"), cex.axis=1.5)
axis(2, at=seq(-0.65,-0.25, length.out=5), labels=c("-0.65","-0.55","-0.45","-0.35","-0.25"), cex.axis=1.5)
mtext ("A", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (4,7, length.out=50)
RDP_me <- sapply (x_seq, function(x) RDP_me_post$a + RDP_me_post$bBM*x) 
RDP_me_phy <- sapply (x_seq, function(x) RDP_me_phy_post$a + RDP_me_phy_post$bBM*x) 
lines(x=x_seq, y=apply(RDP_me,2,mean), col="gray75", lty=2, lwd=3)
lines(x=x_seq, y=apply(RDP_me_phy,2,mean), col="gray75", lwd=3)
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2")
pts <- c(4,6,2,2,1,5)
points (BM_mu,RDP_mu, col=cols, pch=pts, cex=3)

# TJI

plot (NULL, xlab="mean (log (body mass))", ylab="mean (log (TJI))", xlim=c(3.5,7.5), ylim=c(-0.75,-0.35), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(4,7, length.out=7), labels=c("4.0","4.5","5.0","5.5","6.0","6.5","7.0"), cex.axis=1.5)
axis(2, at=seq(-0.75,-0.35, length.out=5), labels=c("-0.75","-0.65","-0.55","-0.45","-0.35"), cex.axis=1.5)
mtext ("B", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (4,7, length.out=50)
TJI_me <- sapply (x_seq, function(x) TJI_me_post$a + TJI_me_post$bBM*x) 
TJI_me_phy <- sapply (x_seq, function(x) TJI_me_phy_post$a + TJI_me_phy_post$bBM*x) 
lines(x=x_seq, y=apply(TJI_me,2,mean), col="gray75", lty=2, lwd=3)
lines(x=x_seq, y=apply(TJI_me_phy,2,mean), col="gray75", lwd=3)
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2")
pts <- c(4,6,2,2,1,5)
points (BM_mu,TJI_mu, col=cols, pch=pts, cex=3)

### IFA

data <- read.csv ("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform (data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

log_IFA <- log (data$IFA)
log_BM <- log (data$BM)
IFA_mu <- tapply (log_IFA, data$species, mean)
BM_mu <- tapply (log_BM, data$species, mean)

plot (NULL, xlab="mean (log (body mass))", ylab="mean (log (IFA))", xlim=c(2.5,7.5), ylim=c(-1.5,-0.5), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3,7, length.out=5), labels=c("3.0","4.0","5.0","6.0","7.0"), cex.axis=1.5)
axis(2, at=seq(-1.5,-0.5, length.out=5), labels=c("-1.50","-1.25","-1.00","-0.75","-0.50"), cex.axis=1.5)
mtext ("C", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq (3,7, length.out=50)
IFA_me <- sapply (x_seq, function(x) IFA_me_post$a + IFA_me_post$bBM*x) 
IFA_me_phy <- sapply (x_seq, function(x) IFA_me_phy_post$a + IFA_me_phy_post$bBM*x) 
lines(x=x_seq, y=apply(IFA_me,2,mean), col="gray75", lty=2, lwd=3)
lines(x=x_seq, y=apply(IFA_me_phy,2,mean), col="gray75", lwd=3)
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2","magenta2")
pts <- c(4,6,2,2,1,5,0)
points (BM_mu,IFA_mu, col=cols, pch=pts, cex=3)

}

par (xpd=NA)
sp_names <- c("B. suillus","C. hottentotus","F. damarensis","F. mechowii","G. capensis","H. argenteocinereus","H. glaber")
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2","magenta2")
pts <- c(4,6,2,2,1,5,0)
legend (9.5,-0.65, col="gray75", lwd=2, lty=c(1,2), legend=c("phylogeny","no phylogeny"), text.font=1, cex=1.5, box.col=NA)
legend (10,-0.85, col=cols, pch=pts, legend=sp_names, text.font=3, cex=1.5, box.col=NA)

dev.off()


############################
### covariance functions ###

{
   
par (mfrow=c(2,2), mai=c(0.7,0.7,0.2,0.05))
   
# RDP

plot (NULL, xlim=c(0,1), ylim=c(0,10), xlab="", ylab="", main="", yaxt="n", xaxt="n", cex.axis=1.5, axes=FALSE)
axis (1, at=seq(from=0, to=1, length.out=5), labels=c("0.0","0.25","0.5","0.75","1.0"), cex.axis=1.5)
axis (2, at=seq(from=0, to=10, length.out=5), labels=c("2","4","6","8","10"), cex.axis=1.5)
mtext ("A", 3, line=-0.5, adj=0.05, cex=1.5)
mtext ("covariance +/- 89% CI", 2, line=2.75, cex=1.5)
mtext ("phylogenetic distance", 1, line=3.15, cex=1.5)

# prior

nsim <- 1000
eta <- abs(rexp(nsim,1))
rho <- abs(rexp(nsim,1))

d_seq <- seq (from=0, to=1, length.out=50)
K <- sapply (d_seq , function(x) eta*exp(-rho*x))
K_mu <- apply (K,2,mean)
K_int <- apply (K,2,PI)
shade (K_int, d_seq, col=col.alpha("gray75",0.55))
lines (d_seq, K_mu, lwd=4, col="gray75")

# posterior 

x_seq <- seq (from=0, to=1, length.out=50)
mcov <- sapply (x_seq, function(x) RDP_me_phy_post$etasq*exp(-RDP_me_phy_post$rhosq*x))
mcov_mu <- apply (mcov,2,mean)
mcov_int <- apply (mcov,2,PI)
shade(mcov_int, x_seq, col=col.alpha("gray25",0.55))
lines (x_seq, mcov_mu, lwd=4, col="gray25")

names <- c("prior","posterior")
cols <- c("gray75","gray25")
legend (0.5,9, col=cols, legend=names, lwd=3, lty=1, title="RDT", title.adj=0.6, cex=1.5, box.col=NA, horiz=FALSE, xpd=NA)


# TJI

plot (NULL, xlim=c(0,1), ylim=c(0,10), xlab="", ylab="", main="", yaxt="n", xaxt="n", cex.axis=1.5, axes=FALSE)
axis (1, at=seq(from=0, to=1, length.out=5), labels=c("0.0","0.25","0.5","0.75","1.0"), cex.axis=1.5)
axis (2, at=seq(from=0, to=10, length.out=5), labels=c("2","4","6","8","10"), cex.axis=1.5)
mtext ("B", 3, line=-0.5, adj=0.05, cex=1.5)
mtext ("covariance +/- 89% CI", 2, line=2.75, cex=1.5)
mtext ("phylogenetic distance", 1, line=3.15, cex=1.5)

# prior

nsim <- 1000
eta <- abs(rexp(nsim,1))
rho <- abs(rexp(nsim,1))

d_seq <- seq (from=0, to=1, length.out=50)
K <- sapply (d_seq , function(x) eta*exp(-rho*x))
K_mu <- apply (K,2,mean)
K_int <- apply (K,2,PI)
shade (K_int, d_seq, col=col.alpha("gray75",0.55))
lines (d_seq, K_mu, lwd=4, col="gray75")

# posterior 

x_seq <- seq (from=0, to=1, length.out=50)
mcov <- sapply (x_seq, function(x) TJI_me_phy_post$etasq*exp(-TJI_me_phy_post$rhosq*x))
mcov_mu <- apply (mcov,2,mean)
mcov_int <- apply (mcov,2,PI)
shade(mcov_int, x_seq, col=col.alpha("gray25",0.55))
lines (x_seq, mcov_mu, lwd=4, col="gray25")

names <- c("prior","posterior")
cols <- c("gray75","gray25")
legend (0.5,9, col=cols, legend=names, lwd=3, lty=1, title="TJI", title.adj=0.55, cex=1.5, box.col=NA, horiz=FALSE, xpd=NA)


# IFA

plot (NULL, xlim=c(0,1), ylim=c(0,10), xlab="", ylab="", main="", yaxt="n", xaxt="n", cex.axis=1.5, axes=FALSE)
axis (1, at=seq(from=0, to=1, length.out=5), labels=c("0.0","0.25","0.5","0.75","1.0"), cex.axis=1.5)
axis (2, at=seq(from=0, to=10, length.out=5), labels=c("2","4","6","8","10"), cex.axis=1.5)
mtext ("C", 3, line=-0.5, adj=0.05, cex=1.5)
mtext ("covariance +/- 89% CI", 2, line=2.75, cex=1.5)
mtext ("phylogenetic distance", 1, line=3.15, cex=1.5)

# prior

nsim <- 1000
eta <- abs(rexp(nsim,1))
rho <- abs(rexp(nsim,1))

d_seq <- seq (from=0, to=1, length.out=50)
K <- sapply (d_seq , function(x) eta*exp(-rho*x))
K_mu <- apply (K,2,mean)
K_int <- apply (K,2,PI)
shade (K_int, d_seq, col=col.alpha("gray75",0.55))
lines (d_seq, K_mu, lwd=4, col="gray75")

# posterior 

x_seq <- seq (from=0, to=1, length.out=50)
mcov <- sapply (x_seq, function(x) IFA_me_phy_post$etasq*exp(-IFA_me_phy_post$rhosq*x))
mcov_mu <- apply (mcov,2,mean)
mcov_int <- apply (mcov,2,PI)
shade(mcov_int, x_seq, col=col.alpha("gray25",0.55))
lines (x_seq, mcov_mu, lwd=4, col="gray25")

names <- c("prior","posterior")
cols <- c("gray75","gray25")
legend (0.5,9, col=cols, legend=names, lwd=3, lty=1, title="IFA", title.adj=0.55, cex=1.5, box.col=NA, horiz=FALSE, xpd=NA)

}

dev.off()


####################################
### stochastic character mapping ###

# load the data and trees

tree_data <- read.csv ("tree_data.csv", header=TRUE, row.names=1)
tree <- read.nexus ("mole-rat_tree-out.nex")
cons_tree <- ls.consensus (tree) # create consensus tree to plot mapped characters on

# replace qualitative statements with integer values

tree_data$DT <- ifelse (tree_data$DT=="Yes",1,0)
tree_data$TT <- ifelse (tree_data$TT=="Yes",1,0) 
tree_data$DFTF <- ifelse (tree_data$DFTF=="Yes",1,0) 
tree_data$SS <- ifelse (tree_data$SS=="SC",1,0) 
oldvals <- c("NO","CT","SD")
newvals <- as.integer (c(0,1,2))
tree_data$DM <- newvals[match(tree_data$DM, oldvals)]


# DT

DT <- setNames (tree_data$DT, rownames(tree_data))
sim_DT <- make.simmap (tree, DT, model="ARD", Q="empirical")
sim_DT_pd <- summary (sim_DT, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(DT)))
sim_DT_ace <- sim_DT_pd$ace
write.csv (sim_DT_ace, file="sim_DT_ace.csv", row.names=FALSE) # save the reconstructed node probabilities

# TT

TT <- setNames (tree_data$TT, rownames(tree_data))
sim_TT <- make.simmap (tree, TT, model="ARD", Q="empirical")
sim_TT_pd <- summary (sim_TT, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(TT)))
sim_TT_ace <- sim_TT_pd$ace
write.csv (sim_TT_ace, file="sim_TT_ace.csv", row.names=FALSE)

# DFTF

DFTF <- setNames (tree_data$DFTF, rownames(tree_data))
sim_DFTF <- make.simmap (tree, DFTF, model="ARD", Q="empirical")
sim_DFTF_pd <- summary (sim_DFTF, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(DFTF)))
sim_DFTF_ace <- sim_DFTF_pd$ace
write.csv (sim_DFTF_ace, file="sim_DFTF_ace.csv", row.names=FALSE)

# SS

SS <- setNames (tree_data$SS, rownames(tree_data))
sim_SS <- make.simmap (tree, SS, model="ARD", Q="empirical")
sim_SS_pd <- summary (sim_SS, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(SS)))
sim_SS_ace <- sim_SS_pd$ace
write.csv (sim_SS_ace, file="sim_SS_ace.csv", row.names=FALSE)

# DM

DM <- setNames (tree_data$DM, rownames(tree_data))
sim_DM <- make.simmap (tree, DM, model="ARD", Q="empirical")
sim_DM_pd <- summary (sim_DM, consensus.tree=cons_tree)
cols <- setNames(c("gray35","gray55","gray85"), levels(as.factor(DM)))
sim_DM_ace <- sim_DM_pd$ace
write.csv (sim_DM_ace, file="sim_DM_ace.csv", row.names=FALSE)


### plot reconstructed probabilities for each trait 
   
plot (sim_DT_pd,
         offset=0.5,
         fsize=1,
         lwd=1,
         cex=c(0.75,0.75),
         colors=cols
)
   
add.simmap.legend (x=0.5,
                   y=1.5,
                   colors=setNames(c("gray35","gray85"), 
                                   c("DT present","DT absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.15
)


plot (sim_TT_pd,
         offset=0.5,
         fsize=1,
         lwd=1,
         cex=c(0.75,0.75),
         colors=cols
)
   
add.simmap.legend (x=0.5,
                   y=1.5,
                   colors=setNames(c("gray35","gray85"), 
                                   c("TT present","TT absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.15
)
   
   
plot (sim_DFTF_pd,
         offset=0.5,
         fsize=1,
         lwd=1,
         cex=c(0.75,0.75),
         colors=cols
)
   
add.simmap.legend (x=0.5,
                   y=1.5,
                   colors=setNames(c("gray35","gray85"), 
                                   c("DFTF present","DFTF absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.15
)
   
   
plot (sim_SS_pd,
         offset=0.5,
         fsize=1,
         lwd=1,
         cex=c(0.75,0.75),
         colors=cols
)
   
add.simmap.legend (x=0.5,
                   y=1.5,
                   colors=setNames(c("gray35","gray85"), 
                                   c("solitary","social")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.15
)
   
   
plot (sim_DM_pd,
         offset=0.5,
         fsize=1,
         lwd=1,
         cex=c(0.75,0.75),
         colors=cols
)
   
add.simmap.legend (x=0.5,
                   y=1.75,
                   colors=setNames(c("gray35","gray55","gray85"), 
                                   c("no digging","chisel-tooth","scratch-digging")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.15
)


################################################################################
################################################################################