###################################################################################################################################
### R code for: Fossorial adaptations in African mole-rats (Bathyergidae), the unique appendicular phenotype of naked mole-rats ###
###################################################################################################################################


library(phytools)
library(rethinking)


# settings recommended by Stan

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

memory.limit(size=2e6) # increase memory


###############################################
### phylogenetic varying effects regression ###

# we will fit two models for each index - with and without body mass 
# then we will compare the two models with WAIC to assess the relative importance of body mass for prediction

data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[,1:4] # index columns for modelling IFA
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

tree <- read.nexus("mole-rat_tree.nex") # read the tree
tree <- ls.consensus(tree) # create a consensus tree 
tree <- drop.tip(tree, tip=c("Bathyergus_janetta")) # drop B. janetta
spp_i <- c("Heliophobius_argenteocinereus","Fukomys_mechowi","Fukomys_damarensis","Cryptomys_hottentotus","Georychus_capensis","Bathyergus_suillus","Heterocephalus_glaber")
tree_trimmed <- keep.tip(tree, spp_i) # combine tips of the phylogeny with societies in the sample
R_OU <- corMartins(1,phy=tree_trimmed,form=~spp_i,fixed=FALSE) # adjust the tree according to the OU model but let the alpha be estimated rather than fixed
V <- vcv(R_OU) # convert to a variance-covariance matrix
V <- V[spp_i,spp_i] # order the species
R <- V/max(V) # convert to a correlation matrix
N <- length(spp_i)

# log-transform the data

log_IFA <- log(data$IFA)
log_BM <- log(data$BM)

# create a data list

d <- list(
   log_IFA=log_IFA,
   log_BM=log_BM,
   species=species,
   R=R,
   N=N
)


# fit a model with body mass

m_IFA_bm <- ulam(
   alist(
      # Gaussian model for the outcome
      log_IFA ~ normal(mu,sigma), 
      mu <- a[species] + bBM[species]*log_BM + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S), # "rep_vector(mu_a,N)" doesn't work :/
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # adaptive prior for slopes 
      bBM[species] ~ multi_normal(c(mu_b,mu_b,mu_b,mu_b,mu_b,mu_b,mu_b),S),
      matrix[N,N]:S <- sigma_phy_b*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept and slope means
      mu_a ~ normal(0,0.5),
      mu_b ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_phy_b,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_IFA_bm_prec <- precis(m_IFA_bm,2) # get model summary
m_IFA_bm_prec <- data.frame(round(m_IFA_bm_prec,2))
write.csv(data.frame(m_IFA_bm_prec), file="m_IFA_bm_prec.csv")

m_IFA_bm_post <- extract.samples(m_IFA_bm) # get posterior samples
write.csv(data.frame(m_IFA_bm_post), file="m_IFA_bm_post.csv")


# fit a model without body mass

m_IFA_int <- ulam(
   alist(
      # Gaussian model for the outcome
      log_IFA ~ normal(mu,sigma), 
      mu <- a[species] + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S),
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept mean
      mu_a ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_IFA_int_prec <- precis(m_IFA_int,2) # get model summary
m_IFA_int_prec <- data.frame(round(m_IFA_int_prec,2))
write.csv(data.frame(m_IFA_int_prec), file="m_IFA_int_prec.csv")

m_IFA_int_post <- extract.samples(m_IFA_int) # get posterior samples
write.csv(data.frame(m_IFA_int_post), file="m_IFA_int_post.csv")


### model comparison ###

IFA_model_comp <- compare(m_IFA_int, m_IFA_bm, func=WAIC)
IFA_model_comp <- data.frame(round(IFA_model_comp,2))
write.csv(IFA_model_comp, file="IFA_model_comp.csv")


### RDT and TJI ###

data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# create correlation matrix again but without H. glaber

tree <- drop.tip(tree, tip=c("Heterocephalus_glaber")) # drop H. glaber
spp_i <- c("Heliophobius_argenteocinereus","Fukomys_mechowi","Fukomys_damarensis","Cryptomys_hottentotus","Georychus_capensis","Bathyergus_suillus")
tree_trimmed <- keep.tip(tree, spp_i) # combine tips of the phylogeny with societies in the sample
R_OU <- corMartins(1,phy=tree_trimmed,form=~spp_i,fixed=FALSE) # adjust the tree according to the OU model but let the alpha be estimated rather than fixed
V <- vcv(R_OU) # convert to a variance-covariance matrix
V <- V[spp_i,spp_i] # order the species
R <- V/max(V) # convert to a correlation matrix
N <- length(spp_i)

# log-transform the data

log_RDT <- log(data$RDT)
log_TJI <- log(data$TJI)
log_BM <- log(data$BM)

# create a data list

d <- list(
   log_RDT=log_RDT,
   log_TJI=log_TJI,
   log_BM=log_BM,
   species=species,
   R=R,
   N=N
)


# RDT
# fit a model with body mass

m_RDT_bm <- ulam(
   alist(
      # Gaussian model for the outcome
      log_RDT ~ normal(mu,sigma), 
      mu <- a[species] + bBM[species]*log_BM + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S),
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # adaptive prior for slopes 
      bBM[species] ~ multi_normal(c(mu_b,mu_b,mu_b,mu_b,mu_b,mu_b),S),
      matrix[N,N]:S <- sigma_phy_b*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept and slope means
      mu_a ~ normal(0,0.5),
      mu_b ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_phy_b,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_RDT_bm_prec <- precis(m_RDT_bm,2) # get model summary
m_RDT_bm_prec <- data.frame(round(m_RDT_bm_prec,2))
write.csv(data.frame(m_RDT_bm_prec), file="m_RDT_bm_prec.csv")

m_RDT_bm_post <- extract.samples(m_RDT_bm) # get posterior samples
write.csv(data.frame(m_RDT_bm_post), file="m_RDT_bm_post.csv")


# fit a model without body mass

m_RDT_int <- ulam(
   alist(
      # Gaussian model for the outcome
      log_RDT ~ normal(mu,sigma), 
      mu <- a[species] + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S),
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept and slope means
      mu_a ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_RDT_int_prec <- precis(m_RDT_int,2) # get model summary
m_RDT_int_prec <- data.frame(round(m_RDT_int_prec,2))
write.csv(data.frame(m_RDT_int_prec), file="m_RDT_int_prec.csv")

m_RDT_int_post <- extract.samples(m_RDT_int) # get posterior samples
write.csv(data.frame(m_RDT_int_post), file="m_RDT_int_post.csv")


### model comparison ###

RDT_model_comp <- compare(m_RDT_int, m_RDT_bm, func=WAIC)
RDT_model_comp <- data.frame(round(RDT_model_comp,2))
write.csv(RDT_model_comp, file="RDT_model_comp.csv")


# TJI
# fit a model with body mass

m_TJI_bm <- ulam(
   alist(
      # Gaussian model for the outcome
      log_TJI ~ normal(mu,sigma), 
      mu <- a[species] + bBM[species]*log_BM + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S),
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # adaptive prior for slopes 
      bBM[species] ~ multi_normal(c(mu_b,mu_b,mu_b,mu_b,mu_b,mu_b),S),
      matrix[N,N]:S <- sigma_phy_b*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept and slope means
      mu_a ~ normal(0,0.5),
      mu_b ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_phy_b,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_TJI_bm_prec <- precis(m_TJI_bm,2) # get model summary
m_TJI_bm_prec <- data.frame(round(m_TJI_bm_prec,2))
write.csv(data.frame(m_TJI_bm_prec), file="m_TJI_bm_prec.csv")

m_TJI_bm_post <- extract.samples(m_TJI_bm) # get posterior samples
write.csv(data.frame(m_TJI_bm_post), file="m_TJI_bm_post.csv")


# fit a model without body mass

m_TJI_int <- ulam(
   alist(
      # Gaussian model for the outcome
      log_TJI ~ normal(mu,sigma), 
      mu <- a[species] + e[species],
      
      # adaptive prior for intercepts 
      a[species] ~ multi_normal(c(mu_a,mu_a,mu_a,mu_a,mu_a,mu_a),S),
      matrix[N,N]:S <- sigma_phy_a*R,
      
      # model for species-specific effects
      e[species] ~ normal(0,sigma_spp),
      
      # hyper-priors for intercept and slope means
      mu_a ~ normal(0,0.5),
      
      # priors for standard deviations
      sigma ~ exponential(2),
      c(sigma_phy_a,sigma_spp) ~ exponential(2)
      
   ), data=d, chains=4, cores=4, iter=5000, log_lik=TRUE,
   control=list(adapt_delta=0.99,
                max_treedepth=15))


m_TJI_int_prec <- precis(m_TJI_int,2) # get model summary
m_TJI_int_prec <- data.frame(round(m_TJI_int_prec,2))
write.csv(data.frame(m_TJI_int_prec), file="m_TJI_int_prec.csv")

m_TJI_int_post <- extract.samples(m_TJI_int) # get posterior samples
write.csv(data.frame(m_TJI_int_post), file="m_TJI_int_post.csv")


### model comparison ###

TJI_model_comp <- compare(m_TJI_int, m_TJI_bm, func=WAIC)
TJI_model_comp <- data.frame(round(TJI_model_comp,2))
write.csv(TJI_model_comp, file="TJI_model_comp.csv")


####################################################
### posterior predictive plots - varying effects ###

RDT_post <- read.csv("m_RDT_post.csv", header=TRUE)
TJI_post <- read.csv("m_TJI_post.csv", header=TRUE)
IFA_post <- read.csv("m_IFA_post.csv", header=TRUE)


{

data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id
   
par(mfrow=c(2,2), mai=c(0.7,0.7,0.2,0.05))

### RDT
# B. suillus

plot(NULL, xlab="log (body mass)", ylab="log (RDT) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.7,-0.3), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.7,-0.3, length.out=5), labels=c("-0.70","-0.60","-0.50","-0.40","-0.30"), cex.axis=1.5)
mtext("A", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq(5.5,7.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.1 + RDT_post$bBM.1*x + RDT_post$e.1) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points(log(data$BM[species==1]), log(data$RDT[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq(3.5,5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.2 + RDT_post$bBM.2*x + RDT_post$e.2) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points(log(data$BM[species==2]), log(data$RDT[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq(4,5.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.3 + RDT_post$bBM.3*x + RDT_post$e.3) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points(log(data$BM[species==3]), log(data$RDT[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq(5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.4 + RDT_post$bBM.4*x + RDT_post$e.4) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points(log(data$BM[species==4]), log(data$RDT[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq(4,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.5 + RDT_post$bBM.5*x + RDT_post$e.5) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points(log(data$BM[species==5]), log(data$RDT[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq(4.5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$a.6 + RDT_post$bBM.6*x + RDT_post$e.6) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points(log(data$BM[species==6]), log(data$RDT[species==6]), pch=5, col="turquoise2")


### TJI
# B. suillus

plot(NULL, xlab="log (body mass)", ylab="log (TJI) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.80,-0.40), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.80,-0.40, length.out=5), labels=c("-0.80","-0.70","-0.60","-0.50","-0.40"), cex.axis=1.5)
mtext("B", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq(5.5,7.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.1 + TJI_post$bBM.1*x + TJI_post$e.1) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points(log(data$BM[species==1]), log(data$TJI[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq(3.5,5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.2 + TJI_post$bBM.2*x + TJI_post$e.2) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points(log(data$BM[species==2]), log(data$TJI[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq(4,5.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.3 + TJI_post$bBM.3*x + TJI_post$e.3) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points(log(data$BM[species==3]), log(data$TJI[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq(5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.4 + TJI_post$bBM.4*x + TJI_post$e.4) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points(log(data$BM[species==4]), log(data$TJI[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq(4,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.5 + TJI_post$bBM.5*x + TJI_post$e.5) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points(log(data$BM[species==5]), log(data$TJI[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq(4.5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) TJI_post$a.6 + TJI_post$bBM.6*x + TJI_post$e.6) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points(log(data$BM[species==6]), log(data$TJI[species==6]), pch=5, col="turquoise2")


### IFA

data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

# B. suillus

plot(NULL, xlab="log (body mass)", ylab="log (IFA) +/- 89% CI", xlim=c(2.5,7.5), ylim=c(-1.6,-1), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(2.5,7.5, length.out=6), labels=c("2.5","3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-1.6,-1, length.out=4), labels=c("-1.60","-1.40","-1.20","-1.00"), cex.axis=1.5)
mtext("C", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq(5.5,7.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.1 + IFA_post$bBM.1*x + IFA_post$e.1) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("springgreen4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="springgreen4", lwd=2)
points(log(data$BM[species==1]), log(data$IFA[species==1]), pch=4, col="springgreen4")

# C. hottentotus

x_seq <- seq(3.5,5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.2 + IFA_post$bBM.2*x+ IFA_post$e.2) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("goldenrod1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="goldenrod1", lwd=2)
points(log(data$BM[species==2]), log(data$IFA[species==2]), pch=6, col="goldenrod1")

# F. damarensis

x_seq <- seq(4,5.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.3 + IFA_post$bBM.3*x+ IFA_post$e.3) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("darkorange1",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="darkorange1", lwd=2)
points(log(data$BM[species==3]), log(data$IFA[species==3]), pch=2, col="darkorange1")

# F. mechowii

x_seq <- seq(5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.4 + IFA_post$bBM.4*x+ IFA_post$e.4) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("red2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="red2", lwd=2)
points(log(data$BM[species==4]), log(data$IFA[species==4]), pch=2, col="red2")

# G. capensis

x_seq <- seq(4,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.5 + IFA_post$bBM.5*x+ IFA_post$e.5) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("royalblue4",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="royalblue4", lwd=2)
points(log(data$BM[species==5]), log(data$IFA[species==5]), pch=1, col="royalblue4")

# Heliophobius argenteocinereus 

x_seq <- seq(4.5,6, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.6 + IFA_post$bBM.6*x+ IFA_post$e.6) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("turquoise2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="turquoise2", lwd=2)
points(log(data$BM[species==6]), log(data$IFA[species==6]), pch=5, col="turquoise2")

# Heterocephalus glaber

x_seq <- seq(2.5,4, length.out=50)
RDT_ml <- sapply(x_seq, function(x) IFA_post$a.7 + IFA_post$bBM.7*x+ IFA_post$e.7) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("magenta2",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="magenta2", lwd=2)
points(log(data$BM[species==7]), log(data$IFA[species==7]), pch=0, col="magenta2")

}

par(xpd=NA)
sp_names <- c("B. suillus","C. hottentotus","F. damarensis","F. mechowii","G. capensis","H. argenteocinereus","H. glaber")
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2","magenta2")
pts <- c(4,6,2,2,1,5,0)
legend(9.5,-1.15, lwd=2, col=cols, pch=pts, legend=sp_names, text.font=3, cex=1.5, box.col=NA)

dev.off()


##########################################################
### posterior predictive plots - population parameters ###

{
   
data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

par(mfrow=c(2,2), mai=c(0.7,0.7,0.2,0.05))

### RDT

plot(NULL, xlab="log (body mass)", ylab="log (RDT) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.90,-0.10), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.90,-0.10, length.out=5), labels=c("-0.90","-0.70","-0.50","-0.30","-0.10"), cex.axis=1.5)
mtext("A", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq(3.5,7.5, length.out=50)
RDT_ml <- sapply(x_seq, function(x) RDT_post$mu_a + RDT_post$mu_b*x) # linear model
shade(apply(RDT_ml,2,PI),x_seq,col=col.alpha("gray55",0.3))
lines(x=x_seq, y=apply(RDT_ml,2,mean), col="gray55", lwd=3)

points(log(data$BM[species==1]), log(data$RDT[species==1]), pch=4, col="springgreen4")
points(log(data$BM[species==2]), log(data$RDT[species==2]), pch=6, col="goldenrod1")
points(log(data$BM[species==3]), log(data$RDT[species==3]), pch=2, col="darkorange1")
points(log(data$BM[species==4]), log(data$RDT[species==4]), pch=2, col="red2")
points(log(data$BM[species==5]), log(data$RDT[species==5]), pch=1, col="royalblue4")
points(log(data$BM[species==6]), log(data$RDT[species==6]), pch=5, col="turquoise2")


### TJI

plot(NULL, xlab="log (body mass)", ylab="log (TJI) +/- 89% CI", xlim=c(3.5,7.5), ylim=c(-0.90,-0.10), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(3.5,7.5, length.out=5), labels=c("3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-0.90,-0.10, length.out=5), labels=c("-0.90","-0.70","-0.50","-0.30","-0.10"), cex.axis=1.5)
mtext("B", 3, line=-0.5, adj=0.05, cex=1.5)

TJI_ml <- sapply(x_seq, function(x) TJI_post$mu_a + TJI_post$mu_b*x) # linear model
shade(apply(TJI_ml,2,PI),x_seq,col=col.alpha("gray55",0.3))
lines(x=x_seq, y=apply(TJI_ml,2,mean), col="gray55", lwd=3)

points(log(data$BM[species==1]), log(data$TJI[species==1]), pch=4, col="springgreen4")
points(log(data$BM[species==2]), log(data$TJI[species==2]), pch=6, col="goldenrod1")
points(log(data$BM[species==3]), log(data$TJI[species==3]), pch=2, col="darkorange1")
points(log(data$BM[species==4]), log(data$TJI[species==4]), pch=2, col="red2")
points(log(data$BM[species==5]), log(data$TJI[species==5]), pch=1, col="royalblue4")
points(log(data$BM[species==6]), log(data$TJI[species==6]), pch=5, col="turquoise2")


### IFA

data <- read.csv("quantitative_data.csv", header=TRUE)
data <- data[,1:4]
data <- data[complete.cases(data),] # drop rows with missing values
data <- transform(data, species_id=as.numeric(factor(data$species))) # assign a unique number to each species
species <- data$species_id

plot(NULL, xlab="log (body mass)", ylab="log (IFA) +/- 89% CI", xlim=c(2.5,7.5), ylim=c(-1.7,-0.5), yaxt="n", xaxt="n", cex.lab=1.5, axes=FALSE)
axis(1, at=seq(2.5,7.5, length.out=6), labels=c("2.5","3.5","4.5","5.5","6.5","7.5"), cex.axis=1.5)
axis(2, at=seq(-1.7,-0.5, length.out=5), labels=c("-1.70","-1.40","-1.10","-0.80","-0.50"), cex.axis=1.5)
mtext("C", 3, line=-0.5, adj=0.05, cex=1.5)

x_seq <- seq(2.5,7.5, length.out=50)
IFA_ml <- sapply(x_seq, function(x) IFA_post$mu_a + IFA_post$mu_b*x) # linear model
shade(apply(IFA_ml,2,PI),x_seq,col=col.alpha("gray55",0.3))
lines(x=x_seq, y=apply(IFA_ml,2,mean), col="gray55", lwd=3)

points(log(data$BM[species==1]), log(data$IFA[species==1]), pch=4, col="springgreen4")
points(log(data$BM[species==2]), log(data$IFA[species==2]), pch=6, col="goldenrod1")
points(log(data$BM[species==3]), log(data$IFA[species==3]), pch=2, col="darkorange1")
points(log(data$BM[species==4]), log(data$IFA[species==4]), pch=2, col="red2")
points(log(data$BM[species==5]), log(data$IFA[species==5]), pch=1, col="royalblue4")
points(log(data$BM[species==6]), log(data$IFA[species==6]), pch=5, col="turquoise2")
points(log(data$BM[species==7]), log(data$IFA[species==7]), pch=0, col="magenta2")
   
}

par(xpd=NA)
sp_names <- c("B. suillus","C. hottentotus","F. damarensis","F. mechowii","G. capensis","H. argenteocinereus","H. glaber")
cols <- c("springgreen4","goldenrod1","darkorange1","red2","royalblue4","turquoise2","magenta2")
pts <- c(4,6,2,2,1,5,0)
legend(9.5,-0.8, lwd=2, col=cols, pch=pts, legend=sp_names, text.font=3, cex=1.5, box.col=NA)

dev.off()


####################################
### stochastic character mapping ###

# load the data and trees

tree_data <- read.csv("tree_data.csv", header=TRUE, row.names=1)
tree <- read.nexus("mole-rat_tree-out.nex")
cons_tree <- ls.consensus(tree) # create consensus tree to plot mapped characters on

# replace qualitative statements with integer values

tree_data$DT <- ifelse(tree_data$DT=="Yes",1,0)
tree_data$TT <- ifelse(tree_data$TT=="Yes",1,0) 
tree_data$DFTF <- ifelse(tree_data$DFTF=="Yes",1,0) 
tree_data$SS <- ifelse(tree_data$SS=="SC",1,0) 
oldvals <- c("NO","CT","SD")
newvals <- as.integer(c(0,1,2))
tree_data$DM <- newvals[match(tree_data$DM, oldvals)]


# DT

DT <- setNames(tree_data$DT, rownames(tree_data))
sim_DT <- make.simmap(tree, DT, model="ARD", Q="empirical")
sim_DT_pd <- summary(sim_DT, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(DT)))
sim_DT_ace <- sim_DT_pd$ace
write.csv(sim_DT_ace, file="sim_DT_ace.csv", row.names=FALSE) # save the reconstructed node probabilities

# TT

TT <- setNames(tree_data$TT, rownames(tree_data))
sim_TT <- make.simmap(tree, TT, model="ARD", Q="empirical")
sim_TT_pd <- summary(sim_TT, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(TT)))
sim_TT_ace <- sim_TT_pd$ace
write.csv(sim_TT_ace, file="sim_TT_ace.csv", row.names=FALSE)

# DFTF

DFTF <- setNames(tree_data$DFTF, rownames(tree_data))
sim_DFTF <- make.simmap(tree, DFTF, model="ARD", Q="empirical")
sim_DFTF_pd <- summary(sim_DFTF, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(DFTF)))
sim_DFTF_ace <- sim_DFTF_pd$ace
write.csv(sim_DFTF_ace, file="sim_DFTF_ace.csv", row.names=FALSE)

# SS

SS <- setNames(tree_data$SS, rownames(tree_data))
sim_SS <- make.simmap(tree, SS, model="ARD", Q="empirical")
sim_SS_pd <- summary(sim_SS, consensus.tree=cons_tree)
cols <- setNames(c("gray85","gray35"), levels(as.factor(SS)))
sim_SS_ace <- sim_SS_pd$ace
write.csv(sim_SS_ace, file="sim_SS_ace.csv", row.names=FALSE)

# DM

DM <- setNames(tree_data$DM, rownames(tree_data))
sim_DM <- make.simmap(tree, DM, model="ARD", Q="empirical")
sim_DM_pd <- summary(sim_DM, consensus.tree=cons_tree)
cols <- setNames(c("gray35","gray55","gray85"), levels(as.factor(DM)))
sim_DM_ace <- sim_DM_pd$ace
write.csv(sim_DM_ace, file="sim_DM_ace.csv", row.names=FALSE)


### plot reconstructed probabilities for each trait 
   
plot(sim_DT_pd,
        offset=0.5,
        fsize=1,
        lwd=1,
        cex=c(0.75,0.75),
        colors=cols
)
   
add.simmap.legend(x=0.5,
                  y=1.5,
                  colors=setNames(c("gray35","gray85"), 
                                c("DT present","DT absent")),
                  prompt=FALSE,
                  vertical=TRUE,
                  shape="circle",
                  cex=0.15
)


plot(sim_TT_pd,
        offset=0.5,
        fsize=1,
        lwd=1,
        cex=c(0.75,0.75),
        colors=cols
)
   
add.simmap.legend(x=0.5,
                  y=1.5,
                  colors=setNames(c("gray35","gray85"), 
                                  c("TT present","TT absent")),
                  prompt=FALSE,
                  vertical=TRUE,
                  shape="circle",
                  cex=0.15
)
   
   
plot(sim_DFTF_pd,
        offset=0.5,
        fsize=1,
        lwd=1,
        cex=c(0.75,0.75),
        colors=cols
)
   
add.simmap.legend(x=0.5,
                  y=1.5,
                  colors=setNames(c("gray35","gray85"), 
                                  c("DFTF present","DFTF absent")),
                  prompt=FALSE,
                  vertical=TRUE,
                  shape="circle",
                  cex=0.15
)
   
   
plot(sim_SS_pd,
        offset=0.5,
        fsize=1,
        lwd=1,
        cex=c(0.75,0.75),
        colors=cols
)
   
add.simmap.legend(x=0.5,
                  y=1.5,
                  colors=setNames(c("gray35","gray85"), 
                                  c("solitary","social")),
                  prompt=FALSE,
                  vertical=TRUE,
                  shape="circle",
                  cex=0.15
)
   
   
plot(sim_DM_pd,
        offset=0.5,
        fsize=1,
        lwd=1,
        cex=c(0.75,0.75),
        colors=cols
)
   
add.simmap.legend(x=0.5,
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