#===============================================================
# Project: Random change point models in Stan
# File:    Fitting the random change point model to a single
#          simulated dataset; this file is included as an 
#          example in the supplementary material for the paper:
#            "Bayesian piecewise linear mixed models with a 
#             random change point: an application to BMI  
#              rebound in childhood"
#===============================================================


#============================================================================
# Preliminaries - load required packages
# Notes: If the following packages aren't already installed, then 
#        they can be  installed from CRAN by typing, for example: 
#        install.packages("package-name"). Note that this requires an 
#        internet connection. Installing RStan for the first time requires
#        a few additional steps, see http://mc-stan.org/interfaces/rstan
#        for details. 
#============================================================================

library(MASS)    # package for mvrnorm
library(doBy)    # package for summaryBy
library(rstan)   # package for running Stan from within R
library(lattice) # package for xyplot
library(latticeExtra)


#============================================================================
# Function to create a single simulated dataset. 
# Notes: The function returns an R list, containing within it two objects: 
#        (i) a vector with the named elements of the vector providing the 
#            true parameter values used in creating the simulated data, and
#        (ii) a list containing the simulated data, to be utilised by RStan.
#        Note that the dataset must be provided as an R list, since this is 
#        the format required by RStan. If you are analysing your own dataset, 
#        then it will need to be converted from its usual format (for 
#        example a data frame) to a named R list.
#======================================================================

createdata <- function(
  seed                     = sample.int(.Machine$integer.max, 1), 
  Npat                     = 100, 
  Npat.forpred             = 16,
  agemin                   = 1, 
  agemax                   = 15,
  intercept.truemean       = 15, 
  intercept.truesd         = 1.2,
  preslope.truemean        = -0.4, 
  preslope.truesd          = 0.2,
  postslope.truemean       = 0.6, 
  postslope.truesd         = 0.3,
  kp.truemean              = 6.5, 
  kp.truesd                = 1.4, 
  betakp.lowerbound        = 3, 
  betakp.upperbound        = 9, 
  error.truesd             = 1,
  corr.intercept.preslope  = 0.3, 
  corr.intercept.postslope = 0.2, 
  corr.intercept.kp        = -0.3,
  corr.preslope.postslope  = -0.2, 
  corr.preslope.kp         = 0.3, 
  corr.postslope.kp        = -0.3) {
  
  set.seed(seed)
  
  # Number of measurements per patient
  Ntimepoints      <- as.integer(runif(Npat, 5, 25))       
  id.long          <- rep(1:Npat, times = Ntimepoints)
  
  # Measurement times
  age              <- runif(n = length(id.long), min = agemin, max = agemax)    
  
  # Measurement errors
  error            <- rnorm(n = length(id.long), mean = 0, sd = error.truesd)
  
  # Long data (measurement times and error terms)
  dat1             <- data.frame(id=id.long, age=age, error=error)        
  
  # True random effects (correlated): intercept, preslope, postslope, knot point
  id.short         <- rep(1:Npat)
  mu               <- c(intercept.truemean, 
                        preslope.truemean, 
                        postslope.truemean, 
                        kp.truemean)
  var1             <- intercept.truesd ^ 2
  var2             <- preslope.truesd ^ 2
  var3             <- postslope.truesd ^ 2
  var4             <- kp.truesd ^ 2
  cov12            <- corr.intercept.preslope * intercept.truesd * preslope.truesd
  cov13            <- corr.intercept.postslope * intercept.truesd * postslope.truesd
  cov14            <- corr.intercept.kp * intercept.truesd * kp.truesd 
  cov23            <- corr.preslope.postslope * preslope.truesd * postslope.truesd
  cov24            <- corr.preslope.kp * preslope.truesd * kp.truesd
  cov34            <- corr.postslope.kp * postslope.truesd * kp.truesd 
  covmatrix        <- matrix(c(var1, cov12, cov13, cov14,  
                        cov12, var2, cov23, cov24,
                        cov13, cov23, var3, cov34,
                        cov14, cov24, cov34, var4), 
                      nrow = 4, ncol = 4)
  cormatrix        <- cov2cor(covmatrix)   
  simparms         <- mvrnorm(n = Npat, mu = mu, Sigma = covmatrix)
  
  # Short data (true random effect parameters: intercept, slopes and knot point)
  dat2             <- data.frame(id = id.short, simparms)       
  
  # Merge long and short data
  simdat           <- merge(dat1, dat2, by = "id")
  simdat           <- simdat[order(simdat$id, simdat$age),]
  
  # Observed outcome
  simdat$y         <- simdat$X1 + 
                      simdat$X2 * pmin(simdat$age - simdat$X4, 0) + 
                      simdat$X3 * pmax(simdat$age - simdat$X4, 0) + 
                      simdat$error        
  
  # Prediction dataset
  temp             <- summaryBy(age ~ id, data = simdat, FUN = c(min, max))
  temp             <- temp[temp$id <= Npat.forpred,]  
  Npred.perpatient <- (agemax - agemin) / 0.1 + 1
  preddata.id      <- rep(temp$id, each = Npred.perpatient)
  preddata.agemin  <- rep(temp$age.min, each = Npred.perpatient)
  preddata.agemax  <- rep(temp$age.max, each = Npred.perpatient)
  preddata.age     <- rep(seq(agemin, agemax, 0.1), times = nrow(temp))
  preddata         <- data.frame(cbind(id = preddata.id, 
                                       agemin = preddata.agemin, 
                                       agemax = preddata.agemax, 
                                       age = preddata.age))
  preddata         <- preddata[(preddata$age > preddata$agemin) & 
                                 (preddata$age < preddata$agemax) , c("id", "age")]
  
  # Return data for Stan
  return(list(
    trueparms = c(
      "beta[1]"       = intercept.truemean, 
      "beta[2]"       = preslope.truemean, 
      "beta[3]"       = postslope.truemean,
      "betakp"        = kp.truemean,
      "u_sd[1]"       = intercept.truesd,
      "u_sd[2]"       = preslope.truesd,
      "u_sd[3]"       = postslope.truesd,
      "u_sd[4]"       = kp.truesd,
      "ukp_sd"        = kp.truesd,
      "u_Corr[1,2]"   = corr.intercept.preslope,
      "u_Corr[1,3]"   = corr.intercept.postslope,
      "u_Corr[1,4]"   = corr.intercept.kp,
      "u_Corr[2,3]"   = corr.preslope.postslope,
      "u_Corr[2,4]"   = corr.preslope.kp,
      "u_Corr[3,4]"   = corr.postslope.kp,  
      "u_Sigma[1,2]"  = corr.intercept.preslope * intercept.truesd * preslope.truesd,
      "u_Sigma[1,3]"  = corr.intercept.postslope * intercept.truesd * postslope.truesd,
      "u_Sigma[1,4]"  = corr.intercept.kp * intercept.truesd * kp.truesd,
      "u_Sigma[2,3]"  = corr.preslope.postslope * preslope.truesd * postslope.truesd,
      "u_Sigma[2,4]"  = corr.preslope.kp * preslope.truesd * kp.truesd,
      "u_Sigma[3,4]"  = corr.postslope.kp * postslope.truesd * kp.truesd,
      "y_sd"          = error.truesd),
    standata = list(
      "N"             = nrow(simdat),
      "Npat"          = Npat,
      "fixedkp"       = kp.truemean,
      "zeros3"        = rep(0,3),
      "zeros4"        = rep(0,4),
      "betakp_lower"  = betakp.lowerbound,
      "betakp_upper"  = betakp.upperbound,
      "id"            = simdat$id,
      "age"           = simdat$age, 
      "y"             = simdat$y,
      "Npred"         = nrow(preddata),
      "Npat_pred" = Npat.forpred,
      "id_pred"       = preddata$id,
      "age_pred"      = preddata$age,
      "agemin"        = agemin,
      "agemax"        = agemax)))
} 


#==============================================================
# Create simulated data, and fit the random change point model
#==============================================================

# Create a single simulated dataset with specified seed
standata <- createdata(seed = 123456)

# Specify the location of the Stan model code file
# Notes: The stancode.filepath object should be changed, so that it contains
#        the pathway to the Stan model code, wherever that file is installed 
#        on your PC.  
stancode.filepath <- "stancode_randomknotcorr.stan"

# Specify which parameters we wish to monitor the MCMC samples for
# Notes: "Monitoring" a parameter means that the MCMC samples related to 
#        that parameter are retained, so they can be used later for post-
#        processing and inference. The default is for RStan to monitor all 
#        model parameters.
#        However to save memory, we wish to only monitor those parameters
#        which will be used for inference or checking model diagnostics.
#        The parameter names are listed in a vector, and the parameters
#        names must align with those that are used in the Stan model code
#        file (ie, the file "stancode_randomknotcorr.stan").
#        Here, we choose to monitor the log posterior (lp__), all key 
#        model parameters, as well as prediction and random effect estimates
#        for a subset of the patients (see the Stan model code and the data  
#        steps above for more detail).
stanmonitor <- c("beta", 
                 "betakp", 
                 "y_sd", 
                 "u_sd", 
                 "u_Sigma", 
                 "u_Corr",
                 "y_mu_pred", 
                 "y_pred", 
                 "alpha_tosave",
                 "lp__")

# Fit the random change point model to the simulated data
# Notes: Here we are fitting 3 chains in parallel, using 3 PC cores. The 
#        number of cores may need to be adjusted according to the PC being
#        used to fit the model. The control option is used to adjust some
#        of the parameters used for the Hamiltonian Monte Carlo sampler; 
#        see the RStan vignette (http://mc-stan.org/interfaces/rstan) for 
#        a discussion of these.
stanfit <- stan(file    = stancode.filepath, 
                data    = standata$standata, 
                pars    = stanmonitor, 
                chains  = 3, 
                cores   = 3, 
                iter    = 3000, 
                warmup  = 1000,
                control = list(adapt_delta = 0.8, max_treedepth = 15))

#=======================================
# Some diagnostics for the fitted model
#=======================================

# Trace plots of the MCMC samples for various model parameters
traceplot(stanfit, 
          pars = c("beta", "betakp"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("u_sd"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("u_Corr"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("lp__"), 
          inc_warmup = TRUE)

# Diagnostics for Hamiltonian Monte Carlo sampler 
# Notes:
#   See page 15 of the RStan vignette (http://mc-stan.org/interfaces/rstan)
#   for some discussion of these diagnostics
summary(do.call(rbind, 
                args = get_sampler_params(stanfit, inc_warmup = FALSE)),
        digits = 2)


#=============================
# Summary of model parameters
#=============================

print(stanfit, pars = c("beta", 
                        "betakp", 
                        "y_sd", 
                        "u_sd", 
                        "u_Sigma", 
                        "u_Corr"))

#=====================================================================
# Plot observed height measurements with fitted trajectories overlaid
#=====================================================================

# Extract MCMC samples for the model
mcmcsamples.pred <- extract(stanfit)

# Calculate the mean trajectory based on the posterior mean for each of the model parameters
pred.data <- data.frame(cbind(id = standata$standata$id_pred, 
                              age = standata$standata$age_pred))
pred.alphas <- data.frame(cbind(id = c(1:length(unique(standata$standata$id_pred))), 
                                apply(mcmcsamples.pred$alpha_tosave, 
                                      MARGIN = c(2,3), 
                                      FUN = mean)))
pred.data <- merge(pred.data, pred.alphas, by = "id")
pred.data[, "y.predtraj"] <- pred.data$V2 + 
                             pred.data$V3 * pmin(pred.data$age - pred.data$V5, 0) + 
                             pred.data$V4 * pmax(pred.data$age - pred.data$V5, 0)

# Calculate percentiles of the posterior predictive distribution
pred.p025 <- data.frame(cbind(id = standata$standata$id_pred, 
                              age = standata$standata$age_pred, 
                              y.predlimit = apply(mcmcsamples.pred$y_pred, 
                                                  MARGIN = 2, 
                                                  FUN = quantile, probs = .025) ))
pred.p975 <- data.frame(cbind(id = standata$standata$id_pred, 
                              age = standata$standata$age_pred, 
                              y.predlimit = apply(mcmcsamples.pred$y_pred, 
                                                  MARGIN = 2, 
                                                  FUN = quantile, probs = .975) ))
pred.limits <- rbind(pred.p975, pred.p025[order(pred.p025$age, decreasing = TRUE), ])

# Observed data
obs.data <- data.frame(cbind(id = standata$standata$id, 
                             age = standata$standata$age, 
                             y = standata$standata$y))

# Create plot
a <- xyplot(y.predlimit ~ age | id, 
            data = pred.limits[pred.limits$id %in% c(1:10), ],
            panel = panel.polygon,
            xlab = "Age (years)", ylab = "BMI (kg/m2)", strip=FALSE, 
            scales = list(alternating = FALSE, relation = "sliced"),
            layout = c(5,2), col = 'lightgray', border = 'lightgray')
b <- xyplot(y ~ age | id, 
            data = obs.data[obs.data$id %in% c(1:10), ],
            panel = panel.xyplot,
            type = 'p', pch = 19, col = 'black')
c <- xyplot(y.predtraj ~ age | id,
            data = pred.data[pred.data$id %in% c(1:10), ],
            panel = panel.xyplot,
            type = 'l', lty = 2, col = 'black')
print(a + as.layer(b) + as.layer(c))







