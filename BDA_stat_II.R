# https://oscarbonilla.com/2009/05/visualizing-bayes-theorem/
# https://michaeldewittjr.com/resources/stan_linear_regression.html
# https://docs.google.com/document/d/1NERhwwyWfYr2ThNnVU74-wKvYAUO8RI0NLbzbdbebu4/edit?ts=5ffc976dhttps://docs.google.com/document/d/1NERhwwyWfYr2ThNnVU74-wKvYAUO8RI0NLbzbdbebu4/edit?ts=5ffc976d
# https://avehtari.github.io/BDA_course_Aalto/demos.html
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# https://mc-stan.org/posterior/
#
# Stan User Guide:       https://mc-stan.org/docs/2_26/stan-users-guide/summing-out-the-responsibility-parameter.html
# Stan Reference Manual: https://mc-stan.org/docs/2_26/reference-manual/increment-log-prob-statement.html
#
# Ch 14 - Bayesian Workflow
# https://mc-stan.org/users/documentation/case-studies/planetary_motion/planetary_motion.html#sec7
#
# Class Questions/Answers
# https://docs.google.com/document/d/1NERhwwyWfYr2ThNnVU74-wKvYAUO8RI0NLbzbdbebu4/edit?ts=600f3c8e
#
# http://www.stat.columbia.edu/~gelman/research/published/entropy-19-00555-v2.pdf

# https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# no pooling, hierarchical model, and complete pooling.  No pooling is hierarchical model 
#     with tau=infinity, complete pooling is with tau=0

colVars <- function(a) {n <- dim(a)[[1]]; c <- dim(a)[[2]]; return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))}


library(devtools)
library(ggplot2)
library(ggrepel)
library("loo")
library("brms")
library("MASS")
library(foreign)
library(reshape2)
library(cowplot)
library("sjPlot")
library("gridExtra")
library ("arm")
library(dplyr)
library(tidyr)
library(binr)
install.packages("latex2exp")
library(latex2exp)
library("cmdstanr")
library(posterior)
library(bayesplot)
library(viridis)
library("rstanarm")
library("rstan")
install.packages("formatR")
library(formatR)
library (stringr)
# example(stan_model, package = "rstan", run.dontrun = TRUE)


#all.equal(test_a, post_bivar$post_d_given_a)
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

invlogit <- plogis
theme_set(theme_bw())
theme_set(theme_cowplot())
options(digits=4)
set.seed (1234)
##########
bivar_norm_dense <- function (x, y, mu_x, mu_y, sigma_x, sigma_y, rho){
  c1 <- (2*(1-rho^2))
  c2 <- 2*pi*sigma_x*sigma_y*sqrt(1-rho^2)
  
  z_x <- (x-mu_x)/sigma_x
  z_y <- (y-mu_y)/sigma_y
  
  dense <- exp(-(z_x^2 - 2*rho *z_x*z_y + z_y^2)/c1)/c2
  return (dense)
}

##### 2-b
#
# a.  Write a Stan program to fit the model, y_n = a + b*x_n + error_n, for n=1,...,N, 
#     with errors that are independent and normally distributed with mean 0 and standard 
#     deviation sigma.  Assume the parameters a, b are restricted to be positive.
# b.  In R, simulate fake data for this model with N=100, x uniformly distributed between 
#     0 and 10, and a, b, sigma taking on the values 2, 3, 0.2.
# c.  Fit the Stan model using your simulated data and check that the true parameter values 
#     are approximately recovered.  Check also that you get approximately the same answer 
#     as from fitting a classical linear regression.
# d.  Make a single graph showing a scatter plot of the simulated data and the fitted model.

N <- 100
a <- 2
b <- 3
sigma <- 0.2

x <- runif (N, 0, 10)
y <- a + b*x + rnorm (N, 0, sigma)
plot (x, y)
data1 <- list (N=N, x=x, y=y)
#data_df <- data.frame(x, y)
#fit_lin <- stan_glm (y ~ x, data = data_df, refresh = 0)
#print (fit_lin)

#linear_fit <- stan(file="3a.stan", data=data1, iter=1000, chains=4)
fit1 <- stan (file = 'linear.stan', data = data1)

linear_model <- stan_model("linear.stan")
fit1 <- sampling (linear_model, data=data1)

print (fit1)
plot (fit1)
traceplot (fit1)
pairs (fit1)

##### 3-a
# Same problem as from previous homework, except the model is 
# y_n = 1/(a + b*x_n) + error_n, for n=1,...,N.  (Stan note:  in writing the program, 
# if your code is vectorized you have to write 1 ./ (a + b*x), not 1 / (a + b*x). 

N <- 100
a <- 2
b <- 3
sigma <- 0.2

x <- runif (N, 0, 10)
y <- 1/(a + b*x) + rnorm (N, 0, sigma)
plot (x, y)
data2 <- list (N=N,x=x, y=y)

inv_linear_model <- cmdstan_model("inv_linear.stan", pedantic=T)
fit2 <- inv_linear_model$sample(data=data2, refresh = 0)
fit2$summary()
fit2$summary("sigma", "mean", "sd")

draws <- fit2$draws()
str(draws)
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)

color_scheme_set("red")
mcmc_intervals()#, pars = c("cyl", "drat", "am", "sigma"))
mcmc_hist(fit2$draws("sigma"))

# convert to regular stan output
stanfit2 <- rstan::read_stan_csv(fit2$output_files())
# regular stan output
print (stanfit2, digits=4)
plot (stanfit2)
traceplot (stanfit2)
pairs (stanfit2)

### Ch 5: Hierarchical Models
# EX 5.1
# n_j - number of rats in j's experiment
# y_j - number of rates with tumors in j's experiment
# assume y_j ∼ Bin (n_j,theta_j),
# theta_j - probability of tumor in control group of rats in experiment j
# theta ~ Beta (alpha, beta)

y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
       2,1,5,2,5,3,2,7,7,3,3,2,9,10,4,4,4,4,4,4,4,10,4,4,4,5,11,12,
       5,5,6,5,6,6,6,6,16,15,15,9,4)
n <- c(20,20,20,20,20,20,20,19,19,19,19,18,18,17,20,20,20,20,19,19,18,18,25,24,
       23,20,20,20,20,20,20,10,49,19,46,27,17,49,47,20,20,13,48,50,20,20,20,20,
       20,20,20,48,19,19,19,22,46,49,20,20,23,19,22,20,20,20,52,46,47,24,14)
prop <- y/n

E_theta <- mean (prop[1:70])   # 0.13
var_theta <- sd (prop[1:70])^2 # 0.0107

# Calculate parmaters of prior Beta(alfa, beta)
#
# alfa + beta = [E(theta)*(1 - E(theta))]/var(theta)  - 1 
# alfa = (alfa + beta) * E(theta)
# beta = (alfa + beta) * (1 - E(theta))

a_plus_b <- E_theta*(1-E_theta)/var_theta - 1
a <- a_plus_b* E_theta     # 1.4
b <- a_plus_b*(1-E_theta)  # 8.6
# Beta_prior (1.4, 8.6)

# Posterior Beta_post given 1 data pair y (success), n (tries):  Beta (a + y, b + n - y)
# Beta (a + y[71], b + n[71] - y[71])
a_post <- a + y[71]
b_post <- b + n[71] - y[71]
# Beta_posterior [5.4, 18.6]
# now find E(theta_71) and var(theta_71) from Beta_posterior[5.4, 18.6]
# E_theta <- alfa/(alfa + beta)
# var_theta <- [E(theta)*(1 - E(theta))]/(alfa + beta + 1)

a_post_plus_b_post <- a_post + b_post

E_theta_post   <- a_post/a_post_plus_b_post                                # 0.223
var_theta_post <- E_theta_post*(1 - E_theta_post)/(a_post_plus_b_post + 1) # 0.0069
sd_theta_post  <- sqrt (var_theta_post)                                    # 0.833

### Full Bayes ####

x <- seq(0.0001, 0.9999, length.out = 1000)
bdens <- function(n, y, x)
  dbeta(x, y+1, n-y+1)

df_sep <- mapply(bdens, n, y, MoreArgs = list(x = x)) %>%
  as.data.frame() %>% cbind(x) %>% gather(ind, p, -x)

labs1 <- paste('posterior of', c('theta_j', 'theta_71'))
plot_sep <- ggplot(data = df_sep) +
  geom_line(aes(x = x, y = p, color = (ind=='V71'), group = ind)) +
  labs(x = expression(theta), y = '', title = 'Separate model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))
# The last one is for emphasize colored red
plot_sep

a <- sum(y)+1
b <- sum(n)-sum(y)+1
a_plus_b <- a + b

E_theta_post_pool <- a/a_plus_b                             # 0.136
var_theta_post_pool <- E_theta*(1 - E_theta)/(a_plus_b + 1) # 0.0000675
sd_theta_post_pool <- sqrt (var_theta_post_pool)
df_pool <- data.frame(x = x, p = dbeta(x, sum(y)+1, sum(n)-sum(y)+1))
plot_pool <- ggplot(data = df_pool) +
  geom_line(aes(x = x, y = p, color = '1')) +
  labs(x = expression(theta), y = '', title = 'Pooled model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = 'red', labels = 'Posterior of common theta') +
  theme(legend.background = element_blank(), legend.position = c(0.7,0.9))
plot_pool           

#
# 2.4 placenta previa 
# data generating process is binomial.. Conjugates are beta functions
#
beta_mean_var <- function (prior, num_prior_data, num_post_data, num_post_success){
  a_prior <- (num_prior_data + 2) * prior
  b_prior <- (num_prior_data + 2) - a_prior
  
  a_post <- a_prior + num_post_success
  b_post <- b_prior + num_post_data - num_post_success
  
  a_plus_b <- a_post + b_post
  E_theta <- a_post/a_plus_b                           
  var_theta <- E_theta*(1 - E_theta)/(a_plus_b + 1)
  sd_theta <- sqrt (var_theta)
  theta <- list (a_prior=a_prior, b_prior=b_prior, E=E_theta, sd=sd_theta )
  return (theta)
}

n_data <- 980 # births
y_data <- 437 # num girls, 543 boys

flat_prior <- c(0.5, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485)
num_prior  <- c(2, 2, 5, 10, 20, 100, 200) - 2
theta <- beta_mean_var(flat_prior, num_prior, n_data, y_data) 

cat ("prior  a    b    a+b       Post mean (theta)        2 sd CI")
for (i in 1:7){
  cat (sprintf("Prior=%.3f a_prior=%6.2f b_prior=%6.2f a+b=%3d E(theta)=%.3f CI=[%.3f , %.3f]\n", 
               prior[i], theta$a_prior[i], theta$b_prior[i], theta$a_prior[i] + theta$b_prior[i], 
               theta$E[i], theta$E[i]-2*theta$sd[i], theta$E[i]+2*theta$sd[i]))
  
}


#
# 2.7 Informative prior distribution for cancer rates
# data generation process is Poisson distribution..  Conjugates are Gamma functions
# Gamma (alfa, beta)    E(Gamma) = alfa/beta, sd (Gamma) = sqrt(alfa)/beta
# if prior=Gamma (alfa, beta) => posterior=Gamma (alfa + successes, beta + tries)
#
# J counties
# n_j = population of county j
# y_j = number of cancer deaths between 1980-89 in county j
# theta_j = true cancer death rate per person per year
# y_j ~ Poisson (10* theta_j * n_j)

# theta_prior_j = Gamma (alfa=20, beta=430000)
alfa<- 20
beta<- 430000
E_theta_prior <- alfa/beta # 4.65E-5
sd_theta_prior <- sqrt (alfa)/beta

# POSTERIOR: P(theta_j | y_j) = theta_post_j <- Gamma (alfa + y_j, beta + 10*n_j)
E_theta_post_j <- (20 + y_j)/(430000 + 10*n_j)
sd_theta_post_j <- sqrt (20 + y_j)/(430000 + 10*n_j)

# if n_j=1000 and n_j=0
E_theta_post <- (20 + 0)/(430000 + 10*1000) # 4.55E-5, close to prior of 4.65E-5
# if n_j=1000 and n_j=1
E_theta_post <- (20 + 1)/(430000 + 10*1000) # 4.77E-5, still close to prior of 4.65E-5
# if n_j=1000 and n_j=2
E_theta_post <- (20 + 2)/(430000 + 10*1000) # 5.00E-5, still close to prior of 4.65E-5
# if n_j large, E_theta_post_j = y_j/(10*n_j)

# How many deaths would we expect in 10 years in a county of N_cnty

theta <- rgamma (500, 20, 430000)


N_cnty <- 2500000
num_deaths <- rep (NA, 500)
for (i in 1:500){
  num_deaths[i] <- rpois (1, 10*N_cnty*theta)
}
num_deaths <- rpois (500, 10*N_cnty*theta)

cat (sprintf("-50pcnt= %d  median= %d +50pcnt= %d", 
             sort(num_deaths)[125], sort (num_deaths)[249], sort(num_deaths)[375]))
cat ("expected deaths per year per 100000" = mean(num_deaths)*100000/N_cnty/10)

num_bars <- max (num_deaths) - min (num_deaths) + 1

title <- paste ("County Size =", N_cnty, " people")
subtitle <- paste("median=", median (num_deaths), "and 50% range")
ggplot (data.frame (num_deaths), aes(x=num_deaths)) + 
  geom_histogram (aes(x=num_deaths, fill=..count..), binwidth=0.5*num_bars/4) +
#  geom_histogram (aes(y = 100*(..count..)/sum(..count..)), binwidth=0.5) +
  geom_vline(xintercept = median (num_deaths), color="red") +
  geom_vline(xintercept = (sort (num_deaths))[125], color="red", linetype="dashed") +
  geom_vline(xintercept = (sort (num_deaths))[375], color="red", linetype="dashed") +
  labs (x="num deaths", y="Occurance num deaths in sample of 500", title = title, subtitle=subtitle)

# How many deaths would we expect in 10 years in a county of population N_city 
N_city <- 0.001*10E6
num_deaths_per_year <- rpois (500, 10*N_city*theta)
mean (num_deaths_per_year/N_city*100000) 
sd (num_deaths_per_year/N_city*100000)  # 105
int_50_pcnt <- 2/3 * sd (num_deaths)
cat (sprintf("Expected deaths in 10 years in county of %4.1fM = %.0f, 50 pct CI = [%.0f , %.0f]", 
             N_city/1000000, mean (num_deaths), mean (num_deaths)-2/3*sd (num_deaths),  mean (num_deaths)+2/3*sd (num_deaths) ))

# Homework 6-6
# Perform a posterior predictive check, using the same test quantity, T = number of switches, 
# but simulating the replications yrep under the new measurement protocol. Display the 
# predictive simulations, T(yrep), and discuss how they differ from Figure 6.5.
     
num_switches <- function (y_rep, num_zeros){
  zeros <- 0
  switch <- 0
  comp <- y_rep[1]
  i <- 1
  while (zeros <= num_zeros -1){
    if (y_rep[i] == 0) zeros <- zeros + 1
    i <- i + 1
    if (y_rep[i] != comp){
      switch <- switch + 1
#      print (switch)
      comp <- y_rep[i]
    }
  }
  res <- list (switch = switch, zeros = zeros, i=i-1)
  return (res)
}

# T(y_rep) - number of swithches until 13 zeros are observed

theta <- rbeta (1, 8, 14)
y_rep <- rbinom (100, 1, theta)
res <- num_switches (y_rep, 13)

cat ("last index =", res$i, "zeros = ", res$zeros, "switches = ", res$switch)
N<- 10000
n_sw <- rep (NA, N)
n_indx <- rep (NA, N)
for (i in 1:N){
  theta <- rbeta (1, 8, 14)
  y_rep <- rbinom (100, 1, theta)
  n_sw[i] <- (num_switches (y_rep, 13))$switch
  n_indx[i] <- (num_switches (y_rep, 13))$i
}

hist (n_sw, breaks = 0:max(n_sw))
hist (n_indx, breaks = 0:max(n_indx))

#### End HW 6-6


# 2-5 Normal Distribution with known variance
#
# y - scalar observation from normal distribution parameterized by mean
# theta and variance sigma^2
# sampling distribution: p(y | theta) = 1/(sqrt(2*pi)*sigma) * exp (- 1/(2*sigma^2)*(y-theta)^2)
#
# prior       p(theta)     ~ exp (- 1/(2*tau0^2)*(theta-mu0)^2)
# posterior   p(theta | y) ~ exp (- 1/(2*tau1^2)*(theta-mu1)^2)
# 1/tau1^2 = 1/tau0^2 + 1/sigma^2
# mu1 = (1/tau0^2*mu0 + 1/sigma^2*y)/(1/tau0^2 + 1/sigma^2)
#
# Posterior predictive distribution y_post
# E (y_post|y) = E(theta|y) = mu1 
# var (y_post|y) = sigma^2 + tau1^2 - sum of predictive variance sigma from model, 
# and posterior uncertainty in theta, tau1

# Normal model with n multiple observations
#
# 1/tau1^2 = 1/tau0^2 + n/sigma^2
# mu1 = (1/tau0^2*mu0 + n/sigma^2*y_bar)/(1/tau0^2 + n/sigma^2)

#
# CH-3 Multi-parameter models
#
# theta = (theta1, theta2)   - we're interested in theta1, theta2 is a nuisance parameter
#
# p (theta1, theta2|y) ~ p(y|theta1, theta2)*p(theta1, theta2)
# integrate over nuisance parameter
#             p (theta1|y) = integral [p(theta1,theta2| y * p theta2)]
# Chain Rule: p (theta1, theta2|y) = p (theta1 |theta2,y) * p(theta2|y)
#             p (theta1|y) = integral [p(theta1 | theta2,y)) * p(theta2|y) * d theta2]  


#
#
# 8 Schools
#
library("rstan") # observe startup messages
school <- c(1:8)
y =     c(28,  8, -3,  7, -1,  1, 18, 12)
sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
J <- length (y)
schools_dat <- list(J=J, y=y, sigma=sigma)

# Model: y ~ normal (theta, sigma) 
eight_schools <- cmdstan_model ("~/Documents/RStanFiles/schools.stan", pedantic = T)
schools_fit <- eight_schools$sample (data = schools_dat, refresh = 0)

# eta ~ normal (0, 1)
# Model: theta[j] = mu + tau * eta[j], 
#        y ~ normal (theta, sigma)
eight_schools1 <- cmdstan_model ("~/Documents/RStanFiles/8schools1.stan", pedantic = T)
schools_fit1 <- eight_schools1$sample(data = schools_dat, refresh = 0)
bayesplot (schools_fit1)
schools_fit1$summary()
schools_fit1$summary("theta")

draws <- schools_fit1$draws()
str(draws)
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)
cat("Chances School A intervention effect is larger than in C =",
#    100*mean (draws_df[,12] > draws_df[,14]),
    100*mean (draws[,,12] > draws[,,14]),
    "%")
###
# tau ~ exponential (1)
eight_schools3 <- cmdstan_model ("~/Documents/RStanFiles/8schools3.stan", pedantic = T)
schools_fit3 <- eight_schools3$sample(data = schools_dat, parallel_chains = 4, refresh = 0)
schools_fit3$summary()
schools_fit3$summary("theta")
# chances school A is better than school C
draws    <- schools_fit3$draws()
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
cat("Chances School A intervention effect is larger than in C =",
    100*mean (draws[,,12] > draws[,,14]), "%")

# Same hieararchical model, change y[1]
y2 =     c(200,  8, -3,  7, -1,  1, 18, 12)
schools_dat2 <- list(J=J, y=y2, sigma=sigma)
schools_fit1_A200 <- eight_schools1$sample (data = schools_dat2, refresh = 0)

schools_fit1_A200$summary("theta")
# chances school A is better than school C
draws    <- schools_fit1_A200$draws()
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
cat("Chances School A intervention effect is larger than in C =",
    100*mean (draws[,,12] > draws[,,14]), "%")
# New model
#
# eta ~ student_t (nu, 0, 1);
# nu ~ normal (7, 1);
# tau ~ cauchy (0, 25);
# Model: theta[j] = mu + tau * eta[j]
#        y ~ normal (theta, sigma)
y =     c(28,  8, -3,  7, -1,  1, 18, 12)
schools_dat <- list(J=J, y=y, sigma=sigma)
eight_schools2 <- stan_model ("~/Documents/RStanFiles/8schools2.stan")
schools_fit2 <- sampling (eight_schools2, data = schools_dat, iter=1000, chains = 4)
print (schools_fit2, pars = "theta")
mean ((as.matrix (schools_fit2))[,11] > (as.matrix (schools_fit2))[,13])

eight_schools <- stan_model ("~/Documents/RStanFiles/schools.stan")
schools_fit <- sampling (eight_schools, data = schools_dat, iter=1000, chains = 4)
# test statistic: diff between best and second best coaching
test <- function (y){
  y_sort <- sort (y, decreasing = TRUE)
  return (y_sort[1] - y_sort[2])
}
# generate output
output <- function (fit, y, sigma){
  print (fit, digits=2)
  #  plot(fit)
  #  print(fit, pars = c("mu", "tau", "theta"))
  
  schools_sim <- extract (fit, permuted=TRUE)
  # posterior inference for theta
  theta_post_mean <- array (NA, J)
  theta_post_sd <- array (NA, J)
  for (j in 1:J){
    theta_post_mean[j] <- mean (schools_sim$theta[,j])
    theta_post_sd[j] <- sd (schools_sim$theta[,j])
  }
  # posterior predictive replicated data (Make up data w/ model theta and data sigma)
  png ("8schools.png", height=500, width=600)
  n_sims <- length (schools_sim$lp__)
  y_rep <- array (NA, c(n_sims, J))
  for (s in 1:n_sims){
    y_rep[s,] <- rnorm (J, schools_sim$theta[s,], sigma)
  }
  par (mfrow=c(5,4), mar=c(2,2,2,2))
  # actual data
  hist (y, xlab="", main="y", breaks=8)
  # 19 sets of simulated 8 schools
  for (s in 1:19){
    hist (y_rep[s,], xlab="", main=paste("y_rep", s), breaks=8)
  }
  dev.off ()
  
  # diff between best and second best coaching
  t_y <- test (y)
  t_rep <- rep (NA, n_sims)
  for (s in 1:n_sims){
    t_rep[s] <- test (y_rep[s,])
  }
  par (mfrow=c(2,2), mar=c(2,2,2,2))
  # posterior inference for tau
  hist (main="Posterior Inference of Tau", xlab = "", schools_sim$tau)
  abline (v=mean(schools_sim$tau), col = "blue")
  if (length (schools_sim$nu) != 0){
    # posterior inference for Nu
    hist (main="Posterior Inference of Nu", xlab = "", schools_sim$nu)
    abline (v=mean(schools_sim$nu), col = "blue")
  }
  # histogram of replicated effect difference between best and second best coaching 
  hist (main="Inference: Best - Second Best Coaching", xlab="", t_rep)
  # actual and replicated difference between best and second best
  abline(v=c(mean(t_rep),t_y), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))
  
  cat ("\nMean Tau =", mean (schools_sim$tau))
  if (length (schools_sim$nu) != 0){
    cat ("\nMean Nu =", mean (schools_sim$nu))
  }
  cat ("\nExpected avg. coaching:  prior= ", mean(y), " posterior =", mean(theta_post_mean))
  cat ("\nPrior sd =", sd (y), " Posterior sd =", sd(theta_post_mean))
  cat ("\nBest -  Second Best coaching: Prior =", t_y, " Posterior = ", mean(t_rep))   
  cat ("\nPosterior Prob School 1 coaching is better than School 3 =", 
       mean (schools_sim$theta[,1] > schools_sim$theta[, 3])*100, "%")
  
  return (list (theta_post_mean=theta_post_mean, theta_post_sd=theta_post_sd))
}

list1 <- output (schools_fit1, y, sigma)
list2 <- output (schools_fit2, y, sigma)

range_y <- c (min(y-sigma), max(y+sigma))
range_model1 <- c (min(list1$theta_post_mean-list1$theta_post_sd), 
                   max(list1$theta_post_mean+list1$theta_post_sd)) 
range_model2 <- c (min(list2$theta_post_mean-list2$theta_post_sd), 
                   max(list2$theta_post_mean+list2$theta_post_sd)) 
range_all <- range (c(range_y, range_model1, range_model2))
range_all[1] <- range_all[1]- 10
par (mfrow=c(1,1), mar=c(4,4,2,2))
plot  (school, ylim=range_all, xlab="School", ylab="Coaching Effect", 
       main="Prior and Posterior Coaching Effect", col ="white")
points(school+.00, y,                     bty='l', pch=20, col = "black")
points(school+.05, list1$theta_post_mean, bty="l", pch=20, col = "red")
points(school+.10, list2$theta_post_mean, bty="l", pch=20, col = "green")
for (j in 1:J){
  lines(rep(school[j],2),       y[j] +                     c(-1,1)*sigma[j], lwd=.5, col="black")
  lines(rep(school[j],2) + .05, list1$theta_post_mean[j] + c(-1,1)*list1$theta_post_sd[j], lwd=.5, col="red")
  lines(rep(school[j],2) + .10, list2$theta_post_mean[j] + c(-1,1)*list2$theta_post_sd[j], lwd=.5, col="green")
  abline (h=0)
}

#########
# MCMC and Bayesplot
#########
color_scheme_set("blue")

mcmc_hist(schools_fit1$draws(c("mu", "tau")))
mcmc_hist(schools_fit1$draws("theta"))

bayesplot_grid(
  mcmc_hist(schools_fit1$draws("theta"), binwidth = 0.025),
  titles = c("Posterior distribution from MCMC"),
  xlim = c(0, 1)
)
bayesplot_grid(
  mcmc_hist(schools_fit1$draws("tau"), binwidth = 0.05),
  titles = c("Posterior distribution from MCMC"),
  xlim = c(0, 4)
)

#
# Bayesian Workflow CH 13
# Golf Putting
# 
golf2 <-  read.table(textConnection("
x n y
0.28 45198 45183
0.97 183020 182899
1.93 169503 168594
2.92 113094 108953
3.93 73855 64740
4.94 53659 41106
5.94 42991 28205
6.95 37050 21334
7.95 33275 16615
8.95 30836 13503
9.95 28637 11060
10.95 26239 9032
11.95 24636 7687
12.95 22876 6432
14.43 41267 9813
16.43 35712 7196
18.44 31573 5290
20.44 28280 4086
21.95 13238 1642
24.39 46570 4767
28.40 38422 2980
32.39 31641 1996
36.39 25604 1327
40.37 20366 834
44.38 15977 559
48.37 11770 311
52.36 8708 231
57.25 8878 204
63.23 5492 103
69.18 3087 35
75.19 1742 24"), header=TRUE) %>% 
  as_tibble () %>%
  mutate (p = y/n,
          se = sqrt (p*(1-p)/n))

golf <-  read.table(textConnection("
x n y
2 1443 1346
3 694 577
4 455 337
5 353 208
6 272 149
7 256 136
8 240 111
9 217 69
10 200 67
11 237 75
12 202 52
13 192 46
14 174 54
15 167 28
16 201 27
17 195 31
18 191 33
19 147 20
20 152 24"), header=TRUE) %>% 
  as_tibble () %>%
  mutate (p = y/n,
          se = sqrt (p*(1-p)/n))
  
ggplot (golf, aes(x=x, y=p)) + 
  geom_pointrange(data=golf, aes(x=x, ymin=p-se, ymax=p+se), color="blue", size = 0.2) +
  geom_pointrange(data=golf2, aes(x=x, ymin=p-se, ymax=p+se), color="red", size = 0.2) +
  labs (x ="Distance of putt", y = "Success rate")

# Logistic regression
# 𝑦𝑗∼binomial𝑛𝑗,logit (𝑎+𝑏𝑥𝑗),for𝑗=1,...,𝐽.
golf_data1 <- list (x=golf$x, y=golf$y, n=golf$n, J=dim(golf)[1])

golf_model <- cmdstan_model ("~/Documents/RStanFiles/golf_logistic.stan", pedantic = T)
fit_logistic<- golf_model$sample (data = golf_data1, refresh = 0)

draws <- fit_logistic$draws()
str(draws)
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)
fit_logistic$summary ()

a <- mean(draws_df$a)
b <- mean(draws_df$b)

ggplot (golf, aes(x=x, y=p)) + 
  xlim (0, 20) + ylim (0, 1) +
  geom_pointrange(aes(x=x, ymin=p-se, ymax=p+se), color="blue", size = 0.2) +
  geom_function(fun = function(x) invlogit(a + b*x), size=1, color = "green") +
  labs (x ="Distance from hole", y = "Success rate")

#
# Geometric Model
# y_j ∼ binomial (n_j, p_j)
# p_j= 2* pbinom (arcsin((R−r)/x_j)* σ)−1, for j=1,…,J.  phi - cumulative normal distribution
r <- 1.68/2/12 # in feet
R <- 4.25/2/12 # in feet

golf_data2 <- c(golf_data1, r=r, R=R)

golf_model <- cmdstan_model ("~/Documents/RStanFiles/golf_angle.stan", pedantic = T)
fit_angle<- golf_model$sample (data = golf_data2, refresh = 0)
fit_angle$summary()

draws <- fit_angle$draws()
draws_df <- as_draws_df(draws)
sigma <- mean (draws_df$sigma) # in radians
sigma_deg <- mean (draws_df$sigma_degrees)
ggplot (golf, aes(x=x, y=p)) + 
  xlim (0.2, 20) + ylim (0, 1) +
  geom_pointrange(aes(x=x, ymin=p-se, ymax=p+se), color="blue", size = 0.2) +
  geom_function(fun = function(x) invlogit(a + b*x), size=1, color = "green") +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sigma))-1, size=1, color = "red") +
  annotate ("text", x=5, y=0.38, label="Geometry Based\nModel", color="red") + 
  annotate ("text", x=12, y=0.5, label="Logistic Regression",  color="green") + 
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nData, Logistic Regression and Geometry Based Model")

sig_deg <- c(0.5, 1, 2, 5, 20)
sig <- sig_deg/180*pi
labels <- paste0("\u03C3=", sig_deg, "\u00B0")
lab_x <- c(2.2, 3.5, 5.5, 7., 11.0)
lab_y <- c(0.2, 0.35, 0.5, 0.7, 0.8)
ggplot () +
  xlim (0.2, 20) + ylim (0, 1) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sigma))-1, size=1, color = "red") +
  annotate("text", x = 12.5, y = 0.3, label = paste("\u03C3=", round(sigma_deg,2), "\u00B0"), color="red")  +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sig[1]))-1) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sig[2]))-1) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sig[3]))-1) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sig[4]))-1) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sig[5]))-1) +
  annotate("text", x = lab_x, y = lab_y, label = labels)  +
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nModeled Pr(\"success\") for different values of \u03C3") 
#
# Fitted model vs new data
#
ggplot (data=golf, aes(x=x, y=p)) + 
  xlim (0.2, 80) + ylim (0, 1) +
  geom_pointrange(data=golf, aes(x=x, ymin=p-se, ymax=p+se, color="blue"), size = 0.2) +
  geom_pointrange(data=golf2, aes(x=x, ymin=p-se, ymax=p+se, color="red"), size = 0.2) +
  geom_function(fun = function(x) 2*pnorm ((asin((R-r)/x)/sigma))-1, color = "red") +
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nOld and New Data, Geometry Based Model") +
  scale_color_identity(name = "",
                        breaks = c("blue", "red"),
                        labels = c("Old Data", "New Data"),
                        guide = "legend") 
#
# In addition to uncertainty as to the angle, 
# There's also uncertainty as to how hard the ball is hit
# golfer hits the ball so it will travel "overshot" feet further than the hole
# u is the shot's potential distance u = (x+1) * (1 + error)
# define sigma_distance such that, u = (x+1, (x+1)*sigma_distance)  - multiplicative
# the distance is acceptable id u = [x, x+distance_tolerance]
# the probability is: pnorm (2/((x+overshot)*sigma_distance)) - 
#                     pnorm (-overshot/((x+overshot)*sigma_distance))
# so, the prob the shot goes in becomes:
#  ( 2 * pnorm (asin((R−r)/x)sigma_angle)−1) * 
#       (pnorm ((distance_tolerance - overshot)/((x+overshot)*sigma_distance)) - 
#        pnorm (-overshot/((x+overshot)*sigma_distance)))

overshot <- 1
distance_tolerance <- 3

golf_data_new <- list (x=golf2$x, y=golf2$y, n=golf2$n, J=dim(golf2)[1])  
golf_data_new <- c(golf_data_new, r=r, R=R, overshot=overshot, distance_tolerance=distance_tolerance)

golf_model <- cmdstan_model ("~/Documents/RStanFiles/golf_angle_distance.stan", pedantic = T)
fit_angle_distance<- golf_model$sample (data = golf_data_new, refresh = 0)
fit_angle_distance$summary()

draws <- fit_angle_distance$draws()
draws_df <- as_draws_df(draws)
sigma_angle <- mean (draws_df$sigma_angle) # in radians
sigma_deg <- mean (draws_df$sigma_degrees)
sigma_distance <- mean (draws_df$sigma_distance)

#sigma_angle <- 0.013

ggplot (data=golf2, aes(x=x, y=p)) + 
  xlim (0.2, 80) + ylim (0, 1) +
  geom_pointrange(data=golf2, aes(x=x, ymin=p-se, ymax=p+se), color="red", size = 0.2) +
  geom_function(fun = function(x) (2*pnorm ((asin((R-r)/x)/sigma_angle))-1) *
                  (pnorm ((distance_tolerance - overshot)/((x+overshot)*sigma_distance)) - 
                   pnorm (-overshot/((x+overshot)*sigma_distance))), 
                color = "red") +
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nGeometry and Distance Based omial BinModel") 

# Model unstable, Rhat's around 2.  more data at short end.  Renormlize n_j same, say 1000
golf3 <- golf2
norm_n <- 10
golf3$y <- round(golf2$y * norm_n/golf2$n)
golf3$n <- norm_n
golf_data4 <- list (x=golf3$x, y=golf3$y, n=golf3$n, J=dim(golf3)[1])  
golf_data4 <- c(golf_data4, r=r, R=R, overshot=overshot, distance_tolerance=distance_tolerance)

golf_model <- cmdstan_model ("~/Documents/RStanFiles/golf_angle_distance.stan", pedantic = T)
fit_angle_distance<- golf_model$sample (data = golf_data4, refresh = 0)
fit_angle_distance$summary()

draws <- fit_angle_distance$draws()
draws_df <- as_draws_df(draws)

sigma_angle <- mean (draws_df$sigma_angle) # in radians
sigma_deg <- mean (draws_df$sigma_degrees)
sigma_distance <- mean (draws_df$sigma_distance)

ggplot (data=golf3, aes(x=x, y=p)) + 
  xlim (0.2, 80) + ylim (0, 1) +
  geom_pointrange(data=golf3, aes(x=x, ymin=p-se, ymax=p+se), color="red", size = 0.2) +
  geom_function(fun = function(x) (2*pnorm ((asin((R-r)/x)/sigma_angle))-1) *
                  (pnorm ((distance_tolerance - overshot)/((x+overshot)*sigma_distance)) - 
                     pnorm (-overshot/((x+overshot)*sigma_distance))), 
                color = "red") +
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nFinal Geometry and Distance Based (normalized) Binomial Model") 

# 
# Change model from binomial to normal.  From:
# y_j ∼ binomial (n_j, p_j), to
# y_j/n_j ∼ normal (p_j, sqrt (p_j * (1−p_j)/n_j + sigma_y^2))
#
golf_data <- list (x=golf2$x, y=golf2$y, n=golf2$n, J=dim(golf2)[1])  
golf_data <- c(golf_data, r=r, R=R, overshot=overshot, distance_tolerance=distance_tolerance)

golf_model <- cmdstan_model ("~/Documents/RStanFiles/golf_angle_distance_f.stan", pedantic = T)
fit_angle_distance_f<- golf_model$sample (data = golf_data_new, refresh = 0)
fit_angle_distance_f$summary()

draws          <- fit_angle_distance_f$draws()
draws_df       <- as_draws_df(draws)
sigma_angle    <- mean (draws_df$sigma_angle) # in radians
sigma_deg      <- mean (draws_df$sigma_degrees)
sigma_distance <- mean (draws_df$sigma_distance)
sigma_y        <- mean (draws_df$sigma_y)

ggplot (data=golf2, aes(x=x, y=p)) + 
  xlim (0.2, 80) + ylim (0, 1) +
  geom_pointrange(data=golf2, aes(x=x, ymin=p-se, ymax=p+se), color="red", size = 0.2) +
  geom_function(fun = function(x) (2*pnorm ((asin((R-r)/x)/sigma_angle))-1) *
                  (pnorm ((distance_tolerance - overshot)/((x+overshot)*sigma_distance)) - 
                     pnorm (-overshot/((x+overshot)*sigma_distance))), 
                color = "red") +
  labs (x ="Distance from hole", y = "Success rate", 
        title = "Putts in pro golf:\nFinal: Geometry and Distance Based Normal Model") 

golf2$residuals <- (fit_angle_distance_f$summary("residuals")[2])$mean

ggplot (data=golf2, aes(x=x, y=residuals)) + geom_line () +
  labs (x="Distance from hole (feet)", y="y/n - fitted E(y/n", 
        title="Residuals from fitted model")

# 
# Lecture on 3/17
#
# needed in designing an experiment:
# 1. Treatment assignment rule
# 2. Number of people who get T and C
# 3. Outcome measure!
# 4. Only 2 treatments?
# 5. What are the treatments?
# 6. How is the treatment applied? Flavor of treatment. Confounding variables.
# 7. Pre-treatment variables ("covariates") - "experimental units" and "observational units"
# 8. Who is the experiment being done on?  What is the "sample"?
# 9. What is the "population" of interest?
# 10. Goals and constraints!

# lecture on 3/28
# http://www.stat.columbia.edu/~gelman/research/published/chung_etal_Pmetrika2013.pdf


# Chapter 15 Movie Ratings
#
# MODEL I: 2 movies.
#
# 102 total ratings, 2 for one movie, 100 ratings for the other
# theta_j, true popularity of two movies, j=1, 2
# j[1] = j[2]= 1, j[3]=...=j[102]=2
# yi - individual rating  i=1,...,102
# yi~ normal (theta_j[i], sigma), i=1,...,102

# ratings_1.stan
y_1 <- c(3, 5)
y_2 <- rep (c(2, 3, 4, 5), c(10, 20, 30, 40))
y <- c(y_1, y_2)
N <- length (y)
movie <- rep (c(1, 2), c(length (y_1), length(y_2)))
movie_data <- list (y=y, N=N, movie=movie)
#
ratings_1 <- cmdstan_model ("~/Documents/RStanFiles/ratings_1.stan", pedantic = T)
ratings_fit1 <- ratings_1$sample(data = movie_data, parallel_chains = 4, refresh = 0)

ratings_fit1$summary()
mcmc_hist(ratings_fit1$draws(c("theta")))

# MODEL II:  J movies
#
J <- 40 # 40 movies
# assign each movie a random number of reviewers from 1 to 100
N_ratings <- sample (1:100, J, replace = TRUE)
# create vector
N <- sum (N_ratings)
movie <- rep (1:J, N_ratings)
theta <- rnorm (J, 3, 0.5)
y <- rnorm (N, theta[movie], 2)
movie_data2 <- list (y=y, N=N, J=J, movie=movie)

ratings_2 <- cmdstan_model ("~/Documents/RStanFiles/ratings_2.stan", pedantic = T)
ratings_fit2 <- ratings_2$sample(data = movie_data2, parallel_chains = 4, refresh = 0)

ratings_fit2$summary()
mcmc_hist(ratings_fit2$draws(c("theta")))

draws <- ratings_fit2$draws()
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)

theta_bar_post <- rep (NA, J)
theta_se_post  <- rep (NA, J)
for (j in 1:J){
  theta_bar_post[j] <- mean(draws[,,j+1])
  theta_se_post[j] <-  sd(draws[,,j+1])
}

J_movies <- data.frame (theta_prior=theta, post=theta_bar_post, se=theta_se_post, N_ratings = N_ratings)

ggplot (data = J_movies) +
  geom_abline(intercept = 0, slope = 1, size=0.5, color = "gray")+
  geom_segment (aes(x=post-2*se, y=theta, xend=post+2*se, yend=theta), size=0.5, color="gray") +
  geom_segment (aes(x=post-0.66*se, y=theta, xend=post+0.66*se, yend=theta), size=0.5, color="black") +
  geom_point (aes (x=post, y=theta)) +
  labs (x="Posterior median, 50%, 95% interval", 
        y="True parameter value", 
        title="Comparing thetas to posterior in ferences") +
  xlim (1, 5) + ylim (1, 5)

ggplot (data = J_movies) +
  geom_point(aes (x=N_ratings, y=4*se/3)) +
  labs (x="Number of reviews", y="Width of 50% posterior interval", title = "More reviews, less uncertainity") + 
  ylim (0, 1.2)

#
# MODEL III:  Item-response
#
# y_i - rating i of movie j[i] by rater k[i]
# y_i ~ normal (a_j[i] - b_k[i], sigma_y)
# a_j - represents quality of movie j
# b_k - respresents toughness of rater k
# Model:
# y_i ~ normal (a_j[i] - b_k[i], sigma_y), i=1, ...N
# a_j ~ normal (mu_a, sigma_a),            j=1, ...J
# b_k ~ normal (0, sigma_k),               k=1, ...K
#
# REPARAMETERIZE:
# y_i ~ normal (mu + sigma_a * alpha_j[i] - sigma_b * beta_k[i], sigma_y), i=1, ...N)
# a_j ~ normal (0, 1),            j=1, ...J
# b_k ~ normal (0, 1),            k=1, ...K
J <- 40   # movies
K <- 100  # raters
N <- J*K  # movies
movie <- rep (1:J, rep (K, J))
rater <- rep (1:K, J)
mu <- 3
sigma_a <- 0.5
sigma_b <- 0.5
sigma_y <- 2
alpha <- rnorm (J, 0, 1)
beta  <- rnorm (K, 0, 1)
y <- rnorm (N, mu + sigma_a * alpha[movie] - sigma_b * beta[rater], sigma_y)
data_3 <- list (N=N, J=J, K=K, movie=movie, rater=rater, y=y)

ratings_3 <- cmdstan_model ("~/Documents/RStanFiles/ratings_3.stan", pedantic = T)
ratings_fit3 <- ratings_3$sample(data = data_3, parallel_chains = 4, refresh = 0)

ratings_fit3$summary(c("mu", "sigma_a", "sigma_b", "sigma_y"))
mcmc_hist(ratings_fit3$draws(c("mu", "sigma_a", "sigma_b", "sigma_y")))

draws <- ratings_fit3$draws()
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)

# Plot alphas
alpha_bar_post <- rep (NA, J)
alpha_se_post  <- rep (NA, J)
for (j in 1:J){
  alpha_bar_post[j] <- mean(draws[,,j+1])
  alpha_se_post[j] <-  sd(draws[,,j+1])
}

response_model <- data.frame (alpha_prior=alpha, post=alpha_bar_post, se=alpha_se_post)

ggplot (data = response_model) +
  geom_abline(intercept = 0, slope = 1, size=0.5, color = "gray")+
  geom_segment (aes(x=post-2*se, y=alpha, xend=post+2*se, yend=alpha), size=0.5, color="gray") +
  geom_segment (aes(x=post-0.66*se, y=alpha, xend=post+0.66*se, yend=alpha), size=0.5, color="black") +
  geom_point (aes (x=post, y=alpha)) +
  labs (x="Posterior median, 50%, 95% interval", 
        y="True parameter value", 
        title="Comparing alphas to posterior inferences\nResponse Model") +
  xlim (-3, 3) + ylim (-3, 3)

# Plot betas
beta_bar_post <- rep (NA, K)
beta_se_post  <- rep (NA, K)
for (k in 1:K){
  beta_bar_post[k] <- mean(draws[,,k+1+J])
  beta_se_post[k] <-  sd(draws[,,k+1+J])
}

response_model <- data.frame (beta_prior=beta, post=beta_bar_post, se=beta_se_post)

ggplot (data = response_model) +
  geom_abline(intercept = 0, slope = 1, size=0.5, color = "gray")+
  geom_segment (aes(x=post-2*se, y=beta, xend=post+2*se, yend=beta), size=0.5, color="gray") +
  geom_segment (aes(x=post-0.66*se, y=beta, xend=post+0.66*se, yend=beta), size=0.5, color="black") +
  geom_point (aes (x=post, y=beta)) +
  labs (x="Posterior median, 50%, 95% interval", 
        y="True parameter value", 
        title="Comparing betas to posterior inferences\nResponse Model") +
  xlim (-3.5, 3.5) + ylim (-3.5, 3.5)

# UNBALANCED DATA
genre <- rep (c("romantic", "crime"), c(round (J/2), J-round (J/2)) )

prob_of_rated <- ifelse (beta[rater] > 0,
                         ifelse (genre[movie] == "romantic", 0.2, 0.7),
                         ifelse (genre[movie] == "crime",    0.7, 0.2))
rated <- rbinom (N, 1, prob_of_rated) == 1 

data_3a <- list (N=sum (rated), J=J, K=K, movie=movie[rated], rater=rater[rated], y=y[rated])
ratings_fit3a <- ratings_3$sample(data = data_3a, parallel_chains = 4, refresh = 0)
ratings_fit3a$summary(c("mu", "sigma_a", "sigma_b", "sigma_y"))
mcmc_hist(ratings_fit3$draws(c("mu", "sigma_a", "sigma_b", "sigma_y")))

draws <- ratings_fit3a$draws()
draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
print(draws_df)
# Plot alphas
alpha_bar_post <- rep (NA, J)
alpha_se_post  <- rep (NA, J)
for (j in 1:J){
  alpha_bar_post[j] <- mean(draws[,,j+1])
  alpha_se_post[j] <-  sd(draws[,,j+1])
}

unbal_model <- data.frame (alpha_prior=alpha, post=alpha_bar_post, se=alpha_se_post, genre=genre)

ggplot (data = unbal_model) +
  geom_abline(intercept = 0, slope = 1, size=0.5, color = "gray")+
  geom_segment (aes(x=post-2*se, y=alpha, xend=post+2*se, yend=alpha), size=0.5, color="gray") +
  geom_segment (aes(x=post-0.66*se, y=alpha, xend=post+0.66*se, yend=alpha), size=0.5, color="black") +
  geom_point (aes (x=post, y=alpha, color=genre)) +
  labs (x="Posterior median, 50%, 95% interval", 
        y="True parameter value", 
        title="Comparing alphas to posterior inferences\nUnbalanced Model") +
  xlim (-3, 3) + ylim (-3, 3)

# Plot betas
beta_bar_post <- rep (NA, K)
beta_se_post  <- rep (NA, K)
for (k in 1:K){
  beta_bar_post[k] <- mean(draws[,,k+1+J])
  beta_se_post[k] <-  sd(draws[,,k+1+J])
}

unbal_model <- data.frame (beta_prior=beta, post=beta_bar_post, se=beta_se_post, 
                              rater=ifelse (beta > 0, "easy", "tough"))

ggplot (data = unbal_model) +
  geom_abline(intercept = 0, slope = 1, size=0.5, color = "gray")+
  geom_segment (aes(x=post-2*se, y=beta, xend=post+2*se, yend=beta), size=0.5, color="gray") +
  geom_segment (aes(x=post-0.66*se, y=beta, xend=post+0.66*se, yend=beta), size=0.5, color="black") +
  geom_point (aes (x=post, y=beta, color=rater)) +
  labs (x="Posterior median, 50%, 95% interval", 
        y="True parameter value", 
        title="Comparing betas to posterior inferences\nUnbalanced Model") +
  xlim (-3.5, 3.5) + ylim (-3.5, 3.5)

# true value of movie j
a <- mu + sigma_a * alpha
# average observed rating for movie j
y_bar <- rep (NA, J)
for (j in 1:J){
  y_bar[j] <- mean (y[movie==j & rated])
}

a_model <- (ratings_fit3a$summary(c("a"))[2])$mean

movie_q <- data.frame (a=a, y_bar=y_bar, a_model=a_model, genre=genre)

p1 <- ggplot (data=movie_q) +
  geom_point (aes(x=y_bar, y=a, color=genre)) +
  geom_abline(intercept = 0, slope = 1) +
  labs (x="Raw average rating for movie j", y="a_j") +
  xlim (1.5, 4.5) + ylim (1.5, 4.5) 

p2 <- ggplot (data=movie_q) +
  geom_point (aes(x=a_model, y=a, color=genre)) +
  geom_abline(intercept = 0, slope = 1) +
  labs (x="Posterior median estimate for movie j", y="a_j") +
  xlim (1.5, 4.5) + ylim (1.5, 4.5)

grid.arrange(p1, p2)

# Bayesian Workflow: Chapter 17
#
series <- matrix(scan("~/Downloads/Series1000.txt"), nrow=1000, ncol=135, byrow=TRUE)

T <- 135
N <- 1000
par(mar=c(3,3,2,0), tck=-.01, mgp=c(1.5,.5,0))
plot(c(1,T), range(series), bty="l", type="n", xlab="Time", ylab="series")
for (n in 1:N) {
  lines(1:T, series[n,], lwd=.5)
}
# IN ggplot
aa <- data.frame (t(rbind (series, x=1:T)))
aa_long <- aa %>%  pivot_longer (-x, values_to = "series")
# OR
aa_long <- data.frame (name=rep (1:N, rep (T, N)), x=rep (1:T, N), series=as.vector(t(series)))
ggplot (data=aa_long) +
  geom_line (aes (x=x, y=series, group=name), size=0.1)

# Calculate slope of each N curve
slope <- rep (NA, N)
se    <- rep (NA, N)
for (n in 1:N){
  time <- 1:T
  data <- series[n, ]
  fit <- lm (data ~ time)
  slope[n] <- 100*coef(fit)[2]
  se[n]    <- 100*se.coef(fit)[2]
}

ggplot (data = data.frame (Slope=slope, SE=se)) +
  geom_point (aes (x=Slope, y=SE))

ggplot (data = data.frame (Slope=slope, SE=se)) +
  geom_histogram (aes(x=Slope), binwidth=0.1, fill="blue", color="white") 
  
# Fit mixture model
y <- slope
K <- 3
N <- 1000
mu <- c(0, -1, 1)
data <- list (y=y, K=K, mu=mu, N=N)

mix <- cmdstan_model ("~/Documents/RStanFiles/mixture.stan", pedantic = T)
mix_fit <- mix$sample(data = data, parallel_chains = 4, refresh = 0)

mix_fit$summary(c("theta", "sigma"))
mcmc_hist(mix_fit$draws(c("theta", "sigma")))

#draws <- mix_fit$draws()
#draws_df <- as_draws_df(draws) # as_draws_matrix() for matrix
#summarize_draws (draws_df) # same as mix_fit$summary()

pp <- mix_fit$summary("p", "mean")
p_hat_v <- pp$mean
# OR
draws_p <- (mix_fit$draws("p"))
p_hat_v <- apply (draws_p, 3, mean)

p_hat <- matrix (p_hat_v, byrow = F, nrow=1000, ncol=3)
print (round(p_hat, 2))
  
max_prob <- apply (p_hat, 1, max)
choice   <- apply (p_hat, 1, which.max)
table (choice)

expected_correct <- sum(max_prob)
sd_correct <- sqrt(sum(max_prob*(1-max_prob)))
pfround (c(expected_correct, sd_correct), 1)

#$$$$$$$$$$$$$$$
# BDA Appendix C.3 Direct simulation
# theta_j ~ normal (mu, tau)
school <- c(1:8)
y =     c(28,  8, -3,  7, -1,  1, 18, 12)
sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
J <- length (y)

# Marginal and conditional simulation for the normal model
# Estimates
mu_hat <- function (tau, y, sigma){
  sum (y/(sigma^2 + tau^2))/sum (1/(sigma^2 + tau^2))
}
V_mu <- function (tau, y, sigma){
  1/sum (1/(sigma^2 + tau^2))
}

n_grid <- 2000
tau_grid <- seq (0.1, 40, length=n_grid)
log_p_tau <- rep (NA, n_grid)  # log posterior density of tau
for (i in 1:n_grid){
  mu <- mu_hat (tau_grid[i], y, sigma)
  V  <- V_mu (tau_grid[i], y, sigma)
  log_p_tau[i] <- 0.5*log (V) - 0.5*sum(log (sigma^2 + tau_grid[i]^2)) -
    0.5 * sum ((y-mu)^2 / (sigma^2 + tau_grid[i]))
}

log_p_tau <- log_p_tau - max (log_p_tau)
p_tau <- exp (log_p_tau)
p_tau <- p_tau/sum (p_tau)
n_sims <- 1000
tau <- sample (tau_grid, n_sims, replace = T, prob=p_tau)

mu <- rep (NA, n_sims)
theta <- array (NA, c(n_sims, J))
for (i in 1:n_sims){
  mu[i] <- rnorm (1, mu_hat (tau[i], y, sigma), sqrt (V_mu (tau[i], y, sigma)))
  theta_mean <- (mu[i]/tau[i]^2 + y/sigma^2)/(1/tau[i]^2 + 1/sigma^2)
  theta_sd <-   sqrt(1/sum (1/(tau[i]^2 + sigma^2)))
  theta[i,] <- rnorm (J, theta_mean, theta_sd)
}

data.frame (param=c(paste ("theta[", 1:8, "]", sep=""), "mu", "tau"), 
            MEAN=c(colMeans(theta), mean(mu), mean(tau)),
            SD=    c(apply(theta, 2, sd), sd(mu), sd(tau)))
            
################
# Gibbs Sampler for normal model
theta_update <- function (){
  theta_hat <- (mu/tau^2 + y/sigma^2)/(1/tau^2 + 1/sigma^2)
  V_theta <- 1/(1/tau^2 + 1/sigma^2)
  rnorm (J, theta_hat, sqrt(V_theta))
}
mu_update <- function (){
  rnorm (1, mean (theta), tau/sqrt(J))
}
tau_update <- function (){
  sqrt (sum ((theta-mu)^2)/rchisq(1, J-1))
}

chains <- 5
iter <- 1000
sims <- array (NA, c(iter, chains, J+2))
dimnames (sims) <- list (NULL, NULL, 
                         c(paste ("theta[", 1:8, "]", sep=""), "mu", "tau"))
for (m in 1:chains){
  mu <- rnorm(1, mean(y), sd(y))
  tau <- runif (1, 0, sd(y))
  for (t in 1:iter){
    theta <- theta_update()
    mu <- mu_update()
    tau <- tau_update()
    sims[t, m, ] <- c(theta, mu, tau)
  }
}
monitor (sims)

#
# Gibbs Sampler: Alternative parameterization


#############
# Hamiltonian Monte Carlo
# 10 parameters, th: mu, tau, theta[1,..8]
log_p_th <- function (th, y, sigma){  # log posterior density
  J <- length (th) - 2
  theta <- th[1:J]
  mu <- th[J+1]
  tau <- th[J+2]
  if (is.nan (tau) | tau <= 0)
    return (-Inf)
  else{
    log_hyperprior <- 1
    log_prior <- sum (dnorm (theta, mu, tau, log=TRUE))
    log_likelihood <- sum (dnorm(y, theta, sigma, log=TRUE))
    return (log_hyperprior + log_prior + log_likelihood)
  }
} 
# Gradients
gradient_th <- function (th, y, sigma){
  J <- length (th) - 2
  theta <- th[1:J]
  mu <- th[J+1]
  tau <- th[J+2]
  if (tau <= 0)
    return (c(0, 0, 0))
  else{
    d_theta <- - (theta-y)/sigma^2 - (theta-mu)/tau^2
    d_mu <- -sum(mu-theta)/tau^2
    d_tau <- -J/tau + sum ((mu-theta)^2)/tau^3
    return (c(d_theta, d_mu, d_tau))
  }   
}
# Numerical derivative - just in case
gradient_th_numerical <- function (th, y, sigma){
  d <- length (th)
  e <- 0.0001
  diff <- rep (NA, d)
  for (k in 1:d){
    th_hi <- th
    th_lo <- th
    th_hi[k] <- th[k] + e
    th_lo[k] <- th[k] - e
    diff[k]  <- (log_p_th (th_hi, y, sigma) - log_p_th (th_lo, y, sigma))/(2*e)
  }
  return (diff)
}
# One HMC iteration
hmc_iteration <- function (th, y, sigma, epsilon, L, M){
  M_inv <- 1/M
  d <- length (th)
  phi <- rnorm (d, 0, sqrt(M)) 
  th_old <- th
  log_p_old <- log_p_th (th, y, sigma) - 0.5*sum(M_inv * phi^2)
  phi <- phi +0.5* epsilon*gradient_th(th, y, sigma)
  for (l in 1:L){
    th <- th + epsilon*M_inv*phi
    phi <- phi + (if (l==L) 0.5 else 1)*epsilon*gradient_th(th, y, sigma)
  }
  phi <- -phi
  log_p_star <- log_p_th(th, y, sigma) - 0.5*sum (M_inv * phi^2)
  r <- exp (log_p_star - log_p_old)
  if (is.nan(r)) r <-0
  p_jump <- min (r, 1)
  th_new <- if (runif (1) < p_jump) th else th_old
  return (list (th=th_new, p_jump=p_jump))
}
# Wrapper
hmc_run <- function (starting_values, iter, epsilon_0, L_0, M){
  chains <- nrow (starting_values)
  d <- ncol (starting_values)
  sims <- array (NA, c(iter, chains, d),
                 dimnames = list(NULL, NULL, colnames (starting_values)))
  warmup <- 0.5*iter
  p_jump <- array (NA, c(iter, chains))
  for (j in 1:chains){
    th <- starting_values[j,]
    for (t in 1:iter){
      epsilon <- runif (1, 0, 2*epsilon_0)
      L <- ceiling (2*L_0*runif(1))
      temp <- hmc_iteration(th, y, sigma, epsilon, L, M)
      p_jump[t, j] <- temp$p_jump
      sims[t, j, ] <- temp$th
      th <- temp$th
    }
  }
  monitor (sims, warmup)
  cat ("Average acceptable probs:",
       fround (colMeans(p_jump[(warmup+1):iter,]), 2), "\n")
  return (list (sims=sims, p_jump=p_jump))
}

parameter_names <- c(paste("theta[", 1:8, "]", sep = ""), "mu", "tau")
d <- length (parameter_names)
chains <- 4

mass_vector <- rep (1/15^2, d)
starts <- array (NA, c(chains, d), dimnames=list(NULL, parameter_names))
for (j in 1:chains){
  starts[j, ] <- rnorm (d, 0, 15)
  starts[j, 10] <- runif (1, 0, 15)
}

M1 <- hmc_run (starting_values = starts, iter = 1000, 
               epsilon_0 = 0.05, L_0 = 20, M = mass_vector )

# Gibbs Sampler: Student t model with fixed degrees of freedom
#
theta_update <- function (){
  theta_hat <- (mu/V + y/sigma^2)/(1/V + 1/sigma^2)
  V_theta <- 1/(1/V + 1/sigma^2)
  rnorm (J, theta_hat, sqrt(V_theta))
}
mu_update <- function (){
  mu_hat <- sum (theta/V)/sum(1/V)
  V_mu <- 1/sum(1/V)
  rnorm (1, mean (mu_hat), sqrt(V_mu))
}
tau_update <- function (){
  sqrt (rgamma (1, J*nu/2 + 1, (nu/2)*sum(1/V)))
}
V_update <- function (){  # individual school variances
  (nu*tau^2 + (theta-mu)^2)/rchisq(J, nu+1)
}

chains <- 5
iter <- 1000
sims <- array (NA, c(iter, chains, J+2))
dimnames (sims) <- list (NULL, NULL, 
                         c(paste ("theta[", 1:8, "]", sep=""), "mu", "tau"))
nu <- 4
for (m in 1:chains){
  mu <- rnorm(1, mean(y), sd(y))
  tau <- runif (1, 0, sd(y))
  for (t in 1:iter){
    theta <- theta_update()
    V <- V_update ()
    mu <- mu_update()
    tau <- tau_update()
    sims[t, m, ] <- c(theta, mu, tau)
  }
}
monitor (sims)

# Gibbs-Metropolis Sampling: Student t model with variable degrees of freedom
#
# log_post - calculate the logarithm of the posterior distribution of 1/nu
# given all the other parameters
log_post <- function (theta, V, mu, tau, nu, y, sigma){
  sum (dnorm (y, theta, sigma, log = TRUE)) +      # data points, y_j
    sum (dnorm (theta, mu, sqrt(V), log = TRUE)) + # school effects, theta_j
    sum (0.5*nu*log(nu/2) + nu*log(tau) -         # variances, V_j
         lgamma (nu/2) - (nu/2+1)*log(V) - 0.5*nu*tau^2/V)
}

nu_update <- function (sigma_jump_nu){
  nu_inv_star <- rnorm (1, 1/nu, sigma_jump_nu)
  if (nu_inv_star <=0 | nu_inv_star >1)
    p_jump <- 0
  else{
    nu_star <- 1/nu_inv_star
    log_post_old  <- log_post(theta, V, mu, tau, nu, y, sigma) 
    log_post_star <- log_post(theta, V, mu, tau, nu_star, y, sigma)
    r <- exp (log_post_star - log_post_old)
    nu <- ifelse (runif (1) < r, nu_star, nu)
    p_jump <- min (r, 1)
  }
  return (list(nu=nu, p_jump=p_jump))
}

chains <- 5
iter <- 10000
sigma_jump_nu <- 1
p_jump_nu <- array (NA, c(iter, chains))
sims <- array (NA, c(iter, chains, J+3))
dimnames (sims) <- list (NULL, NULL, 
                         c(paste ("theta[", 1:8, "]", sep=""), "mu", "tau", "nu"))
for (m in 1:chains){
  mu <- rnorm(1, mean(y), sd(y))
  tau <- runif (1, 0, sd(y))
  V <- runif (J, 0, sd(y))^2
  nu < 1/runif (1, 0, 1)
  for (t in 1:iter){
    theta <- theta_update()
    V <- V_update ()
    mu <- mu_update()
    tau <- tau_update()
    temp <- nu_update (sigma_jump_nu)
    nu <- temp$nu
    p_jump_nu[t,m] <- temp$p_jump
    sims[t, m, ] <- c(theta, mu, tau, nu)
  }
}
monitor (sims)

good_sims <- sims[((iter/2)+1):iter, ,]

mix_sims <- data.frame(apply(good_sims, c(1,3), mean))
means_sims <- data.frame (means=colMeans(mix_sims))
sims_long <- mix_sims %>% pivot_longer (everything(), values_to = "vals") %>%
  mutate (means=rep(means_sims$means, nrow(mix_sims)))

means_df <- data.frame(name=colnames(mix_sims), means=means_sims$means, 
                    means_lab=paste0("mean=", round(means_sims$means,1)))

sims_long_theta <- sims_long %>% filter (name %in% c("mu", "nu", "tau")==F)
means_df_theta  <- means_df  %>% filter (name %in% c("mu", "nu", "tau")==F)
x_min <- floor(min(sims_long_theta$vals))
x_max <- ceiling(max(sims_long_theta$vals))  
ggplot (data=sims_long_theta, aes(x=vals, y=..density..)) +
  facet_wrap(~ name, ncol=4) +
  geom_histogram(binwidth=2, fill="blue", color="lightblue") +
  geom_segment (data=means_df_theta, aes(x=means, xend=means, y=0, yend = Inf), color="red") +
  geom_text (data=means_df_theta, aes(label=means_lab, x=x_min, y=0.1, hjust=0, vjust=0)) +
  xlim (x_min, x_max)  

sims_long_not_theta <- sims_long %>% filter (name %in% c("mu", "nu", "tau"))
means_df_not_theta  <- means_df  %>% filter (name %in% c("mu", "nu", "tau"))
x_min <- floor(min(sims_long_not_theta$vals))
x_max <- ceiling(max(sims_long_not_theta$vals))  
ggplot (data=sims_long_not_theta, aes(x=vals, y=..density..)) +
  facet_wrap(~ name, ncol=4) +
  geom_histogram(binwidth=1, fill="blue", color="lightblue") +
  geom_segment (data=means_df_not_theta, aes(x=means, xend=means, y=0, yend = Inf), color="red") +
  geom_text (data=means_df_not_theta, aes(label=means_lab, x=x_min, y=0.18, hjust=0, vjust=0)) +
  xlim (x_min, x_max)  


#############
# Lecture 4/14
#From Charles Margossian to Everyone: (10:36 AM)
# some recommended supplementary reading: chapter 1 - 7 from Hastie et al, https://web.stanford.edu/~hastie/ElemStatLearn/ and case studies by Betancourt https://betanalpha.github.io/writing/
#  From Andrew Gelman to Everyone: (10:46 AM)
# model is Pr(y=1) = invlogit(a + b*x) to get my guess:(i) suppose x=300 is the population avg, guess:  invlogit(a + 300b) = 0.5 thus:  a + 300b = 0 (ii) guess for b suppose that prob is 50% for a 3rd-generation immigrant but only 40% for a 4th+-generation invlogit(a + 400b) = 0.4 
#From Andrew Gelman to Everyone: (10:49 AM)
# a + 400b = 4*(0.4 - 0.5)
# a + 400b = -0.4
# so 100 b = 0.4
# b = +0.004
# a + 300 b = 0
# a = 300*0.004 = -1.2
# assume prior independence of (a +300b) and b.
# a + 300 b ~ normal(0, 1)
# b ~ normal(-0.004, 0.004)
# From Ling Chen to Everyone: (10:54 AM)
# Can we standardize X in the first place?
#  From Andrew Gelman to Everyone: (10:54 AM)
# x = (X - mean(X))/sd(X)
# From Andrew Gelman to Everyone: (10:55 AM)
# in stan:  a_centered = a + 300 b; (in transformed parameters)
# a_centered ~ normal(0, 1); (in model block)
# y ~ bernoulli_inv_logit(a + b*x);
# x_obs ~ normal(x, sigma_x_obs); // measurement error model
# x ~ normal(0, 1); // prior
#From Andrew Gelman to Everyone: (10:57 AM)
# Here, y and x_obs are in the data block
# and x is in the parameters block
#From Charles Margossian to Everyone: (11:05 AM)
# More on autodiff than you’d want to know: https://onlinelibrary.wiley.com/doi/10.1002/widm.1305 And this only scratches the surface.
# Alright, I can’t just refer you to what I wrote. Here’s another very good and popular review: https://arxiv.org/abs/1502.05767
# From Charles Margossian to Everyone: (11:06 AM)
# An if you want to enter a world of pain and suffering, the seminal reference is the book by Griewank and Walter 2008
# The Jax autodiff cookbook: https://jax.readthedocs.io/en/latest/notebooks/autodiff_cookbook.html
#From xintian to Everyone: (11:17 AM)
# This is exactly the answer I have in my solution!
#  From Andrew Gelman to Everyone: (11:21 AM)
# regarding latex and word, see the Technical Note of this paper:  http://www.stat.columbia.edu/~gelman/research/published/zombies.pdf
################
# BIOASSAY CH 3.7
#
# x: toxin amount
# num_c: number of cows per bin
# y: number of deaths per bin

log_x <- c(-0.86, -0.3, -0.05, 0.73)
num_c <- c(5, 5, 5, 5)
y     <- c(0, 1, 3, 5)

# https://cran.r-project.org/web/packages/rstanarm/vignettes/binomial.html
# create granular dead/alive data per cow per x

######
# bivariate normal density
mu_a <- 0
mu_b <- 10
sigma_a <- 2
sigma_b <- 10 
rho <- 0.5
Mu <- c(mu_a, mu_b)
Sigma <- matrix (c(sigma_a^2, sigma_a*sigma_b*rho, sigma_a*sigma_b*rho, sigma_b^2),
                 ncol=2)
###########
logxs <- rep (log_x, c(rep (5, 4)))
deads <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1)
data.df <- data.frame (logxs=logxs, deads=deads)
####### 
# aplha, beta estimates
fit.glm <- glm (formula = deads ~ logxs, data = data.df, family = binomial)
print (fit.glm, digits=2)
a_opt <- coef(fit.glm)[1]
b_opt <- coef(fit.glm)[2]

ggplot (data.frame(log_x=log_x, num_c=num_c, y=y)) + 
  geom_point (aes(x=log_x, y=y/num_c), color="red") + 
  geom_function(fun = function (x) invlogit (a_opt + b_opt*x)) +
  labs (x="log (treatment)", y="Pr(cow dies)", title="Averaged data vs Invlogit fit")

################ 
# Bayesian Inference

unif_d <- function (a_gid, b_grid){
  density <- array (1, c(length (a_grid), length (b_grid)))
  return (density/sum(density))
}

point_likelihood <- function (a, b){
  return (prod ((invlogit (a + b*log_x))^y *
          (1 - invlogit (a + b*log_x))^(num_c - y)))
}

bivar_d <- function (a, b){
  return (mvtnorm::dmvnorm(c(a, b), mean=Mu, sigma=Sigma))
}
#################

n_a <- 1000 # num grid points in a
n_b <- 1000 # num grid points in b
a_grid <- seq (-5, 10,  length.out=n_a)
b_grid <- seq (-10, 40, length.out=n_b)

# Likelihood density
start_time <- Sys.time()
likelihood <- outer(a_grid, b_grid, Vectorize(point_likelihood))
end_time <- Sys.time()
(elapsed <- end_time - start_time)
likelihood <- likelihood/sum (likelihood)

# 1) Uniform prior -> posterior
d_unif  <- unif_d  (a_grid, b_grid)
post_unif_d  <- d_unif*likelihood/sum (d_unif*likelihood)
all.equal(likelihood, post_unif_d)

# 2) Bivariate normal prior -> posterior
start_time <- Sys.time()
d_bivar <- outer (a_grid, b_grid, Vectorize(bivar_d))
end_time <- Sys.time()
(elapsed <- end_time - start_time)

d_bivar <- d_bivar/sum (d_bivar)
post_bivar_d <- d_bivar*likelihood/sum (d_bivar*likelihood)
##################
# Visualizations
# Sample from posterior density
n_sims <- 1000
a_sim <- rep (NA, n_sims)
b_sim <- rep (NA, n_sims)

# uniform prior
for (s in 1:n_sims){
  idx <- sample (1:n_a, size=1, replace = T, prob = rowSums (post_unif_d))
  a_sim[s] <- a_grid[idx]
  b_sim[s] <- sample (b_grid, size=1, replace = T, prob = post_unif_d[idx, ])
}

p_cont_unif <- ggplot (data.frame (alpha=a_sim, beta=b_sim)) + 
  geom_jitter (aes (x=alpha, y=beta), width=0.1, height=0.3, size=0.4, color="blue") + 
  geom_point (aes(x=a_opt, y=b_opt),color="red") + 
  xlim (-4, 10) + ylim (-10, 40) +
  labs (title="Posterior density of alpha and beta - uniform prior") +
  geom_density_2d_filled (aes (x=alpha, y=beta), alpha=0.5, size=0.5, color="black")

#bivariate prior  
for (s in 1:n_sims){
  idx <- sample (1:n_a, size=1, replace = T, prob = rowSums (post_bivar_d))
  a_sim[s] <- a_grid[idx]
  b_sim[s] <- sample (b_grid, size=1, replace = T, prob = post_bivar_d[idx, ])
}

p_cont_bivar <- ggplot (data.frame (alpha=a_sim, beta=b_sim)) + 
  geom_jitter (aes (x=alpha, y=beta), width=0.1, height=0.3, size=0.4, color="blue") + 
  geom_point (aes(x=a_opt, y=b_opt),color="red") + 
  xlim (-4, 10) + ylim (-10, 40) +
  labs (title="Posterior density of alpha and beta - bivariate prior") +
  geom_density_2d_filled (aes (x=alpha, y=beta), alpha=0.5, size=0.5, color="black")

grid.arrange(p_cont_unif, p_cont_bivar)

# Countour graphs
par(mfrow=c(2,1))
contour (x=a_grid, y=b_grid, z=post_unif_d, lwd=2, col="red", 
         vfont = c("sans serif", "bold"),
         xlim = c(-2, 4), ylim =c(-0, 30))
contour (x=a_grid, y=b_grid, z=post_bivar_d, lwd=2, col="red", 
         vfont = c("sans serif", "bold"),
         xlim = c(-2, 4), ylim =c(-0, 30))
par(mfrow=c(1,1))

xx <- a_grid[rep (1:n_a, rep (n_b,n_a))]
yy <- b_grid[rep((1:n_b), n_a)]

p1 <- ggplot (data.frame (alpha=xx, beta=yy, zz=as.vector ((d_unif))), 
              aes(x=alpha, y=beta, z=zz)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Prior density: Uniform")

p2 <- ggplot (data.frame (alpha=xx, beta=yy, zz=as.vector ((d_bivar))), 
              aes(x=alpha, y=beta, z=zz)) +
  xlim (-5, 10) + ylim (-10, 30) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Prior density: bivariate, rho=0.5")

p3 <- ggplot (data.frame (alpha=xx, beta=yy, likelihood=as.vector (t(likelihood))), 
              aes(x=alpha, y=beta, z=likelihood)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Likelihood of alpha and beta")

p4 <- ggplot (data.frame (alpha=xx, beta=yy, post_dense=as.vector(t(post_bivar_d))), 
              aes(x=alpha, y=beta, z=post_dense)) +
    xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Posterior density: bivariate normal Prior")

p5 <- ggplot (data.frame (alpha=xx, beta=yy, post_dense=as.vector(t(post_unif_d))), 
              aes(x=alpha, y=beta, z=post_dense)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Posterior density: uniform prior")

grid.arrange(p1, p3, p5)
grid.arrange(p2, p3, p4)
##################
# 3D plot
col.pal <- colorRampPalette(c("red", "orange", "yellow"))
colors <- col.pal(100)
z <- post_bivar_d
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
z.facet.range <- cut(z.facet.center, 100)

persp (x=a_grid, y=b_grid, z=post_bivar_d,
       xlab="alpha", ylab="beta", zlab="Posterior density", 
       main="Posterior density of alpha, beta",
       phi = 5,
       theta = 30,
       col=colors[z.facet.range],
       border="grey20",
       ticktype = "detailed",
       axes = T)
###########
persp (a_grid, b_grid, d_unif, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, d_bivar, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, likelihood, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, post_unif_d, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, likelihood - post_unif_d, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, post_bivar_d, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, likelihood - post_bivar_d, theta=30, ticktype = "detailed")
persp (a_grid, b_grid, post_bivar_d - post_unif_d, theta=30, ticktype = "detailed")

######### 
#Bivariate PDFs, CDFs, plots
# https://cran.r-project.org/web/packages/bivariate/vignettes/bivariate.pdf
install.packages("barsurf")
install.packages("bivariate")
library (barsurf)
library (bivariate)
f <- nbvpdf (0, 0, 1, 1, 0) # mux, muy, sigmax, sigmay, corr
(f(0, 0))

F <- nbvcdf (0, 0, 1, 1, 0)
F(0, 0)
plot (f, theme="blue")
plot (f, FALSE, theme="blue")
plot (F, theme="blue")
################

################
# CORONAVIRUS MODEL CH 18
# https://www.medrxiv.org/content/10.1101/2020.05.22.20108944v3.full
# https://bob-carpenter.github.io/diagnostic-testing/reports/case-study.html
################
# Shortest Posterior Interval (SPIN)
spin <- function(x, lower=NULL, upper=NULL, conf=0.95){
  x <- sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n <- length(x)
  gap <- round(conf*n)
  width <- x[(gap+1):n] - x[1:(n-gap)]
  index <- min(which(width==min(width)))
  x[c(index, index + gap)]
}
##########
var_names <- c("pi","spec_unk", "sens_unk", "mu_logit_spec", "mu_logit_sens", 
                "sigma_logit_spec", "sigma_logit_sens")

print_spin_results <- function (draws_df, var_names){

  medians <- rep (NA, length (var_names))
  results <- array (NA, c(length (var_names), 2))

  for (i in 1:num_vars){
      medians[i] <- round(median(pull (draws_df, var_names[i])), 3)
      results[i, ] <- round (spin (pull (draws_df, var_names[i]), conf=0.95), 3)
  }
  intervals <- sprintf ("%5.3f - %5.3f", results[, 1], results[, 2])
  
  print (data_frame (variable=var_names, median=medians, "95% interval"=intervals),
         row.names = F)
  return (medians %>% cbind (results))
}


# y - test result positive/negative
# sens - sensitivity 
# spec - specificity
# pi- prevalence
# p - probability of positive
# 
# p = pi*sens + (1-pi)*(1-spec)
# pi = (p + spec - 1)/(sens + spec -1)
#
n_sample <- 3330 # number of tests in sample
y_sample <- 50    # number of positives
n_spec <- 401 # number of known negatives
y_spec <- 399 # negative results out of known negatives
n_sens <- 122 # number of known positives
y_sens <- 103 # positive results out of known positives

covid_dat <- list(n_sample=n_sample, y_sample=y_sample, 
                  n_spec=n_spec, y_spec=y_spec, 
                  n_sens=n_sens, y_sens=y_sens)
covid_model1 <- cmdstan_model("covid.stan", pedantic=F)
fit1 <- covid_model1$sample (data=covid_dat, refresh = 0)
fit1$summary()
fit1$summary(c("pi", "spec", "sens"))

draws <- fit1$draws()
draws_df <- as_draws_df(draws) 

p1 <- ggplot (draws_df, aes(x=spec, y=pi)) +
  geom_point (size=0.5) +
  labs (x="Specificity", y="prevalence")

p2 <- ggplot (draws_df) + 
  geom_histogram (aes(x=pi, y = (..count..)/sum(..count..)), 
                  fill="gray", color="black", binwidth =0.002, boundary=0 ) +
  labs (x="prevalence, pi", y="Density")

grid.arrange(p1, p2, ncol=2)
########
# Non-informative prior
covid2_dat <- list(K_pos=4,
                   N_pos=c(0, 85, 37, 35),
                   n_pos=c(0, 78, 27, 25),
                   K_neg=14,
                   N_neg=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
                   n_neg=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),                   
                   N_unk=3330,
                   n_unk=50)
covid2_dat$sigma_sigma_logit_sens <- 1
covid2_dat$sigma_sigma_logit_spec <- 1

covid_model2 <- cmdstan_model("covid2.stan", pedantic=F)
fit2 <- covid_model2$sample (data=covid2_dat, 
                             iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                             refresh = 0, adapt_delta=0.9) 

print (fit2, digits=3, var_names, c("median"))
#print (fit2, digits=3, c("p", "spec", "sens","sigma_logit_spec", "sigma_logit_sens"), max_rows = 30)

draws2 <- fit2$draws()
draws2_df <- as_draws_df(draws2) 
non_inf_result <- print_spin_results (draws2_df, var_names)
########
# Informative prior
covid2b_dat <- covid2_dat
covid2b_dat$sigma_sigma_logit_sens <- 0.3
covid2b_dat$sigma_sigma_logit_spec <- 0.3

fit2b <- covid_model2$sample (data=covid2b_dat, 
                              iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                              refresh = 0, adapt_delta=0.9) 
                              
fit2b$summary (var_names)
print (fit2b, var_names, digits=3)

draws2b <- fit2b$draws()
draws2b_df <- as_draws_df(draws2b) 
inf_result <- print_spin_results (draws2b_df, var_names)
################
spec_sens_plots <- function (data, draws_df, title){
  spec_sens_df <- draws_df %>% select (starts_with (c("sens[", "spec[")))
  
  lo_spin <- rep (NA, dim (spec_sens_df)[2])
  medians <- rep (NA, dim (spec_sens_df)[2])
  hi_spin <- rep (NA, dim (spec_sens_df)[2])
  for (i in 1:dim(spec_sens_df)[2]){
    aa <- spin (as.matrix(spec_sens_df)[,i], lower=0, upper=1, conf=0.95)
    lo_spin[i] <- aa[1]
    hi_spin[i] <- aa[2]
    medians[i] <- median (as.matrix(spec_sens_df)[,i])
  }
  names <- colnames(spec_sens_df)

  spec_sens_dat <- data.frame (orig=c(data$n_pos/data$N_pos, data$n_neg/data$N_neg), names, lo_spin, medians, hi_spin)

  spec_sens_dat$names <- factor (spec_sens_dat$names, levels = spec_sens_dat$names)

  p1 <- ggplot (spec_sens_dat %>% filter(str_detect(names, "spec"))) +
    geom_pointrange(aes(x=names, y=medians, ymin=lo_spin, ymax=hi_spin, color="blue")) +
    geom_point (aes(x=names, y=orig, color="red")) +
    labs (title=title,
          subtitle="Specificity: point and posterior estimates",
          x="test", y="estimate") +
    theme(axis.text.x = element_text(angle = 45)) +
    scale_color_identity(name = "",
                         breaks = c("red", "blue"),
                         labels = c("Prior point estimate", "Posterior median and spin"),
                         guide = "legend") 
  p2 <- ggplot (spec_sens_dat %>% filter(str_detect(names, "sens"))) +
    geom_pointrange(aes(x=names, y=medians, ymin=lo_spin, ymax=hi_spin, color="blue")) +
    geom_point (aes(x=names, y=orig, color="red")) +
    labs (title=title,
          subtitle="Sensitivity: point and posterior estimates",
          x="test", y="estimate") +
    theme(axis.text.x = element_text(angle = 45)) +
    scale_color_identity(name = "",
                         breaks = c("red", "blue"),
                         labels = c("Prior point estimate", "Posterior median and spin"),
                         guide = "legend") 
  
  
  spec_long <- spec_sens_df %>% 
    pivot_longer (everything(), values_to="vals") %>% 
    filter (str_detect(name, "spec"))
  sens_long <- spec_sens_df %>% 
    pivot_longer (everything(), values_to="vals") %>% 
    filter (str_detect(name, "sens"))

  spec_long$name <- factor (spec_long$name, levels = unique(spec_long$name))
  p3 <- ggplot (data=spec_long, aes(x=vals, y=..density..)) +
    facet_wrap(~ name, ncol=4) +
    geom_histogram(binwidth=0.01, fill="blue", color="lightblue") +
    labs (title=title,
          subtitle="Posterior specificity of test", x="Specifity")

  p4 <- ggplot (data=sens_long, aes(x=vals, y=..density..)) +
    facet_wrap(~ name, ncol=4) +
    geom_histogram(binwidth=0.05, fill="blue", color="lightblue") +
    labs (title=title,
          subtitle="Posterior sensitivity of test", x="Sensitivity")
  
  return (list (p1, p2, p3, p4))
}

non_inf_plots <- spec_sens_plots (covid2_dat, draws2_df, "Non-Informative Prior")
inf_plots     <- spec_sens_plots (covid2b_dat, draws2b_df, "Informative Prior")
################
# prevalence Visualizations
pi_df <- data.frame(noniformative=draws2_df$pi, informative=draws2b_df$pi) %>%
  pivot_longer (everything(), names_to="prior", values_to = "values")

ggplot (pi_df) +
  geom_density(aes(x=values, color=prior, fill=prior), size=1, alpha=0.5) +
  labs (title="COVID prevalence, pi,  non-infomative and informative priors",
        x="pi", y="Density(pi)") +
  xlim (0, 0.1) 
##################
# Explicit - K_UNK MODEL
# https://bob-carpenter.github.io/diagnostic-testing/reports/case-study.html
###################
covid3_dat <- list(K_pos=3,
                   N_pos=c(85, 37, 35),
                   n_pos=c(78, 27, 25),
                   K_neg=13,
                   N_neg=c(371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
                   n_neg=c(368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
                   N_unk=3330,
                   n_unk=50)
covid3_dat$sigma_sigma_logit_sens <- 1
covid3_dat$sigma_sigma_logit_spec <- 1

model <- cmdstan_model("prior2_sensitivity.stan")
fit3 <-  model$sample(data = covid3_dat, 
                     iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                     refresh = 0, adapt_delta=0.9)

print (fit3, max_rows=42, digits=3)
draws3 <- fit3$draws()
draws3_df <- as_draws_df(draws3) 

model3_result <- print_spin_results (draws3_df, var_names)
model3_plots  <- spec_sens_plots (covid3_dat, draws3_df, "Non-Informative")
#########
# Informative Explicit model
#########
covid3b_dat <- covid3_dat
covid3b_dat$sigma_sigma_logit_sens <- 0.3
covid3b_dat$sigma_sigma_logit_spec <- 0.3

model <- cmdstan_model("prior2_sensitivity.stan")
fit3b <-  model$sample(data = covid3b_dat, 
                      iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                      refresh = 0, adapt_delta=0.9)

print (fit3b, max_rows=42, digits=3)

print(fit3b, digits = 3, c("pi"))

draws3b <- fit3b$draws()
draws3b_df <- as_draws_df(draws3b) 

model3b_result <- print_spin_results (draws3b_df, var_names)
model3b_plots  <- spec_sens_plots (covid3b_dat, draws3_df, "Informative")
##########
# visualizations
pi_df <- data.frame(noniformative=draws3_df$pi, informative=draws3b_df$pi) %>%
  pivot_longer (everything(), names_to="prior", values_to = "values")

ggplot (pi_df) +
  geom_density(aes(x=values, color=prior, fill=prior), size=1, alpha=0.5) +
  labs (title="COVID prevalence, pi,  non-infomative and informative priors",
        x="pi", y="density(pi)") +
  xlim (0, 0.1) 

n3_unk <- rbinom (length (draws3_df$pi), 3330, 
                 draws3_df$pi * draws3_df$sens_unk + 
                   (1 - draws3_df$pi) * (1 - draws3_df$spec_unk))
n3b_unk <- rbinom (length (draws3b_df$pi), 3330, 
                   draws3b_df$pi * draws3b_df$sens_unk + 
                     (1 - draws3b_df$pi) * (1 - draws3b_df$spec_unk))

n_unk <- data.frame (n3_unk, n3b_unk)
ggplot(data=n_unk) +
  geom_density (aes(x=n3_unk, y=..density..), binwidth = 1, 
                  fill="blue", color="blue") +
  geom_segment (aes(x=median (n3_unk), y=0, xend=median (n3_unk), yend=Inf), 
                size=0.5, color="red", alpha = 0.5) +
  geom_density (aes(x=n3b_unk, y=..density..), binwidth = 1, 
                  fill="red", color="red", alpha = 0.5) +
  geom_segment (aes(x=median (n3b_unk), y=0, xend=median (n3b_unk), yend=Inf), size=0.5, color="blue") +
  labs (title = "COVID infection distribution in a population of 3330 people")

ggplot () +
  geom_density (data=draws3b_df, aes(x=(pi*sens_unk + (1 - pi) * (1 - spec_unk))), color="blue", alpha=0.5) +
  geom_density (data=draws3_df, aes(x=(pi*sens_unk + (1 - pi) * (1 - spec_unk))), color="red", alpha=0.5) +
  xlim (0, 0.1)

##################
# prior sensitivity
###################
ribbon_df <- data.frame(sigma_sens = c(), sigma_spec = c(),
                        prev05 = c(), prev50 = c(), prev95 = c())

model <- cmdstan_model("prior2_sensitivity.stan")

sigma_senss <- c(0.01, 0.25, 0.5, 0.75, 1)
sigma_specss <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for (sigma_sens in sigma_senss) {
  for (sigma_spec in sigma_specss) {
    print(c(sigma_sens, sigma_spec))
    covid3_dat$sigma_sigma_logit_sens <- sigma_sens
    covid3_dat$sigma_sigma_logit_spec <- sigma_spec
    fit <-  model$sample(data = covid3_dat, 
                         iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                         refresh = 0, adapt_delta=0.9)
    pis <- as.vector(fit$draws()[,,"pi"])
    ribbon_df <- rbind(ribbon_df,
                       data.frame(sigma_sens = paste("sensitivity hyperprior", "=", sigma_sens),
                                  sigma_spec = sigma_spec,
                                  prev05 = quantile(pis, 0.05),
                                  prev50 = quantile(pis, 0.5),
                                  prev95 = quantile(pis, 0.95)))
  }
}

# ribbon version
plot_ribbon <- ggplot(ribbon_df, aes(x = sigma_spec)) +
  facet_wrap(~ sigma_sens, nrow = 1) +
  geom_ribbon(aes(ymin = prev05, ymax = prev95), fill = "gray95") +
  geom_line(aes(y = prev50), size = 0.5) +
  geom_line(aes(y = prev05), color = "darkgray", size = 0.25) +
  geom_line(aes(y = prev95), color = "darkgray", size = 0.25) +
  scale_y_log10(limits = c(0.0009, 1.1), breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 1),
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  ylab("prevalence") +
  xlab("specificity hyperprior") +
  theme_bw() +
  theme(panel.spacing = unit(0.25, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_ribbon

###############
# Multilevel Regression Poststratification MRP Model
##############
N <- 3330
y <- sample (rep (c(0, 1), c(3330-50, 50)))
male <- sample (rep (c(0, 1), c(2101, 1229)))
eth <- sample(rep(1:4, c(2118, 623, 266, 306+17)))
age <- sample(rep(1:4, c(71, 550, 2542, 167)))
N_zip <- 58 # number of different zip codes
zip <- sample (1:N_zip, N, replace = T)
x_zip <- rnorm (N_zip, 50, 20) # Zip code level predictor

J <- 2*4*4*N_zip # 1856 cells
N_pop <- rep (1000, J) # assume each cell has 1000 people

mrp_data <- list(N=3330,
                 y=y,
                 male=male,
                 eth=eth,
                 age=age,
                 zip=zip,
                 N_zip=N_zip,
                 x_zip=x_zip,
                 J_spec=14,
                 y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
                 n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
                 J_sens=4,
                 y_sens=c(0, 78, 27, 25),
                 n_sens=c(0, 85, 37, 35),
                 logit_spec_prior_scale=0.3,
                 logit_sens_prior_scale=0.3,
                 coef_prior_scale=0.5,
                 J=J,
                 N_pop=N_pop)

model_MRP <- cmdstan_model("Hierarchical_MRP.stan")

fit_MRP <- model_MRP$sample(data = mrp_data, 
             iter_warmup = 10000, iter_sampling = 10000, seed = 1234, parallel_chains=4,
             refresh = 0, adapt_delta=0.9)
print (fit_MRP, c("p_avg"), digits=4)

draws_MRP <- fit_MRP$draws()
draws_MRP_df <- as_draws_df(draws_MRP) 

aa <- print(spin(draws_MRP_df$p_avg, lower=0, upper=1, conf=0.95))

var_MRP_names <- c("p_avg", "b", "a_age", "a_eth", 
                   "sigma_eth", "sigma_age", "sigma_zip", 
                   "mu_logit_spec", "sigma_logit_spec",  "mu_logit_sens", 
                   "sigma_logit_sens", "p_pop[1]", "p_pop[2]", "p_pop[3]")
print (fit_MRP, var_MRP_names, digits=3, max_rows = 23)

ggplot (draws_MRP_df) +
  geom_histogram (aes(x=p_avg, y=..density..), fill="blue", color="lightblue", 
                  binwidth=0.002, boundary=0, alpha=0.5) +
  geom_segment (aes(x=median (p_avg), y=0, xend=median (p_avg), yend=Inf), size=0.5, color="red") +
  geom_segment (aes(x=aa[1], y=0, xend=aa[1], yend=Inf), size=0.2, 
                color="red", linetype="dashed") +
  geom_segment (aes(x=aa[2], y=0, xend=aa[2], yend=Inf), size=0.2, 
                color="red", linetype="dashed") +
  labs (title="p_avg, COVID prevalence averaged over MRP cells\nMedian and 95% shortest posterior interval")
  
##############
# World Cup 2014
###############
# i - game index
# j - team index
# a_j - ability of team j
# z_j_i - goals scored by team j in game i
# prior_score_j - score for team j before tournament
# a_j ~ normal (a + b * prior_score_j, sigma_a)
# or
# a_j <- a + b * prior_score_j + sigma_a * alpha_j
# where alpha_j ~ normal (0, 1) for j=1, ..., J


# team rankings before games
teams <- as.vector(unlist(read.table("soccerpowerindex.txt", header=FALSE)))

N_teams <- length(teams)
prior_score <- rev(1:N_teams)
prior_score <- (prior_score - mean(prior_score)) / (2 * sd(prior_score))

data2014 <- read.table ("worldcup2014.txt", header=FALSE)
N_games <- nrow(data2014)
team_1 <- match(as.vector(data2014[[1]]), teams)
score_1 <- as.vector(data2014[[2]])
team_2 <- match(as.vector(data2014[[3]]), teams)
score_2 <- as.vector(data2014[[4]])
df <- 7
data_worldcup <- list (N_teams=N_teams, N_games=N_games, 
                       prior_score=prior_score, 
                       team_1=team_1, team_2=team_2, 
                       score_1=score_1, score_2=score_2, df=df)

model_0 <- cmdstan_model("world_cup_0.stan")
fit_0 <- model_0$sample(data = data_worldcup, refresh = 0)
print (fit_0, max_rows=130)

model_1 <- cmdstan_model("world_cup_1.stan")
fit_1 <- model_1$sample(data = data_worldcup, refresh = 0)
print (fit_1, max_rows=130)

draws <- fit_1$draws()
# draws <- fit_0$draws()
draws_y_rep_df <- as_draws_df(draws) %>% select(starts_with("y_rep["))
draws_a_df <- as_draws_df(draws) %>% select(starts_with("a["))

median <- apply(draws_a_df, 2, quantile, 0.5)
q25_a <- apply(draws_a_df, 2, quantile, 0.16)
q75_a <- apply(draws_a_df, 2, quantile, 0.84)

q25_y_rep <- apply(draws_y_rep_df, 2, quantile, 0.025)
q75_y_rep <- apply(draws_y_rep_df, 2, quantile, 0.975)

par(mar = rep(2, 4))

png("teams.png", 750, 500)
coefplot (rev(median), sds=rep(0, N_teams),
          lower.conf.bounds=rev(q25_a), upper.conf.bounds=rev(q75_a),
          varnames=rev(paste(teams)),
          main="Team quality (estimate +/- 1 s.e.)\n(model with no square root)",
          mar=c(2,7,10,2), xlim=c(-1.25, 1.25))
dev.off()

png("games.png", 750, 1000)
coefplot (rev(score_1 - score_2), sds=rep(0, N_games),
          lower.conf.bounds=rev(q25_y_rep), upper.conf.bounds=rev(q75_y_rep),
          varnames=varnames,
          main="Game score differentials\ncompared to 95% predictive interval from model\n(model with no square root)",
          mar=c(2,10,10,2))
dev.off()

#
# BIOASSAY 3-7 redo
# 
# x_i     - toxin level
# num_c_i - number of cows exposed t toxin level x_i
# y_i     - number of cows exposed to x_i and died
# theta_i - probability of death for cow exposed to x_i
# y_i | theta_i ~ Bin (n_i, theta_i)
# logit(theta_i) = alpha + beta * x_i
#
log_x <- c(-0.86, -0.30, -0.05, 0.73)
num_c <- c(5, 5, 5, 5)
y     <- c(0, 1, 3, 5)

data_agg.df <- data.frame (log_x, num_c, y)

# make individualized, per cow, data
logxs <- c( -0.86, -0.86, -0.86, -0.86, -0.86,
             -0.30, -0.30, -0.30, -0.30, -0.30,
             -0.05, -0.05, -0.05, -0.05, -0.05,
              0.73,  0.73,  0.73,  0.73,  0.73)
deads <- c(0, 0, 0, 0, 0,
           1, 0, 0, 0, 0,
           1, 1, 1, 0, 0,
           1, 1, 1, 1, 1)

# Get point estimated for alpha, beta, to create grid for Bayesian posterior
data.df <- data.frame (deads, logxs)
fit.glm <- glm (formula = deads ~ logxs, data = data.df, family = binomial)

alpha <- fit.glm$coefficients[1]  # 0.847
beta  <- fit.glm$coefficients[2]  # 7.75

ggplot (data = data_agg.df) +
  geom_jitter(aes(x=log_x, y=y/num_c), width=0.01, height = 0.01) +
  geom_function(fun = function (x) invlogit (alpha + beta * x), color="red") +
  labs (title="Pr (Cow dying given toxin=log_x)", 
        x="log of toxin", y="Proportion of cows dead")

# Given the point estimate, generate grid for Bayesian posterior
# use a for alpha, b for beta
n_a <- 1000
n_b <- 1000

a_min <- -5
a_max <-  10
b_min <- -10
b_max <-  40

flat_prior <- array (1, c(n_a, n_b))
likelihood <- array (NA, c(n_a, n_b))

a_grid <- seq (a_min, a_max, length.out=n_a)
b_grid <- seq (b_min, b_max, length.out=n_b)

start_time <- Sys.time()
for (a in 1:n_a){
  for (b in 1:n_b){
    likelihood[a, b] <- prod(invlogit(a_grid[a] + b_grid[b] * log_x)^y *
                                    (1 - invlogit(a_grid[a] + b_grid[b] * log_x))^(num_c-y))
  }
}
end_time <- Sys.time()
(elapsed <- end_time - start_time)

flat_prior <- flat_prior / sum (flat_prior)
likelihood <- likelihood / sum (likelihood)
flat_posterior <- flat_prior * likelihood
flat_posterior <- flat_posterior / sum (flat_posterior)  
######
# prior Bi-variate Normal
mu_a <- 0
mu_b <- 10
sigma_a <- 2
sigma_b <- 10 
rho <- 0.5
Mu <- c(mu_a, mu_b)
Sigma <- matrix (c(sigma_a^2, sigma_a*sigma_b*rho, sigma_a*sigma_b*rho, sigma_b^2),
                 ncol=2)

bivar_prior <- array (NA, c(n_a, n_b))

start_time <- Sys.time()
for (a in 1:n_a){
  for (b in 1:n_b){
    bivar_prior[a, b] <- (mvtnorm::dmvnorm(c(a_grid[a], b_grid[b]), mean=Mu, sigma=Sigma))
  }
}
end_time <- Sys.time()
(elapsed <- end_time - start_time)

bivar_prior     <- bivar_prior / sum (bivar_prior)
bivar_posterior <- bivar_prior * likelihood
bivar_posterior <- bivar_posterior/sum(bivar_posterior)

# Sample from the posterior using density, generate contour graphs
#
# 
n_sims <- 1000
a_sims_flat  <- rep (NA, n_sims)
b_sims_flat  <- rep (NA, n_sims)
a_sims_bivar <- rep (NA, n_sims)
b_sims_bivar <- rep (NA, n_sims)
for (s in 1:n_sims){
# pick from *marginal* distribution of alpha (integrate joint distribution over beta): 
  indx_flat  <- sample (1:n_a, 1, T, prob = rowSums (flat_posterior))
  indx_bivar <- sample (1:n_a, 1, T, prob = rowSums (bivar_posterior))
  a_sims_flat[s]  <- a_grid[indx_flat]
  a_sims_bivar[s] <- a_grid[indx_bivar]
  # sample beta, given a
  b_sims_flat[s]  <- sample (b_grid, 1, T, prob = flat_posterior[indx_flat, ])
  b_sims_bivar[s] <- sample (b_grid, 1, T, prob = bivar_posterior[indx_bivar, ])
}

distr_flat.df  <- data.frame (a_sims_flat, b_sims_flat)
flat_post <- ggplot (data = distr_flat.df) +
  geom_point(aes(x=a_sims_flat, y=b_sims_flat), size=0.1) +
  geom_point(aes(x=alpha, y=beta), color="red") +
  geom_density_2d_filled (aes (x=a_sims_flat, y=b_sims_flat), alpha=0.5, size=0.5, color="black") +
  labs (x="alpha", y="beta", title="Posterior alpha, beta - flat prior") +
  xlim (-5, 10) + ylim (-10, 40)

distr_bivar.df <- data.frame (a_sims_bivar, b_sims_bivar)
bivar_post <- ggplot (data = distr_bivar.df) +
  geom_point(aes(x=a_sims_bivar, y=b_sims_bivar), size=0.1) +
  geom_point(aes(x=alpha, y=beta), color="red") +
  geom_density_2d_filled (aes (x=a_sims_bivar, y=b_sims_bivar), alpha=0.5, size=0.5, color="black") +
  labs (x="alpha", y="beta", title="Posterior alpha, beta - bivariate normal prior") +
  xlim (-5, 10) + ylim (-10, 40)

grid.arrange(flat_post, bivar_post)

#
# Contour graphs
# Countour graphs
par(mfrow=c(2,1))
contour (x=a_grid, y=b_grid, z=flat_posterior, lwd=2, col="red", 
         vfont = c("sans serif", "bold"),
         xlim = c(-2, 4), ylim =c(-0, 30))
contour (x=a_grid, y=b_grid, z=bivar_posterior, lwd=2, col="red", 
         vfont = c("sans serif", "bold"),
         xlim = c(-2, 4), ylim =c(-0, 30))
par(mfrow=c(1,1))

persp (a_grid, b_grid, flat_posterior, theta=30, ticktype = "detailed")

xx <- a_grid[rep (1:n_a, rep (n_b,n_a))]
yy <- b_grid[rep((1:n_b), n_a)]

p1 <- ggplot (data.frame (alpha=xx, beta=yy, zz=as.vector (t(likelihood))), 
              aes(x=alpha, y=beta, z=zz)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Likelihood")
p2 <- ggplot (data.frame (alpha=xx, beta=yy, zz=as.vector ((flat_posterior))), 
              aes(x=alpha, y=beta, z=zz)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="FlatPrior")
p3 <- ggplot (data.frame (alpha=xx, beta=yy, zz=as.vector ((bivar_posterior))), 
              aes(x=alpha, y=beta, z=zz)) +
  xlim (-5, 10) + ylim (0, 25) +
  geom_contour_filled (alpha=0.6) +
  labs (title ="Bivariate Normal Prior")

grid.arrange(p1, p2, p3)
