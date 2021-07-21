#
# Regression Anrew Gelman
# http://www.stat.columbia.edu/~gelman/regression/
#
# When signal is low and noise is high, statistically significant patterns in data
# are likely to be wrong... the results are unlikely to replicate

# Cool Trix:
#
# df <- melt(data.frame(fit_1=br2_1,fit_3=br2_3))  
#
# sims  <- as.matrix(fit_bayes, pars = c("(Intercept)", "x"))
#
# attach("/Users/nevinaltaras/Downloads/nes.rda")
#
# the package sjplot has a nice function, plot_model()
#
# grid.arrange (p0, p1, ncol=2)  
#
# summary(fit$fitted, digits=2)
#
#  scale_color_manual(name = "Pr(y=1)", labels = c("0", "1"), values = c("red", "blue")) 
#
# 2   sd covers 95% of normal distr:  pnorm (2, lower.tail = T) - pnorm (2, lower.tail = F)
# 1   sd covers 68% of normal distr
# 2/3 sd covers 50% of normal distr
#
# refresh=0 results in no output in stan_glm
#
# +/- is \U00B1
#
#scale_color_identity(name = "Label",
#                     breaks = c("red", "blue","black"),
#                     labels = c("Prior", "Data", "Bayes"),
#                     guide = "legend") 
#
# ppc_dens_overlay (earnings$log_earn, yrep_1[subset,])

library(devtools)
library(ggplot2)
library(ggrepel)
library("bayesplot")
library("loo")
library("rstanarm")
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

invlogit <- plogis
theme_set(theme_bw())
theme_set(theme_cowplot())
options(digits=4)
set.seed (1234)

############################
# Restart Sep/2020
############################

# Ex 3-1
n1 <- 200
n2 <- 250
n3 <- 300
n4 <- 250

N <- n1+n2+n3+n4

w1 <- n1/N
w2 <- n2/N
w3 <- n3/N
w4 <- n4/N

0.5*w1 + 0.6*w2 +0.4*w2 + 0.3*w4 # = 0.425

#Ex 3-3
a <- rnorm (1000, 0, 1) 
b <- rnorm (1000, 5, 1) 
c <- rnorm (1000, -4, 3)
d <- rnorm (1000, -1, 7)

df <- data.frame(cbind (a, b, c, d))

pb <- ggplot(df) + geom_histogram(aes(x=b), binwidth=5/15) + xlim (-10, 10)
pc <- ggplot(df) + geom_histogram(aes(x=c), binwidth=1.5)  + xlim (-10, 10)
pd <- ggplot(df) + geom_histogram(aes(x=d), binwidth=1)    + xlim (-10, 10)

grid.arrange (pb, pc, pd, ncol=2)  

par (mfrow=c(2,2))
hist (rnorm (1000, 0, 1)) 
hist (rnorm (1000, 5, 1)) 
hist (rnorm (1000, -4, 3))
hist (rnorm (1000, -1, 7))
par (mfrow=c(1,1))

# Ex 3-5
par (mfrow=c(2,2))
hist (rbinom (10, 20, .3)) 
hist (rbinom (100, 20, .3)) 
hist (rbinom (1000, 20, .3))
hist (rbinom (10000, 10, .3))
par (mfrow=c(1,1))

#
# Ch 4
#
# The standard error is a measure of the variation in an estimate of the mean of
# the POPULATION as inferred from the SAMPLE, and gets smaller  
# as sample size gets larger, converging on zero as the sample increases in size. 

# When estimating the mean of an infinite population, given a simple random sample 
# of size n, the  standard error is σ/√ n, where σ is the standard deviation of 
# the measurements in the population. 

# SE for proportions: se = sqrt (estimate*(1-estimate)/n)

# SE of a difference: se = sqrt (se1^2 + se2^2)


# Ex 4-1
mu_treated = 0.50 
mu_control = 0.40 
se_treated = sqrt(0.5*0.5/500)  # 2.2%
se_control = sqrt(0.4*0.6/500)  # 2.2%

mean_diff = mu_treated - mu_control          # 0.1
se_diff <- sqrt (se_treated^2+se_control^2)  # 3.1

# Ex 4.2
#se <- sqrt(se_men^2 + se_women^2)
# assume se_men = se_women
# 0.05 <- sqrt (2*se_men^2)
se_men <- sqrt(0.05^2/2)  # 0.0353
# se_men <- 0.5/sqrt(n_men)
n_men <- (0.5/se_men)^2   # 200 
n_total <- 2*n_men        # 400

# Ex 4.3

N <- 100000
t <-0
for (i in 1:N){
  a <- sum (rbinom (20, 1, 0.3))
  b <- sum (rbinom (20, 1, 0.4))
  if (b > a) t <- t+1
}
cat ("\nN=", N, "prob =", t/N*100)   # 69.3%

# Ex 4.4
se <- 0.05
mu_1 <- 0.4 #, se_1 = sqrt(0.4*0.6/N_1)
mu_1 <- 0.3 #, se_2 = sqrt(0.3*0.7/N_2)
# (0.05) <- sqrt (se_1^2 + se_2^2)
# assume se_1 = se_2 = se, N_1=N_2=N
#se <- 0.5/sqrt(N)
#0.05 = sqrt(2 * 0.5^2/N) = sqrt (0.5/N)
N <- 0.5/0.05^2  # 200
N_total <- 400

# Alternate solution
N <- 10000
J <- 200
probs <- rep (NA, J)
for (j in 1:J){
  t <-0
  for (i in 1:N){
    a <- sum (rbinom (j, 1, 0.3))
    b <- sum (rbinom (j, 1, 0.4))
    if (b > a) t <- t+1
  }
  probs[j] <- t/N
  cat ("\nj=", j, "prob =", t/N*100)   # 69.3%
}

df <- data.frame (N=1:J, probs)
ggplot(df) + geom_point (aes(x=N, y=probs)) +
  ggtitle("Probability player A is better than player B") +
  xlab ("Number of shots taken by each player") +
  ylab ("Probability")
  
  
#Ex 5-1

ctr <- rep (NA, 1000)

for (c in 1:1000){
  ctr[c]=0
  prev_shot <- 1
  shot      <- 1
  while (!(shot==0 & prev_shot==0)){
    prev_shot <- shot
    shot <-  rbinom (1, 1, 0.6)
#    print (shot)
    ctr[c] <- ctr[c]+1
  }
}

# alternate
ctr  <- rep (NA, 1000)
pctg <- rep (NA, 1000)

for (c in 1:1000){
  shot <- rbinom (200, 1, 0.6)
  for (i in 1:200){
    if (shot[i]==0 & shot[i+1]==0) {
      ctr[c] <- i+1
      pctg[c] <- sum (shot[1:(ctr[c])])/ctr[c]
      break
    }
  }
}
  
hist (ctr, breaks=1000)
df <- data.frame (Shots=ctr, Percentage=pctg)
ggplot(df) + geom_point (aes(x=Shots, y=Percentage)) +
  ggtitle("Shooting percentage increases with Streak") +
  xlab ("Shots taken before two misses in a row") +
  ylab ("Shooting Percentage")

# Ex 5-2

J <- 1000
w <- rep (NA, J)
for (j in 1:J){
  w[j] <- 0
  sex <- sample(c("M", "F"), 10, prob=c(0.48,0.52), replace=T)
  for (i in 1:10){
    if (sex[i]=="M") w[j] <- w[j] + exp(rnorm (1, 5.13, 0.17))
    if (sex[i]=="F") w[j] <- w[j] + exp(rnorm (1, 4.96, 0.20))
  }
}
sum (w > 1750)/1000 # 5.5%

# Ex 5-3 a)
dbinom(3, size=10, prob=0.4)  # 21.5%
# Ex 5-3 b)
ctr <- 0
for (i in 1:10000){
  if (sum (rbinom (10, 1, 0.4)) == 3) ctr <- ctr+1
}
ctr/10000 # 21.2%

#################
# Chapter 5
#################
n_sims <- 1000
N <- 400 
p_g <- 0.488 # prob of girl

n_girls <- rep (NA, n_sims)
for (i in 1:n_sims){
  n_girls[i] <- rbinom (1, N, p_g)  # 192
}
# alternatively
n_girls <- rbinom (n_sims, N, p_g)  
hist (n_girls, breaks = 100)

# Account for twins
# FT = fraternal twin
# IT = identical twin
# SB = single birth

birth_type <- sample (c("FT", "IT", "NB"), replace=T, 400, prob=c(1/125, 1/300, 1-1/125-1/300 ))
girls <- rep (NA, 400)
for (i in 1:400){
  if(birth_type[i]=="NB") girls[i] <-   rbinom (1, 1, 0.488)
  if(birth_type[i]=="IT") girls[i] <- 2*rbinom (1, 1, 0.495)
  if(birth_type[i]=="FT") girls[i] <-   rbinom (1, 2, 0.495)
}
# alternatively
girls <- ifelse ((birth_type=="NB"),   rbinom (400, 1, 0.488), 
  ifelse((birth_type=="IT"),         2*rbinom (400, 1, 0.495),
                                       rbinom (400, 2, 0.495)))

n_girls <- sum(girls)

N <- 1000
girls <- rep (NA, N)
for (n in 1:1000){
  girls[n] <- sum(ifelse ((birth_type=="NB"),   rbinom (400, 1, 0.488), 
                  ifelse((birth_type=="IT"),  2*rbinom (400, 1, 0.495),
                                                rbinom (400, 2, 0.495))))
}
# alternatively
girls <- replicate (N, sum(ifelse ((birth_type=="NB"),   rbinom (400, 1, 0.488), 
                           ifelse((birth_type=="IT"),  2*rbinom (400, 1, 0.495),
                                                         rbinom (400, 2, 0.495)))))
hist (girls, breaks=100)

# Men's heights
# Men comprise 0.48% of population
# Men's height   ~ N(69.1, 2.9)
# Women's height ~ N(63.7, 2.7)

N_sim <- 1000
N_ppl <- 10
male <- rbinom (N_ppl, 1, 0.48)
avg_height <- mean (ifelse ((male == 1), rnorm (N_ppl, 69.1, 2.9), rnorm (N_ppl, 63.7, 2.7)))

avg_height_v <- replicate (N_sim, mean (ifelse ((male == 1), rnorm (N_ppl, 69.1, 2.9), rnorm (N_ppl, 63.7, 2.7))))
max_height_v <- replicate (N_sim, max (ifelse ((male == 1), rnorm (N_ppl, 69.1, 2.9), rnorm (N_ppl, 63.7, 2.7))))
hist (avg_height_v, breaks=100)
hist (max_height_v, breaks=100)

###################
# Bootstrapping
###################
earnings <- data.frame(read.csv('earnings.csv', header=T)) %>% filter (earn > 0)

earn <- earnings$earn
male <- earnings$male
ratio <- median (earn[male==0])/ median (earn[male==1]) # 60%

nrow <- nrow (earnings)
boot <- sample (nrow, replace=TRUE)
earn_boot <- earn[boot]
male_boot <- male[boot]
ratio_boot <- median (earn_boot[male_boot==0])/ median (earn_boot[male_boot==1]) # 64%

ratio_boot <- function (data){
  nrow <- nrow (data)
  boot <- sample (nrow, replace=TRUE)
  earn_boot <- data$earn[boot]
  male_boot <- data$male[boot]
  ratio_boot <- median (earn_boot[male_boot==0])/ median (earn_boot[male_boot==1]) 
  
  return (ratio_boot)
}

n_sims <- 10000
ratios_boot <- replicate (n_sims, ratio_boot(data=earnings))
mean (ratios_boot) # 60
sd (ratios_boot)   # 0.025

hist (ratios_boot, breaks=50)

#################
# Chapter 6:  Regression
##################

# make fake data
x <- 1:20
n <- length (x)
a <- 0.2
b <- 0.3
sigma <- 0.5
y <- a + b*x + sigma * rnorm (n)

df <- data.frame (x, y)
fit_1 <- stan_glm (y ~ x, data=df, refresh=0)
print (fit_1, digits=2)
b_hat <- fit_1$coefficients

df <- cbind (df, fitted=fit_1$fitted.values)

ggplot (df) + geom_point(aes(x=x, y=y)) +
#  geom_line(aes(x=x, y=fitted), color="red", size=2) +
  geom_abline(intercept=b_hat[1], slope=b_hat[2], color="blue", size=1) +
  ggtitle("Data and Fitted Regression Line") +
  geom_text(x= mean(x), y=b_hat[1] + b_hat[2]*mean(x) - 1, 
            label = paste ("y=", round(b_hat[1], 2), "+", round(b_hat[2], 2), "* x"))

# Use earnings_1990 data
earnings_1990 <- as_tibble (read.csv ('earnings_1990', header = T))

fit_2 <- stan_glm (earnk ~ height + male, data = earnings_1990, refresh=0)
print (fit_2, digits=2)

# Heights
heights <- as_tibble(read.table("Heights.txt", header=TRUE))

fit_3 <- stan_glm(daughter_height ~ mother_height, data = heights, refresh=0)
print (fit_3, digits=2)

b_hat <- fit_3$coefficients

df <- data.frame (daughter_height=heights$daughter_height, mother_height=heights$mother_height)
ggplot (df) + 
#  xlim (50, 70) +
  geom_jitter (aes(x=mother_height, y=daughter_height), width=0.5, size=0.1) +
  geom_abline(intercept=b_hat[1], slope=b_hat[2], color="blue", size=1) +
  ggtitle("Data and Fitted Regression Line") +
  geom_text(x= mean(heights$mother_height) -5, y=b_hat[1] + b_hat[2]*mean(heights$mother_height) -10, 
            label = paste ("y=", round(b_hat[1], 2), "+", round(b_hat[2], 2), "* x"), size=5, color="blue")

# Fake data on final grade vs midterm grade
n <- 1000
true_ability <- rnorm (n, 50, 10)
noise_1      <- rnorm (n,  0, 10)
noise_2      <- rnorm (n,  0, 10)
midterm <- true_ability + noise_1
final   <- true_ability + noise_2
exams <- data.frame(midterm, final)

fit_1 <- stan_glm(final ~ midterm, data=exams)
print (fit_1, digits=2)

b_hat    <- fit_1$coefficients 
b_hat_se <- fit_1$ses
# alternatively 
b_hat    <- coef(fit)
b_hat_se <- se(fit)

ggplot (exams) + 
  geom_point (aes(x=midterm, y=final)) +
  geom_abline(intercept=b_hat[1], slope=b_hat[2], color="blue", size=1) +
  ggtitle("Data and Fitted Regression Line") +
  geom_text(x= mean(exams$midterm) -5, y=b_hat[1] + b_hat[2]*mean(exams$midterm) -30, 
            label = paste ("y=", round(b_hat[1], 2), "+", round(b_hat[2], 2), "* x"), size=5, color="blue")


# Ex: 6-1
n <- 50
a <- 30
b <- 1.2
sigma <- 10

midterm <- runif (n, 0, 50)
final   <- a + b * midterm + rnorm (n, 0, sigma)
df <- data.frame(midterm, final)
fit <- stan_glm (final ~ midterm, data=df, refresh=0)
print (fit)

ggplot (df) + geom_point (aes(x=midterm, y=final)) +
  geom_abline(intercept = a, slope = b, color="red") +
  geom_abline(intercept = coef(fit)[1]+2*sigma(fit), slope = coef(fit)[2], color="grey",  linetype="dashed") +
  geom_abline(intercept = coef(fit)[1]+  sigma(fit), slope = coef(fit)[2], color="black", linetype="dashed") +
  geom_abline(intercept = coef(fit)[1],              slope = coef(fit)[2], color="blue") +
  geom_abline(intercept = coef(fit)[1]-  sigma(fit), slope = coef(fit)[2], color="black", linetype="dashed") +
  geom_abline(intercept = coef(fit)[1]-2*sigma(fit), slope = coef(fit)[2], color="grey",  linetype="dashed") +
  ggtitle ("Final grade projected from midterm") +
  geom_text (x=0, y=100, hjust=0, label ="True effect",     color="red") +
  geom_text (x=0, y= 95, hjust=0, label ="Regression line", color="blue") +
  geom_text (x=0, y= 90, hjust=0, label ="Fit \U00B1 1 sd", color="black") +
  geom_text (x=0, y= 85, hjust=0, label ="Fit \U00B1 2 sd", color="grey") 

  
# EX: 6-2
n <- 100
a <- 30
b <- 1.2
sigma <- 10

stat <- function (n, a, b, sigma){
  x <- runif (n, 0, 100)
  y   <- a + b * x + rnorm (n, 0, sigma)
  df <- data.frame(x, y)
  fit <- stan_glm (y ~ x, data=df, refresh=0)
  print (fit)
  b_hat <- fit$coefficients
  ggplot (df) + geom_point (aes(x=x, y=y)) +
  geom_abline(intercept=b_hat[1], slope=b_hat[2], color="blue", size=1) 
}
stat (n, a, b, sigma)

################
# Chapter 7
################
hibbs <- read.table ("hibbs.dat", header = T) %>% rename (vote=inc.party.vote)

M1 <- stan_glm (vote ~ growth, data=hibbs, refresh=0)
print (M1) # b_hat[1]=, 46.2, b_hat[2]=3.1, sigma=3.9
b_hat <- M1$coefficients
avg_growth  <- mean (hibbs$growth)
ggplot (hibbs) + geom_point (aes(x=growth, y=vote)) +
  geom_abline (intercept = b_hat[1], slope = b_hat[2]) +
  geom_text(x= mean(avg_growth) + 1.5, y=b_hat[1] + b_hat[2]*avg_growth + 2, 
            label = paste ("y=", round(b_hat[1], 1), "+", round(b_hat[2], 1), "* x"), size=4) +
  xlab ("Economic growth") + 
  ylab ("Incumbent's party vote share")

# growth in 2016 was 2%
vote_hillary <- 46.2 + 3.1 * 2 + 3.9 * rnorm (10000, 0, 1)
hist (vote_hillary, breaks=100)

prob_hillary_win <- mean (vote_hillary > 50) # 74%
# alternatively
prob_hillary_win <- 1-pnorm (50, 52.3, 3.9)  # 72%

##################
# Create regenerative model - fake data: assume a=46.3, b=3.0, sigma=3.9
##################
b     <- c(46.3, 3.0)
sigma <- 3.9
x <- hibbs$growth
n <- length (x)  # 16 - not very large
# simulate fake votes, y
y <- b[1]+b[2]*x+rnorm (n, 0, sigma)
df <- data.frame (x, y)
fit <- stan_glm (y ~ x, data=df, refresh=0)
print (fit) # a_hat=46.3 +/- 1.5, , b_hat=2. 7+/- 0.6, sigma=3.4 +/- 0.7
# alternatively
b_hat    <- fit$coefficients
b_hat_se <- fit$ses
cat("b=", b[2], ", b_hat=", round(b_hat[2],1), ", se_b=", round(b_hat_se[2], 1))

n_sims <- 1000

covers_68 <- rep (NA, n_sims)
covers_95 <- rep (NA, n_sims)

for (s in 1:n_sims){
  y <- b[1]+b[2]*x+rnorm (n, 0, sigma)
  df <- data.frame (x, y)
  fit <- stan_glm (y ~ x, data=df, refresh=0)
  
  b_hat    <- fit$coefficients
  b_hat_se <- fit$ses

  covers_68[s] <- abs (b[2]-b_hat[2])/b_hat_se[2] 
  covers_95[s] <- abs (b[2]-b_hat[2])/b_hat_se[2] 
  
 cat (s, covers_68[s] < 1, covers_95[s] < 2, "\n")
}

mean (covers_68 <1) # 68%
mean (covers_95 <2) # 95%

# n=16 so normal distribution may not be best
# try t distribution

t_68 <- qt (0.84,  n-2) # 1.03
t_95 <- qt (0.975, n-2) # 2.14

mean (covers_68/t_68 < 1) # 69%  
mean (covers_95/t_95 < 1) # 97%

# Indicator Variables
n_0 <- 20
y_0 <- rnorm (n_0, 2, 5)
fake_0 <- data.frame(y_0)
fit_0 <- stan_glm (y_0 ~ 1, data = fake_0, refresh=0)
print (fit_0)

n_1 <- 30
y_1 <- rnorm (n_1, 8, 5)

diff <- mean (y_1) - mean (y_0)   # 4.1
se_0 <- sd (y_0)/sqrt (n_0)       # 1.1
se_1 <- sd (y_1)/sqrt (n_1)       # 0.9
se   <- sqrt(se_1^2 + se_1^2)     # 1.3

# let x be indicator variable
x <- c(rep (0, n_0), rep (1, n_1))
y <- c(y_0, y_1)
df <- data.frame (x, y)
fit <- stan_glm (y ~ x, data =df, refresh=0 )
print (fit)
print (fit$coefficients[2])       # 4.2

# Ex 7-2
n <- 100
a <- 5
b <- 7
sigma <- 3

x <- runif (n, 0, 50)
y   <- a + b*x + rnorm (n, 0, sigma)
df <- data.frame(x, y)
fit <- stan_glm (y ~ x, data=df, refresh=0)
print (fit)

ggplot (df) + geom_point (aes(x=x, y=y)) +
  geom_abline(intercept = a, slope = b, color="red") +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color="blue") +
#  ggtitle ("Final grade projected from midterm") +
  geom_text (x=0, y=340, hjust=0, label = paste("Tue line=", a, "+ ", b, "* x"), color="red") +
  geom_text (x=0, y=320, hjust=0, label =
               paste("Fitted line=", round (coef(fit)[1], 1), "+ ", round (coef(fit)[2], 1), "* x"),
                                               color="blue") 

# Ex 7-6
hibbs$gr_indic <- ifelse (hibbs$growth < 2, 0, 1 )
# a)
mean (hibbs$vote[hibbs$gr_indic==1]) - mean (hibbs$vote[hibbs$gr_indic==0])       # b_hat=5.5%
se_indic0 <- sd (hibbs$vote[hibbs$gr_indic==0])/sqrt(sum (hibbs$gr_indic==0))
se_indic1 <- sd (hibbs$vote[hibbs$gr_indic==1])/sqrt(sum (hibbs$gr_indic==1))
sqrt (se_indic0^2 + se_indic1^2)                                                  # =2.5

# b)
fit <- stan_glm (vote ~ gr_indic, data = hibbs, refresh=0)
print (fit)  # b_hat==5.5 +/- b_hat_se=2.6


##############
# Chapter 8
##############
x <- 1:10
y <- c(1, 1, 2, 3, 5, 8, 13, 21, 34, 55)
df <- data.frame(x, y)
fit <- stan_glm (y ~ x, data = df, refresh=0)
print (fit)
sims <- as.matrix (fit)  # x=5.1 ,x_se=1.1 -> 2 sd CI = [2.9, 7.3]

quantile (sims[,2], c(0.025, 0.975))  #                 [2.6, 7.4]

# back to Hibbs data
hibbs <- read.table ("hibbs.dat", header = T) %>% rename (vote=inc.party.vote)

M1 <- stan_glm (vote ~ growth, data=hibbs, refresh=0)
print (M1) 
sims <- as.matrix(M1)
b_hat <- M1$coefficients

avg_growth  <- mean (hibbs$growth)

n_sim <- 50 # sample n_sim simulations from sims matrix
s_sample <- sample(nrow (sims), n_sim, replace = F)
ggplot (hibbs) + 
  geom_abline (intercept = sims[s_sample,1], slope = sims[s_sample,2], color="gray") +
  geom_abline (intercept = b_hat[1], slope = b_hat[2], color="red", size=2) +
  geom_text(x= mean(avg_growth) + 1.5, y=b_hat[1] + b_hat[2]*avg_growth + 2, 
            label = paste ("y=", round(b_hat[1], 1), "+", round(b_hat[2], 1), "*x"), size=4) +
  geom_point (aes(x=growth, y=vote)) +
  xlab ("Economic growth") + 
  ylab ("Incumbent's party vote share")

##################
# Chapter-9
##################
# back to Hibbs data
hibbs <- read.table ("hibbs.dat", header = T) %>% rename (vote=inc.party.vote)

M1 <- stan_glm (vote ~ growth, data=hibbs, refresh=0)
print (M1) 

b_hat <- M1$coefficients # 46, 3.0
b_hat <- coef(M1)  # Alternative
b_hat_se <- M1$ses       # 1.8, 0.7
b_hat_se <- se(M1) # Alternative
sigma <- sigma (M1)      # 3.9

sims <- as.matrix(M1)
n_sims <- nrow (sims)

Median <- apply (sims, 2, median) # 46.3, 3.0, 3.9
MAD_SD <- apply (sims, 2, mad)    # 1.8, 0.7, 0.8

a <- sims[,1]
b <- sims[,2]
df <- data.frame(a, b)

ggplot(df) + geom_histogram(aes(x=a), binwidth=1, breaks=seq(round(min(a)), round(max(a)), by=1), color="gray") +
  geom_vline (xintercept = Median[1], size=1) +
  geom_segment (x=Median[1]-MAD_SD[1], xend=Median[1]+MAD_SD[1], y=600, yend=600, 
                arrow = arrow()) +
  geom_segment (xend=Median[1]-MAD_SD[1], x=Median[1]+MAD_SD[1], y=600, yend=600, 
                arrow = arrow()) +
  geom_segment (x=Median[1]-2*MAD_SD[1], xend=Median[1]+2*MAD_SD[1], y=250, yend=250, 
                arrow = arrow()) +
  geom_segment (xend=Median[1]-2*MAD_SD[1], x=Median[1]+2*MAD_SD[1], y=250, yend=250, 
                arrow = arrow()) +
  ggtitle ("Posterior simulations of the intercept, a, \nand posterior median \U00B1 1 and 2 std err")

ggplot(df) + geom_histogram(aes(x=b), breaks=seq(round(min(b)), round(max(b)), by=0.5), color="gray") +
  geom_vline (xintercept = Median[2], size=1) +
  geom_segment (x=Median[2]-MAD_SD[2], xend=Median[2]+MAD_SD[2], y=600, yend=600, 
                arrow = arrow()) +
  geom_segment (xend=Median[2]-MAD_SD[2], x=Median[2]+MAD_SD[2], y=600, yend=600, 
                arrow = arrow()) +
  geom_segment (x=Median[2]-2*MAD_SD[2], xend=Median[2]+2*MAD_SD[2], y=250, yend=250, 
                arrow = arrow()) +
  geom_segment (xend=Median[2]-2*MAD_SD[2], x=Median[2]+2*MAD_SD[2], y=250, yend=250, 
                arrow = arrow()) +
  ggtitle ("Posterior simulations of the slope, b, \nand posterior median \U00B1 1 and 2 std err")

ggplot (df) + 
  geom_point(aes(x=a, y=b), size=.1) +
  ggtitle ("Posterior draws on regression coefficients, a, b")

s_sample <- sample (n_sims, 100)
ggplot (hibbs) +
  geom_abline (intercept = sims[s_sample,1], slope = sims[s_sample,2], size=0.1) +
  geom_point (aes(x=growth, y=vote)) +
  ggtitle("Data and 100 posterior draws of th eline a + b*x") +
  xlab ("Growth rate") + 
  ylab ("Incumbent party's vote share")

#################################
# Prediction and Uncertainty: 
#################################
# y_new, given new data x_new
# let x: medicine dose, 
#     y: blood pressure, 
#
# 1) Point prediction: Scalar,
#          Best estimate of average blood pressure in population conditional on dose=x_new 
#          y_new <- a_hat + b_hat*x_new  - no uncertainty
# 2) Linear Prediction w/ Uncertainty:  Vector
#          Modeled Average blood pressure conditional on dose=x_new in the population.  
#          Uncertainty due to estimation of a and b 
#          y_new <- a     + b    *x_new  (a, b are vectors from simulations matrix )
# 3) Predictive Distribution: Vector of distributions
#          Blood pressure of a single individual drawn at random from population
#          y_new <- a     + b    *x_new  + rnorm (0, sigma)
#  As sample size approached infinity, a and b are estimated more precisely, 
#  uncertainty in Linear Prediction approaches zero
#  uncertainty in Predictive Distribution approaches sigma
####################################
posterior_graph_2 <- function (x, y, x_new, x_lab, y_lab){
  df <- data.frame (x, y)
  fit <- stan_glm(y ~ x, data = df, refresh=0)
#  print (fit) 
  new <- data.frame (x=x_new)

  point_predict <- predict            (fit, newdata = new)
  post_linpred  <- posterior_linpred  (fit, newdata = new)
  post_predict  <- posterior_predict  (fit, newdata = new)
  
  se_post_linpred <- mad (post_linpred)  
  se_post_predict <- mad (post_predict) 
  
  df_predict      <- data.frame (x1=x_new[1], y1=point_predict[1],
                                 x2=x_new[2], y2=point_predict[2])
  df_post_linpred <- data.frame (x1=rep(x_new[1], 4000), y1=post_linpred[,1],
                                 x2=rep(x_new[2], 4000), y2=post_linpred[,2])
  df_post_predict <- data.frame (x1=rep(x_new[1], 4000), y1=post_predict[,1],
                                 x2=rep(x_new[2], 4000), y2=post_predict[,2])
  
  p <-  ggplot () + 
    geom_point(data=df_post_predict, aes(x=x1, y=y1, color="yellow")) +
    geom_point(data=df_post_predict, aes(x=x2, y=y2, color="yellow")) +
    geom_point(data=df_post_linpred, aes(x=x1, y=y1, color="red"),  size=0.6) +
    geom_point(data=df_post_linpred, aes(x=x2, y=y2, color="red"),  size=0.6) +
    geom_point(data=df_predict,      aes(x=x1, y=y1, color="blue"), size=3) +
    geom_point(data=df_predict,      aes(x=x2, y=y2, color="blue"), size=3) +
    geom_point(data=df, aes(x=x, y=y), size=0.5) +
    geom_abline(intercept = coef(fit)[1], slope=coef(fit)[2]) +
    xlab (x_lab) + ylab (y_lab) +
    ggtitle (paste(y_lab,  "vs", x_lab, " - data, best linear fit", 
                   "\nPredictions for", x_lab, "= {", toString(x_new), "}")) +
    scale_color_identity(name = "Predictions",
                         breaks = c("blue", "red","yellow"),
                         labels = c("Point", "LinPred", "PostPred"),
                         guide = "legend") 
  return (p)
}

# Earnings Ex
earnings_1 <- as_tibble (read.csv ("Earnings", header=T)) %>%
  filter (is.na (height) != T & is.na (weight)  !=T)

earnings_1$c_height <- earnings_1$height - mean (earnings_1$height)
fit_2 <- stan_glm (weight ~ c_height, data = earnings_1, refresh=0)
print (fit_2)
sims <- as.matrix(fit_2)

sigma <- sigma (fit_2) # 29.0
x_new <- c(4, 4)

y_new_point <- coef(fit_2)[1] + 4*coef(fit_2)[2]  # 176 lb
y_new_line  <- sims[,1] + 4*sims[,2]                                 # linear predictor
y_new_pred  <- sims[,1] + 4*sims[,2] + rnorm (nrow(sims), 0, sigma)  # predicted value

n <- nrow (earnings_1)  # 1816
sigma_hat_linpred    <- sigma (fit_2) * 
  sqrt (1/n + x_new/sum (earnings_1$c_height^2) )     #  0.8 (x_bar =0)
sigma_hat_prediction <- sigma (fit_2) * 
  sqrt (1 + 1/n + x_new/sum (earnings_1$c_height^2) ) # 29.0 (x_bar =0)

p <- posterior_graph_2 (x=earnings_1$c_height, y=earnings_1$weight, x_new=x_new, 
                        x_lab="Height", y_lab="Weight")
  

###############
# Election Ex
new <- data.frame (growth=2.0)
# 1)
y_point_predict <- predict (M1, newdata=new)     # 52.4
# alternatively
a_hat <- coef(M1)[1]
b_hat <- coef(M1)[2]
y_point_predict <- a_hat + 2.0*b_hat # 52.4, Hilary's predicted vote %age

#2)
y_linpred <- posterior_linpred(M1, newdata = new)

# alternatively
sims <- as.matrix (M1)
a <- sims[, 1]
b <- sims[, 2]
y_linpred <- a + 2.0*b

#3)
y_pred <- posterior_predict (M1, newdata = new)

# alternatively
sigma <- sims[, 3]
y_pred <- a + 2.0*b + rnorm (n_sims, 0, sigma)

# Calculate statistics
y_pred_median <- median (y_pred) # 52.3
y_pred_mad    <- mad (y_pred)    #  3.9
win_prob <- mean (as.numeric(y_pred > 50))  # 72.4  Hilary's win probability
##############
# Histograms of y_new
###############
y_pp_inf_data <- rnorm (4000, y_point_predict, median(sigma))

df_point_pred  <- data.frame (y_new=y_point_predict)
df_linpred     <- data.frame (y_new=y_linpred)
df_post_pred   <- data.frame (y_new=y_pred) %>% rename (y_new = X1)
df_pp_inf_data <- data.frame (y_new=y_pp_inf_data)

p1 <- ggplot(df_post_pred)   + geom_histogram (aes(x=y_new), breaks=seq(30, 75, by=1)) +
  ggtitle ("Posterior preditcion")
p2 <- ggplot(df_linpred)     + geom_histogram (aes(x=y_new), breaks=seq(30, 75, by=1))  +
    ggtitle ("Linear preditcion")
p3 <- ggplot(df_point_pred)  + geom_histogram (aes(x=y_new), breaks=seq(30, 75, by=1)) +
  ggtitle ("Point preditcion")
p4 <- ggplot(df_pp_inf_data) + geom_histogram (aes(x=y_new), breaks=seq(30, 75, by=1)) +
  ggtitle ("Post. pred. with infinite data")

grid.arrange (p3, p2, p1, p4, ncol=2)  

###################
# Begin make up example
###################
n <- 100
x <- runif(n, 0, 10)
a <- 50
b <- 2
sigma <- 5
y <- a + b*x + rnorm (n, 0, sigma)
df <- data.frame (x, y)
fit <- stan_glm (y ~ x, data = df, refresh =0)
print (fit)

x_new <- c(0, 5)

posterior_graph_2 (x=x, y=y, x_new=x_new, x_lab="x", y_lab="y")

new <- data.frame (x=x_new)

y_point_pred <- predict           (fit, newdata = new)
y_line_pred  <- posterior_linpred (fit, newdata = new)
y_post_pred  <- posterior_predict (fit, newdata = new)

mad (y_line_pred[,1]) # 1.0
mad (y_line_pred[,2]) # 0.8
mad (y_post_pred[,1]) # 5.0
mad (y_post_pred[,2]) # 5.0

df_y_point_pred <- data.frame(x1=x_new[1], y1=y_point_pred[1], 
                             x2=x_new[2], y2=y_point_pred[2])
df_y_line_pred <- data.frame(x1=rep(x_new[1], n_sims), y1=y_line_pred[,1], 
                             x2=rep(x_new[2], n_sims), y2=y_line_pred[,2])
df_y_post_pred <- data.frame(x1=rep(x_new[1], n_sims), y1=y_post_pred[,1], 
                             x2=rep(x_new[2], n_sims), y2=y_post_pred[,2])

# Differences
point_diff    <- rep(y_point_pred[2] - y_point_pred[1], 1000)  # 17.5
linepred_diff <- y_line_pred[,2] - y_line_pred[,1]  
postpred_diff <- y_post_pred[,2] - y_post_pred[,1]  

ggplot () +
  geom_histogram(aes(x=postpred_diff, color="yellow"), binwidth = 0.5) +
  geom_histogram(aes(x=linepred_diff, color="blue"),   binwidth = 0.5) +
  geom_histogram(aes(x=point_diff,    color="red"),    binwidth = 0.1) +
  ggtitle ("Predictions for difference between x=5 and x=0")
  
# Vary number of data
pred_plots <- function (n){
  x <- runif(n, 0, 10)
  a <- 50
  b <- 2
  sigma <- 5
  y <- a + b*x + rnorm (n, 0, sigma)
  df <- data.frame (x, y)
  fit <- stan_glm (y ~ x, data = df, refresh =0)
  
  x_new <- c(0, 5)
  new <- data.frame (x=x_new)

  y_line_pred  <- posterior_linpred (fit, newdata = new)
  y_post_pred  <- posterior_predict (fit, newdata = new)
  
  linepred_diff <- y_line_pred[,2] - y_line_pred[,1]  
  postpred_diff <- y_post_pred[,2] - y_post_pred[,1]  
  
  p_line_pred <- ggplot () + geom_histogram(aes(x=linepred_diff), color="red",    
                                            breaks=seq( 0, 25, .5)) +
    ggtitle (paste (n, "data points"))
  p_post_pred <- ggplot () + geom_histogram(aes(x=postpred_diff), color="yellow", 
                                            breaks=seq(-20, 40, 1.2)) +
        ggtitle (paste (n, "data points"))
  
  plots <- list(p1=p_line_pred, p2=p_post_pred)
  
  return (plots)  
}

plots <- pred_plots (n=50)
p1 <- plots$p1
p2 <- plots$p2

plots <- pred_plots (n=200)
p3 <- plots$p1
p4 <- plots$p2

plots <- pred_plots (n=10000)
p5 <- plots$p1
p6 <- plots$p2

grid.arrange (p1, p2, p3, p4, p5, p6, ncol=2)  
###################
# End make up example
###################

# Prediction given a range of input values
# Assume we want to predict Hilary's chances for grid of 
# growth rates from -2 to +4
new_grid <- data.frame (growth=seq(from=-2, to=4, by=0.5))

# 1) Point prediction
y_point_predict_grid <-     predict (M1, newdata=new_grid) # 13 scalars (one for each growth rate)
# 2) Linear prediction
y_linpred_grid <- posterior_linpred (M1, newdata=new_grid) # 13 X n_sims
# 3) Posterior prediction
y_pred_grid <-    posterior_predict (M1, newdata=new_grid) # 13 X n_sims

# Propagating uncertainty:  What if  gdp measurement was uncertain?

# assume gdp <- N (2, 0.3)
n_sims <- 4000
x_new <- rnorm (n_sims, 2, 0.3)
y_pred <- rnorm (n_sims, a + b*x_new, sigma)
hist (y_pred)
median(y_pred) # 52.3
mad(y_pred)    #  4.1
mean (y_pred > 50) # 71
# Now compare to when no uncertainty in gdp estimation
x_new <- rnorm (n_sims, 2, 0)
y_pred <- rnorm (n_sims, a + b*x_new, sigma)
hist (y_pred)
median(y_pred) # 52.3
mad(y_pred)    #  3.9
mean (y_pred > 50) # 73.4
# with less gdp uncertainty, Hillary's vote has lower se, and her changes of winning is higher.

# Earnings
earnings_1 <- as_tibble (read.csv ("Earnings", header=T))

fit_1 <- stan_glm (weight ~ height, data = earnings_1, refresh=0)
print (fit_1)  # a=-173 +/- 11.8, b=5.0 +/- 0.2, sigma=29
# a=-173 implies, height=0in implies weight=-173lb
# recast with centering height
earnings_1$c_height <- earnings_1$height - mean (earnings_1$height)
fit_2 <- stan_glm (weight ~ c_height, data = earnings_1, refresh=0)
print (fit_2)  # a=156 +/- 0.7, b=5.0 +/- 0.2, sigma=29
# People of average height 66.5in are expected to weigh 156 lbs
# people of 67.5 in are expected to weigh 156 + 5=161 lbs

# Let's do predictions for a person 4 in taller than average, i.e. 70.5in tall
new <- data.frame (c_height=4)
y_point_pred <- predict (fit_2, newdata = new)  # 176

y_linpred   <- posterior_linpred  (fit_2, newdata = new)  
median (y_linpred) # 176
mad (y_linpred)    #   1

y_postpred   <- posterior_predict  (fit_2, newdata = new)  
median (y_postpred) # 176
mad (y_postpred)    #  29

####################
# BAYESIAN SYNTHESIS
####################
# theta - parameter we want to estimate
# theta_hat_prior - prior estimate,    with se=se_prior
# theta_hat_data  - data estimate,     with se=se_data
# theta_hat_bayes - bayesian estimate, with se=se_bayes

# FOR NORMALLY DISTRIBUTED theta
# define: precision = 1/se^2 = 1/var

# 1) 
# theta_hat_bayes = precision weighted average {theta_hat_prior, theta_hat_data}
#                 = (prec_prior*theta_hat_prior + prec_data*theta_hat_data)/(prec_prior + prec_data)

# 2)
# prec_bayes = prec_prior + prec_data
#
# Combining prior and data information increases precision!

# Example: Voting in 2016: 
theta_prior <- 0.524 
se_prior    <- 0.041
# data: n=400, 190 will vote for Hilary
n <- 400
theta_data <- 190/400                # 0.475
se_data <- sqrt (0.475*(1-0.475)/n)  # 0.025

prec_prior <- 1/se_prior^2            #  594
prec_data  <- 1/se_data^2             # 1604
prec_bayes <- prec_prior + prec_data  # 2199

se_bayes    <- sqrt (1/prec_bayes)       # 0.21 - lower than se_prior and se_data!
theta_bayes <- (theta_prior*prec_prior + theta_data*prec_data)/(prec_prior + prec_data) 
# 48.8% - between 47.5% and 52.4%, closer to data, 47.5, 
# which has higher precision than prior - 1604 vs 594

# graphs
ggplot () +
  geom_function(fun = dnorm, args = list(mean = theta_prior, sd = se_prior), aes(color="red")) +
  geom_function(fun = dnorm, args = list(mean = theta_data,  sd = se_data),  aes(color='blue')) +
  geom_function(fun = dnorm, args = list(mean = theta_bayes, sd = se_bayes), aes(color="black")) +
  scale_color_identity(name = "Distr's",
                       breaks = c("red", "blue","black"),
                       labels = c("Prior", "Data", "Bayes"),
                       guide = "legend") +
  ggtitle ("theta_hat: Prior, Data, and Bayesian Estimates") +
  xlim (0.35, 0.7) 

# now if 
se_data <- 0.075 # (instead) of 0.25

prec_data  <- 1/se_data^2             #  178
prec_bayes <- prec_prior + prec_data  #  772

se_bayes    <- sqrt (1/prec_bayes)       # 0.36 - still lower than se_prior and se_data!
theta_bayes <- (theta_prior*prec_prior + thata_data*prec_data)/(prec_prior + prec_data) 
# 51.3% - between 47.5% and 52.4%, closer to prior, 47.5, 
# which has higher precision than data - 594 vs 178
#################
# Bayesian Priors
#################

# Hibbs: uniform prior
M3 <- stan_glm (vote ~ growth, data=hibbs, refresh=0,
                prior = NULL, prior_intercept = NULL, prior_aux = NULL)
print (M3)

# Hibbs: default prior
# prior for intercept   a   is N (mean_y, 2.5*sd(y))
# prior for coefficient b_k is N (0, 2.5*sd(y)/sd(x))
# prior for sigma           is exp (1/sd(y))
sd_x <- sd (hibbs$growth)
sd_y <- sd (hibbs$vote)
mean_y <- mean (hibbs$vote)

M1 <- stan_glm (vote ~ growth, data = hibbs, refresh=0)
print (M1) # a=46.3 +/- 1.7, b=3.0 +/- 0.7, sigma=3.9 +/- 0.8
M1a <- stan_glm (vote ~ growth, data = hibbs, refresh=0,
                 prior           = normal (0, 2.5*sd(y)/sd(x)),
                 prior_intercept = normal (mean_y, 2.5*sd(y)), 
                 prior_aux       = exponential(1/sd_y))

# Weakly informative prior
# prior intercept is a + x_bar*b - what can we say about vote percentage on average?
# 50 +/- 10  => N (50, 10)
# prior slope - should be positive, 1% growth should give rise to less than 5 change in vote
# 5, but not muck more than 10.  => 5 +/- 5    => N (5, 5)

M4 <- stan_glm (vote ~ growth, data = hibbs, refresh=0,
                prior           = normal ( 5,  5),
                prior_intercept = normal (50, 10))
print (M4) # a=46.1 +/- 1.7, b=3.1 +/- 0.7, sigma=3.9 +/- 0.8   SAME AS M1!
# Weak informative prior makes no difference!!!

#################
# Strong prior
# Example: Sex and beauty: Difference between percentage of girls to beautiful parents 
#          vs ugly parents
n <- 3000
n_beauty <- n/10
n_ugly   <- n - n_beauty

theta_beauty <- 0.56  # proportion of girls among beautiful parents
se_beauty <- sqrt (theta_beauty*(1-theta_beauty)/n_beauty) # 0.03
theta_ugly   <- 0.46  # proportion of girls among ugly parents
theta_data <- 0.08
se_ugly <- sqrt (theta_ugly*(1-theta_ugly)/(n-n/10))       # 0.01
se_data <- sqrt (se_beauty^2 + se_ugly^2)                  # 0.03      
# 2 std CI [2, 14], doesnot contain 0

# Prior to seeing the data we assume there is no difference between 
# percentage of girls born to beautiful or ugly parents.  
# Assume mean of difference
theta_prior = 0
# Given that the proportion of girls is pretty constant in all populations
# assume small se_prior
# Assume 
se_prior <- 0.0025 # quarter of 1%
# theta_prior=N(0, 0.0025)

precision_prior <- 1/se_prior^2                      # 160000
precision_data  <- 1/se_data^2                       #   1095
precision_bayes <- precision_prior + precision_data  # 161095

se_bayes <- sqrt (1/precision_bayes) # 0.25% !!
theta_bayes <- (theta_prior*precision_prior + theta_data*precision_data)/
  (precision_prior+precision_data)   # 0.0005, 0.05%!!!

# The Bayesian estimate of difference in sexes is 0.05% with se 0.25% => 2sd CI [-0.2, 0.3]
# Includes zero: not statistically significant!
x <- seq(-2,2,1)
y <- c(50, 44, 50, 47, 56)
sexratio <- data.frame(x, y)

# first no (weak prior)
fit_weak_prior <- stan_glm(y ~ x, data = sexratio, refresh = 0) 
print (fit_weak_prior) # a=49.4    2.0, b=1.4    1.4, sigma=4.6    1.7

sims <- as.matrix(fit_weak_prior) 
s_sample <- sample (4000, 100)
p1 <- ggplot () +
  geom_abline(intercept = sims[s_sample, 1], slope = sims[s_sample, 2], color="gray") +
  geom_abline(intercept = coef(fit)[1],      slope = coef(fit)[2]) +
  geom_point (aes(x=x, y=y)) +
  xlab ("Attractiveness of parent") +
  ylab ("Percentage of girl babies") +
  ggtitle ("Regression line and \nposterior stimulations \ngiven weak prior")

df <- data.frame (a = sims[,1], b = sims[,2])
p2 <- ggplot (df) +
  geom_point (aes (x=a, y=b), size=0.01) +
  ggtitle("Posterior simulations \nunder default prior") +
  xlim (30, 70) + ylim (-15, 15)
# now stronger prior
# Prior intercept: we know that %girl births are very stable around 48.8 =>
#  ~ N (48.8, 0.5)
# Slope intercept: Before we do the study we don't know which direction parents' 
# attractiveness will influence girl births.  If anything we'll assume small amount
#  ~ N (0, 0.2)

fit_strong_prior <- stan_glm(y ~ x, data = sexratio, refresh = 0,
                             prior = normal (0, 0.2),
                             prior_intercept = normal (48.8, 0.5))

print (fit_strong_prior) #48.8 +/- 0.5, b=0.0 +/- 0.2, sigma=43 +/- 1.3

sims <- as.matrix(fit_strong_prior) 
s_sample <- sample (4000, 100)
p3 <- ggplot () +
  geom_abline(intercept = sims[s_sample, 1], slope = sims[s_sample, 2], color="gray") +
  geom_abline(intercept = coef(fit)[1],      slope = coef(fit)[2]) +
  geom_point (aes(x=x, y=y)) +
  xlab ("Attractiveness of parent") +
  ylab ("Percentage of girl babies") +
  ggtitle ("Regression line and \nposterior stimulations \ngiven strong prior")

df <- data.frame (a = sims[,1], b = sims[,2])
p4 <- ggplot (df) +
  geom_point (aes (x=a, y=b), size=0.01) +
  ggtitle("Posterior simulations \nunder strong prior") +
  xlim (30, 70) + ylim (-15, 15)

grid.arrange (p2, p1, p4, p3, ncol=2)  

# Ex 9-1
n <- 50
income <- sample (0:100, n, replace=T)
score  <- 20 + 0.5*income + rnorm (n, 0, 5)
x_new <- c(40, 80)
df_x_new <- data.frame (income=x_new)

df <- data.frame (income, score)

fit <- stan_glm(score ~ income, data = df, refresh=0)
print (fit)

point_est <- predict          (fit, newdata = df_x_new)
post_pred <- posterior_predict(fit, newdata = df_x_new)

mu_diff           <- point_est[2]  - point_est[1]
post_pred_diff    <- post_pred[,2] - post_pred[,1] 
quantile (post_pred_diff, c(0.05, 0.95)) 

# alternatively, brute force
aa                <- sort (abs (post_pred_diff - mu_diff))
ninety_percent_se <- aa[round (0.9*4000)]  # 36 
cat ("90% Predictive interval= [", 
     round(mu_diff-ninety_percent_se, 2), ",",
     round(mu_diff+ninety_percent_se, 2), "]")

posterior_graph_2 (x=income, y=score, x_new = x_new, x_lab = "Income(K)", y_lab = "Test score")

# Ex 9.3
hibbs <- read.table ("hibbs.dat", header = T) %>% rename (vote=inc.party.vote)
M1 <- stan_glm (vote ~ growth, data=hibbs, refresh=0)
print (M1)

growth_new <- 2
point_pred_a <- coef(M1)[1] + growth_new*coef(M1)[2] # 52.4

# a) 
sims <- as.matrix (M1)
sigma <- median (sims[,3]) # 3.9

sigma_linepred_a <- sigma *sqrt (1/nrow(hibbs) + (growth_new - mean (hibbs$growth))/
                                sum ((hibbs$growth - mean(hibbs$growth))^2)) # 1.0
sigma_postpred_a <- sigma *sqrt (1 + 1/nrow(hibbs) + (growth_new - mean (hibbs$growth))/
                                sum ((hibbs$growth - mean(hibbs$growth))^2)) # 4.0

point_pred_a <- coef(M1)[1] + 2*coef(M1)[2] # 52.4
# b)
new <- data.frame (growth=growth_new)
point_pred_b <- predict         (M1, newdata = new) # 52.4
sigma_linepred_b  <- mad (posterior_linpred (M1, newdata = new)) 1.0
sigma_postpred_b  <- mad (posterior_predict (M1, newdata = new)) 4.0

# Ex 9-4
theta_prior  <- 0.42
theta_prior_se <- 0.05
theta_data <- 0.54
theta_posterior <- 0.49

prec_prior <- 1/theta_prior_se^2
# theta_posterior <- (prec_prior*theta_prior + prec_data*theta_data)/(prec_prior+prec_data)
prec_data <- prec_prior * (theta_prior - theta_posterior)/(theta_posterior - theta_data)

data_se <- sqrt(1/prec_data) # 0.04
#data_se <- 0.5/sqrt(n)
n <- (0.5/data_se)^2 # 140

# Ex 9-5
theta_prior    <- -0.02
theta_data     <-  0.16
theta_prior_se <-  0.05
theta_data_se  <-  0.08

prec_prior <- 1/theta_prior_se^2
prec_data  <- 1/theta_data_se^2

theta_posterior <- (prec_prior*theta_prior + prec_data*theta_data)/(prec_prior+prec_data) # 0.03
prec_posterior <- prec_prior+prec_data
theta_posterior_se <- sqrt (1/prec_posterior) # 0.04

# Ex 9-9
n <- 50
a <- 10
b <- 0.85
N <- 1000

# model 1
midterm <- rnorm (n, 70, 20)
final   <- a + b*midterm + rnorm (n, 0, 20)
df <- data.frame(midterm, final)
fit <- stan_glm(final ~ midterm, data = df, refresh=0)
avg_diff <- mean (final) - mean (midterm)
ggplot () + geom_point(data = df, aes(x=midterm, y=final), size=0.5) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2]) +
  ggtitle (paste("avg(final) - avg(midterm) =", round (avg_diff, 1), 
                 "\nCorr (midterm, final) =", round(cor(midterm, final), 2)))

# alternative model
true_ability <- rnorm (n, 50, 20)
noise_1 <- rnorm (n, 0, 10)
noise_2 <- rnorm (n, 0, 10)
midterm <- true_ability + noise_1
final   <- true_ability + noise_2
df <- data.frame(midterm, final)
fit_alt <- stan_glm(final ~ midterm, data = df, refresh=0)
avg_diff <- mean (final) - mean (midterm)
ggplot () + geom_point(data = df, aes(x=midterm, y=final), size=0.5) +
  geom_abline(intercept = coef(fit_alt)[1], slope = coef(fit_alt)[2]) +
  ggtitle (paste("avg(final) - avg(midterm) =", round (avg_diff, 1), 
                 "\nCorr (midterm, final) =", round(cor(midterm, final), 2)))

find_prior_info <- function (){
  midterm <- rnorm (n, 70, 20)
  final   <- a + b*midterm + rnorm (n, 0, 20)
  prior_info <- list (diff_avg=abs(mean (final) - mean (midterm)), cor_m_f=100*cor (midterm, final))
  return (prior_info)
}
prior_info <- replicate (N, find_prior_info(), simplify = "array")
prior_info_unlisted <- array (unlist(prior_info), c(2, 1000) )
cat ("average (avg_final-avg_midterm) =",  mean (prior_info_unlisted[1,]))
cat ("average cor of midterm and final =", mean (prior_info_unlisted[2,]))

find_prior_info_alt <- function (){
  true_ability <- rnorm (n, 50, 20)
  noise_1 <- rnorm (n, 0, 10)
  noise_2 <- rnorm (n, 0, 10)
  midterm <- true_ability + noise_1
  final   <- true_ability + noise_2
  prior_info <- list (diff_avg=abs(mean (final) - mean (midterm)), cor_m_f=100*cor (midterm, final))
  return (prior_info)
}

prior_info_alt <- replicate (N, find_prior_info_alt())
prior_info_alt_unlisted <- array (unlist(prior_info_alt), c(2, 1000) )
cat ("average (avg_final-avg_midterm) =",  mean (prior_info_alt_unlisted[1,]))
cat ("average cor of midterm and final =", mean (prior_info_alt_unlisted[2,]))

max (prior_info_unlisted[1,])
max (prior_info_alt_unlisted[1,])


#########
# Chapter-10
kidiq <- as_tibble(read.dta(file="kidiq.dta"))
kidiq$mom_hs <- as.factor(kidiq$mom_hs)

fit <- stan_glm (kid_score ~ mom_hs, data = kidiq, refresh=0)
print (fit)

ggplot (kidiq) +
  geom_jitter(aes (x=mom_hs, y=mom_hs), width = 0.05) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2]) 
               
fit2 <- stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit2)

ggplot (kidiq) +
  geom_point(aes (x=mom_iq, y=kid_score)) +
  geom_abline(intercept = coef(fit2)[1], slope = coef(fit2)[2]) 

fit3 <- stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq, refresh=0)
print (fit3)

ggplot () +
  geom_point (data = kidiq, aes (x=mom_iq, y=kid_score, color = mom_hs)) +
  geom_abline(intercept = coef(fit3)[1], slope = coef(fit3)[3], color = "red") +
  geom_abline(intercept = coef(fit3)[1] + coef(fit3)[2], slope = coef(fit3)[3], color = "blue") +
  scale_color_manual(breaks = c(0, 1), values = c("red", "blue")) +
  ggtitle ("Kid score vs mom's IQ and mom hs graduation status\nData and fitted lines") +
  xlab ("Mother IQ score") +
  ylab ("Child test score")
  
fit4 <- stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, data = kidiq, refresh = 0)
print (fit4)

ggplot () +
  geom_point (data = kidiq, aes (x=mom_iq, y=kid_score, color = mom_hs)) +
  geom_abline(intercept = coef(fit4)[1], slope = coef(fit4)[3], color = "red") +
  geom_abline(intercept = coef(fit4)[1] + coef(fit4)[2], slope = coef(fit4)[3] + coef(fit4)[4], color = "blue") +
  scale_color_manual(breaks = c(0, 1), values = c("red", "blue")) +
  ggtitle ("Kid score vs mom's IQ and mom hs graduation status\nData and fitted lines") +
  xlab ("Mother IQ score") +
  ylab ("Child test score")

# Earnings_1 data
colnames (earnings_1)

fit <- stan_glm (weight ~ c_height + male, data = earnings_1, refresh = 0)
print (fit)
 
# weight of a 70 in tall woman
woman <- data.frame (c_height=(70-mean(earnings_1$height)), male=0)

woman_weight <- posterior_predict(fit, woman)
cat ("Average weight for a woman of 70in =", median (woman_weight), "with sd =", mad (woman_weight))

# weight of a 70 in tall man
man <- data.frame (c_height=(70-mean(earnings_1$height)), male=1)

man_weight <- posterior_predict(fit, man)
cat ("Average weight for a man of 70in =", median (man_weight), "with sd =", mad (man_weight))

fit2 <- stan_glm (weight ~ c_height + male + factor(ethnicity), data = earnings_1, refresh = 0)
print (fit2)

# baseline category in ethnicity is "Black"
#Median MAD_SD
#(Intercept)               156.6    2.3 
#c_height                    3.9    0.3 
#male                       12.0    2.1 
#factor(ethnicity)Hispanic  -6.1    3.6 
#factor(ethnicity)Other    -12.3    5.3 
#factor(ethnicity)White     -5.2    2.3 

# Change baseline to White 
earnings_1$eth <- factor ((earnings_1$ethnicity), levels =c("White", "Black", "Hispanic", "Other"))

fit3 <- stan_glm (weight ~ c_height + male + eth, data = earnings_1, refresh = 0)
print (fit3)

# baseline category is now White
#Median MAD_SD
#(Intercept) 151.4    1.1 
#c_height      3.9    0.2 
#male         12.1    1.9 
#ethBlack      5.1    2.2 
#ethHispanic  -1.0    3.0 
#ethOther     -6.8    4.8 

# Congressional elections
c_86 <- as_tibble(read.table("congress_86", header=F)) %>%
  mutate (percent_dem_86 = V4/(V4+V5), 
  percent_dem_86 = if_else(percent_dem_86 > 0.99, 0.75, percent_dem_86),
  percent_dem_86 = if_else(percent_dem_86 < 0.01, 0.25, percent_dem_86)) %>%
  rename  (inc_86 = V3) %>% 
  filter (V4 != -9 | V5 != -9)  %>%
  select (-V4, -V5)
c_88 <- as_tibble(read.table("congress_88", header=F)) %>%
  mutate (percent_dem_88 = V4/(V4+V5)) %>% 
  filter (percent_dem_88 < 0.99 & percent_dem_88 > 0.01) %>%
  rename  (inc_88 = V3) %>% 
  filter (V4 != -9 | V5 != -9) %>%
  select (-V4, -V5)
c_90 <- as_tibble(read.table("congress_90", header=F)) %>%
  mutate (percent_dem_90 = V4/(V4+V5)) %>% 
  rename  (inc_90 = V3) %>% 
  filter (V4 != -9 | V5 != -9) %>%
  select (-V4, -V5)

congr_dat <- c_86 %>% 
  left_join(c_88, by=c("V1", "V2")) %>%
  left_join(c_90, by=c("V1", "V2")) %>%
  na.omit()

ggplot (congr_dat) +
  geom_point(aes(x=percent_dem_86, y=percent_dem_88, color=factor(inc_88))) +
  scale_color_manual(breaks = c(-1, 0, 1), 
                     values = c("red", "green", "blue"),
                     labels = c("R", "D", "open")) +
  geom_abline(intercept = 0, slope = 1) +
  xlim (0,1) + ylim (0,1)
  
df_88 <- data.frame (vote      = congr_dat$percent_dem_88, 
                     past_vote = congr_dat$percent_dem_86, 
                     inc       = congr_dat$inc_88)
fit_88 <- stan_glm (vote ~ past_vote + inc, data = df_88, refresh=0)  
print (fit_88, digits=2)

# predict 1990 results from 1988
df_90 <- data.frame (past_vote = congr_dat$percent_dem_88, 
                     inc       = congr_dat$inc_90)
pred_90 <- posterior_predict (fit_88, newdata = df_90)
dim (pred_90)
# 4000 simulations for the vote Dem percentage of each congressional seat
# each column represents 4000 sims for a congressional seat
# each row represents a simulation for each seat
dems_pred <- rowSums (pred_90 > 0.5) # 
median (dems_pred) # 197 on average Dems were expected to win 197 of 347 seats
mad (dems_pred)

pre_90_sd <- apply (pred_90, 2, sd) 
each_distr_avg <- apply (pred_90, 2, mean)
hist (each_distr_avg)
sum (each_distr_avg > 0.5) # 198, number of seats Dems expected to win out od 347 in 1990 
each_distr_se  <- apply (pred_90, 2, sd) #  sd of outcomes foe each district
mean (each_distr_se)  # 7%, average se in each congressional district 

# Ex 10-1
N <- 100
sigma <- 3
 
z <- rbinom (N, 1, 0.5)
x <- rnorm (N, z, 1)
error <- rnorm (N, 0, sigma)
b <- c(1, 2, -1, -1)
y <- b[1] + b[2]*x + b[3]*z + b[4]*x*z + error

df <- data.frame (x, y, z)

fit_1 <- stan_glm (y ~ x, data=df, subset=(z==0), refresh=0)
print (fit_1)
fit_2 <- stan_glm (y ~ x, data=df, subset=(z==1), refresh=0)
print (fit_2)

ggplot (df) + 
  geom_point (aes(x=x, y=y, color=factor(z))) +
  geom_abline(intercept = coef(fit_1)[1], slope = coef(fit_1)[2], color="red") +
  geom_abline(intercept = coef(fit_2)[1], slope = coef(fit_2)[2], color="blue") +
  scale_color_manual(breaks = c(0, 1), values = c("red", "blue")) +
  ggtitle ("Model without x:y interaction")
  
fit_3 <- stan_glm (y ~ x + x:y, data=df, subset=(z==0), refresh=0)
print (fit_3)
fit_4 <- stan_glm (y ~ x + x:y, data=df, subset=(z==1), refresh=0)
print (fit_4)

ggplot (df) + 
  geom_point (aes(x=x, y=y, color=factor(z))) +
  geom_abline(intercept = coef(fit_3)[1], slope = coef(fit_3)[2], color="red") +
  geom_abline(intercept = coef(fit_4)[1], slope = coef(fit_4)[2], color="blue") +
  scale_color_manual(breaks = c(0, 1), values = c("red", "blue")) +
  ggtitle ("Model with x:y interaction")

# Ex 10-2
# y <- 1.2 + 1.6*x + 2.7*z + 0.7*x*z + normal (0, 0.5)
# control, z=0
# y_control <- 1.2 + 1.6*x + error
# treat, z=1
# y_treat   <- 3.9 + 2.3*x + error

N <- 100
z <- rbinom (N, 1, 0.5)
x <- runif (N, 0, 10)
y <-1.2 + 1.6*x + 2.7*z + 0.7*x*z + rnorm (N, 0, 0.5)

df <- data.frame(x, y, z)
ggplot (df) +
  geom_point(aes(x=x, y=y, color=factor(z))) +
  geom_abline(intercept = 1.2, slope=1.6, color="red") +
  geom_abline(intercept = 3.9, slope=2.3, color="blue") +
  scale_color_manual(breaks = c(0, 1), values = c("red", "blue")) +
  ggtitle("y = 1.2 + 1.6*x + 2.7*z + 0.7*x*z + normal (0, 0.5)")
  

# Chapter 11
kidiq <- as_tibble(read.dta(file="kidiq.dta")) %>% 
  mutate (mom_hs = as.factor (mom_hs))

mom_hs_bar <- mean (as.numeric(as.character(kidiq$mom_hs)))
mom_iq_bar <- mean (kidiq$mom_iq)

kidiq <- kidiq %>% 
# kid_score v mom_iq
fit_1 <- stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit_1) 
sims <- as.matrix (fit_1)
pick_sims <- sample(4000, 10)
ggplot (kidiq) +
  geom_abline(intercept = sims[pick_sims, 1], slope = sims[pick_sims, 2], color="gray") +
  geom_point (aes(x=mom_iq, y=kid_score)) +
  geom_abline(intercept = coef (fit_1)[1], slope = coef (fit_1)[2]) +
  xlab ("IQ of mom") + ylab ("Kid's test score") +
  ggtitle ("Kid's test score vs Mom's IQ, Best linear fit, 10 draws from simulations")
# kid_score v mom_hs + mom_iq - no interction
fit_3 <- stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq, refresh = 0)
print (fit_3)  
p3 <- ggplot (kidiq) +
  geom_point (aes(x=mom_iq, y=kid_score, color = mom_hs)) +
  geom_abline(intercept = coef(fit_3)[1], slope = coef(fit_3)[3], color = "red") +
  geom_abline(intercept = coef(fit_3)[1]+coef(fit_3)[2], slope = coef(fit_3)[3], color = "blue") +
  scale_color_manual(breaks = c(0,1), values = c("red", "blue")) +
  xlab ("Mom's IQ") + ylab ("Kid's Score") +
  ggtitle ("No interaction")
# kid_score v mom_hs + mom_iq - no interction
fit_4 <- stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, data = kidiq, refresh = 0)
print (fit_4)  
p4 <- ggplot (kidiq) +
  geom_point (aes(x=mom_iq, y=kid_score, color = mom_hs)) +
  geom_abline(intercept = coef(fit_4)[1], slope = coef(fit_4)[3], color = "red") +
  geom_abline(intercept = coef(fit_4)[1]+coef(fit_4)[2], 
              slope = coef(fit_4)[3] +coef(fit_4)[4], color = "blue") +
  scale_color_manual(breaks = c(0,1), values = c("red", "blue")) +
  xlab ("Mom's IQ") + ylab ("Kid's Score") +
  ggtitle ("With interaction")

grid.arrange (p3, p4, ncol=1)  
# Display one plot per input variable, use fit_3
# kid_score = b1 + b2*mom_hs + b3*mom_iq
fit_3 <- stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq, refresh = 0)
print (fit_3) 
sims_3 <- as.matrix (fit_3)
b_hat <- apply (sims_3, 2, median)
# kid_score = (b1 + b2*mom_hs_bar) + b3*mom_iq
p1 <- ggplot (kidiq) +
  geom_abline(intercept = sims_3[pick_sims, 1] + mom_hs_bar*sims_3[pick_sims, 2], 
              slope = sims_3[pick_sims, 3], color="gray") +
  geom_point (aes(x=mom_iq, y=kid_score)) +
  geom_abline(intercept = b_hat[1] + mom_hs_bar*b_hat[2], slope = b_hat[3]) +
  xlab ("IQ of mom") + ylab ("Kid's test score") +
  ggtitle ("Kid's test score vs Mom's IQ | Mom'a HS = average (0.78)\nBest linear fit, 10 draws from simulations")
# kid_score = (b1+ b3*mom_iq_bar) + b2*mom_hs 
p2 <- ggplot (kidiq) +
  geom_abline(intercept = sims_3[pick_sims, 1] + mom_iq_bar*sims_3[pick_sims, 3], 
              slope = sims_3[pick_sims, 2], color="gray") +
  geom_jitter (aes(x=mom_hs, y=kid_score), width = 0.05) +
  geom_abline(intercept = b_hat[1] + mom_iq_bar*b_hat[3], slope = b_hat[2])+
  xlab ("HS of mom") + ylab ("Kid's test score") +
  ggtitle ("Kid's test score vs Mom's HS | Mom's IQ = average (100)\nBest linear fit, 10 draws from simulations")

grid.arrange (p1, p2, ncol=2)

# Multiple Regression
# y=a +b1*x1 + b2*x2 + ... + bk*xk + theta*z + error

N <- 100
K <- 10
X <- array (runif (N*K, 0, 1), c(N,K))
z <- sample (c(0, 1), N, replace = TRUE)
a <- 1
b <- 1:K
theta <- 5
sigma <- 2

y <- a + X %*% b + theta*z + rnorm (N, 0, sigma)

df <- data.frame (X=X, y=y, z=z)
fit <- stan_glm (y ~ X + z, data = df, refresh = 0)
print (fit)

new <- data.frame (X=X)
y_hat <- predict (fit, newdata = new)

df <- cbind (df, y_hat)

p1 <- ggplot (subset (df, z==0)) +
  geom_point(aes(x=y, y=y_hat)) +
  geom_abline(intercept = 0, slope=1) + 
  xlim (10, 50) + ylim (10, 50) +
  xlab ("linear predictor, y_hat") +
  ylab ("Outcome, y") +
  ggtitle ("z=0")
  
p2 <- ggplot (subset (df, z==1)) +
  geom_point(aes(x=y, y=y_hat)) +
  geom_abline(intercept = 0, slope=1) + 
  xlim (10, 50) + ylim (10, 50) +
  xlab ("linear predictor, y_hat") +
  ylab ("Outcome, y") +
  ggtitle ("z=1")

max_y <- max (y-y_hat)
min_y <- min (y-y_hat)
p3 <- ggplot (subset (df, z==0)) +
  geom_point(aes(x=y, y=y-y_hat)) +
  geom_abline(intercept = 0, slope=0) + 
  xlim (10, 50) + ylim (min_y, max_y) +
  xlab ("linear predictor, y_hat") +
  ylab ("Residual, y-y_hat") +
  ggtitle ("z=0")

p4 <- ggplot (subset (df, z==1)) +
  geom_point(aes(x=y, y=y-y_hat)) +
  geom_abline(intercept = 0, slope=0) + 
  xlim (10, 50) + ylim (min_y, max_y) +
  xlab ("linear predictor, y_hat") +
  ylab ("Residual, y-y_hat") +
  ggtitle ("z=1")

grid.arrange (p1, p2, p3, p4, ncol=2)

# More residual plots
introclass <- read.table("introclass", header=TRUE)
fit <- stan_glm (Final ~ Midterm, data = introclass, refresh = 0)
print (fit)

introclass$Final_pred <- predict (fit)

introclass$resid <- introclass$Final - introclass$Final_pred

p1 <- ggplot (introclass) +
  geom_point (aes(x=Final_pred, y=resid)) +
  geom_abline(intercept = 0, slope=0) +
  xlab ("Final - Predicted") + ylab ("Residual") +
  ggtitle("Residual vs Predicted Final")
p2 <- ggplot (introclass) +
  geom_point (aes(x=Final, y=resid)) +
  geom_abline(intercept = 0, slope=0) +
  xlab ("Final - Actual") + ylab ("Residual") +
  ggtitle("Residual vs Actual Final")
grid.arrange (p1, p2, ncol=2)
#
# Plot residuals vs Predicted points!!!
# Not residuals vs data!!!

# Generate fake data from output of fit
a <- coef (fit)[1]
b <- coef (fit)[2]
sigma <- sigma (fit)
n <- nrow (introclass)
introclass$final_fake <- a + b*introclass$Midterm + rnorm (n, 0, sigma)
# fit model to data
fit_fake <- stan_glm (final_fake ~ Midterm, data = introclass, refresh = 0) 

sims <- as.matrix (fit_fake)
introclass$predicted_fake <- colMeans(sims[,1] + sims[,2] %*% t(introclass$Midterm))
# or, alternatively
predicted_fake_2 <- predict (fit_fake)
introclass$resid_fake <- introclass$predicted_fake - introclass$final_fake

p1 <- ggplot (introclass) +
  geom_point (aes(x=predicted_fake, y=resid_fake)) +
  geom_abline(intercept = 0, slope=0) +
  xlab ("Final - Predicted") + ylab ("Residual") +
  ggtitle("Fake Residual vs Predicted Final")
p2 <- ggplot (introclass) +
  geom_point (aes(x=final_fake, y=resid_fake)) +
  geom_abline(intercept = 0, slope=0) +
  xlab ("Final - Actual") + ylab ("Residual") +
  ggtitle("Fake Residual vs Actual Final")
grid.arrange (p1, p2, ncol=2)
#
# Plot residuals vs Predicted points!!!
# Not residuals vs data!!!

# Predictive checking!
newcomb <- data.frame(y=c(28,26,33,24,34,-44,27,16,40,-2,29,22,24,21,25,30,
                          23,29,31,19,24,20,36,32,36,28,25,21,28,29,37,25,28,
                          26,30,32,36,26,30,22, 36,23,27,27,28,27,31,27,26,
                          33,26,32,32,24,39,28,24,25,32,25,29,27,28,29,16,23)) 
fit <- stan_glm (y ~ 1, data = newcomb, refresh = 0)
print (fit)
hist (newcomb$y, breaks = 100)
# generate data from fit
n <- dim (newcomb)[1]
sigma <- sigma (fit)
sims  <- as.matrix (fit)
n_sims <- dim (sims)[1]

y_rep <- array (NA, c(n_sims, n))
for (s in 1:n_sims) 
  y_rep[s,] <- rnorm (n, sims[s,1], sims[s,2]) 
colMeans(y_rep)
# alternatively
y_rep_3 <- array(rnorm (n_sims*n, sims[,1], sims[,2]), c(n_sims, n))
colMeans(y_rep_3)
# alternatively
y_rep_2  <- posterior_predict (fit)
colMeans(y_rep_2)

par(mfrow=c(5,4))  
for (s in sample(n_sims, 20)) hist(y_rep[s,])   
par(mfrow=c(1,1))  

# pick 20 samples
samples <- sample(n_sims, 20)

# ggplot
dft <- data.frame (t(y_rep[samples,])) 
colnames (dft) <- samples
dft_long <- dft %>%  pivot_longer (everything(), names_to = "sims", values_to = "y_rep")
ggplot (dft_long) + 
  geom_histogram(aes(x=y_rep), bins=20, color="black", fill="gray") + 
  facet_wrap (~ as.numeric(sims), ncol=5) +
  ggtitle ("Normal fit to Newcomb data: 20 simulated sets of 66 observations")

# Test metric
dft_all <- data.frame (t(y_rep)) 
y_rep_mins <- data.frame(mins=apply (dft_all, 2, min))
ggplot (y_rep_mins) + 
  geom_histogram(aes(x=mins), bins=30, color="black", fill="gray") +
  geom_vline(xintercept = min(newcomb$y), size=1) +
  xlab ("Mimimum of Newcomb data and minima of replicated values from model") +
  ggtitle ("Test statistic for normal model fitted to Newcomb data")
  
#####
# time series model - autoregression
# this year's unemployment is related to last year's
#####
unemp <- as_tibble (read.table("unemployment", header=TRUE))
n <- dim (unemp)[1]
unemp$y_lag <- c(NA, unemp$y[1:(n-1)])

fit_lag <- stan_glm (y ~ y_lag, data = unemp, refresh = 0)
print (fit_lag)

sim_lag <- as.matrix(fit_lag)

sigma <- median (sim_lag[, 3])
# or alternatively
sigma <- sigma (fit_lag)

n_sims <- dim (sim_lag)[1]

y_rep <- array (NA, c(n_sims, n))
for (s in 1:n_sims){
  y_rep[s, 1] <- unemp$y[1]
  for (t in 2:n){
    y_rep[s,t] <- sim_lag[s,1] + sim_lag[s,2]*y_rep[s, t-1] + rnorm (1, 0, sigma)
  }
}
samples <- sample(4000, 20)
y_rep_samples <- y_rep[samples,]
dft <- data.frame(t(y_rep_samples))
colnames (dft) <- samples

dft_long <- cbind (year=unemp$year, dft) %>% 
  pivot_longer(-year, names_to = "sims", values_to = "y_reps")

ggplot (dft_long) +
  geom_line(aes(x=year, y=y_reps)) + 
  facet_wrap( ~ sims, ncol=5)

# test statistic
# figure out the number of "turns" data and each simulation
# i.e.: deriv at t+1 != deriv at t
# i.e.: unemp(t)-unemp(t-1) != unemp(t-1)-unemp(t-2)

num_switches <- function (y){
 n <- length (y) 
  y_lag   <- c(NA, y[1:(n-1)])
  y_lag_2 <- c(NA, NA, y[1:(n-2)])
  n_switch <- sum (sign (y-y_lag) != sign (y_lag-y_lag_2), na.rm = TRUE)
  return (n_switch)
}

n_switches_data <- num_switches (unemp$y)
cat (n_switches_data, "direction switches in unemployment data")

# do same for simulations
ss <- rep (NA, n_sims)
for (s in 1:n_sims){
  ss[s] <- num_switches(y_rep[s,])
}

df_ss <- data.frame (switches=ss)
ggplot (df_ss) + 
  geom_histogram(aes(x=switches), bins=100, color="black", fill="gray") +
  geom_vline(xintercept = n_switches_data, size=1) +
  xlab ("Number of derivative changes") +
  ggtitle ("Number of direction changes\nData, (vertical line) and 4000 simulations")

cat ("!!!", round (mean(ss > n_switches_data)*100, 0), "of the simulations have more derivative changes than data!!!\n")
aa <- quantile (ss, c(0.10, 0.90), names=FALSE)
cat ("80% of the data is between", aa[1], "and", aa[2])

##########
# Residuals R2 = var(y_hat)/var(y)
##########
x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3

df <- data.frame(x, y)
fit <- stan_glm (y ~ x, data=df, refresh=0)
print (fit)
sims <- as.matrix(fit)
sigma <- median(sims[,3])

ggplot () + geom_point(aes(x=x, y=y)) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])

y_hat <- coef(fit)[1] + coef(fit)[2]*x
r <- y-y_hat
sd_y     <- sd (y)
sd_y_hat <- sd (y_hat)
sd_r     <- sd (r)
R2       <- sd_y_hat^2/sd_y^2
cat ("sd(y)=", sd_y, ", sd (y_hat)=", sd_y_hat, ", sd (r)=", sd_r, "\nR2=", R2)

# Refit with STRONG priors:
# beta_0_hat ~ N (0, 0.2)
# beta_1_hat ~ N (1, 0.2)

a_hat_prior <- 0
b_hat_prior <- 1
fit_strong_prior <- stan_glm(y ~ x, data = df, refresh = 0,
                             prior_intercept = normal (a_hat_prior, 0.2),
                             prior = normal (b_hat_prior, 0.2))
print (fit_strong_prior)
sims <- as.matrix(fit_strong_prior)
a_hat_post <- median (sims[,1])
b_hat_post <- median (sims[,2])
sigma_hat_post <- median (sims[,3])


y_hat_post <-  a_hat_post + b_hat_post*x

sd_y_hat_post <- sd (y_hat_post)
r_post <- y-y_hat_post
sd_r_post     <- sd (r_post)
R2_post       <- sd_y_hat_post^2/sd_y^2
cat ("sd(y)=", sd_y, ", sd (y_hat_post)=", sd_y_hat_post, ", sd (r_post)=", sd_r_post, 
     "\nR2_post=", R2_post, "!!!! R > 1 !!!!!")

s_20 <- sample(4000, 20)

cbind (y_hat_post, df)
ggplot (df) + 
  geom_abline(intercept = sims[s_20,1],   slope = sims[s_20,2],   color="blue",   size=0.2) +
  geom_point(aes(x=x, y=y,          color="black")) +
  geom_point(aes(x=x, y=y_hat_post, color="blue"), size=2) +
  geom_abline(aes(intercept = coef(fit)[1], slope = coef(fit)[2], color="black"), size=1) +
  geom_abline(aes(intercept = a_hat_prior,  slope = b_hat_prior,  color="red"),   size=1, linetype="dashed") +
  geom_abline(aes(intercept = a_hat_post,   slope = b_hat_post,   color="blue"),  size=1) +
  scale_color_identity(name = "Data and lines",
                       breaks = c("black", "red","blue"),
                       labels = c("Least SquaresFit", "Prior Line", "PosteriorFit"),
                       guide = "legend") +
  ggtitle ("Least Squares Fit, Prior, Posterior Fit and 20 posterior simulations")

# Alternative definition: Bayesian R2 for each simulation
# R2B_s <- var(y_hat_s)/(var(y_hat_s) + sigma_s^2

# point estimate 
sd (y_hat_post)^2/(sd (y_hat_post)^2 + sigma_hat_post^2) # 78

# general case
R2B_1_bar <- round (median (bayes_R2(fit_strong_prior)), 2) # 78
R2B_1_sd  <- round (sd (bayes_R2(fit_strong_prior)), 2) # 13
R2B_1     <- bayes_R2(fit_strong_prior)

ggplot (data.frame(R2B_1)) +
  geom_histogram(aes(x=R2B_1), bins=50, color ="black", fill="gray") +
  xlab("Bayesian R2") +
  ggtitle ("Bayesian R2 for 4000 simulations - mean R2 =", R2B_1_bar)

# alternatively - DIY
y_reps <- posterior_predict (fit_strong_prior)

R2B      <- rep (NA, 4000)
for (s in 1:4000){
  R2B[s]  <- sd(y_reps[s,])^2/(sd(y_reps[s,])^2 + sims[s,3]^2)
}

R2B_bar <- round(mean(R2B), 2)
R2B_sd <- sd (R2B) 
cat ("Bayesian R2 has mean=", R2B_bar, "sd=", round (R2B_sd, 2) )

ggplot (data.frame(R2B)) +
  geom_histogram(aes(x=R2B), bins=50, color ="black", fill="gray") +
  xlab("Bayesian R2") +
  ggtitle ("Bayesian R2 for 4000 simulations - mean R2 =", R2B_bar)

# Ex 11-5
pyth <- read.table("pyth", header=TRUE)
pyth_data <- na.omit(pyth)

fit <- stan_glm (y ~ x1 + x2, data = pyth_data, refresh = 0)
print (fit)  
aa <- range (pyth_data$x1)
bb <- range (pyth_data$x2)
x2_bar <- mean (pyth_data$x2)
pyth_data$y_fit <- coef(fit)[1] + coef(fit)[2]*pyth_data$x1 + coef(fit)[3]*pyth_data$x2
# alternatively
pyth_data$y_fit <- predict (fit)
pyth_data$fit_error <- pyth_data$y - pyth_data$y_fit

ggplot (subset (pyth, is.na(y)==FALSE)) +
  geom_point (aes(x=x1, y=y)) +
  geom_abline(aes(intercept = coef(fit)[1] + coef(fit)[3]*bb[1],  slope = coef(fit)[2], color="red")) +
  geom_abline(aes(intercept = coef(fit)[1] + coef(fit)[3]*x2_bar, slope = coef(fit)[2], color="black")) +
  geom_abline(aes(intercept = coef(fit)[1] + coef(fit)[3]*bb[2],  slope = coef(fit)[2], color="blue")) +
  scale_color_identity(name = "",
                     breaks = c("blue","black", "red"),
                     labels = c("max x2", "avg x2", "min x2"),
                     guide = "legend") +
  ggtitle ("Model : y = b0 + b1*x1 + b2*x2")
  
ggplot (pyth_data) +
  geom_point (aes(x=y, y=y_fit)) + 
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("Predicted y vs actual y")

ggplot (pyth_data) +
  geom_point (aes(x=y_fit, y=fit_error)) + 
  geom_abline(intercept = sigma (fit), slope = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 0) +
  geom_abline(intercept = -sigma (fit), slope = 0, linetype = "dashed") +
  ggtitle("Prediction error \U00B1 1 se vs Predicted y") +
  xlab ("Predicted y") + ylab ("Prediction error")

new_pyth <- subset (pyth, is.na(y)==TRUE)  %>% select (x1, x2)

predict (fit, newdata = new_pyth)
  
ggplot () +
  geom_point (data=pyth_data, aes(x=x1, y=x2, color="black")) +
  geom_point (data=new_pyth,  aes(x=x1, y=x2, color="red")) +
  scale_color_identity (name = "",
                       breaks = c("black", "red"),
                       labels = c("Sample", "Inference"),
                       guide = "legend") +
  ggtitle("Range of predictors for sample and predictors for inference\nGood overlap")

  
########
# Ch 12
#######
earnings <- as_tibble (read.csv ("Earnings", header=T))  
fit <- stan_glm (earn ~ height, data = earnings, refresh = 0 )
print (fit)  
cat ("intercept = $", coef(fit)[1], ", slope =", coef(fit)[2], "$/in , sigma =", sigma (fit))
sims <- as.matrix (fit)
select_20 <- sample (4000, 20)
ggplot (earnings) +
  geom_abline(intercept = sims[select_20,1], slope = sims[select_20,2], color="gray") +
  geom_abline(intercept = median(sims[,1]), slope = median(sims[,2]), color="black") +
  geom_jitter (aes(x=height, y=earn), width=0.1) +
 ylim (0, 200001)

# Center kidiq
kidiq <- as_tibble(read.dta(file="kidiq.dta")) %>%
  mutate (c_mom_iq = mom_iq - mean (mom_iq),
          c_mom_hs = mom_hs - mean (mom_hs))
fitc <- stan_glm (kid_score ~ c_mom_hs + c_mom_iq + c_mom_hs:c_mom_iq, data = kidiq, refresh = 0)
print (fitc)

# Standardize kiqiq:  z_var = (var - mean(bar))/(2*sd(var))
kidiq <- kidiq %>%
  mutate (z_mom_iq = (mom_iq - mean (mom_iq))/(2*sd(mom_iq)),
          z_mom_hs = (mom_hs - mean (mom_hs))/(2*sd(mom_hs)))
fitz <- stan_glm (kid_score ~ z_mom_hs + z_mom_iq + z_mom_hs:z_mom_iq, data = kidiq, refresh = 0)
print (fitz)

# log transformation: log (y) = b0 + b1*x1 + b2*x2 + ....
fit_log1 <- stan_glm (log(earn) ~ height, data = subset(earnings, earn > 0) , refresh = 0)
print (fit_log1, digits=2)

sims <- as.matrix (fit_log1)
select_20 <- sample (4000, 20)
# log scale
ggplot (subset(earnings, earn > 0)) +
  geom_abline(intercept = sims[select_20,1], slope = sims[select_20,2], color="gray") +
  geom_abline(intercept = median(sims[,1]), slope = median(sims[,2]), color="black") +
  geom_jitter (aes(x=height, y=log(earn)), width=0.1) +
  ylim (5, log(200001))

# true scale
ggplot (subset(earnings, earn > 0)) +
#  geom_function(fun = function(earn) exp(sims[select_20,1] + sims[select_20,2]*earn)) +
  geom_function(fun = function(earn) exp(median(sims[,1]) + median(sims[,2])*earn)) +
  geom_jitter (aes(x=height, y=earn), width=0.1) +
  ylim (0, 200001)

# POSTERIOR PREDICTIVE CHECKING: Did the model do a good job fitting data?
# first linear fit
dim (y_rep_log1) # 4000 X 1629
pick_100 <- sample (dim(y_rep_log1)[1], 100)

y_rep_1 <- posterior_predict (fit)
ppc_dens_overlay(earnings$earn, y_rep_1[pick_100,])
# not great.. predictions have negative earnings

# try the log(earn fit)
y_rep_log1 <- posterior_predict (fit_log1)
dim (y_rep_log1) # 4000 X 1629
pick_100 <- sample (dim(y_rep_log1)[1], 100)

ppc_dens_overlay(log(subset(earnings, earn > 0)$earn), y_rep_log1[pick_100,])

## Add another predictor: male
fit_log2 <- stan_glm (log(earn) ~ height + male, data = subset(earnings, earn > 0) , refresh = 0)
print (fit_log2, digits=3)
# coefficient of male is 0.37..  Means a man of same height expected to make
# exp(0.37)=1.45, about 45% times more!
# coefficient of height falls to 0.02, Means a woman 1 in taller  expect to make
# exp(0.02)=1.02, about 2% more
#
# log earning:  a 70 in person expected to have 7.99 + 70*0.024 = 9.4
# log CI.. sigma =0.87, 1 sd CI -> 9.4 +/- 0.9 = [8.5 , 10.3] =>
# linear CI = [5k , 29k]  Huge!!!
R2 <- median (bayes_R2 (fit_log2)) #*****
ci_70_lo <- exp(coef(fit_log2)[1]+coef(fit_log2)[2]*70 - sigma(fit_log2))
ci_70_hi <- exp(coef(fit_log2)[1]+coef(fit_log2)[2]*70 + sigma(fit_log2))
cat ("R2=", round (R2*100, 0), "%  - 1 sd CI= [", round(ci_70_lo,0),",",round(ci_70_hi,0),"]")

# include interaction height:male
fit_log3 <- stan_glm (log(earn) ~ height + male + height:male, data = subset(earnings, earn > 0) , refresh = 0)
print (fit_log3, digits=3)
# hard to interpret intercept and coefficients when height is set to zero, 
# so center and standardize height
mean(earnings$height) # 66.6in
sd(earnings$height)   # 3.8in
earnings$z_height <- (earnings$height - mean(earnings$height))/sd(earnings$height)
fit_log4 <- stan_glm (log(earn) ~ z_height + male + z_height:male, data = subset(earnings, earn > 0) , refresh = 0)
print (fit_log4, digits=3)

#  INTERPRET coefficients
#              Median MAD_SD
# (Intercept)   9.545  0.037  when z_height and male=0: ie. 66.6in tall woman exp to make exp(9.54) = 14k
# z_height      0.057  0.043  when male is zero: a woman 3.8 in taller will make 6% more
# male          0.350  0.067  66.6 in male exp to make more than 66.6 woman by exp(0.35), 42% more
# z_height:male 0.076  0.063  a 3.8 in taller male is expected to make 8% more than a 3.8 in taller woman
#                             a 66.6+3.8 in man exp to make 6+8=14% more than a 66.6 in woman

####
# Ex: mesquite bushes
####
mesquite <- as_tibble (read.table("mesquite", header=TRUE))

fit_1 <- stan_glm (weight ~ diam1 + diam2 + canopy_height + total_height + density + group, 
                   data = mesquite, refresh = 0)
print (fit_1)
# evaluate model with LeaveOneOut
(loo_1 <- loo (fit_1))  # elpd_loo -335     ...error messages, use k-fold cv
kfold_1 <- kfold (fit_1, K=10)

# refit on log scale
fit_2 <- stan_glm (log(weight) ~ log(diam1) + log(diam2) + log(canopy_height) + 
                     log(total_height) + log(density) + group, 
                   data = mesquite, refresh = 0)
print (fit_2)
(loo_2 <- loo (fit_2)) # elpd_loo -19

# COMPARE USING PREDICTIVE CHECKS - ppc_dens_overlay!!!!
y_rep_1 <- posterior_predict (fit_1)
n_sims <- nrow (y_rep_1)
subset <- sample (n_sims, 100)
p_1 <- ppc_dens_overlay(mesquite$weight, y_rep_1[subset,]) + # weight goes negative, not good
  ggtitle ("Model 1: 6 predictors", subtitle="weight ~ diam1 + diam2 + height + total_height + density + group")

predict (fit_2)
y_rep_2 <- posterior_predict (fit_2)
n_sims <- nrow (y_rep_2)
subset <- sample (n_sims, 100)
p_2 <- ppc_dens_overlay(log(mesquite$weight), y_rep_2[subset,]) + # much better
  ggtitle ("Model 2: 6 predictors", subtitle="log(weight) ~ log(diam1) + log(diam2) + log(height) + 
                     log(total_height) + log(density) + group")

# try simpler model
mesquite <- mesquite %>%
  mutate (canopy_volume = diam1*diam2*canopy_height)
fit_3 <- stan_glm (log(weight) ~ log(canopy_volume), data = mesquite, refresh = 0)
print (fit_3)
(loo_3 <- loo (fit_3)) # edlp -26.5
(loo_compare(loo_2, loo_3)) # model2 is better has higher elpd

median (loo_R2(fit_2))  # 85%
median (loo_R2(fit_3))  # 78%   model2 better... but far fewer parameters!!

y_rep_3 <- posterior_predict (fit_3)
n_sims <- nrow (y_rep_3)
subset <- sample (n_sims, 100)
p_3 <- ppc_dens_overlay(log(mesquite$weight), y_rep_3[subset,]) + # still decent compare to model 2
  ggtitle ("Model 3: 1 predictor", subtitle="log(weight) ~ log(volume)")

# try some more stuff
mesquite <- mesquite %>%
  mutate (canopy_area  = diam1*diam2,
          canopy_shape = diam1/diam2)

fit_4 <- stan_glm (log(weight) ~ log(canopy_volume) + log(canopy_area) + log(canopy_shape) +
                     log(total_height) + log(density) + group, 
                   data = mesquite, refresh = 0)
print (fit_4)
(loo_4 <- loo (fit_4)) # edlp -19.5 similar to model2
(loo_compare(loo_2, loo_4)) # model2 is better has higher elpd
median (loo_R2(fit_4))      # 85% same as model 2 

# another
fit_5 <- stan_glm (log(weight) ~ log(canopy_volume) + log(canopy_area) + log(canopy_shape) + group, 
                   data = mesquite, refresh = 0)
print (fit_5)
(loo_5 <- loo (fit_5)) # edlp -18.7 a bit better than model2
(loo_compare(loo_2, loo_5)) # model5 is better has higher elpd
median (loo_R2(fit_5))  # 85% same as model2 but fewer predictors!!!

y_rep_5 <- posterior_predict (fit_5)
n_sims <- nrow (y_rep_5)
subset <- sample (n_sims, 100)
p_5 <- ppc_dens_overlay(log(mesquite$weight), y_rep_5[subset,]) +
    ggtitle ("Model 5: 4 predictors", subtitle="log(weight) ~ log(volume) + log(area) + log(shape) + group") 
  
# still decent compare to model 2, with 4 predictors versus 6.

grid.arrange(p_1, p_2, p_3, p_5, ncol=2)

######
# Discrete predictors
#######
kidiq <- as_tibble(read.dta(file="kidiq.dta"))
fit <- stan_glm (kid_score ~ as.factor(mom_work), data = kidiq, refresh = 0)
print (fit)
y1 <-    coef(fit)[1]
y2 <- y1+coef(fit)[2]
y3 <- y1+coef(fit)[3]
y4 <- y1+coef(fit)[4]

df_bar <- data.frame(x=c(1, 2, 3, 4, 5), y=c(y1, y2, y3, y4, y4))
kidiq <- kidiq %>%
  mutate (
    level <- case_when(
      mom_work == 1 ~ coef(fit)[1],
      mom_work == 2 ~ y1 + coef(fit)[2],
      mom_work == 3 ~ y2 + coef(fit)[3],
      mom_work == 4 ~ y3 + coef(fit)[4]
))
fit2 <- stan_glm (kid_score ~ mom_work, data = kidiq, refresh = 0)
print (fit2)

ggplot (kidiq) +
  geom_jitter(aes(x=mom_work, y=kid_score), width=0.1, size=0.5) +
  geom_step (data=df_bar, aes(x=x-.5, y=y), size=1) +
#  geom_segment(x=0.5, xend=1.5, y=y1, yend=y1) +
#  geom_segment(x=1.5, xend=2.5, y=y2, yend=y2) +
#  geom_segment(x=2.5, xend=3.5, y=y3, yend=y3) +
#  geom_segment(x=3.5, xend=4.5, y=y4, yend=y4) +
  geom_abline(intercept = coef(fit2)[1], slope = coef(fit2)[2], color="red") +
  labs (title = "Kid score vs mom work status - as continuous and categorical variable",
        x="Mom's work status", y="Kid's score")
  
#####
# students (Portugal)
#####
students <- as_tibble(read.csv('https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Student/data/student-merged.csv',
                               header=T)) 
colnames (students)
predictors <- c( "school",     "sex",        "age",        "address",    "famsize",    "Pstatus",   
                "Medu",       "Fedu",      
                "traveltime", "studytime",  "failures",   "schoolsup",  "famsup",     "paid",      
                "activities", "nursery",    "higher",     "internet",   "romantic",   "famrel",    
                "freetime",   "goout",      "Dalc",       "Walc",       "health",     "absences")
data_G3mat <- subset (students, G3mat > 0, select =c("G3mat", predictors))

fit0 <- stan_glm (G3mat ~ ., data = data_G3mat, refresh=0)
print (fit0)
median(bayes_R2(fit0)) # 30
median(loo_R2(fit0))   # 17
y_rep_0 <- posterior_predict (fit0)

subset <- sample (4000, 100)
p0 <- ppc_dens_overlay(data_G3mat$G3mat, y_rep_0[subset,]) +
  ggtitle ("Model0")

mcmc_areas(as.matrix(fit0), pars=vars(-'(Intercept)',-sigma),
           prob_outer=0.95, area_method = "scaled height") +
  xlim(c(-3.2,2.4))

# "Standardize" predictors
datastd_G3mat <- data_G3mat
datastd_G3mat[, predictors] <- scale(data_G3mat[, predictors])

fit1 <- stan_glm (G3mat ~ ., data = datastd_G3mat, refresh=0)
median(bayes_R2(fit1)) # 31
median(loo_R2(fit1))   # 17

y_rep_1 <- posterior_predict (fit1)
p1 <- ppc_dens_overlay(data_G3mat$G3mat, y_rep_1[subset,])  +
  ggtitle ("26 predictors, default prior") # same?
mcmc_areas(as.matrix(fit1), pars=vars(-'(Intercept)',-sigma),
           prob_outer=0.95, area_method = "scaled height") +
  xlim(c(-1.2,0.8))

# "Standardize" predictors, use prior standardized by num of predictors
fit2 <- stan_glm (G3mat ~ ., data = datastd_G3mat, refresh=0,
                   prior=normal(scale=sd(datastd_G3mat$G3mat)/sqrt(0.3*26)))
median(bayes_R2(fit2)) # 30
median(loo_R2(fit2))   # 17

y_rep_2 <- posterior_predict (fit2)
p2 <- ppc_dens_overlay(datastd_G3mat$G3mat, y_rep_2[subset,]) +
  ggtitle ("26 predictors, standardized prior") # same
mcmc_areas(as.matrix(fit2), pars=vars(-'(Intercept)',-sigma),
           prob_outer=0.95, area_method = "scaled height") +
  xlim(c(-1.2,0.8))

# "Standardize" predictors, use regularized "horseshoe" prior 
p <- length (predictors)  # 26
n <- nrow (data_G3mat)    # 357
p0 <- 6
slab_scale   <- sqrt(0.3/p0)*sd(data_G3mat$G3mat) # 0.73
global_scale <- (p0/(p-p0))/sqrt(n)               # 0.016

fit3 <- stan_glm (G3mat ~ ., data = datastd_G3mat, refresh=0,
                  prior=hs(global_scale=global_scale, slab_scale=slab_scale))
median(bayes_R2(fit3)) # 23
median(loo_R2(fit3))   # 19

y_rep_3 <- posterior_predict (fit3)
p3 <- ppc_dens_overlay(datastd_G3mat$G3mat, y_rep_3[subset,]) +
  ggtitle ("26 predictors, horseshoe prior") # same
mcmc_areas(as.matrix(fit3), pars=vars(-'(Intercept)',-sigma),
           prob_outer=0.95, area_method = "scaled height") +
  xlim(c(-1.2,0.8))
# predictors that stand out are: failures, schoolsup, goout, absences

# Use only above 4 predictors, with default weak priors!!!
fit4 <- stan_glm (G3mat ~ failures + schoolsup + goout + absences, 
                  data = datastd_G3mat, refresh=0)
median(bayes_R2(fit4)) # 20
median(loo_R2(fit4))   # 17

y_rep_4 <- posterior_predict (fit4)
p4 <- ppc_dens_overlay(datastd_G3mat$G3mat, y_rep_4[subset,]) +
  ggtitle ("4 predictors, default prior")
mcmc_areas(as.matrix(fit4), pars=vars(-'(Intercept)',-sigma),
           prob_outer=0.95, area_method = "scaled height") +
  xlim(c(-1.2,0.8))

grid.arrange(p1, p2, p3, p4, ncol=2)


# Ex 12-1
# wt_1 <- 150 + 1.8*age_10               sigma = 34.5
# wt_2 <- 108 + 21*age_10 -2*age_10_sq   sigma = 33.9

age    <- c(0, 50, 100)
age_10 <- c(0, 5,  10) 
wt_1   <- c(150, 159, 168)
wt_2   <- c(110, 163, 118)  # max wt_2 at age_10=5.25, wt_2(5.25)=165

sd_1 <- 34.5
sd_2 <- 33.9
age <- 19:90
age_10 <- age/10
age_10_sq <- age_10^2
n <- length (age)
wt_1 <- 148.7 + 1.8*age_10 + rnorm (n, 0, sd_1)
wt_2 <- 108 + 21.3*age_10 - 2*age_10_sq + rnorm (n, 0, sd_2)

108 + 21.3*age_10 - 2*age_10^2
df <- data.frame(age, age_10, age_10_sq, wt_1, wt_2)
ggplot (df) +
  geom_point (aes(x=age_10, y=wt_1), color = "red") +
  geom_abline(intercept = 148.7,        slope = 1.8, color="red") +
  geom_abline(intercept = 148.7 + sd_1, slope = 1.8, color="red", linetype = "dashed") +
  geom_abline(intercept = 148.7 - sd_1, slope = 1.8, color="red", linetype = "dashed") +
  geom_point (aes(x=age_10, y=wt_2), color = "blue") +
  geom_function(fun = function(x) 108 +       21.3*x - 2*x^2, color="blue") +
  geom_function(fun = function(x) 108 + sd_2+ 21.3*x - 2*x^2, color="blue", linetype = "dashed") +
  geom_function(fun = function(x) 108 - sd_2+ 21.3*x - 2*x^2, color="blue", linetype = "dashed") 

# Ex 12-2
earn12 <- as_tibble (read.csv('https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Earnings/data/earnings.csv',header=T))

age_divs <- c(17.9, 29, 44, 64, 91.1) 
earn12$age_factor <- cut (earn12$age, breaks = age_divs)

fit <- stan_glm (weight ~ age_factor, data = earn12, refresh = 0)
print (fit)
fit2 <- stan_glm (weight ~ age, data = earn12, refresh = 0)
print (fit2)

y1 <-    coef(fit)[1]
y2 <- y1+coef(fit)[2]
y3 <- y1+coef(fit)[3]
y4 <- y1+coef(fit)[4]

df_step <- data.frame(age_divs, levels=c(y1, y2, y3, y4, y4))

ggplot (earn12) +
  geom_point(aes(x=age, y=weight), size=0.5) +
  geom_step(data=df_step, aes(x=age_divs, y=levels), size=1) +
  geom_abline(intercept = coef(fit2)[1], slope = coef(fit2)[2], color="red", size=1) +
  labs(x="Age", y="Weight", title="Weight vs Age Categories") +
  xlim (min(age_divs), max(age_divs))


########
# Chapter 13
########
# logit(x) = log (x/(1-x)), invlogit (x) = exp(x)/(1+exp(x))

########
attach("/Users/nevinaltaras/Downloads/nes.rda")
nes92 <- nes %>% 
  filter (year == 1992) %>% 
  as_tibble () %>% 
  select (income, rvote) %>% 
  na.omit()

fit1 <- stan_glm (rvote ~ income, family = binomial(link = logit), data = nes92, refresh = 0)
print (fit1, digits=2)
sims <- as.matrix(fit1)
f <- function(x) invlogit(median(sims[,1]) + median(sims[,2])*x)
sample20 <- sample(4000, 20)

y <- array (NA, c(nrow(nes92), 20))
for (s in 1:20){
  y[,s] <- invlogit(sims[sample20[s],1]+ sims[sample20[s],2]*nes92$income)
}
y <- as_tibble (cbind(income=nes92$income,y))
df_long <- y %>% pivot_longer(-income, names_to = "sims", values_to = "prob")
ggplot(nes92) +
  geom_jitter (aes(x=income, y=rvote), width=0.04, height = 0.02, size=0.1) +
  geom_line (data=df_long, aes(x=income, y=prob, group=sims), color="gray") +
  geom_function(fun = function(x) invlogit(median(sims[,1]) + median(sims[,2])*x), size = 1) +
  labs (x="Income", y="Probability voting R") + 
  xlim (0, 6)

# predictions:  Probability in the population of voters with income=5 of voting R 
new <- data.frame(income=5)
pred <- predict (fit1, newdata = new, type = "response") # 47%
#alternatively
pred <- invlogit (predict (fit1, newdata = new)) # 47%
#
# get simulation draws for liner predictors, a + b*5  
linpred <- posterior_linpred (fit1, newdata = new)
invlogit(linpred)
# alternatively
epred <- posterior_epred(fit1, newdata = new)

median (invlogit(linpred)) # 47
median (epred)             # 47

cat ("mean epred=", round(median (epred), 2), "sd epred=", round(sd (epred), 3))
#
# get the vote of a single voter with income=5: 1=R, 0=D
postpred <- posterior_predict (fit1, newdata = new)
mean(postpred) # 47%  - E(ynew | xnew)

# instead of new <- 5, try array of new <- c(1, 2, 3, 4, 5)
new <- data.frame(income=1:5)
pred     <- predict          (fit1, newdata = new, type="response")
linpred  <- posterior_linpred(fit1, newdata = new)
epred    <- posterior_epred  (fit1, newdata = new)
postpred <- posterior_predict(fit1, newdata = new)

head (epred)
head (postpred)

mean (epred[,5] > epred[,4]) # 100% !!!!!!  

cat ("mean difference in support for Bush with income 5 vs 4", 
     round(mean(epred[,5] - epred[,4]),3),
     "sd =", 
     round(sd(epred[,5] - epred[,4]),3))
quantile (epred[,5] - epred[,4], c(0.025, 0.975)) # [3.8%, 9.3%]
  
total <- apply (postpred, 1, sum)

mean (total >= 1)  # probability 1 income level supports Bush = 88%
mean (total >= 3)  # probability 3 income levels support Bush = 23%
mean (total >= 4)  # probability 4 income levels support Bush = 5%
mean (total == 5)  # probability 4 income levels support Bush = 0.2%

# Bayesian inference

bayes_sim <- function (n, a=-2, b=0.8){
  x <- runif (n, -1, 1)
  z <- rlogis (n, a + b*x, 1)
  y <- ifelse (z > 0, 1, 0)
  fake <- data.frame(x, y, z)
  glm_fit <- glm (y ~ x, family=binomial (link="logit"), data=fake)
  stan_glm_fit <- stan_glm (y ~ x, family=binomial(link="logit"), data=fake,
                            prior=normal (0.5, 0.5), refresh=0)
  display (glm_fit, digits=1)
    
  stan_glm_fit <- stan_glm(y ~ x, family=binomial(link="logit"), data=fake, 
                           prior=normal(0.5, 0.5), refresh = 0) 
  print (stan_glm_fit, digits=1)
}

bayes_sim (n=10)   # coeffs: glm (-2.4, 2.3), stan_glm (-1.5, 0.6)
bayes_sim (n=100)  # coeffs: glm (-2.0, 0.5), stan_glm (-2.0, 0.5)
bayes_sim (n=1000) # coeffs: glm (-2.0, 1.1), stan_glm (-1.9, 1.1)
# as we increase number of observations, coef of b in stan_glm
# moves from prior of 0.5 to same as the coeff of b in glm 

#############
# Bangladesh
#############
wells <- read.table ("arsenic", header = T, sep = ",") %>% as_tibble ()
wells$label <- as.factor(ifelse (wells$switch==0, "Switch: No", "Switch: Yes"))

ggplot(wells) + 
  geom_histogram(aes(x=dist), bins = 30, color="black", fill="gray") +
  labs (x="Distance t nearest well")

fit1 <- stan_glm (switch ~ dist100, family = binomial(link = "logit"), 
                  data = wells, refresh = 0)
print (fit1, digits=3)

# switch prob = invlogit (0.61 -0.63*dist100)
# Interpret:
# if dist100 is zero, i.e. there is a safe well next to your house, the switch
# probability is invlogit(0.61) = 65%
# probability of switching at 100m is invlogit (0.61 - 0.63*1) = 50%
dist_bar <- mean (wells$dist100) # 48
# probability of switching at avg distance is invlogit (0.61 - 0.63*0.48) = 57%
# Marginal probability of switching at average dist of 48m:
# DIVIDE BY 4 RULE:  -0.63/4 = -0.16: At the average distance, 48m, the probabilty  
# of switching to a safe well at 148m decreases by 16% : 57 - 16 = 41%

ggplot (wells) +
  geom_jitter(aes(x=dist, y=switch), height=0.02, size=0.05) +
  geom_function(fun = function (x) invlogit(coef(fit1)[1] + coef(fit1)[2]*x/100), size=1) +
  labs (x="Distance to safe well", y="Switch probability",
  title = "Probability of switching vs distance to safe well") 

ggplot(wells) + 
  geom_histogram(aes(x=dist), bins = 30, color="black", fill="gray") +
  facet_wrap (~ label, ncol=2) +
  labs (x="Distance to nearest well")

# add arsenic to fit
fit2 <- stan_glm (switch ~ dist100 + arsenic, family = binomial(link = "logit"), 
                  data = wells, refresh = 0)
print (fit2, digits=2)

# switch probability = invlogit (-0.9*dist100 + 0.46*arsenic)
# interpret
# if arsenic level in current well=0 and dist to nearest safe level=0,
# switch prob = invlogit(0) = 50%
# DIVIDE BY 4 RULE:
# -0.9/4 = -0.22 : 100m father of safe well reduces switch prob by 22%
# 0.46/4 = 0.11  : increase of 1 unit in arsenic content in current well
# increases switch prob by 11%

# compare models (useless here)
subset <- sample (4000, 20) 
yrep1 <- posterior_predict(fit1)
ppc_dens_overlay (wells$switch, yrep1[subset,]) # good fit
yrep2 <- posterior_predict(fit2)
ppc_dens_overlay (wells$switch, yrep2[subset,]) # indistinguishable - useless

p1 <- ggplot (wells) +
  geom_jitter (aes(x=arsenic, y=switch), height=0.02, size=0.1) +
  geom_function(fun = function (x) invlogit(cbind (1, 0, x) %*% coef(fit2))) +
  geom_function(fun = function (x) invlogit(cbind (1, 0.5, x) %*% coef(fit2)),linetype="dashed") +
  xlim (0, ceiling(max(wells$arsenic))) +
  labs (x="Arsenic", y="Switch probability",
        title ="Model w/ no Interactions\nSwitch probability to safe well vs Arsenic level in current well",
        subtitle ="Safe well distances=0(solid) and 50m(dashed)")

p2 <- ggplot (wells) +
  geom_jitter (aes(x=dist100, y=switch), height=0.02, size=0.1) +
  geom_function(fun = function (x) invlogit(cbind (1, x, 0.5) %*% coef(fit2))) +
  geom_function(fun = function (x) invlogit(cbind (1, x, 1) %*% coef(fit2)), linetype="dashed") +
  xlim (0, ceiling(max(wells$dist100))) +
  labs (x="Distance", y="Switch probability",
        title ="Model w/ no Interactions\nSwitch probability vs Distance to safe well",
        subtitle ="Arsenic levels=0.5(solid) and 1(dashed)")

# Ex 13-1
attach("/Users/nevinaltaras/Downloads/nes.rda")

nes92_play <- nes %>% 
  filter (year == 1992) %>% 
  select (rvote, income, female, age, white, religion, educ1) %>%
  as_tibble () %>% 
  na.omit()

fit01 <- stan_glm (rvote ~ income,
#                 + factor(female) + age + factor(white) + factor(religion) + factor(educ1),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit01, digits=2)
loo01 <- loo (fit01)  # -867.8

fit011 <- stan_glm (rvote ~ factor(white),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit011, digits=2)
(loo011 <- loo (fit011))  # -857.8

fit012 <- stan_glm (rvote ~ factor(religion),
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit012, digits=2)
(loo012 <- loo (fit012))  # -855.6  Best single predictor

fit013 <- stan_glm (rvote ~ factor(white) + factor(religion),
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit013, digits=2)
(loo013 <- loo (fit013))  # -828.4 GOOD!

fit014 <- stan_glm (rvote ~ factor(female),
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit014, digits=2)
(loo014 <- loo (fit014))  # -879  BAD!!

fit015 <- stan_glm (rvote ~ factor(religion) + age,
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit015, digits=2)
(loo015 <- loo (fit015))  # -856.2 BAD!!

fit016 <- stan_glm (rvote ~ factor(white) + factor(religion) + educ1,
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit016, digits=2)
(loo016 <- loo (fit016))  # -818.9 GOOD!

fit017 <- stan_glm (rvote ~ educ1,
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit017, digits=2)
(loo017 <- loo (fit017))  # -871.6 

fit02 <- stan_glm (rvote ~ income + factor(religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit02, digits=2)
(loo02 <- loo (fit02)) # -842.0

fit03 <- stan_glm (rvote ~ income + factor(female),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit03, digits=2)
(loo03 <- loo (fit03)) # -868.7

fit04 <- stan_glm (rvote ~ income + factor(white),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit04, digits=2)
loo04 <- loo (fit04) # -850.8  GOOD!

fit05 <- stan_glm (rvote ~ income + factor(white) + age,
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit05, digits=2)
loo05 <- loo (fit05) # -851.2  

fit06 <- stan_glm (rvote ~ income + factor(white) + factor (educ1),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit06, digits=2)
(loo06 <- loo (fit06)) # -850.9  

fit07 <- stan_glm (rvote ~ income + factor(white) + factor (educ1) + factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit07, digits=2)
(loo07 <- loo (fit07)) # -817.4 GOOD!  

fit08 <- stan_glm (rvote ~ income + factor(white)  + factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit08, digits=2)
loo08 <- loo (fit08) # -820.0 

# add some interactions to fit07
fit09 <- stan_glm (rvote ~ income + factor(white) + educ1 + factor (religion) +
                   income:factor(white) + factor(white):factor (religion), 
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit09, digits=2)
(loo09 <- loo (fit09)) # -810.5 GOOD! 

fit091 <- stan_glm (rvote ~ income + factor(white) + educ1 + factor (religion) +
                    factor(white):factor (religion), 
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit091, digits=2)
(loo091 <- loo (fit091)) # -810.5 GOOD! 

fit092 <- stan_glm (rvote ~ factor(white) + educ1 + factor (religion) +
                      factor(white):factor (religion), 
                    family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit092, digits=2)
(loo092 <- loo (fit092)) # -812.9 Income is not a great predictor 

fit10 <- stan_glm (rvote ~ income + factor(white) + factor (educ1) + factor (religion) +
                   income:factor(white) + factor(white):factor (religion) +
                   income:factor(educ1) + income:factor(religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit10, digits=2)
loo10 <- loo (fit10) # -814.5 

# take education out
fit11 <- stan_glm (rvote ~ income + factor(white) + factor (religion) +
                   income:factor(white) + factor(white):factor (religion) +
                   income:factor(religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit11, digits=2)
(loo11 <- loo (fit11, k_threshold = 0.7)) # -814.8  

fit12 <- stan_glm (rvote ~ income + factor(white) + factor (religion) +
                     income:factor(white) + factor(white):factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit12, digits=2)
(loo12 <- loo (fit12)) # -814.0  GOOD!  

fit13 <- stan_glm (rvote ~ income + factor(white) + factor (religion) +
                     income:factor(white),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit13, digits=2)
(loo13 <- loo (fit13)) # 819

fit14 <- stan_glm (rvote ~ income + factor(white) + factor (religion) +
                     factor(white):factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit14, digits=2)
(loo14 <- loo (fit14)) # -814.1

fit15 <- stan_glm (rvote ~ income + factor(white) + factor (religion) + factor(educ1) +
                     factor(white):factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit15, digits=2)
(loo15 <- loo (fit15)) # -811.1

fit16 <- stan_glm (rvote ~ factor(income) + factor(white) + factor (religion) + factor(educ1) +
                     factor(white):factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit16, digits=2)
(loo16 <- loo (fit16)) # -813.9

fit17 <- stan_glm (rvote ~ income + factor(white) + factor(religion) + educ1 +
                     factor(white):factor (religion),
                   family = binomial(link="logit"), data = nes92_play, refresh = 0)
print (fit17, digits=2)
(loo17 <- loo (fit17)) # 810.5

# Worst predictors: female, age, income
# Best predictors:  religion, race and educ1 

# Ex 13-2 
# a) regular invlogit
# b) shifted to left by 2
# c) compressed in x by a factor of 2
# d) compressed in x by a factor of 2 and then shifted to left by 2
# e) compressed in x by a factor of 2 and transposed around x=0

# Ex 13-3

n <- 4000
a <- 46.3
b <- 3.1
sigma <- 3.8
gdp <- seq(-2, 4.5, length.out = 100)

win_prob <- rep (NA, length(gdp))
for (i in 1:length(gdp)){
  win_prob[i] <- mean (a + b*gdp[i] + rnorm (n, 0, sigma) > 50)
}

df <- cbind (gdp, win_prob, logit_win_prob = logit (win_prob)) %>% data.frame ()
df_long <- df %>%  pivot_longer (-gdp, names_to = "prob_type", values_to = "probability")
df_long$prob_type <- factor(df_long$prob_type,levels=unique(df_long$prob_type))
ggplot (df_long) +
  geom_point(aes(x=gdp, y=probability)) + 
  facet_wrap (~ prob_type, scales = "free") +
  labs (x="GDP", title="Pr(win), logit (Pr(win)")
  
fit <- stan_glm (logit_win_prob ~ gdp,  data = df, refresh = 0)
print (fit)

label_fit <- paste0("invlogit(", round (coef(fit)[1], 1), " + ", round (coef(fit)[2], 1), "*GDP)")

df <- cbind (df, fitted = invlogit (coef(fit)[1] + coef(fit)[2]*gdp))
ggplot (df) +
  geom_point(aes(x=gdp, y=win_prob)) +
  geom_line (aes(x=gdp, y=fitted), color='red', size=1) +
  annotate("text", x = 0.50, y = 0.75, label = label_fit, color = 'red') +
  labs (x="GDP", y="Pr(Vote>50%)",
  title="Pr (Incumbent party wins) vs GDP, Invlogit fit")
  
cat ("For a 1% higher GDP, expect:", round (coef(fit)[2]/4*100, 0), 
     "% change in probability of winning")

########
# Ch 14
########
# model: Pr(y=1) = invlogit (a + b*x)
n <- 50
a <- 2
b <- 3
x_mean <- -a/b
x_sd   <- b/4
x <- rnorm (n, x_mean, x_sd)
y <- rbinom (n, 1, invlogit (a+b*x))
df <- data.frame(x,y)
fit1 <- stan_glm (y ~ x, family = binomial(link = "logit"), data = df, refresh = 0)
print (fit1)
ggplot (df) + 
  geom_point(aes(x=x, y=y)) +
  geom_function(fun = function (x) invlogit (a + b*x), color = 'red') +
  annotate("text", x = -0.6, y = 0.15, label = "True Curve:\ninvlogit(2 + 3*x)", color = 'red') +
  geom_function(fun = function (x) invlogit (cbind (1, x) %*% coef(fit1)), , color = 'blue') +
  annotate("text", x = -1.5, y = 0.65, label = "Fitted Curve:\ninvlogit(1.9 + 2.3*x)", color = 'blue') +
  xlim (floor(min(x)), ceiling(max(x))) +
  labs (x="x", y="y", title="Data, Regenerative model, Fitted Curve")

# bin data
n_bins = 10
bins <- as.numeric (cut (x, n_bins))
x_bar <- rep (NA, n_bins)
y_bar <- rep (NA, n_bins)
for (k in 1:n_bins){
  x_bar[k] <- mean(x[bins==k])
  y_bar[k] <- mean(y[bins==k])
}

df_bar <- data.frame(x_bar, y_bar)
ggplot () + 
  geom_point(data=df,     aes(x=x,     y=y)) +
  geom_point(data=df_bar, aes(x=x_bar, y=y_bar), color="red", size=3) +
  labs (title="Binning: data aggregated to 5 bins")

# model: Pr(y=1) = invlogit (a + b1*x1 + b2*x2)
n <- 100
beta <- c(2, 3, 4)
x1 <- rnorm (n, 0, 0.4)
x2 <- rnorm (n, -0.5, 0.4)
y  <- as.factor(rbinom (n, 1, invlogit (cbind(1, x1, x2) %*% beta)))
df2 <- data.frame(x1, x2, y)

fit2 <- stan_glm (y ~ x1 + x2, family = binomial(link = "logit"), data = df2, refresh = 0)
print (fit2)
b0 <- coef(fit2)[1]
b1 <- coef(fit2)[2]
b2 <- coef(fit2)[3]
# Pr(y=1) = 50% line when invlogit (b0 + b1*x1 + b2*x2) = 0.5 =>
#                        b0 + b1*x1 + b2*x2 = logit (0.5) = 0
# x1 = - (b0 + b1*x2)/b2       = -b0/b2 - b1/b2 * x2
# Pr(y=1) = 90% line when b0 + b1*x1 + b2*x2 = logit (0.9) = 2.2
# x2 = (2.2 - (b0 + b1*x2))/b2 = (2.2 - b0)/b2 - b1/b2 * x2
# Pr(y=1) = 10% line when b0 + b1*x1 + b2*x2 = logit (0.1) = -2.2
# x2 = (-2.2 - (b0 + b1*x2))/b2 = (-2.2 - b0)/b2 - b1/b2 * x2

ggplot (df2) +
  geom_point (aes(x=x1, y=x2, color=y)) +
  geom_abline(intercept = (2.2-b0)/b2, slope = - b1/b2, linetype="dashed") +
  geom_abline(intercept = -b0/b2, slope = - b1/b2) +
  geom_abline(intercept = (-2.2-b0)/b2, slope = - b1/b2, linetype="dashed") +
  labs (title = "Data and Pr(y=1)=10%, 50%, 90% (discrimination) lines",
        "\nfrom fitted logistic regression")
########
# back to arsenic
########
# interaction after centering
wells <- wells %>%
  mutate (dist100_c = dist100 - mean (dist100),
          arsenic_c = arsenic - mean (arsenic))

fit5 <- stan_glm (switch ~ dist100_c + arsenic_c + dist100_c:arsenic_c, family = binomial(link="logit"), 
                  data = wells, refresh = 0)
print (fit5, digits=2)
(loo5 <- loo (fit5)) -1968

# Median MAD_SD
# (Intercept)          0.35   0.04 
# dist100_c           -0.88   0.10 
# arsenic_c            0.47   0.04 
# dist100_c:arsenic_c -0.18   0.10 
#
# Expected prob of switching when dist to safe well is 48m and the arsenic level in current well is 1.66
# is invlogit (coef(fit5)[1]) = invlogit (0.35) = 59%
# Expected prob of switching when arsenic level in current well is 1.66 decreases by 22% if the safe well 
# is 100m further
# Expected prob of switching when dist to safe well is 48m, increases by 12% if the current well arsenic 
# increases by 1
# 
# interaction: 
# a) For each additional unit of arsenic, the coefficient of distance100 changes by -0.18.
# at average arsenic level of 1.66, the coefficient of distance is -88.  So, this coefficient 
# becomes -0.88 - 0.18 = -1.06
# b) For each additional 100m to next safe well, the coefficient of arsenic changes by -0.18.
# at average distance of 48m, the coefficient of distance is 0.47.  So, this coefficient 
# becomes 47 - 0.18 = 0.29.  

p3 <- ggplot (wells) +
  geom_jitter (aes(x=arsenic, y=switch), height=0.02, size=0.1) +
  geom_function(fun = function (x) invlogit(cbind (1, 0,   x, 0)     %*% coef(fit5))) +
  geom_function(fun = function (x) invlogit(cbind (1, 0.5, x, 0.5*x) %*% coef(fit5)),linetype="dashed") +
  xlim (0, ceiling(max(wells$arsenic))) +
  labs (x="Arsenic", y="Switch probability",
        title ="Model w/ Interactions\nSwitch probability to safe well vs Arsenic level in current well",
        subtitle ="Safe well distances=0(solid) and 50m(dashed)")

p4 <- ggplot (wells) +
  geom_jitter (aes(x=dist100, y=switch), height=0.02, size=0.1) +
  geom_function(fun = function (x) invlogit(cbind (1, x, 0.5, 0.5*x) %*% coef(fit5))) +
  geom_function(fun = function (x) invlogit(cbind (1, x, 1,       x) %*% coef(fit5)), linetype="dashed") +
  xlim (0, ceiling(max(wells$dist100))) +
  labs (x="Distance", y="Switch probability",
        title ="Model w/ Interactions\nSwitch probability vs Distance to safe well",
        subtitle ="Arsenic levels=0.5(solid) and 1(dashed)")
grid.arrange(p1, p2, p3, p4, ncol=2)

# add predictors
wells <- wells %>%
  mutate (educ4 = educ/4,
          educ4_c =educ4 - mean (educ4))

fit6 <- stan_glm (switch ~ dist100 + arsenic + educ4 + assoc, 
                  family = binomial(link = "logit"), data = wells, refresh = 0)
print (fit6, digits=2)
# coeff of assoc is negative (weird) and has high se.  Drop it
# center, add interactions
fit7 <- stan_glm (switch ~ dist100_c + arsenic_c + educ4_c +
                    dist100_c:educ4_c +
                    arsenic_c:educ4_c, 
                  family = binomial(link = "logit"), data = wells, refresh = 0)
print (fit7, digits=2)
(loo6 <- loo (fit6)) # -1959
(loo7 <- loo (fit7)) # -1953
(loo_compare(loo7, loo5))
#loo7 elpd is higher by 15.4

# Predictive simulation
# simple version:
fit <- stan_glm (switch ~ dist100, family = binomial(link = "logit"), data = wells, refresh = 0)
print (fit, digits=2)
sims <- as.matrix (fit)
n_sims <- nrow (sims)

df <- data.frame (beta0 = sims[,1], beta1 = sims[,2])
ggplot (df) + geom_point(aes(x=beta0, y=beta1), size = 0.5)

sample20 <- sample (4000, 20)
y <- array (NA, c(nrow(wells), 20))
for (s in 1:20){
  y[,s] <- invlogit(sims[sample20[s],1] + sims[sample20[s],2]*wells$dist100)
}
y <- as_tibble (cbind(dist100=wells$dist100,y))
df_long <- y %>% pivot_longer(-dist100, names_to = "sims", values_to = "prob")

ggplot (wells) +
  geom_jitter (aes(x=dist100*100, y=switch), height = 0.02, size = 0.15) +
  geom_line (data=df_long, aes (x=dist100*100, y=prob, group=sims), color="gray") +
  geom_function(fun = function (x) invlogit (coef(fit)[1] + coef(fit)[2]*x/100), size=1) +
  labs (x="Distance to nearest safe well (m)", y="Pr(Switch)", 
        title="Switch probability: Data, logistical fit, 20 simulated draws",
        subtitle="logit (0.60 - 0.62 * dist100)")

# Predictive simulation with binomial distribution.  n_new households
n_new <- 100
X_new <- cbind (1, sort(as.numeric(sample (300, n_new)/100)))
y_new <- array (NA, c(n_sims, n_new))
for (s in 1:n_sims){
  p_new <- invlogit (X_new %*% sims[s,])
  y_new[s, ] <- rbinom (n_new, 1, p_new)
}

expected_switch_prob <- colMeans(y_new)

df <- data.frame(cbind (X_new, expected_switch_prob))

ggplot (df) +
  geom_point (aes(x=V2, expected_switch_prob)) +
  geom_function(fun = function (x) invlogit (coef(fit)[1] + coef(fit)[2]*x), size=1) +
  ylim (0,1) +
  labs (x="Distance from safe wells", y="Expected Pr(Switch)", 
        title = "100 simulations",
        subtitle="Expected Pr(Switch) of simulated homes, Model fit to actual homes")

# Residuals
predicted <- invlogit (predict (fit7))
error <- wells$switch - predicted
df <- data.frame (arsenic = wells$arsenic, dist = wells$dist, predicted, error)
ggplot (df) + geom_point(aes (x=predicted, y=error)) 
# bin errors
binnedplot(predicted, error, nclass = 40,
           main="Residuals: Estimated Pr (switching)")

binnedplot(wells$dist, error, nclass = 40,
           main="Residuals: Distance from nearest safe well")

binnedplot(wells$arsenic, error, nclass = 40,
           main="Residuals: Arsenic at current well")

# Redo fit7 with taking log of arsenic
wells$log_arsenic <- log (wells$arsenic)
wells$log_arsenic_c <- log (wells$arsenic) - mean (log (wells$arsenic))
fit7p <- stan_glm (switch ~ dist100_c + log_arsenic_c + educ4_c +
                    dist100_c:educ4_c +
                    log_arsenic_c:educ4_c, 
                  family = binomial(link = "logit"), data = wells, refresh = 0)
print (fit7p, digits=2)
(loo7p <- loo (fit7p)) # -1938  20 better!!

predicted_p <- invlogit (predict (fit7p))
error_p <- wells$switch - predicted_p

binnedplot(wells$arsenic, error_p, nclass = 40,
           main="Residuals: Arsenic at current well")
###
nes_play <- nes %>% 
  select (year, rvote, income, female, black) %>%
  as_tibble () %>% 
  na.omit()

years <- unique (nes_play$year)
n_years <- length (years)
fit_coefs       <- array (NA, c(n_years, 4), 
                          dimnames = list (years, c("rvote", "income", "female", "black")))
fit_ses         <- array (NA, c(n_years, 4), 
                          dimnames = list (years, c("servote", "seincome", "sefemale", "seblack")))
fit_bayes_coefs <- array (NA, c(n_years, 4), 
                          dimnames = list (years, c("rvote", "income", "female", "black")))
fit_bayes_ses   <- array (NA, c(n_years, 4), 
                          dimnames = list (years, c("servote", "seincome", "sefemale", "seblack")))

for (i in 1:n_years){
  cat ("\n", years[i])
  fit       <-      glm (rvote ~ income + female + black, data = nes_play,
                             family = binomial(link = "logit"), subset = (year == years[i]))
  fit_coefs[i,] <- coef(fit)
  fit_ses[i,]   <- se.coef(fit)
  
  fit_bayes <- stan_glm (rvote ~ income + female + black,data = nes_play,
                             family = binomial(link = "logit"), subset = (year == years[i]),
                             refresh = 0 )
  fit_bayes_coefs[i,] <- coef(fit_bayes)
  fit_bayes_ses[i,]   <- se(fit_bayes)
}

df <- cbind (years=years, fit_coefs, fit_ses) %>% data.frame() 

df_bayes <-  cbind (years=years, fit_bayes_coefs, fit_bayes_ses) %>% data.frame()

df_long <-       df %>%  pivot_longer (-years, names_to = "predictor", values_to = "vals")
df_bayes_long <- df_bayes %>%  pivot_longer (-years, names_to = "predictor", values_to = "vals")

p1 <- ggplot (df) +
  geom_point (aes(x=years, y=black)) + 
  geom_segment(aes(x=years, xend=years, y=black-seblack, yend=black+seblack))


p2 <- ggplot (df_bayes) +
  geom_point (aes(x=years, y=black)) +
  geom_segment(aes(x=years, xend=years, y=black-seblack, yend=black+seblack))

grid.arrange(p1, p2, ncol = 2)


# Ex 14-1  
n <- 50
A <- -10
B <- 15
x = runif(n, A, B) 
prob <- invlogit(0.4 - 0.3*x)
y <- rbinom (n, 1, prob)
df <- cbind (x, y, prob) %>% data.frame()

ggplot (df) +
  geom_jitter (aes(x=x, y=y), height = 0.01)

fit <- stan_glm(y ~ x, family = binomial(link="logit"), data = df, refresh = 0)
print (fit)

ggplot (df) +
  geom_jitter (aes(x=x, y=y), height = 0.01) +
  geom_function(fun = function (x) invlogit (coef(fit)[1] + coef(fit)[2]*x)) +
  geom_function(fun = function (x) invlogit (0.4 - 0.3*x), color="red") +
  labs (x="x", y="rbinom (invlogit(0.4 - 0.3*x))")
  
# Ex 14-2
n <- 50
A1 <- -7.5
B1 <- 7.5
A2 <- -12.5
B2 <- 15
x1 <- runif (n, A1, B1)
x2 <- runif (n, A2, B2)
prob <- invlogit (0.4 - 0.3*x1 - 0.2*x2)
beta <- c(0.4, -0.3, -0.2)
y <- rbinom (n, 1, invlogit (cbind (1, x1, x2) %*% beta))
sum (y)

df <- cbind (x1, x2, y) %>% data.frame()
fit <- stan_glm (y ~ x1 + x2, data = df, family = binomial(link = "logit"), refresh = 0)
print (fit)
beta_hat <- coef (fit)

# Pr(y=1)=0.5 line: invlogit (beta_hat[1] + beta_hat[2]*x1 + beta_hat[3]*x2) = 0.5
# beta_hat[1] + beta_hat[2]*x1 + beta_hat[3]*x2) = 0
# x2 = - beta_hat[1]/beta_hat[3] - beta_hat[2]/beta_hat[3]*x1

# Pr(y=1)=0.1 line: invlogit (beta_hat[1] + beta_hat[2]*x1 + beta_hat[3]*x2) = 0.1
# beta_hat[1] + beta_hat[2]*x1 + beta_hat[3]*x2) = -2.2
# x2 = -(2.2 + beta_hat[1])/beta_hat[3] - beta_hat[2]/beta_hat[3]*x1

ggplot (df) +
  geom_point (aes(x=x1, y=x2, color=factor(y))) +
  geom_abline(intercept = (logit(0.5) - beta_hat[1])/beta_hat[3], slope = -beta_hat[2]/beta_hat[3]) +
  geom_abline(intercept = (logit(0.9) - beta_hat[1])/beta_hat[3], slope = -beta_hat[2]/beta_hat[3], linetype="dashed") +
  geom_abline(intercept = (logit(0.1) - beta_hat[1])/beta_hat[3], slope = -beta_hat[2]/beta_hat[3], linetype="dashed") +
  xlim (A1, B1) + ylim (A2, B2) +
  labs (x="x1", y="x2", title = "Data and 10%, 50%, and 90% dicrimination lines\nfrom fittrd logistic regression")

# Ex 14-6
n <- 20
x <- 1:n
y <- sample (c(0, 1), n, replace = TRUE)
a <- -20
b <-  2
yp <- rbinom (n, 1, invlogit (a + b*x))
df <- data.frame (x, y, yp)

fit <- stan_glm (y ~ x, data = df, family = binomial(link = "logit"), refresh = 0)
print (fit)
fitp <- stan_glm (yp ~ x, data = df, family = binomial(link = "logit"), refresh = 0)
print (fitp)

ggplot (df) +
  geom_jitter (aes(x=x, y=y), color="black", height = 0.02) +
  geom_jitter (aes(x=x, y=yp), color="red" , height = 0.02) +
  geom_function(fun = function (x) invlogit (coef(fit)[1] + coef(fit)[2]*x)) +
  geom_function(fun = function (x) invlogit (coef(fitp)[1] + coef(fitp)[2]*x), color="red") +
  annotate("text", x = 6,  y = 0.75,    label = "y=rbinom(invlogit(a + b*x))", color="red") +
  annotate("text", x = 13, y = 0.50,    label = "y=sample(0,1)") 

# Test
ggplot () +
  geom_function(fun = function (x) invlogit (40*(9.5 - x))) +
  xlim (1, 20)
  
#######
# Ch 16
#######
# Design calculations
#
p <- 0.6 # proportion who support death penalty
# the se for estimating p_dp in a sample of n people is
# se = sqrt (p_dp*(1-p_dp)/n) = sqrt (0.6*.0.4/n) = 0.49/sqrt(n)
# se = 0.5/sqrt(n)
# if need se=5%
se= 0.05
n <- (sqrt(p*(1-p))/se)^2  #  96
# or approximately:
n <- (0.5/se)^2            # 100

# POWER
# A - demonstrate: more than 0.5 of population supports death penalty
# we "know" true proportion is p=0.6
# what sample size is needed
se <- 0.5/sqrt(n)
# p_bar: estimate for true p
# 95% CI for p = [p_bar +/- 1.96 * 0.5/sqrt(n)]
# we have demonstrated A, if whole interval is above 0.5:
# p_bar - 1.96*0.5/sqrt(n) > 0.5 =>
# p_bar > 0.5 + 1.96*0.5/sqrt(n)

# Choose n such that 80% of the 95% intervals will not include p_0=0.5
0.5 + 1.96*se = 0.6 - qnorm (0.8)*se
0.5 + 1.96*se = 0.6 - 0.84*se
# se = (p-p_0)/2.8
se <- 0.1/2.8 # 0.0357

p   <- 0.6
# sqrt(p*(1-p)/n) = (p-p_0)/2.8
n <- p*(1-p)*(2.8/(p-p_0))^2 #or
n <- (0.5*2.8/(p-p_0))^2  # 196

n_sim <- 100
sims <- rnorm (n_sim, 0.6, 0.0357142)
# percentage of 2*se ranges that are greater than 0.5
mean(sims - 2*se > 0.5) # 79.8

x_lo <- sims - 2*se
x_hi <- sims + 2*se
y_lo <- seq (0, 10, length.out = 100)
y_hi <- y_lo

df <- data.frame (x_lo, x_hi, y_lo, y_hi)

ggplot (df) +
  geom_function(fun = dnorm, args = list(mean = 0.6, sd = 0.03571)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_segment(aes(x=x_lo, xend=x_hi, y=y_lo, yend=y_hi), color='red') +
  xlim (0.4, 0.8) 

# In summary, to have 80% power, the true value of the parameter must be 
# 2.8 standard errors  away from the comparison point:
# the value 2.8 is 1.96 from the 95% interval, 
# plus 0.84 to reach the  80th percentile of the normal distribution.

# Design with fake data simulation
n <- 100
y_if_control <- rnorm (n, 60, 20) # post-test control
y_if_treated <- y_if_control + 5  # post_test treated
  
z <- sample (rep(c(0, 1), n/2))

y <- ifelse (z==1, y_if_treated, y_if_control)
df <- data.frame(y, z)

diff <- mean (y[z==1]) - mean (y[z==0])  # -1.2, 5.97, 3.6
# results all over the place: True value = 5
se_diff <- sqrt (sd(y[z==0])^2/sum(z==0) + sd(y[z==1])^2/sum(z==1))  # 3.5, 4.0, 3.9
# relatively stable
fit_1a <- stan_glm (y ~ z, data = df, refresh = 0)
print (fit_1a) # z: 5.2+/-3.8, 1.8+/-4.4, 10.2+/-3.9
# slope of z unstable, se ok

df$x <- rnorm (n, 60, 20)  # pre-test score
fit_1b <- stan_glm (y ~ z + x, data = df, refresh = 0)
print (fit_1b)

# new approach
n <- 100
true_ability <- rnorm (n, 50, 16)
x <- true_ability + rnorm (n, 0, 12)                  # pre_test
y_if_control <- true_ability + rnorm (n, 0, 12)  + 10 # post_test_control
y_if_treated <-  y_if_control + 5                     # post_test_treated
z <- sample(rep(c(0,1), n/2)) # replace=FALSE, ergo as many 1's (0's) as in the vector
                              # vector to sample from!!!
y <- ifelse (z==1, y_if_treated, y_if_control)
df2 <- data.frame (x, y, z)
fit_2a <- stan_glm (y ~ z, data = df2, refresh = 0)
print (fit_2a)
fit_2b <- stan_glm (y ~ z + x, data = df2, refresh = 0)
print (fit_2b)

# Simulate selection BIAS
z <- rbinom (n, 1, invlogit (-(x-50)/20))
df <- data.fram (x, z)
ggplot (df) +
  geom_point(aes(x=x, y=z)) +
  geom_function(fun = function (x) invlogit (-(x-50)/20)) +
  annotate("text", x = 40 , y = 0.95,    label = "Assigned to control group, z=0") +
  annotate("text", x = 40 , y = 0.05,    label = "Assigned to treatment group, z=1") +
  labs (x="Pretest score", y="Pr(z=1)")
y <- ifelse (z==1, y_if_control, y_if_control)

df <- data.frame (x, y, z)
fit_3a <- stan_glm (y ~ z, data = df, refresh = 0)
print (fit_3a)
fit_3b <- stan_glm (y ~ z + x, data = df, refresh = 0)
print (fit_3b)
rbind (c(coef(fit_3a)[2], se(fit_3a)[2]), c(coef(fit_3b)[2], se(fit_3b)[2]))

experiment <- function (n) {
  true_ability <- rnorm (n, 50, 16)
  x <- true_ability + rnorm (n, 0, 12)
  y_if_control <- true_ability + rnorm (n, 0, 12) + 10
  y_if_treated <- y_if_control + 5
  z <- rbinom (n, 1, invlogit (-(x-50)/20))   # bias
#  z <- sample(rep(c(0,1), n/2))              # no bias
  y <- ifelse (z==1, y_if_treated, y_if_control)
  df <- data.frame (x, y, z)
  fit_3a <- stan_glm (y ~ z, data = df, refresh = 0)
  fit_3b <- stan_glm (y ~ z + x, data = df, refresh = 0)
  rbind (c(coef(fit_3a)[2], se(fit_3a)[2]), c(coef(fit_3b)[2], se(fit_3b)[2]))
}

n < - 100
n_loop <- 50
results <- array (NA, c(n_loop, 2, 2), 
                  dimnames = list (1:n_loop, c("Simple", "Adjusted"), c("Estimate", "se")))
for (loop in 1:n_loop){
  print (loop)
  results[loop,,] <- experiment(n)
}

results_avg <- apply (results, c(2, 3), mean)
results_avg

#########
# Ch 19
########

electric_plots <- function (df_dat, identifier) {
  n_grades <- length (unique(df_dat$grade))
  results <- array (NA, c(n_grades, 2, 2, 3), 
                    dimnames = list (1:n_grades, 
                                     c("Simple", "Adjusted"), 
                                     c("Estimate", "se"), 
                                     c("intercept", "treatment", "pretreat")))
  for (k in 1:n_grades){
    fit_1 <- stan_glm (posttreat ~ treatment, 
                       data = df_dat, subset = (grade==k), refresh = 0)
    fit_2 <- stan_glm (posttreat ~ treatment + pretreat, 
                       data = df_dat, subset = (grade==k), refresh = 0)
    results[k,,,1] <- rbind (c(coef(fit_1)[1], se(fit_1)[1]), c(coef(fit_2)[1], se(fit_2)[1]))
    results[k,,,2] <- rbind (c(coef(fit_1)[2], se(fit_1)[2]), c(coef(fit_2)[2], se(fit_2)[2]))
    results[k,,,3] <- rbind (c(0,              0           ), c(coef(fit_2)[3], se(fit_2)[3]))
  }
  
  df <- data.frame (grade=(1:4), 
                    a_simple=results[, 1, 1, 1], a_adjusted=results[, 2, 1, 1],
                    b_simple=results[, 1, 1, 2], b_adjusted=results[, 2, 1, 2],
                    c_simple=results[, 1, 1, 3], c_adjusted=results[, 2, 1, 3],
                    a_simple_se=results[, 1, 2, 1], a_adjusted_se=results[, 2, 2, 1],
                    b_simple_se=results[, 1, 2, 2], b_adjusted_se=results[, 2, 2, 2],
                    c_simple_se=results[, 1, 2, 3], c_adjusted_se=results[, 2, 2, 3])
  
  p1 <- ggplot (df) +
    geom_pointrange(aes(x=b_simple, y=grade+0.1, xmin=b_simple-2*b_simple_se, xmax=b_simple+2*b_simple_se,
                    color='red'), size=0.5) +
    geom_linerange(aes(x=b_simple, y=grade+0.1, xmin=b_simple-0.67*b_simple_se, xmax=b_simple+0.67*b_simple_se,
                    color='red'), size=1) +
    geom_pointrange (aes(x=b_adjusted, y=grade-0.1, xmin=b_adjusted-2*b_adjusted_se, xmax=b_adjusted+2*b_adjusted_se, 
                     color="blue"), size=0.5) +
    geom_linerange (aes (x=b_adjusted, y=grade-0.1, xmin=b_adjusted-0.67*b_adjusted_se, xmax=b_adjusted+0.67*b_adjusted_se, 
                    color="blue"), size=1) +
    geom_vline(xintercept = 0, linetype= "dashed") +
    labs (title=identifier,
          subtitle="Regression with and without adjusting for pre-test as predictor - 50% and 95% intervals",
          x="Coefficient of treatment indicator", y="Grades") +
    scale_color_identity(name = "Label",
                         breaks = c("red", "blue"),
                         labels = c("No adjustment for pre-test", "Adjusted for pre-test"),
                         guide = "legend") 
    
  p2 <- ggplot () +
    geom_point(data=df_dat, aes(x=pretreat, y=posttreat, color=factor(treatment))) +
    geom_abline(data=df, aes(intercept=a_adjusted+b_adjusted, slope=c_adjusted), color="blue") +
    geom_abline(data=df, aes(intercept=a_adjusted, slope=c_adjusted), color="red") +
    facet_wrap (~ as.numeric(grade), scales="free", ncol=4) +
    labs (title=identifier,
          subtitle="Regression without interaction", 
          x="pre-treat, xi", y="post_treat, yi")+
    xlim (0, 130) + ylim (0, 130)

  plots <- list (p1, p2)
  return (plots)
}

electric<- read.table("electric", header=TRUE) %>% as_tibble() 
  
electric1 <- electric %>%
  select (Grade, treated.Posttest, control.Posttest) 

electric_long <- pivot_longer(electric1, !Grade) 

avgs <- electric_long %>%
  group_by(Grade, name) %>%
  summarise(avg = mean(value))  %>% data.frame()
sds <- electric_long %>%
  group_by(Grade, name) %>%
  summarise(sd = sd(value))  %>% data.frame()
sds$label <- paste0("avg=", round(avgs$avg, 0), "\nsd=", round(sds$sd, 0))

ggplot(electric_long) +
  geom_histogram(aes(x=value), bins = 15, color="black", fill="gray") +
#  facet_grid(Grade ~ name)
  facet_grid (name ~ Grade) +
  labs (x="Grades") +
  geom_vline (data=avgs, aes(xintercept = avg)) + #, inherit.aes=FALSE) +
  geom_text (data=sds,  aes(label=label), x = 40, y = 7, hjust=0, vjust=0) 

x <- c(electric$treated.Pretest,  electric$control.Pretest)
y <- c(electric$treated.Posttest, electric$control.Posttest)
z <- c(rep(1, dim(electric)[1]), rep(0, dim(electric)[1]))
grade <- rep (electric$Grade, 2)

df_dat <- data.frame (pretreat=x, posttreat=y, treatment=z, grade)

plots <- electric_plots(df_dat, "Estimated effect of Treatment vs Control")
plots[[1]]
plots[[2]]

n_grades <- 4
results <- array (NA, c(n_grades, 4), 
                  dimnames = list (1:n_grades, c("intercept", "treatment", "pretreat",
                                                 "interaction")))
for (k in 1:n_grades){
  fit_3 <- stan_glm (posttreat ~ treatment + pretreat + treatment:pretreat, 
                     data = df_dat, subset = (grade==k), refresh = 0)
  results[k,] <- coef(fit_3)
}

results <- data.frame (cbind(grade=1:4, results))
ggplot (data=df_dat) +
  geom_point(data=df_dat, aes(x=pretreat, y=posttreat, color=factor(treatment))) +
  geom_abline(data=results, aes(intercept=intercept+treatment, slope=pretreat+interaction), color="blue") +
  geom_abline(data=results, aes(intercept=intercept, slope=pretreat), color="red") +
  facet_wrap (~ as.numeric(grade), ncol=4) +
  labs (title="Regression with interaction", x="pre-treat, xi", y="post_treat, yi") +
  xlim (0, 130) + ylim (0, 130)

# plot uncertainty for Grade 4
fit_4 <- stan_glm (posttreat ~ treatment + pretreat + treatment:pretreat, 
                   data = df_dat, subset = (grade==4), refresh = 0)
print (fit_4)
coef(fit_4)
sims <- as.matrix (fit_4)

select_20 <- sample (1:4000, 20)

ggplot (subset(df_dat, grade==4), x=pretreat) +
  geom_abline (intercept = sims[select_20, 2], slope=sims[select_20, 4], color="darkgray", size=0.4) +
  geom_abline (intercept = coef(fit_4)[2], slope=coef(fit_4)[4]) +
  geom_hline (yintercept = 0, linetype = "dashed") +
  xlim (80, 120) + ylim (-5, 10) +
  labs (title="Grade 4", x="pre-test", y="treatment effect")

# Mean treatment for grade 4
n_sims <- nrow (sims)
effect <- array (NA, c(n_sims, sum (grade==4)))
for (s in 1:n_sims){
  effect[s,] <- sims[s,2] + sims[s,4]*df_dat$pretreat[grade==4]
}
avg_effect <- (rowMeans(effect))
mean (avg_effect) # 1.8
sd (avg_effect)   # 0.70  SAME AS MODEL WITH NO INTERCTIONS!!!
  
# Ex 19-4
n <- 100
b <- 0.7
theta <- 10
# a)

pre_test <- rnorm (n, 40, 15)
z <- sample (rep (c(0,1), n/2), n)
error <- rnorm (n, 0, 10)

post_test_control_err <- function(a){
  post_test <- a + b*pre_test + theta*z[z==0] + error
  err <- mean (post_test) - 50
}
a <- round (uniroot(post_test_control_err, interval=c(0,30), tol= 0.000001)$root)

# for control, z=0
# post_test[z==0] = a + b*pre_test + error
# mean(post_test[z==0]) <- a + b*mean(pre_test)
 50 = a + b * 40
a <- 50 - 40*b # 22

# b)
a <- 22
b <- 0.7
pre_test <- rnorm (n, 40, 15)
z <- sample (rep (c(0,1), n/2), n)
error <- rnorm (n, 0, 10)
post_test         <- a + b*pre_test + theta*z + error
post_test_control <- a + b*pre_test + error
post_test_treated <- a + b*pre_test + theta + error

sd (post_test_control)   # 13.8, 14.5, 17
# or
sqrt ( (b*15)^2 + 10^2)   # 14.5

# c)
sd (post_test_treated)   # 14.5
mean (post_test_treated) # 60

# Ex 19-5

y <- invlogit (coef(fit)[1] + coef(fit)[2]*z + coef(fit)[3]*age + coef(fit)[4]*z*age)
y_z0 <- invlogit (coef(fit)[1] + coef(fit)[3]*age)
y_z1 <- invlogit (coef(fit)[1] + coef(fit)[2] + (coef(fit)[3] + coef(fit)[4])*age)

avg_treat_effect <- (sum (y_z1*n_pop) - sum (y_z0*n_pop))/sum (n_pop)

# Ex 19-8
cows <- read.table("cows", header=TRUE) %>% 
  rename (initial_wt = initial.weight,
          final_wt   = final.weight) %>%
  as_tibble()

cows_avgs <- cows %>% group_by (level) %>%
  summarize (avg_fat = mean(fat))

cows$avg_fat <- cows_avgs$avg_fat[match(cows$level, cows_avgs$level)]

fit_1 <- stan_glm (fat ~ level, data = cows, refresh = 0)
print (fit_1)

ggplot (cows) +
  geom_point (aes (x=level, y=fat)) +
  geom_point (aes (x=level, y=avg_fat), color="red") +
  geom_abline(intercept = coef(fit_1)[1], slope = coef(fit_1)[2]) +
  labs (title="Linear regression with Treatment level as only predictor", 
        x="Treatment level", y="Fat produced")

fit_2 <- stan_glm (fat ~ level + lactation + age + initial_wt, data = cows, refresh = 0)
print (fit_2, digits=2)
coefs2 <- coef(fit_2)

fit_3 <- stan_glm (fat ~ factor(level) + lactation + age + initial_wt, data = cows, refresh = 0)
print (fit_3, digits=2)
coefs3 <- coef(fit_3)
df_step3 <- data.frame (level_divs=c(unique(cows$level), 0.4), 
                       fat_divs=c(coefs3[1], 
                                  coefs3[1] + coefs3[2], 
                                  coefs3[1] + coefs3[3], 
                                  coefs3[1] + coefs3[4],
                                  coefs3[1] + coefs3[4]))
 ggplot (cows) +
  geom_point (aes (x=level, y=fat)) +
  geom_point (aes (x=level, y=avg_fat), color="red") +
  geom_step(data=df_step3, aes(x=level_divs-0.05, y=fat_divs), size=1) +
  labs (title="Regression with teatment level as categorical variable", 
        x="Treatment level", y="Fat produced")

fit_4 <- stan_glm (fat ~ factor(level) + lactation, data = cows, refresh = 0)
print (fit_4, digits=2)

fit_5 <- stan_glm (fat ~ factor(level) + lactation + age + initial_wt, data = cows, refresh = 0)
print (fit_5, digits=2)
coefs5 <- coef(fit_5)
ses5   <- se(fit_5)
f_treat_5 <- data.frame (level_divs=c(unique(cows$level)), 
                        fat_divs=c(0, coefs5[2], coefs5[3], coefs5[4]), 
                        se_divs =c(0, ses5[2],   ses5[3],   ses5[4]))
df_treat_step5 <- data.frame (level_divs=c(unique(cows$level), 0.4), 
                              fat_divs=c(0, coefs5[2], coefs5[3], coefs5[4], coefs5[4]))

ggplot (df_treat_5) +
  geom_pointrange(aes(x=level_divs, y=fat_divs, 
                      ymin = fat_divs-se_divs, ymax = fat_divs+se_divs), color='red', size=0.7) +
  geom_step(data=df_treat_step5, aes(x=level_divs-0.05, y=fat_divs), size=1) +
  geom_abline(aes(intercept = 0, slope = coefs2[2]), color='blue') +
  labs (title="Cows - Regression with level as categorical variable", 
        subtitle="Predictors: level, lactation, age, initial weight",
        x="level", y="treatment effect")

# Ch 20
electric<- read.table("electric", header=TRUE) %>% as_tibble() 
# don't use controls
x <- electric$treated.Pretest
y <- electric$treated.Posttest
z <- ifelse (electric$Supplement. == "R", 0, 1 )
grade <- electric$Grade

df_dat <- data.frame (pretreat=x, posttreat=y, treatment=z, grade)

plots <- electric_plots (df_dat, "Estimated effect of Supplement vs Replacement")
plots[[1]]
plots[[2]]

# Ch 21

sesame <- read.csv('https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Sesame/data/sesame.csv', header=T) %>%
  as_tibble() %>% 
  select (watched, encouraged, prelet, postlet, site, setting)

fit <- stan_glm (postlet ~ watched, data = sesame, refresh = 0)
print (fit) # 12.6
cat ("The effect of watching, if watching were randomized, is", round(coef(fit)[2]), "points",
"\nBUT watching IS NOT RANDOMIZED!!!")

itt_zt <- stan_glm (watched ~ encouraged, data = sesame, refresh = 0)
print (itt_zt, digits = 2)
paste0("Percentage of kids who ended up watching because they were encouraged and (compliers) = ", 100*round(coef (itt_zt)["encouraged"], 2),"%")
# another way, directly from data
aa <- sesame %>% filter (encouraged == 1) %>% 
  summarize (pct_enc_and_watched = sum(watched == 1)/n())
bb <- sesame %>% filter (encouraged == 0) %>% 
  summarize (pct_notenc_and_watched = sum(watched == 1)/n())
paste0("Percentage of kids who would not have watched the show if they were not encouraged =",
  round ((aa$pct_enc_and_watched - bb$pct_notenc_and_watched)*100), " %")


itt_zy <- stan_glm (postlet ~ encouraged, data = sesame, refresh = 0)
print (itt_zy, digits = 2)
paste0("Improvement for kids who were encouraged = ", round(coef (itt_zy)["encouraged"], 2), " points")
# if only encouraged kids improved by 2.83 points, and only 36% of kids who were encouraged
# actually watched, 
paste0 ("Improvement for kids who watched = ", 
        round(coef(itt_zy)["encouraged"]/coef (itt_zt)["encouraged"], 1), " points")

# Another approach
fit_2a <- stan_glm (watched ~ encouraged,  data = sesame, refresh = 0)
sesame$watched_hat <- fit_2a$fitted
fit_2b <- stan_glm (postlet ~ watched_hat, data = sesame, refresh = 0)
print (fit_2b, digits = 2)
paste0 ("The effect of watching on those who were encouraged is = ",
        round( coef(fit_2b)["watched_hat"], 1), " points")

# Adjust for covariates
fit_3a <- stan_glm (watched ~ encouraged + prelet + factor(site) + setting,
                     data = sesame, refresh = 0)
print (fit_3a, digits=2)
sesame$watched_hat <- fit_3a$fitted
fit_3b <- stan_glm (postlet ~ watched_hat + prelet + factor(site) + setting, 
                    data = sesame, refresh = 0)
print (fit_3b, digits = 2)
paste0 ("The estimated effect of watching Sesame Street on the letter recognition test is ", round(coef(fit_3b)[2], 0), " points")

###
# brms - properly estimates standard errors
###
f1 <- bf (watched ~ encouraged)
f2 <- bf (postlet ~ watched)
IV_brm_a <- brm (f1 + f2, data = sesame)
print (IV_brm_a)  # postlet_watched        8   +/-   4.7

f3 <- bf (watched ~ encouraged + prelet + setting + factor(site))
f4 <- bf (postlet ~ watched    + prelet + setting + factor(site))
IV_brm_b <- brm (f3 + f4, data = sesame)
print (IV_brm_b) # postlet_watched        14   +/-   3.9 

# Regression Discontinuity





############################
# End September/2020  $$$$$$$$$$
############################
par (mar=c(2, 3, 2, 1))
#
n_girls <- rbinom (1, 400, 0.488)
print(n_girls)

#
# loop
#
n_sims  <- 1000
n_girls <- rep(NA, n_sims)
for (s in 1:n_sims){
  n_girls[s] <- rbinom(1, 400, 0.488)
}
list(n_girls)
#
# or inside rbinom function
#
n_girls <- rbinom (1000, 400, 0.488)
#
# or replicate
#
n_girls <- replicate (1000, rbinom(1, 400, 0.488))

hist (n_girls)
sd (n_girls)
mean (n_girls)
#
# Complicate model
#
birth_type <- sample (c("fraternal twin", "identical twin", "single birth"),
                      size=400, replace = TRUE, 
                      prob=c(1/125, 1/300, 1-1/125-1/300))
#
# Loop
#
girls <- rep(NA, 400)
for (i in 1:400){
  if(birth_type[i]=="single birth"){
    girls[i] <- rbinom (1, 1, 0.488)
  } else if (birth_type[i]== "identical twin"){
    girls[i] <- 2*rbinom (1, 1, 0.495)
  } else if (birth_type[i]== "fraternal twin"){
    girls[i] <- rbinom (1, 2, 0.495)
  }
}
n_girls <- sum (girls)
print (n_girls)

# or
#
girls <- ifelse (birth_type == "single birth", rbinom (400, 1, 0.488), 
                 ifelse (birth_type == "identical twin", 2*rbinom (400, 1, 0.495),
                         rbinom (400, 2, 0.495)))
n_girls <- sum (girls)
print (n_girls)
#
# Simulate n_girls 1000 times
#
n_sims <- 1000
n_girls <- rep (NA, n_sims)
for (s in 1:n_sims){
  birth_type <- sample (c("fraternal twin", "identical twin", "single birth"),
                        size=400, replace = TRUE, 
                        prob=c(1/125, 1/300, 1-1/125-1/300))
  girls <- ifelse (birth_type == "single birth", rbinom (400, 1, 0.488), 
                   ifelse (birth_type == "identical twin", 2*rbinom (400, 1, 0.495),
                           rbinom (400, 2, 0.495)))
  n_girls[s] <- sum (girls)
}

hist (n_girls)
#
# Distributions: normal, lognormal, binomial, Poisson
#
n_sims <- 1000
y1 <- rnorm (n_sims, 3, 0.5)
y2 <- exp (y1)
y3 <- rbinom (n_sims, 20, 0.6)
y4 <- rpois (n_sims, 5)
par (mfrow=c(2, 2))
hist(y1, breaks = seq(floor(min(y1)), max (y1)+0.2, 0.2),
     main="1000 draws from normal dist. wiht mean 3, sd 0.5")
hist(y2, breaks = seq(0, max(y2)+ 5, 5),
     main="Corresponding lognormal dist")
hist(y3, breaks = seq(-0.5, 20.5, 1),
     main="1000 draws from binmial dist. with 20 tries, probablity of success 0.6")
hist(y4, breaks = seq(-0.5, max (y4)+1, 1),
     main="1000 draws from Poisson dist. wiht mean 5")
#
#  Continuous + Discrete simulations
#
male <- rbinom (1, 1, 0.48)
height <- ifelse(male==1, rnorm (1, 69.1, 2.9), rnorm (1, 63.7, 2.7))

N <- 10
male <- rbinom (N, 1, 0.48)
# should we use N to sample for men and women?
height <- ifelse(male==1, rnorm (N, 69.1, 2.9), rnorm (N, 63.7, 2.7))
avg_height <- mean (height)
#
# simulate average height 1000 times
#
par (mfrow=c(1, 1))
n_sims <- 1000
avg_height <- rep (NA, n_sims)
max_height <- rep (NA, n_sims)
for (s in 1:n_sims){
  N <- 10
  male <- rbinom (N, 1, 0.48)
  height <- ifelse(male==1, rnorm (N, 69.1, 2.9), rnorm (N, 63.7, 2.7))
  avg_height[s] <- mean (height) 
  max_height[s] <- max (height) 
}
#
# or use a function
#
height_sim <- function (N){
  male <- rbinom (N, 1, 0.48)
  height <- ifelse(male==1, rnorm (N, 69.1, 2.9), rnorm (N, 63.7, 2.7))
  return(mean(height))
}

avg_height <- replicate (1000, height_sim(N=10))
hist (avg_height)
# 
# mean, median, sd, mad sd
#
z <- rnorm (10000, 5, 2)
cat ("mean=", mean(z), ", median=", median(z), ", sd=", sd(z), ", mad sd=", mad(z))
#
# The sd of the population is 2, => the sample mean's sample error (se) is 2/sqrt(10000) = 0.02
# 
# Quantile  ???
#
quantile(z, c(0.25, 0.75))    # central 50% interval
quantile(z, c(0.025, 0.975))  # central 95% interval
quantile(z, c(0.49, 0.51))    # central 2% interval
#
# Bootstrapping
#
ID     <- seq(1, 10, 1)
earn   <- c(50000, 60000, 30000, 51000, 9000, 29000, 32000, 2000, 27000, 6530)
height <- c(74, 66, 64, 63, 64, 62, 73, 72, 72, 70)
male   <- c(1, 0, 0, 0, 0, 0, 1, 1, 1, 1)
earnings <- cbind (data.frame(ID, earn, height, male))

earnings <- read.csv('earnings.csv',header=T)
ratio <- median(earnings$earn[earnings$male==0]) / median (earnings$earn[earnings$male==1])
cat ("Ratio of Female/Male earnings", ratio)

n <- nrow (earnings)       
boot <- sample(n, replace = TRUE)
earn_boot  <- earnings$earn[boot]
male_boot  <- earnings$male[boot]
ratio_boot <- median(earn_boot[male_boot==0]) / median(earn_boot[male_boot==1])
print (ratio_boot)

Boot_ratio <- function (data){
  n <- nrow (data)       
  boot <- sample(n, replace = TRUE)
  earn_boot  <- earnings$earn[boot]
  male_boot  <- earnings$male[boot]
  ratio_boot <- median(earn_boot[male_boot==0]) / median(earn_boot[male_boot==1])
  return(ratio_boot)
}

n_sims <- 10000
output <- replicate (n_sims, Boot_ratio(data=earnings))
hist (output)
cat ("median =", median(output), "sd =", sd(output))

#
# Chapter 6 - Regression
#

#
# Simulate 20 points yi= a + b * xi + ei
# yi - observed value
# xi - predictor
# ei - error, normally distributed with mean=0, sd=sigma
#

x <- 1:20
n <- length (x)
a <- 0.2
b <- 0.3
sigma <- 0.5
y <- a + b * x + sigma * rnorm (n)
# Make data frame
fake <- data.frame (x, y)
# Fit generalized linear model
fit_1 <- stan_glm (y ~ x, data=fake)
print (fit_1, digits = 2)

a_hat <- coef(fit_1)[1]
b_hat <- coef(fit_1)[2]
y_hat <- a_hat + b_hat * x

y_post <- rep (NA, n)
for (i in 1:n){
  y_post[i] <- y_hat[i] + rnorm (1, 0, 1)
}

cat ("mean (y)       sd (y)       =>", round(mean(y), 2),       ", ", round(sd(y),       2),
   "\nmean (y_hat)   sd (y_hat)   =>", round(mean(y_hat), 2),   ", ", round(sd(y_hat),   2),
   "\nmean (y_post)  sd (y_post)  =>", round(mean(y_post), 2),  ", ", round(sd(y_post),  2),
   "\nmean (y-y_hat) sd (y-y_hat) =>", round(mean(y-y_hat), 2), ", ", round(sd(y-y_hat), 2))


xbar <- mean (fake$x)

p1 <- ggplot (fake, aes(x=x, y=y)) + geom_point() +
  geom_abline(intercept=a_hat, slope=b_hat, color="blue", size=1) +
  ggtitle("Data and Fitted Regression Line") +
  geom_text(x= xbar, y=a_hat + b_hat*xbar - 1, 
            label = paste ("y=", round(a_hat, 2), "+", round(b_hat, 2), "* x"))
p1

#
# Earnings 2D linear fit
#
earnings$earnk <- earnings$earn/1000

fit_2 <- stan_glm(earnk ~ height + male, data=earnings)
print(fit_2)

#
# Heights
#
heights <- read.table("Heights.txt", header=TRUE)
heights

fit_1 <- stan_glm(daughter_height ~ mother_height, data = heights)
print (fit_1)

n <- nrow (heights)
heights$mother_height_jitt   <- heights$mother_height   + runif (n, -0.5, 0.5)
heights$daughter_height_jitt <- heights$daughter_height + runif (n, -0.5, 0.5)

ab_hat <- coef (fit_1)

sigma_hat <- median(as.matrix (fit_1)[,3])

ggplot (heights, aes(x=mother_height_jitt, y=daughter_height_jitt)) + geom_point() +
  xlab("Mother's heights in inches") + ylab("Adult Daughter's height in inches") +
  geom_abline(intercept = ab_hat[1], slope = ab_hat[2]) 

#
# Simulatiing Midterm and Finals Exam data
#
n <- 1000
true_ability <- rnorm (n, 50, 10)
noise_1 <- rnorm (n, 0, 10)
noise_2 <- rnorm (n, 0, 10)
midterm <- true_ability + noise_1
final   <- true_ability + noise_2
exams <- data.frame (midterm, final)

fit_1 <- stan_glm(final ~ midterm, data = exams)
print (fit_1)

sims <- as.matrix (fit_1)
a_hat <- median(sims[,1])
b_hat <- median(sims[,2])

ggplot (exams, aes(x=midterm, y=final)) + geom_point() +
  xlab("Midterm Score") + ylab("Final Score") +
  geom_abline(intercept = a_hat, slope = b_hat)
#
# Chapter 7; Linear Regression w/ SinglePredictor
#
hibbs <- read.table ("hibbs.dat", header = T)

M1 <- stan_glm(inc.party.vote ~ growth, data = hibbs)
print (M1)

sims <- as.matrix (M1)
a_hat <- median(sims[,1])
b_hat <- median(sims[,2])

xbar <- mean (hibbs$growth)

xoff <-  0.5
yoff <- -0.5

ggplot (hibbs, aes(x=growth, y=inc.party.vote)) + geom_point() +
  xlab("Economic Growth") + ylab("Incumbent Party's Vote Share") +
  geom_abline(intercept = a_hat, slope = b_hat) +
  geom_text(x= xbar +xoff, y=a_hat + b_hat*xbar + yoff, 
            label = paste ("y=", round(a_hat, 1), "+", round(b_hat, 1), "* x")) +
  theme_bw()
#
# Hillary v Trump
# growth = 2%
# "prior": vote_share = 46.2 + growth_rate * 3.1
# Dem expected vote share in 2016 = 46.2 + 2*3.1 = 52.4
# with std = 3.9%
# What were chances of Dems getting > 50% of the vote?
#
1-pnorm (50, 52.4, 3.9)  # = 73%

# 
# Simulate fake data
#
 
a <- 46.3
b <- 3.1
sigma <- 3.9
x <- hibbs$growth
n <- length (x)

y <- a + b*x + rnorm (n, 0, sigma) 
fake <- data.frame(x, y)

fit <- stan_glm (y ~ x, data = fake)
print (fit)

b_hat <- coef (fit)["x"]  # => slope of fit, 4.1
b_se  <- se (fit)["x"]    # => standard error of slope of fit, 0.8

# b 68% pobability between 3.3 and 4.9, 95% probability between 2.5 and 5.7

cover_68 <- abs (b-b_hat) < b_se
cover_95 <- abs (b-b_hat) < 2*b_se
cat (paste("68% coverage: ", cover_68, "\n"))
cat (paste("95% coverage: ", cover_95, "\n"))

# Do this many times

n_fake <- 1000
cover_68 <- rep (NA, n_fake)
cover_95 <- rep (NA, n_fake)  

for (s in 1:n_fake){
  y <- a + b*x + rnorm (n, 0, sigma) 
  fake <- data.frame(x, y)
  fit <- stan_glm (y ~ x, data = fake, refresh = 0)
  b_hat <- coef (fit)["x"] 
  b_se  <- se (fit)["x"]
  cover_68[s] <- abs (b-b_hat) < b_se
  cover_95[s] <- abs (b-b_hat) < 2*b_se
}  
cat (paste("68% coverage: ", mean(cover_68), "\n"))
cat (paste("95% coverage: ", mean(cover_95), "\n"))

# Try t-distrubtion with 14 df..  RESULTS DONT REPLICATE ****
#
t_68 <- qt (0.84, n-2)  # 84 = 100- (100-68)/2
t_95 <- qt (0.975, n-2) # 97.5 = 100- (100-95)/2
for (s in 1:n_fake){
  y <- a + b*x + rnorm (n, 0, sigma) 
  fake <- data.frame(x, y)
  fit <- stan_glm (y ~ x, data = fake, refresh = 0)
  b_hat <- coef (fit)["x"] 
  b_se  <- se (fit)["x"]
  cover_68[s] <- abs (b-b_hat) < t_68*b_se
  cover_95[s] <- abs (b-b_hat) < t_95*b_se
}  
cat (paste("68% coverage: ", mean(cover_68), "\n"))
cat (paste("95% coverage: ", mean(cover_95), "\n"))

# 
# 7-3 
#
n_0 <- 20
y0 <- 2
y_0 <- rnorm (n_0, y0, 5)
y_0

se_0 <-  sd (y_0)/sqrt (n_0)

cat ("mean y_0=", round(mean(y_0), 2), 
     "\nsd   y_0=", round (sd(y_0), 2),
     "\nse   y_0=", round (se_0, 2),
     "\nTrue y0=", y0,", 2 se CI = [", round (mean(y_0)-2*se_0, 2), ",", round (mean(y_0)+2*se_0, 2), "]")
     
fake_0 <- data.frame(y_0)     
fit_0 <- stan_glm(y_0 ~ 1, data = fake_0)
print (fit_0, digits = 2)

# Another data set
n_1 <- 30
y1 <- 8
y_1 <- rnorm (n_1, y1, 5)
fake_1 <- data.frame(y_1)
y_1
se_1 <-  sd (y_1)/sqrt (n_1)
cat ("mean y_1=", round(mean(y_1), 2), 
     "\nsd   y_1=", round (sd(y_1), 2),
     "\nse   y_1=", round (se_1, 2),
     "\nTrue y1=", y1,", 2 se CI = [", round (mean(y_1)-2*se_1, 2), ",", round (mean(y_1)+2*se_1, 2), "]")

fake_1 <- data.frame(y_1)
fit_1 <- stan_glm(y_1 ~ 1, data = fake_1)
print (fit_1, digits = 2)

# mean and sd of difference between y_1 and y_0
diff <- mean (y_1) - mean (y_0)
se_1 <- sd (y_1)/sqrt (n_1)
se   <- sqrt (se_0^2 + se_1^2)

cat ("mean y_1-y_0=", round(diff, 2), "\nse   y_1-y_0=", round (se, 2))

#
# Combine data sets
n <- n_0 + n_1
y <- c(y_0, y_1)
x <- c(rep(0, n_0), rep(1, n_1))
fake <- data.frame (x, y)
fit <- stan_glm (y ~ x, data = fake)
print (fit, digits = 2)

sims <- as.matrix (fit)
a_hat <- median (sims[,1]) # 1.8
b_hat <- median (sims[,2]) # 7.4

y_0bar <- mean(y_0)  # 1.74
y_1bar <- mean(y_1)  # 9.25
ggplot (fake, aes(x=x, y=y)) + geom_point() +
  ggtitle ("Least squares regression on an indicator is the same as computing the difference in means") +
  geom_abline(intercept=a_hat, slope=b_hat) +  
  geom_abline(intercept=y_0bar, slope=0, linetype = "dashed") +  
  geom_abline(intercept=y_1bar, slope=0, linetype = "dashed") +  
  annotate("text", x = 0.1, y = 1,    label = "y0_bar=1.74") +
  annotate("text", x = 0.9, y = 10.5, label = "y1_bar=9.25") +
  annotate("text", x = 0.55, y = 7.,  label = "y=1.8 + 7.4x") +
  xlab("Indicator, x") + ylab("y")
#
# Chapter 8 - Fitting Regression Models
# RSS - Residual Sum of Squares
#

rss <- function (x, y, a, b){
  resid <- y - (a + b*x)
  return (sum(resid^2))
}

rss (hibbs$growth, hibbs$inc.party.vote, a=46.3, b=3.1)  # minimum rss at optimal a and b

#
x <- 1:10
y <- c(1, 1, 2, 3, 5, 8, 13, 21, 34, 55)
fake <- data.frame(x, y)
fit <- stan_glm(y ~ x, data = fake)
print (fit) # x: median: 5.1,  sd: 1.1 => 2 * sd CI is: [2.9, 7.3]

sims <- as.matrix (fit)
sims

quantile (sims[,"x"], c(0.025, 0.975)) # [2.6 7.4 ]

#
# Simulation 2/10/2020
#
# PlayerA:, 30% shooter, PlayerB:, 40% shooter, take 20 shorts each.
# What is the probability PlayerB outshoots PlayerA, (yB > yA) or (yB - YA > 0)
#
#
# Method1 - Exact
#

# Wide Format
p <- array (NA, c(21, 21))  # Matrix of probabilities
sum <- 0
for (i in 1:21){
  for(j in 1:21){
    p[i,j] <- dbinom(i-1, 20, 0.3)*dbinom(j-1, 20, 0.4)
    if (j>i){
      sum <- sum + p[i,j]
    }
  }
}
print (sum)

# Long Format
yA <- rep(NA, 21^2)
yB <- rep(NA, 21^2)
count <- 0
for (i in 1:21){
  for(j in 1:21){
    count <- count+1
    yA[count] <- i-1
    yB[count] <- j-1
  }
}
p <- dbinom (yA, 20, 0.3)*dbinom (yB, 20, 0.4)

subset <- yB > yA
sum (p[subset])

#
# Method2 - Normal Approximation
#
mean_of_yA <- 20 * 0.3
mean_of_yB <- 20 * 0.4
sd_of_yA <- sqrt (0.3*(1-0.3)*20)
sd_of_yB <- sqrt (0.4*(1-0.4)*20)
mean_of_diff <- mean_of_yB - mean_of_yA
sd_of_diff   <- sqrt (sd_of_yA^2 + sd_of_yB^2)

est_prob_yB_gt_yA <- 1 - pnorm (0.5, mean_of_diff, sd_of_diff)
print (est_prob_yB_gt_yA)

#
# Method3 - Simulation
#
#set.seed(1)
n_loop <- 1000
out <- rep (NA, n_loop)
for (loop in 1:n_loop){
  shotA <- rbinom (1, 20, 0.3)
  shotB <- rbinom (1, 20, 0.4)
  out[loop] <- shotB > shotA
}
print (mean(out))

# Or more concise
mean(replicate (n_loop, rbinom (1, 20, 0.4) > rbinom (1, 20, 0.3)))

#
# Chapter9 -  Prediction and Bayesian Inferecne
#

hibbs <- read.table ("hibbs.dat", header = T)

M1 <- stan_glm(inc.party.vote ~ growth, data = hibbs)
print (M1)

sims <- as.matrix (M1)

Median <- apply (sims, 2, median)
MAD_SD <- apply (sims, 2, mad)
print (cbind(Median, MAD_SD))

a     <- sims[, 1]
b     <- sims[, 2]
sigma <- sims[, 3]

hist (a) # histogram of intercept, a
hist (b) # histogram of splope, b
hist (sigma) # histogram of residual std, sigma

data1 <- data.frame(a, b)
  
ggplot (data1, aes(x=a, y=b)) + geom_point() +
  ggtitle ("Posterior Draws of Regression Coefficients a and b") +
  xlab("Intercept, a") + ylab("Slope, b")

#
# 9-2 Prediction and Uncertainity:
# -predict
# -posterior_linpred
# -posterior_predict
#
new_growth <- 2
new <- data.frame (growth=new_growth)
#
# predict - Point prediction 
#
y_point_pred <- predict (M1, newdata = new) # 52.4
# Same as
a_hat <- coef(M1)[1] # 46.3
b_hat <- coef(M1)[2] # 3.04
y_point_pred <- a_hat + b_hat*new_growth  # 52.3

Median <- median (y_point_pred) # 52.4
MAD_SD <- mad    (y_point_pred) # 0

win_prob <- mean (y_point_pred > 50) # 1 - no uncertainity
#
# posterior_linpred- Linear Predictor 
#
y_lin_pred <- posterior_linpred(M1, newdata = new)
# Same as
a <- sims[, 1]
b <- sims[, 2]
y_lin_pred <- a + b*new_growth 

Median <- median (y_lin_pred)  # 52.4
MAD_SD <- mad    (y_lin_pred)  # 0.97

win_prob <- mean (y_lin_pred > 50) # 98.5%
#
# posterior_predict -  Predictive Distribution
#
y_pred <- posterior_predict(M1, newdata = new)
# Same as
sigma <- sims[, 3]
n_sims <- nrow (sims)

y_pred <- y_lin_pred + rnorm (n_sims, 0, sigma)

Median <- median (y_pred)  # 52.4
MAD_SD <- mad    (y_pred)  # 4.0

win_prob <- mean (y_pred > 50)  # 73%

cat ("Predicted Clinton Win of 2 Pary Vote: ", round(Median, 1),
     "With S.E ", round(MAD_SD, 1), "\nPr (Cinton Win) = ", round(win_prob, 2))

#
# Predict Prob (Clinton Win) For a Range of GDP Values, from -2 t0 +4
#
new_grid <- data.frame(growth=seq (-2.0, 4.0, 0.5))
y_point_pred_grid <- predict (M1, newdata = new_grid)
y_lin_pred_grid   <- posterior_linpred (M1, newdata = new_grid)
y_pred_grid       <- posterior_predict (M1, newdata = new_grid)                                        
#
# Add more uncertainity by using uncertainlty in GDP prediction, x
#
gdp <- rnorm (n_sims, 2.0, 0.3)
y_hat <- a + b*gdp
y_pred <- rnorm (n_sims, y_hat, sigma)

Median <- median (y_pred)  # 52.4
MAD_SD <- mad    (y_pred)  # 4.1

win_prob <- mean (y_pred > 50)  # 71%

cat ("Predicted Clinton Win of 2 Pary Vote: ", round(Median, 1),
     "With S.E ", round(MAD_SD, 1), "\nPr (Cinton Win) = ", round(win_prob, 2))

#
# Simulating Uncertainity for linear predictors and predictoed values
#
vitals_full <- read.table ("vitals.dat", header = T, sep =',')
# Clean up data
minidata <- vitals_full[, c("weight","height","female","ethnicity","exercise","smokenow")]
minidata$weight[minidata$weight>990]=NA;
ok <- apply(is.na(minidata), 1, sum) == 0 # only rows that do not have NAs

vitals <- minidata[ok,]
fit_1 <- stan_glm(weight ~ height, data = vitals, refresh = 0)
print (fit_1)
#
# Center height at 66 inches
#
vitals$c_height <- vitals$height - 66
fit_2 <- stan_glm(weight ~ c_height, data = vitals, refresh = 0)
print (fit_2)

#
# Predict the weight of someone 70 inches tall (c_height = 4.0 in)
#
new <- data.frame (c_height=4.0)

point_pred_2 <- predict (fit_2, newdata = new) # 173 lb
Median <- median (point_pred_2)
MAD_SD <- mad    (point_pred_2)
cat ("Predicted Weight of 70 in. person: ", round(Median, 0),
     "lbs, With S.E ", round(MAD_SD, 1))


linpred_2    <- posterior_predict(fit_2, newdata = new)
Median <- median (linpred_2)
MAD_SD <- mad    (linpred_2)
cat ("Predicted average Weight of 70 in. people: ", round(Median, 0),
     "lbs, With S.E ", round(MAD_SD, 1))
#
# Bayesian Information Aggregation
#
theta_hat_prior <- 0.524
se_prior        <- 0.041

n <- 400
y <- 190
theta_hat_data <- y/n            # 0.475
se_data <- sqrt (y/n*(1-y/n)/n)  # 0.025

theta_hat_bayes <- (theta_hat_prior/se_prior^2 + theta_hat_data/se_data^2)/
  (1/se_prior^2 + 1/se_data^2)                     # 0.48
se_bayes <- sqrt (1/(1/se_prior^2 + 1/se_data^2))  # 0.21

#
# Ex 8-1
#
n_data <- 101
a_n <- rep (NA, n_data)
b_n <- rep (NA, n_data)
rss_const_a <- rep (NA, n_data)
rss_const_b <- rep (NA, n_data)
for (n in 1:n_data){
  a_n[n] <- 46.3 + (n-51)
  b_n[n] <-  3.1 + (n-51)
  rss_const_a[n] <- rss (hibbs$growth, hibbs$inc.party.vote, a=46.3, b=b_n[n])
  rss_const_b[n] <- rss (hibbs$growth, hibbs$inc.party.vote, a=a_n[n], b=3.1)
}

hibbs_rss <- data.frame(a_n, rss_const_b, b_n, rss_const_a)


indx = which (rss_const_b ==min (rss_const_b))
x_lab <- paste ("Intercept a, minimum at", min (a_n[indx]))

ggplot (hibbs_rss, aes(x=a_n, y=rss_const_b)) + geom_point() +
  ggtitle ("RSS - Residual Sum of Squares - Slope=3.1") +
  xlab(x_lab) + ylab("RSS") + 
  theme_bw()

indx = which (rss_const_a ==min (rss_const_a))
x_lab <- paste ("Slope b, minimum at", min (b_n[indx]))

ggplot (hibbs_rss, aes(x=b_n, y=rss_const_a)) + geom_point() +
  ggtitle ("RSS - Residual Sum of Squares - Intercept=46.3") +
  xlab(x_lab) + ylab("RSS") + 
  theme_bw()
#
# Ex 8-5
#
n_data <- nrow (hibbs)
a_hat_noi      <- rep (NA, n_data)
b_hat_noi      <- rep (NA, n_data)

avg_growth <- mean(hibbs$growth)
avg_vote   <- mean(hibbs$inc.party.vote)

avg_growth_noi <- (avg_growth*n_data   - hibbs$growth        )/(n_data -1)
avg_vote_noi   <- (avg_vote  *n_data   - hibbs$inc.party.vote)/(n_data -1)

for (n in 1:n_data){
  sum_num <- 0
  sum_den <- 0
  for(cnt in 1:n_data){
    if (cnt !=  n){
      sum_num <- sum_num + (hibbs$growth[cnt] - avg_growth_noi[cnt]) * 
        hibbs$inc.party.vote[cnt]
      sum_den <- sum_den + (hibbs$growth[cnt] - avg_growth_noi[cnt])^2 
    }
  }
  b_hat_noi[n] <- sum_num/sum_den
}

a_hat_noi <- avg_vote_noi - b_hat*avg_growth_noi

b_hat_influence <- abs((b_hat_noi - b_hat)*100/b_hat)
a_hat_influence <- abs((a_hat_noi - a_hat)*100/a_hat)
growth_fm_mean <- hibbs$growth - mean (hibbs$growth)

influence <- data.frame (growth_fm_mean, a_hat_influence, b_hat_influence)

ggplot (influence, aes(x=growth_fm_mean, y=a_hat_influence)) + geom_point() +
  ggtitle ("% Influence of point i on intercept, a_hat") +
  xlab("Growth from mean") + ylab("%Influence of point i") + 
  theme_bw()

ggplot (influence, aes(x=growth_fm_mean, y=b_hat_influence)) + geom_point() +
  ggtitle ("% Influence of point i on slope, b_hat") +
  xlab("Growth from mean") + ylab("%Influence of point i") + 
  theme_bw()
#
# Just kidding 
#
# a_hat <- sum (yi) - b_hat * sum (xi)
# b_hat <- sum (xi - xavg) * yi/(sum(xi-xavg)^2)
#
# Change of yi by 1 changes a_hat by 1/n (add 1 to numerator)
# Change of yi by 1 changes b_hat by sum (xi - xavg)* yi /(sum(xi-xavg)^2)
hibbs <- read.table ("hibbs.dat", header = T)

hibbs <- hibbs[order(hibbs$growth),]
M1 <- stan_glm(inc.party.vote ~ growth, data = hibbs)
b_hat <- coef(M1)[2]  # 3.02

hibbs$growth_fm_mean <- hibbs$growth - mean(hibbs$growth)

sum_num <- rep (NA, n_data)
for (n in 1:n_data){
  sum_num[n] <- 0
  for(cnt in 1:n_data){
    sum_num[n] <- sum_num[n] + 
      hibbs$growth_fm_mean[cnt] * ((hibbs$inc.party.vote[cnt]) + (n==cnt))
  }
}

hibbs$b_hat_prime     <- sum_num / sum(growth_fm_mean^2)
hibbs$b_hat_influence <- abs(b_hat_prime - mean(b_hat_prime))*100/mean(b_hat_prime)

ggplot (hibbs, aes(x=growth_fm_mean, y=b_hat_influence)) + geom_point() +
  ggtitle ("% Influence of point i on slope, b_hat") +
  xlab("Growth - Mean (Growth)") + ylab("%Influence of point i") + 
  theme_bw()

#
# Refit data n_data times, changing vote% by one for each one
#
b_hat <- rep (NA, n_data)

for (n in 1:n_data){
  hibbs$inc.party.vote[n] <- hibbs$inc.party.vote[n] + 1
  M1 <- stan_glm(inc.party.vote ~ growth, data = hibbs, refresh=0)
  b_hat[n] <-  coef(M1)[2]       
  hibbs$inc.party.vote[n] <- hibbs$inc.party.vote[n] - 1
}

hibbs$influence_glm <- abs((b_hat - mean (b_hat))*100/mean(b_hat))

ggplot (hibbs, aes(x=growth_fm_mean, y=influence_glm)) + geom_point() +
  ggtitle ("% Influence of point i on slope, b_hat - GLM") +
  xlab("Growth - Mean (Growth)") + ylab("%Influence of point i") + 
  theme_bw()


x <- rep (NA, 2*n_data)
y <- rep (NA, 2*n_data)
Analysis <- rep (NA, 2*n_data)
for (n in 1:n_data){
  x[n] <- hibbs$growth_fm_mean[n]
  x[n_data+n] <- hibbs$growth_fm_mean[n]
  y[n] <- hibbs$b_hat_influence[n]
  y[n_data+n] <- hibbs$influence_glm[n]
  Analysis[n] <- "By_Hand"
  Analysis[n_data+n] <- "GLM"  
}

testd <- data.frame (x, y, Analysis)

ggplot (testd, aes(x=x, y=y, color = Analysis)) + geom_point() +
  ggtitle ("% Influence of point i on slope, b_hat - GLM") +
  xlab("Growth - Mean (Growth)") + ylab("%Influence of point i") + 
  theme_bw()
#
# Class Ex
#
sd <- 0.5
a <- 2
b <- 0.5

x <- seq(1,5)
y <- rnorm(length (x), a + b*x, sd)
df <- data.frame (x, y)
fit <- stan_glm (y~x, data=df, refresh=0)
ab_hat <- coef (fit)
sims <- as.matrix (fit)
a_s     <- sims[, 1]
b_s     <- sims[, 2]
a_hat   <- median(sims[,1])
b_hat   <- median(sims[,2])
x_bar <- mean(x)

ggplot (df, aes(x=x, y=y)) + ylim (0, 7) +
  geom_abline(intercept=a_s, slope=b_s, color="black", size=.05) +
  geom_point(color = "white") +
  geom_abline(intercept=a, slope=b, color="red", size=1) +
  geom_abline(intercept=a_hat, slope=b_hat, color="blue", size=1) +
  ggtitle("Data/Regression Line")  +
  geom_text(x= x_bar, y=ab_hat[1] + ab_hat[2]*x_bar -0.5, color="blue",label = "a_hat + b_hat * x") +
  geom_text(x= x_bar, y=a+b*x_bar +0.5, color="red", label = "a + b * x")
#
# Ex9-4
#
hibbs <- read.table ("hibbs.dat", header = T)

M1 <- stan_glm(inc.party.vote ~ growth, data = hibbs)
print (M1)

sims  <- as.matrix (M1)
a_hat <- median(sims[,1])
b_hat <- median(sims[,2])
ab_hat <- coef (M1)
growth_new = 2.0

sigma_hat <- 3.9
n_data    <- nrow (hibbs)

growth_avg <- mean (hibbs$growth)

sum_dist_sq <- sum ((hibbs$growth - growth_avg)^2)

sigma_hat_linpred <- sigma_hat * 
  sqrt ((1/n_data + (growth_new - growth_avg)^2/sum_dist_sq))      # 0.98
sigma_hat_prediction <- sigma_hat * 
  sqrt ((1 + 1/n_data + (growth_new - growth_avg)^2/sum_dist_sq))  # 4.02

new <- data.frame(growth=growth_new)
point_pred <- predict           (M1, newdata = new)  # 52.3, se=0
lin_pred   <- posterior_linpred (M1, newdata = new)  # 52.3, se=1.0
post_pred  <- posterior_predict (M1, newdata = new)  # 52.3, se=4.25

cat ("sigma_hat_linpred=", round (sigma_hat_linpred, 2), 
     "sd (lin_pred) = ", round(sd(lin_pred), 2))
cat ("sigma_hat_prediction=", round(sigma_hat_prediction, 2), 
     "sd (post_pred) = ", round(sd(post_pred), 2))
#
# Ex 9-6
#
sigma <- 0.2
a     <- 0.0
b     <- 0.5
cost  <- 0.3

x_new <- 20
y_new_point <- a + b*x_new # 4

n <- 1000000
y <- a + (rnorm (n, 0.5, 0.2)-0.3)* x_new 
sum (y<0)/n*100 # 16%

#
# Bayesion Inference: beauty and sex ratio
#
n_data <- 3000
theta_data  <- 8 # %girls - %boys in sample
se_data     <- 3
theta_prior <- 0 # %girls - %boys in population
se_prior    <- 0.25

theta_bayes <- (theta_prior/se_prior^2 + theta_data/se_data^2)/
  (1/se_prior^2 + 1/se_data^2)                         # 0.06
se_bayes    <- 1/sqrt((1/se_prior^2 + 1/se_data^2))    # 0.25                                       # 0.25

cat ("theta_bayes =", round(theta_bayes,2), "\nse_bayes =", round(se_bayes, 2))

# Why?  What's se in a sample of 3000, with p ~ 0.5
se <- sqrt (0.5*(1-0.5)/3000) # 0.9%
# what is the se of the difference of two 1500 people samples?
se <- sqrt (0.5*(1-0.5)/1500 + 0.5*(1-0.5)/1500) # 1.8%
# in this example se ~ 3, i.e. the sample of beautiful people was ONLY 10% of 3000
se <- sqrt (0.5*(1-0.5)/300 + 0.5*(1-0.5)/2700) # 3%

#
# Priors in stan_glm
#
#
# Uniform Prior
#
hibbs <- read.table ("hibbs.dat", header = T)

M3 <- stan_glm(inc.party.vote ~ growth, data = hibbs, refresh = 0,
               prior_intercept = NULL, prior = NULL, prior_aux = NULL)
sims <- as.data.frame(M3)
a <- sims[,1]
b <- sims[,2]

df <- data.frame (a, b)
ggplot (df, aes(x=a, y=b)) + 
  geom_point(color = "black")  +  
  ggtitle("4000 posterior draws (a,b)") 
#
# Weakly informative prior - similar result to uniform prior
#
M4 <- stan_glm(inc.party.vote ~ growth, data = hibbs, refresh = 0,
               prior = normal(5, 5, autoscale = FALSE),
               prior_intercept = normal (50, 50, autoscale = FALSE))
print (M4)
sims <- as.data.frame(M4)
a <- sims[,1]
b <- sims[,2]

df <- data.frame (a, b)
ggplot (df, aes(x=a, y=b)) + 
  geom_point(color = "black")  +  
  ggtitle("4000 posterior draws (a,b)") 
#
# Informative prior
#
x <- seq(-2,2,1)
y <- c(50, 44, 50, 47, 56)
sexratio <- data.frame(x, y)
fit <- lm(y ~ x, data = sexratio)
print(fit)
# Linear Fit
fit1 <- stan_glm(y ~ x, data = sexratio, refresh=0)
print(fit1)

sims <- as.data.frame(fit1)
a_l <- sims[,1]
b_l <- sims[,2]

df_l <- data.frame (a_l, b_l) 
ggplot (df_l, aes(x=a_l, y=b_l)) + xlim (30,70) + ylim (-12, 15) +
  geom_point(color = "black")  +  
  ggtitle("Flat Prior: 4000 posterior draws (a,b)") 

ggplot (sexratio, aes(x=x, y=y)) +# ylim (0, 7) +
  geom_abline(intercept=a_l, slope=b_l, color="black", size=0.05) +
  geom_point(color = "red") +
  ggtitle("Flat Prior: Data/Regression Line") 

# Bayesian Informative Prior

fit_post <- stan_glm (y~x, data = sexratio, refresh=0,
                      prior = normal (0, 0.2, autoscale = FALSE),
                      prior_intercept = normal (48.8, 0.5, autoscale = FALSE))
print (fit_post)
sims <- as.data.frame(fit_post)
a_b <- sims[,1]
b_b <- sims[,2]

df_b <- data.frame (a_b, b_b)
ggplot (df_b, aes(x=a_b, y=b_b)) + xlim (30,70) + ylim (-12,15) +
  geom_point(color = "black")  +  
  ggtitle("Informative Prior: 4000 posterior draws (a,b)") 

ggplot (sexratio, aes(x=x, y=y)) +# ylim (0, 7) +
  geom_abline(intercept=a_b, slope=b_b, color="black", size=0.05) +
  geom_point(color = "red") +
  ggtitle("Informative Prior: Data/Regression Line") 

#
# Chapter 10
#

kidiq <- read.dta(file="kidiq.dta")

# One binary predictor: mom's high school status
fit1 = stan_glm (kid_score ~ mom_hs, data = kidiq, refresh=0)
print (fit1)
a_hat <- coef (fit1)[1] # 77.5
b_hat <- coef (fit1)[2] # 11.8

cat ("kid_score =", round (a_hat, 0), "+", round (b_hat, 0), "* mom_hs + error")

# One continuous predictor: mom's IQ
fit2 = stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit1)
a_hat <- coef (fit2)[1]
b_hat <- coef (fit2)[2]

cat ("kid_score =", round (a_hat, 0), "+", round (b_hat, 1), "* mom_IQ + error")


# Two predictors: mom's high school status + mom's IQ
fit3 = stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq, refresh=0)
print (fit3)
a_hat <- coef (fit3)[1]  # 25.9
b_hat <- coef (fit3)[2]  # 6.0
c_hat <- coef (fit3)[3]  # 0.56

cat ("kid_score =", round (a_hat, 0), "+", round (b_hat, 0), "* mom_hs +",
     round (c_hat, 1), "* mom_iq + error")
#
# Interactions
#
fit4 = stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, 
                 data = kidiq, refresh=0)
print (fit4)
a_hat <- coef (fit4)[1]
b_hat <- coef (fit4)[2]
c_hat <- coef (fit4)[3]
d_hat <- coef (fit4)[4]

cat ("kid_score =", round (a_hat, 1), "+", round (b_hat, 1), "* mom_hs +",
     round (c_hat, 1), "* mom_iq + ", round (d_hat, 1), " mom_hs * mom_iq + error")

cat ("If mom_hs = 0, kid_score =", round (a_hat, 1),       "+ ", 
     round (c_hat, 1),       "* mom_iq")
cat ("If mom_hs = 1, kid_score =", round (a_hat+b_hat, 1), "+ ", 
     round (c_hat+d_hat, 1), "* mom_iq")

kidiq$mom_hs <- as.factor (kidiq$mom_hs)
ggplot (kidiq, aes(x=mom_iq, y=kid_score, color=mom_hs) )+ 
  scale_color_manual(values=c("blue", "darkgreen")) +
  xlim (0, 150) + ylim (0, 150) +
  geom_point() +
  geom_abline(intercept=a_hat, slope=c_hat, color="blue", size=1) +
  geom_abline(intercept=a_hat+b_hat, slope=c_hat+d_hat, color="darkgreen", size=1) +
  ggtitle("Child's Test Score Predictors")  +
  xlab("Mother's IQ Score") + ylab("Childs's Test Score") +
  theme_bw()
#
# Indicator Variables
#
fit_1 <- stan_glm (weight ~ height, data = vitals, refresh = 0)
print (fit_1)

coefs_1 <- coef (fit_1)
predicted_1 <- coefs_1[1] + coefs_1[2]*66.0  # 153

new <- data.frame(height=66)
pred <- posterior_linpred (fit_1, newdata = new)
cat ("The average weight of all 66-inch tall people is", round(mean (pred), 1), 
     "pounds, with a se of", round (sd(pred), 1))
pred <- posterior_predict (fit_1, newdata = new)
cat ("Predicted weight of a 66-inch tall person is", round(mean (pred), 1), 
     "pounds, with a se of", round (sd(pred), 1))
#
# Center height by subtracting 66 inches
#
fit_2 <- stan_glm (weight ~ c_height, data = vitals, refresh = 0)
print (fit_2)
#
# Expand model, add indicator variable, sex
#
fit_3 <- stan_glm (weight ~ c_height + female, data = vitals, refresh = 0)
print (fit_3)
coefs_3 <- coef (fit_3)
# Weight of a 70 inch woman
predicted <- coefs_3[1] + coefs_3[2]*4.0 + coefs_3[3]*1.0  # 164
# or
new <- data.frame(c_height=4, female=1)
pred <- posterior_predict(fit_3, newdata = new)
cat ("Predicted weight of a 70-inch tall woman is", round(median (pred), 1), 
     "pounds, with a sd of", round (sd(pred), 1))
# Weight of a 70 inch man
new <- data.frame(c_height=4, female=0)
pred <- posterior_predict(fit_3, newdata = new)
cat ("Predicted weight of a 70-inch tall man is", round(median (pred), 1), 
     "pounds, with a sd of", round (sd(pred), 1))
#
# Multiple levels of categorical predictor, add ethnicity
#
fit_4 <- stan_glm (weight ~ c_height + female + factor (ethnicity), 
                   data = vitals, refresh = 0)
print (fit_4)
#
# Change baseline factor level
#
vitals$eth <- factor (vitals$ethnicity,
                      levels=c("white", "black", "hispanic", "other"))
fit_5 <- stan_glm (weight ~ c_height + female + eth, 
                   data = vitals, refresh = 0)
print (fit_5)
# or
vitals$eth_white    <- ifelse (vitals$ethnicity=="white",    1, 0)
vitals$eth_black    <- ifelse (vitals$ethnicity=="black",    1, 0)
vitals$eth_hispanic <- ifelse (vitals$ethnicity=="hispanic", 1, 0)
vitals$eth_other    <- ifelse (vitals$ethnicity=="other",    1, 0)

fit_6 <- stan_glm (weight ~ c_height + female + eth_black + eth_hispanic + eth_other, 
                    data = vitals, refresh = 0)
print (fit_6)
#
# 10-5 Uncertainity predicting congressional elections
#
# Use 86-88 elections to generate prior
# Using prior, and 90 data as likelihood, project Bayesian posterior
#
c_86 <- read.table("congress_86", header=F)
c_86$percent_dem <- c_86$V4/(c_86$V4+c_86$V5)
c_86$percent_dem <- ifelse (c_86$percent_dem < 0.1, 0.25, c_86$percent_dem)
c_86$percent_dem <- ifelse (c_86$percent_dem > 0.9, 0.75, c_86$percent_dem)

bad_86 <- c_86$V4==-9 | c_86$V5==-9
c_86$percent_dem[bad_86] <- NA

c_88 <- read.table("congress_88", header=F)
c_88$percent_dem <- c_88$V4/(c_88$V4+c_88$V5)
c_88$percent_dem <- ifelse (c_88$percent_dem < 0.1, 0.25, c_88$percent_dem)
c_88$percent_dem <- ifelse (c_88$percent_dem > 0.9, 0.75, c_88$percent_dem)

bad_88 <- c_88$V4==-9 | c_88$V5==-9
c_88$percent_dem[bad_88] <- NA
inc_88 <- ifelse (c_88$V3 == -9, 0, c_88$V3)

c_90 <- read.table("congress_90", header=F)
c_90$percent_dem <- c_90$V4/(c_90$V4+c_90$V5)
c_90$percent_dem <- ifelse (c_90$percent_dem < 0.1, 0.25, c_90$percent_dem)
c_90$percent_dem <- ifelse (c_90$percent_dem > 0.9, 0.75, c_90$percent_dem)

bad_90 <- c_90$V4==-9 | c_90$V5==-9
c_90$percent_dem[bad_90] <- NA
inc_90 <- ifelse (c_90$V3 == -9, 0, c_90$V3)

hist (c_86$percent_dem, breaks = seq(0, 1, 0.05))
hist (c_88$percent_dem, breaks = seq(0, 1, 0.05))
hist (c_90$percent_dem, breaks = seq(0, 1, 0.05))

data_88 <- data.frame(past_vote=c_86$percent_dem, 
                      vote     =c_88$percent_dem,
                      inc      =inc_88)

ggplot (data_88, aes(x=past_vote, y=vote, color=factor(inc))) + 
  geom_point() + xlim (0, 1) + ylim (0,1) +
  ggtitle ("Adjusted Data") +
  xlab("Democratic vote share in 1986") + ylab("Democratic vote share in 1988") + 
  geom_abline(intercept=0, slope=1, color="black", size=0.5) +
  theme_bw()

# Fit model
fit_88 <- stan_glm (vote ~ past_vote + inc, data = data_88)
print (fit_88, digits=3)

a_hat <- coef (fit_88)[1]  # 0.2 - Even if dems got 0% in a district in 86, 
                           #       they expect to ge 20% in 88
b_hat <- coef (fit_88)[2]  # 0.6 - all else same, a district that got 10% more than another
                           #       in 86 is expected to get 6% more in 88
c_hat <- coef (fit_88)[3]  # 0.08- all else same, a district has an incumbent in 88 
                           #       is expected to get 8% more than one that hasn't

ggplot (data_88, aes(x=past_vote, y=vote, color=factor(inc))) + 
  geom_point() + xlim (0, 1) + ylim (0,1) +
  ggtitle ("Adjusted Data") +
  xlab("Democratic vote share in 1986") + ylab("Democratic vote share in 1988") + 
  theme_bw() +
  geom_abline(intercept=a_hat, slope=b_hat, color="black", size=0.5) 
  
sims_88 <- as.matrix(fit_88)
#
# Predict win probability given incumbency for '90
#
data_90 <- data.frame(past_vote=c_88$percent_dem, inc =inc_90)
data_90 <- na.omit (data_90)

pred_90 <- posterior_predict(fit_88, newdata = data_90)

# standard error of posterior prediction of each district
pre_90_sd <- apply (pred_90, 2, sd)  # 0.65-0.69
hist (pre_90_sd)

dems_pred <- rep (NA, n_sims)
for (s in 1:n_sims){
  dems_pred[s] <- sum (pred_90[s,] > 0.5)
}

cat ("Number of seats Dems expected to win in '90 =", round (mean(dems_pred)),
     ", se =", round (sd(dems_pred), 1) )

#
# Ex 9-1a
#
theta_hat_prior <- 0.42
theta_hat_data  <- 0.54
theta_hat_bayes <- 0.49
se_prior        <- 0.05                      

theta_hat_bayes_err <- function(n, theta_hat_prior, theta_hat_data, theta_hat_bayes, se_prior) {
  se_data <- sqrt (theta_hat_data * (1-theta_hat_data)/n)
  err <- (theta_hat_prior/se_prior^2 + theta_hat_data/se_data^2)/
    (1/se_prior^2 + 1/se_data^2) - theta_hat_bayes
  return(err)
  }

root <- uniroot(theta_hat_bayes_err, 
                interval=c(0,1000),
                theta_hat_prior = 0.42, 
                theta_hat_data  = 0.54, 
                theta_hat_bayes = 0.49, 
                se_prior        = 0.05,
                tol= 0.000000000000000001)$root 
n <- round (root)  # 139

se_data  <- sqrt (theta_hat_data * (1-theta_hat_data)/n)     # 4.2%
se_bayes <- sqrt (1/(1/se_prior^2 + 1/se_data^2))            # 3.2%

#se_prior = 5.0%
#se_data  = 4.2%
#se_bayes = 3.1%

# 
# Check
#
theta_hat_prior <- 0.42
se_prior        <- 0.05
theta_hat_data  <- 0.54 

n <- 144

se_data <- sqrt (theta_hat_data*(1-theta_hat_data)/n)      # 0.042
se_bayes <- sqrt (1/(1/se_prior^2 + 1/se_data^2))          # 0.032
theta_hat_bayes <- (theta_hat_prior/se_prior^2 + theta_hat_data/se_data^2) * se_bayes^2  # 0.49
 
#
# 9-1b
#
n_sim <- 100000

theta_hat_bayes <- 0.49
se_bayes        <- 0.032

sum(rnorm (n_sim, theta_hat_bayes, se_bayes) > 0.5)/n_sim # 38%
#
# 9-7
#

# final_score ~ N (70, 25)
# midterm_score ~ N (70, 25)
# slope ~ N (0.65, 0.15)
# intercept ~ N (70, 20) - after centering scores around 70

fit <- stan_glm(final_score ~ midterm_score, data = scores, refresh = 0,
                prior           = normal (0.65, 0.15, autoscale = T),
                prior_intercept = normal (75,   20,   autoscale = T))
#
# 10-1
#
var1 <- rnorm (1000, 0, 1)
var2 <- rnorm (1000, 0, 1)
cor (var1, var2) # 0.29

vars <- data.frame(var1, var2)
fit <- stan_glm (var2 ~ var1, data = vars, refresh=0)
coef(fit) [1]
coef(fit) [2]

ggplot (vars, aes(x=var2, y=var1)) + 
  geom_point() + 
  theme_bw() +
  geom_abline(intercept=0, slope=1, color="black", size=0.5) 

#
# 10-2
#
n <- 100
z_scores <- rep (NA, n)
for (k in 1:n){
  var1 <- rnorm (1000, 0, 1)
  var2 <- rnorm (1000, 0, 1)  
  fake <- data.frame(var1, var2)
  fit <- stan_glm (var2 ~ var1, data = fake, refresh=0)
  z_scores[k] <- coef(fit)[2]/se(fit)[2]
}
z_scores
sum (abs(z_scores)>2)

#
# Chapter 11 
#

# One predictor - mom_iq
fit_2 <- stan_glm (kid_score ~ mom_iq, data = kidiq)
print (fit_2)

ggplot (kidiq, aes(x=mom_iq, y=kid_score)) + geom_point() + 
  xlab ("Mother IQ Score") + ylab ("Child test score") +
  theme_bw() +
  geom_abline(intercept=coef(fit_2)[1], slope=coef(fit_2)[2], color="black", size=0.5) 

# Two predictors, mom_iq, mom_hs
fit_3 <- stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq)
print (fit_3)
b_hat <- coef(fit_3)
ggplot (kidiq, aes(x=mom_iq, y=kid_score, color = mom_hs)) + geom_point() + 
  xlab ("Mother IQ Score") + ylab ("Child test score") +
  geom_abline(intercept=b_hat[1] + b_hat[2], slope=b_hat[3], color=3, size=0.5) +
  geom_abline(intercept=b_hat[1], slope=b_hat[3], color=2, size=0.5) 

# Two predictors, with interaction
fit_4 <- stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, data = kidiq)
print (fit_4)
b_hat <- coef(fit_4)
ggplot (kidiq, aes(x=mom_iq, y=kid_score, color = mom_hs)) + geom_point() + 
  xlab ("Mother IQ Score") + ylab ("Child test score") +
  theme_bw() +
  geom_abline(intercept=b_hat[1] + b_hat[2], slope=b_hat[3]+ b_hat[4], color=3, size=0.5) +
  geom_abline(intercept=b_hat[1], slope=b_hat[3], color=2, size=0.5) 

# Display uncertainity in one predictor model
fit_2 <- stan_glm (kid_score ~ mom_iq, data = kidiq)
print (fit_2)

sims_2 <- as.matrix (fit_2)
beta_hat_2 <- apply (sims_2, 2, median)

ggplot (kidiq, aes(x=mom_iq, y=kid_score)) + 
  theme_bw() +
  geom_abline(intercept=sims_2[,1], slope=sims_2[,2], color="gray", size=0.1) +
  geom_abline(intercept=beta_hat_2[1], slope=beta_hat_2[2], color="black", size=0.5) +
  geom_point () + 
  xlab ("Mother IQ Score") + ylab ("Child test score") 

# Display one plot per input variable, holding other variable constant at its average
fit_3 <- stan_glm (kid_score ~ mom_hs + mom_iq, data = kidiq)
print (fit_3)

sims_3 <- as.matrix(fit_3)
n_sims_3 <- nrow (sims_3)

avg_hs <- mean(as.numeric(as.character(kidiq$mom_hs)))
ggplot (kidiq, aes(x=mom_iq, y=kid_score)) +
  ggtitle ("Child Test Score - Mean mother graduation rate") +
  xlab ("Mother IQ Score") + ylab ("Child test score") +
  theme_bw() +
  geom_abline(intercept=sims_3[,1] + avg_hs*sims_3[,3], slope=sims_3[,3], color="gray", size=0.1) +
  geom_abline(intercept=median(sims_3[,1]) + avg_hs*median(sims_3[,3]), 
              slope=median(sims_3[,3]), color="black", size=0.5) +
  geom_point() 

avg_iq <- mean(kidiq$mom_iq)
ggplot (kidiq, aes(x=mom_hs, y=kid_score))  + ylim (0, 150) +
  ggtitle ("Child Test Score - Mean mother IQ") +
  xlab ("Mother High School status") + ylab ("Child test score") +
  geom_abline(intercept=sims_3[,1]+avg_iq*sims_3[,3], 
              slope=sims_3[,3], color="gray", size=0.5) +
  geom_abline(intercept=median(sims_3[,1])+avg_iq*median(sims_3[,3]), 
              slope=median(sims_3[,3]), color="black", size=0.5) +
  geom_point() 
#
# Outcome vs Continuous Predictor + Binary treatment
#
# y = a + b*x + theta*z + error
# x - pretreatment predictor
# z - treatment indicator
#
N <- 100
x <- runif (N, 0, 1)
z <- sample (c(0, 1), N, replace = TRUE)
a <- 1
b <- 2
theta <- 5
sigma <- 2

y <- a + b*x + theta*z + rnorm (N, 0, sigma)

y_z0 <- a + b*x + rnorm (N, 0, sigma)
y_z1 <- a + b*x + theta + rnorm (N, 0, sigma)

yz <- data.frame(x, y, y_z0, y_z1)
ggplot (yz, aes(x=x, y=y_z0))  + ylim (-5, 15) +
  ggtitle ("z=0: No Treatment") + 
  xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  theme_bw() +
  geom_abline(intercept=a, slope=b, color="black", size=0.1) +
  geom_point() 

ggplot (yz, aes(x=x, y=y_z1))  + ylim (-5, 15) +
  ggtitle ("z=1: Treatment") + xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  theme_bw() +
  geom_abline(intercept=a+theta, slope=b, color="black", size=0.1) +
  geom_point() 

#
# Outcome vs Multiple Continuous Predictors + Binary treatment
#
# y = a + sum(bi*xi) + theta*z + error
# xi - pretreatment predictor, i
# z - treatment indicator
N <- 100
K <- 10
X <- array(runif (N*K, 0, 1), c(N,K))
z <- sample (c(0, 1), N, replace = TRUE)
a <- 1
b <- 1:K
theta <- 5
sigma <- 2

y <- a + X %*% b + theta*z + rnorm (N, 0, sigma)

fake <- data.frame(X=X, y=y, z=z)
fit  <- stan_glm(y ~ X + z, data = fake)
print (fit)

y_hat <- predict (fit)

y_y_hat_z_0 <- data.frame(y=y[z==0], y_hat=y_hat[z==0])
y_y_hat_z_1 <- data.frame(y=y[z==1], y_hat=y_hat[z==1])

p0 <- ggplot (y_y_hat_z_0, aes(x=y_hat, y=y))  + xlim (10, 60) +ylim (10, 60) +
  ggtitle ("z=0: No Treatment") + 
  xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  geom_abline(intercept=0, slope=1, color="black", size=0.5) + geom_point(shape=19) 

p1 <- ggplot (y_y_hat_z_1, aes(x=y_hat, y=y))  + xlim (10, 60) +ylim (10, 60) +
  ggtitle ("z=1: Treatment") + 
  xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  geom_abline(intercept=0, slope=1, color="black", size=0.5) + geom_point(shape=1) 

grid.arrange (p0, p1, nrow=1, ncol=2)

# Plot Residuals, y-y_bar
#
p0 <- ggplot (y_y_hat_z_0, aes(x=y_hat, y=y-y_hat))  + xlim (10, 50) +ylim (-15, 15) +
  ggtitle ("z=0: No Treatment") + 
  xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="black", size=0.5) +
  geom_point(shape=19) 

p1 <- ggplot (y_y_hat_z_1, aes(x=y, y=y-y_hat))  + xlim (10, 50) +ylim (-15, 15) +
  ggtitle ("z=1: Treatment") + 
  xlab ("Pre-treatment predictor, x") + ylab ("Outcome, y") +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="black", size=0.5) +
  geom_point(shape=1) 

grid.arrange(p0, p1, ncol=2)

# Residuals for kidiq example: kid_score ~ mom_iq
#
fit2 = stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit1)

kid_score_hat <- predict (fit2)

res_mean <- mean (kidiq$kid_score-kid_score_hat)
res_sd   <- sd   (kidiq$kid_score-kid_score_hat)

fake <- data.frame(y = kidiq$kid_score, y_hat = kid_score_hat)

ggplot (fake, aes(x=y_hat, y=y-y_hat))+
  xlab ("Mother's IQ score") + ylab ("Residuals, y") 
  geom_abline(intercept=-res_sd, slope=0, color="black", linetype = "dashed", size=0.5) +
  geom_abline(intercept=0,       slope=0, color="black",                      size=0.5) +
  geom_abline(intercept= res_sd, slope=0, color="black", linetype = "dashed", size=0.5) +
  geom_point() 
#
# Introclass
#
introclass <- read.table("introclass", header=TRUE)

fit_1 <- stan_glm (Final ~ Midterm, data = introclass)
print (fit_1)

sims <- as.matrix(fit_1)

y      <- introclass$Final
y_pred <- colMeans(sims[,1] + sims[,2] %*% t(introclass$Midterm))

r <- y-y_pred  

fake <- data.frame(y=y, y_pred=y_pred, r=r)

p0 <- ggplot (fake, aes(x=y_pred, y=r))+
    ggtitle ("Predicted Value vs Residual") +
    xlab ("Predicted Value") + ylab ("Residual") +
    theme_bw() +
    geom_abline(intercept=0, slope=0, color="black", size=0.5) +
    geom_point() 
p1 <- ggplot (fake, aes(x=y, y=r))+
  ggtitle ("Observed Value vs Residual") +
  xlab ("Observed Value") + ylab ("Residual") +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="black", size=0.5) +
  geom_point() 

grid.arrange (p0, p1, ncol=2)  

# Fake data simulation: Model is by definition "Correct"
#
a <- 64.5
b <- 0.7
sigma <- 14.8
n <- nrow (introclass)
introclass$final_fake <- a + b*introclass$Midterm + rnorm (n, 0, sigma)

fit_fake <- stan_glm(final_fake ~ Midterm, data = introclass)  
print (fit_fake)  

sims <- as.matrix (fit_fake)  
predicted_fake <- colMeans((sims[,1] + sims[,2] %*% t(introclass$Midterm)))  
  
fake <- data.frame(y=introclass$final_fake, y_pred=predicted_fake, 
                   r=introclass$final_fake-predicted_fake)  
  
p0 <- ggplot (fake, aes(x=y_pred, y=r))+
  ggtitle ("Predicted Value vs Residual") +
  xlab ("Predicted Value") + ylab ("Residual") +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="black", size=0.5) +
  geom_point() 
p1 <- ggplot (fake, aes(x=y, y=r))+
  ggtitle ("Observed Value vs Residual") +
  xlab ("Observed Value") + ylab ("Residual") +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="black", size=0.5) +
  geom_point() 
grid.arrange (p0, p1, ncol=2)    

#
# Section 11-4
#
newcomb <- data.frame(y=c(28,26,33,24,34,-44,27,16,40,-2,29,22,24,21,25,30,
                          23,29,31,19,24,20,36,32,36,28,25,21,28,29,37,25,28,
                          26,30,32,36,26,30,22, 36,23,27,27,28,27,31,27,26,
                          33,26,32,32,24,39,28,24,25,32,25,29,27,28,29,16,23))  

fit <- stan_glm(y ~ 1, data = newcomb)  
print (fit)  

sims <- as.matrix (fit)
n_sims <- nrow (sims)

n <- length(newcomb$y)
y_rep <- array (NA, c(n_sims, n))
for (s in 1:n_sims){
  y_rep[s,] <- rnorm(n, sims[s,1], sims[s,2])
}

# OR!
y_rep <- posterior_predict (fit)
  
par (mfrow=c(5,4))
for (s in sample (n_sims, 20)){
  hist(y_rep[s,])
}
par (mfrow=c(1,1))
#
# Create test function
#
Test <- function (y){
  min (y)
}
test_rep <- apply (y_rep, 1, Test)  

hist (test_rep, xlim=range(Test(newcomb$y), test_rep))
abline (v=Test(newcomb$y), lwd = 3)

#
# 11-5 Predictive Simulation
#
unemp <- read.table("unemp", header=TRUE)
head (unemp)

n <- nrow (unemp)
y <- unemp$unemployed.pct
unemp$y_lag <- c(NA, y[1:(n-1)])
fit_lag <- stan_glm(unemp$unemployed.pct ~ y_lag, data = unemp)

print (fit_lag, digits=2)

sims <- as.matrix(fit_lag)
n_sims <- nrow (sims)

# Replicate datasets
y_rep <- array(NA, c(n_sims, n))
for (s in 1:n_sims){
  y_rep[s,1] <- y[1]
  for (t in 2:n){
    y_rep[s,t] <- sims[s, "(Intercept)"] + sims[s, "y_lag"]*y_rep[s, t-1] +
      rnorm(1, 0, sims[s, "sigma"])
  }
}

# Can NOT use posterior predict
#
idx <- sample (c(1:4000), size=15, replace = TRUE)

plot_list = list()
for (s in 1:15){
  i<- idx[s]
  y_rep_i <- y_rep[i, ]
  fake <- data.frame(y=y_rep_i, x=unemp$year)
  title <- paste ("Sim:", i) 
  plot_list[[s]] <- ggplot (fake, aes(x=x, y=y)) + 
    ggtitle(title) + geom_line() 
}

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
          plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], 
          plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], 
          ncol=5, nrow=3)    
  
#
# Test Function
#
Test <- function (y){
  n <- length(y)
  y_lag   <- c(NA,     y[1:(n-1)])
  y_lag_2 <- c(NA, NA, y[1:(n-2)])
  return (sum(sign(y-y_lag) != sign(y_lag-y_lag_2), na.rm = TRUE))
}

test_y <- Test (y)
test_rep <- apply (y_rep, 1, Test)

round (sum (test_rep > test_y)/length (test_rep) * 100 )  # 97% of simulations change direction
                                                          # more than actual data

#
# Residual standard deviation
# R^2 = 1 - (sigma_hat)^2/(s_y)^2  
# If model fit uses least squares:
##### R^2 = Var(y_hat)/Var(y) #####
# R^2 = 1/(n-1) * sum (y_hat_i - y_hat_bar)^2/s_y^2
#
idiq <- read.dta(file="kidiq.dta")
# One predictor: mom's high school status
fit1 = stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
sims <- as.matrix (fit1)

a_hat     <- median(sims[,1]) 
b_hat     <- median(sims[,2]) 
sigma_hat <- median(sims[,3])  

n <- nrow (kidiq)

y         <- kidiq$kid_score
y_hat     <- a_hat + b_hat * kidiq$mom_iq
y_hat_bar <- mean (y_hat)

#sigma_hat    18.3%
sd (y)      # 20.0%
sd (y_hat)  #  9.1%

R2 <- 1 - ((sigma_hat)/sd(y))^2                      # 19.7 %
# OR - IF LINEAR REGRESSION, i.e. y_hat_i = a_hat + x_i * beta_hat 
R2 <- 1/(n-1) * sum ((y_hat - y_hat_bar)^2)/sd(y)^2  # 20.0 %

# 
# Bayesian R^2
#
# Bayesian R^2 = Var(y_hat)/[(Var(y_hat) + Var(r)]
# r - y - y_hat (residual)
#

x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3
xy <- data.frame(x,y)

fit <- stan_glm(y ~ x, data = xy)
print (fit)

sims <- as.matrix (fit)

a_hat     <- median(sims[,1])  # 0.0083
b_hat     <- median(sims[,2])  # 0.589
sigma_hat <- median(sims[,3])  # 0.699
y_hat <- a_hat + b_hat * x

r <- y - y_hat

fake <- data.frame(y, y_hat)

ggplot (fake, aes(x=y, y=y_hat)) + geom_point() +
  geom_abline(intercept=a_hat, slope=b_hat, color="black", size=1) +
  xlab ("y") + ylab("y_bar") 
  
sd (r)      #  0.525
sd (y)      #  1.084
sd (y_hat)  #  0.931

rsq_1 <- var(y_hat)/(var(y))               # 74%
rsq_2 <- var(y_hat)/(var(y_hat) + var(r))  # 76%

# Bayesian fit with prior
#
fit_bayes <- stan_glm(y ~ x, data = xy,
                      prior_intercept = normal(0, 0.2, autoscale = FALSE),
                      prior           = normal(1, 0.2, autoscale = FALSE),
                      prior_aux = NULL, refresh = 0)

posterior  <- as.matrix(fit_bayes, pars = c("(Intercept)", "x"))
post_means <- colMeans(posterior)  # -0.002, 0.837

median (bayesR2<-bayes_R2(fit_bayes))  # 80%

ggplot (xy, aes(x=x, y=y)) + geom_point() +
  geom_abline(intercept=a_hat,         slope=b_hat, color="black", size=1) +
  geom_abline(intercept=0,             slope=1,     color="black", linetype = "dashed", size=1) +
  geom_abline(intercept=post_means[1], slope=post_means[2], color="black", linetype = "dotted", size=1) +
  ggtitle ("Data, Linear fit (solid), Beyesian prior (dashed), posterior (dotted)")
  xlab ("x") + ylab("y") 

a_bayes <- rep (NA, 20)
b_bayes <- rep (NA, 20)
idx <- sample (c(1:4000), size=20, replace = TRUE)
for (s in 1:20){
  i<- idx[s]
  a_bayes[s] <- posterior[i, 1]
  b_bayes[s] <- posterior[i, 2]
}

ggplot (xy, aes(x=x, y=y)) + geom_point() +
  geom_abline(intercept=a_bayes,      slope=b_bayes, color="black", size=0.1) +
  geom_abline(intercept=post_means[1], slope=post_means[2], color="black", linetype = "solid", size=1) +
  ggtitle ("Bayesian Posterior Mean and 20 draws from posterior distribution") +
  xlab ("x") + ylab("posterior y (Bayes)") 
#
# 11-8 Cross Validation and Deviance
#
fit_3 <- stan_glm(kid_score ~ mom_hs + mom_iq, data=kidiq, refresh = 0)
print(fit_3, digits = 2)

# Calculate Log Score
sigmas <- as.matrix(fit_3)[,'sigma']
preds <- posterior_linpred(fit_3)
nsims <- nrow(preds)

logscore_3 <- sum(log(rowMeans(sapply(1:nsims, FUN = function(i) dnorm(kidiq$kid_score, preds[i,], sigmas[i], log=FALSE)))))
round(logscore_3, 1)  # -1871.2

loo_3 <- loo (fit_3)
print (loo_3)   # elpd_loo = -1876, ~ (logscore_3 - 4), 4 is p_loo, num_parameters
#
# Add 5 noise predictors
#
n <- nrow (kidiq)
kidiqr <- kidiq
kidiqr$noise <- array (rnorm(5*n), c(n, 5))

fit_3n <- stan_glm(kid_score ~ mom_hs + mom_iq + noise, data=kidiqr, refresh = 0) 
print(fit_3n, digits = 2)

loo_3n <- loo (fit_3n)
print (loo_3n)   # elpd_loo = -1880

# Calculate Logscore
sigmas <- as.matrix(fit_3n)[,'sigma']
preds <- posterior_linpred(fit_3n)
nsims <- nrow(preds)

logscore_3n <- sum(log(rowMeans(sapply(1:nsims, FUN = function(i) dnorm(kidiq$kid_score, preds[i,], sigmas[i], log=FALSE)))))
round(logscore_3n, 1)  # 1871

# Go back to simple model
fit_1 <- stan_glm (kid_score ~ mom_hs, data = kidiq)
loo_1 <- loo (fit_1)
print (loo_1)  # loo_1 elpd = -1915

#
# Compare models fit_1 anf fit_3
#
# Compute elpds
compare_models (loo_3, loo_1)   # elpd_diff = -38.83, se=8.3
# Compute Bayes R2
br2_1<-bayes_R2(fit_1)  
median (br2_1)           # 0.056
br2_3<-bayes_R2(fit_3)
median (br2_3)           # 0.214

df <- melt(data.frame(fit_1=br2_1,fit_3=br2_3))  #************#
ggplot(df, aes(x=value, linetype=variable)) +
  geom_density(alpha=0.25, show.legend=FALSE) +
  ggtitle("Bayes-R2 posteriors") +
  labs(x="Bayes-R^2", y="Density") +
  scale_y_continuous(breaks=NULL) + 
  annotate("text", x = 0.107, y = 16.2, label = "kid_score ~ mom_hs") +
  annotate("text", x = 0.28, y = 13.2, label = "kid_score ~ mom_hs + mom_iq")
#
# Compare models fit_3 anf fit_3n - adding arbitrary noise
#
# Compute elpds
compare_models (loo_3n, loo_3)   # elpd_diff = 4.1, se=1.4
# Compute Bayes R2
br2_3n<-bayes_R2(fit_3n)  
median (br2_3n)           # 0.224
br2_3<-bayes_R2(fit_3)
median (br2_3)            # 0.214

df <- melt(data.frame(fit_3=br2_3,fit_3n=br2_3n))
ggplot(df, aes(x=value, linetype=variable)) +
  geom_density(alpha=0.25, show.legend=FALSE) +
  ggtitle("Bayes-R2 posteriors") +
  labs(x="Bayes-R^2", y="Density") +
  scale_y_continuous(breaks=NULL) + 
  annotate("text", x = 0.15, y = 11.5, label = "kid_score ~ mom_hs + mom_iq") +
  annotate("text", x = 0.285, y = 12, label = "kid_score ~ mom_hs + mom_iq +\n noise")

#
# Define Loo R2
#
looR2 <- function(fit) {
  y <- rstanarm::get_y(fit)
  ypred <- posterior_linpred(fit)
  ll <- log_lik(fit)
  r_eff <- relative_eff(exp(ll), chain_id = rep(1:4, each = 1000))
  psis_object <- psis(log_ratios = -ll, r_eff = r_eff)
  ypredloo <- E_loo(ypred, psis_object, log_ratios = -ll)$value
  eloo <- ypredloo-y
  return(1-var(eloo)/var(y))
}
####
round(looR2(fit_1),3)  # 0.047  mom_hs
round(looR2(fit_3),3)  # 0.203  mom_hs + mom_iq
round(looR2(fit_3n),3) # 0.187  added noise

# Try fit_4
fit_4 <- stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, data = kidiq)
loo_4 <- loo (fit_4)
print (loo_4)  # loo_1 elpd = -1872

compare_models(loo_3, loo_4)  # elpd_diff=3.6, se=2.5 .. not substantially better

# Calculate Loo R2 fordifferent versions of model
round(looR2(fit_1), 3)  # 0.047  mom_hs
round(looR2(fit_3), 3)  # 0.203  mom_hs + mom_iq
round(looR2(fit_3n),3)  # 0.187  added noise
round(looR2(fit_4), 3)  # 0.216  mom_hs:mom_iq interaction
#
# K-fold cross validation
#
n   <- 60
k   <- 30
rho <- 0.8
Sigma <- rho * array(1, c(k,k)) + (1-rho)* diag(k)
X     <- mvrnorm (n, rep(0,k), Sigma)

b <- c(c(-1,1,2), rep(0, k-3))
y <- X %*% b + 2*rnorm(n)

fake <- data.frame(X,y)

# Use model with WEAK prior
fit_1 <- stan_glm (y ~ ., prior=normal(0, 10, autoscale = FALSE), data = fake)
loo_1 <- loo (fit_1)
print (loo_1)  # get nasty messages, use kfold

# Try kfold approach
kfold_1 <- kfold (fit_1, K=10)
print (kfold_1)

# Now use model with more INFORMTIVE prior, horseshoe, hs_prior

k0 <- 2 # prior guess for the number of relevant variables
tau0 <- k0/(k-k0) * 1/sqrt(n)
hs_prior <- hs(df=1, global_df=1, global_scale=tau0, slab_scale=3, slab_df=7)

fit_2 <- stan_glm(y ~ ., prior=hs_prior, data=fake)
print (fit_2)

loo_2 <- loo (fit_2)  
print (loo_2)  # get nasty messages, use kfold

# Try kfold approach
kfold_2 <- kfold (fit_2, K=10)
print (kfold_2)

compare_models(kfold_1, kfold_2)  # elpd_diff=26.2, se=6.1  => model2 is better    

#
# 11-3
#
n <- 100

a <- 2
b <- 3

n_loop <- 1000
test <- rep (NA, n_loop)
for (i in 1:n_loop){
  x <- runif (n, 0, 10) 
#  y <- a + b * x + rnorm (n, 0, 1)
    y <- a + b * x + sample (0:1, n, replace = T)
  xy <- data.frame(x, y)

  fit <- stan_glm(y~x, data = xy, refresh=0)

  sims <- as.matrix (fit)
  b_hat   <- median(sims[,2])
  test[i] <- abs(b_hat - b) < 2*sd (sims[,2])
}
  
sum (test)/n_loop

#
# CH 12
#
# Logarithmic transformations
#
earnings <- read.csv('earnings.csv',header=T)

# Eliminate earnings = 0
earnings <- earnings[!earnings$earn == 0, ]
earnings$log_earn <- log (earnings$earn)

# Linear fit
fit_1 <- stan_glm(earn ~ height, data = earnings)
print (fit_1, digits=2)

sims <- as.matrix(fit_1)
a_hat <- median(sims[,1])
b_hat <- median(sims[,2])

fake <- data.frame(x=earnings$height, y=earnings$earn)
ggplot (fake, aes(x=x, y=y)) + geom_point() +
  xlab("Height") + ylab("Earnings") +
  geom_abline(intercept = a_hat, slope = b_hat) 

# Log Fit
logmodel_1 <- stan_glm(log_earn ~ height, data = earnings)
print (logmodel_1, digits=2)

sims <- as.matrix(logmodel_1)
a_hat <- median(sims[,1])
b_hat <- median(sims[,2])

pick_10 <- sample (4000, 10)
fake <- data.frame(x=earnings$height, y=earnings$log_earn)

ggplot (fake, aes(x=x, y=y)) + geom_point() +
  xlab("Height") + ylab("Log (earnings)") +
  ggtitle ("Log Fit with 10 draws") +
  geom_abline(intercept = a_hat, slope = b_hat) +
  geom_abline(intercept = sims[pick_10,1 ], slope = sims[pick_10, 2], color="darkgray", size = 0.2) 

# Replicte data for linear fit
yrep_1 <- posterior_predict(fit_1)
n_sims <- nrow (yrep_1)

subset <- sample (n_sims, 100)
ppc_dens_overlay (earnings$earn, yrep_1[subset,])

# Replicate data for log fit
yrep_1 <- posterior_predict(logmodel_1)
n_sims <- nrow (yrep_1)

subset <- sample (n_sims, 100)
ppc_dens_overlay (earnings$log_earn, yrep_1[subset,])

# Add another predictor, male
logmodel_2 <- stan_glm(log_earn ~ height + male, data = earnings)
print (logmodel_2, digits=2)

sims <- as.matrix(logmodel_2)
b0_hat <- median(sims[,1])
b1_hat <- median(sims[,2])
b2_hat <- median(sims[,3])

y     <- earnings$log_earn
y_hat <- b0_hat + b1_hat * earnings$height + b2_hat * earnings$male
r     <- y_hat - y

new1 <- data.frame("height" = 70, "male" = 0)
y_point_pred <- predict (logmodel_2, newdata = new1)  # 9.59
y_lin_pred   <- posterior_linpred (logmodel_2, newdata = new1)
y_pred       <- posterior_predict (logmodel_2, newdata = new1)     

sd (r)  # 0.88: 68% of log earnings will be 0.88 of predicted value.  
        # For ex, woman of height=70, y_predict = 8.14 + 0.02 * 70 + 0 = 9.54, with sd 0.88
        # 68% chance log earnings will be [9.54 +/- 0.88] = [8.66, 10.42]
        # earnings ~ [6,000, 34,000]

rsq_2 <- var(y_hat)/(var(y_hat) + var(r)) # 0.09%: Only 9% of the varience of the
                                          # log transformed data is explained by model

yrep_2 <- posterior_predict(logmodel_2)
n_sims <- nrow (yrep_2)

subset <- sample (n_sims, 100)
ppc_dens_overlay (earnings$log_earn, yrep_2[subset,])

# Include Interaction
logmodel_3 <- stan_glm(log_earn ~ height + male + height:male, data = earnings)
print (logmodel_3, digits=3)

sims <- as.matrix(logmodel_3)

b0_hat <- median(sims[,1])
b1_hat <- median(sims[,2])
b2_hat <- median(sims[,3])
b3_hat <- median(sims[,4])

# Z score
earnings$z_height <- (earnings$height - mean (earnings$height))/sd(earnings$height)

logmodel_4 <- stan_glm(log_earn ~ z_height + male + z_height:male, data = earnings)
print (logmodel_4, digits=2)

sims <- as.matrix(logmodel_4)
b0_hat <- median(sims[,1])
b1_hat <- median(sims[,2])
b2_hat <- median(sims[,3])
b3_hat <- median(sims[,4])

# Log-Log transforms
earnings$log_height <- log (earnings$height)
logmodel_5 <- stan_glm(log_earn ~ log_height + male, data = earnings)
print (logmodel_5, digits=2)

y <- earnings$log_earn
y_bar <- coef (logmodel_5)[1] + 
  coef (logmodel_5)[2] * earnings$log_height + 
  coef (logmodel_5)[3] * earnings$male
r <- y-y_bar
#
# Ex 11-4
#
pyth <- read.table("pyth", header=TRUE)

pyth_no_na <- na.omit (pyth)


fit <- stan_glm (y ~ x1 + x2, data = pyth_no_na)
print (fit)

sims <- as.matrix(fit)
b0_hat <- median(sims[,1])
b1_hat <- median(sims[,2])
b2_hat <- median(sims[,3])

x1bar <- mean (pyth_no_na$x1)
x2bar <- mean (pyth_no_na$x2)

ggplot (pyth_no_na, aes(x=x1, y=y)) + geom_point() +
  ggtitle ("y=b0+b1*x1 + b2*x2 vs x1, x2 = x2bar") +
  xlab("x1") + ylab("y") +
  geom_abline(intercept = b0_hat + b2_hat * x2bar, slope = b_hat) 

ggplot (pyth_no_na, aes(x=x2, y=y)) + geom_point() +
  ggtitle ("y=b0+b1*x1 + b2*x2 vs x2, x1 = x1bar") +
  xlab("x2") + ylab("y") +
  geom_abline(intercept = b0_hat + b1_hat * x1bar, slope = b2_hat) 

r <- pyth_no_na$y - (b0_hat + b1_hat*pyth_no_na$x1 + b2_hat*pyth_no_na$x2)

pyth_no_na <- cbind(pyth_no_na, r)

ggplot (pyth_no_na, aes(x=x1, y=r)) + geom_point() +
  ggtitle ("residual vs x1") + xlab("x1") + ylab("residual") +
  geom_abline(intercept = 0, slope = 0) 
ggplot (pyth_no_na, aes(x=x2, y=r)) + geom_point() +
  ggtitle ("residual vs x2") + xlab("x2") + ylab("y") +
  geom_abline(intercept = 0, slope = 0) 
#
# Posterior predict
#
new1 <- data.frame("x1" = 9.63, "x2" = 12.16)

b0_hat + b1_hat * new1$x1 + b2_hat * new1$x2

y_point_pred <- predict (fit, newdata = new1)
y_lin_pred   <- posterior_linpred (fit, newdata = new1)
y_pred       <- posterior_predict (fit, newdata = new1)     

pyth_na <-  pyth[41:60,]
pyth_na$y <- NULL
new <- pyth_na

y_point_pred <- predict (fit, newdata = new)
y_lin_pred   <- posterior_linpred (fit, newdata = new)
y_pred       <- posterior_predict (fit, newdata = new)  

cbind (pyth_na, y_point_pred)

x1bar <- mean (pyth_na$x1)
x2bar <- mean (pyth_na$x2)

ggplot (pyth_na, aes(x=x1, y=y_point_pred)) + geom_point() +
  ggtitle ("y=b0+b1*x1 + b2*x2 vs x2, x1 = x1bar") +
  xlab("x2") + ylab("y") +
  geom_abline(intercept = b0_hat + b2_hat * x2bar, slope = b1_hat) 
ggplot (pyth_na, aes(x=x2, y=y_point_pred)) + geom_point() +
  ggtitle ("y=b0+b1*x1 + b2*x2 vs x2, x1 = x1bar") +
  xlab("x2") + ylab("y") +
  geom_abline(intercept = b0_hat + b1_hat * x1bar, slope = b2_hat)

#
# 11-5
#
# y = 3 + 0.1 * x1 + 0.5 * x2 + error
#
n <- 1000
x_1 <- 1:n
x_2 <- rbinom (n, 1, 0.5)

y = 3 + 0.1 * x_1 + 0.5 * x_2 + 1/5 * rt (n, 4)
#y = 3 + 0.1 * x_1 + 0.5 * x_2 + 1/5 * rnorm (n, 0, 1)

xy <- data.frame(x_1, x_2, y)

x_2_bar <- mean (x_2)

ggplot (xy, aes(x=x_1, y=y)) + geom_point() 
xy <- data.frame(x_1, x_2, y)

fit <- stan_glm (y ~ x_1 + x_2, data = xy)
print (fit, digits = 2)

sims <- as.matrix(fit)
b0_hat <- median(sims[,1])
b1_hat <- median(sims[,2])
b2_hat <- median(sims[,3])
sigma_hat <- median(sims[,4])

y_hat <- b0_hat + b1_hat * xy$x_1 + b2_hat * xy$x_2  

(n - sum (abs (y - y_hat) > sigma_hat))/n

# 
# Mesquite 
#
mesquite <- read.table("mesquite", header=TRUE)

fit_1 <- stan_glm(formula = weight ~ diam1 + diam2 + canopy_height + total_height + density + group,
                  data = mesquite)
print (fit_1)

(loo_1 <- loo (fit_1))
(kfold_1 <- kfold (fit_1, K=10))

# Try log of parameters
fit_2 <- stan_glm(formula = log(weight) ~ log(diam1) + log(diam2) + 
                  log(canopy_height) + log(total_height) + log(density) + group,
                  data = mesquite)
print (fit_2)

(loo_2 <- loo (fit_2))

#
#  REDO: Chapter 11
#

# one predictor - mom_iq
fit_2 <- stan_glm (kid_score ~ mom_iq, data = kidiq)
print (fit_2, digits = 2)

sims <- as.matrix (fit_2)
a_hat <- median (sims[,1])
b_hat <- median (sims[,2])

ggplot (kidiq, aes(x=mom_iq, y=kid_score)) + geom_point() +
  xlab("mom iq") + ylab("child score") +
  geom_abline(intercept = a_hat, slope = b_hat) 

# two predictors - mom_iq, mom_hs
fit_3 <- stan_glm (kid_score ~ mom_hs + mom_iq , data = kidiq)
print (fit_3, digits = 2)

b_hat <- coef (fit_3)

# child_score = b_hat[1] + b_hat[2] * mom_hs + b_hat[3]  * mom_iq
ggplot (kidiq, aes(x=mom_iq, y=kid_score, color=mom_hs)) + geom_point() +
  xlab("mom iq") + ylab("child score") +
  geom_abline(intercept = b_hat[1] + b_hat[2], slope = b_hat[3], color = "green") +
  geom_abline(intercept = b_hat[1], slope = b_hat[3], color = "red") 

# model with interaction 
fit_4 <- stan_glm (kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq , data = kidiq)
print (fit_4, digits = 2)

b_hat <- coef (fit_4)
# child_score = b_hat[1] + b_hat[2]*mom_hs + b_hat[3]*mom_iq + b_hat[4]*mom_hs*mom_iq
ggplot (kidiq, aes(x=mom_iq, y=kid_score, color=mom_hs)) + geom_point() +
  xlab("mom iq") + ylab("child score") +
  geom_abline(intercept = b_hat[1] + b_hat[2], slope = b_hat[3]+ b_hat[4], color = "green") +
  geom_abline(intercept = b_hat[1], slope = b_hat[3], color = "red") 

# simulate n=10 line simulations with fit_2
sims_2 <- as.matrix (fit_2)
n_sims_2 <- nrow (sims_2)
beta_hat_2 <- apply (sims_2, 2, median) 

sims_display <- sample (n_sims_2, 10)

# show uncertainity in betas
ggplot (kidiq, aes(x=mom_iq, y=kid_score)) + geom_point() +
  ggtitle("Child's score v. Mom's IQ - fit with b_hats and 10 simulations") +
  xlab("Mother's iq score") + ylab("Child's  test score") +
  geom_abline(intercept = beta_hat_2[1], slope = beta_hat_2[2]) +
  geom_abline(intercept = sims_2[sims_display,1], slope = sims_2[sims_display,2], 
              size=0.1) 
#
# back to model3
sims_3 <- as.matrix (fit_3)
n_sims_3 <- nrow (sims_3)
beta_hat_3 <- apply (sims_3, 2, median) 
sims_display <- sample (n_sims_2, 10)

mom_hs_bar <- mean(as.numeric(as.character(kidiq$mom_hs)))
# child_score = b_hat[1] + b_hat[2] * mom_hs + b_hat[3]  * mom_iq
ggplot (kidiq, aes(x=mom_iq, y=kid_score)) + geom_point() +
  ggtitle("Child's score v. Mom's IQ - At average Mom High School status") +
  xlab("Mother's iq score") + ylab("Child's  test score") +
  geom_abline(intercept = beta_hat_3[1] + beta_hat_3[2] * mom_hs_bar, 
              slope = beta_hat_3[3]) +
  geom_abline(intercept = sims_3[sims_display,1] + sims_3[sims_display,2] * mom_hs_bar, 
              slope = sims_3[sims_display,3], 
              size=0.1) 

mom_iq_bar <- mean(as.numeric(as.character(kidiq$mom_iq)))

# child_score = b_hat[1] + b_hat[2] * mom_hs + b_hat[3]  * mom_iq
ggplot (kidiq, aes(x=mom_hs, y=kid_score)) + geom_point() + 
  ggtitle("Child's score v. Mom's HS staus - At average Mom IQ") +
  xlab("Mother's iq score") + ylab("Child's  test score") +
  geom_abline(intercept = beta_hat_3[1] + beta_hat_3[3] * mom_iq_bar, 
              slope = beta_hat_3[2]) +
  geom_abline(intercept = sims_3[sims_display,1] + sims_3[sims_display,3] * mom_iq_bar, 
              slope = sims_3[sims_display,2], 
              size=0.1) 
#
# Continuous predictor
# y = a + b*x + theta*z + error  
# z - binary variable
#
N <- 100
x <- runif (N, 0, 1)
z <- sample (c(0,1), N, replace=TRUE)
a     <- 1
b     <- 2
theta <- 5
sigma <- 2

y <- a + b*x + theta*z + rnorm (N, 0, sigma)

fake <- data.frame(x, y, z)
fit <- stan_glm (y~x + z, data = fake)
print (fit)

ggplot (fake, aes(x=x, y=y, color=factor(z))) + geom_point() +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color="red") +
  geom_abline(intercept = coef(fit)[1] + coef(fit)[3], slope = coef(fit)[2], color="green") +
  xlab ("Pre-treatment predictor - x") + ylab ("Outcome - y") +
  ggtitle("y for different levels of binary variable z")

# or 2 searate plots
y_z0 <- coef(fit)[1] +                coef(fit)[2]*x + rnorm (N, 0, sigma)
y_z1 <- coef(fit)[1] + coef(fit)[3] + coef(fit)[2]*x + rnorm (N, 0, sigma)

fake_1 <- cbind (fake, y_z0, y_z1)

p0 <-ggplot (fake_1, aes(x=x, y=y_z0)) + geom_point(shape=19) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2]) +
  xlab ("Pre-treatment predictor - x") + ylab ("Outcome - y") +
  ggtitle("Treament, z = 0") + ylim (-2, 15)
p1 <-ggplot (fake_1, aes(x=x, y=y_z1)) + geom_point(shape=1) +
  geom_abline(intercept = coef(fit)[1] + coef(fit)[3], slope = coef(fit)[2]) +
  xlab ("Pre-treatment predictor - x") + ylab ("Outcome - y") +
  ggtitle("Treatment, z =1") + ylim (-2, 15)

grid.arrange(p0, p1, ncol=2)
#
# Multiple continuous predictors
# y = b0 + b1*x1 + b2*x2 + ... + bk*xk + theta*z + error  
# z - binary variable
#  
N <- 100
K <- 10
X <- array(runif (N*K, 0, 1), c(N,K))
z <- sample (c(0,1), N, replace=TRUE)
a     <- 1
b     <- 1:K
theta <- 5
sigma <- 2

y    <- a + X %*% b + theta*z + rnorm (N, 0, sigma)

fake <- data.frame(X=X, y=y, z=z)

fit <- stan_glm (y ~X + z, data=fake)
print (fit)

y_hat <- predict (fit)

y0 <- y[z==0]
y_hat0 <- y_hat[z==0]
fake_20 <- data.frame (y=y0, y_hat=y_hat0, r=y_hat0-y0)

y1 <- y[z==1]
y_hat1 <- y_hat[z==1]
fake_21 <- data.frame (y=y1, y_hat=y_hat1, r=y_hat1-y1)


p0 <-ggplot (fake_20, aes(x=y_hat, y=y)) + geom_point(shape=19) +
  geom_abline(intercept = 0, slope = 1) +
  xlab ("Linear predictor, y_hat") + ylab ("Outcome - y") +
  ggtitle("Treament, z = 0") + xlim (14,55) + ylim (14, 55)
p1 <-ggplot (fake_21, aes(x=y_hat, y=y1)) + geom_point(shape=1) +
  geom_abline(intercept = 0, slope = 1) +
  xlab ("Linear predictor, y_hat") + ylab ("Outcome - y") +
  ggtitle("Treament, z = 0") + xlim (14,55) + ylim (14, 55)
grid.arrange (p0, p1, ncol=2)

# Residuals
p0 <-ggplot (fake_20, aes(x=y_hat, y=r)) + geom_point(shape=19) +
  xlab ("Linear predictor, y_hat") + ylab ("Residual") +
  ggtitle("Treament, z = 0") + xlim (10,50) + ylim (-14, 23) +
  geom_abline(intercept = 0, slope = 0)
p1 <-ggplot (fake_21, aes(x=y_hat, y=r)) + geom_point(shape=1) +
  xlab ("Linear predictor, y_hat") + ylab ("Residual") +
  ggtitle("Treament, z = 0")+ xlim (10,50) + ylim (-14, 23) +
  geom_abline(intercept = 0, slope = 0)
grid.arrange (p0, p1, ncol=2)

# 
# Introclass
#
introclass <- read.table("introclass", header=TRUE)

fit_1 <- stan_glm (Final ~ Midterm, data = introclass)
print (fit_1)

sims <- as.matrix (fit_1)
y_pred1 <- median (sims[,1]) + median (sims[,2]) * introclass$Midterm
# Or
y_pred2 <- t(median (sims[,1]) + median (sims[,2]) %*% t (introclass$Midterm))
# Or
y_pred3 <- predict (fit_1)

ff1 <- data.frame (x=introclass$Midterm, y_pred1)
ff2 <- cbind (ff1, y_pred2=y_pred2)
ff3 <- cbind (ff2, y_pred3=y_pred3)
p0 <- ggplot (ff1, aes(x=x, y=y_pred1)) + geom_point()
p1 <- ggplot (ff2, aes(x=x, y=y_pred2)) + geom_point()
p2 <- ggplot (ff3, aes(x=x, y=y_pred3)) + geom_point()
grid.arrange(p0, p1, p2, ncol=3)
#
# Plot residuals
#
introclass$y_pred <- y_pred3
introclass$resid1 <- introclass$Final - introclass$y_pred
p0 <- ggplot (introclass, aes(x=y_pred, y=resid1 )) + geom_point() + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Residuals vs. Predicted Values") +
  xlab ("Predicted Value") + ylab ("residual")
p1 <- ggplot (introclass, aes(x=Final, y=resid1 )) + geom_point() + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Residuals vs. Observed Values") +
  xlab ("Observeded Value") + ylab ("residual")
grid.arrange(p0, p1, ncol=2)

# Simulate Fake final data
n <- nrow (introclass)

err <- rnorm (n , 0, median (sims[,3]))
introclass$final_fake1 <- median (sims[,1]) + median (sims[,2]) * introclass$Midterm + err
 # Or
introclass$final_fake2 <- y_pred1 + err
# Or????
pred_final_fake3 <- posterior_predict (fit_1)
apply (pred_final_fake3, 2, sd)
introclass$final_fake3 <- apply (pred_final_fake3, 2, median)

p0 <- ggplot (introclass, aes(x=Midterm, y=final_fake1 )) + geom_point()
p1 <- ggplot (introclass, aes(x=Midterm, y=final_fake2 )) + geom_point()
p2 <- ggplot (introclass, aes(x=Midterm, y=final_fake3 )) + geom_point()
grid.arrange(p0, p1, p2, ncol=3)
#
# Residual Plots
#
# ri <- yi - Xi %*% Beta_hat

# kidiq example
kidiq <- read.dta(file="kidiq.dta")

# One binary predictor: mom's high school status
fit1 = stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit1)

kidiq$kid_score_hat <- coef (fit1)[1] + coef (fit1)[2] * kidiq$mom_iq
kidiq$res <- kidiq$kid_score - kidiq$kid_score_hat
sd (kidiq$res)
ggplot (kidiq, aes(x=kid_score_hat, y=res )) + geom_point() +
  geom_abline(intercept = sd (kidiq$res),  slope = 0, linetype = "dashed") +  
  geom_abline(intercept = 0,               slope = 0) +  
  geom_abline(intercept = -sd (kidiq$res), slope = 0, linetype = "dashed") +
  xlab("Mother IQ Score") + ylab("Residuals") + ylim (-60, 60)
#
# Comparing data to replications from fitted model
#
fit <- stan_glm (y ~ 1, data = newcomb)
print (fit)

sims <- as.matrix (fit)

n_sims <- nrow (sims)
n      <- length (newcomb$y)

y_rep <- array (NA, c(n_sims, n))
for (s in 1:n_sims){
  y_rep[s,] <- rnorm (n, sims[s,1], sims[s,2])
}
# Or ##########
y_rep  <- posterior_predict (fit)  

mean_y_rep <- rep (NA, n)
sd_y_rep   <- rep (NA, n)
for (i in 1:n){
  mean_y_rep[i] <- mean(y_rep[,i])
  sd_y_rep[i]   <- sd(y_rep[,i])  
}
#
# Plot 10 simulations of 60 y_pred each
#

# actural
par(mfrow=c(1,1))
hist (newcomb$y, breaks=seq(-50, 40, 2))
# simulated
par(mfrow=c(2,5))
for (s in sample (n_sims, 10)){
  hist (y_rep[s,])
}

# define test statitic as the minima of each data set
Test <- function (y){
  min(y)
}
test_rep <- apply (y_rep, 1, Test)
hist (test_rep, breaks=20, xlim=range(Test(newcomb$y), test_rep))
abline (v=Test(newcomb$y), lwd = 3)

hist (test_rep, xlim=range(Test(newcomb$y), test_rep))
abline (v=Test(newcomb$y), lwd = 3)

#
# 11-5  Predictive simulation to check fit of time-series model
#
unemp <- read.table("unemp", header=TRUE)
head (unemp)

n <- nrow (unemp)
unemp$y <- unemp$unemployed.pct
unemp$y_lag <- c(NA, unemp$y[1:(n-1)])

fit_lag <- stan_glm (y ~ y_lag, data = unemp)
print (fit_lag, digits=2)                     

sims <- as.matrix(fit_lag)                     
n_sims <- nrow (sims)                     
                     
y_rep <- array (NA, c(n_sims, n))

for (s in 1:n_sims){
  y_rep[s, 1] <- unemp$y[1]
  for (t in 2:n){
    y_rep[s,t] <- sims[s,"(Intercept)"] + sims[s,"y_lag"] * y_rep[s,t-1] *
      rnorm(1,0, sims[s,"sigma"])
  }
}                  

# plot actual data 
ggplot (unemp, aes(x=year, y=y)) + ggtitle("Historical UR") + geom_line() 

# make and plot fake data 
plot_list = list()
j<-1
for (i in sample(n_sims, 15)){
  y_rep_i <- y_rep[i, ]
  fake <- data.frame(y=y_rep_i, x=unemp$year)
  title <- paste ("Sim:", i) 
  plot_list[[j]] <- ggplot (fake, aes(x=x, y=y)) + ggtitle(title) + geom_line() 
  j <- j+1
}

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
          plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], 
          plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], 
          ncol=5, nrow=3)    
#
# Simulations look more jagged than actual data, so test
#
Test <- function (y){
  n <- length (y)
  y_lag   <- c(NA, y[1:(n-1)])
  y_lag_2 <- c(NA, NA, y[1:(n-2)])
  return (sum(sign(y-y_lag) != sign(y_lag-y_lag_2), na.rm = TRUE))
}
 
test_y <- Test (unemp$y)          
test_rep <- apply (y_rep, 1, Test)

hist (test_rep, breaks=seq (7, 45, 2))
abline (v=test_y, lwd = 3)

cat (round(sum (test_rep > test_y)/4000*100, 1), 
     "percent of simultions have more jaggedness than real data.",
     "\nSo model is not satisfactory")

#
# R^2 = 1 - (sigma_hat)^2/(s_y)^2   
# If model fit uses least squares:
##### R^2 = Var(y_hat)/Var(y) #####
# R^2 = 1/(n-1) * sum (y_hat_i - y_hat_bar)^2/s_y^2
#idiq <- read.dta(file="kidiq.dta")
# One predictor: mom's high school status

# Ex
idiq <- read.dta(file="kidiq.dta")
fit1 = stan_glm (kid_score ~ mom_iq, data = kidiq, refresh=0)
print (fit1)
sims <- as.matrix (fit1)

y     <- kidiq$kid_score
y_hat <- median ((sims[,1])) + median ((sims[,2])) * kidiq$mom_iq

n     <- length (y)
y_hat_bar <- mean (y_hat)        # 86.8


s_y       <- sd (y)              # 10.4%
s_y_hat   <- sd (y_hat)          #  9.1%

R2 <- 1 - (sigma_hat)^2/(s_y)^2                     # 19.7%  
# OR - IF LINEAR REGRESSION, i.e. y_hat_i = a_hat + x_i * beta_hat 
R2 <- 1/(n-1) * sum ((y_hat-y_hat_bar)^2)/s_y^2   # 20.0%

# 
# Bayesian R^2
#
##### Bayesian R^2 = Var(y_hat)/[(Var(y_hat) + Var(r)]  #####
# r - y - y_hat (residual)
#
hist(bayes_R2(fit1), breaks =20)
mean(bayes_R2(fit1))

#
# Ex
#
x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3
xy <- data.frame(x,y)

fit <- stan_glm(y ~ x, data = xy)
print (fit)
sims <- as.matrix (fit)

a_hat <- median ((sims[,1]))     # 0.002
b_hat <- median ((sims[,2]))     # 0.59
sigma_hat <- median ((sims[,3])) # 18.3%

y_hat <- a_hat + b_hat * x
r     <- y - y_hat

s_y       <- sd (y)               # 1.08
s_y_hat   <- sd (y_hat)           # 0.93
s_r       <- sd (r)               # 0.52

ggplot (xy, aes(x=x, y=y)) + geom_point() +
  ggtitle("Least Squares and Bayes Fits") +
    geom_abline(intercept = a_hat, slope = b_hat, linetype="solid", size=1)
  
R2 <- var (y_hat)/var(y)                       # 0.74
1 - (sigma_hat/s_y)^2
R2_Bayes <- var (y_hat)/(var(y_hat) + var(r))  # 0.76
#
# Calculate posterior usig informative prior
# prior line : 0 + 1* x
# 
ggplot (xy, aes(x=x, y=y)) + geom_point() +
  ggtitle("Least Squares and Bayes Fits") +
  geom_abline(intercept = a_hat, slope = b_hat, linetype="solid",  size=1) +
  geom_abline(intercept = 0,     slope = 1,     linetype="dotted", size=1) +
  xlim (-2, 2) + ylim (-2, 2)
  
fit_bayes <- stan_glm(y ~ x, data = xy,
                      prior_intercept = normal(0, 0.2, autoscale = FALSE),
                      prior           = normal(1, 0.2, autoscale = FALSE),
                      prior_aux = NULL, refresh = 0)
print (fit_bayes)
sims_bayes <- as.matrix(fit_bayes)
a_hat_bayes <- median (sims_bayes[,1]) # 0.00
b_hat_bayes <- median (sims_bayes[,2]) # 0.82

ggplot (xy, aes(x=x, y=y)) + geom_point() +
  ggtitle("Least Squares and Bayes Fits") +
  geom_abline(intercept = a_hat,           slope = b_hat,            color="blue",  size=1) +
  geom_abline(intercept = 0,               slope = 1,                color="red", size=1) +
  geom_abline(intercept = a_hat_bayes,     slope = b_hat_bayes,      color="green", size=1) +
  annotate("text", x = 1, y = 0.25,    label = "Least-squares fit",  color="blue") +
  annotate("text", x = -0.8, y = -1.3, label = "Posterior mean fit", color="red") +
  annotate("text", x = 1.5, y = 1.2, label = "(Prior regression line)", color="green") +
  xlim (-2, 2) + ylim (-2, 2)

# Draw 20 lines from posterior fit
pick20 <- sample(4000, 20)

ggplot (xy, aes(x=x, y=y)) + geom_point() +
  ggtitle("Bayes Posterior Simulations") +
  geom_abline(intercept = a_hat_bayes,          slope = b_hat_bayes, size=1) +
  geom_abline(intercept = sims_bayes[pick20,1], slope = sims_bayes[pick20, 2], 
              color="darkgray", size = 0.2) +
  xlim (-2, 2) + ylim (-2, 2) 

median(bayes_R2(fit_bayes)) # 80
mean(bayes_R2(fit_bayes))   # 79
sd(bayes_R2(fit_bayes))     # 0.03
hist(bayes_R2(fit_bayes), breaks=40) 

#
# WHICH MODEL IS BETTER?
# 11-8 Cross Validation and Deviance
#
fit3 <- stan_glm(kid_score ~ mom_hs + mom_iq, data=kidiq, refresh = 0)
print(fit3, digits = 2)

sigma_hat3 <- median((as.matrix(fit3, pars="sigma")))  # 18.1

loo3 <- loo (fit3)
print (loo3)   # elpd_loo = -1876
# add five pute noise predictors
n <- nrow (kidiq)
kidiqr <- kidiq  
kidiqr$noise <- array (rnorm (5*n), c(n, 5))  

fit3n <- stan_glm(kid_score ~ mom_hs + mom_iq + noise, data=kidiqr, refresh = 0)
print(fit3n, digits = 2)

sigma_hat3n <- median(as.matrix(fit3n)[,"sigma"])  # 18.2
loo3n <- loo (fit3n)
print (loo3n)   # elpd_loo = -1881
# adding 5 noise terms makes sigma a bit worse from 18.1 t0 18.2, 
# Posterior Predictive log store also worse from -1876 to -1881 with addition of 5 parameters

# simpler model 
fit1 <- stan_glm(kid_score ~ mom_hs, data=kidiq, refresh = 0)
print(fit1, digits = 2)
sigma_hat1 <- median(as.matrix(fit1)[,"sigma"])  # 19.9 - a bit wors

loo1 <- loo (fit1)
print (loo1) # elpd -1915, a decline of 39 from model3 - not good

compare_models (loo3, loo1)  # elpd diff -39, se 8.3

# go to more sophisticated model
fit4 <- stan_glm(kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq, data=kidiq, refresh = 0)
print(fit4, digits = 2)

sigma_hat4 <- median((as.matrix(fit4, pars="sigma")))  # 18.0

loo4 <- loo (fit4)
print (loo4)                # elpd 1873

compare_models (loo3, loo4)  # elpd diff 3.4, se 2.6 - pretty indistinguishable from zero

# compute Bayesian R2
median(bayes_R2(fit1))  # 0.056
median(bayes_R2(fit3))  # 0.21
median(bayes_R2(fit3n)) # 0.22
median(bayes_R2(fit4))  # 0.23
# compute loo R2
round(looR2(fit1),3)    # 0.047  mom_hs
round(looR2(fit3),3)    # 0.203  mom_hs + mom_iq
round(looR2(fit3n),3)   # 0.183  added noise
round(looR2(fit4),3)    # 0.187  added mom_hs:mom_iq interaction

#
# K-fold cross validation
#
n   <- 60
k   <- 30
rho <- 0.8
Sigma <- rho * array(1, c(k,k)) + (1-rho)* diag(k)
X     <- mvrnorm (n, rep(0,k), Sigma)

b <- c(c(-1,1,2), rep(0, k-3))
y <- X %*% b + 2*rnorm(n)

fake <- data.frame(X,y)

# Use model with WEAK prior
fit_1 <- stan_glm (y ~ ., prior=normal(0, 10, autoscale = FALSE), data = fake)
loo_1 <- loo (fit_1)
print (loo_1)  # get nasty messages, use kfold

# Try kfold approach
kfold_1 <- kfold (fit_1, K=10)
print (kfold_1)

# Now use model with more INFORMTIVE prior, horseshoe, hs_prior

k0 <- 2 # prior guess for the number of relevant variables
tau0 <- k0/(k-k0) * 1/sqrt(n)
hs_prior <- hs(df=1, global_df=1, global_scale=tau0, slab_scale=3, slab_df=7)

fit_2 <- stan_glm(y ~ ., prior=hs_prior, data=fake)
print (fit_2)

loo_2 <- loo (fit_2)  
print (loo_2)  # get nasty messages, use kfold

# Try kfold approach
kfold_2 <- kfold (fit_2, K=10)
print (kfold_2)

compare_models(kfold_1, kfold_2)  # elpd_diff=14.9, se=6.1  => horseshoe prior is better    

#
# CH 12
#
earnings <- read.csv('earnings.csv',header=T)

# Eliminate earnings = 0
earnings <- earnings[earnings$earn != 0, ]
earnings$earn <- earnings$earn/1000
earnings$log_earn <- log (earnings$earn)

# simple
fit <- stan_glm (earn ~ height, data = earnings)
print (fit)
b <- coef (fit)  # b[1] = median (as.matrix (fit, pars="(Intercept)"))

ggplot (earnings, aes(x=height, y=earn)) + geom_point() +
  geom_abline(intercept = b[1], slope = b[2])

# Replicte data for linear fit
yrep1 <- posterior_predict(fit)
n_sims <- nrow (yrep1)

subset <- sample (n_sims, 100)
ppc_dens_overlay (earnings$earn, yrep1[subset,])
# a bunch of replicated negative earnings even tho data has none!

# log fit
fit2 <- stan_glm (log_earn ~ height, data = earnings)
print (fit2)
b <- coef (fit2)  
sims <- as.matrix (fit2)
pick_sims <- sample (4000, 20)

ggplot (earnings, aes(x=height, y=log_earn)) + geom_point() +
  geom_abline(intercept = b[1], slope = b[2]) +
  geom_abline(intercept = sims[pick_sims,1], slope = sims[pick_sims,2], size=0.1) +
  xlab ("Height") + ylab("Log (earnings - in K)")

# Replicte data for log fit
yrep2 <- posterior_predict(fit2)
ppc_dens_overlay (earnings$log_earn, yrep2[subset,])
br2_2 <- median(bayes_R2(fit2)) # 0.06

#
# log fit add male
fit3 <- stan_glm (log_earn ~ height + male, data = earnings)
print (fit3)
b <- coef (fit3)  
sims <- as.matrix (fit3)

exp (b[1] + b[2] * 70 + b[3] * 0)

ggplot (earnings, aes(x=height, y=log_earn)) + geom_point() +
  geom_abline(intercept = b[1], slope = b[2]) +
  geom_abline(intercept = sims[pick_sims,1], slope = sims[pick_sims,2], size=0.1) +
  xlab ("Height") + ylab("Log (earnings - in K)")

# Replicte data for log fit + male
yrep3 <- posterior_predict(fit3)
ppc_dens_overlay (earnings$log_earn, yrep3[subset,])
br2_3 <- median(bayes_R2(fit3)) # 0.09

#
#
new <- data.frame("height" = 70)
y_point_pred <- exp(predict (fit2, newdata = new))            # 19.8k
y_lin_pred   <- exp(posterior_linpred (fit2, newdata = new))  # 19.8k, se=640
y_pred       <- exp(posterior_predict (fit2, newdata = new) )

# include interaction
fit4 <- stan_glm(log_earn ~ height + male + height:male, data = earnings)
print (fit4, digits=3)

yrep4 <- posterior_predict(fit4)
ppc_dens_overlay (earnings$log_earn, yrep4[subset,])
br2_4 <- median(bayes_R2(fit4)) # 0.09

# log earning in thousands = 1.4 + 0.02*height + 0.07*male + 0.005* height * male

#
# use z-score for height
#
earnings$z_height <- (earnings$height - mean (earnings$height))/sd (earnings$height)

fit5 <- stan_glm(log_earn ~ z_height + male + z_height:male, data = earnings)
print (fit5, digits=3)

yrep5 <- posterior_predict(fit5)
ppc_dens_overlay (earnings$log_earn, yrep5[subset,])
br2_5 <- median(bayes_R2(fit5)) # 0.09

# Log-Log transforms
earnings$log_height <- log (earnings$height)
  
fit6 <- stan_glm(log_earn ~ log_height + male, data = earnings)
print (fit6, digits=3)
b <- coef (fit6)

# for females, male=0
ggplot (earnings, aes(x=log_height, y=log_earn)) + geom_point() +
  geom_abline(intercept =b[1], slope =b[2])
  
yrep6 <- posterior_predict(fit6)
ppc_dens_overlay (earnings$log_earn, yrep6[subset,])
br2_6 <- median(bayes_R2(fit6)) # 0.09

#
# Mesquite Bushes
#
mesquite <- read.table("mesquite", header=TRUE)

# all predictors
fit_1 <- stan_glm (weight ~ diam1 + diam2 + canopy_height + 
                     total_height + density + group, data = mesquite)
print (fit_1)
yrep_1 <- posterior_predict(fit_1)
ppc_dens_overlay (mesquite$weight, yrep_1[subset,])
br2_1 <- median(bayes_R2(fit_1)) # 0.83
# evaluate fit
loo_1 <- loo (fit_1)  # elpd = -334.5    not stable, use K-fold

kfold_1 <- kfold (fit_1, K=10)  # elpd = -346, se 25
kfold_1

# logarithmic scale
fit_2 <- stan_glm (log(weight) ~ log(diam1) + log(diam2) + log(canopy_height) + 
                     log(total_height) + log(density) + group, data = mesquite)
print (fit_2)

yrep_2<- posterior_predict(fit_2)
ppc_dens_overlay (log(mesquite$weight), yrep_2[subset,])

br2_2 <- median(bayes_R2(fit_2)) # 0.87
loo_2 <- loo (fit_2) # elpd = -19.2, se 5.3
loo_2

# simpler model
mesquite$canopy_volume <- mesquite$diam1 * mesquite$diam2 * mesquite$canopy_height

fit_3 <- stan_glm (log(weight) ~ log(canopy_volume), data = mesquite)
print (fit_3)

yrep_3<- posterior_predict(fit_3)
ppc_dens_overlay (log(mesquite$weight), yrep_3[subset,])

br2_3 <- median(bayes_R2(fit_3)) # 0.80
loo_3 <- loo (fit_3) # elpd = -26.7, se 5.1
loo_3

compare_models(loo_2, loo_3)  # elpd diff=-7.5, se=5  model_2 better

# iterate model
mesquite$canopy_area  <- mesquite$diam1 * mesquite$diam2 
mesquite$canopy_shape <- mesquite$diam1 / mesquite$diam2 

fit_4 <- stan_glm (log(weight) ~ log(canopy_volume) + log(canopy_area) +
                     log(canopy_shape) + log(total_height) + log(density) + group, data = mesquite)
print (fit_4)

yrep_4<- posterior_predict(fit_4)
ppc_dens_overlay (log(mesquite$weight), yrep_4[subset,])

br2_4 <- median(bayes_R2(fit_4)) # 0.88
loo_4 <- loo (fit_4) # elpd = -19.6, se 5.4, some error messages
loo_4

compare_models(loo_2, loo_4)  # elpd diff=-0.3, se=0.1   ~identical

# levae predictors in that have small se's
fit_5 <- stan_glm (log(weight) ~ log(canopy_volume) + log(canopy_shape) + group, data = mesquite)
print (fit_5)

yrep_5<- posterior_predict(fit_5)
ppc_dens_overlay (log(mesquite$weight), yrep_5[subset,])

br2_5 <- median(bayes_R2(fit_5)) # 0.87
loo_5 <- loo (fit_5) # elpd = -18.2, se 5.4
loo_5

compare_models(loo_4, loo_5)  # elpd diff=1.5, se=1.5   ~identical but fewer predictive prameters.
                              # Winner!
#
# Chapter 13 - Logistic Regression
#
# logit (x)    = log (x/(1-x))
# invlogit (x) = e^x/(1+e^x)
# invlogit (-inf) = 0
# invlogit (0)    = 1/2
# invlogit (inf)  = 1

# Pr (yi = 1) = invlogit (Xi * beta)
# for linear model, a + b * x, slope = d(invlogit(x))/dx
# = b * [((1+e^x)e^x - e^x(e^x))/(1+e^x)^2]
# = b * e^x/(1+e^x)^2
# slope maximum at x=0,  b * e^0/(1+e^0)^2 
# max slope = b/4
attach("/Users/nevinaltaras/Downloads/nes.rda")
head (nes)

nes92 <- subset (nes, year== 1992)
#nes92$income_jitt <- nes92$income + runif (length(nes92$income), -0.25, 0.25)

fit_1 <- stan_glm (rvote ~ income, family = binomial(link = "logit"), data = nes92)
print (fit_1, digits = 2)
sims <- as.matrix(fit_1)

pr_vote <- invlogit (coef(fit_1)[1] + coef(fit_1)[2] * nes92$income)

sample_20 <- sample(4000,20)
pr_vote_s <- array(NA, c(20,nrow(nes92)))
k <- 1
for (i in sample_20){
  for(j in 1:nrow (nes92)){
    pr_vote_s[k,j] <- invlogit (sims[i,1]+sims[i,2]*nes92$income[j])
  }
  k<-k+1
}

pr_vote_s <- data.frame(t(pr_vote_s), income=nes92$income, rvote=nes92$rvote, pr_vote)

pr_vote_s_melt <- melt (pr_vote_s, id=c("income", "rvote", "pr_vote"))  #******#

ggplot (pr_vote_s_melt) + 
  geom_line   (aes(x=income, y=value, color = variable)) +
  geom_jitter (aes(x=income, y=rvote), width = 0.25, height = 0.05) +
  geom_line(aes(x=income, y=pr_vote))

# Predictions
# what is the probability that a voter in income category 5 will vote for Bush?
new <- data.frame(income=5)
linpred <- posterior_linpred(fit_1, transform = TRUE, newdata = new)
cat ("Prob income=5 voting for Bush=", median (linpred), "(", sd(linpred), ")")

#
# Wells in Bangladesh
wells <- read.table ("arsenic", header = T, sep = ",")
wells$y <- wells$switch
wells <- na.omit (wells)

fit_1 <- stan_glm (y ~ dist100, family = binomial(link = "logit"), data = wells)
print (fit_1, digits = 2)
(loo1 <- loo(fit_1)) 

hist (wells$dist, breaks = seq(0, 340, 5))
mean (wells$dist) # 48.3

wells$fit <- invlogit (coef(fit_1)[1] + coef(fit_1)[2] * wells$dist100)
wells$fitted <- fit_1$fitted

ggplot (wells) + 
  geom_jitter (aes(x=dist, y=y), height = 0.01) +
  geom_line (aes(x=dist, y=fit)) +
  geom_line (aes(x=dist, y=fitted), color="red") +
  ggtitle("Switch Probability vs Safe Well Distance") +
  xlab ("Distance to nearest safe well") + ylab ("Probability of switching")

# what is a) the Pr(switch) at average distance?
#         b) the sensitivity of same

# a)
mean (wells$dist) # 48
pr_switch_mean_dist  <-invlogit (coef(fit_1)[1] + coef(fit_1)[2] * mean(wells$dist)) # 60%
# or
new_data <- data.frame(dist = mean (wells$dist))
pr_switch_mean_dist <- invlogit (predict(fit_1, newdata = new_data)) 
# b)  invlogit (x) = e^x/(1+e^x) where x=a+b*dist
# take derivative  ((1+e^x)*e^x- e^x*e^x)/(1+e^x)^2 *dx/d dist
#  = b * e^x/(1+e^x)^2
coef(fit_1)[2]* exp(pr_switch_mean_dist)/(1+exp(pr_switch_mean_dist))^2
# -0.14/meter..  at  meters,  if the well is 1 meter 
# farther(nearer), ie at 49(47) meters, the prob of switching falls(rises) by 0.14%

# OR, rule of four:  
coef(fit_1)[2]/4 # = =-0.16/meter

#
# Add another predictor, arsenic
hist(wells$arsenic, breaks=seq(0,.25+max(wells$arsenic),.25), freq=TRUE,
     xlab="Arsenic concentration in well water", ylab="", main="", mgp=c(2,.5,0))

fit_2 <- stan_glm (y ~ dist100 + arsenic, family = binomial(link = "logit"), data = wells)
print (fit_2, digits = 2) # dist100:-.90(0.11), arsenic:0.46(0.04)
(loo2 <- loo(fit_2))

compare_models(loo2, loo1) # elpd diff= -72(12) favoring fit_2
# Rule of 4: 
# For every 100m distance in a well, prob of switching goes down by .90/4=22%
# For every 1 in arsenic concentration, prob of switching goes down up .46/4=11%
#
# Plot
# arsenic = {0.5, 1.0}
as_1.0 <- invlogit (coef(fit_2)[1] + 
                    coef(fit_2)[2] * wells$dist100 + 
                    coef(fit_2)[3] * 0.5)
as_0.5 <- invlogit (coef(fit_2)[1] + 
                    coef(fit_2)[2] * wells$dist100 + 
                    coef(fit_2)[3] * 1.0)

fake <- data.frame (y=wells$switch, dist=wells$dist100, as_0.5, as_1.0)

ggplot (fake) + 
  geom_jitter (aes(x=dist, y=y), height = 0.01) +
  geom_line (aes(x=dist, y=as_0.5)) +
  geom_line (aes(x=dist, y=as_1.0), color="red") +
  annotate("text", x = 0.80, y = 0.30, label = "As=0.5") +
  annotate("text", x = 1.2,  y = 0.45, label = "As=1.0", color="red") +
  ggtitle("Switch Probability vs Safe Well Distance") +
  xlab ("Distance to nearest safe well") + ylab ("Probability of switching")

# distance = {0, 50}
dist_0  <- invlogit (coef(fit_2)[1] + 
#                     coef(fit_2)[2] * wells$dist100 + 
                     coef(fit_2)[3] * wells$arsenic)
dist_50 <- invlogit (coef(fit_2)[1] + 
                      coef(fit_2)[2] * 0.5 + 
                      coef(fit_2)[3] * wells$arsenic)

fake <- data.frame (y=wells$switch, arsenic=wells$arsenic, dist_0, dist_50)

ggplot (fake) + 
  geom_jitter (aes(x=arsenic, y=y), height = 0.01) +
  geom_line (aes(x=arsenic, y=dist_0)) +
  geom_line (aes(x=arsenic, y=dist_50), color="red") +
  annotate("text", x = 1.0, y = 0.70, label = "Dist=0") +
  annotate("text", x = 4.1, y = 0.75, label = "Dist=50", color="red") +
  ggtitle("Switch Probability v Arsenic Concentration") +
  xlab ("Arsenic Concentration in Well") + ylab ("Probability of switching")

#
# Chapter 14: Working with Logistic Regression
#
n <- 50
a <- 2
b <- 3
x_mean <- -a/b
x_sd <- 4/b
x <- rnorm(n, x_mean, x_sd)
y <- rbinom(n, 1, invlogit(a + b*x))
fake_1 <- data.frame(x, y)

fit <- stan_glm(y ~ x, family=binomial(link="logit"), data=fake_1)
print (fit, digits=3)

orig_y <- invlogit (a + b*x)
fake_2 <- cbind (fake_1, orig_y, fitted_y=fit$fitted)

ggplot (fake_2) + geom_point (aes(x=x, y=y)) +
  geom_line (aes(x=x, y=orig_y)) +
  geom_line (aes(x=x, y=fitted_y), linetype="dashed") +
  annotate("text", x = -1, y = 0.70, label = "Original Curve") +
  annotate("text", x = 0.5, y = 0.80, label = "Fitted Curve") 

#Binning
n_bins <- 5
bins <- as.numeric(cut (x, n_bins))

x_bar <- rep (NA, n_bins)
y_bar <- rep (NA, n_bins)
for (i in 1:n_bins){
  x_bar[i] <- mean (x[bins==i])
  y_bar[i] <- mean (y[bins==i])
}

x_bins <- x_bar[bins]
y_bins <- y_bar[bins]
fake_3 <- cbind (fake_2, x_bins, y_bins)

ggplot (fake_3) + geom_point (aes(x=x, y=y)) +
  geom_line (aes(x=x, y=orig_y)) +
  geom_line (aes(x=x, y=fitted_y), linetype="dashed") +
  annotate("text", x = -1, y = 0.70, label = "Original Curve") +
  annotate("text", x = 0.5, y = 0.80, label = "Fitted Curve") +
  geom_point (aes(x=x_bins, y=y_bins), color='red', size=3) 

# Logistic regression with interactions
#
fit_4 <- stan_glm (y ~ dist100 + arsenic + dist100:arsenic, family = binomial(link = "logit"), data = wells)
print (fit_4, digits = 2) 
(loo4 <- loo(fit_4))
loo_compare (loo4, loo2)
# Use centering
wells$c_dist100 <- wells$dist100 - mean (wells$dist100)
wells$c_arsenic <- wells$arsenic - mean (wells$arsenic)

fit_5 <- stan_glm (y ~ c_dist100 + c_arsenic + c_dist100:c_arsenic, family = binomial(link = "logit"), data = wells)
print (fit_5, digits = 2) 
(loo5 <- loo(fit_5))
loo_compare (loo5, loo2)
#(Intercept)          0.35 -invlogit(0.35)=59% Pr(switch) at avg dist and avg arsenic
#c_dist100           -0.87 -Pr(switch) decreases 22% for 100 more dist for good well
#c_arsenic            0.47 -Pr(switch) increases 11% for 1 more unit As at curr well
#c_dist100:c_arsenic -0.18
invlogit (0.35) # 59% prob of switching at avg dist and avg arsenic

# Predictive Simulation Using Binomial Dsitribution
fit_1 <- stan_glm (y ~ dist100, family = binomial(link = "logit"), data = wells)
print (fit_1, digits = 2)

sims <- as.matrix(fit_1)
n_sims <- nrow (sims)

dists <- c (0.2, 0.5, 0.8, 1.3, 1.8, 2.2, 2.3, 3.1, 3.3, 3.4)
n_new <- length (dists)

dist_new <- data.frame (dist100 = dists)

probs <- posterior_linpred (fit_1, transform = TRUE, newdata = dist_new)

y_new <- array (NA, c(n_sims, n_new))
for (s in 1:n_sims){
  y_new[s, ] <- rbinom (n_new, 1, prob=probs[s,])
}
print (colMeans (y_new), digits=2)

# OR
one <- rep (1, n_new)
X_new <- cbind (one, dist_new)

y_new_p <- array (NA, c(n_sims, n_new))
for (s in 1:n_sims){
  pr <- invlogit (X_new %*% sims[s,])
  y_new_p[s, ] <- rbinom (n_new, 1, prob=pr)
}
print (colMeans (y_new_p), digits=2)

# Average Predictive comparisons
wells$educ4 <- wells$educ/4
wells$c_educ4 <- wells$educ4 - mean (wells$educ4)
fit_7 <- stan_glm (y ~ dist100 + arsenic + educ4, 
                   family = binomial(link = "logit"), data = wells)
print (fit_7, digits=2)
b <- coef (fit_7)
# Ps (switch) = invlogit (b[1] + b[2]*dist100 + b[3]*arsenic + b[4]*educ4)

# Vary dist100: Let
dist_lo <- 0
dist_hi <- 1

delta <- invlogit (b[1] + b[2]*dist_hi + b[3]*wells$arsenic + b[4]*wells$educ4) - 
         invlogit (b[1] + b[2]*dist_lo + b[3]*wells$arsenic + b[4]*wells$educ4)

mean (delta)  # -20% lower prob of switching if house is 100m farhter 
              # for same arsenic, educ level

# Vary arsenic: Let
arsenic_lo <- 0.5
arsenic_hi <- 1.0

delta <- invlogit (b[1] + b[2]*wells$dist100 + b[3]*arsenic_hi + b[4]*wells$educ4) - 
         invlogit (b[1] + b[2]*wells$dist100 + b[3]*arsenic_lo + b[4]*wells$educ4)

mean (delta) # 6% higher prob of switching if well has 0.5 higher arsenic 
             # for same distance, educ level

# Average predictive comparisons in the presence of interavctions
fit_8 <- stan_glm (y ~ c_dist100 + c_arsenic + c_educ4 + 
                   c_dist100:c_educ4 + c_arsenic:c_educ4,
                   family = binomial(link = "logit"), data = wells)
print (fit_8, digits=2)

# Vary dist100: Let
dist_lo <- 0
dist_hi <- 1

b <- coef (fit_8)

delta <- invlogit (b[1] + b[2]*dist_hi + b[3]*wells$c_arsenic + b[4]*wells$c_educ4 +
                   b[5]*wells$c_dist100*wells$c_educ4 + b[6]*wells$c_arsenic*wells$c_educ4) - 
         invlogit (b[1] + b[2]*dist_lo + b[3]*wells$c_arsenic + b[4]*wells$c_educ4 +
                   b[5]*wells$c_dist100*wells$c_educ4 + b[6]*wells$c_arsenic*wells$c_educ4)  
mean (delta) # -21%

#
# Residuals for Logistic Regression- ***BINNING***
# ri = yi - E (yi | Xi) = yi - invlogit (Xi * beta)

wells$fitted <- fit_8$fitted
r <- wells$switch - wells$fitted
fake <- data.frame(estimated = wells$fitted, r)

ggplot (fake) + geom_point (aes(x=estimated, y=r)) +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Residuals") +
  xlab ("Estimated Pr(switching)") + ylab ("Observed - Estimated")

# Not helpful.  Do binned residuals
n_bins <- 40 
bins <- as.numeric(cut_number (fake$estimated, n_bins))

estimated_bar <- rep (NA, n_bins)
r_bar         <- rep (NA, n_bins)
se_bin        <- rep (NA, n_bins)
for (i in 1:n_bins){
  estimated_bar[i] <- mean (fake$estimated[bins==i])
  r_bar[i]         <- mean (fake$r[bins==i])
  se_bin[i]        <- sd (fake$r[bins==i])/sqrt(sum (bins==i))
}

fake1 <- data.frame(estimated_bar, r_bar, se_bin)

ggplot (fake1) + geom_point (aes(x=estimated_bar, y=r_bar)) +
  geom_line (aes(x=estimated_bar, y=-2*se_bin)) +
  geom_line (aes(x=estimated_bar, y= 2*se_bin)) +
  ggtitle("Binned Residual Plot") +
  xlab ("Estimated probability of switching") + ylab ("Average residual")

# OR
binnedplot(fake$estimated, fake$r, nclass = 40)

#
# Bin residuals for only arsenic as predictor
#
wells$fitted <- fit_8$fitted
r <- wells$switch - wells$fitted
fake <- data.frame(arsenic = wells$arsenic, r)

n_bins <- 40

bins <- as.numeric(cut_number (wells$arsenic, n_bins))

r_bar         <- rep (NA, n_bins)
as_bar        <- rep (NA, n_bins)
se_bin        <- rep (NA, n_bins)
for (i in 1:n_bins){
  r_bar[i]    <- mean (fake$r[bins==i])
  as_bar[i]   <- mean (fake$arsenic[bins==i])
  se_bin[i]   <- sd (fake$r[bins==i])/sqrt(sum (bins==i))
}

fake1 <- data.frame(as_bar, r_bar, se_bin)

ggplot (fake1) + geom_point (aes(x=as_bar, y=r_bar)) +
  geom_line (aes(x=as_bar, y=-2*se_bin), color="gray") +
  geom_line (aes(x=as_bar, y= 2*se_bin), color="gray") +
  ggtitle("Binned Residual Plot") +
  xlab ("Arsenic Level") + ylab ("Average residual") 

# OR
binnedplot(fake$arsenic, fake$r, nclass=40)
# Transform arsenic as predictor to log(arsenic)
#
wells$c_log_arsenic <- log (wells$arsenic) - mean(log (wells$arsenic))
fit_9 <- stan_glm (y ~ c_dist100 + c_log_arsenic + c_educ4 + 
                     c_dist100:c_educ4 + c_log_arsenic:c_educ4,
                   family = binomial(link = "logit"), data = wells)
print (fit_9, digits=2)

(loo8 <- loo(fit_8))
(loo9 <- loo(fit_9))
loo_compare (loo9, loo8) # 14.7(4.3)  log transforming arsenic improves predictive power


#
#
nes64 <- subset (nes, year== 1964)
fit_glm <- glm (rvote ~ female + black + income, family = binomial(link = "logit"), data = nes64)
display(fit_glm) # black: -16 (420)
coefplot (fit_glm)
summary(fit_glm)$coefficients[, 1]  # coefs
summary(fit_glm)$coefficients[, 2]  # se's

fit_stan <- stan_glm (rvote ~ female + black + income, family = binomial(link = "logit"), data = nes64)
print (fit_stan, digits = 2) # black: -4.5 (1.1)

x <- seq (-15, 5, 0.1)
likely <- (invlogit(-1.2 - 4.5*x))
black <- data.frame (x, likely)

ggplot (black) + geom_line (aes(x=x, y=likely))
plot (x, invlogit(-1.2 - 4.5*x))
#
# Ex 5-1
# 
N <- 50
x <- runif (N, -10, 15)
y <- invlogit (0.4 - 0.3*x)
y_yn <- rbinom (N, 1, y)

fake <- data.frame (x, y, y_yn)
ggplot (fake) + geom_jitter (aes(x=x, y=y_yn), size=0.7, height = 0.01) +
  geom_line (aes(x=x, y=y)) +
  ggtitle ("X yes or no") +
  xlab ("Runif (50, -10, 10)") + ylab ("Pr(x=1)")
#
# Ex 15-2
#
N <- 100
x1 <- runif (N, -20, 0)
x2 <- runif (N, -20, 0)

mean (x1)
mean (x2)
pr_true <- invlogit (0.4 - 0.3*x1 + 0.2*x2)
y <- rbinom (N, 1, pr_true)

plot (x1, y)
plot (x2, y)

fake <- data.frame(x1, x2, y)

fit <- stan_glm (y ~ x1 + x2, family = binomial(link = "logit"), data = fake)
print (fit, digits = 2)
beta <- coef (fit)

# pr(y=1) = 0.5 when  beta[1] + beta[2]*x1 + beta[3]*x2 = logit(0.5) = 0
# => x2 = (- beta[1])/beta[3] - beta[2]/beta[3] * x1

# pr(y=1) = 0.1 when  beta[1] + beta[2]*x1 + beta[3]*x2 = logit(0.1) = 2.2
# => x2 = (2.2 - beta[1])/beta[3] - beta[2]/beta[3] * x1
intercept50 <- (- beta[1])/beta[3]
slope50     <- -beta[2]/beta[3]
intercept10 <- (logit(0.1)- beta[1])/beta[3]
slope10     <- -beta[2]/beta[3]
intercept90 <- (logit(0.9)- beta[1])/beta[3]
slope90     <- -beta[2]/beta[3]
ggplot (fake) + geom_point (aes(x=x1, y=x2, color=factor(y))) +
  geom_abline(intercept=intercept10, slope=slope10, color='red') +
  geom_abline(intercept=intercept50, slope=slope50) +
  geom_abline(intercept=intercept90, slope=slope90, color='blue') +
  annotate("text", x = -12, y = -7,  label = "90%", color="blue") +
  annotate("text", x = -10, y = -12, label = "50%", color="black") +
  annotate("text", x = -5, y = -17,  label = "10%", color="red") +
  ggtitle ("Data and 10%, 50%, 90% discrimination lines\nfrom fitted logistic regression") +
  scale_color_manual(name = "Pr(y=1)", labels = c("0", "1"), values = c("red", "blue"))  
#
# Ex 15-3
#
wells$log_dist <- log (wells$dist)
fit_15_3 <- stan_glm (y ~ log_dist, family = binomial(link = "logit"), data = wells)
print (fit_15_3, digits = 2)

beta <- coef (fit_15_3)
X <- cbind (1, wells$log_dist)
prob_swp <- invlogit(X %*% beta)
prob_sw  <- invlogit (beta[1] + beta[2]*wells$log_dist)

fake <- data.frame (log_dist = wells$log_dist, y=wells$y, prob_sw, prob_swp)

ggplot (fake) + geom_jitter (aes(x=log_dist, y=y), size=0.7, height = 0.01) +
  geom_line (aes(x=log_dist, y=prob_sw)) +
  geom_line (aes(x=log_dist, y=prob_swp), color="red") +
  ggtitle ("Switch Well: yes or no") +
  xlab ("Log (distance from well)") + ylab ("Pr(switch well")

resid <- wells$y - prob_sw

binnedplot(prob_sw, resid, nclass = 40)
# 
# Ex 15-5
# Transform  x=N (mu, sigma) to unit_norm=N (0, 1)
# unit_n0rm = (x - mu)/sigma
# x = unit_norm*sigma + mu

N<-50
set.seed (1234)
midterm  <- (rnorm (N, 60, 15))
set.seed (1234)
unit_norm <- rnorm (N, 0, 1)

pr_pass  <- invlogit (-24 + 0.4 * midterm)

midterm_p <- unit_norm*15+60
pr_pass_p  <- invlogit (-24 + 0.4 * midterm_p)
pass_yn <- rbinom (N, 1, y)

fake <- data.frame(midterm, midterm_p, pr_pass, pr_pass_p, pass_yn)

ggplot (fake) + geom_line (aes(x=midterm, y=pr_pass)) + 
  geom_line (aes(x=midterm_p, y=pr_pass_p), color='red') + 
  xlab ("Midterm Grade") + ylab ("Probability of Passing Course")

ggplot (fake) + geom_jitter (aes(x=midterm, y=pass_yn), size=0.7, height = 0.01) +
  geom_line (aes(x=midterm, y=pr_pass)) + 
  ggtitle ("Pass/Fail vs Midterm Grade") +
  xlab ("Midterm Grade") + ylab ("Pr(Pass Course)")

fit <- stan_glm (pass_yn ~ midterm, family = binomial(link = "logit"), data = fake)
print (fit)

set.seed (1234)
midterm_diffuse  <- rnorm (N, 60, 15) + rnorm (N, 0, 1)
fake_diff <- cbind (fake, midterm_diffuse)

fit_diff <- stan_glm (pass_yn ~ midterm_diffuse, family = binomial(link = "logit"), data = fake_diff)
print (fit_diff)

(loo_fit      <- loo(fit))
(loo_fit_diff <- loo(fit_diff))
loo_compare (loo_fit, loo_fit_diff)

# Chapter 16 - Design
#

# Sample size to achieve a desired standard error, se

# Views on death penaly

# n -  number in sample
p_hat <- 0.6 # estimate of death penalty support in population 60%

se <- sqrt (p_hat * (1-p_hat)/n) # standard error of sample

# need se          < 0.05, 5%
# sqrt (p_hat * (1-p_hat)/n) < 0.05
# sqrt (p_hat * (1-p_hat)) / sqrt(n)   < 0.05
# n > p_hat * (1-p_hat)/(0.05)^2
n_min <- p_hat*(1-p_hat)/(0.05)^2  # 96

# if we dont have estimate for p_hat, conservatively use p_hat = 0.5
# n_min <- (sqrt (0.5*0.5)/0.05)^2           # 100
# n_min <- (0.5/0.05)^2

# want to prove p > 0.5, i.e more than half favor death penalty
# Estimate p_hat = 0.5
# se (p_hat) = 0.5/sqrt(n)  [ sqrt (0.5*(1-0.5)/n) ]  *********
# 95% CI for p is [p_hat +/- 1.96 * 0.5/sqrt(n)] 
# we can say p > 0.5 if whole CI lies above 0.5, i.e lower bound of CI > 0.5
# p_hat - 1.96 * 0.5 / sqrt(n) > 0.5, i.e
# p_hat > 0.5  + 1.96 * 0.5 / sqrt (n)

# POWER: The probability that the 95% CI will be entirely above the 0.5 comparison point
# if p = 0.6 
# To find the value of n such that exactly 80% of the estimates will be at least 1.96 standard errors
# from 0.5, we need:
# 0.5 + 1.96 * s.e. = 0.6 - 0.84 * se    [ pnorm (0.8) = 0.84 ], [ se = 0.5/sqrt(n)]
# 0.1 = 2.8 * se
# 0.1 = 2.8 * 0.5/sqrt(n)
# n = 196

# IN GENERAL,       n = p*(1-p)*(2.8/(p-p_hat))^2
# and if p ~ 0.5,   n = (0.5*2.8/(p-p_hat))^2       


# To have 80% power, the true value of the paremeter must be 2.8 std away from the comparison point.
# 1.96 std for the 95% interval, and additional 0.84 to reach the 80th 
# percentile of the normal distribution

# Comparisons of Proportions - standard error
# se = sqrt (p1*(1-p1)/n1 + p2*(1-p2)/n2)
# upper bound when p1=p2=0.5, se = 0.5 * sqrt (1/n1 + 1/n2)
# if n1 = n2 = n/2,           se = 1/sqrt(n)   ************
# A specified se can be attained by a sample size of n > 1/(se)^2

# Ex: suspect death penalty 10% more popular in US than in Canada.  
# Will sample n/2 people in each country
# How large should n be to achieve statistical significance?
# p1_hat - favorability in US
# p2_hat - favorability in Canada
# se of p1_hat - p2_hat = sqrt (p1_hat*(1-p1_hat)/n/2 +p2_hat*(1-p2_hat/n/2))
# assume p1 ~ p2 ~ 0.5
# se = sqrt (0.5*0.5*2/n+0.5*0.5/n)
# se = 1/sqrt (n), n = 1/(se)^2
# For 10% to be 2.8 std from 0, we need n > (2.8/0.1)^2
# n > 784 
  
# Ex: Epidiemiology.  80% control, 20% treatment, 
# se(diff) = sqrt (0.5*0.5/0.8n + 0.5*0.5/0.2n)
#          = 0.5/sqrt(n)*sqrt(1/0.8 + 1/0.2)
#          = 1.25/sqrt(n) => n > (1.25/se)^2
# sample size, n, to achieve 80% power to detect 10% difference
# n > (2.8*1.25/0.1)^2
# n > 1225

#
# Continuous Outcomes - population std, sigma, has to be specified
#
# theta   - true population mean
# thata_0 - population mean post treatment
# y1, y2, ..., yn - sample
# y_bar - sample mean
# sigma - population standard deviation
# se = standard error in sample

# se = sigma/sqrt(n)
# To achieve a certain se, n > (sigma/se)^2
# If goal is to have 80% power to distinguish between theta and theta_0
# n > (2.8*sigma/(theta-theta_0))^2

# Comparison of means
se_diff_y1_bar_y2_bar <-  sqrt (sigma_1^2/n1 + sigma_2^2/n2)
# if sigma_1=sigma_2=sigma, 
se_diff_y1_bar_y2_bar  <- sigma* sqrt (1/n1 + 1/n2)
# if n1=n2=n, 
se_diff_y1_bar_y2_bar <-  sqrt (sigma_1^2 + sigma_2^2) / sqrt(n)
n <- (sigma_1^2 + sigma_2^2)/se^2
# if goal is to detect a difference of delta, ie se=delta, and 
# want to achieve 80% power:
n <- (sigma_1^2 + sigma_2^2) * (2.8/delta)^2
# if sigma_1=sigma_2=sigma,
n <- (5.6*sigma/delta)^2 

 # 
#Estimating standard deviatons from previous studies
#
# Rosato etal.
#
# delta - Effect size
# sigma - population standard deviation
#
# delta <- avg (wo/zinc) - avg (w/zinc)
# se    <- sqrt (sum_i (se_i^2)/ sum_i(1))

delta    <- (1.1+1.4)/2 - (0.7+0.8)/2                    # 0.5 episodes/year
se_delta <- sqrt ((0.2^2 +0.2^2)/4 + (0.1^2 +0.1^2)/4)   # 0.15

est_treatment_effect  <- delta/(1.1+1.4)/2               # 0.40

# in-group sigmas:  se = sigma/sqrt(n) => sigma = se*sqrt(n)
sigma_placebo <- 0.2 * sqrt (56) # 1.5
sigma_fe      <- 0.2 * sqrt (54) # 1.5
sigma_z       <- 0.1 * sqrt (54) # 0.7
sigma_fe_z    <- 0.1 * sqrt (55) # 0.7

# delta <- sqrt (1.5^2/n/2 + 0.7^2/n/2) = 2.3/sqrt(n)
n <- (2.3/delta)^2
# 80% power
n_80 <- n * 2.8^2 # 166 people in 2 groups one with one without zinc

# Ruel etal
8.1 * 365/100 # placebo: 30 episodes/year
6.3 * 365/100 # treatment: 23 episodes/year
delta <- 8.1-6.3    # 1.8

est_treatment_effect  <- delta/8.1                       # 0.22 - n needs to be higher than
                                                         # Rosato, which has ETE 0.40
se_placebo <- (10.2-5.8)/4  # 1.1
se_treated <- (8.9-4.2)/4   # 1.2

se_delta <- sqrt ((se_placebo^2 + se_treated^2))  # 1.6   ?????

sigma_placebo <- 1.1 * sqrt (44)
sigma_treated <- 1.2 * sqrt (45)

# Lira etal, logarithmic scale (consider 5mg treatment only)
delta <- log (0.68)                        # -0.4

# convert uncertainity [0.49, 0.95] from logarithmic to additive scale
#  => [-0.71, -0.05]
se_delta <- (0.71 - 0.05)/4  # 0.17

# delta is 0.4/0.17= 2.4 se's away from zero, need 2.8 se's for 80% power
n = 71
n_80 <- 71 * (2.8/2.4)^2  # 97




# Sample size
#
n <- 1000
sigma <- 10
y <- rnorm (n, 0, sigma)

x1 <- sample (c(-0.5, 0.5), n, replace = TRUE)
x2 <- sample (c(-0.5, 0.5), n, replace = TRUE)
fake <- data.frame(y, x1, x2)
fit1 <- stan_glm(y ~ x1, data=fake)
print (fit1)
fit2 <- stan_glm(y ~ x1 + x2 + x1:x2, data=fake)
print (fit2) # se(x1)=se(x2)= 0.6 = 2*sigma/sqrt(n) = 2*10/sqrt(1000) 

x1 <- sample (c(0, 1), n, replace = TRUE)
x2 <- sample (c(0, 1), n, replace = TRUE)
fake <- data.frame(y, x1, x2)
fit1 <- stan_glm(y ~ x1, data=fake)
print (fit1)
fit2 <- stan_glm(y ~ x1 + x2 + x1:x2, data=fake)
print (fit2) # se(x1)=se(x1)= 0.9 Edge-condition
#
# 16-6 Design with Fake Data  ****************
#
n <- 100
y_if_control <- rnorm (n, 60, 20)
y_if_treated <- y_if_control + 5  # true treatment effect = 5
#z <- sample (rep (c(0,1), n/2))
z <- sample (c(0,1), n, replace = TRUE)
y <- ifelse (z==1, y_if_treated, y_if_control)

diff <- mean (y[z==1]) - mean (y[z==0])
se_diff <- sqrt ((sd (y[z==0]))^2/sum(z==0) + (sd (y[z==1]))^2/sum(z==1))  # 4.1
# or
fake <- data.frame (y, z)
fit_1a <- stan_glm (y ~ z, data=fake, refresh=0)
print (fit_1a, digits=2) # z: 7.1 (4.0)
# inclease n by 4, n=400, se_diff: 2.0, z:3.6 (2), se decreases by sqrt(4)

#
# add pre-test score, *DONT* link post test control+treated to pre_test
n <- 100
fake$x <- rnorm(n, 0, 12)
fit_1b <- stan_glm (y ~ z + x, data=fake, refresh=0)
print (fit_1b, digits=2) # z: 3.8 (4.0) # s.e doesn't change

#
# add pre-test score, *DO* link post-test control+treated to pre-test
n <- 100
true_ability <- rnorm(n, 50, 16)
x <- true_ability + rnorm(n, 0, 12)
y_if_control <- true_ability + rnorm(n, 0, 12) + 10
y_if_treated <- y_if_control + 5     # true treatment effect = 5
z <- sample(rep(c(0,1), n/2))
#z <- sample(c(0,1), n, replace=TRUE)
y <- ifelse (z==1, y_if_treated, y_if_control)

fake2 <- data.frame (x, y, z)
fit_2a <- stan_glm (y ~ z,     data=fake2, refresh=0)
fit_2b <- stan_glm (y ~ z + x, data=fake2, refresh=0)
print (fit_2a, digits=2) # z: 6.5 (4.0)
print (fit_2b, digits=2) # z: 7.9 (3.2)  *** s.e. decreases
#
# Now assume, treatment is not given out arbitrarily, but to students
# who performed poorly on the pre-test, with probability:
# PR (zi=1) = invlogit (-(xi-50)/20)
# This will lead to biased estimate of the treatment effect, which is 5 ***

true_ability <- rnorm(n, 50, 16)
x <- true_ability + rnorm(n, 0, 12)
z <- rbinom (n, 1, invlogit (-(x-50)/20))
prob <- invlogit (-(x-50)/20)

plot (x, z)
curve (invlogit (-(x-50)/20), from=0, to=100, add=TRUE)

fake <- data.frame(x, prob, z)
ggplot (fake) + geom_point (aes(x=x, y=z)) + 
  geom_line (aes(x=x, y=prob)) +
  xlab ("Pre-test score") + ylab ("Probability of getting tratment") +
  ggtitle ("Treatment status v Pre-test score")

Experiment <- function (n){
  true_ability <- rnorm(n, 50, 16)
  x <- true_ability + rnorm(n, 0, 12)
  y_if_control <- true_ability + rnorm(n, 0, 12) + 10
  y_if_treated <- y_if_control + 5
  z <- rbinom (n, 1, invlogit (-(x-50)/20))
  y <- ifelse (z==1, y_if_treated, y_if_control)
  fake3 <- data.frame (x, y, z)
  fit_3a <- stan_glm (y ~ z,     data=fake3, refresh=0)
  fit_3b <- stan_glm (y ~ z + x, data=fake3, refresh=0)
  inferences <- rbind (c(coef(fit_3a)["z"], se(fit_3a)["z"]),
                       c(coef(fit_3b)["z"], se(fit_3b)["z"]))
  return (inferences)
}

n <- 100
n_loop <- 50
results <- array (NA, c(n_loop, 2, 2), 
                  dimnames = list (1:n_loop, c("Simple", "Adjusted"), c("Estimate", "SE")))
for (loop in 1:n_loop){
  results[loop, , ] <- Experiment(n)
}

results_avg <- apply (results, c(2,3), mean)
results_avg
# Ex 16-1
# a)
# se = sqrt (p*(p-1)/n) => 
p <- 0.5
se <- 0.03
n = round(p*(1-p)/se^2) # 278
# b)
p <- 0.14
se <- 0.03
n = round(p*(1-p)/se^2) # 134
#c)
se <- 0.01
n = round(p*(1-p)/se^2) # 1204

#
# Ex 16-2
#

#se = sqrt (p1*(1-p1)/n1 + p2*(1-p2)/n2)

se_prop_diff <- function (N, target_se){
p1 <- 0.45
p2 <- 0.35
n1 <- p1 * N # 225
n2 <- p2 * N # 175
se = sqrt (p1*(1-p1)/n1 + p2*(1-p2)/n2)
return (se - target_se)
}

root <- uniroot(se_prop_diff, interval=c(0,100000), 
                target_se = 0.05,
                tol= 0.000001)$root 
root # 480

#
# Ex 16-4
#
n <- 500
p_bar <- 0.5

# p = p_bar +/- 1.96 * se                - 5% error
# p = p_bar +/- 1.96 * sqrt (p_bar*(1-p_bar)/n)
# p = p_bar +/- 1.96 * 0.5 / sqrt (n)    - 5% error

# 50% power
p <- p_bar+ 1.96 * 0.5/sqrt(n) 
(p-p_bar) * 100   # 4.4%

# 80% power
p <- p_bar+ 2.8 * 0.5/sqrt(n) 
(p-p_bar) * 100   # 6.2%
###

n<- 96
# power=0.5
p <- p_bar+ 1.96 * 0.5/sqrt(n) 
(p-p_bar) * 100   # 10%
# power=0.8
p <- p_bar+ 2.8 * 0.5/sqrt(n) 
(p-p_bar) * 100   # 14.3%

#
# Ex 16-6
# a)
# p = p_bar +/- 2.8 * se   - 5% error, 80% power
# p-pbar <- 2.8 * sqrt (pbar*(1-pbar)/n)
p<- 0

pbar<- 0.0001
n <- 2.8*sqrt(pbar*(1-pbar))/(p-pbar)^2  # 2.8m
pbar<- .01
n <- 2.8*sqrt(pbar*(1-pbar))/(p-pbar)^2  # 2.8k
pbar<- 1
n <- 2.8*sqrt(pbar*(1-pbar))/(p-pbar)^2  # 0

# b)
# log(p) = log(p_bar) +/- 2.8 * se   - 5% error, 80% power
# p-pbar <- 2.8 * sqrt (pbar*(1-pbar)/n)
p<- 0

pbar<- log(100) * 0.0001
n <- 2.8*sqrt(pbar*(1-pbar))/(p-pbar)^2  # 283k
pbar<- log(10000) * 0.0001
n <- 2.8*sqrt(pbar*(1-pbar))/(p-pbar)^2  # 100k

#
# Ex 16-7
#
# a)
p <- 0.5
n <- 1000

n_i <- n/J

se_i <- sqrt (p*(1-p)/n_i)
se_i <- 0.5/sqrt(n_i)     

# b)
delta <- 0.02
n_delta_2p <- (0.5/delta)^2 # 625

cost <- n_delta_2p*50 + J*500

# 16-8
electric<- read.table("electric", header=TRUE)
a <- subset (electric, electric$Grade==2)
dim (a) # 34 classes

sigma_treated <- sd ((a$treated.Pretest+a$treated.Posttest)/2)/100 # 0.10
sigma_control <- sd ((a$control.Pretest+a$control.Posttest)/2)/100 # 0.12

avg_improvement <- mean ((a$treated.Posttest+a$treated.Pretest)/2) - 
                   mean ((a$control.Posttest+a$control.Pretest)/2)  # 6.77

# 0.05 <- sqrt (sigma_treated^2/n + sigma_control^2/n)
n <- (sigma_treated^2 + sigma_control^2)*(2.8/0.05)^2    #  80 classes 
# Have data on 34 classes
n_34 <- 34
se_34 = 0.077
# 34 <- (sigma_treated^2 + sigma_control^2)*(2.8/se_34)^2   #  34 classes 
#34 classes can provid 80% power on se of 7.7

sigma_treated <- sd ((a$treated.Posttest-a$treated.Pretest)/2)/100 # 0.0302
sigma_control <- sd ((a$control.Posttest-a$control.Pretest)/2)/100 # 0.0300
# 0.05 <- sqrt (sigma_treated^2/n + sigma_control^2/n)
n <- (sigma_treated^2 + sigma_control^2)*(2.8/0.05)^2     # 6 classes


# use regression
x_all <- c(a$treated.Pretest,  a$control.Pretest)
y_all <- c(a$treated.Posttest, a$control.Posttest)
t_all <- rep (c(1,0), rep (nrow (a),each=2))

fake_all <- data.frame(x_all, y_all, t_all)
fit_all <- stan_glm(y_all ~ x_all + t_all, data = fake_all)
print (fit_all, digits=2)
# t_all: 4.3 (1.35): Avg improvement of treatment is 4.3(1.35) points 
se <- 0.0135


# Ex 16-9
mean_control <- (1.1+1.4)/2  # 1.25
mean_treated <- (0.7+0.8)/2  # 0.75

sd_placebo  <- 0.2*sqrt(56)  # 1.5
sd_fe       <- 0.2*sqrt(54)  # 1.5
# (110 people in control)
sd_fe_z     <- 0.1*sqrt(54)  # 0.7
sd_z        <- 0.1*sqrt(55)  # 0.7
# (109 people in treated)
sd_control  <- 1.5
sd_treated  <- 0.7

diff    <- mean_control - mean_treated  # 0.5 
se_diff <- sqrt ((0.2^2 + 0.2^2)/4 + (0.1^2 + 0.1^2)/4)  # 0.16
se_diff <- sqrt (sd_control^2 + sd_treated^2)/sqrt(110)  # 0.16

n <- 100
control <- rnorm (n, mean_control, sd_control) 
zinc    <- rnorm (n, mean_treated, sd_treated) 

# let n_control=n_treated=n/2
se <- sqrt (sd_control^2 + sd_treated^2)/sqrt(110)  # 1.6

se_5 <- 0.05
n_80 = (sd_control^2 + sd_treated^2)/(se_5/2.8)^2  # 17186

n_control <- 8592
n_treated <- 8592

se = 2.8 * sqrt (sd_control^2/n_control + sd_treated^2/n_treated)

cost_exp <- function (n_control){
  n_treated <- sd_treated^2/((0.05/2.8)^2 - sd_control^2/n_control)
  cost = 100 * n_control + 150 * n_treated
  return(c(n_treated, cost))
}

n_control <- rep (NA, 200)
n_treated <- rep (NA, 200)
n_control <- seq (8000, 8000+10000, 50)
n_treated <- sd_treated^2/((0.05/2.8)^2 - sd_control^2/n_control)

n_total  <- n_control + n_treated
cost = 100 * n_control + 150 * n_treated

exp_param <- data.frame (n_control = n_control/1000, 
                         n_treated = n_treated/1000, 
                         n_total   = n_total/1000, 
                         cost      = cost/1000000)
plot (n_control, cost)
plot (n_control, n_treated)
plot (n_control, n_total)

x_min_total <- exp_param$n_control[which (n_total == min (n_total))]
y_min_total <- exp_param$n_total[which (n_total == min (n_total))]
x_min_cost  <- exp_param$n_control[which (cost == min (cost))]
y_min_cost  <- exp_param$cost[which (cost == min (cost))]

p1 <- ggplot (exp_param) + 
  geom_line(aes(x=n_control, y=n_treated)) + 
  xlab ("n control(k)") + ylab ("n treated(k)") +
  ggtitle ("n control v. n treated") 
p2 <- ggplot (exp_param) + 
  geom_line(aes(x=n_control, y=n_total)) + 
  annotate("text", x=14, y=15.18, label="min:(10.4k,15.2k)") +
  geom_segment(x=10.15, xend=10.55, y=15.18, yend=15.18, color='red') +
  geom_segment(x=10.35, xend=10.35, y=15.08, yend=15.28, color='red') +
  xlab ("n control(k)") + ylab ("n total(k)") +
  ggtitle ("n total v. n control") 
p3 <- ggplot (exp_param) + 
  geom_line(aes(x=n_control, y=cost)) + 
  annotate("text", x = 15, y = 1.74,   label="min:(11.1k,1.74m)") +
  geom_segment(x=10.9, xend=11.3, y=1.74, yend=1.74, color='red') +
  geom_segment(x=11.1, xend=11.1, y=1.72, yend=1.76, color='red') +
  xlab ("n control(k))") + ylab ("Cost($m)") +
  ggtitle ("cost($m) v. n control") 

grid.arrange (p1, p2, p3, ncol=3)  

indx <- which (exp_param$cost==min (exp_param$cost))
cat ("Minimum Cost=", exp_param$cost[indx]/1000000, "m, ", 
     "Control size=", round(exp_param$n_control[indx]/1000,2), "k, ", 
     "Treated size=", round(exp_param$n_treated[indx]/1000,1), "k")
#****************************

#
# Ch 17:  Poststratification and Missing-Data Imputation
#

#  **Create some fake data**
n_pid <- c(254, 282, 242)
n <- sum(n_pid)
pid_names <- c("Republican", "Democrat", "Independent")
pid <- rep(pid_names, n_pid)
n_vote <- as.list(rep(NA, 3))
n_vote[[1]] <- round(c(0.77, 0.08)*n_pid[1])
n_vote[[2]] <- round(c(0.05, 0.89)*n_pid[2])
n_vote[[3]] <- round(c(0.36, 0.38)*n_pid[3])
vote <- NULL
y_bar_cells <- rep(NA, 3)
for (j in 1:3){
  n_vote[[j]]<- c(n_vote[[j]], n_pid[j] - sum(n_vote[[j]]))
  vote <- c(vote, rep(c(1, 0, NA), n_vote[[j]]))
  y_bar_cells[j] <- mean(vote[pid==pid_names[j]], na.rm=TRUE)
  round(y_bar_cells[j], 3)
}
poll <- data.frame(vote, pid)

fit1 <- stan_glm (vote ~ pid, data = poll)
print (fit1, digits=2)

poststrat_data <- data.frame(pid = c("Republican", "Democrat", "Independent"), 
                             N = c(0.33, 0.36,0.31))

# posterior uncertainty  in the support for Trump in each of the 3 categories:
predict1 <- posterior_linpred(fit1, newdata = poststrat_data) 

# uncertainty in population average
poststrat_pred1 <- predict1 %*% poststrat_data$N/sum (poststrat_data$N)
print (c(mean(poststrat_pred1), mad(poststrat_pred1)), digits=2)   # 0.47, 0.013

# If we assume there is 2% uncertainty in polls
n_sims <- nrow (predict1)

poststrat_pred2 <- poststrat_pred1 + rnorm (n_sims, 0, 0.02)
print (c(mean(poststrat_pred2), mad(poststrat_pred2)), digits=2)   # 0.47, 0.024

#
# 17.2 Fake data simulation for regression and post-stratification
#
# y - Yes/No
# x - sex, age (18–29, 30–44, 45–64, 65+),  
#     ethnicity (non-hispanic white, black, hispanic, other)
#     2 X 4 X 4  categories

J <- c(2, 4, 4)
poststrat <- as.data.frame(array(NA, c(prod(J), length(J)+1)))  # 32 by (3+1)
colnames (poststrat) <- c("sex", "age", "eth", "N")
count <- 0
for (i1 in 1:J[1]){
  for (i2 in 1:J[2]){
    for (i3 in 1:J[3]){
     count <- count +1
     poststrat[count, 1:3] <- c(i1, i2, i3)
    }
  }
}
# Sub-populations, N, for a 250mm country
p_sex <- c(0.52, 0.48)
p_age <- c(0.20, 0.25, 0.30, 0.25)
p_eth <- c(0.7,  0.1,  0.1,  0.1 )

for (j in 1:prod(J)){
  poststrat$N[j] <- 250e6 * p_sex[poststrat[j,1]] * p_age[poststrat[j,2]] * p_eth[poststrat[j,3]] 
}

# Non-response probabilitiies
p_response_baseline <- 0.1
p_response_sex <- c(1, 0.8)
p_response_age <- c(1, 1.2, 1.6, 2.5)
p_response_eth <- c(1, 0.8, 0.7, 0.6)

p_response <- rep(NA, prod(J))
for (j in 1:prod(J)){
  p_response[j] <- p_response_baseline * p_response_sex[poststrat[j,1]] *
    p_response_age[poststrat[j,2]] * p_response_eth[poststrat[j,3]]
}

# Now sample
#
n <- 1000
# people[i] is the i'th respondents poststrat cell (from 1 to 32)
people <- sample (prod(J), n, replace = TRUE, prob = poststrat$N*p_response)

n_cell <- rep (NA, prod(J))
for (j in 1:prod(J)){
  n_cell[j] <- sum (people==j)
}

print (cbind (poststrat, n_cell/n, poststrat$N/sum(poststrat$N)))

# ASSUME survey responses, fitted with logistic regression
coef_intercept <- 0.6
coef_sex <- c(0, -0.2)
coef_age <- c(0, -0.2, -0.3, -0.4)
coef_eth <- c(0,  0.6,  0.3,  0.3)

# YES probabilities
prob_yes <- rep (NA, prod(J))
for(j in 1:prod(J)) {
  prob_yes[j] <- invlogit(coef_intercept + coef_sex[poststrat[j,1]] +
                            coef_age[poststrat[j,2]] + coef_eth[poststrat[j,3]])
}

y <- rbinom (n, 1, prob_yes[people])

# Perform regression and poststratification
sex <- poststrat[people, 1]
age <- poststrat[people, 2]
eth <- poststrat[people, 3]

fake <- data.frame(sex, age, eth)

fit <-stan_glm (y ~ factor (sex) + factor (age) + factor(eth), 
                family=binomial(link="logit"), data = fake)

print (fit, digits=2)

# Estimate proportion of Yes responses
pred_sim <- posterior_linpred(fit, transform=TRUE,newdata=as.data.frame(poststrat))
pred_est <- colMeans(pred_sim)

print(cbind(poststrat, prob_yes, pred_est))

poststrat_sim <- pred_sim %*% poststrat$N / sum(poststrat$N)
round(c(mean(poststrat_sim), sd(poststrat_sim)), 3)

# 
# Missingness
#


# Exercise 17.1
#
vitals_full <- read.table ("vitals.dat", header = T, sep =',')
# Clean up data
minidata <- vitals_full[, c("weight","height","female","ethnicity","exercise","smokenow")]
minidata$weight[minidata$weight>990]=NA;
ok <- apply(is.na(minidata), 1, sum) == 0 # only rows that do not have NAs

vitals <- minidata[ok,]

vitals$c_height <- vitals$height - 66

ggplot (vitals, aes(x=c_height, y=weight, color=factor(female))) + geom_point()

# predictor: c_height
fit_1 <- stan_glm(weight ~ c_height, data = vitals, refresh = 0)  # coef (height) = 4.9
print (fit_1)

sum (vitals$female)/nrow (vitals)  # 62%
# Female ratio in population = 52%, in sample, 62%:  Need to poststratify
sum(vitals$height)/nrow(vitals) # avg height in sample 66.6 in 

sum(vitals$height*(vitals$female==1))/sum(vitals$female==1)   # avg woman's height = 64.5 in 
sum(vitals$height*(vitals$female==0))/sum(vitals$female==0)   # avg man's height 70 in

avg_height_in_sample <- sum(vitals$height)/nrow(vitals) # 66.6 in 
avg_height_in_pop    <- 0.48 * 70 + 0.52 * 64.5         # 67.1 in
height_adj <- avg_height_in_pop - avg_height_in_sample  # 0.55 in

# Fit using female as additional predictor
fit_2 <- stan_glm(weight ~ c_height + female, data = vitals, refresh = 0)
print (fit_2, digits = 2) # intercept=161.6, coef(female)=-12, coef(height)=3.8, sigma = 28

# POST_STRATIFY: generate population where 52% is women and obey rules above
#
post_stratting <- function (fit_in, height, female){
  fake <- data.frame (c_height=round(height), female=female)
  qq <- posterior_linpred (fit_in, newdata = fake)
  weight <- apply (qq, 2, median)

  sigma <- median(as.matrix(fit_in, pars = "sigma"))
  
  weight <- weight + rnorm (n, 0, sigma)
  data <- data.frame (weight = weight, c_height = height, female = female)
  
  return (data)
}

n       <- nrow (vitals)   # 1983
n_men   <- round (n*0.48)  #  952
n_women <- round (n*0.52)  # 1031
min_c_height <- min(vitals$c_height)  # -9
max_c_height <- max(vitals$c_height)  #  16
post_strat_height <- sort(sample (min_c_height:max_c_height, n, replace=TRUE))

# 1) Females at arbitrary height
post_strat_height <- sample (c(0,1),n, replace=TRUE, prob=c(0.48, 0.52))

post_strat_data <- post_stratting (fit_in=fit_2, height=post_strat_height, female=post_strat_female )

ggplot (post_strat_data, aes(x=c_height, y=weight, color=factor(female))) + geom_point()

fit_3 <- stan_glm(weight ~ c_height + female, data = post_strat_data)
print (fit_3, digits = 2) # intercept=160.9, coef(female)=-12, coef(height)=3.8, sigma = 29
#
# 2) Females at lower height
post_strat_female <- c(rep(1, n_women),rep(0, n_men))

post_strat_data <- post_stratting (fit_in=fit_2, height=post_strat_height, female=post_strat_female )

ggplot (post_strat_data, aes(x=c_height, y=weight, color=factor(female))) + geom_point()

fit_4 <- stan_glm(weight ~ c_height + female, data = post_strat_data)
print (fit_4, digits = 2) # intercept=160.5, coef(female)=-12, coef(height)=3.7, sigma = 27.5
#
# 3) Females women at higher height
post_strat_female <- c(rep(0, n_men), rep(1, n_women))

post_strat_data <- post_stratting (fit_in=fit_2, height=post_strat_height, female=post_strat_female )

ggplot (post_strat_data, aes(x=c_height, y=weight, color=factor(female))) + geom_point()

fit_5 <- stan_glm(weight ~ c_height + female, data = post_strat_data)
print (fit_5, digits = 2)  # intercept=162.6, coef(female)=-15, coef(height)=4.0. sigma=27.5

#
# 7.1-b Intercation between female:height
#
fit_6 <- stan_glm(weight ~ c_height + female + female:c_height, data = poststrat_data)
print (fit_6, digits = 2) # intercept=162.0, coef(female)=-13, coef(height)=3.7. sigma=28
# c_height:female median:0.32   sd:0.36 - not a good predictor 
# CI [-0.40, 1.04] contains zero

# 1) Females at arbitrary height
post_strat_female <- sample (c(0,1), n, replace=TRUE, prob=c(0.48, 0.52))
post_strat_data <- post_stratting (fit_in=fit_6, height=post_strat_height, female=post_strat_female )
ggplot (post_strat_data, aes(x=c_height, y=weight, color=factor(female))) + geom_point()
fit_7 <- stan_glm(weight ~ c_height + female + female:c_height, data = post_strat_data)
print (fit_7, digits = 2) # intercept=162.7, coef(female)=-13, coef(height)=3.6. sigma=28
# c_height:female median:0.12   sd:0.16 - not a good predictor 
# CI [-0.04, 0.28] contains zero

# 2) Females at lower height
post_strat_female <- c(rep(1, n_women),rep(0, n_men))
post_strat_data <- post_stratting (fit_in=fit_6, height=post_strat_height, female=post_strat_female )
ggplot (post_strat_data, aes(x=c_height, y=weight, color=factor(female))) + geom_point()

fit_8 <- stan_glm(weight ~ c_height + female + female:c_height, data = post_strat_data)
print (fit_8, digits = 2)  # intercept=160, coef(female)=-10, coef(height)=3.8. sigma=28
# c_height:female median:0.22   sd:0.34 - not a good predictor 
# CI [-0.12, 0.56] contains zero

#
# CH 19 Causal Inference
#
a <- read.table("electric", header=TRUE)

par (mfrow=c(4, 2))
for (i in 1:4){
  hist(a$control.Posttest[a$Grade==i], breaks = seq(0, 125, 5))
  abline (v=mean(a$control.Posttest[a$Grade==i]), lwd = 3)
  hist(a$treated.Posttest[a$Grade==i], breaks = seq(0, 125, 5))
  abline (v=mean(a$treated.Posttest[a$Grade==i]), lwd = 3)
}

par (mfrow=c(2, 4))
for (i in 1:4){
  hist(a$control.Posttest[a$Grade==i], breaks = seq(40, 125, 5))
  abline (v=mean(a$control.Posttest[a$Grade==i]), lwd = 3)
}
for (i in 1:4){
  hist(a$treated.Posttest[a$Grade==i], breaks = seq(40, 125, 5))
  abline (v=mean(a$treated.Posttest[a$Grade==i]), lwd = 3)
}

#
# 
x_all <- c(a$treated.Pretest,  a$control.Pretest)
y_all <- c(a$treated.Posttest, a$control.Posttest)
t_all <- rep (c(1,0), rep (nrow (a),each=2))
grade_all <- c(a$Grade, a$Grade)

fake_all <- data.frame(x_all, y_all, t_all, grade_all)
fit_all <- stan_glm(y_all ~ t_all, data = fake_all)
print (fit_all, digits=2)

ggplot (fake_all, aes(x=t_all, y=y_all, color=factor(t_all)))+ geom_point() +
  ggtitle("Posttest results as a function of Treatment") +
  geom_abline(intercept=mean (y_all[t_all==0]), slope=0, color="red", size=1) +
  geom_abline(intercept=mean (y_all[t_all==1]), slope=0, color="green", size=1) +
  annotate("text", x = 0.10, y = 90,   label = "treat=no\navg=94.3") +
  annotate("text", x = 0.90, y = 105,  label = "treat=yes\navg=100") +
  xlab ("Treatment Yes/No") + ylab ("Posttest")

mean (y_all[t_all==0])  # 94.3
mean (y_all[t_all==1])  #100.0

t_coef <- rep (NA, 4)
t_se <- rep (NA, 4)
grade <- rep (NA, 4)
fit_1 <- as.list (rep(NA, 4))
for (i in 1:4){
  fit_1[[i]] <- stan_glm(y_all ~ t_all, data=fake_all, subset=(grade_all==i))
  sims <- as.matrix(fit_1[[i]])
  t_coef[i] <- median (sims[,2])
  t_se[i]   <- sd (sims[,2])
  grade[i] <- i
}  


for (i in 1:4){
  cat ("Grade", grade[i], "Treatment=", round(t_coef[i], 2),  
       " - CI=[", round(t_coef[i]-2*t_se[i], 2),",", 
       round(t_coef[i]+2*t_se[i], 2), "]\n") 
  
}

treat <- data.frame(grade, t_coef, t_se)  
ggplot (treat, aes(x=t_coef, y=grade)) + #
  geom_segment (x=t_coef-2*t_se, y=grade, xend=t_coef+2*t_se, yend=grade, size=1, color="gray") +
  geom_segment (x=t_coef-t_se, y=grade, xend=t_coef+t_se, yend=grade, size=1) +
  geom_point() + xlim(-5, 20) + geom_vline (xintercept = 0, linetype="dotted") +
  xlab ("Treatment") + ylab ("Grade")

#
# Add pre_test to analysys
# y = int + coef_t * t + coef_x * x 
#
t_coef <- rep (NA, 4)
t_se   <- rep (NA, 4)
grade  <- rep (NA, 4)

fit_1 <- as.list (rep(NA, 4))
for (i in 1:4){
  fit_1[[i]] <- stan_glm(y_all ~ t_all + x_all, data=fake_all, subset=(grade_all==i))
  sims <- as.matrix(fit_1[[i]])
  t_coef[i] <- median (sims[,2])
  t_se[i]   <- sd (sims[,2])
  grade[i] <- i
}   

plot_list = list()
for (k in 1:4){
  data_k <- subset(fake_all, grade_all==k)
  fit <- stan_glm(y_all ~ t_all + x_all, data=data_k)
  print (fit, digits=2)
  sims <- as.matrix(fit)
  int     <- median (sims[,1])
  coef_t  <- median (sims[,2])
  coef_x  <- median (sims[,3])
  
  title <- paste ("Grade", k, "- Model: (y=a + b*x + theta*z)")
  plot_list[[k]] <- ggplot (data_k, aes(x=x_all, y=y_all, color=factor(t_all))) + geom_point() +
    geom_abline(intercept = int, slope = coef_x, color = "red")+
    geom_abline(intercept = int + coef_t, slope = coef_x, color = "green") +
    ggtitle (title) + 
    xlab("Pre test, xi") + ylab("Post test, yi") +
    xlim (0, 125) + ylim (0,125)
}


for (i in 1:4){
  cat ("Grade", grade[i], "Treatment=", round(t_coef[i], 2),  
       " - CI=[", round(t_coef[i]-2*t_se[i], 2),",", 
       round(t_coef[i]+2*t_se[i], 2), "]\n") 
}

treat <- data.frame(grade, t_coef, t_se)  
ggplot (treat, aes(x=t_coef, y=grade)) + #
  geom_segment (x=t_coef-2*t_se, y=grade, xend=t_coef+2*t_se, yend=grade, size=1, color="gray") +
  geom_segment (x=t_coef-t_se, y=grade, xend=t_coef+t_se, yend=grade, size=1) +
  geom_point() + xlim(0, 15) + geom_vline (xintercept = 0, linetype="dotted") +
  xlab ("Treatment") + ylab ("Grade")
 
# 
# interction between pretest and treatment
# y = int + coef_t * t + coef_x * x + coef_tx * x * t
#
plot_list = list()
for (k in 1:4){
  data_k <- subset(fake_all, grade_all==k)
  fit <- stan_glm(y_all ~ t_all + x_all + t_all:x_all, data=data_k)
  print (fit, digits=2)
  sims <- as.matrix(fit)
  int     <- median (sims[,1])
  coef_t  <- median (sims[,2])
  coef_x  <- median (sims[,3])
  coef_tx <- median (sims[,4])

  data_k$t_all <- factor(data_k$t_all)
  title <- paste ("Grade", k)#, "- Model: (y=a + b*x + theta*z + c*z*x)")
  plot_list[[k]] <- ggplot (data_k, aes(x=x_all, y=y_all, color=t_all)) + geom_point() +
    geom_abline(intercept = int, slope = coef_x, color = "red")+
    geom_abline(intercept = int+coef_t, slope = coef_x + coef_tx, color = "green") +
    ggtitle (title) + 
    xlab("Pre test, xi") + ylab("Post test, yi") + 
    xlim (0, 125) + ylim (0,125) +
    theme(legend.position = c(0.65, 0.25)) 
}

grid.arrange (plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], ncol = 4)

# 
# 19-4 Varying treatment effects, interactions, poststratification
#

# predictor: treatment only, grade=4
fit1 <- stan_glm(y_all ~ t_all, data=fake_all, subset=(grade_all==4))
print(fit1) # coeff_t= 3.6 +/- 1.8, sigma=6.0

# predictor: treatment + pretest, grade=4
fit2 <- stan_glm(y_all ~ t_all + x_all, data=fake_all, subset=(grade_all==4))
print(fit2) # coeff_t= 1.7 +/- 0.7, sigma=2.2

# predictor: treatment + pretest + treatment:pretest, grade=4
fit3 <- stan_glm(y_all ~ t_all + x_all + t_all:x_all, data=fake_all, subset=(grade_all==4))
print(fit3, digits=2) # coeff_t= 12 +/- 8, coeff_tx = 0.1 +/- 0.07, sigma=2.2

# Treatment effect = 11.95 + 0.1* x_all
# for grade4  78 < x_all < 120  => Treatment effect ~ [4.15, -0.05]

# or center x_all for grade 4
fake_all$x_all_c4 <- x_all - mean ((subset(fake_all, grade_all==4)$x_all))

fit4 <- stan_glm(y_all ~ t_all + x_all_c4 + t_all:x_all_c4, data=fake_all, subset=(grade_all==4))
print(fit4, digits=2) # coeff_t= 1.8 +/- 0.7, coeff_tx = -0.15 +/- 0.09, sigma=2.2

data_4 <- subset(fake_all, grade_all==4)
sims <- as.matrix (fit3)
pick_sims <- sample (4000, 20)

ggplot (data_4, aes(x=x_all, y=y_all))  +
  geom_abline(intercept = median(sims[,2]),  slope = median(sims[,4])) +
  geom_abline(intercept = sims[pick_sims,2], slope = sims[pick_sims,4], size=0.1) +
  xlim (80, 120) + ylim (-5, 10) +
  xlab ("Pre-test") + ylab ("Treatment Effect") +
  ggtitle("20 Draws exhibiting uncertainty in treatment effect vs pre-test")

n_sims <- nrow(sims)
n_4    <- nrow (data_4) 
effect <- array(NA, c(n_sims, n_4))
for (i in 1:n_4){
  for (s in 1:n_sims){
    effect[s,i] <- sims[s,2] + sims[s,4]*data_4$x_all[i]
  }
}
avg_effect <- rowMeans(effect)
median (avg_effect) # 1.73
sd (avg_effect)     # 0.70

#
# Ex:19-4
#

# post_test <- a + b * pre_test + theta * z + error   
n <- 100
b <- 0.7
theta <- 10
pre_test  <- rnorm (n, 40, 15) 

# Control, z=0, 
# E(post_test_control) = 50 
# E(post_test_control) <- a + b * E(pre_test) + E(error)
# 50 = a + 0.7 * 40
a = 50 - 0.7 * 40  # 22

#sd (post_test_control) = sqrt (b*sd(pre_test)^2 + sd(error)^2)
# let error ~ N(0, err_sd)
b <- 0.7
err_sd <- 10
n_control <- 50
pre_test <- rnorm (n_control, 40, 15) 
error    <- rnorm (n_control,  0, err_sd)

post_test_control <- a + b * pre_test + error

fake <- data.frame(pre_test, post_test_control)
ggplot (fake, aes(x=pre_test, y=post_test_control)) + geom_point() 

est_sd_post_test <- sqrt ((b*sd(pre_test))^2 + sd(error)^2)

cat ("mean post_test control=", round(mean (post_test_control), 0), 
     "ESTIMATED: mean post_test control=", 50, "\n")
cat ("sd post_test control=",   round(sd (post_test_control), 0),
     "ESTIMATED: sd post_test control=", round (est_sd_post_test, 0))

#
# Regression

n <- 100
post_test <- rep (NA, n)
b <- 0.7
theta <- 10
pre_test  <- rnorm (n, 40, 15) 
z <-  c(rep (1,50), rep (0,50))


err_sd <- 9.5 
error    <- rnorm (n,  0, err_sd)
post_test_control <- a + b * pre_test + error

post_test <- a + b * pre_test + error
post_test[z==1] <- post_test[z==1] + theta

cat ("mean post-test score:  No Treatment=", round (mean(post_test[z==0])), 
     "   With Treatment=", round (mean(post_test[z==1])))

fake <- data.frame(pre_test, post_test)
fit <- stan_glm(post_test ~ pre_test + z, data = fake)
print (fit) # Residual sd=10  when error ~ N(0, 9.5)

#
# Ex 19-6 Cows
#
cows <- read.table ("cows", header = T)

ggplot (cows, aes(x=level, y=fat)) + geom_point () 
ggplot (cows, aes(x=initial.weight, y=final.weight)) + geom_point () 

# final.weight vs initial.weight + treatment as predictor
fit_w <- stan_glm (final.weight ~ factor(level) + initial.weight, data = cows)
print (fit_w, digits=2) # level=2.00 (0.55), sigma=0.44

# fat vs level
fit_1 <- stan_glm (fat ~ level, data = cows)
print (fit_1, digits=2) # level=2.00 (0.55), sigma=0.44
# fat vs level + lactation, age, initial.weight
fit_2 <- stan_glm (fat ~ level + lactation + age + initial.weight, data = cows)
print (fit_2, digits=2) # level=2.87 (0.53), sigma=0.42
# fat vs level + lactation
fit_3 <- stan_glm (fat ~ level + lactation, data = cows)
print (fit_3, digits=2) # level=2.05 (0.54), sigma=0.42
# fat v level + lactation, initial.weight, lactation*initial.weight
fit_4 <- stan_glm (fat ~ level + lactation + initial.weight + lactation:initial.weight, data = cows)
print (fit_4, digits=2) # level=2.04 (0.55), sigma=0.43

# fat v factor(level)
fit_1l <- stan_glm (fat ~ factor(level), data = cows)
print (fit_1l, digits=2) # sigma=0.44
# fat vs factor(level) + lactation, age, initial.weight
fit_2l <- stan_glm (fat ~ factor(level) + lactation + age + initial.weight, data = cows)
print (fit_2l, digits=2) # sigma=0.43

plot (fit_1l)
plot_model (fit_1l)

sims  <- as.matrix(fit_1l)
level <- c(0.1, 0.2, 0.3)
coef  <- c(0.09, 0.25, 0.60)
se    <- c(0.18, 0.18, 0.18)

for (i in 1:3){
  cat ("Level", level[i],":", coef[i], "(", se[i], ")\n")
}

treat <- data.frame(level, coef, se)  
ggplot (treat, aes(x=coef, y=level)) + #
  geom_segment (x=coef-2*se, y=level, xend=coef+2*se, yend=level, size=1, color="gray") +
  geom_segment (x=coef-se,   y=level, xend=coef+se,   yend=level, size=1) +
  geom_point() + xlim(-0.5, 1) + geom_vline (xintercept = 0, linetype="dotted") +
  xlab ("Treatment Level") + ylab ("Fat")

#
# Ex 19-9
c_86 <- read.table("congress_86", header=F)
c_86$percent_dem <- c_86$V4/(c_86$V4+c_86$V5)

bad_86 <- c_86$V4==-9 | c_86$V5==-9
c_86$percent_dem[bad_86] <- NA

c_88 <- read.table("congress_88", header=F)
c_88$percent_dem <- c_88$V4/(c_88$V4+c_88$V5)

bad_88 <- c_88$V4==-9 | c_88$V5==-9
c_88$percent_dem[bad_88] <- NA
inc_88 <- ifelse (c_88$V3 == -9, 0, c_88$V3)

data_88 <- data.frame(inc_party = ifelse (c_86$percent_dem > 0.5, 1, -1),
                      past_vote = c_86$percent_dem,
                      vote      = c_88$percent_dem,
                      win_88    = ifelse (c_88$percent_dem > .5, 1, -1),
                      inc_88,
                      inc         = factor(ifelse (inc_88 !=0, 1, 0)),
                      past_vote_margin = ifelse (c_86$percent_dem > 0.5, 
                                                 c_86$percent_dem-0.5,
                                                 0.5-c_86$percent_dem),
                      vote_margin      = ifelse (c_88$percent_dem > 0.5, 
                                                 c_88$percent_dem-0.5,
                                                 0.5-c_88$percent_dem)
                      )

data_88_contested <- subset (data_88, past_vote > .25 & past_vote < .75 & vote > .25 & vote < .75)

cat (length (data_88_contested$inc), "contested elections,", 
     sum (data_88_contested$inc), "incumbents,", 
     round ((sum (data_88_contested$inc)/length(data_88_contested$inc))*100, 0), 
     "% incumbents")

ggplot (data_88_contested, aes(x=past_vote, y=vote, color=factor(inc_88))) + geom_point() +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("1986 Dem Vote") + ylab("1988 Dem Vote")

ggplot (data_88_contested, aes(x=past_vote, y=win_88, color=factor(inc_88))) + geom_point() +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("1986 Dem Vote") + ylab("1988 Dem Win")

ggplot (data_88_contested, aes(x=inc_88, y=vote, color=factor(inc_88))) + 
  geom_jitter (width = 0.05, height = 0.0) +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("Incumbency") + ylab("1988 Dem Vote")

ggplot (data_88_contested, aes(x=inc_88, y=win_88, color=factor(inc_88))) + 
  geom_jitter(width = 0.2, height = 0.2) +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("Incumbency") + ylab("1988 Dem Win")

ggplot (data_88_contested, aes(x=inc, y=vote_margin, color=factor(inc_88))) + 
  geom_jitter(width = 0.05, height = 0.0) +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("Incumbency") + ylab("1988 Dem Win")

ggplot (data_88_contested, aes(x=past_vote_margin, y=vote_margin, color=factor(inc))) + 
  geom_point() +
  ggtitle("1988 Congressional Democratic Win") + 
  xlab("Incumbency") + ylab("1988 Dem Win")


fit1 <- stan_glm (vote ~ past_vote, data = data_88_contested)
print (fit1, digits=2) # past_vote=91(3)
fit2 <- stan_glm (vote ~ inc_88, data = data_88_contested)
print (fit2, digits=2) # inc_88=15(1)
fit3 <- stan_glm (win_88 ~ inc_88, data = data_88_contested)
print (fit3, digits=2) # inc_88=94(3)
fit4 <- stan_glm (vote ~ past_vote + inc_88, data = data_88_contested)
print (fit4, digits=2) # past_vote=50(5), inc_88=8(1)
fit5 <- stan_glm (vote ~ past_vote + inc_88 + past_vote:inc_88, data = data_88_contested)
print (fit5, digits=2) # past_vote=49(5), inc_88=4(3), past_vote:inc=7(5)
fit6 <- stan_glm (vote ~ past_vote + inc_88 + past_vote:inc_88 + inc_party, data = data_88_contested)
print (fit6, digits=2) # past_vote=45(6), inc_88=4(3), inc_party=2(2), past_vote:inc=7(5)
fit7 <- stan_glm (vote_margin ~ inc, data = data_88_contested)
print (fit7, digits=2) # inc=7(2)
fit7 <- stan_glm (vote_margin ~ past_vote_margin, data = data_88_contested)
print (fit7, digits=2) # mast_vote_margin=39(5)
fit8 <- stan_glm (vote_margin ~ past_vote_margin + inc, data = data_88_contested)
print (fit8, digits=2) # mast_vote_margin=42(5), inc=9(1)
fit9 <- stan_glm (vote_margin ~ past_vote_margin + inc + inc_88, data = data_88_contested)
print (fit9, digits=2) # mast_vote_margin=43(5), inc=9(1), inc_88=-1(0)

coef (fit8)
beta0 <- coef (fit8)[1]
beta1 <- coef (fit8)[2]
beta2 <- coef (fit8)[3]

ggplot (data_88_contested, aes(x=past_vote_margin, y=vote_margin, color=factor(inc_88))) + #****
  geom_point() +
  geom_abline(intercept = beta0, slope = beta1, color="green") +
  geom_abline(intercept = beta0+beta2, slope = beta1, color="brown") +
  ggtitle("1988 Congressional Winning Margin For Competitive Districts") + 
  xlab("86 Winning Margin") + ylab("88 Winning Margin") +
  scale_color_manual(name = "Incumbent", labels = c("Dem", "None", "GOP"),
                     values = c("blue", "green", "red"))  

#
# Poststratification Example
#
fill_strat_treat_cells <- function (fit, cells){
  pointpred_treated <- predict(fit, newdata=data.frame(cells, z=1))
  pointpred_control <- predict(fit, newdata=data.frame(cells, z=0))
  
  epred_treated <- posterior_linpred(fit, newdata=data.frame(cells, z=1)) 
  epred_control <- posterior_linpred(fit, newdata=data.frame(cells, z=0)) 
  
  pt_treat_eff_in_cells <- pointpred_treated - pointpred_control
  treat_eff_in_cells <- epred_treated - epred_control
  
  treat_list <- list(pt_treat_eff_in_cells, treat_eff_in_cells)
  return (treat_list)
}  
# 1. Simple poststrat for a survey with 1 categorical predictor

# Create the fake world and simulate fake data


age_cat <- (1:4)
sample_N     <- c(10,     20,   30,   40)
population_N <- c(2000, 3000, 2000, 3000)

n <- 100
age_cats <- rep(age_cat, sample_N)
y <- 2 + 3*age_cats + rnorm(n, 0, 1)
fake <- data.frame(age_cats, y)

# Fit a regression model to data
fit_1 <- stan_glm(y ~ age_cats, data=fake, refresh=0)
print(fit_1) # and compare to true values of parameters

ggplot (fake, aes(x=age_cats, y=y)) + geom_jitter(width=0.05, height=0, color="blue") +
  geom_abline(intercept = coef(fit_1)[1], slope=coef(fit_1)[2], color="blue") +
  xlab ("Age Category") + ylab ("Outcome")

# Define stratification cells: age_cat  1 x 4
cells <- data.frame(age_cat=1:4)

# Prediction for each cell
pointpred <- predict          (fit_1, newdata=cells) # point prediction
epred     <- posterior_linpred(fit_1, newdata=cells) # prediction with uncertainity

# Calculate average in sample
prestrat_epred <- epred %*% sample_N / sum(sample_N) # This gives a vector of length N_sims
prestrat_pointpred <- sum(pointpred * sample_N) / sum(sample_N)

cat ("posterior_linpred: sample avg=", median(prestrat_epred), ", sd=", mad(prestrat_epred),
     " --  predict: sample average=", prestrat_pointpred)

# Poststratify!
poststrat_epred <- epred %*% population_N / sum(population_N) # This gives a vector of length N_sims
poststrat_pointpred <- sum(pointpred * population_N) / sum(population_N)

cat ("posterior_linpred: population avg=", median(poststrat_epred), ", sd=", mad(poststrat_epred),
     " --  predict: sample average=", poststrat_pointpred)

# 2. Possible generalizations:
# - causal inference
# - nonlinear model & interactions
# - continuous predictor
# - adding other predictors in the poststrat
# - goal is to estimate a regression, not just a population avg

# 3. Causal inference in simplest case

# 3a. Potential outcomes are implicit
age_cat <- (1:4)
sample_N     <- c(10, 20, 30, 40)
population_N <- c(2000, 3000, 2000, 3000)

n <- 100
age_cats <- rep(age_cat, sample_N)
z <- sample(c(0,1), n, replace=TRUE)
y <- 2 + 3*age_cats + 2*z + rnorm(n, 0, 1)
fake <- data.frame(age_cats, z, y)

fit_2 <- stan_glm(y ~ age_cats + z, data=fake, refresh=0)
print(fit_2) # and compare to true values of parameters

ggplot (fake, aes(x=age_cats, y=y, color=factor(z)))+ 
  geom_jitter(width=0.05, height=0) +
  geom_abline(intercept = coef(fit_2)[1], slope=coef(fit_2)[2], color="blue") +
  geom_abline(intercept = coef(fit_2)[1] + coef(fit_2)[3], slope=coef(fit_2)[2], color="red") +
  scale_color_manual(name = "Treatment", values = c("blue", "red")) +
  xlab ("Age Category") + ylab ("Outcome")

# Define cells
cells_2 <- data.frame(age_cats=1:4)

# Prediction for each cell.  Trated and Control *Separately*
treat_list <- fill_strat_treat_cells (fit_2, cells_2)

pt_treat_eff_in_cells <- treat_list[[1]]
treat_eff_in_cells    <- treat_list[[2]]

# Prestratify
prestrat_treat_eff <- (treat_eff_in_cells %*% sample_N) / sum(sample_N)
cat("Sample average:", median(prestrat_treat_eff), ", sd:", mad(prestrat_treat_eff))

# Poststratify!
poststrat_treat_eff <- (treat_eff_in_cells %*% population_N)/sum(population_N)
cat("Population average:", median(poststrat_treat_eff), ", sd:", mad(poststrat_treat_eff))

pt_poststrat_treat_eff <- (pt_treat_eff_in_cells %*% population_N)/sum(population_N)
cat ("Population average:", pt_poststrat_treat_eff)

# 4. Casual inference, estimating an avg regression relationship
age_cat  <- (1:4)
male_cat <- c(0,1)
treat    <- c(0,1)

sample_N     <- c(  10,   10,   20,   20,   30,   30,  40,   40)
population_N <- c(1000, 1000, 1500, 1500, 1100, 900, 2000, 1000)

n <- 100
age_cats <- rep(rep(age_cat, each=2), sample_N)
males    <- sample(male_cat, n, prob=c(0.7, 0.3), replace=TRUE)
z        <- sample(treat, n, replace=TRUE)

y <- 2 + 3*age_cats + 2*z + 1*males + rnorm(n, 0, 1)
fake <- data.frame(age_cats, z, males, y)

fit_3 <- stan_glm(y ~ age_cats + z + males, data=fake, refresh=0)
print(fit_3) # and compare to true values of parameters

ggplot (fake, aes(x=age_cats, y=y, color=factor(z)))+ 
  geom_jitter(width=0.05, height=0) +
  geom_abline(intercept = coef(fit_3)[1], slope=coef(fit_3)[2], color="blue") +
  geom_abline(intercept = coef(fit_3)[1] + coef(fit_3)[3], slope=coef(fit_2)[2], color="red") +
  geom_abline(intercept = coef(fit_3)[1] + coef(fit_3)[4], slope=coef(fit_3)[2], color="blue", linetype = "dashed") +
  geom_abline(intercept = coef(fit_3)[1] + coef(fit_3)[3] + coef(fit_3)[4], slope=coef(fit_2)[2], color="red", linetype = "dashed") +
  scale_color_manual(name = "Treatment", values = c("blue", "red")) +
  xlab ("Age Category") + ylab ("Outcome") + ggtitle ("Outcome -  Women: solid, Men: dashed")

cells_3 <- data.frame(age_cats=c(1,1, 2,2, 3,3, 4,4), males=c(0, 1, 0, 1, 0, 1, 0, 1))

# Prediction for each cell.  Trated and Control *Separately*
treat_list <- fill_strat_treat_cells (fit_3, cells_3)

pt_treat_eff_in_cells <- treat_list[[1]]
treat_eff_in_cells    <- treat_list[[2]]

# Prestratify
prestrat_treat_eff <- (treat_eff_in_cells %*% sample_N) / sum(sample_N)
cat("Sample average:", median(prestrat_treat_eff), ", sd:", mad(prestrat_treat_eff))

# Poststratify!
poststrat_treat_eff <- (treat_eff_in_cells %*% population_N) / sum(population_N)
cat("Population average:", round (median(poststrat_treat_eff),2), ", sd:", mad(poststrat_treat_eff))

pt_poststrat_treat_eff <- (pt_treat_eff_in_cells %*% population_N)/sum(population_N)
cat ("Population average:", pt_poststrat_treat_eff)

# constructing cells if z weren't treated separately: 4 x 2 x 2 = 16 rows
age_cats <- rep (age_cat, each=(length(male_cat))*(length(treat)))
males    <- rep (rep(male_cat, each=(length(treat)), length(age_cat)))
z        <- rep (treat, length(male_cat)*length(age_cat))

cells <- data.frame (age_cat=age_cats, males=males, z=z)

# suppose we wanted to estimate the regression of y on age_cat and z, averaging over sex?

#
# Chapeter 21
# Estimatimg causal effects indirectly
#
sesame <- read.dta("sesame.dta")
#' **Load data**

sesame$watched <- ifelse(sesame$viewcat==1, 0, 1)
sesame$encouraged <- ifelse(sesame$viewenc==2, 0, 1)
sesame$y       <- sesame$postlet   # letter recogintion outcome   
sesame$pretest <- sesame$prelet    # letter recognition pretest   

(sesame_tab <- table(sesame[,c('watched','encouraged')]))
round(prop.table(sesame_tab, margin=2), digits=2)

#' Estimate the intent-to-treat (ITT) effect of the instrument
#' (encouragement) on the treatment (regular watching), that is,
#' percentage of children actually induced to watch Sesame Street by
#' the intervention**
itt_zt <- stan_glm (watched ~ encouraged, data=sesame)
print (itt_zt, digits=2) # encouraged = 0.36(0.05)

#' Estimate the intent-to-treat (ITT) estimate on the outcome
#' (post-treatment letter identification)
itt_zy <- stan_glm(postlet ~ encouraged, data=sesame)
print(itt_zy, digits=2) # encouraged = 2.9(1.8)

# Wald estimator
wald_est <- coef (itt_zy)["encouraged"]/coef (itt_zt)["encouraged"] # 8 point
# The estimated effect of regularly watching sesame street 
# is 8 points on the letter recognition test

#
# TWO -STAGE least squares
#

# start same as before:
fit_2a <- stan_glm (watched ~ encouraged, data=sesame)
print (fit_2a, digits=2) # encouraged = 0.36(0.05)
summary(fit_2a$fitted, digits=2)

sesame$watched_hat <- fit_2a$fitted
         
fit_2b <- stan_glm (postlet ~ watched_hat, data=sesame)
print (fit_2b, digits=2) # watched_hat 8!! (4.8)  - se is not correct!!!

# Add other predictors
fit_3a <- stan_glm (watched ~ encouraged + prelet + as.factor(site) + setting, data=sesame)
print (fit_3a, digits=2) # encouraged = 0.34(0.1)

sesame$watched_hat_3 <- fit_3a$fitted

fit_3b <- stan_glm (postlet ~ watched_hat_3 + prelet + as.factor(site) + setting, data=sesame)
print (fit_3b, digits=2) # watched_hat 14 (3.7)  - se is not correct!!!

# USE BAYESIAN BRMS
f1 <- bf(watched ~ encour)
f2 <- bf(postlet ~ watched)
IV_brm_a <- brm(f1 + f2, data=sesame, rescor = T)
print(IV_brm_a, digits=1) # postlet_watched 8(4.8)

# Add other predictors
f1 <- bf(watched ~ encour + prelet + setting + factor(site))
f2 <- bf(postlet ~ watched + prelet + setting + factor(site))
IV_brm_b <- brm(f1 + f2, data=sesame)
print(IV_brm_a, digits=1) # postlet_watched 8(4.8)




