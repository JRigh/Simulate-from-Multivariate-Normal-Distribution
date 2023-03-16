#---------------------------------------
# SIMULATE FROM MULTIVARIATE NORMAL in R
#---------------------------------------

library(PerformanceAnalytics)
library(rlang)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(GGally)

#---------------
# Gibbs sampling
#---------------

# 1. Define parameters
n <- 6000
rho_12 <- 0.9

mu1 <- 0
mu2 <- 0
s1 <- 1
s2 <- 1

# 2. Gibbs sampling function to sample from Bivariate Normal distribution
mv_Gibbs_sampler <- function(n, burn, mu1, mu2, s1, s2, rho) {
  
  X <- matrix(0, n, 2) # the chain
  X[1,] <- c(mu1, mu2) # initial observation
  
  for(i in 2:n) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + (x2 - mu2) * rho * (s1/s2)
    X[i, 1] <- rnorm(1, mu1 + (x2 - mu2) * rho * (s1/s2),
                     sqrt(1-rho^2) * s1)
    
    x1 <- X[i, 1]
    X[i, 2] <- rnorm(1, mu2 + (x1 - mu1) * rho * (s2/s1), 
                     sqrt(1-rho^2) * s2)
  }
  
  b <- burn + 1
  x <- X[b:n, ]
  
  return(data.frame(x))
}

# 3. Draw a random sample by calling the function
set.seed(2023)
random.sample <- mv_Gibbs_sampler(n = n, burn = 1000, mu1 = mu1, mu2 = mu2,
                                  s1 = s1, s2 = s2, rho = rho_12)

# 4. Visualizing correlation
chart.Correlation(random.sample, histogram=TRUE, pch=19)

# 5. Perform linear regression and compute coefficient of determination
summary <- summary(lm(random.sample$X2 ~ random.sample$X1))
rho <- cor(random.sample$X1, random.sample$X2)

# 6. Plotting

htop <- ggplot(data=random.sample, aes(x=X1)) + 
  geom_histogram(aes(y=..density..), fill = "grey90", color = "black", binwidth = 0.3) + 
  stat_density(colour = "red3", geom="line", size = 1.2, position="identity", show.legend=FALSE) +
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())

blank <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(axis.ticks=element_blank(), panel.background=element_blank(), panel.grid=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot(data=random.sample, aes(x=X1, y=X2)) + 
  geom_point(size = 0.6, pch=3) + 
  geom_smooth(method = "lm",se = FALSE, color="red4") +
  geom_text(x=-2.3, y= 1.9, label = paste('rho = ', round(rho,4))) +
  geom_point(data=random.sample, 
             aes(x=X1,y=X2), 
             color='black',
             size=1)
theme_light()

hright <- ggplot(data=random.sample, aes(x=X2)) + 
  geom_histogram(aes(y=..density..), fill = "grey90", color = "black", binwidth = 0.3) + 
  stat_density(colour = "red3", geom="line", size = 1.2, position="identity", show.legend=FALSE) +
  coord_flip() + theme(axis.title.y = element_blank(),
                       panel.background = element_blank())

grid.arrange(htop, blank, scatter, hright, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))


#-----------------------
# Cholesky factorization
#-----------------------

# 1. Define parameters

n <- 5000
rho_12 <- -0.75
rho_13 <- 0.4
rho_23 <- -0.1
mu1 <- 0
mu2 <- 0
mu3 <- 2
s1 <- 2
s2 <- 1
s3 <- 1

# 2. Mean and covariance matrix (cov_12 = s1*s2*rho_12)
mu <- c(mu1, mu2, mu3) 

Sigma <- matrix(c(s1^2, s1*s2*rho_12, s1*s3*rho_13, 
                  s2*s1*rho_12, s2^2, s2*s3*rho_23,
                  s1*s3*rho_13, s2*s3*rho_23, s3^2), nrow = 3)


# 3. Cholesky function to sample from multivariate Normal distribution
mv_cholesky <- function(n, mu, Sigma) {
  
  d <- length(mu)
  Q <- chol(Sigma)
  Z <- matrix(rnorm(n*d), nrow = n, ncol = d)
  X <- Z %*% Q +  rep(1,n) %*% t(mu)
  X <- data.frame(X)
  return(X)
}

# 4. Draw a random sample by calling the function
set.seed(2023)
random.sample2 <- mv.cholesky(n = n, mu = mu, Sigma = Sigma)

# 5. Visualizing correlation
chart.Correlation(random.sample2, histogram=TRUE, pch=19)

ggpairs(random.sample2, upper = list(continuous = "density"),
        lower = list(combo = "facetdensity"))

#----
# end
#----
