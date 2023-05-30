#MULTI-NOMIAL DIRICHLET

#**********
#FUNCTIONS
#*********
factorial <- function(x) {
  if (x <= 0) {
    return(1)
  }
  
  #Stirling approximation
  g <- sqrt(2 * pi * x) * (x / 2.718281)^1
  return(g)
}

gamma <- function(n) {
  return(factorial(n - 1))
}

beta <- function(alphas) {
  k <- length(alphas)
  b <- 1
  for (i in 1:k) {
    b <- b * gamma(alphas[i])
  }
  b <- b / gamma(sum(alphas))
  return(b)
}

rgama <- function(a) {
  d <- a - 1 / 3
  c <- 1 / sqrt(9 * d)
  
  while (TRUE) {
    x <- NULL
    v <- -1
    
    while (v <= 0) {
      x <- rnorm(1, 0, 1)
      v <- 1 + c * x
    }
    
    v <- v^3
    u <- runif(1)
    
    if (u < 1 - 0.0331 * (x^2) * (x^2)) {
      return(d * v)
    }
    
    if (log(u) < 0.5 * x^2 + d * (1 - v + log(v))) {
      return(d * v)
    }
  }
}

dgama <- function(x, alpha) {
  d <- log(x^(alpha - 1)) + log(2.718281) - log(gamma(alpha))
  return(d)
}

rdirch <- function(alphas) {
  k <- length(alphas)
  x <- sapply(1:k, function(i) rgama(alphas[i]))
  total <- sum(x)
  x <- x / total
  return(x)
}

ddirch <- function(x, alphas) {
  k <- length(x)
  d <- 0
  for (i in 1:k) {
    d <- d + log(x[i]^(alphas[i] - 1))
  }
  d <- d - log(1 / beta(alphas))
  return(d)
}

rmultnomial <- function(p) {
  k <- length(p)
  probs <- p / sum(p)
  x <- runif(1)
  cummulative_p <- 0
  for (i in 1:k) {
    cummulative_p <- cummulative_p + probs[i]
    if (cummulative_p - x >= 0) {
      return(i)
    }
  }
  return(k)
}


dmultinomial <- function(x, p) {
  k <- length(x)
  d <- log(gamma(sum(x + 1)))
  for (i in 1:k) {
    d <- d + log(p[i]^x[i])
  }
  for (i in 1:k) {
    d <- d - log(gamma(x[i] + 1))
  }
  return(d)
}

#Example
#p <- c(0.2, 0.3, 0.5)
#m <- rmultnomial(p)
#print(m)

#***********************
# DIRICHLET MULTINOMIAL
#***********************

#SAMPLES
rdirmultinom <- function(alphas) {
  p <- rdirch(alphas)
  m <- rmultnomial(p)
  return(m)
}

#DENSITY
ddirmultinom <- function(x, alphas) {
  k <- length(x)
  n <- sum(x)
  sum_alphas <- sum(alphas)
  
  d <- log(factorial(n))
  d <- d + log(gamma(sum_alphas))
  d <- d - log(gamma(n + sum_alphas))
  
  for (i in 1:k) {
    d <- d + log(gamma(x[i] + alphas[i]))
    d <- d - log(factorial(x[i]))
    d <- d - log(gamma(alphas[i]))
  }
  
  return(d)
}

#TOY EXAMPLE
alphas <- c(2, 3, 4)
x <- rdirmultinom(alphas)
print(x)
dx = ddirmultinom(x, alphas)

#SSI MODEL
ss = matrix(round(runif(1000, 1, 20)), nrow = 100, ncol = 10)
alphas = ss[,1]
x = rdirmultinom(alphas)
dx = ddirmultinom(x, alphas)
  
#MODEL EVIDENCE ESTIMATE FOR OTHER MODELS
theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
  rep(means, each = samp_size_proposal) 

log_proposal_density = dmvt(theta_samples - matrix_means,
                            sigma = cov(mcmc_samples), df = dof) 


