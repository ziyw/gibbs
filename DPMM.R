
# Data parameters
D <- 2
N <- 100
K_true <- 4 

# Generate data 
mu_scale <- 4.0 
cov_scale <- 0.7

z_true <- sample(1:K_true, N, replace=TRUE)

mu <- matrix(rnorm(D * K_true), D, K_true) * mu_scale

X <- mu[,z_true] + matrix(rnorm( D * N ), D, N) * cov_scale 
X <- t(X)

plot(X[,1], X[,2], col = z_true)




# gibbs sampling 

library("MCMCpack")
library("mvtnorm")
library("LaplacesDemon")

get_normalized <- function(x) {x / sum(x)}
get_index <- function(x){which(x == 1)[1]}

# Iteration times 
iter_time <- 100

# size of samples and dimension 
N <- nrow(X)
D <- 2
K <- 4

# parameters for IW: m_0, lambda_0, nu_0,S_0
m_0 <- rep(0,D)
lambda_0 <- cov_scale^2/(mu_scale^2)
nu_0 <- D + 3
S_0 <- cov_scale^2*nu_0*diag(D)

# mu, mu_k = mu[k,]
mu <- array( rep(0, K*D), dim = c(K,D))
# sigma, sigma_k = sigma[,,k]
sigma <- array( rep(0, K*D*D), dim = c(D,D,K))

for (k in c(1:K)){
  sigma[,,k] <- diag(D)
}

# initilize 
z <- sample(1:K, N, replace=TRUE)

# init mu_k
for(k in c(1:K)){
  
  N_k <- sum(z == k)
  
  if(N_k == 1){
    x_bar <- X[z == k,]
  } else if(N_k == 0){
    x_bar <- rep(0,D)
  } else {
    x_bar <- colSums(X[z == k,]) / N_k    
  }
  mu[k,] <- x_bar 
}

# init sigma_k
for( k in c(1:K)){
  
  N_k <- sum(z == k)
  
  if(N_k == 1){
    x_bar <- X[z == k,]
    sigma[,,k] <- var(x_bar)
  } else if(N_k == 0){
    next
  } else {
    x_bar <- X[z == k,]
    sigma[,,k] <- var(x_bar)
  }
  
}

alpha <- rep(1.0, K)
pi <- rep(1.0/K, K)

# gibbs sampler 

for (iter in c(1:iter_time)) {
  print("iteration times:")
  print(iter)
  
  for (i in c(1:N)){
    
    # sample process
    # parameters: z, pi, mu, sigma
    
    # sample z 
    prob <- rep(0, K)
    
    for (k in c(1:K)){
      prob[k] <- pi[k] * dmvnorm(X[i,],mu[k,],sigma[,,k])
    }
    
    new_k <- rmultinom(1, 1, prob = prob)
    z[i] <- get_index(new_k)
    # sample pi 
    
    for( k in c(1:K)){
      alpha[k] <- alpha[k] + sum(z == k)
    }
    
    pi <- rdirichlet(1, alpha)
    
    for (k in c(1:K)){
      
      n_k <- sum(z == k)
      x_bar <- colMeans(X[z == k,])
      
      mu_k  <- ( lambda_0 * m_0 + n_k * x_bar ) / (lambda_0 + n_k)
      lambda_k <- lambda_0 + n_k
      nu_k <- nu_0 + n_k
      S_k <- S_0 + (lambda_0 * n_k) / (lambda_0 + n_k) * (x_bar - m_0) %*% t(x_bar - m_0)
      
      S <- array(rep(0,D*D),dim = c(D,D))
      if(n_k !=0 ){
        x_sample <- X[z == k,]
        for(i in c(1:n_k)){
          S <- S + (x_sample[i,] - x_bar) %*% t(x_sample[i,] - x_bar) 
        }
      }
      
      S_k <- S_k + S
      sample <- rnorminvwishart(n=1, mu_k, lambda_k, S_k, nu_k)
      mu[k,] <- sample$mu
      sigma[,,k] <- sample$Sigma 
    }
  }
}


plot(X[,1],X[,2],col = z)

