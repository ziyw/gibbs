# generate data 
#######
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
############








# library("mvtnorm")
library("mvnfast")
library("matrixStats")
get_normalized <- function(x) {x / sum(x)}
get_index <- function(x){which(x == 1)[1]}

draw <- function(prob) {
  k_uni <- runif(1)
  
  for (i in c(1: length(prob))){
    
    k_uni <- k_uni - prob[i]
    if (k_uni < 0) { 
      return(i)
    } 
  }
  return(length(prob) - 1)
}


get_log_p1 <- function(z,index,k, alpha, K){
  N <- length(z)
  z <- z[-index]
  p <- log(sum(z == k) + (alpha / K) ) - log(N + alpha - 1)

  return(p)
  
}


get_log_p2<- function(my_x, my_z, index, k, m_0, k_0, v_0, S_0){
  x_k <- my_x[index,]
  
  my_x <- my_x[-index, ]
  my_z <- my_z[-index]
  
  my_x <- my_x[my_z == k,]
  my_z <- my_z[my_z == k]
  
  N <- nrow(my_x)
  if(is.null(N)){
    N = 0 
  }
  k_n <- k_0 + N
  v_n <- v_0 + N
  
  if (N == 0){
    m_n <- k_0 * m_0
  } else {
    m_n <- k_0 * m_0 + N * colMeans(my_x)  
  }
  
  
  m_n <- m_n / k_n
  
  S_n <- S_0 + k_0 * m_0 %*% t(m_0) - k_n * m_n %*% t(m_n)
  
  S <- matrix(rep(0,D*D),D,D)
  
  if( N != 0 ){
    for ( i in c(1:N)) {
      S = S + my_x[i,] %*% t(my_x[i,])
    }
  }
  
  
  S_n <- S_n + S
  
  corr <- ((k_n + 1) / (k_n * (v_n - D + 1)) ) * S_n
  p <- dmvt(x_k, mu = m_n, sigma = corr, df = v_n - D + 1, log = TRUE)
  return(p)
  
}

# Model Parameters 
alpha <- 1.0 
K <- 4 
iteration <- 20

# init prior 
m_0 <- rep(0,D)
k_0 <- cov_scale^2/(mu_scale^2)
v_0 <- D + 3
S_0 <- cov_scale^2*v_0*diag(D)

# gibbs sampler 

# init z (randomly)
z <- sample(1:K, N, replace=TRUE)

for (it in c(1:iteration)) {
  
  for ( i in c(1:N)) {
    prob <- rep(0,K)
    
    for ( k in c(1:K)) {

      prob[k] <- get_log_p1(z, k,i, alpha, K) + get_log_p2(X, z,index = i , k = k , m_0, k_0, v_0, S_0)
      #prob[k] <- exp(prob[k])
    }
    prob <- exp(prob - logSumExp(prob))
    #prob <- get_normalized(prob)
    new_k <- rmultinom(1, 1, prob = prob)
    z[i] <- get_index(new_k)
  }
}

plot(X[,1],X[,2],col = z)