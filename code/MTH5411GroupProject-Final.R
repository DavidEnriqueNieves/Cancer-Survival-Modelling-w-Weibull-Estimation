# Working with maxLik package
# If this is an error, then the package needs to be installed
library(maxLik)


# Importing sample data for optimization
t<- c(0)
filename <- "C:\\Users\\SplashFreeze\\Downloads\\projdata.csv"
data <- read.csv(filename, header=TRUE)

df_male <- data[data$sex == '1',]
df_female <- data[data$sex == '2',]
df_male
t_male <- df_male$time
t_female <- df_female$time
t_male
t_female
t <- data$time

# Gathering initial set of values for each parameter
# Not all of these may be used
init_lambda = 40.1
init_omega =  0.035
init_phi = max(data$time) + 0.0005
init_nu = 6

# Set this to true/false if you want printed function values for each iteration
print_results = TRUE

log_likelihood_maxlik <- function(p,x) {
  # Initialize values
  x <- t
  lam <- p[1]
  omega <- p[2]
  phi <- init_phi
  n <- length(x)
  
  # Calculate nu
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  # Calculate log likelood function
  lnL <- n * (log(nu) + log(omega) + log(lam) + lam * log(phi)) + sum(((omega* lam - 1) * log(x) - ((omega + 1) * log(phi^lam - x^lam) - nu * (x^lam / (phi^lam - x^lam))^omega)))
  
  # Print results
  if (print_results){
  cat("lam=", lam , "," )
  cat("nu=", nu , ",")
  cat("phi=", phi  , ",")
  cat("omega=", omega  , ",")
  cat("f=", lnL , ",")
  cat("\n")
  }
  
  return(lnL)
}

gradient_maxlik <- function(p,x) {
  # Initialize values
  x <- t
  lam <- p[1]
  omega <- p[2]
  phi <- init_phi
  n <- length(x)
  
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  
  grad <- c(0,0)
  differ = phi^lam-x^lam
  log_differ = (phi^lam)*log(phi)-(x^lam)*log(x)
  
  # Lam Derivative
  sub_value = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  # Calculate log likelihood gradient for lambda
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*sub_value  )
  
  # Omega Derivative
  # Calculate log likelihood gradient for omega
  grad[2]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(grad)
}

log_likelihood_optim <- function(p,t) {
  x<- t
  result <- log_likelihood_maxlik(p,x)
  return(-1 * result)
}

gradient_optim <- function(p,t) {
  # Initialize values
  x <- t
  lam <- p[1]
  omega <- p[2]
  phi <- init_phi
  n <- length(x)
  
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  
  grad <- c(0,0)
  differ = phi^lam-x^lam
  log_differ = (phi^lam)*log(phi)-(x^lam)*log(x)
  
  # Lam Derivative
  sub_value = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  # Calculate log likelihood gradient for lambda
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*sub_value  )
  
  # Omega Derivative
  # Calculate log likelihood gradient for omega
  grad[2]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(grad)
}

ml <- maxLik(logLik =log_likelihood_maxlik, grad=gradient_maxlik, start=c(init_lambda, init_omega), control=list(printLevel=2))
ml

optim(c(init_lambda, init_omega),log_likelihood_optim, gr=gradient_optim, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

