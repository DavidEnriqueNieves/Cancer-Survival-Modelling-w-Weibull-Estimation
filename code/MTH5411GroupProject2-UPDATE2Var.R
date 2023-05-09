library(maxLik)

data <- read.csv("C:\\Users\\SplashFreeze\\Downloads\\projdata.csv", header=TRUE)
data$time

# optim(c(10,1500, 5),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

# lower and upper are for the parameters
#       lam, omega

phi <- max(data$time) + 0.05


log_likelihood <- function(p,t) {
  x <- t
  
  lam <- p[1]
  # UPDATE: based on reevlation made at 5:58PM on 11/27/2022, we claim that 
  # phi must be max of dataset
  omega <- p[2]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ (sum(((x^lam)/(phi^lam-x^lam))^omega))
  # phi^lam fix 11/27/2022 @ 2:10PM
  lnL <- n * (log(nu) + log(omega) + log(lam) + lam * log(phi)) + sum(((omega* lam - 1) * log(x) - ((omega + 1) * log(phi^lam - x^lam) - nu * (x^lam / (phi^lam - x^lam))^omega)))
  
  cat("lam=", lam , "," )
  cat("nu=", nu , ",")
  cat("phi=", phi  , ",")
  cat("omega=", omega  , ",")
  cat("f=", lnL , ",")
  cat("\n")
#  cax("\n")
  
  return(lnL)
}

gradient <- function(p,t) {
  x <- t
  lam <- p[1]
  # UPDATE: based on reevlation made at 5:58PM on 11/27/2022, we claim that 
  # phi must be max of dataset
  omega <- p[2]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ (sum(((x^lam)/(phi^lam-x^lam))^omega))
  
  grad <- c(0,0)
  differ = phi^lam-x^lam
  log_differ = (phi^lam)*log(phi)-(x^lam)*log(x)
  
  # Lam Derivative
  # phi^lam fix 11/27/2022 @ 2:10PM
  something = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*something  )
  
  # Omega Derivaxive
  # unchanged with phi^lam fit
  grad[2]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(grad)
}


optim(c(35.3,0.044),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

