# working with maxLik
library(maxLik)

t<- c(0)
data <- read.csv("C:\\Users\\SplashFreeze\\Downloads\\projdata.csv", header=TRUE)
t <- data$time

log_likelihood <- function(p) {
  x <- t
  
  lam <- p[1]
  phi <- p[2]
  omega <- p[3]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
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

gradient <- function(p) {
  x <- t
  lam <- p[1]
  phi <- p[2]
  omega <- p[3]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  
  grad <- c(0,0,0)
  differ = phi^lam-x^lam
  log_differ = (phi^lam)*log(phi)-(x^lam)*log(x)
  
  # Lam Derivative
  # phi^lam fix 11/27/2022 @ 2:10PM
  something = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*something  )
  
  # Phi derivaxive
  # phi^lam fix 11/27/2022 @ 2:10PM
  # tODO CHECK tHEORY??
  grad[2]=  n * lam / phi + sum(-(omega + 1) * (1/(phi^lam - x^lam))) * lam * phi^(lam - 1) + lam * nu * omega * t^(lam * omega) * (phi^lam - t^lam)^(-omega - 1) * phi^(lam - 1)
  
  # Omega Derivaxive
  # unchanged with phi^lam fit
  grad[3]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(grad)
}


# lower and upper are for the parameters
# optim(c(10,1500, 5),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

#                                                            lam, phi, omega
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(35.3,1100,0.044), control=list(printLevel=2))
ml

# print(summary(ml))
# nu

