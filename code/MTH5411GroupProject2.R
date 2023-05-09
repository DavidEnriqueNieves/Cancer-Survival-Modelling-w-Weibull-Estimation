library(maxLik)
log_likelihood <- function(p,x) {
  
  lam <- p[1]
  phi <- p[2]
  omega <- p[3]
  n <- length(t)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((t^lam)/(phi^lam-t^lam))^omega)
  # phi^lam fit 11/27/2022 @ 2:10PM
  lnL <- n * (log(nu) + log(omega) + log(lam) + lam * log(phi)) + sum(((omega* lam - 1) * log(t) - ((omega + 1) * log(phi^lam - t^lam) - nu * (t^lam / (phi^lam - t^lam))^omega)))
  
  cat("lam=", lam , "," )
  cat("nu=", nu , ",")
  cat("phi=", phi  , ",")
  cat("omega=", omega  , ",")
  cat("f=", lnL , ",")
  cat("\n")
#  cat("\n")
  
  return(-lnL)
}

log_likelihood_maxLik <- function(p,x) {
  lam <- p[1]
  phi <- p[2]
  omega <- p[3]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  # phi^lam fit 11/27/2022 @ 2:10PM
  lnL <- n * (log(nu) + log(omega) + log(lam) + lam * log(phi)) + sum(((omega* lam - 1) * log(x) - ((omega + 1) * log(phi^lam - x^lam) - nu * (x^lam / (phi^lam - x^lam))^omega)))
  
  cat("lam=", lam , "," )
  cat("nu=", nu , ",")
  cat("phi=", phi  , ",")
  cat("omega=", omega  , ",")
  cat("f=", lnL , ",")
  cat("\n")
#  cat("\n")
  
  return(-lnL)
}
gradient <- function(p,t) {
  x <- t
  lam <- p[1]
  phi <- p[2]
  omega <- p[3]
  n <- length(t)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((t^lam)/(phi^lam-t^lam))^omega)
  
  grad <- c(0,0,0)
  differ = phi^lam-t^lam
  log_differ = (phi^lam)*log(phi)-(t^lam)*log(t)
  
  # Lam Derivative
  # phi^lam fit 11/27/2022 @ 2:10PM
  something = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  grad[1]= n/lam + n * log(phi) + sum(omega*log(t) - (omega+1)*(log_differ)/(differ) - nu*something  )
  
  # Phi derivative
  # phi^lam fit 11/27/2022 @ 2:10PM
  # tODO CHECK tHEORY??
  grad[2]=  n * lam / phi + sum(-(omega + 1) * (1/(phi^lam - t^lam))) * lam * phi^(lam - 1) + lam * nu * omega * t^(lam * omega) * (phi^lam - t^lam)^(-omega - 1) * phi^(lam - 1)
  
  # Omega Derivative
  # unchanged with phi^lam fit
  grad[3]=  n/omega + sum(lam*log(t) - log(differ) - nu*((t^lam)/(differ))^omega * log((t^lam)/(differ)))
  
  return(-grad)
}

gradient_maxLik <- function(p,x) {
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
  # phi^lam fit 11/27/2022 @ 2:10PM
  something = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (x^lam * log(x) * differ^(-1) - x^lam * log_differ/(differ)^2)
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*something  )
  
  # Phi derivative
  # phi^lam fit 11/27/2022 @ 2:10PM
  # tODO CHECK tHEORY??
  grad[2]=  n * lam / phi + sum(-(omega + 1) * (1/(phi^lam - x^lam))) * lam * phi^(lam - 1) + lam * nu * omega * x^(lam * omega) * (phi^lam - x^lam)^(-omega - 1) * phi^(lam - 1)
  
  # Omega Derivative
  # unchanged with phi^lam fit
  grad[3]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(-grad)
}


data <- read.csv("C:\\Users\\SplashFreeze\\Downloads\\projdata.csv", header=TRUE)
data$time
max(data$time)

# optim(c(10,1500, 5),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")
maxLik(start=c(10, 1500, 5), logLik=log_likelihood_maxLik, grad=gradient_maxLik, control=list(printLevel=2))
# ml

# print(summary(ml))

