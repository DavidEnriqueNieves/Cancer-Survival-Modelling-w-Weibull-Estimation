# working with maxLik
library(maxLik)

t<- c(0)
data <- read.csv("C:\\Users\\SplashFreeze\\Downloads\\projdata.csv", header=TRUE)
t <- data$time

df_male <- data[data$sex == '1',]
df_female <- data[data$sex == '2',]
df_male
t_male <- df_male$time
t_female <- df_female$time
t_male
t_female
# t <- data$time
max(t_female)
t <- t_male

# optim(c(10,1500, 5),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

# lower and upper are for the parameters
#       lam, omega:w

phi <- max(data$time) + 0.0005


log_likelihood <- function(p,x) {
  x <- t
  
  lam <- p[1]
  # UPDATE: based on reevlation made at 5:58PM on 11/27/2022, we claim that 
  # phi must be max of dataset
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

gradient <- function(p,x) {
  x <- t
  lam <- p[1]
  # UPDATE: based on reevlation made at 5:58PM on 11/27/2022, we claim that 
  # phi must be max of dataset
  phi <- p[2]
  omega <- p[3]
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  
  grad <- c(0,0)
  differ = phi^lam-x^lam
  log_differ = (phi^lam)*log(phi)-(x^lam)*log(x)
  
  # Lam Derivative
  # phi^lam fix 11/27/2022 @ 2:10PM
  something = omega*(x^lam/(phi^lam - x^lam))^(omega- 1) * (t^lam * log(t) * differ^(-1) - t^lam * log_differ/(differ)^2)
  grad[1]= n/lam + n * log(phi) + sum(omega*log(x) - (omega+1)*(log_differ)/(differ) - nu*something  )
  
  # Phi derivative
  # phi^lam fit 11/27/2022 @ 2:10PM
  # tODO CHECK tHEORY??
  grad[2]=  n * lam / phi + sum(-(omega + 1) * (1/(phi^lam - x^lam))) * lam * phi^(lam - 1) + lam * nu * omega * x^(lam * omega) * (phi^lam - x^lam)^(-omega - 1) * phi^(lam - 1)
  
  # Omega Derivative
  # unchanged with phi^lam fit
  grad[3]=  n/omega + sum(lam*log(x) - log(differ) - nu*((x^lam)/(differ))^omega * log((x^lam)/(differ)))
  
  return(grad)
}

#                                                            lam, phi     omega
# # gathers -618.241
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(35.5,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gathers -632
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(40,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gathrs -618.4594
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gathers -625.6156
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(25,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# 
# # gathrs -642.5122
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(20,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gathrs -1074
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(0.5,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gathers -619
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(29,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # get NaNs due to overflows
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(100,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # also get NaNs 
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(50,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gets -632 before going into NaN
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(45,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # NaNs
# ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(75,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# 
# # varying initial guesses for omega with lambda = 30
# 
# # gers NaN but finishes with 0.0217
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500,  1),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # NaNs
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 10),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gets 0.0204 but NaN later
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 2),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gets 2.351152 * 10^-2
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 1.5),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # gets to 0.07520
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5, 1.5),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # 0.07459
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.5),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # 0.0751, 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.25),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # 0.0791
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.10),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # 0.06871
# 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.0721),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.06871),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(42.97,1500,0.023),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(45.1,1500,0.023),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.0245),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# # -616.8963
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.02451914),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# 
# # goes back to -632
# ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.98582,1500,0.0272111),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# 
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.0245191),iterlim=20,tol=0.1, control=list(printLevel=2))

# -709
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(15,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))



# -615.9547
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(48,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -614
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -615
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.85,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -615.9
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.899,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -616.656
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.901,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -616
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.8,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -616
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.8,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))


# conclusion: -614 loglik, with phi=1022, lam=48.43294,nu = 4.754, omega=0.02718752