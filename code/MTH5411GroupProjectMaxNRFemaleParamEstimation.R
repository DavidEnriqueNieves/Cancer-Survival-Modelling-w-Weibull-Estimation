# working with maxLik
library(maxLik)
library(ramify)

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
t <- t_female

# optim(c(10,1500, 5),log_likelihood, gr=gradient, control=list(maxit=100), t=data$time, method="CG") # , lower=c(0,0.1,1022,0.1),upper=c(1000,1000, 1200,1000), method="L-BFGS-B")

# lower and upper are for the parameters
#       lam, omega:w

phi <- max(data$time) + 0.0005

last_params <<- c(0,0,0)

log_likelihood <- function(p,x) {
  x <- t
  
  lam <- p[1]
  # UPDATE: based on reevlation made at 5:58PM on 11/27/2022, we claim that 
  # phi must be max of dataset
  phi <- p[2]
  omega <- p[3]
  # cat("last_params = " , last_params , "\n")
  last_params <<- p
  n <- length(x)
  # tODO CHeck that log evaluates to ln
  nu <- n/ sum(((x^lam)/(phi^lam-x^lam))^omega)
  # phi^lam fix 11/27/2022 @ 2:10PM
  lnL <- n * (log(nu) + log(omega) + log(lam) + lam * log(phi)) + sum(((omega* lam - 1) * log(x) - ((omega + 1) * log(phi^lam - x^lam) - nu * (x^lam / (phi^lam - x^lam))^omega)))
  last_loglikelihood <<- lnL
  
   cat("lam=", lam , "," )
   cat("nu=", nu , ",")
   cat("phi=", phi  , ",")
   cat("omega=", omega  , ",")
   cat("f=", lnL , ",")
   cat("\n")
    cat("\n")
  
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
n <- 10000

bounds = c(1,500)
lam_space <- seq(from=bounds[1], to=bounds[2], by=((bounds[2]-bounds[1]) / n))

bounds = c(1023,1500)
phi_space <- seq(from=bounds[1], to=bounds[2], by=((bounds[2]-bounds[1]) / n))

bounds = c(0, 1)
omega_space <- seq(from=bounds[1], to=bounds[2], by=((bounds[2]-bounds[1]) / n))
last_params
omega_hist_vector <- c(0)


last_params
lambda_hist_vector <- c(0)
loglik_hist_vector <- c(0)

successes <- 0
errors <- 0
warnings <- 0

last_loglikelihood_hist_vector <- 0
for(lambda in lam_space){
    out <- tryCatch(
    {
      ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(lambda,1500,  0.02),iterlim=50,tol=0.0001, control=list(printLevel=2))
      last_params
      print("Successful return!")
      successes <<- successes +1
      return(last_params[1])
    },
    error=function(cond){
      print("Error detected!")
      errors <<- errors + 1
      return(last_params[1])
    },
    warning=function(cond){
      print("Warning emitted!")
      print('last params are ')
      print(last_params)
      print(cond)
      warnings <<- warnings + 1
      return(last_params[1])
    }
  )
    cat("last_params is \n")
    cat(out)
  lambda_hist_vector <- c(lambda_hist_vector, out)
  last_loglikelihood_hist_vector <- c(last_loglikelihood_hist_vector, last_loglikelihood)
}

successes
errors
warnings
lambda_hist_vector
hist(lambda_hist_vector)
hist(last_loglikelihood_hist_vector)

omega_hist_vector <- c(0)
for(omega in omega_space){
    out <- tryCatch(
    {
      ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(600,1500,omega),iterlim=50,tol=0.0001, control=list(printLevel=2))
      last_params
      print("Successful return!")
      successes <<- successes +1
      return(last_params[3])
    },
    error=function(cond){
      print("Error detected!")
      errors <<- errors + 1
      return(last_params[3])
    },
    warning=function(cond){
      print("Warning emitted!")
      print('last params are ')
      print(last_params)
      print(cond)
      warnings <<- warnings + 1
      return(last_params[3])
    }
  )
    cat("last_params is \n")
    cat(out)
  omega_hist_vector <- c(omega_hist_vector, last_params[3])
  last_loglikelihood_hist_vector <- c(last_loglikelihood_hist_vector, last_loglikelihood)
}

successes
errors
warnings
lambda_hist_vector
hist(omega_hist_vector)
hist(last_loglikelihood_hist_vector)


ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(lambda,1500,  0.02),iterlim=50,tol=0.0001, control=list(printLevel=2))
last_params
last_loglikelihood



mle_lambda_hist <- c(0)


lambda_hist_vector <- c(0)
for(lambda in lam_space){
  lnL <- log_likelihood( c(lambda,1022.01,0.023 ),t)
  if(!(is.na(lnL))){
    mle_lambda_hist <<- c(mle_lambda_hist, lnL)
    lambda_hist_vector<- c(lambda_hist_vector, lambda)
  }
}


hist(mle_lambda_hist)
hist(lambda_hist_vector)
# plot(lambda_hist_vector, mle_lambda_hist, ylim = c(-1700, -1050))
# plot(lambda_hist_vector, mle_lambda_hist, ylim=c(-1000, -600), xlim = c(0,90))
plot(lambda_hist_vector, mle_lambda_hist, ylim=c(-1500, -1050), xlim = c(0,90))

# plot(lambda_hist_vector, mle_lambda_hist, ylim=c(-600, -415), xlim = c(0,95))

max_loglik <- -9999
max_loglik_index <- 0
index <- 1
for(loglik in mle_lambda_hist){
  cat(loglik, "\n")
  if(loglik > max_loglik && loglik != 0){
    max_loglik <<- loglik
    max_loglik_index <- index
  }
  index <- index + 1
}
max_loglik
max_loglik_index
lambda_hist_vector[max_loglik_index]



max_loglik[1]
max_loglik_index <- which(mle_lambda_hist == max_loglik)
max_loglik_index[1]
max_lambda <- mle_lambda_hist[max_loglik_index]
max_lambda


bounds = c(0,5)
omega_space <- seq(from=bounds[1], to=bounds[2], by=((bounds[2]-bounds[1]) / n))
last_params
omega_hist_vector <- c(0)
omega_hist_vector <- c(0)
mle_omega_hist <- c(0)
omega_hist_vector <- c(0)
for(omega in omega_space){
  lnL <- log_likelihood( c(48.3,1022.01,omega),t)
  if(!(is.na(lnL))){
    mle_omega_hist <<- c(mle_omega_hist, lnL)
    omega_hist_vector<- c(omega_hist_vector, omega)
  }
}
hist(mle_omega_hist)
hist(omega_hist_vector)
plot(omega_hist_vector, mle_omega_hist, ylim=c(-2000, -1000))
plot(omega_hist_vector, mle_omega_hist)
# plot(omega_hist_vector, mle_omega_hist, ylim=c(-1400,-600))
# plot(omega_hist_vector, mle_omega_hist, ylim=c(-900,-420))

max_loglik <- -9999
max_loglik_index <- 0
index <- 1
for(loglik in mle_omega_hist){
  cat(loglik, "\n")
  if(loglik > max_loglik && loglik != 0){
    max_loglik <<- loglik
    max_loglik_index <- index
  }
  index <- index + 1
}
max_loglik
max_loglik_index
omega_hist_vector[max_loglik_index]


 # gathers -404
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(35.5,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
last_params

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(35.5,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))


log_likelihood(last_params,t)
# gathers -402.4
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(38.4,965.5,  0.0382 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(38.4,965.001,  0.0382 ),iterlim=20,tol=0.0001, control=list(printLevel=2))
# -402.55
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(38.4,965.001,  0.03819),iterlim=20,tol=0.0001, control=list(printLevel=2))
# -404.21
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(38.41,965.001,  0.03819),iterlim=20,tol=0.0001, control=list(printLevel=2))
# gathres -404
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(38.40062,965.5,  0.0379 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 
# # gathers -403.6
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(40,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))


# gathrs -409
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# gathers -417.9
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(25,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))


# gathrs -432
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(20,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# gathrs -417
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(0.5,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# gathers -410.9
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(29,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# get NaNs due to overflows
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(100,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# also get NaNs 
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(50,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# gets -402.7
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(45,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))

# NaNs
ml <- maxLik(logLik =log_likelihood, grad=gradient, start=c(75,1500,  0.035 ),iterlim=20,tol=0.0001, control=list(printLevel=2))


# varying initial guesses for omega with lambda = 30

# -423
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500,  1),iterlim=20,tol=0.0001, control=list(printLevel=2))

# NaNs
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 10),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -424.7866
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 2),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# -476
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1500, 1.5),iterlim=20,tol=0.0001, control=list(printLevel=2))

# gets to 0.07520
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5, 1.5),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 0.07459
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.5),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 0.0751, 
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.25),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 0.0791
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.10),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 0.06871

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.0721),iterlim=20,tol=0.0001, control=list(printLevel=2))

# 
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(30,1022.5,0.06871),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(42.97,1500,0.023),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(45.1,1500,0.023),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.0245),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -616.8963
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.02451914),iterlim=20,tol=0.0001, control=list(printLevel=2))


# goes back to -632
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.98582,1500,0.0272111),iterlim=20,tol=0.0001, control=list(printLevel=2))
# 
# 
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.5122,1500,0.0245191),iterlim=20,tol=0.1, control=list(printLevel=2))

# -709
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(15,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))



# -615.9547
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(48,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

# -614
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,1500,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,1000,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,965.5,0.0245191),iterlim=20,tol=0.0001, control=list(printLevel=2))
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,965.5,0.02444),iterlim=20,tol=0.0001, control=list(printLevel=2))

ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,965.1,0.02439),iterlim=20,tol=0.0001, control=list(printLevel=2))
# starts increasing, -410
ml <- maxLik( logLik =log_likelihood, grad=gradient, start=c(47.9,965.1,0.02435),iterlim=20,tol=0.0001, control=list(printLevel=2))

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

