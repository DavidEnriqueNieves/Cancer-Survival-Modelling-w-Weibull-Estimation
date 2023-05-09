

import math

var_lambda = 0
nu = 0
phi = 0
omega = 0
t_samples=  [0, 1, 2, 3,5]
a=0

def partial_wrt_nu(var_lambda, nu, phi, omega,t_samples):
    print(f"partial_wrt_nu")
    n = len(t_samples)
    result = n / nu 

    sum = 0
    for t in t_samples:
        print(t**var_lambda / (phi**var_lambda - t**var_lambda))
        sum+= (

            (t**var_lambda / (phi**var_lambda - t**var_lambda))**omega
        )
    result = result - sum

    return result


def partial_wrt_omega(var_lambda, nu, phi, omega,t_samples):
    print("partial_wrt_omega")
    n = len(t_samples)
    result = n / omega 

    sum = 0
    for t in t_samples:
        print(f"({phi ** var_lambda - t ** var_lambda})=")
        sum+= (
            var_lambda + math.log(t) - math.log(phi ** var_lambda - t ** var_lambda)
            - 
            nu * ((t**var_lambda / (phi**var_lambda - t**var_lambda))** omega) 
            * math.log(t**var_lambda / (phi**var_lambda - t**var_lambda))

        )
    result = result + sum
    return result

def partial_wrt_lambda(var_lambda, nu, phi, omega,t_samples):
    print("partial_wrt_lambda")
    n = len(t_samples)
    result = n / var_lambda 

    sum = 0
    for t in t_samples:
        sum+= (
            omega * math.log(t) - (omega + 1) * (1/(phi**var_lambda - t**var_lambda + 0.001)) * (phi**var_lambda * math.log(phi) - t**var_lambda * math.log(t))
            - nu * (omega * t**(var_lambda * omega) * math.log(t) * (phi**var_lambda - t**var_lambda)**(-omega) 
             + t**(var_lambda * omega) * (-omega) *
              (phi**var_lambda - t**var_lambda)**(-omega - 1) * 
              (phi**var_lambda * math.log(phi) - t**var_lambda * math.log(t))
            )

        )
    result = result + sum
    return result

def partial_wrt_phi(var_lambda, nu, phi, omega,t_samples):
    print("partial_wrt_phi")
    n = len(t_samples)
    result = 0

    sum = 0
    for t in t_samples:
        sum+= (
            var_lambda * -(omega + 1)/(phi**var_lambda - t**var_lambda) * phi**(var_lambda - 1) +
            nu * omega * (t**(var_lambda * omega) * (phi**var_lambda - t**var_lambda)**(-omega-1) * phi**(var_lambda - 1)) * var_lambda
        )
    result = result + sum
    return result

def step(steps, alpha=0.0001, epsilon = 0.01):
    var_lambda = 0.01
    nu = 5
    phi = 6
    omega = 1
    t_samples=  [2, 1, 2, 3,5]
    params = [var_lambda, nu, omega, phi]
    for i in range(steps):
        gradient = [
            partial_wrt_lambda(var_lambda, nu, phi, omega,t_samples),
            partial_wrt_nu(var_lambda, nu, phi, omega,t_samples),
            partial_wrt_omega(var_lambda, nu, phi, omega,t_samples),
            partial_wrt_phi(var_lambda, nu, phi, omega,t_samples)]
        for j in range(len(params)):
            params[j] += alpha * gradient[j]

        var_lambda = params[0]
        nu = params[1]
        omega = params[2]
        phi = params[3]
        # calculate new gradient and see if it is close to 0
        print(f"{gradient=}")
    
step(100)
        


        


