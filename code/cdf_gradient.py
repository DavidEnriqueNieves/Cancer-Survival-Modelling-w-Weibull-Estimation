
import math
var_lambda = 0
nu = 0
phi = 0
omega = 0
t_samples=  [0, 1, 2, 3,5]
a=0

def calculate_a(var_lambda, nu, phi, omega, t_samples):
    a = 0
    for t in t_samples:
        a += (1 /( ((phi**var_lambda)/(t**var_lambda)) - 1))**omega
    a = -nu * a
a = calculate_a(var_lambda, nu, phi, omega, t_samples)

def calculate_prefix(var_lambda=var_lambda, nu=nu, phi=phi, omega=omega, t_samples=t_samples, a=a):
    n = len(t_samples)
    result = 1/(1 + (((-1)**n)/n) * math.e**a) 
    result  *= (((-1)**n)/n) * math.e**a

def partial_wrt_lambda(var_lambda, nu, phi, omega,t_samples, a=a):
    prefix = calculate_prefix()
    result = prefix 

    sum = 0
    for t in t_samples:
        sum+= (
            (
            -1 * (omega * (((phi ** var_lambda )/(t**var_lambda)) - 1)**(omega-1)) * 
            ( t**var_lambda * math.log(phi) * phi **var_lambda - math.log(t) * (t**var_lambda)  * phi**var_lambda )/(t**var_lambda)**2
            )
        
        /

        (((phi ** var_lambda )/(t**var_lambda)) - 1)**(2 * omega)

        )
    sum = -nu * sum

    result = result * sum

    return result


def partial_wrt_nu(var_lambda, nu, phi, omega,t_samples, a=a):
    prefix = calculate_prefix()
    result = prefix * -(a / nu)
    return result

def partial_wrt_phi(var_lambda, nu, phi, omega,t_samples, a=a):
    prefix = calculate_prefix()
    result = prefix 

    sum = 0
    for t in t_samples:
        sum+= (
        (1 / ((phi ** var_lambda ) / (t**var_lambda)) - 1)**(omega - 1) 
        * 
        -1 * (var_lambda * phi**(var_lambda - 1) /t**var_lambda)/(((phi ** var_lambda ) / (t**var_lambda)) - 1)**2
        )
    sum = -nu * omega   * sum

    result = result * sum

    return result

def partial_wrt_omega(var_lambda, nu, phi, omega,t_samples, a=a):
    prefix = calculate_prefix()
    result = prefix * -(a / nu)

    sum = 0
    for t in t_samples:
        sum+= (math.log(1 / ((phi ** var_lambda ) / (t**var_lambda)) - 1) 
        * 
        (1 / ((phi ** var_lambda ) / (t**var_lambda)) - 1)**omega)
    sum = -nu   * sum

    result = result * sum

    return result