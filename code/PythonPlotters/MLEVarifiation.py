import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_csv(r'..\\projdata.csv')
tlist = []


def ln(x):
    n = float(x)
    return math.log(float(n))

def calculate_optimal_mu(phi, omega, lam, t):
    n = len(t)
    sum = 0
    for t_i in t:
        sum+=((t_i**lam) / (phi**lam - t_i**lam))**omega
    mu = n / sum
    return mu

def lnMLE(lam, mu, phi, omega, tlist):
    sum = 0
    for Ti in tlist:
        print(f"{sum=}")
        if phi < Ti:
            return None
        first = ln(mu) + ln(omega) + ln(lam) + lam * ln(phi) + (lam * omega - 1) * ln(Ti) - (omega + 1) * (
            ln(((phi ** lam) - (Ti ** lam))))
        second = first - (mu * (Ti ** lam) / ((phi ** lam) - (Ti ** lam))) ** omega
        sum +=second
    return sum

if __name__ == "__main__":
    phi = 1022.01

    mu = 4.9

    # lambda_space = np.arange(1, 75,0.01)
    # lambda_space = list(lambda_space)
    # omega_space = np.arange(0, 5, 0.01)
    # omega_space = list(omega_space)

    n = 10
    lam_bounds=(0,70)
    lambda_space = np.arange(lam_bounds[0], lam_bounds[1],(lam_bounds[1] - lam_bounds[0])/n)
    lambda_space = list(lambda_space)

    om_bounds = (0, 7)
    omega_space = np.arange(om_bounds[0],om_bounds[1], (om_bounds[1] - om_bounds[0]) / n)
    omega_space = list(omega_space)
    
    lnMLE_Space = []
    print(f"{lambda_space=}")
    print(f"{omega_space=}")

    for i, value in enumerate(lambda_space):
        print(f"{value=}")
        omega_vars = []
        for j, value in enumerate(omega_space):
            print(f"{value=}")
            mle = lnMLE(lambda_space[i], mu, phi, omega_space[j], tlist)
            print(f"{mle=}")
            omega_vars.append(mle)
        lnMLE_Space.append(omega_vars)


    lam_bounds=(0,70)
    lambda_space = np.arange(lam_bounds[0], lam_bounds[1],(lam_bounds[1] - lam_bounds[0])/n)
    lambda_space = list(lambda_space)

    om_bounds = (0, 7)
    omega_space = np.arange(om_bounds[0],om_bounds[1], (om_bounds[1] - om_bounds[0]) / n)
    omega_space = list(omega_space)

    final_triples = [[omega_space[i], lambda_space[i], lnMLE_Space[i]] for i, val in enumerate(lambda_space)]
    print(f"{final_triples=}")

    # plt.xlabel("MLE")
    # plt.ylabel("lambda")
    ax = plt.axes(projection='3d')
    lnMLE_Space = np.asarray(lnMLE_Space)
    omega_space = np.asarray(omega_space)
    lambda_space = np.asarray(lambda_space)
    
    print(len(lambda_space))
    print(len(omega_space))
    print(len(lnMLE_Space))
    plotx,ploty, = np.meshgrid(np.linspace(np.min(lambda_space),np.max(lambda_space),10),\
                           np.linspace(np.min(omega_space),np.max(omega_space),10))
    # plotz = interp.griddata((lambda_space, omega_space),lnMLE_Space,(plotx,ploty),method='linear')


    #  ax.plot_surface(omega_space, lambda_space, lnMLE_Space)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(plotx, ploty, lnMLE_Space, alpha=0.5, rstride=1, cstride=1)


    plt.title('MLE wrt lambda')
    plt.show()


    # run_phi()
    # run_omega()
    # run_lam()
    # run_mu()

