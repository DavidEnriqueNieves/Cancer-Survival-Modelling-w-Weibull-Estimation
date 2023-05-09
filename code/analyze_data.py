# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:17:27 2022

@author: james
"""

import seaborn
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

def analyze_data(data, png_location = 'dataplot.png'):
    ''' 
    Plots the nxn plots of values for each of the data points to visualize the 
    correlation between each variable. A correlogram.
    data = pd.Dataframe
    png_location = str
    
    '''
    grid = seaborn.PairGrid(data)
    grid.map_diag(seaborn.displot)
    grid.map_lower(seaborn.regplot)
    grid.map_upper(corr)
    grid.savefig(png_location)
    
    # plots = seaborn.pairplot(data)
    # plots.savefig(png_location)
    
def corr(x, y, label=None,color=None,**kwargs):
    """
    Allows the upper traingle of matrix plot to be string values instead of a 
    duplication of the lower half.
    """
    data = pd.DataFrame({'x':x, 'y':y})
    r = data.corr()
    ax = plt.gca()
    r = r['x']['y']
    ax.annotate('r = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center', fontsize=25)
    ax.set_axis_off()
    
# data = pd.read_csv(r'C:\Users\james\Documents\projdata.csv')
# analyze_data(data, r'C:\Users\james\Documents\corr_plot.png' )

def f(t, lam, mu, phi, omega):
    # Calculate F(t)
    if phi < t:
        return None
    p = ((pow(t,lam) / ((pow(phi,lam) - pow(t,lam)))))
    pp = -mu * pow(p, omega)
    x = pow(math.e, pp)
    return(1 - x)

def h(t, lam, mu, phi, omega): 
    # Hazard function
    if phi < t:
        return None
    h = ( mu * omega * lam * pow(phi, lam) * pow(t, omega * lam - 1))
    h = h/(pow( pow(phi,lam) - pow(t,lam) , omega + 1))
    h = h * pow(pow(t,lam) / (pow(phi,lam) - pow(t,lam)) , omega)
    return h

t_values = np.linspace(0,10,num=1000, retstep=True)[0]

# Constants
mu=1
lam=1
omega=1
all_lines = {}
for phi in range(0,10):
    # iterate through fi
    f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
    plt.plot(t_values, f_vals, label=str(phi))

plt.xlabel("t")
plt.ylabel("F")
plt.title('F(t)')
plt.show()



phispace = np.linspace(0, 4, 50)
for phi in phispace:
    # print(Fspace)
    hazard = [h(t,lam,mu,phi,omega) if t < phi else 0 for t in t_values]
    plt.plot(t_values, hazard)
# plt.title(f" Plot of H(t) versus t with λ : {l}, ν : {v}, φ : {omega}, ω : {w}")
plt.xlabel("t")
plt.xlim([min(t_values), max(phispace)]) # evaluates to ax.xlim([0, 15]) with example data above
plt.ylabel("H")
plt.title('H(t)')
plt.show()




