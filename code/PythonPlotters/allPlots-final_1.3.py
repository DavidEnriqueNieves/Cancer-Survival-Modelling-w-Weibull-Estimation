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

def run_corr_plot():
    """
    Essentially the 'main' for generating the correlelogram
    """
    data = pd.read_csv(r'C:\Users\james\Documents\projdata.csv')
    data.drop(columns=['Unnamed: 0', 'inst', 'status'], inplace=True)
    
    # filter the data to a specific group (male, female, etc)
    filt = data[data.sex == 2]
    filt.drop(columns=['sex'], inplace=True)
    
    # run analysis
    analyze_data(filt, r'C:\Users\james\Documents\corr_plot.png' )



## Functions for problem functions
def f(t, lam, mu, phi, omega):
    # Calculate F(t)
    if phi < t:
        return None
    f_t = ((pow(t,lam) / ((pow(phi,lam) - pow(t,lam)))))
    f_t = -mu * pow(f_t, omega)
    f_t = pow(math.e, f_t)
    return(1 - f_t)

def h(t, lam, mu, phi, omega):
    # Hazard function
    if phi < t:
        return None
    try:
        h_t = ( mu * omega * lam * pow(phi, lam) * pow(t, lam * omega - 1))
    except:
        print(mu, omega, lam, phi, t)
    h_t = h_t/(pow( pow(phi,lam) - pow(t,lam) , omega + 1))
    return h_t

def pdf(t, lam, mu, phi, omega):
    # Calculate f(t)
    # Since this is a PDF, each individual value should be close to 0. F(t) matters more
    prob = h(t, lam, mu, phi, omega) * ( 1 - f(t, lam, mu, phi, omega))
    return prob

# Defined values for t and interative space for variables
t_values = np.linspace(0,1020,num=1000, retstep=True)[0]
iter_space = np.linspace(0, 1, 10)

# #------------------------------------------------------------------------------------------
# # w.r.t phi
def run_phi():
    # # Constants
    mu=1
    lam=1
    omega=1
    
    # Graph F(t)
    legend = []
    for phi in range(0,10):
        f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, f_vals, label=str(phi))
        legend.append(str(phi))
    
    plt.xlabel("t")
    plt.ylabel("F")
    plt.legend(legend)
    plt.title('F(t) wrt phi')
    plt.show()
    
    
    # Graph h(t)
    legend=[]
    for phi in iter_space:
        hazard = [h(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, hazard)
        legend.append(str(phi))
    
    plt.xlabel("t")
    plt.ylabel("H")
    plt.legend(legend)
    plt.title('H(t) w.r.t. phi')
    plt.show()

# #------------------------------------------------------------------------------------------
# # w.r.t mu
def run_mu():
    phi=50
    lam=1
    omega=1
    
    legend=[]
    for mu in range(0,50, 10):
        f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, f_vals, label=str(mu))
        legend.append(str(mu))
    
    plt.xlabel("t")
    plt.ylabel("F")
    plt.legend(legend)
    plt.title('F(t) w.r.t. mu')
    plt.show()
    
    legend=[]
    for mu in iter_space:
        hazard = [h(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, hazard)
        legend.append(str(mu))
    
    plt.xlabel("t")
    plt.ylabel("H")
    plt.legend(legend)
    plt.title('H(t) w.r.t. mu')
    plt.show()


# #------------------------------------------------------------------------------------------
# # w.r.t lam
def run_lam():
    phi=5
    mu=1
    omega=1
    
    legend=[]
    for lam in range(40,50):
        f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, f_vals, label=str(phi))
        legend.append(str(lam))
    
    plt.xlabel("t")
    plt.ylabel("F")
    plt.legend(legend)
    plt.title('F(t) w.r.t. lam')
    plt.show()
    
    legend=[]
    for lam in iter_space:
        hazard = [h(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, hazard)
        legend.append(str(lam))
    
    plt.xlabel("t")
    plt.ylabel("H")
    plt.legend(legend)
    plt.title('H(t) w.r.t. lam')
    plt.show()


# #------------------------------------------------------------------------------------------
# # w.r.t omega
def run_omega():
    phi=1
    lam=1
    mu=1
    
    legend=[]
    for omega in range(0,10):
        f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, f_vals, label=str(phi))
        legend.append(str(omega))
    
    plt.xlabel("t")
    plt.ylabel("F")
    plt.legend(legend)
    plt.title('F(t) w.r.t. omega')
    plt.show()
    
    legend=[]
    for omega in iter_space:
        hazard = [h(t,lam,mu,phi,omega) for t in t_values]
        plt.plot(t_values, hazard)
        legend.append(str(omega))
    
    plt.xlabel("t")
    plt.ylabel("H")
    plt.legend(legend)
    plt.title('H(t) w.r.t. omega')
    plt.show()


def plot_f():
    data = pd.read_csv(r'C:\Users\james\Documents\projdata.csv')
    
    x = []
    y = []
    itersize = 50
    for num in range(10,1030,itersize):
        count = len( data[(data.time <= num) & (data.time > num - itersize)])
        x.append(num)
        y.append(count)
        
    # Plot PDF
    plt.plot(x, np.array(y)/sum(y))
    plt.title('PDF f(t) based on data')
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.show()
    
    x = []
    y = []
    itersize = 10
    for num in range(10,1030,itersize):
        count = len( data[(data.time <= num) & (data.time > num - itersize)])
        x.append(num)
        y.append(count)
        
    
    # Plot CDF
    plt.plot(x, np.cumsum(y)/sum(y))
    plt.xlim(0,max(x))
    plt.ylim(0,1)
    plt.title('CDF F(t) based on data')
    plt.xlabel("t")
    plt.ylabel("F(t)")
    plt.show()
    
    
def plot_true_vs_cal_F(lam,mu,phi,omega):
    data = pd.read_csv(r'C:\Users\james\Documents\projdata.csv')
    
    f_vals = [f(t,lam,mu,phi,omega) for t in t_values]
    
    x = []
    y = []
    itersize = 10
    for num in range(10,1030,itersize):
        count = len( data[(data.time <= num) & (data.time > num - itersize)])
        x.append(num)
        y.append(count)
    
    
    # Plot CDF
    plt.plot(x, np.cumsum(y)/sum(y), label="Data")
    plt.plot(t_values, f_vals, label="Calculated")
    plt.legend(['Data', 'Calculated'])
    plt.xlim(0,max(x))
    plt.ylim(0,1)
    plt.title('CDF F(t) based on data')
    plt.xlabel("t")
    plt.ylabel("F(t)")
    plt.show()


def plot_true_vs_cal_f(lam,mu,phi,omega):
    data = pd.read_csv(r'C:\Users\james\Documents\projdata.csv')
    
    f_vals = [pdf(t,lam,mu,phi,omega) for t in t_values]
    
    x = []
    y = []
    f_vals = []
    itersize = 10
    for num in range(10,1030,itersize):
        count = len( data[(data.time <= num) & (data.time > num - itersize)])
        x.append(num)
        y.append(count)
        f_vals.append(f(num,lam,mu,phi,omega) - f(num-itersize,lam,mu,phi,omega))
        
    # Plot PDF
    plt.plot(x, np.array(y)/sum(y), label="Data")
    plt.plot(x, f_vals, label="Calculated")
    plt.legend(['Data', 'Calculated'])
    plt.xlim(0,max(x))
    plt.ylim(0,0.3)
    plt.title(f'PDF f(t) based on data, bin={itersize}')
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.show()
    

if __name__ == "__main__":
    corr_plots = False
    param_plots= False
    plot_dist = False
    plot_true=True
    
    if param_plots:
        print("Generating plots for a delta in each parameter.")
        run_mu()
        run_phi()
        run_lam()
        run_omega()
    if corr_plots:
        print('Generating correlation plots for data analysis.')
        run_corr_plot()
    if plot_dist:
        plot_f()
    if plot_true:
        phi=1022.005
        omega=0.0522538
        mu=4.799849
        lam=27.02045
        
        # plot_true_vs_cal_F(lam,mu,phi,omega)
        
        plot_true_vs_cal_f(lam,mu,phi,omega)