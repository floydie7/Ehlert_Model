'''
ehlert.py
@author: Benjamin Floyd

This script is designed to take in the appropriate data and preform a Markov Chain Monte Carlo analysis in order to
determine the posterior distributions of a modified model from Ehlert et al. 2015.
'''

#%matplotlib inline
from __future__ import print_function
from scipy.integrate import quadrature
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scipy.optimize as op
import emcee
import corner
import tools # Bandon Decker's tools package

np.random.seed(100)

# Define constants
M_sun =

# Define calculated quantities
D_A = tools.angDiaDist(z)
r500 = tools.m500Tor500(M500)

# Define modified Ehlert Model
N_obs = N * D_A**2 * r500 * Phi * ( r / r500 )**beta + C
N = N_0 * ( M500 / 1e15 * M_sun)**zeta
beta = beta_0 + beta_m * ( M500 / 1e15 * M_sun)


# Luminosity Function from Assef+2011
PhidL = lambda L : Phi_star * 1.0 / ( (L / L_star)**gamma_1 + (L / L_star)**gamma_2 )
epsilon_L = 1e-5 # Epsilon to create narrow interval in L for our LF
Phi = quadrature(PhidL, L - epsilon_L, L + epsilon_L) # Compute the density from the LF at a given luminosity

# Define our likelihood function, our priors, and the probability
def lnlike(param, N_obs, M500, r500, D_A, ): # What are the independent variables for this?
    gamma_1, gamma_2, L_star, Phi_star, zeta, beta_0, beta_m, C = param
    model = N * D_A**2 * r500 * Phi * (r / r500)**beta +C
    inv_sigma2 = 1.0/() # What are the errors for this?
    return -0.5 * (np.sum( N_obs - model))**2 * inv_sigma2)) # What are the data points here? N_obs?

def lnprior(param):
    gamma_1, gamma_2, L_star, Phi_star, zeta, beta_0, beta_m, C = param
    if gamma_1 in np.random.normal(3.35,0.05) and gamma_2 in np.random.normal(0.37,0.27) \
            and L_star in np.random.normal() and Phi_star in np.random.normal() \
            and C in np.random.normal(330,33) and -5.0 < zeta < 5.0 and -5.0 < beta_0 < 5.0 and -5.0 < beta_m < 5.0
        return 0.0
    return -np.inf

def lnprob(param, N_obs, ):
    lp = lnprior(param)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(param, N_obs, )

