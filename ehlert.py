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

np.random.seed(100)

# Define modified Ehlert Model
N_obs = N * D_A**2 * r500 * Phi * ( r / r500 )**beta + C
N = N_0 * ( M500 / 1e15 * M_sun)**zeta
beta = beta_0 + beta_m * ( M500 / 1e15 * M_sun)


# Luminosity Function from Assef+2011
PhidL = lambda L : Phi_star * 1.0 / ( (L / L_star)**gamma_1 (L / L_star)**gamma_2 )
epsilon_L = 1e-5 # Epsilon to create narrow interval in L for our LF
Phi = quadrature(PhidL, L - epsilon_L, L + epsilon_L) # Compute the density from the LF at a given luminosity





