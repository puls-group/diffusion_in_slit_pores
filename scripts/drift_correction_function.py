#!/usr/bin/env python3
import csv
import sys
import math
import matplotlib.pyplot as plt
import argparse

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# from graph_tools import *
# setPlotHalfA4Portrait()


def drift_lifetime_coeff_gamma(gamma, accuracy=100):
    # This function calculates K_B(gamma) with the denoted numerical accuracy
    # Higher accuracy includes more terms of the series
    drift_ratio = gamma

    sum = 0.0

    for n in range(1, accuracy+1):
        if n % 2 == 0:
            alternate_coeff = 1
        else:
            alternate_coeff = -1

        upper = (1-alternate_coeff*np.cos(drift_ratio/2.))
        offset = (drift_ratio/math.pi)**2/4.
        lower_first = n**2 + offset

        sum += upper/(lower_first**2)

    sum *= 24./(math.pi**4)*drift_ratio**2 / (np.cosh(drift_ratio)-1)

    return sum


def drift_lifetime_coeff(L, D, mu, accuracy=100):
    # Helper function to calculate K_B(gamma), where gamma is calculated
    # for the user from the values of L, D and \mu
    gamma = L*mu/D
    return drift_lifetime_coeff_gamma(gamma, accuracy=accuracy)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Script to calculate the value of $K_B(\gamma)$.')
    parser.add_argument('gamma', required=True,
                        help='argument of the function $K_B$ to be calculated')

    args = parser.parse_args()

    gamma = args.gamma

    # Calculate the function $K_B(gamma)$:
    drift_coeff = drift_lifetime_coeff_gamma(gamma)

    print("K_B(", gamma, ") \t=\t", drift_coeff)

    """ If you want the function to be calculated for a range of gamma instead, use this:
    
    gamma_range = np.linspace(0, 8, num=400)

    drift_coeff = drift_lifetime_coeff_gamma(gamma_range)
    """
