import numpy as n
import pyspec as p
import time as ti
import random as r
import numpy.linalg as lin
from scipy.optimize import minimize, curve_fit
from scipy.special import erf
import matplotlib.pyplot as plt
from sympy.physics.wigner import wigner_6j


filename = '4497fit.txt'
def spec_fit(guess):
	vec[12] = guess[0]
	vec[10] = guess[1]
	vec[14] = guess[2]
	vec[15] = filename
	execfile('OPV4.py')
	return chi_2_red

guess = [1.0,300,0.4]
min = minimize(spec_fit,guess)