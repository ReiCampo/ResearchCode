#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:37:01 2019

@author: RachelCampo
"""

import math
import fileinput
import numpy as np
from numpy import exp as e
import scipy as sp
import matplotlib.pyplot as plt
import astropy as ap
from astropy import units as u
from astropy import constants as c
from astropy.cosmology import Planck15 as p15

from scipy.integrate import quad
from scipy.misc import derivative
from scipy.integrate import odeint

n = 6.2 * 10 ** -10
m_e = 511 * 10 ** 3 # in eV/c^2
Q = 13.57 # in eV
k = 8.617 * 10 ** -5 #in ev/K
T = np.linspace(0, 1e4, 1000)


#Part 1a

S = 3.84 * n * ((k*T/m_e))**(3/2)*(e(Q/(k*T)))

def saha(T):
    X = (-1+np.sqrt(1+4*S))/(2*S)
    return X

def test(T):
    return np.sin(T)

'''
fig = plt.figure(figsize=(7,7))
axis = fig.add_subplot(1,1,1)

axis.plot(T, saha(T), 'g', linewidth = 1.5)
axis.set_xlabel('Temperature (K)', fontsize = 10)
axis.set_ylabel('Saha Equation', fontsize = 10)
axis.set_title('Saha Equation as a Function of Temperature', fontsize = 20)

plt.show()
'''

#Part 1b

H_0 = (p15.H(0)).to('1/yr').value


#def Friedman(a):
#    Or = 9.03 * 10**-5
#    Om = 0.306
#    Ol = 0.692
#    FrEq = 1/(H_0*(Or*a**-2 + Om*a**-1 + Ol*a**2)**(1/2))
#    return FrEq
#
#a = np.linspace(1e-3, 1e-4, 100)
#t_space = np.array([])

###############################################################################
Or = 9.03e-5

def Friedman(a, Or):
    Om = 0.306
    Ol = 0.692
    FrEq = 1/((H_0) * np.sqrt((Or*a**-2 + Om*a**-1 + Ol*a**2)))
    return FrEq

all_a_list = np.logspace(-5,10,1000)
all_t_list = np.array([quad(Friedman, 0, all_a_list[x], args=(Or))[0] for x in range(len(all_a_list))])
###############################################################################

def a_f(t):
    return np.interp(t, all_t_list, all_a_list)

t_list = np.linspace(1.e4, 3.e5, 1000)
a_list = a_f(t_list)

def X_t(t):
    a=a_f(t)
    T0=2.7
    T=T0/a
    X=saha(T)
    return X

    
fig = plt.figure(figsize=(7,7))
axis = fig.add_subplot(1,1,1)
axis.plot(t_list, X_t(t_list), 'b')
plt.show()

'''
fig=plt.figure(figsize=(5,4))
ax=fig.add_subplot(111)
ax.plot(all_t_list,all_a_list)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('t')
ax.set_ylabel('a')
plt.show()
plt.close()
'''

fig=plt.figure(figsize=(5,4))
ax=fig.add_subplot(111)
ax.plot(all_t_list,a_f(all_t_list))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('t')
ax.set_ylabel('a')
plt.show()
