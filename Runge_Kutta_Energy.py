#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 22:02:13 2023

@author: zhaoyilin
"""

import numpy as np
import matplotlib.pyplot as plt
from DoublePendulum_Funcs import w1_dot, w2_dot, g, m
from Runge_Kutta_method import runge_kutta



#Initial Conditions
time = 100

#The time steps should be less than 0.02s
timestep = 0.01

#Input the theta degrees in radians
initial_theta_1 = np.pi/2
initial_theta_2 = np.pi/2



# Equation for theta1
def dp3(w1):
    return w1

# Equation for theta2
def dp4(w2):
    return w2

# Method used to calculate total energy
def energy(w_1, w_2, theta_1, theta_2):
    total_en_list = []
    T = 0
    V = 0
    
    # Number of Timesteps
    interval = int(time/timestep)
    
    for i in np.arange(0, interval):
        T = m*(L**2)*(w_1[i]**2 + 0.5*w_2[i]**2 + w_1[i]*w_2[i]*np.cos(theta_1[i] - theta_2[i]))
        V = -m*g*L*(2*np.cos(theta_1[i]) + np.cos(theta_2[i]))
        total_en = T + V
        total_en_list.append(total_en)
        
    time_list = np.linspace(0, time, interval)
    
    # Make the plot
    
    fig = plt.figure(figsize = (7,4))
    #plt.figure(figsize = (16,6))
    plt.plot(time_list, total_en_list, label = r"Total Energy")
    plt.title("Total Energy under RK method for 1000s")
    plt.xlabel("time(s)")
    plt.ylabel("Total Energy(J)")
    plt.legend(loc = 'upper right')
    txt=r"Under Condition of part (b): initial $\theta_1 = \theta_2= \pi/2$, initial $\omega_1 = \theta_2 = 0$, L = 0.4m, $g = 9.8m~s^{-2}$, m = 1kg "
    fig.text(.5, .0001, txt, ha='center')
    plt.tight_layout()
    plt.show()

    return total_en_list

theta1, theta2, w1, w2,tpoints = runge_kutta(L, g, initial_theta_1, initial_theta_2, time, timestep)
total_energy = energy(w1, w2, theta1, theta2)
