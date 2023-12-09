#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 22:02:13 2023

@author: zhaoyilin
"""

import numpy as np
import matplotlib.pyplot as plt

#Initial Conditions
L = 0.4
g = 9.8
time = 100

#The time steps should be less than 0.02s
timestep = 0.01

#Input the theta degrees in radians
initial_theta_1 = np.pi/2
initial_theta_2 = np.pi/2

def runge_kutta(L, g, theta_1, theta_2, t, N):
    a = 0.0
    b = t 
    # h is the time interval
    h = N


    tpoints = np.arange(a,b,h)
    
    # Create two lists to store the angle degrees
    theta1_points = []
    theta2_points = []
    w1_points = []
    w2_points = []
    
    # initial velocity
    w1 = 0
    w2 = 0
    for t in tpoints:
        
        k11 = h*dp1(w1, w2 ,t, theta_1, theta_2)
        k12 = h*dp2(w1, w2 ,t, theta_1, theta_2)
        k13 = h*dp3(w1)
        k14 = h*dp4(w2)
        
        
        k21 = h*dp1(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1 + 0.5*k13, theta_2 + 0.5*k14)
        k22 = h*dp2(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1 + 0.5*k13, theta_2 + 0.5*k14 )
        k23 = h*dp3(w1 + 0.5*k11)
        k24 = h*dp4(w2 + 0.5*k12)
    
        k31 = h*dp1(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1 + 0.5*k23, theta_2 + 0.5*k24)
        k32 = h*dp2(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1 + 0.5*k23, theta_2 + 0.5*k24)
        k33 = h*dp3(w1 + 0.5*k21)
        k34 = h*dp4(w2 + 0.5*k22)
    
        k41 = h*dp1(w1 + k31, w2 + k32, t + h, theta_1 + k33, theta_2 + k34)
        k42 = h*dp2(w1 + k31, w2 + k32, t + h, theta_1 + k33, theta_2 + k34)
        k43 = h*dp3(w1 + k31)
        k44 = h*dp4(w2 + k32)
        
        #Update w1, w2, theta1 and theta2
        w1 += (k11 + 2*k21 + 2*k31 + k41)/6
        w2 += (k12 + 2*k22 + 2*k32 + k42)/6
        theta_1 += (k13 + 2*k23 + 2*k33 + k43)/6
        theta_2 += (k14 + 2*k24 + 2*k34 + k44)/6
        
        #Record w1, w2, theta1 and theta2
        w1_points.append(w1)
        w2_points.append(w2)
        theta1_points.append(theta_1)
        theta2_points.append(theta_2)
        
    
    return theta1_points,theta2_points, w1_points, w2_points, tpoints

# define the two first-order differentiation equation
# Equation for w1
def dp1(w1, w2, t, theta_1, theta_2):
    g = 9.8
    L = 0.4
    w1_term = w1**2*(np.sin(2*theta_1 - 2*theta_2))
    w2_term = 2*w2**2*np.sin(theta_1 - theta_2)
    ###
    numerator = w1_term*(-1) - w2_term - (g/L)*(np.sin(theta_1 - 2*theta_2) + 3*np.sin(theta_1))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom

# Equation for w2
def dp2(w1, w2, t, theta_1, theta_2):
    g = 9.8
    L = 0.4
    w1_term = 4*w1**2*(np.sin(theta_1 - theta_2))
    w2_term = w2**2*(np.sin(2*theta_1 - 2*theta_2))
    numerator = w1_term + w2_term + (2*g/L)*(np.sin(2*theta_1 - theta_2) - np.sin(theta_2))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom

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
    m = 1
    
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
