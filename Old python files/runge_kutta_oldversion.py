#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:22:32 2023

@author: zhaoyilin
"""
import numpy as np
import matplotlib.pyplot as plt

L = 1
g = 9.8


def position(L, g, theta_1, theta_2, t, N):
    a = 0.0
    b = t 
    # h is the time interval
    h = (b-a)/N 


    tpoints = np.arange(a,b,h)
    
    # Create two lists to store the angle degrees
    theta1_points = []
    theta2_points = []
    
    # initial velocity
    w1 = 0
    w2 = 0
    for t in tpoints:
        theta1_points.append(theta_1)
        theta2_points.append(theta_2)
        k11 = h*dp1(w1, w2 ,t, theta_1, theta_2)
        k12 = h*dp2(w1, w2 ,t, theta_1, theta_2)
    
        k21 = h*dp1(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1, theta_2)
        k22 = h*dp2(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1, theta_2)
    
        k31 = h*dp1(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1, theta_2)
        k32 = h*dp2(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1, theta_2)
    
        k41 = h*dp1(w1 + k31, w2 + k32, t + h, theta_1, theta_2)
        k42 = h*dp2(w1 + k31, w2 + k32, t + h, theta_1, theta_2)
    
        w1 += (k11 + 2*k21 + 2*k31 + k41)/6
        w2 += (k12 + 2*k22 + 2*k32 + k42)/6
    
        theta_1 += w1*h
        theta_2 += w2*h

    
    plt.figure(figsize=(18, 6))    
    plt.plot(tpoints, theta1_points, label = "theta_1")
    plt.plot(tpoints, theta2_points, label = "theta_2")
    plt.xlabel("t")
    plt.ylabel("degree of angle theta")
    plt.legend()
  
    txt1="The initial theta_1 is: " + str(theta_1)
    txt2 ="  The initial theta_2 is: " + str(theta_2)
    plt.figtext(0.5, 0.01, txt1 + txt2, wrap=True, horizontalalignment='center', fontsize=12)
    plt.show()
    
    return theta1_points, theta2_points


# define the two first-order differentiation equation
# Equation for w1
def dp1(w1, w2, t, theta_1, theta_2):
    w1_term = w1**2*(np.sin(2*theta_1 - 2*theta_2))
    w2_term = 2*w2**2*np.sin(theta_1 - theta_2)
    ###
    numerator = w1_term*(-1) - w2_term - (g/L)*(np.sin(theta_1 - 2*theta_2) + 3*np.sin(theta_1))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom

# Equation for w2
def dp2(w1, w2, t, theta_1, theta_2):
    w1_term = 4*w1**2*(np.sin(theta_1 - theta_2))
    w2_term = w2**2*(np.sin(2*theta_1 - 2*theta_2))
    numerator = w1_term + w2_term + (2*g/L)*(np.sin(2*theta_1 - theta_2) - np.sin(theta_2))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom
