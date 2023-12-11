#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:36:12 2023

@author: zhaoyilin
"""

import numpy as np
import matplotlib.pyplot as plt
from DoublePendulum_Funcs import w1_dot, w2_dot, g, m
from Runge_Kutta_method import runge_kutta



#Initial Conditions
time = 10

#The time steps should be less than 0.02s
timestep = 0.05

# theta should be in radians instead of degress
initial_theta_1 = np.pi/4
initial_theta_2 = np.pi/4


# Equation for theta1
def dp3(w1):
    return w1

# 
def dp4(w2):
    return w2



# Create lists to contain stepsizes and errors
theta1_error = []
theta2_error = [0]
theta1_points = [0]
step_points = []

stepsize = timestep

for i in np.arange(1,10,1):
    # For each loops, the stepsize decreases by a factor of 2
    stepsize = stepsize/2
    
    theta1_l,theta2_l, w1_points, w2_points, tpoints = runge_kutta(L, g, initial_theta_1, initial_theta_2, time, stepsize)

    # Get the error by subtract the theta from the previous theta degree
    error1 = np.abs(theta1_l[-1] - theta1_points[i-1])
    
    # Collect stepsize for each loop
    step_points.append(stepsize)
    
    theta1_points.append(theta1_l[-1])
    
    # Collect the error for each loop
    theta1_error.append(error1)
    

    
# Simulation of the error changes
y_points = []
for i in step_points[1:]:
    # We know that the error changes in the order of 4.
    y_points.append((i**4*2000))

# Make the plot
plt.plot(step_points[1:], theta1_error[1:], label = "Convergence for RK method")
plt.scatter(step_points[1:], theta1_error[1:])
plt.plot(step_points[1:], y_points, label = r"$y = 2000x^4$")

plt.title("Test the Convergence of Runge_Kutta method")
plt.xlabel("time step(s)")
plt.ylabel(r"error of $\theta$(radians)")
plt.yscale('log')
plt.xscale('log')
plt.legend(loc = "upper left")
plt.show()
