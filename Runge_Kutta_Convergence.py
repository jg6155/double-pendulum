#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:36:12 2023

@author: zhaoyilin
"""

import numpy as np
import matplotlib.pyplot as plt
from DoublePendulum_Funcs import w1_dot, w2_dot

L = 0.4
g = 9.8

time = 10

#The time steps should be less than 0.02s
timestep = 0.05

# theta should be in radians instead of degress
initial_theta_1 = np.pi/4
initial_theta_2 = np.pi/4

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
        
        k11 = h*w1_dot(w1, w2 ,theta_1, theta_2)
        k12 = h*w2_dot(w1, w2 ,theta_1, theta_2)
        k13 = h*dp3(w1)
        k14 = h*dp4(w2)
        
        
        k21 = h*w1_dot(w1 + 0.5*k11, w2 + 0.5*k12, theta_1 + 0.5*k13, theta_2 + 0.5*k14)
        k22 = h*w2_dot(w1 + 0.5*k11, w2 + 0.5*k12, theta_1 + 0.5*k13, theta_2 + 0.5*k14 )
        k23 = h*dp3(w1 + 0.5*k11)
        k24 = h*dp4(w2 + 0.5*k12)
    
        k31 = h*w1_dot(w1 + 0.5*k21, w2 + 0.5*k22, theta_1 + 0.5*k23, theta_2 + 0.5*k24)
        k32 = h*w2_dot(w1 + 0.5*k21, w2 + 0.5*k22, theta_1 + 0.5*k23, theta_2 + 0.5*k24)
        k33 = h*dp3(w1 + 0.5*k21)
        k34 = h*dp4(w2 + 0.5*k22)
    
        k41 = h*w1_dot(w1 + k31, w2 + k32, theta_1 + k33, theta_2 + k34)
        k42 = h*w2_dot(w1 + k31, w2 + k32, theta_1 + k33, theta_2 + k34)
        k43 = h*dp3(w1 + k31)
        k44 = h*dp4(w2 + k32)
        
        
        w1 += (k11 + 2*k21 + 2*k31 + k41)/6
        w2 += (k12 + 2*k22 + 2*k32 + k42)/6
        theta_1 += (k13 + 2*k23 + 2*k33 + k43)/6
        theta_2 += (k14 + 2*k24 + 2*k34 + k44)/6
        
        
        w1_points.append(w1)
        w2_points.append(w2)
        theta1_points.append(theta_1)
        theta2_points.append(theta_2)
        
    
    return theta1_points,theta2_points, w1_points, w2_points, tpoints


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
