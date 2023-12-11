#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 21:15:28 2023

@author: zhaoyilin
"""

import numpy as np
import matplotlib.pyplot as plt
from animator import PendulumAnimator
from DoublePendulum_Funcs import w1_dot, w2_dot, L, g, m
from Runge_Kutta_method import runge_kutta

#Initial Conditions
time = 100

#The time steps should be less than 0.02s
timestep = 0.02

#Input the theta degrees in radians
initial_theta_1 = np.pi/2
initial_theta_2 = np.pi/2


# Equation for theta1
def dp3(w1):
    return w1

# 
def dp4(w2):
    return w2

theta1, theta2, w1_points, w2_points, tpoints = runge_kutta(L, g, initial_theta_1, initial_theta_2, time, timestep)

plt.figure(figsize = (10,5))
plt.plot(tpoints, theta1, label = r"$\theta_1$")
plt.plot(tpoints, theta2, label = r"$\theta_2$")
plt.title("Test of the RK integrator under the condition of part (b)")
plt.xlabel("time(s)")
plt.ylabel("Angle degree(Radian)")
plt.legend(loc = 'upper right')

animator = PendulumAnimator()
animator.set_data((theta1, theta2))
animator.animate()
