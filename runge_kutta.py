#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 19:33:14 2023

@author: zhaoyilin
"""

# Input: length of string, g, initial theta_1, initial theta_2, time, intervals

import numpy as np

def runge_kutta(L, g, initial_theta_1, initial_theta_2, time, interval):
    def dp1(w1, w2, theta_1, theta_2):
        # Differential equation for w1
        w1_term = w1**2 * np.sin(2*theta_1 - 2*theta_2)
        w2_term = 2*w2**2 * np.sin(theta_1 - theta_2)
        numerator = -w1_term - w2_term - (g/L) * (np.sin(theta_1 - 2*theta_2) + 3*np.sin(theta_1))
        denom = 3 - np.cos(2*theta_1 - 2*theta_2)
        return numerator/denom

    def dp2(w1, w2, theta_1, theta_2):
        # Differential equation for w2
        w1_term = 4*w1**2 * np.sin(theta_1 - theta_2)
        w2_term = w2**2 * np.sin(2*theta_1 - 2*theta_2)
        numerator = w1_term + w2_term + (2*g/L) * (np.sin(2*theta_1 - theta_2) - np.sin(theta_2))
        denom = 3 - np.cos(2*theta_1 - 2*theta_2)
        return numerator/denom

    a = 0.0
    b = time
    h = (b - a) / interval

    tpoints = np.arange(a, b, h)

    theta1 = initial_theta_1
    theta2 = initial_theta_2
    w1 = 0.0
    w2 = 0.0

    theta1_points, theta2_points, w1_points, w2_points = [], [], [], []

    for t in tpoints:
        theta1_points.append(theta1)
        theta2_points.append(theta2)
        w1_points.append(w1)
        w2_points.append(w2)

        k11 = h * dp1(w1, w2, theta1, theta2)
        k12 = h * dp2(w1, w2, theta1, theta2)

        k21 = h * dp1(w1 + 0.5*k11, w2 + 0.5*k12, theta1, theta2)
        k22 = h * dp2(w1 + 0.5*k11, w2 + 0.5*k12, theta1, theta2)

        k31 = h * dp1(w1 + 0.5*k21, w2 + 0.5*k22, theta1, theta2)
        k32 = h * dp2(w1 + 0.5*k21, w2 + 0.5*k22, theta1, theta2)

        k41 = h * dp1(w1 + k31, w2 + k32, theta1, theta2)
        k42 = h * dp2(w1 + k31, w2 + k32, theta1, theta2)

        w1 += (k11 + 2*k21 + 2*k31 + k41) / 6
        w2 += (k12 + 2*k22 + 2*k32 + k42) / 6

        theta1 += w1 * h
        theta2 += w2 * h

    return theta1_points, theta2_points, w1_points, w2_points
