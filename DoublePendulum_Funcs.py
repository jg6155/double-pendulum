"""
Author: Daniel E. Tanagho
-------------------------
This file contains all the double pendulum constants:
-----------------------------------------------------
m (mass of pendulum bobs) (kg)
g (gravitational acceleration) (m/s^2)
l (length of pendulum rods) (m)

As well as the functions for the two 1-st order ODEs and the Total Energy:
--------------------------------------------------------------------------
w1_dot: Time derivative of the upper bob's angular velocity (rad/s^2)
w2_dot: Time derivative of the lower bob's angular velocity (rad/s^2)
Energy: Total Energy of the system (Joules)

Arguments:
----------
w1 = array of values for the angular velocity of upper bob (rad/s)
w2 = array of values for the angular velocity of lower bob (rad/s)
theta1 = array of values for the position of the upper bob (rad)
theta2 = array of values for the position of the lower bob (rad)
--------------------------------------------------------------------------
"""
import numpy as np
import matplotlib.pyplot as plt

m = 1.0
g = 9.8
L = 0.4
l = 0.4

def w1_dot(w1, w2, theta1, theta2):
	return ((w1**2*np.sin(2*theta1 - 2*theta2) + 2*w2**2*np.sin(theta1 - theta2) + (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))) / (np.cos(2*theta1 - 2*theta2) - 3))

def w2_dot(w1, w2, theta1, theta2):
	return ((4*w1**2*np.sin(theta1 - theta2) + w2**2*np.sin(2*theta1 - 2*theta2) + 2*(g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))) / (3 - np.cos(2*theta1 - 2*theta2)))

def Energy(w1, w2, theta1, theta2):
	return (m*l**2*(w1**2 + 0.5*w2**2 + w1*w2*np.cos(theta1-theta2)) - m*g*l*(2*np.cos(theta1) + np.cos(theta2)))
