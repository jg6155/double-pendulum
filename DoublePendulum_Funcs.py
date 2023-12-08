import numpy as np
import matplotlib.pyplot as plt

m = 1.0
g = 9.8
l = 0.4

def w1_dot(w1, w2, theta1, theta2):
	return ((w1**2*np.sin(2*theta1 - 2*theta2) + 2*w2**2*np.sin(theta1 - theta2) + (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))) / (np.cos(2*theta1 - 2*theta2) - 3))

def w2_dot(w1, w2, theta1, theta2):
	return ((4*w1**2*np.sin(theta1 - theta2) + w2**2*np.sin(2*theta1 - 2*theta2) + 2*(g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))) / (3 - np.cos(2*theta1 - 2*theta2)))

def Energy(w1, w2, theta1, theta2):
	return (m*l**2*(w1**2 + 0.5*w2**2 + w1*w2*np.cos(theta1-theta2)) - m*g*l*(2*np.cos(theta1) + np.cos(theta2)))
