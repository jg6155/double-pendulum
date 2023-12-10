"""
Author: Daniel E. Tanagho
-------------------------
This LeapFrog function integrates the four 1-st Order simultaneous ODEs
for the double pendulum problem.
----------------------------------------------------------------------
Units:
------
Angles measured in radians
Time measured in seconds
----------------------------------------------------------------------
Arguments:
----------
w1_init : initial angular velocity for upper pendulum bob
w2_init : initial angular velocity for lower pendulum bob
theta1_init : initial position for upper pendulum bob
theta2_init : initial position for lower pendulum bob
ti = initial time
tf = final time
N = number of intervals between ti and tf
----------------------------------------------------------------------
Output:
-------
w1 = array of values for the angular velocity of upper bob
w2 = array of values for the angular velocity of lower bob
theta1 = array of values for the position of the upper bob
theta2 = array of values for the position of the lower bob
tpoints = array of time points of equal size corresponding to 
		  the coordiantes
----------------------------------------------------------------------

"""
import numpy as np 
import matplotlib.pyplot as plt
from DoublePendulum_Funcs import w1_dot, w2_dot, Energy


def LeapFrog_Int(w1_init, w2_init, theta1_init, theta2_init, ti, tf, N):
	## Time Step
	h = (tf - ti)/N
	tpoints = np.arange(ti, tf, h)

	## Arrays to store the values
	w1 = np.zeros(N)
	w2 = np.zeros(N)
	theta1 = np.zeros(N)
	theta2 = np.zeros(N)
	w1_halves = np.zeros(N)
	w2_halves = np.zeros(N)
	theta1_halves = np.zeros(N)
	theta2_halves = np.zeros(N)

	## Starting Points
	w1[0] = w1_init
	w2[0] = w2_init
	theta1[0] = theta1_init
	theta2[0] = theta2_init

	## Initial Half Steps with Euler's Method
	w1_halves[0] = w1_init + 0.5*h*w1_dot(w1_init, w2_init, theta1_init, theta2_init)
	w2_halves[0] = w2_init + 0.5*h*w2_dot(w1_init, w2_init, theta1_init, theta2_init)

	theta1_halves[0] = theta1_init ##Initial velocity is zero
	theta2_halves[0] = theta2_init ##Initial velocity is zero

	## Integration Loop
	for i in range(N-1):
		w1[i+1] = w1[i] + h*w1_dot(w1_halves[i], w2_halves[i], theta1_halves[i], theta2_halves[i])
		w2[i+1] = w2[i] + h*w2_dot(w1_halves[i], w2_halves[i], theta1_halves[i], theta2_halves[i])
		theta1[i+1] = theta1[i] + h*w1_halves[i]
		theta2[i+1] = theta2[i] + h*w2_halves[i]

		w1_halves[i+1] = w1_halves[i] + h*w1_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		w2_halves[i+1] = w2_halves[i] + h*w2_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		theta1_halves[i+1] = theta1_halves[i] + h*w1[i+1]
		theta2_halves[i+1] = theta2_halves[i] + h*w2[i+1]

	return (w1, w2, theta1, theta2, tpoints)



