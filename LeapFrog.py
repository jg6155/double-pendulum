import numpy as np 
import matplotlib.pyplot as plt
from DoublePendulum_Funcs import w1_dot, w2_dot, Energy


def LeapFrog(w1_init, w2_init, theta1_init, theta2_init, h):
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

	return (w1, w2, theta1, theta2)








