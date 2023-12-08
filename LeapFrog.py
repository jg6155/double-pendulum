import numpy as np 
import matplotlib.pyplot as plt
# from animator import PendulumAnimator

m = 1.0
g = 9.8
L = 0.4
ti = 0.0
tf = 5
N = 100
h = (tf - ti)/N
tpoints = np.arange(ti, tf, h)
# w1_dot = (w1**2*np.sin(2*theta1 - 2*theta2) + 2*w2**2*np.sin(theta1 - theta2) + (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))) / (np.cos(2*theta1 - 2*theta2) - 3)
# return w1_dot

# w2_dot = (4*w1**2*np.sin(theta1 - theta2) + w2**2*np.sin(2*theta1 - 2*theta2) + 2*(g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))) / (3 - np.cos(2*theta1 - 2*theta2))
# return w2_dot


def w1_dot(w1, w2, theta1, theta2):
	w1_term = w1**2*(np.sin(2*theta1 - 2*theta2))
	w2_term = 2*w2**2*np.sin(theta1 - theta2)
	numerator = w1_term*(-1) - w2_term - (g/L)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))
	denom = 3 - np.cos(2*theta1 - 2*theta2)
	return numerator/denom


def w2_dot(w1, w2, theta1, theta2):
	w1_term = 4*w1**2*(np.sin(theta1 - theta2))
	w2_term = w2**2*(np.sin(2*theta1 - 2*theta2))
	numerator = w1_term + w2_term + (2*g/L)*(np.sin(2*theta1 - theta2) - np.sin(theta2))
	denom = 3 - np.cos(2*theta1 - 2*theta2)
	return numerator/denom

def LeapFrog_Int():
	"""Initial Conditions"""
	theta1_init = np.pi/2
	theta2_init = np.pi/2
	w1_init = 0.0
	w2_init = 0.0 

	"""Arrays to store my values"""
	w1 = np.zeros(N)
	w2 = np.zeros(N)
	theta1 = np.zeros(N)
	theta2 = np.zeros(N)
	w1_halves = np.zeros(N)
	w2_halves = np.zeros(N)
	theta1_halves = np.zeros(N)
	theta2_halves = np.zeros(N)

	"""Initial Half Steps with Euler's Method"""
	w1_FirstHalf = w1_init + 0.5*h*w1_dot(w1_init, w2_init, theta1_init, theta2_init)
	w2_FirstHalf = w2_init + 0.5*h*w2_dot(w1_init, w2_init, theta1_init, theta2_init)

	theta1_FirstHalf = theta1_init ##Initial velocity is zero
	theta2_FirstHalf = theta2_init ##Initial velocity is zero

	"""Initial Full Steps with Euler's Method"""	
	w1_FirstFull = w1_init + h*w1_dot(w1_FirstHalf, w2_FirstHalf, theta1_FirstHalf, theta2_FirstHalf)
	w2_FirstFull = w2_init + h*w2_dot(w1_FirstHalf, w2_FirstHalf, theta1_FirstHalf, theta2_FirstHalf)

	theta1_FirstFull = theta1_init + h*w1_FirstHalf
	theta2_FirstFull = theta2_init + h*w2_FirstHalf

	w1[0] = w1_init
	w1[1] = w1_FirstFull
	w2[0] = w2_init
	w2[1] = w2_FirstFull
	theta1[0] = theta1_init
	theta1[1] = theta1_FirstFull
	theta2[0] = theta2_init
	theta2[1] = theta2_FirstFull

	w1_halves[0] = w1_FirstHalf
	w2_halves[0] = w2_FirstHalf
	theta1_halves[0] = theta1_FirstHalf
	theta2_halves[0] = theta2_FirstHalf

	"""Integration Loop"""
	for i in range(N-2):
		w1_halves[i+1] = w1_halves[i] + h*w1_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		w2_halves[i+1] = w2_halves[i] + h*w2_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		theta1_halves[i+1] = theta1_halves[i] + h*w1[i+1]
		theta2_halves[i+1] = theta2_halves[i] + h*w2[i+1]

		w1[i+2] = w1[i+1] + h*w1_dot(w1_halves[i+1], w2_halves[i+1], theta1_halves[i+1], theta2_halves[i+1])
		w2[i+2] = w2[i+1] + h*w2_dot(w1_halves[i+1], w2_halves[i+1], theta1_halves[i+1], theta2_halves[i+1])
		theta1[i+2] = theta1[i+1] + h*w1_halves[i+1]
		theta2[i+2] = theta2[i+1] + h*w2_halves[i+1]

	"""Plotting the Test Data"""
	plt.plot(tpoints, theta1, 'b-', label = r'$\theta_1$')
	plt.plot(tpoints, theta2, 'r-', label = r'$\theta_2$')
	plt.plot(tpoints, w1, 'k-', label = r'$\omega_1$')
	plt.plot(tpoints, w2, 'g-', label = r'$\omega_2$')
	plt.xlabel('Time (s)')
	plt.ylabel('Time-Dependent Variable')
	plt.legend(loc = 'best')
	plt.show()

	return (w1, w2, theta1, theta2)

w1, w2, theta1, theta2 = LeapFrog_Int()

def EnergyPlotter(w1, w2, theta1, theta2):
	"""Total Energy"""
	Energy = m*L**2*(w1**2 + 0.5*w2**2 + w1*w2*np.cos(theta1-theta2)) - m*g*L*(2*np.cos(theta1) + np.cos(theta2))

	"""Total Energy Plot as a function of Time"""
	plt.plot(tpoints, Energy, 'k-')
	plt.xlabel('Time (s)')
	plt.ylabel('Energy (Joules)')
	plt.show()

EnergyPlotter(w1, w2, theta1, theta2)
# animator = PendulumAnimator()
# animator.set_data((theta1, theta2))
# animator.animate()






