import numpy as np 
import matplotlib.pyplot as plt
from animator import PendulumAnimator
from DoublePendulum_Funcs import w1_dot, w2_dot, Energy

ti = 0.0
tf = 1000.0
N = 10000
h = (tf - ti)/N
tpoints = np.arange(ti, tf, h)
"""Initial Conditions"""
theta1_init = np.pi/2
theta2_init = -np.pi/2
w1_init = 0.0
w2_init = 0.0 

def LeapFrog_Int(w1_init, w2_init, theta1_init, theta2_init, h):
	"""Arrays to store my values"""
	w1 = np.zeros(N)
	w2 = np.zeros(N)
	theta1 = np.zeros(N)
	theta2 = np.zeros(N)
	w1_halves = np.zeros(N)
	w2_halves = np.zeros(N)
	theta1_halves = np.zeros(N)
	theta2_halves = np.zeros(N)

	"""Starting Points"""
	w1[0] = w1_init
	w2[0] = w2_init
	theta1[0] = theta1_init
	theta2[0] = theta2_init

	"""Initial Half Steps with Euler's Method"""
	w1_halves[0] = w1_init + 0.5*h*w1_dot(w1_init, w2_init, theta1_init, theta2_init)
	w2_halves[0] = w2_init + 0.5*h*w2_dot(w1_init, w2_init, theta1_init, theta2_init)

	theta1_halves[0] = theta1_init ##Initial velocity is zero
	theta2_halves[0] = theta2_init ##Initial velocity is zero

	"""Integration Loop"""
	for i in range(N-1):
		w1[i+1] = w1[i] + h*w1_dot(w1_halves[i], w2_halves[i], theta1_halves[i], theta2_halves[i])
		w2[i+1] = w2[i] + h*w2_dot(w1_halves[i], w2_halves[i], theta1_halves[i], theta2_halves[i])
		theta1[i+1] = theta1[i] + h*w1_halves[i]
		theta2[i+1] = theta2[i] + h*w2_halves[i]

		w1_halves[i+1] = w1_halves[i] + h*w1_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		w2_halves[i+1] = w2_halves[i] + h*w2_dot(w1[i+1], w2[i+1], theta1[i+1], theta2[i+1])
		theta1_halves[i+1] = theta1_halves[i] + h*w1[i+1]
		theta2_halves[i+1] = theta2_halves[i] + h*w2[i+1]

	"""Plotting the Test Data"""
	# plt.plot(tpoints, theta1, 'b-', label = r'$\theta_1$')
	# plt.plot(tpoints, theta2, 'r-', label = r'$\theta_2$')
	# plt.plot(tpoints, w1, 'k-', label = r'$\omega_1$')
	# plt.plot(tpoints, w2, 'g-', label = r'$\omega_2$')
	# plt.xlabel('Time (s)')
	# plt.ylabel('Time-Dependent Variable')
	# plt.legend(loc = 'best')
	# plt.show()

	return (w1, w2, theta1, theta2)

w1, w2, theta1, theta2 = LeapFrog_Int(w1_init, w2_init, theta1_init, theta2_init, h)

# E = Energy(w1, w2, theta1, theta2)

# print(theta1)

# def Plotter(w1, w2, theta1, theta2, E):
# 	plt.figure(1)
# 	plt.plot(tpoints, theta1, 'b-', label = r'$\theta_1$')
# 	plt.plot(tpoints, theta2, 'r-', label = r'$\theta_2$')
# 	plt.plot(tpoints, w1, 'k-', label = r'$\omega_1$')
# 	plt.plot(tpoints, w2, 'g-', label = r'$\omega_2$')
# 	plt.xlabel('Time (s)')
# 	plt.ylabel('Time-Dependent Variable')
# 	plt.legend(loc = 'best')
# 	plt.show()

# 	plt.figure(2)
# 	plt.plot(tpoints, E, 'k-')
# 	plt.xlabel('Time (s)')
# 	plt.ylabel('Energy (Joules)')
# 	plt.show()

# Plotter(w1, w2, theta1, theta2, E)

"""Saving Data Into CSV files"""
# np.savetxt('w1_double.csv', w1, delimiter=',')
# np.savetxt('w2_double.csv', w2, delimiter=',')
# np.savetxt('theta1_double.csv', theta1, delimiter=',')
# np.savetxt('theta2_double.csv', theta2, delimiter=',')
# np.savetxt('E_double.csv', E, delimiter=',')

# animator = PendulumAnimator()
# animator.set_data((theta1, theta2))
# animator.animate()






