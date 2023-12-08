import numpy as np 
import matplotlib.pyplot as plt
from animator import PendulumAnimator

m = 1.0
g = 9.8
l = 0.4
ti = 0.0
tf = 10.0
N = 10000
h = (tf - ti)/N
tpoints = np.arange(ti, tf, h)

def w1_dot(w1, w2, theta1, theta2):
	return ((w1**2*np.sin(2*theta1 - 2*theta2) + 2*w2**2*np.sin(theta1 - theta2) + (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))) / (np.cos(2*theta1 - 2*theta2) - 3))


def w2_dot(w1, w2, theta1, theta2):
	return ((4*w1**2*np.sin(theta1 - theta2) + w2**2*np.sin(2*theta1 - 2*theta2) + 2*(g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))) / (3 - np.cos(2*theta1 - 2*theta2)))


def LeapFrog_Int():
	"""Initial Conditions"""
	theta1_init = np.pi/10
	theta2_init = np.pi/10
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

w1, w2, theta1, theta2 = LeapFrog_Int()

# np.savetxt('w1_double.csv', w1, delimiter=',')
# np.savetxt('w2_double.csv', w2, delimiter=',')
# np.savetxt('theta1_double.csv', theta1, delimiter=',')
# np.savetxt('theta2_double.csv', theta2, delimiter=',')

# print(w1, w2, theta1, theta2, sep='\n')

def Energy(w1, w2, theta1, theta2):
	return (m*l**2*(w1**2 + 0.5*w2**2 + w1*w2*np.cos(theta1-theta2)) - m*g*l*(2*np.cos(theta1) + np.cos(theta2)))

E = Energy(w1, w1, theta1, theta2)

# np.savetxt('E_double.csv', E, delimiter=',')

plt.plot(tpoints, E, 'k-')
plt.xlabel('Time (s)')
plt.ylabel('Energy (Joules)')
plt.show()

animator = PendulumAnimator()
animator.set_data((theta1, theta2))
animator.animate()






