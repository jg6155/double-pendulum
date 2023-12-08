import numpy as np 
import matplotlib.pyplot as plt

m = 1.0
g = 9.8
l = 0.4
ti = 0.0
tf = 10.0
N = 10000
h = (tf - ti)/N
tpoints = np.arange(ti, tf, h)

def w_dot(theta):
	return -(g/l)*np.sin(theta)

def LeapFrogInt():
	w_init = 0.0
	theta_init = np.pi/2

	w = np.zeros(N)
	theta = np.zeros(N)

	w_halves = np.zeros(N)
	theta_halves = np.zeros(N)

	w[0] = w_init
	theta[0] = theta_init

	w_halves[0] = w_init + 0.5*h*w_dot(theta_init)
	theta_halves[0] = theta_init

	for i in range(N-1):
		w[i+1] = w[i] + h*w_dot(theta_halves[i])
		theta[i+1] = theta[i] + h*w[i]

		w_halves[i+1] = w_halves[i] + h*w_dot(theta[i+1])
		theta_halves[i+1] = theta_halves[i] + h*w[i+1]

	return w, theta

w, theta = LeapFrogInt()

def Energy(w, theta):
	return (0.5*m*l**2*w**2 + m*g*l*(1-np.cos(theta)))

E = Energy(w, theta)

plt.plot(tpoints, E)
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.show()




