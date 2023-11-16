import numpy as np 
import matplotlib.pyplot as plt

def W_Dot(W, Theta, t):
	"""
	Calculating each successive value of Angular Acceleration (W_Dot) 
	using the previous value of W and Theta
	"""
	g = 9.8
	l = 0.4
	w1 = W[0]
	w2 = W[1] 
	theta1 = Theta[0]
	theta2 = Theta[1]
	w1_dot = (w1**2*np.sin(2*theta1 - 2*theta2) + 2*w2**2*np.sin(theta1 - theta2) \
			+ (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))) / (np.cos(2*theta1 - 2*theta2) - 3)
	w2_dot = (4*w1**2*np.sin(theta1 - theta2) + w2**2*np.sin(2*theta1 - 2*theta2) \
			+ 2*(g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))) / (3 - np.cos(2*theta1 - 2*theta2))
	return np.array([w1_dot, w2_dot],float)

def Theta_Dot(W, Theta, t):
	"""
	Calculating each successive value of Angular Velocity (Theta_Dot, equivalently W)
	using the previous value of W
	"""
	g = 9.8
	l = 0.4
	w1 = W[0]
	w2 = W[1] 
	theta1 = Theta[0]
	theta2 = Theta[1]
	theta1_dot = w1
	theta2_dot = w2
	return np.array([theta1_dot, theta2_dot], float)

def LeapFrog_Int():
	"""Constants of Integration and Time Interval"""
	g = 9.8
	l = 0.4
	ti = 0.0
	tf = 10.0
	N = 1000
	h = (tf - ti)/N

	"""Initial Conditions"""
	theta1_i = np.pi/2
	theta2_i = np.pi/2
	w1_i = 0.0
	w2_i = 0.0 
	Theta_i = np.array([theta1_i, theta2_i], float)
	W_i = np.array([w1_i, w2_i], float)

	"""Lists to store desired values. The Thetas will be used to animate, and both the Thetas
	and the Ws will be used to test the energy conservation
	"""
	tpoints = np.arange(ti, tf, h)
	theta1_points = []
	theta2_points = []
	w1_points = []
	w2_points = []

	"""Calculating the first HALF step of W & Theta using Euler's method"""
	W_first_half = W_i + 0.5*h*W_Dot(W_i, Theta_i, 0)
	Theta_first_half = Theta_i + 0.5*h*Theta_Dot(W_i, Theta_i, 0) 

	"""Calculating the first FULL step of W & Theta using Euler's method"""
	W_first_full = W_i + h*W_Dot(W_first_half, Theta_first_half, 0.5*h)
	Theta_first_full = Theta_i + 0.5*h*Theta_Dot(W_first_half, Theta_first_half, 0.5*h)

	for t in tpoints:
		W_first_full += h*W_Dot(W_first_half, Theta_first_half, (t+0.5*h)) 
		Theta_first_full += h*Theta_Dot(W_first_half, Theta_first_half, (t+0.5*h))

		W_first_half += h*W_Dot(W_first_full, Theta_first_full, (t+h))
		Theta_first_half += h*Theta_Dot(W_first_full, Theta_first_full, (t+h))

		theta1_points.append(Theta_first_full[0])
		theta2_points.append(Theta_first_full[1])
		w1_points.append(W_first_full[0])
		w2_points.append(W_first_full[1])

	plt.plot(tpoints, theta1_points, 'b-', label = 'Theta 1')
	plt.plot(tpoints, theta2_points, 'r-', label = 'Theta 2')
	plt.plot(tpoints, w1_points, 'k-', label = 'W 1')
	plt.plot(tpoints, w2_points, 'g-', label = 'W 2')
	plt.legend(loc = 'best')
	plt.show()



if __name__ == "__main__":
	LeapFrog_Int()











