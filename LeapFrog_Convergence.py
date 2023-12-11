"""
Author: Daniel E. Tanagho
-------------------------
This piece of code tests the convergence of the LeapFrog integrator
as a function of time step.

Parameters:
-----------
w1_init, w2_init, theta1_init, theta2_init : 
initial angular velocities (rad/s)
and initial positions (rad) of the pendulum bobs

M : the power of 2 that divides the original time step of 0.01s
"""

import numpy as np 
import matplotlib.pyplot as plt
from LeapFrog import LeapFrog_Int
from DoublePendulum_Funcs import Energy
from animator import PendulumAnimator

w1_init = 0.0
w2_init = 0.0
theta1_init = -np.pi/100
theta2_init = np.pi/100
ti = 0.0
tf = 10.0
N = 1000
M = 15
theta2 = np.zeros([M+1, (2**M*N)])
tpoints = np.zeros([M+1, (2**M*N)])
theta2_diff = np.zeros(M)


for i in range(M+1):
	w1, w2, theta1, theta2[i, 0:((2**i)*N)], tpoints[i, 0:((2**i)*N)] = LeapFrog_Int(w1_init, w2_init, theta1_init, theta2_init, ti, tf, ((2**i)*N))

for j in range(M):
	theta2_diff[j] = np.abs(theta2[j+1, ((2**(j+1)*N) - 1)] - theta2[j, (2**j*N - 1)])


np.savetxt('LeapFrog_theta2_Conv.csv', theta2_diff, delimiter=',')

delta_t = np.zeros(15)

for i in range(15):
	delta_t[i] = 1.e-2 / (2**i)

log10_theta2_diff = np.log10(theta2_diff)
log10_delta_t = np.log10(delta_t)

plt.plot(log10_delta_t, log10_theta2_diff, 'bo', label=r'$\Delta \theta_2$')
plt.xlim(-1.5, -6.5)
plt.xlabel(r'$log_{10}(\Delta t)$')
plt.ylabel(r'$log_{10}(\Delta \theta_2)$')
plt.legend(loc='best')
plt.savefig('LeapFrog_Conv_NoChaos.png')

