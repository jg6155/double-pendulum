from animator import PendulumAnimator

import numpy as np
import matplotlib.pyplot as plt

# Initial Conditions for comparing the answeres of two integrations
# Attention! The codes may not run unless the time interval is higher than 50 per second.
L = 0.4
g = 9.8
m = 1
time = 100
interval = time*50
initial_theta_1 = 30
initial_theta_2 = 30

def position(L, g, theta_1, theta_2, t, N):
    a = 0.0
    b = t 
    # h is the time interval
    h = (b-a)/N 


    tpoints = np.arange(a,b,h)
    
    # Create two lists to store the angle degrees
    theta1_points = []
    theta2_points = []
    
    # initial velocity
    w1 = 0
    w2 = 0
    for t in tpoints:
        theta1_points.append(theta_1)
        theta2_points.append(theta_2)
        k11 = h*dp1(w1, w2 ,t, theta_1, theta_2)
        k12 = h*dp2(w1, w2 ,t, theta_1, theta_2)
    
        k21 = h*dp1(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1, theta_2)
        k22 = h*dp2(w1 + 0.5*k11, w2 + 0.5*k12, t + 0.5*h, theta_1, theta_2)
    
        k31 = h*dp1(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1, theta_2)
        k32 = h*dp2(w1 + 0.5*k21, w2 + 0.5*k22, t + 0.5*h, theta_1, theta_2)
    
        k41 = h*dp1(w1 + k31, w2 + k32, t + h, theta_1, theta_2)
        k42 = h*dp2(w1 + k31, w2 + k32, t + h, theta_1, theta_2)
    
        w1 += (k11 + 2*k21 + 2*k31 + k41)/6
        w2 += (k12 + 2*k22 + 2*k32 + k42)/6
    
        theta_1 += w1*h
        theta_2 += w2*h

    
    plt.figure()    
    plt.plot(tpoints, theta1_points, label = "theta_1")
    plt.plot(tpoints, theta2_points, label = "theta_2")
    plt.title("Position of double pendulum versusu time")
    plt.xlabel("t")
    plt.ylabel("degree of angle theta")
    plt.legend()
   
    plt.show()
    
    return theta1_points, theta2_points


# define the two first-order differentiation equation
# Equation for w1
def dp1(w1, w2, t, theta_1, theta_2):
    w1_term = w1**2*(np.sin(2*theta_1 - 2*theta_2))
    w2_term = 2*w2**2*np.sin(theta_1 - theta_2)
    ###
    numerator = w1_term*(-1) - w2_term - (g/L)*(np.sin(theta_1 - 2*theta_2) + 3*np.sin(theta_1))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom

# Equation for w2
def dp2(w1, w2, t, theta_1, theta_2):
    w1_term = 4*w1**2*(np.sin(theta_1 - theta_2))
    w2_term = w2**2*(np.sin(2*theta_1 - 2*theta_2))
    numerator = w1_term + w2_term + (2*g/L)*(np.sin(2*theta_1 - theta_2) - np.sin(theta_2))
    denom = 3 - np.cos(2*theta_1 - 2*theta_2)
    return numerator/denom

theta1, theta2 = position(L, g, initial_theta_1*np.pi/180, initial_theta_2*np.pi/180, time, interval)

animator = PendulumAnimator()
animator.set_data((theta1, theta2))
animator.animate()

