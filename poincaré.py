from runge_kutta import runge_kutta
import matplotlib.pyplot as plt

L = 0.4
g = 9.8
m = 1
time = 100
interval = time * 50
initial_theta_1 = 30
initial_theta_2 = 30

theta1_vals, theta2_vals, w1_vals, w2_vals = runge_kutta(L, g, initial_theta_1, initial_theta_2, time, interval)
phase_plane = []
for theta1, theta2, w1, w2 in zip(theta1_vals, theta2_vals, w1_vals, w2_vals):
    if w1 > 0 and theta1 == 0:
        phase_plane.append((theta2, w2))

# Extract theta2 and w2 values for plotting
theta2_points = [point[0] for point in phase_plane]
w2_points = [point[1] for point in phase_plane]

# Plotting
plt.scatter(theta2_points, w2_points)
plt.xlabel('Theta 2')
plt.ylabel('Angular Velocity (w2)')
plt.title('Phase Plane Plot')
plt.show()
