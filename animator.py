import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PendulumAnimator:
    def __init__(self, l=1.0):
        self.l = l
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], 'o-', lw=2)

    def set_data(self, results):
        self.theta1, self.theta2 = results

    def init_animation(self):
        self.line.set_data([], [])
        return self.line,

    def animate_frame(self, i):
        x = [0, self.l * np.sin(self.theta1[i]), self.l * (np.sin(self.theta1[i]) + np.sin(self.theta2[i]))]
        y = [0, -self.l * np.cos(self.theta1[i]), -self.l * (np.cos(self.theta1[i]) + np.cos(self.theta2[i]))]
        self.line.set_data(x, y)
        return self.line,

    def animate(self):
        self.ax.set_xlim(-2.5, 2.5)
        self.ax.set_ylim(-2.5, 2.5)
        ani = animation.FuncAnimation(self.fig, self.animate_frame, frames=len(self.theta1), init_func=self.init_animation, blit=True)
        plt.show()

# Generating test data
t = np.linspace(0, 10, 200)
theta1_vals = np.pi/4 * np.sin(t)
theta2_vals = np.pi/4 * np.sin(t + np.pi/8)  # slightly out of phase with the first pendulum

# Using the animator
animator = PendulumAnimator()
animator.set_data((theta1_vals, theta2_vals))
animator.animate()