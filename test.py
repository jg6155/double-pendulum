from animator import PendulumAnimator
from runge_kutta import position
import numpy as np

theta_1, theta_2 = position(1, 9.8, 30.0*np.pi/180, 15.0*np.pi/180, 250, 2000)
animator = PendulumAnimator()
animator.set_data(theta_1, theta_2)
animator.animate()

