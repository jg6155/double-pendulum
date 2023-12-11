import numpy as np
import matplotlib.pyplot as plt
from runge_kutta import runge_kutta

class Poincare:
    def __init__(self, L, g, time, interval, num_nearby_points, epsilon):
        self.L = L
        self.g = g
        self.time = time
        self.interval = interval
        self.num_nearby_points = num_nearby_points
        self.epsilon = epsilon

    def generate_nearby_points(self, fiducial):
        """Generate nearby points around the fiducial point."""
        nearby_points = []
        for _ in range(self.num_nearby_points):
            perturbation = np.random.normal(0, self.epsilon, 4)
            point = np.array(fiducial) + perturbation
            nearby_points.append(point)
        return nearby_points

    def euclidean_distance(self, p1, p2):
        """Calculate Euclidean distance between two points in 4D space."""
        return np.linalg.norm(np.array(p1) - np.array(p2))

    def run_simulation(self, initial_conditions):
        """Run the simulation for given initial conditions."""
        # Ensure initial_conditions is a tuple of length 4 (theta1, theta2, w1, w2)
        if len(initial_conditions) != 4:
            raise ValueError("Initial conditions should be a tuple of four values: theta1, theta2, w1, w2")

        initial_theta_1, initial_theta_2, initial_w1, initial_w2 = initial_conditions
        # Call the runge_kutta function with the correct number of arguments
        return runge_kutta(self.L, self.g, initial_theta_1, initial_theta_2, self.time, self.interval)
    
    def generate_poincare_section(self, initial_conditions):
        """Generate and plot a Poincaré section for given initial conditions."""
        theta1_vals, theta2_vals, w1_vals, w2_vals = self.run_simulation(initial_conditions)

        poincare_theta2 = []
        poincare_w2 = []

        # Map theta values to range [-π, π]
        theta1_vals = np.arctan2(np.sin(theta1_vals), np.cos(theta1_vals))
        theta2_vals = np.arctan2(np.sin(theta2_vals), np.cos(theta2_vals))
        # Define the Poincaré section conditions (θ1 = 0 and θ1_dot > 0)
        for i in range(1, len(theta1_vals)):
            if abs(theta1_vals[i]) < .1 and w1_vals[i] > 0:
                poincare_theta2.append(theta2_vals[i])
                poincare_w2.append(w2_vals[i])

        plt.figure(figsize=(8, 6))
        plt.scatter(poincare_theta2, poincare_w2, s=10)
        plt.xlabel('Theta2')
        plt.ylabel('Angular Velocity w2')
        plt.title(f'Poincaré Section for Initial Conditions: {initial_conditions}')
        plt.show()

    
    def analyze_trajectories(self, fiducial_points):
        for idx, fiducial in enumerate(fiducial_points):
            nearby_points = self.generate_nearby_points(fiducial)
            distances_over_time = []

            # Simulate fiducial trajectory
            theta1_vals_fiducial, theta2_vals_fiducial, w1_vals_fiducial, w2_vals_fiducial = self.run_simulation(fiducial)

            # Simulate nearby trajectories and calculate distances
            for nearby in nearby_points:
                theta1_vals, theta2_vals, w1_vals, w2_vals = self.run_simulation(nearby)
                distances = [self.euclidean_distance((t1, t2, w1, w2), (t1f, t2f, w1f, w2f)) 
                             for t1, t2, w1, w2, t1f, t2f, w1f, w2f in zip(theta1_vals, theta2_vals, w1_vals, w2_vals, theta1_vals_fiducial, theta2_vals_fiducial, w1_vals_fiducial, w2_vals_fiducial)]
                distances_over_time.append(distances)

            # Calculate cumulative differences and select trajectories
            cumulative_differences = [np.sum(distances) for distances in distances_over_time]
            max_diff_idx = np.argmax(cumulative_differences)
            min_diff_idx = np.argmin(cumulative_differences)
            median_diff_idx = np.argsort(cumulative_differences)[len(cumulative_differences) // 2]

            # Create a new figure for each fiducial point
            plt.figure(figsize=(10, 5))
            for selected_idx in [max_diff_idx, median_diff_idx, min_diff_idx]:
                initial_conditions = nearby_points[selected_idx]
                label = f'IC: {", ".join(f"{val:.2f}" for val in initial_conditions)}'
                plt.plot(range(len(distances_over_time[selected_idx])), distances_over_time[selected_idx], label=label)

            # Set titles and labels for the plot
            fiducial_str = ', '.join(f"{val:.4g}" for val in fiducial) if isinstance(fiducial, tuple) else f"{fiducial:.4g}"
            plt.title(f'Distances for Energy Starting at ({fiducial_str})')
            plt.xlabel('Time Steps')
            plt.ylabel('Distance from Fiducial Trajectory')
            plt.legend()

            # Display the plot
            plt.show()

# Example usage
poincare = Poincare(L=0.4, g=9.8, time=100, interval=5000, num_nearby_points=10, epsilon=0.1)
epsilon = .05
stationary_points = [
    (0+epsilon, 0, 0, 0),
    (np.pi +epsilon, 0, 0, 0),
    ( 0, np.pi+epsilon, 0, 0)
]
poincare.analyze_trajectories(stationary_points)

initial_conditions_sets = [
    (0 + epsilon, 0, 0, 0),  # Both pendulums slightly perturbed from hanging down
    (-0.19, -0.06, 0.1, 0.1),
    (-0.19, 0.06, 0.1, 0.1),
    (0.06, -.19, 0.1, 0.1),
    (np.pi + epsilon, 0, 0, 0),  # One pendulum inverted, the other hanging down
    (-3.14, -0.32, 0.1, 0.1),
    (-3.14, -0.19, 0.1, 0.1),
    (-3.14, -0.06, 0.1, 0.1),
    (0, np.pi + epsilon, 0, 0),  # Both pendulums slightly perturbed from being inverted
    (0.32, 3.14, 0.1, 0.1),
    (0.19, 3.14, 0.1, 0.1),
    (0.06, 3.14, 0.1, 0.1)
]
from animator import PendulumAnimator



# Generate Poincaré sections for each set of initial conditions
for ic in initial_conditions_sets:
    #poincare.generate_poincare_section(ic)
    # theta1, theta2, _, _ = poincare.run_simulation(ic)
    # animate = PendulumAnimator()
    # animate.set_data((theta1,theta2))
    # animate.animate()
    pass