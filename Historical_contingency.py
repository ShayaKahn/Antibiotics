import numpy as np
import random
from GLV_model import Glv

class HC:
    """
    ....
    """
    def __init__(self, num_samples, pool_size, num_survived, mean,
                 sigma, p_vector, delta, final_time, max_step, epsilon):
        """
        :param num_samples: The number of samples.
        :param pool_size: The total size of the population.
        :param num_survived: The number of survived species.
        """
        self.num_samples = num_samples
        self.pool_size = pool_size
        self.num_survived = num_survived
        self.mean = mean
        self.sigma = sigma
        self.p_vector = p_vector
        self.y0, self.survived_matrix = self.set_initial_conditions()
        self.A = self.set_interaction_matrix()
        self.r = self.set_growth_rate()
        self.s = self.set_logistic_growth()
        self.steady_state = np.linalg.solve(self.A+np.diag(self.s), -self.r)
        self.delta = delta
        self.final_time = final_time
        self.max_step = max_step
        self.perturbed_state = self.apply_GLV(self.y0, norm=False)
        self.epsilon = epsilon
        self.new_perturbed_state = self.insert_total_pool()
        self.post_pertubed_state = self.apply_GLV(self.new_perturbed_state, norm=True)

    def set_initial_conditions(self):
        y0 = np.zeros((self.num_samples, self.pool_size))
        survived_matrix = np.array([random.sample(range(0, self.pool_size-1),
                                                  self.num_survived) for _ in range(self.num_samples)])
        for y, s in zip(y0, survived_matrix):
            y[s] = np.random.rand(1, self.num_survived)
        return y0, survived_matrix

    def set_interaction_matrix(self):
        # Create a num_species x num_species matrix filled with zeros
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))
        # Generate random numbers
        random_numbers = np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size))
        # Create a mask based on the probability p
        mask = np.random.rand(self.pool_size, self.pool_size) < self.p_vector
        # Apply the mask to the random numbers
        interaction_matrix[mask] = random_numbers[mask]
        return interaction_matrix

    def set_growth_rate(self):
        return np.random.uniform(0, 1, self.pool_size)

    def set_logistic_growth(self):
        return np.ones(self.pool_size)

    def apply_GLV(self, init_cond, norm):
        glv_object = Glv(self.num_samples, self.pool_size, self.delta, self.r, self.s, self.A, init_cond, self.final_time,
                         self.max_step, normalize=norm)
        final_abundances = glv_object.solve()
        return final_abundances

    def insert_total_pool(self):
        new_perturbed_state = self.perturbed_state.copy()
        for p, s in zip(new_perturbed_state, self.survived_matrix):
            mask = np.ones(p.shape[0], dtype=bool)
            mask[s] = False
            p[mask] = self.epsilon
        return new_perturbed_state