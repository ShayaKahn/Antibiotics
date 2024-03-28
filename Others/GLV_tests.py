import numpy as np
from GLV_model import Glv

n_samples = 100
n_species = 100
sigma = 0.1
mean = 0.0
tf = 1000
max_step = 0.1
delta = 1e-3
eps = 1e-8
method = 'BDF'#'RK45'#
p_vector = np.random.uniform(0.1, 0.5, (1, n_species))
def set_initial_conditions(n_samples, n_species):
    y0 = np.random.uniform(0, 1, (n_samples, n_species))
    return y0

def set_interaction_matrix(n_species, mean, sigma, p_vector):
    # Create a num_species x num_species matrix filled with zeros
    interaction_matrix = np.zeros((n_species, n_species))
    # Generate random numbers
    random_numbers = np.random.normal(mean, sigma, size=(n_species, n_species))*3
    # Create a mask based on the probability p
    mask = np.random.rand(n_species, n_species) < p_vector
    # Apply the mask to the random numbers
    interaction_matrix[mask] = random_numbers[mask]
    return interaction_matrix

def set_growth_rate(n_species):
    return np.random.uniform(0, 1, n_species)

def set_logistic_growth(n_species):
    return np.ones(n_species)

initial_cond = set_initial_conditions(n_samples, n_species)
interaction_matrix = set_interaction_matrix(n_species, mean, sigma, p_vector)
r = set_growth_rate(n_species)
s = set_logistic_growth(n_species)

glv_object = Glv(n_samples, n_species, delta, r, s, interaction_matrix, initial_cond,
                 tf, max_step, normalize=True, method=method)
final_abundance = glv_object.solve()

import matplotlib.pyplot as plt

plt.imshow(final_abundance, aspect='auto')

plt.show()

a=1