from GLV_model import Glv
from Historical_contingency import *
import numpy as np
from Historical_contingency import HC

n_samples = 100
n_species = 50
n_survived = 30
sigma = 0.1
mean = 0.0
tf = 1000
max_step = 0.5
delta = 1e-3
eps = 1e-8
#p_vector = np.random.uniform(0.1, 0.5, (1, n_species))
p_vector = np.ones((1, n_species))
HC_object = HC(n_samples, n_species, n_survived, mean,
               sigma, p_vector, delta, tf, max_step, eps)
