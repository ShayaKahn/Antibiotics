import numpy as np
import random
from Classes.GLV_model import Glv
from Classes.overlap import Overlap
from scipy.spatial.distance import braycurtis
from scipy.stats import powerlaw
import networkx as nx


class HC2:
    """
    This class is responsible for the simulation of the historical contingency model.
    """
    def __init__(self, num_samples, pool_size, num_survived, mean, sigma, c, delta, final_time, max_step, epsilon,
                 threshold, max_growth, min_growth, alpha=None, method='RK45', measure='Jaccard', multiprocess=True):
        """
        :param num_samples: The number of samples.
        :param pool_size: The total size of the population.
        :param num_survived: The number of survived species.
        :param mean: For the interaction matrix generation, the mean of the normal distribution.
        :param sigma: For the interaction matrix generation, The standard deviation of the normal distribution.
        :param c: The Connectance.
        :param delta: The stop condition for the steady state.
        :param final_time: The final time of the integration.
        :param max_step: The maximal allowed step size.
        :param epsilon: The value to insert in the non-survived species.
        :param threshold: The threshold to remove low abundances.
        :param method: The method to solve the GLV model.
        :param measure: The measure to calculate the similarity between the samples.
        :param max_growth: The maximum growth rate.
        :param min_growth: The minimum growth rate.
        :param alpha: The power-law exponent for the interaction matrix strength.
        """
        if not (isinstance(num_samples, int) and isinstance(pool_size, int) and
                isinstance(num_survived, int)):
            raise ValueError("num_samples, pool_size, and num_survived must be integers.")
        if num_samples <= 0 or pool_size <= 0 or num_survived <= 0:
            raise ValueError("num_samples, pool_size, and num_survived must be greater than 0.")
        if num_survived > pool_size:
            raise ValueError("num_survived must be smaller then pool_size.")
        self.num_samples = num_samples
        self.pool_size = pool_size
        self.num_survived = num_survived
        if not (isinstance(mean, float) or isinstance(mean, int)):
            raise ValueError("mean must be integer or float.")
        if not (isinstance(sigma, float) or isinstance(sigma, int)):
            raise ValueError("sigma must be integer or float.")
        self.mean = mean
        self.sigma = sigma
        if method not in ['RK45', 'BDF', 'RK23', 'Radau', 'LSODA', 'DOP853']:
            raise ValueError("method must be one of the following methods: RK45, BDF, RK23, Radau, LSODA, and DOP853")
        self.method = method
        if not (isinstance(c, float)):
            raise ValueError("c must be float.")
        self.c = c
        if not (isinstance(threshold, float) or isinstance(threshold, int)):
            raise ValueError("threshold must be a number.")
        if threshold <= 0:
            raise ValueError("threshold must be greater than 0.")
        self.threshold = threshold
        # Check if delta is a number between 0 and 1
        if not (0 < delta < 1):
            raise ValueError("delta must be a number between 0 and 1.")
        self.delta = delta
        # Check if final_time is a number greater than zero
        if not (isinstance(final_time, (int, float)) and final_time > 0):
            raise ValueError("final_time must be a number greater than zero.")
        # Check if max_step is a number greater than zero and smaller than final_time
        if not (isinstance(max_step, (int, float)) and 0 < max_step < final_time):
            raise ValueError("max_step must be a number greater than zero and smaller than final_time.")
        self.final_time = final_time
        self.max_step = max_step
        if not (isinstance(epsilon, float) and 0 < epsilon < 1):
            raise ValueError("epsilon must be a number between 0 and 1.")
        self.epsilon = epsilon
        if measure not in ['Jaccard', 'braycurtis']:
            raise ValueError("measure must be one of the following measures: Jaccard and braycurtis.")
        self.measure = measure
        if not (isinstance(max_growth, (float, int)) and 0 < max_growth <= 1):
            raise ValueError("max_growth must be a number between 0 and 1.")
        if not (isinstance(min_growth, (float, int)) and 0 <= min_growth < 1):
            raise ValueError("min_growth must be a number between 0 and 1.")
        if min_growth > max_growth:
            raise ValueError("min_growth must be smaller than max_growth.")
        self.max_growth = max_growth
        self.min_growth = min_growth
        if not ((isinstance(alpha, (float, int)) and alpha > 0) or alpha is None):
            raise ValueError("alpha must be a number greater than 0.")
        self.alpha = alpha
        if not (isinstance(multiprocess, bool)):
            raise ValueError("multiprocess must be of type bool.")
        self.multiprocess = multiprocess
        self.r = self._set_growth_rate()
        self.s = self._set_logistic_growth()
        self.y0, self.survived_matrix = self._set_initial_conditions()
        self.N = self._set_symmetric_interaction_matrix()#self._set_interaction_matrix()#
        self.H = self._set_int_strength_matrix()
        self.A = np.dot(self.N, self.H)
        self.perturbed_state = self.y0#self.perturbed_state, event_not_satisfied_ind = self._apply_GLV(self.y0, norm=True)
        self.new_perturbed_state = self._normalize_cohort(self._insert_total_pool())
        self.post_perturbed_state, event_not_satisfied_ind_post = self._apply_GLV(self.new_perturbed_state, norm=True)
        self.filtered_post_perturbed_state = self._remove_low_abundances()
        unique_event_not_satisfied_ind = event_not_satisfied_ind_post#unique_event_not_satisfied_ind = sorted(list(set(event_not_satisfied_ind + event_not_satisfied_ind_post)))
        self.perturbed_state = np.delete(self.perturbed_state, unique_event_not_satisfied_ind, axis=0)
        self.filtered_post_perturbed_state = np.delete(self.filtered_post_perturbed_state,
                                                       unique_event_not_satisfied_ind, axis=0)
        self.survived_matrix = np.delete(self.survived_matrix, unique_event_not_satisfied_ind, axis=0)
        self.shuffled_state = self._generate_shuffled_state()

        self.jaccard_perturbed, self.vals_real, self.vals_shuffled = self._calculate_Jaccard_similarity()

    def _set_initial_conditions(self):
        y0 = np.zeros((self.num_samples, self.pool_size))
        survived_matrix = np.vstack([np.sort(random.sample(range(0, self.pool_size),
                                     self.num_survived)) for _ in range(self.num_samples)])
        for y, s in zip(y0, survived_matrix):
            y[s] = np.random.rand(1, self.num_survived)
        return y0, survived_matrix

    def _set_interaction_matrix(self):
        # Create a num_species x num_species matrix filled with zeros
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))
        random_numbers = np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size))
        # Create a mask based on the probability p
        mask = np.random.rand(self.pool_size, self.pool_size) < self.c
        # Apply the mask to the random numbers
        interaction_matrix[mask] = random_numbers[mask]
        #fp = np.linalg.solve(interaction_matrix + np.diag(self.s), -self.r).reshape(self.pool_size, 1)
        #jacobian = np.multiply(fp, interaction_matrix + np.diag(self.s))
        #eigenvalues, _ = np.linalg.eig(jacobian)
        #try:
        #    assert np.all(eigenvalues.real < 0)
        #except AssertionError:
        #    print("Not all real parts of the eigenvalues are negative.")
        return -np.abs(interaction_matrix)

    def _set_int_strength_matrix(self):
        if self.alpha is None:
            return np.diag(np.ones(self.pool_size))
        else:
            diag = powerlaw.rvs(self.alpha, size=self.pool_size)
            mean = np.mean(diag)
            normalized_diag = diag / mean
            return np.diag(normalized_diag)

    def _set_symmetric_interaction_matrix(self):
        G = nx.binomial_graph(self.pool_size, self.c)
        mask = nx.to_numpy_array(G).astype(bool)
        random_numbers = np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size))
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))
        interaction_matrix[mask] = random_numbers[mask]
        return -np.abs(interaction_matrix)

    def _set_growth_rate(self):
        return np.random.uniform(self.min_growth, self.max_growth, self.pool_size)

    def _set_logistic_growth(self):
        return np.ones(self.pool_size)

    def _apply_GLV(self, init_cond, norm):
        glv_object = Glv(self.num_samples, self.pool_size, self.delta, self.r, self.s, self.A, init_cond,
                         self.final_time, self.max_step, normalize=norm, method=self.method,
                         multiprocess=self.multiprocess)
        final_abundances = glv_object.solve()
        glv_object.close_client()
        return final_abundances

    def _insert_total_pool(self):
        new_perturbed_state = self.perturbed_state.copy()
        for p, s in zip(new_perturbed_state, self.survived_matrix):
            mask = np.ones(p.shape[0], dtype=bool)
            mask[s] = False
            p[mask] = self.epsilon
        return new_perturbed_state

    def _remove_low_abundances(self):
        post_perturbed_state_copy = self.post_perturbed_state.copy()
        zero_ind = np.where(post_perturbed_state_copy < self.threshold)
        post_perturbed_state_copy[zero_ind] = 0.0
        return self._normalize_cohort(post_perturbed_state_copy)

    def _generate_shuffled_state(self):
        if self.measure == "Jaccard":
            shuffled_state = self.filtered_post_perturbed_state.copy()
            random_vals_jaccard = self.filtered_post_perturbed_state.copy()
            np.random.shuffle(random_vals_jaccard)
            for j, (ind, rand) in enumerate(zip(self.survived_matrix, random_vals_jaccard)):
                mask = np.ones(len(shuffled_state[j, :]), dtype=bool)
                mask[ind] = False
                shuffled_state[j, :][mask] = rand[mask]
            return self._normalize_cohort(shuffled_state)
        elif self.measure == "braycurtis":
            shuffled_state = self.post_perturbed_state.copy()
            random_vals_bc = self.post_perturbed_state.copy()
            np.random.shuffle(random_vals_bc)
            for j, (ind, rand) in enumerate(zip(self.survived_matrix, random_vals_bc)):
                mask = np.ones(len(shuffled_state[j, :]), dtype=bool)
                mask[ind] = False
                shuffled_state[j, :][mask] = rand[mask]
            return self._normalize_cohort(shuffled_state)

    def _calculate_Jaccard_similarity(self):
        jaccard_perturbed = []
        vals_real = []
        vals_shuffled = []
        for i in range(self.perturbed_state.shape[0] - 1):
            for j in range(i + 1, self.perturbed_state.shape[0]):
                overlap_object_perturbed = Overlap(self.perturbed_state[i, :],
                                                   self.perturbed_state[j, :],
                                                   overlap_type="Jaccard")

                jaccard_perturbed.append(overlap_object_perturbed.calculate_overlap())
                intrsection_ind = overlap_object_perturbed.intersection_bool
                mask = np.ones(len(intrsection_ind), dtype=bool)
                mask[intrsection_ind] = False

                if self.measure == "Jaccard":
                    try:
                        overlap_object_post_perturbed = Overlap(self.filtered_post_perturbed_state[i, :][mask],
                                                                self.filtered_post_perturbed_state[j, :][mask],
                                                                overlap_type="Jaccard")
                        vals_real.append(overlap_object_post_perturbed.calculate_overlap())
                        overlap_object_shuffled = Overlap(self.shuffled_state[i, :][mask],
                                                          self.shuffled_state[j, :][mask],
                                                          overlap_type="Jaccard")
                        vals_shuffled.append(overlap_object_shuffled.calculate_overlap())

                        size = np.shape(self.filtered_post_perturbed_state[i, :][mask])[0]
                        assert size > 25
                        assert np.shape(np.where(self.shuffled_state[j, :][mask] != 0))[1] / size > 0.1
                        assert np.shape(np.where(self.shuffled_state[i, :][mask] != 0))[1] / size > 0.1
                        assert np.shape(np.where(self.filtered_post_perturbed_state[i, :][mask] != 0))[1] / size > 0.1
                        assert np.shape(np.where(self.filtered_post_perturbed_state[j, :][mask] != 0))[1] / size > 0.1
                    except AssertionError:
                        print("Conditions to get unbiased Jaccard measure not satisfied.")
                elif self.measure == "braycurtis":
                    vals_real.append(braycurtis(self._normalize_cohort(self.post_perturbed_state[i, :][mask]),
                                                self._normalize_cohort(self.post_perturbed_state[j, :][mask])))
                    vals_shuffled.append(braycurtis(self._normalize_cohort(self.shuffled_state[i, :][mask]),
                                                    self._normalize_cohort(self.shuffled_state[j, :][mask])))
        return jaccard_perturbed, vals_real, vals_shuffled

    def get_results(self):
        results = {"growth_rate": self.r, "logistic_growth": self.s,
                   "interaction_matrix": self.A, "initial_conditions": self.y0,
                   "perturbed_state": self.perturbed_state,
                   "new_perturbed_state": self.new_perturbed_state,
                   "post_perturbed_state": self.post_perturbed_state,
                   "filtered_post_perturbed_state": self.filtered_post_perturbed_state,
                   "shuffled_state": self.shuffled_state,
                   "similarity_perturbed": self.jaccard_perturbed,
                   "real": self.vals_real,
                   "shuffled": self.vals_shuffled,
                   "survived_matrix": self.survived_matrix}
        return results

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized