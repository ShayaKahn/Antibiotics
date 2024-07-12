import numpy as np
import random
from Classes.GLV_model import Glv
from scipy.stats import skewnorm
from Classes.overlap import Overlap
import networkx as nx

class HC:
    """
    This class is responsible for the simulation of the historical contingency model.
    """
    def __init__(self, num_samples, pool_size, num_survived_min, num_survived_max, mean, sigma, lower, higher, delta,
                 final_time, max_step, epsilon, threshold, min_skew, max_skew, interaction_strength, max_growth,
                 min_growth, method='RK45', measure='Jaccard'):
        """
        :param num_samples: The number of samples.
        :param pool_size: The total size of the population.
        :param num_survived_max: The maximal number of survived species.
        :param num_survived_min: The minimal number of survived species.
        :param mean: For the interaction matrix generation, the mean of the normal distribution.
        :param sigma: For the interaction matrix generation, The standard deviation of the normal distribution.
        :param lower: For the probability vector generation(sparsity of the interaction matrix),
                      the lower bound of the uniform distribution.
        :param higher: For the probability vector generation(sparsity of the interaction matrix),
                       the higher bound of the uniform distribution.
        :param delta: The stop condition for the steady state.
        :param final_time: The final time of the integration.
        :param max_step: The maximal allowed step size.
        :param epsilon: The value to insert in the non-survived species.
        :param threshold: The threshold to remove low abundances.
        :param min_skew: The minimum value for the skewness of the interaction matrix.
        :param max_skew: The maximum value for the skewness of the interaction matrix.
        :param interaction_strength: The strength of the interaction.
        :param method: The method to solve the GLV model.
        :param measure: The measure to calculate the similarity between the samples.
        :param max_growth: The maximum growth rate.
        :param min_growth: The minimum growth rate.
        """
        if not (isinstance(num_samples, int) and isinstance(pool_size, int) and
                isinstance(num_survived_min, int) and isinstance(num_survived_max, int)):
            raise ValueError("num_samples, pool_size, num_survived_min, and num_survived_max must be integers.")
        if num_samples <= 0 or pool_size <= 0 or num_survived_min <= 0 or num_survived_max <= 0:
            raise ValueError("num_samples, pool_size, num_survived_min and num_survived_max must be greater than 0.")
        if num_survived_max > pool_size or num_survived_min > num_survived_max:
            raise ValueError("num_survived_max must be smaller then pool_size and num_survived_min must be smaller then"
                             " num_survived_max.")
        self.num_samples = num_samples
        self.pool_size = pool_size
        self.num_survived_min = num_survived_min
        self.num_survived_max = num_survived_max
        self.num_survived_list = self._create_num_survived_list()
        if not (isinstance(mean, float) or isinstance(mean, int)):
            raise ValueError("mean must be integer or float.")
        if not (isinstance(sigma, float) or isinstance(sigma, int)):
            raise ValueError("sigma must be integer or float.")
        self.mean = mean
        self.sigma = sigma
        if method not in ['RK45', 'BDF', 'RK23', 'Radau', 'LSODA', 'DOP853']:
            raise ValueError("method must be one of the following methods: RK45, BDF, RK23, Radau, LSODA, and DOP853")
        self.method = method
        if not (isinstance(lower, float) and isinstance(higher, float)):
            raise ValueError("lower and higher must be floats.")
        if (lower < 0 or lower >= 1) or (higher <= 0 or higher > 1):
            raise ValueError("lower must be in the range (0, 1] and higher must be in the range [0, 1).")
        if lower >= higher:
            raise ValueError("lower must be smaller than higher.")
        self.lower = lower
        self.higher = higher
        if not (isinstance(threshold, float) or isinstance(threshold, int)):
            raise ValueError("threshold must be a number.")
        if threshold <= 0:
            raise ValueError("threshold must be greater than 0.")
        self.threshold = threshold
        if not (isinstance(min_skew, float) or isinstance(min_skew, int)):
            raise ValueError("min_skew must be a number.")
        if not (isinstance(max_skew, float) or isinstance(max_skew, int)):
            raise ValueError("max_skew must be a number.")
        if min_skew >= max_skew:
            raise ValueError("min_skew must be smaller than max_skew.")
        self.min_skew = min_skew
        self.max_skew = max_skew
        if not (isinstance(interaction_strength, float) or isinstance(interaction_strength, int)):
            raise ValueError("interaction_strength must be a number.")
        if interaction_strength <= 0:
            raise ValueError("interaction_strength must be greater than 0.")
        self.interaction_strength = interaction_strength
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
        if measure not in ['Jaccard', 'intersection']:
            raise ValueError("measure must be one of the following measures: Jaccard and intersection.")
        self.measure = measure
        if not (isinstance(max_growth, (float, int)) and 0 < max_growth <= 1):
            raise ValueError("max_growth must be a number between 0 and 1.")
        if not (isinstance(min_growth, (float, int)) and 0 <= min_growth < 1):
            raise ValueError("min_growth must be a number between 0 and 1.")
        if min_growth > max_growth:
            raise ValueError("min_growth must be smaller than max_growth.")
        self.max_growth = max_growth
        self.min_growth = min_growth
        self.r = self._set_growth_rate()
        self.s = self._set_logistic_growth()
        self.p_vector = self._set_p_vector()
        self.y0, self.survived_matrix = self._set_initial_conditions()
        self.A = self._set_interaction_matrix()#self.A = self._set_interaction_matrix()
        self.steady_state = np.linalg.solve(self.A+np.diag(self.s), -self.r)
        self.perturbed_state, event_not_satisfied_ind = self._apply_GLV(self.y0, norm=True)
        self.new_perturbed_state = self._normalize_cohort(self._insert_total_pool())
        self.post_perturbed_state, event_not_satisfied_ind_post = self._apply_GLV(self.new_perturbed_state, norm=True)
        self.filtered_post_perturbed_state = self._remove_low_abundances()
        unique_event_not_satisfied_ind = sorted(list(set(event_not_satisfied_ind + event_not_satisfied_ind_post)))
        self.perturbed_state = np.delete(self.perturbed_state, unique_event_not_satisfied_ind, axis=0)
        self.filtered_post_perturbed_state = np.delete(self.filtered_post_perturbed_state,
                                                       unique_event_not_satisfied_ind, axis=0)
        self.survived_matrix = [item for idx, item in enumerate(
            self.survived_matrix) if idx not in unique_event_not_satisfied_ind]#np.delete(self.survived_matrix, unique_event_not_satisfied_ind, axis=0)
        ####################################
        self.shuffled_state = self._generate_shuffled_state()
        self.jaccard_perturbed, self.jaccard_real, self.jaccard_shuffled = self._calculate_Jaccard_similarity()

    def _set_initial_conditions(self):
        y0 = np.zeros((self.num_samples, self.pool_size))
        survived_matrix = [random.sample(range(0, self.pool_size),
                                         self.num_survived_list[i]) for i in range(self.num_samples)]
        for index, (y, s) in enumerate(zip(y0, survived_matrix)):
            y[s] = np.random.rand(1, self.num_survived_list[index])
        return y0, survived_matrix

    def _set_interaction_matrix(self):
        # Create a num_species x num_species matrix filled with zeros
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))

        ###############################################
        #random_numbers = np.zeros((self.pool_size, self.pool_size))
        # Generate random numbers
        #skew_vals = np.random.uniform(self.min_skew, self.max_skew, self.pool_size)
        #for i in range(self.pool_size):
            #random_numbers[i, :] = skewnorm.rvs(a=skew_vals[i], loc=self.mean,
                                                #scale=self.sigma, size=self.pool_size)*self.interaction_strength
        ###############################################

        ###############################################
        #upper_triangle = -np.abs(np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size)))*self.interaction_strength
        # Ensure the matrix is symmetric by copying the upper triangle to the lower triangle
        #random_numbers = np.triu(upper_triangle) + np.triu(upper_triangle, 1).T
        # Apply negative absolute value to all elements
        #upper_triangle_indices = np.triu_indices(self.pool_size)
        #upper_triangle_flattened = random_numbers[upper_triangle_indices]

        # Calculate how many elements to set to zero (50% of the upper triangle, including the diagonal)
        #num_zeros = len(upper_triangle_flattened) // 2

        # Randomly choose indices to set to zero
        #zero_indices = np.random.choice(range(len(upper_triangle_flattened)), size=num_zeros, replace=False)

        # Set chosen elements to zero
        #upper_triangle_flattened[zero_indices] = 0

        # Update the matrix with the modified upper triangle
        #random_numbers[upper_triangle_indices] = upper_triangle_flattened

        # Mirror the upper triangle to the lower triangle to maintain symmetry
        #interaction_matrix = np.triu(random_numbers) + np.triu(random_numbers, 1).T
        ###############################################


        ###############################################
        random_numbers = -np.abs(np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size)))*self.interaction_strength
        # Create a mask based on the probability p
        mask = np.random.rand(self.pool_size, self.pool_size) < self.p_vector
        # Apply the mask to the random numbers
        interaction_matrix[mask] = random_numbers[mask]
        print(np.linalg.eigvals(interaction_matrix))
        ###############################################
        return interaction_matrix

    def _set_symmetric_interaction_matrix(self):
        G = nx.binomial_graph(self.pool_size, 0.5)
        mask = nx.to_numpy_array(G).astype(bool)
        random_numbers = -np.abs(np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size))) * self.interaction_strength
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))
        interaction_matrix[mask] = random_numbers[mask]
        return interaction_matrix

    def _set_growth_rate(self):
        return np.random.uniform(self.min_growth, self.max_growth, self.pool_size)#np.random.lognormal(mean=0.0, sigma=0.1, size=self.pool_size)

    def _set_logistic_growth(self):
        return np.ones(self.pool_size)

    def _set_p_vector(self):
        return np.random.uniform(self.lower, self.higher, (1, self.pool_size))

    def _create_num_survived_list(self):
        return np.random.randint(self.num_survived_min, self.num_survived_max, self.num_samples)

    def _apply_GLV(self, init_cond, norm):
        glv_object = Glv(self.num_samples, self.pool_size, self.delta, self.r, self.s, self.A, init_cond,
                         self.final_time, self.max_step, normalize=norm, method=self.method)
        final_abundances = glv_object.solve()
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

    #def _generate_shuffled_state(self):
    #    shuffled_state = self.filtered_post_perturbed_state.copy()
    #    for sample, ind in zip(shuffled_state, self.survived_matrix):
    #        random_vals = np.array([np.random.choice(self.filtered_post_perturbed_state[:, col]) for
    #                                col in range(self.filtered_post_perturbed_state.shape[1])])
    #        mask = np.ones(len(sample), dtype=bool)
    #        mask[ind] = False
    #        sample[mask] = random_vals[mask]
    #    return self._normalize_cohort(shuffled_state)

    def _generate_shuffled_state(self):
        shuffled_state = self.post_perturbed_state.copy()
        for sample, ind in zip(shuffled_state, self.survived_matrix):
            random_vals = np.array([np.random.choice(self.post_perturbed_state[:, col]) for
                                    col in range(self.post_perturbed_state.shape[1])])
            mask = np.ones(len(sample), dtype=bool)
            mask[ind] = False
            sample[mask] = random_vals[mask]
        return self._normalize_cohort(shuffled_state)

    def _calculate_Jaccard_similarity(self):
        jaccard_perturbed = []
        jaccard_real = []
        jaccard_shuffled = []
        for i in range(self.perturbed_state.shape[0] - 1):
            for j in range(i + 1, self.perturbed_state.shape[0]):
                overlap_object_perturbed = Overlap(self.perturbed_state[i, :],
                                                   self.perturbed_state[j, :],
                                                   overlap_type="Jaccard")
                if self.measure == "Jaccard":
                    jaccard_perturbed.append(overlap_object_perturbed.calculate_overlap())
                elif self.measure == "intersection":
                    overlap_object_perturbed.calculate_overlap()
                    jaccard_perturbed.append(np.sum(overlap_object_perturbed.intersection_bool) / len(
                        overlap_object_perturbed.intersection_bool))


                intrsection_ind = overlap_object_perturbed.intersection_bool

                mask = np.ones(len(intrsection_ind), dtype=bool)
                mask[intrsection_ind] = False

                #overlap_object_post_perturbed = Overlap(self.filtered_post_perturbed_state[i, :][mask],
                #                                        self.filtered_post_perturbed_state[j, :][mask],
                #                                        overlap_type="Overlap")#overlap_type="Jaccard")

                #overlap_object_post_perturbed = Overlap(self.post_perturbed_state[i, :],
                #                                        self.post_perturbed_state[j, :],
                #                                        overlap_type="Overlap")  # overlap_type="Jaccard")

                overlap_object_post_perturbed = Overlap(self.post_perturbed_state[i, :][mask],
                                                        self.post_perturbed_state[j, :][mask],
                                                        overlap_type="Overlap")#overlap_type="Jaccard")

                if self.measure == "Jaccard":
                    from scipy.spatial.distance import braycurtis
                    #jaccard_real.append(overlap_object_post_perturbed.calculate_overlap())
                    jaccard_real.append(braycurtis(self.post_perturbed_state[i, :],
                                                   self.post_perturbed_state[j, :]))
                elif self.measure == "intersection":
                    overlap_object_post_perturbed.calculate_overlap()
                    jaccard_real.append(np.sum(overlap_object_post_perturbed.intersection_bool) / len(
                        overlap_object_post_perturbed.intersection_bool))

                overlap_object_shuffled = Overlap(self._normalize_cohort(self.shuffled_state[i, :][mask]),
                                                  self._normalize_cohort(self.shuffled_state[j, :][mask]),
                                                  overlap_type="Overlap")#overlap_type="Jaccard")

                #overlap_object_shuffled = Overlap(self.shuffled_state[i, :],
                #                                  self.shuffled_state[j, :],
                #                                  overlap_type="Overlap")  # overlap_type="Jaccard")

                if self.measure == "Jaccard":
                    #jaccard_shuffled.append(overlap_object_shuffled.calculate_overlap())
                    jaccard_shuffled.append(braycurtis(self._normalize_cohort(self.shuffled_state[i, :]),
                                                       self._normalize_cohort(self.shuffled_state[j, :])))
                elif self.measure == "intersection":
                    overlap_object_shuffled.calculate_overlap()
                    jaccard_shuffled.append(np.sum(overlap_object_shuffled.intersection_bool) / len(
                        overlap_object_shuffled.intersection_bool))

        return jaccard_perturbed, jaccard_real, jaccard_shuffled

    def get_results(self):
        results = {"growth_rate": self.r, "logistic_growth": self.s,
                   "probability_vector": self.p_vector,
                   "interaction_matrix": self.A, "initial_conditions": self.y0,
                   "steady_state": self.steady_state,
                   "perturbed_state": self.perturbed_state,
                   "new_perturbed_state": self.new_perturbed_state,
                   "post_perturbed_state": self.post_perturbed_state,
                   "filtered_post_perturbed_state": self.filtered_post_perturbed_state,
                   "shuffled_state": self.shuffled_state,
                   "similarity_perturbed": self.jaccard_perturbed,
                   "similarity_real": self.jaccard_real,
                   "similarity_shuffled": self.jaccard_shuffled,
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
