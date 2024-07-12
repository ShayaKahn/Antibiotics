import numpy as np
import random
from Classes.GLV_model import Glv
from Classes.overlap import Overlap
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import braycurtis
from scipy.stats import powerlaw
import networkx as nx

class HC:
    """
    This class is responsible for the simulation of the historical contingency model.
    """
    def __init__(self, num_samples, pool_size, num_survived_min, num_survived_max, mean, sigma, c, delta, final_time,
                 max_step, epsilon, threshold, min_growth, max_growth, rep, symmetric=True, alpha=None, method='RK45',
                 multiprocess=True, weighted=False):
        """
        :param num_samples: The number of samples.
        :param pool_size: The total size of the population.
        :param num_survived_max: The maximal number of survived species.
        :param num_survived_min: The minimal number of survived species.
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
        :param multiprocess: If True, the class will use the multiprocessing module.
        :param rep: The number of repetitions to the shuffling procedure.
        :param symmetric: If True, the interaction matrix will be symmetric.
        #:param weighted: If True, the Jaccard similarity will be weighted.
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
        try:
            assert min(self.num_survived_list) >= self.num_survived_min, \
                "The minimal number of survived species is less than the expected minimum!"
        except AssertionError as e:
            print(f"Assertion Error: {e}")
        try:
            assert max(self.num_survived_list) <= self.num_survived_max, \
                "The maximal number of survived species is greater than the expected maximum!"
        except AssertionError as e:
            print(f"Assertion Error: {e}")
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
        if not (isinstance(rep, int) and rep > 0):
            raise ValueError("rep must be integer and greater then 0.")
        self.rep = rep
        if not isinstance(symmetric, bool):
            raise ValueError("symmetric must be of type bool.")
        self.symmetric = symmetric
        if not isinstance(weighted, bool):
            raise ValueError("weighted must be of type bool.")
        self.weighted = weighted
        self.multiprocess = multiprocess
        self.r = self._set_growth_rate()
        self.s = self._set_logistic_growth()
        self.y0, self.survived_matrix = self._set_initial_conditions()
        if self.symmetric:
            self.N = self._set_symmetric_interaction_matrix()
            try:
                mat = self.N.copy()
                index = np.where(mat != 0)
                mat[index] = 1
                assert np.array_equal(mat, mat.T)
            except AssertionError:
                print("The interaction matrix is not symmetric.")
        else:
            self.N = self._set_interaction_matrix()
        self.H = self._set_int_strength_matrix()
        self.A = np.dot(self.N, self.H)
        self.perturbed_state, event_not_satisfied_ind = self._apply_GLV(self.y0, norm=True)
        self.perturbed_state[self.perturbed_state < self.epsilon] = 0
        self.perturbed_state = self._normalize_cohort(self.perturbed_state)
        for j, sample in enumerate(self.perturbed_state):
            index = list(np.where(sample != 0)[0])
            self.survived_matrix[j] = [x for x in self.survived_matrix[j] if x in index]
        self.new_perturbed_state = self._normalize_cohort(self._insert_total_pool())
        self.post_perturbed_state, event_not_satisfied_ind_post = self._apply_GLV(self.new_perturbed_state, norm=True)
        self.filtered_post_perturbed_state = self._remove_low_abundances()
        unique_event_not_satisfied_ind = sorted(list(set(event_not_satisfied_ind + event_not_satisfied_ind_post)))
        self.perturbed_state = np.delete(self.perturbed_state, unique_event_not_satisfied_ind, axis=0)
        self.filtered_post_perturbed_state = np.delete(self.filtered_post_perturbed_state,
                                                       unique_event_not_satisfied_ind, axis=0)
        for index in sorted(unique_event_not_satisfied_ind, reverse=True):
            del self.survived_matrix[index]
        self.shuffled_state = self._generate_shuffled_state()
        self.weights = self._generate_weights()
        (self.jaccard_perturbed, self.jaccard_post_perturbed,
         self.jaccard_shuffled_matrices) = self._calculate_Jaccard_similarity()
        #self.bc_post_perturbed, self.bc_shuffled_matrices = self._calculate_bc_distance()

    def _set_initial_conditions(self):
        y0 = np.zeros((self.num_samples, self.pool_size))
        survived_matrix = [random.sample(range(0, self.pool_size),
                                         self.num_survived_list[i]) for i in range(self.num_samples)]
        for index, (y, s) in enumerate(zip(y0, survived_matrix)):
            y[s] = np.random.rand(1, self.num_survived_list[index])
        return y0, survived_matrix

    def _set_interaction_matrix(self):
        interaction_matrix = np.zeros((self.pool_size, self.pool_size))
        random_numbers = np.random.normal(self.mean, self.sigma, size=(self.pool_size, self.pool_size))
        mask = np.random.rand(self.pool_size, self.pool_size) < self.c
        interaction_matrix[mask] = random_numbers[mask]
        np.fill_diagonal(interaction_matrix, 0)
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
        np.fill_diagonal(interaction_matrix, 0)
        return -np.abs(interaction_matrix)

    def _set_growth_rate(self):
        return np.random.uniform(self.min_growth, self.max_growth, self.pool_size)

    def _set_logistic_growth(self):
        return np.ones(self.pool_size)

    def _create_num_survived_list(self):
        return np.random.randint(self.num_survived_min, self.num_survived_max, self.num_samples)

    def _generate_weights(self):
        init_cond = np.random.rand(1, self.pool_size)
        glv_object = Glv(1, self.pool_size, self.delta, self.r, self.s, self.A, init_cond,
                         self.final_time, self.max_step, normalize=True, method=self.method,
                         multiprocess=self.multiprocess)
        weights, _ = glv_object.solve()
        return weights.squeeze()

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
        shuffled_state = self.filtered_post_perturbed_state.copy()
        extended_shuffled_state = np.tile(shuffled_state, (self.rep, 1))
        np.apply_along_axis(np.random.shuffle, 0, extended_shuffled_state)
        return extended_shuffled_state

    def _calculate_Jaccard_similarity(self):
        jaccard_perturbed = []
        jaccard_post_perturbed = []
        jaccard_shuffled_matrices = []
        for i in range(self.perturbed_state.shape[0] - 1):
            for j in range(i + 1, self.perturbed_state.shape[0]):
                if self.weighted:
                    weighted_jaccard, intrsection_ind = self._weighted_jaccard(self.perturbed_state[i, :],
                                                                               self.perturbed_state[j, :],
                                                                               self.weights)
                    jaccard_perturbed.append(weighted_jaccard)
                else:
                    overlap_object_perturbed = Overlap(self.perturbed_state[i, :],
                                                       self.perturbed_state[j, :],
                                                       overlap_type="Jaccard")

                    jaccard_perturbed.append(overlap_object_perturbed.calculate_overlap())
                    intrsection_ind = overlap_object_perturbed.intersection_bool

                    mask = np.ones(len(intrsection_ind), dtype=bool)
                    mask[intrsection_ind] = False
                    overlap_object_post_perturbed = Overlap(self.filtered_post_perturbed_state[i, :][mask],
                                                            self.filtered_post_perturbed_state[j, :][mask],
                                                            overlap_type="Jaccard")
                    jaccard_post_perturbed.append(overlap_object_post_perturbed.calculate_overlap())
                    jaccard_shuffled_matrices.append(self._calculate_jaccard_similarity_matrix(self.shuffled_state,
                                                                                               mask))
        return jaccard_perturbed, jaccard_post_perturbed, jaccard_shuffled_matrices

    #def _calculate_bc_distance(self):
    #    bc_post_perturbed = []
    #    bc_shuffled_matrices = []
    #    for i in range(self.perturbed_state.shape[0] - 1):
    #        for j in range(i + 1, self.perturbed_state.shape[0]):
    #            overlap_object_perturbed = Overlap(self.perturbed_state[i, :],
    #                                               self.perturbed_state[j, :],
    #                                               overlap_type="Jaccard")
    #            overlap_object_perturbed.calculate_overlap()
    #            intrsection_ind = overlap_object_perturbed.intersection_bool
    #
    #            mask = np.ones(len(intrsection_ind), dtype=bool)
    #            mask[intrsection_ind] = False
    #            bc = braycurtis(self._normalize_cohort(self.filtered_post_perturbed_state[i, :][mask]),
    #                            self._normalize_cohort(self.filtered_post_perturbed_state[j, :][mask]))
    #            bc_post_perturbed.append(bc)
    #            bc_shuffled_matrices.append(self._calculate_bc_matrix(self.shuffled_state, mask))
    #    return bc_post_perturbed, bc_shuffled_matrices

    def get_results(self):
        results = {"growth rate": self.r, "logistic growth": self.s,
                   "interaction matrix": self.A, "initial conditions": self.y0,
                   "perturbed state": self.perturbed_state,
                   "new perturbed state": self.new_perturbed_state,
                   "post perturbed state": self.post_perturbed_state,
                   "filtered post perturbed state": self.filtered_post_perturbed_state,
                   "shuffled state": self.shuffled_state,
                   "survived matrix": self.survived_matrix,
                   "jaccard perturbed": self.jaccard_perturbed,
                   "jaccard post perturbed": self.jaccard_post_perturbed,
                   "jaccard shuffled matrices": self.jaccard_shuffled_matrices,}
                   #"braycurtis post perturbed": self.bc_post_perturbed,
                   #"braycurtis shuffled matrices": self.bc_shuffled_matrices}
        return results

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized

    @staticmethod
    def _weighted_jaccard(first_sample, second_sample, weights):
        smaple_first_boolean = np.where(first_sample != 0, 1, 0)
        smaple_second_boolean = np.where(second_sample != 0, 1, 0)
        intersection_bool = np.logical_and(smaple_first_boolean, smaple_second_boolean)
        intersection_binary = intersection_bool.astype(int)
        union_binary = np.logical_or(smaple_first_boolean, smaple_second_boolean).astype(int)
        weighted_jaccard = np.divide(np.dot(weights, intersection_binary), np.dot(weights, union_binary))
        return weighted_jaccard, intersection_bool

    @staticmethod
    def _calculate_jaccard_similarity_matrix(cohort, mask):
        binary_matrix = cohort[:, mask].copy()
        ind = np.where(binary_matrix != 0)
        binary_matrix[ind] = 1
        jaccard_similarity_mat = 1 - pairwise_distances(binary_matrix.astype(bool), metric='jaccard')
        return jaccard_similarity_mat

    def _calculate_bc_matrix(self, cohort, mask):
        mat = cohort[:, mask].copy()
        mat = self._normalize_cohort(mat)
        bc_mat = pairwise_distances(mat, metric='braycurtis')
        return bc_mat