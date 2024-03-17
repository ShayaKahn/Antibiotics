import numpy as np
import random
from scipy.spatial.distance import braycurtis, jensenshannon
from Null_model_functions import create_shuffled_samples

class NullModel:
    """
    This class is responsible to create the null model for a given subject,
     given a baseline sample, antibiotics sample and test_sample (post-antibiotics sample).
    """

    def __init__(self, baseline_sample, ABX_sample, test_sample, baseline_cohort, num_reals=100,
                 baseline_sample_sum=None, ABX_sample_sum=None, exclude_test=False):
        """
        :param baseline_sample: baseline sample of the subject, numpy array of shape (1, # species)
        :param ABX_sample: antibiotics sample of the subject, numpy array of shape (1, # species)
        :param test_sample: after antibiotics sample of the subject, numpy array of shape (1, # species)
        :param baseline_cohort: baseline cohort, numpy array of shape (# samples, # species)
        :param baseline_sample_sum: sum of the baseline samples, numpy array of shape (1, # species), if not provided,
               the find_ARS function will use the baseline_sample
        :param ABX_sample_sum: sum of the antibiotics samples, numpy array of shape (1, # species), if not provided,
               the find_ARS function will use the ARS_sample
        :param num_reals: number of realizations of the null model (number of shuffled samples)
        :param sparse: if True, use sparse matrices
        """
        self.exclude_test = exclude_test
        if not all(isinstance(arr, np.ndarray) and arr.ndim == 1 for arr in [baseline_sample, ABX_sample, test_sample]):
            raise TypeError("baseline_sample, ABX_sample, and test_sample must be 1D numpy arrays.")
        if not all(len(arr) == len(baseline_sample) for arr in [ABX_sample, test_sample]):
            raise ValueError("ABX_sample and test_sample must have the same size as baseline_sample.")
        if not (isinstance(baseline_cohort, np.ndarray) and baseline_cohort.ndim == 2 and
                baseline_cohort.shape[1] == len(baseline_sample)):
            raise TypeError("baseline_cohort must be a 2D numpy array with the same number of columns as"
                            "baseline_sample, ABX_sample, and test_sample.")
        self.baseline_sample = baseline_sample
        self.ABX_sample = ABX_sample
        if baseline_sample_sum is None:
            self.baseline_sample_sum = self.baseline_sample
            self.ABX_sample_sum = self.ABX_sample
        else:
            self.baseline_sample_sum = baseline_sample_sum
            self.ABX_sample_sum = ABX_sample_sum
        self.test_sample = test_sample
        self.baseline_cohort = baseline_cohort
        self.num_species = np.size(np.nonzero(self.test_sample))
        if not (isinstance(num_reals, int) and num_reals > 0):
            raise ValueError("num_reals must be an integer greater than 0.")
        self.num_reals = num_reals
        self.ARS, self.non_ARS = self._find_ARS()
        self.shuffled_samples = self._create_shuffled_samples()

    def _find_ARS(self):
        """
        Find the ARS and non-ARS species
        :return: Indexes of the ARS species, i.e, the species that present in the baseline, ABX and future samples.
         Indexes of the non-ARS species, i.e, the species that are not ARS species.
        """
        nonzero_ABX = np.nonzero(self.ABX_sample_sum)  # find nonzero elements of ABX sample
        nonzero_base = np.nonzero(self.baseline_sample_sum)  # find nonzero elements of baseline sample
        nonzero_test = np.nonzero(self.test_sample)  # find nonzero elements of the test sample
        intersect_ABX_base = np.intersect1d(nonzero_ABX, nonzero_base)  # find species present in ABX and baseline
        intersect_ABX_test = np.intersect1d(nonzero_ABX, nonzero_test)  # find species present in ABX and test
        if self.exclude_test:
            ARS = intersect_ABX_base
            non_ARS = np.setdiff1d(np.arange(np.size(self.baseline_sample_sum)), ARS)
        else:
            ARS = np.intersect1d(intersect_ABX_base, intersect_ABX_test)  # find ARS species by definition
            non_ARS = np.setdiff1d(np.arange(np.size(self.baseline_sample_sum)), ARS)  # find non-ARS species
        column_sums = self.baseline_cohort.sum(axis=0)  # sum of each column of the baseline cohort
        zero_sum_column_indexes = np.where(column_sums == 0)[0]  # find species that are not present in the baseline
        # remove species that are not present in the baseline cohort from the non-ARS species
        non_ARS = np.setdiff1d(non_ARS, zero_sum_column_indexes)
        return ARS, non_ARS

    def _create_shuffled_samples(self):
        return create_shuffled_samples(self.num_reals, self.test_sample.astype(np.float64), self.non_ARS.astype(np.int64),
                                       self.baseline_cohort.astype(np.float64), self.ARS.astype(np.int64))

    #def _create_shuffled_samples(self):
        # initialize shuffled samples
    #    shuffled_samples = np.zeros((self.num_reals, np.size(self.test_sample)))

        # set the constraint on the number of species
    #    stop = np.size(np.nonzero(self.test_sample[self.non_ARS]))

        # define the pool of species to be shuffled not including the ARS species
    #    pool = self.baseline_cohort[self.baseline_cohort != 0]
    #    pool_indices = np.where(self.baseline_cohort != 0)[1]
    #    index_to_remove = [i for i, ind in enumerate(pool_indices) if ind in self.ARS]
    #    pool_indices = np.array([ind for ind in pool_indices if ind not in self.ARS])
    #    pool = np.array([pool[i] for i in range(len(pool)) if i not in index_to_remove])

    #    for shuffled_sample in shuffled_samples:
    #        pool_copy = pool.copy()
    #        pool_indices_copy = pool_indices.copy()
    #        counter = 0
    #        while counter < stop:
    #            index = random.randint(0, len(pool_indices_copy) - 1)
    #            shuffled_sample[pool_indices_copy[index]] = pool_copy[index]
    #            similar_specie = np.where(pool_indices_copy == pool_indices_copy[index])[0]
    #            pool_copy = np.delete(pool_copy, similar_specie)
    #            pool_indices_copy = np.delete(pool_indices_copy, similar_specie)
    #            counter += 1

        # normalize the shuffled_samples
    #    shuffled_samples = self._normalize_cohort(shuffled_samples)

    #    return shuffled_samples

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized

    def distance(self, method="Bray Curtis"):
        """
        :param method: distance method, either "Bray Curtis" or "rJSD".
        :return: distance between the test sample and the baseline sample and the distances between the shuffled samples
                 and the baseline sample.
        """
        if method not in ["Bray Curtis", "rJSD"]:
            raise ValueError("method must be either 'Bray Curtis' or 'rJSD'.")
        if method == "Bray Curtis":
            bc_real = braycurtis(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
                                 self._normalize_cohort(self.test_sample[self.non_ARS]))
            bc_shuffled = np.zeros(self.num_reals)
            for i in range(self.num_reals):
                bc_shuffled[i] = braycurtis(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
                                            self.shuffled_samples[:, self.non_ARS][i, :])
            return bc_real, bc_shuffled
        elif method == "rJSD":
            bc_real = jensenshannon(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
                                    self._normalize_cohort(self.test_sample[self.non_ARS]))
            bc_shuffled = np.zeros(self.num_reals)
            for i in range(self.num_reals):
                bc_shuffled[i] = jensenshannon(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
                                               self.shuffled_samples[:, self.non_ARS][i, :])
            return bc_real, bc_shuffled