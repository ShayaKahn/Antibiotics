from Classes.overlap import Overlap
import numpy as np
import copy

class StochasticSurrogate:
    """
    This class implements stochastic surrogate data analysis
    """
    def __init__(self, test_base_sample, base_samples_collections, test_post_ABX_sample, test_ABX_sample, reals,
                 weighted=False):
        """
        :param test_base_sample: dictionary that contain a numpy vector of shape (# species, ) that represents the
         baseline sample of a test subject. The key is subjects identifier.
        :param base_samples_collections: dictionary of numpy vectors of shape (# species,) that represent the baseline
         of different subjects in the experiment. The keys are the identifiers of the subjects.
        :param test_post_ABX_sample: numpy array of shape (# species, ) that contains the post antibiotics sample of the
         test subject.
        :param test_ABX_sample: numpy array of shape (# species, ) that contains the antibiotics sample of the
         test subject.
        :param reals: the number of realizations of the stochastic surrogate data.
        :param weighted: boolean that indicates if the weighted version of the similarity is used.
        """
        # separate the key and the matrix from the dictionary
        if not (isinstance(test_base_sample, dict)):
            raise ValueError("test_base_sample must be dictionary")
        self.test_key, self.test_base_sample = next(iter(test_base_sample.items()))
        if not (isinstance(self.test_key, str) and isinstance(self.test_base_sample, np.ndarray)):
            raise ValueError("test_base_sample must be dictionary with string as a key and numpy array of"
                             " shape (# species,) as value.")
        if not (isinstance(base_samples_collections, dict)):
            raise ValueError("base_samples_collections must be dictionary")
        self.base_samples_collections = base_samples_collections
        if not (isinstance(test_post_ABX_sample, np.ndarray)):
            raise ValueError("test_samples must be numpy array of shape (# species,)")
        if not (test_post_ABX_sample.shape == test_ABX_sample.shape):
            raise ValueError("The dimensions are not consistent")
        if not (isinstance(reals, int) and reals > 0):
            raise ValueError("reals must be an integer greater than 0")
        self.reals = reals
        if not (isinstance(weighted, bool)):
            raise ValueError("weighted must be boolean")
        self.weighted = weighted
        self.test_post_ABX_sample = test_post_ABX_sample
        self.test_ABX_sample = test_ABX_sample
        # find the indexes of the ARS and the non-ARS species
        self.non_ARS, self.ARS = self._find_non_ARS()
        # Generate the stochastic baseline sampes
        self.stochastic_base_samples_collections = self._generate_stochastic_surrogate()

    def _find_non_ARS(self):
        # find the indexes of the ARS species, which are the species that are present in both the baseline and ABX
        # samples
        cols_with_zeros = np.any(np.vstack([self.test_ABX_sample, self.test_base_sample]) == 0, axis=0)
        indexes = np.where(cols_with_zeros)[0]
        indexes_comp = np.where(~cols_with_zeros)[0]
        return indexes, indexes_comp

    def _generate_stochastic_surrogate(self):
        base_samples_collections_copy = copy.deepcopy(self.base_samples_collections)
        stochastic_base_samples_collections = {}
        # find Delta
        size_test = np.size(np.nonzero(self.test_base_sample)[0])
        size_ARS_test = np.size(self.ARS)
        delta = size_ARS_test / size_test

        test_base_sample_copy = self.test_base_sample.copy()
        test_base_sample_copy[self.ARS] = 0
        for key in base_samples_collections_copy:
            # find Gamma
            size = np.size(np.nonzero(base_samples_collections_copy[key].squeeze())[0])
            size_ARS = np.size(np.nonzero(base_samples_collections_copy[key].squeeze()[self.ARS]))
            gamma = int(size * delta - size_ARS)
            if gamma < 0:
                gamma = 0
            # set Gamma species to zero
            base_samples_collections_copy[key][self.ARS] = 0
            stochastic_base_samples = np.tile(base_samples_collections_copy[key], (self.reals, 1))
            indexes = np.where(base_samples_collections_copy[key] != 0)[0]
            for j in range(stochastic_base_samples.shape[0]):
                random_values = np.random.choice(indexes, gamma, replace=False)
                stochastic_base_samples[j, random_values] = 0
            # normalize the stochastic baseline samples
            stochastic_base_samples_collections[key] = self._normalize_cohort(stochastic_base_samples)
        return stochastic_base_samples_collections

    def apply_surrogate_data_analysis(self, method="Szymkiewicz Simpson first"):
        # calculate results using overlap class
        results = {}
        if self.weighted:
            if method == "Szymkiewicz Simpson first":
                # calculate weighted specificity
                measure = self._specificity_weighted(self.test_base_sample)
                results[self.test_key] = measure
                for key in self.stochastic_base_samples_collections.keys():
                    measure_surrogate = np.mean([self._specificity_weighted(smp) for smp in
                                                 self.stochastic_base_samples_collections[key]])
                    results[key] = measure_surrogate
            elif method == "Szymkiewicz Simpson second":
                # calculate weighted recovery
                measure = self._recovery_weighted(self.test_base_sample)
                results[self.test_key] = measure
                for key in self.stochastic_base_samples_collections.keys():
                    measure_surrogate = np.mean([self._recovery_weighted(smp) for smp in self.stochastic_base_samples_collections[key]])
                    results[key] = measure_surrogate
        else:
            # Calculate unweighted specificity/ recovery
            measure = Overlap(self.test_base_sample[self.non_ARS], self.test_post_ABX_sample[self.non_ARS],
                              overlap_type=method).calculate_overlap()
            # store the results of the test subject
            results[self.test_key] = measure
            # calculate the unweighted specificity/ recovery of the surrogate data
            for key in self.stochastic_base_samples_collections.keys():
                measure_surrogate = np.mean([Overlap(smp[self.non_ARS], self.test_post_ABX_sample[self.non_ARS],
                                              overlap_type=method).calculate_overlap() for smp in
                                              self.stochastic_base_samples_collections[key]])
                results[key] = measure_surrogate
        return results

    def _specificity_weighted(self, sample_base):
        # find the sum of the post ABX sample
        sum_post_non_ARS = np.sum(self.test_post_ABX_sample[self.non_ARS])
        # find the intersection to the baseline
        intersect = np.where(np.logical_and(sample_base[self.non_ARS] != 0,
                                            self.test_post_ABX_sample[self.non_ARS] != 0))
        sum_post = np.sum(self.test_post_ABX_sample[self.non_ARS][intersect])
        return sum_post / sum_post_non_ARS

    def _recovery_weighted(self, sample_base):
        # find the sum of the baseline sample
        sum_base_non_ARS = np.sum(sample_base[self.non_ARS])
        # find the intersection to the baseline
        intersect = np.where(np.logical_and(sample_base[self.non_ARS] != 0,
                                            self.test_post_ABX_sample[self.non_ARS] != 0))
        sum_base = np.sum(sample_base[self.non_ARS][intersect])
        return sum_base / sum_base_non_ARS

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized