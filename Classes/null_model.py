import numpy as np
import random
from scipy.spatial.distance import braycurtis, jensenshannon
#from Null_model_functions import create_shuffled_samples
from Classes.overlap import Overlap
from skbio.diversity import beta_diversity

class NullModel:
    """
    This class is responsible to create the null model for a given subject,
     given a baseline sample, antibiotics sample and test_sample (post-antibiotics sample).
    """
    def __init__(self, baseline_sample, ABX_sample, test_sample, baseline_cohort, num_reals=100):
        """
        :param baseline_sample: baseline sample of the subject, numpy array of shape (1, # species)
        :param ABX_sample: antibiotics sample of the subject, numpy array of shape (1, # species)
        :param test_sample: after antibiotics sample of the subject, numpy array of shape (1, # species)
        :param baseline_cohort: baseline cohort, numpy array of shape (# samples, # species)
        :param num_reals: number of realizations of the null model (number of shuffled samples)
        """
        if not all(isinstance(arr, np.ndarray) and arr.ndim == 1 for arr in [baseline_sample, ABX_sample, test_sample]):
            raise TypeError("baseline_sample, ABX_sample, and test_sample must be 1D numpy arrays.")
        if not all(len(arr) == len(baseline_sample) for arr in [ABX_sample, test_sample]):
            raise ValueError("ABX_sample and test_sample must have the same size as baseline_sample.")
        if not (isinstance(baseline_cohort, np.ndarray) and baseline_cohort.ndim == 2 and
                baseline_cohort.shape[1] == len(baseline_sample)):
            raise TypeError("baseline_cohort must be a 2D numpy array with the same number of columns as"
                            "baseline_sample, ABX_sample, and test_sample.")
        try:
            if np.isclose(np.sum(baseline_sample), 1.0):
                self.baseline_sample = baseline_sample
            else:
                self.baseline_sample = self._normalize_cohort(baseline_sample)
        except Exception as e:
            print(f"An error occurred: {e}")
        try:
            if np.isclose(np.sum(ABX_sample), 1.0):
                self.ABX_sample = ABX_sample
            else:
                self.ABX_sample = self._normalize_cohort(ABX_sample)
        except Exception as e:
            print(f"An error occurred: {e}")
        try:
            if np.isclose(np.sum(test_sample), 1.0):
                self.test_sample = test_sample
            else:
                self.test_sample = self._normalize_cohort(test_sample)
        except Exception as e:
            print(f"An error occurred: {e}")
        try:
            if np.isclose(np.sum(baseline_cohort.sum(axis=1)), baseline_cohort.shape[0]):
                self.baseline_cohort = baseline_cohort
            else:
                self.baseline_cohort = self._normalize_cohort(baseline_cohort)
        except Exception as e:
            print(f"An error occurred: {e}")
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
        nonzero_ABX = np.nonzero(self.ABX_sample)  # find nonzero elements of ABX sample
        nonzero_base = np.nonzero(self.baseline_sample)  # find nonzero elements of baseline sample
        nonzero_test = np.nonzero(self.test_sample)  # find nonzero elements of the test sample
        intersect_ABX_base = np.intersect1d(nonzero_ABX, nonzero_base)  # find species present in ABX and baseline
        intersect_ABX_test = np.intersect1d(nonzero_ABX, nonzero_test)  # find species present in ABX and test
        ARS = np.intersect1d(intersect_ABX_base, intersect_ABX_test)  # find ARS species by definition
        non_ARS = np.setdiff1d(np.arange(np.size(self.baseline_sample)), ARS)  # find non-ARS species
        column_sums = self.baseline_cohort.sum(axis=0)  # sum of each column of the baseline cohort
        zero_sum_column_indexes = np.where(column_sums == 0)[0]  # find species that are not present in the baseline
        # remove species that are not present in the baseline cohort from the non-ARS species
        non_ARS = np.setdiff1d(non_ARS, zero_sum_column_indexes)
        return ARS, non_ARS

    #def _find_ARS(self):
    #    """
    #    Find the ARS and non-ARS species
    #    :return: Indexes of the ARS species, i.e, the species that present in the baseline, ABX and future samples.
    #     Indexes of the non-ARS species, i.e, the species that are not ARS species.
    #    """
    #    nonzero_ABX = np.nonzero(self.ABX_sample)  # find nonzero elements of ABX sample
    #    nonzero_base = np.nonzero(self.baseline_sample)  # find nonzero elements of baseline sample
    #    ARS = np.intersect1d(nonzero_ABX, nonzero_base)  # find species present in ABX and baseline
    #    non_ARS = np.setdiff1d(np.arange(np.size(self.baseline_sample)), ARS)  # find non-ARS species
    #    column_sums = self.baseline_cohort.sum(axis=0)  # sum of each column of the baseline cohort
    #    zero_sum_column_indexes = np.where(column_sums == 0)[0]  # find species that are not present in the baseline
    #    # remove species that are not present in the baseline cohort from the non-ARS species
    #    non_ARS = np.setdiff1d(non_ARS, zero_sum_column_indexes)
    #    return ARS, non_ARS

    #def _create_shuffled_samples(self):
    #    return create_shuffled_samples(self.num_reals, self.test_sample.astype(np.float64), self.non_ARS.astype(np.int64),
    #                                   self.baseline_cohort.astype(np.float64), self.ARS.astype(np.int64))

    def _create_shuffled_samples(self):
        # initialize shuffled samples
        shuffled_samples = np.zeros((self.num_reals, np.size(self.test_sample)))
        #########################
        shuffled_samples[:, self.ARS] = self.test_sample[self.ARS]
        #shuffled_samples[:, self.ARS] = self.baseline_sample[self.ARS]
        #########################

        # set the constraint on the number of species
        stop = np.size(np.nonzero(self.test_sample[self.non_ARS]))

        # define the pool of species to be shuffled not including the ARS species
        pool = self.baseline_cohort[self.baseline_cohort != 0]
        pool_indices = np.where(self.baseline_cohort != 0)[1]
        index_to_remove = [i for i, ind in enumerate(pool_indices) if ind in self.ARS]
        pool_indices = np.array([ind for ind in pool_indices if ind not in self.ARS])
        pool = np.array([pool[i] for i in range(len(pool)) if i not in index_to_remove])

        for shuffled_sample in shuffled_samples:
            pool_copy = pool.copy()
            pool_indices_copy = pool_indices.copy()
            counter = 0
            while counter < stop:
                index = random.randint(0, len(pool_indices_copy) - 1)
                shuffled_sample[pool_indices_copy[index]] = pool_copy[index]
                similar_specie = np.where(pool_indices_copy == pool_indices_copy[index])[0]
                pool_copy = np.delete(pool_copy, similar_specie)
                pool_indices_copy = np.delete(pool_indices_copy, similar_specie)
                counter += 1

        # normalize the shuffled_samples
        shuffled_samples = self._normalize_cohort(shuffled_samples)
        return shuffled_samples

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized

    #def distance(self, method="Bray Curtis"):
    #    """
    #    :param method: distance method, either "Bray Curtis" or "rJSD".
    #    :return: distance between the test sample and the baseline sample and the distances between the shuffled samples
    #             and the baseline sample.
    #    """
    #    if method not in ["Bray Curtis", "rJSD"]:
    #        raise ValueError("method must be either 'Bray Curtis' or 'rJSD'.")
    #    if method == "Bray Curtis":
    #        bc_real = braycurtis(self.baseline_sample, self.test_sample)
    #        bc_shuffled = np.zeros(self.num_reals)
    #        for i in range(self.num_reals):
    #            bc_shuffled[i] = braycurtis(self.baseline_sample, self.shuffled_samples[i, :])
    #        return bc_real, bc_shuffled
    #    elif method == "rJSD":
    #        bc_real = jensenshannon(self.baseline_sample, self.test_sample)
    #        bc_shuffled = np.zeros(self.num_reals)
    #        for i in range(self.num_reals):
    #            bc_shuffled[i] = jensenshannon(self.baseline_sample, self.shuffled_samples[i, :])
    #        return bc_real, bc_shuffled

    def distance(self, method="Bray Curtis", otu_ids=None, tree=None):
        """
        :param method: distance method, either "Bray Curtis" or "rJSD".
        :return: distance between the test sample and the baseline sample and the distances between the shuffled samples
                 and the baseline sample.
        """
        if method not in ["Bray Curtis", "rJSD", "Jaccard", "Weighted unifrac", "Unweighted unifrac",
                          "Szymkiewicz Simpson second", "Szymkiewicz Simpson first"]:
            raise ValueError("method must be one of the following: Bray Curtis, rJSD, Jaccard, Weighted unifrac, "
                             "Unweighted unifrac, Szymkiewicz Simpson first, and Szymkiewicz Simpson second.")
        if method == "Bray Curtis":
            #dist_real = braycurtis(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
            #                       self._normalize_cohort(self.test_sample[self.non_ARS]))
            dist_real = braycurtis(self.baseline_sample, self.test_sample)#
            dist_shuffled = np.zeros(self.num_reals)
            #for i in range(self.num_reals):
                #dist_shuffled[i] = braycurtis(self.baseline_sample[self.non_ARS],
                #                            self.shuffled_samples[:, self.non_ARS][i, :])
            #    dist_shuffled[i] = braycurtis(self.baseline_sample,
            #                                  self.shuffled_samples[i, :])
            #dist_shuffled = np.apply_along_axis(lambda x: braycurtis(
            #    self._normalize_cohort(self.baseline_sample[self.non_ARS]), x[self.non_ARS]), 1,
            #                           self.shuffled_samples)
            dist_shuffled = np.apply_along_axis(lambda x: braycurtis(self.baseline_sample, x), 1,
                                                self.shuffled_samples)
            return dist_real, dist_shuffled
            return dist_real, dist_shuffled
        elif method == "rJSD":
            #dist_real = jensenshannon(self._normalize_cohort(self.baseline_sample[self.non_ARS]),
            #                        self._normalize_cohort(self.test_sample[self.non_ARS]))
            dist_real = jensenshannon(self.baseline_sample, self.test_sample)
            #dist_shuffled = np.apply_along_axis(lambda x: jensenshannon(
            #    self._normalize_cohort(self.baseline_sample[self.non_ARS]), x[self.non_ARS]), 1,
            #                                    self.shuffled_samples)
            dist_shuffled = np.apply_along_axis(lambda x: jensenshannon(self.baseline_sample, x), 1,
                                                self.shuffled_samples)
            return dist_real, dist_shuffled
        elif method == "Jaccard":
            #dist_real = 1 - Overlap(self.baseline_sample[self.non_ARS], self.test_sample[self.non_ARS],
            #                        overlap_type="Jaccard").calculate_overlap()
            dist_real = 1 - Overlap(self.baseline_sample, self.test_sample, overlap_type="Jaccard").calculate_overlap()
            #dist_shuffled = 1 - np.apply_along_axis(lambda x: Overlap(self.baseline_sample[self.non_ARS], x[self.non_ARS],
            #                                                          overlap_type="Jaccard").calculate_overlap(),
            #                                                          1, self.shuffled_samples)
            dist_shuffled = 1 - np.apply_along_axis(lambda x: Overlap(self.baseline_sample, x,
                                                                      overlap_type="Jaccard").calculate_overlap(),
                                                    1, self.shuffled_samples)
            return dist_real, dist_shuffled
        elif method == "Weighted unifrac":
            dist_real, dist_shuffled = self._unifrac(method, otu_ids, tree)
            return dist_real, dist_shuffled
        elif method == "Unweighted unifrac":
            dist_real, dist_shuffled = self._unifrac(method, otu_ids, tree)
            return dist_real, dist_shuffled
        elif method == "Szymkiewicz Simpson second":
            #dist_real = 1 - Overlap(self.baseline_sample[self.non_ARS], self.test_sample[self.non_ARS],
            #                        overlap_type="Szymkiewicz Simpson second").calculate_overlap()
            dist_real = 1 - Overlap(self.baseline_sample, self.test_sample,
                                    overlap_type="Szymkiewicz Simpson second").calculate_overlap()
            #dist_shuffled = 1 - np.apply_along_axis(
            #    lambda x: Overlap(self.baseline_sample[self.non_ARS], x[self.non_ARS],
            #                      overlap_type="Szymkiewicz Simpson second").calculate_overlap(),
            #    1, self.shuffled_samples)
            dist_shuffled = 1 - np.apply_along_axis(lambda x: Overlap(self.baseline_sample, x,
                                                                      overlap_type="Szymkiewicz Simpson second"
                                                                      ).calculate_overlap(),
                                                1, self.shuffled_samples)
            return dist_real, dist_shuffled
        elif method == "Szymkiewicz Simpson first":
            #dist_real = 1 - Overlap(self.baseline_sample[self.non_ARS], self.test_sample[self.non_ARS],
            #                        overlap_type="Szymkiewicz Simpson first").calculate_overlap()
            dist_real = 1 - Overlap(self.baseline_sample, self.test_sample,
                                    overlap_type="Szymkiewicz Simpson first").calculate_overlap()
            #dist_shuffled = 1 - np.apply_along_axis(
            #    lambda x: Overlap(self.baseline_sample[self.non_ARS], x[self.non_ARS],
            #                      overlap_type="Szymkiewicz Simpson first").calculate_overlap(),
            #    1, self.shuffled_samples)
            dist_shuffled = 1 - np.apply_along_axis(lambda x: Overlap(self.baseline_sample, x,
                                                                      overlap_type="Szymkiewicz Simpson first"
                                                                      ).calculate_overlap(),
                                                    1, self.shuffled_samples)
            return dist_real, dist_shuffled

    def _construct_filtered_data(self, first_sample, second_sample):
        data = np.vstack([first_sample, second_sample])
        filtered_data = self._normalize_cohort(data[:, self.non_ARS].squeeze())
        return filtered_data

    def _construct_pruned_tree(self, otu_ids, tree):
        remove_otu_ids = [otu_ids[i] for i in range(len(otu_ids)) if i in self.ARS]
        pruned_tree = tree.copy()
        for node in remove_otu_ids:
            to_delete = pruned_tree.find(node)
            pruned_tree.remove_deleted(lambda x: x == to_delete)
            pruned_tree.prune()
        return pruned_tree

    def _filter_otu_ids(self, otu_ids):
        return [otu_ids[i] for i in range(len(otu_ids)) if i in self.non_ARS]

    @staticmethod
    def _create_sample_ids(data):
        return np.arange(0, data.shape[0], 1).tolist()

    def _unifrac(self, method, otu_ids=None, tree=None):
        filtered_data = self._construct_filtered_data(self.baseline_sample, self.test_sample)
        sample_ids = self._create_sample_ids(filtered_data)
        filtered_otu_ids = self._filter_otu_ids(otu_ids)
        pruned_tree = self._construct_pruned_tree(otu_ids, tree)
        if method == "Weighted unifrac":
            dist_real = beta_diversity(metric='weighted_unifrac', counts=filtered_data,
                                       ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree, validate=False)[1, 0]
            dist_shuffled = []
            for smp in self.shuffled_samples:
                filtered_data = self._construct_filtered_data(self.baseline_sample, smp)
                dist_shuffled.append(beta_diversity(metric='weighted_unifrac', counts=filtered_data,
                                                    ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree,
                                                    validate=False)[1, 0])
            return dist_real, dist_shuffled
        elif method == "Unweighted unifrac":
            dist_real = beta_diversity(metric='unweighted_unifrac', counts=filtered_data,
                                       ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree, validate=False)[1, 0]
            dist_shuffled = []
            for smp in self.shuffled_samples:
                filtered_data = self._construct_filtered_data(self.baseline_sample, smp)
                dist_shuffled.append(beta_diversity(metric='unweighted_unifrac', counts=filtered_data,
                                                    ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree,
                                                    validate=False)[1, 0])
            return dist_real, dist_shuffled