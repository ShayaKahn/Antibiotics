from Classes.overlap import Overlap
from skbio.diversity import beta_diversity
from skbio.tree._tree import TreeNode
import numpy as np

class Surrogate:
    """
    This class implements surrogate data analysis
    """
    def __init__(self, test_base_samples_collection, base_samples_collections, test_samples, upper, lower=0):
        """
        :param test_base_samples_collection: dictionary that contain a numpy matrix of shape (# baseline_samples,
         # species) that represents the OTU table of the baseline samples of a test subject. The key is subjects
         identifier.
        :param base_samples_collections: dictionary of numpy matrices of shape (# baseline_samples, # species) that
         represent the baseline collections of different subjects in the experiment. The keys are the identifiers of the
         subjects.
        :param test_samples: numpy matrix of shape (# samples, # species) that contains the time series samples of the
         test subject
        :param upper: the upper threshold of calculating the ARS species
        """
        # separate the key and the matrix from the dictionary
        if not (isinstance(test_base_samples_collection, dict)):
            raise ValueError("test_base_samples_collection must be dictionary")
        self.test_key, self.test_base_samples_collection = next(iter(test_base_samples_collection.items()))
        if not (isinstance(self.test_key, str) and isinstance(self.test_base_samples_collection, np.ndarray)):
            raise ValueError("test_base_samples_collection must be dictionary with string as a key and numpy matrix of"
                             " shape (# baseline_samples, # species) as value.")
        if not (isinstance(base_samples_collections, dict)):
            raise ValueError("base_samples_collections must be dictionary")
        self.base_samples_collections = base_samples_collections
        if not (isinstance(test_samples, np.ndarray)):
            raise ValueError("test_samples must be numpy matrix of shape (# samples, # species)")
        if not (test_samples.shape[1] == self.test_base_samples_collection.shape[1]):
            raise ValueError("The dimensions are not consistent")
        self.test_samples = test_samples
        if not (isinstance(upper, int) and (2 <= upper)):
            raise ValueError("upper must be integer greater than 1")
        if not (isinstance(lower, int) and (0 <= lower < upper)):
            raise ValueError("lower must be integer greater or equal to 0 and smaller than upper")
        self.upper = upper
        self.lower = lower
        # find the indexes of the ARS and the non-ARS species
        self.non_ARS, self.ARS = self._find_non_ARS()

    def _find_non_ARS(self):
        # find the indexes of the ARS species, which are the species that are present in the baseline and ABX samples
        cols_with_zeros = np.any(self.test_samples[self.lower:self.upper, :] == 0, axis=0)
        indexes = np.where(cols_with_zeros)[0]
        indexes_comp = np.where(~cols_with_zeros)[0]
        return indexes, indexes_comp

    def apply_surrogate_data_analysis(self, method="Szymkiewicz Simpson first", tree=None, otu_ids=None):
        if method == "unweighted_unifrac":
            if tree is None or otu_ids is None:
                raise ValueError("tree and otu_ids must be provided")
            if not (isinstance(tree, TreeNode) and isinstance(otu_ids, list)):
                raise ValueError("tree must be an instance of NullModel and otu_ids must be a list")
            # calculate results using unweighted unifrac
            results = self._unweighted_unifrac_surrogate(tree, otu_ids)
        else:
            # calculate results using overlap class
            results = {}
            if self.test_base_samples_collection.shape[0] != 1:
                 measures = [np.mean([Overlap(base[self.non_ARS], sample[self.non_ARS],
                                      overlap_type=method).calculate_overlap() for base in
                                      self.test_base_samples_collection if not np.array_equal(
                                      base, sample)]) for sample in self.test_samples]
            else:
                print(f"Subject {self.test_key}")
                print(np.size(np.nonzero(self.test_base_samples_collection.squeeze()[self.ARS])[0]))
                measures = [Overlap(self.test_base_samples_collection.squeeze()[self.non_ARS], sample[self.non_ARS],
                                    overlap_type=method).calculate_overlap() for sample in self.test_samples]
            # store the results of the test subject
            results[self.test_key] = measures

            samples_copy = self.test_samples.copy()
            samples_copy[:, self.ARS] = 0

            # calculate the results of the surrogate data
            for key in self.base_samples_collections.keys():
                if self.base_samples_collections[key].shape[0] != 1:
                    measures_surrogate = [np.mean([Overlap(base[self.non_ARS], sample[self.non_ARS],
                                                           overlap_type=method).calculate_overlap() for base in
                                                   self.base_samples_collections[key]]) for sample in self.test_samples]
                else:
                    print(np.size(np.nonzero(self.base_samples_collections[key].squeeze()[self.ARS])[0]))
                    measures_surrogate = [Overlap(self.base_samples_collections[key].squeeze()[self.non_ARS],
                                                  sample[self.non_ARS],
                                                  overlap_type=method).calculate_overlap(
                                            ) for sample in self.test_samples]
                results[key] = measures_surrogate

        return results

    def _construct_filtered_data(self, first_sample, second_sample):
        data = np.vstack([first_sample, second_sample])
        filtered_data = self._normalize_cohort(data[:, self.non_ARS].squeeze())
        return filtered_data

    def _construct_pruned_tree(self, otu_ids, tree):
        # This method constructs a pruned tree by filtering the ARS species
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
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized

    @staticmethod
    def _create_sample_ids(data):
        return np.arange(0, data.shape[0], 1).tolist()

    def _unweighted_unifrac_surrogate(self, tree, otu_ids):
        """
        This method implements the surrogate data analysis using unweighted unifrac
        """
        results = {}
        # build a pruned tree by filtering the ARS species
        pruned_tree = self._construct_pruned_tree(otu_ids, tree)
        # filter the ARS species from the OTU ids
        filtered_otu_ids = self._filter_otu_ids(otu_ids)
        measures = []
        for sample in self.test_samples:
            # filter the ARS species from the data
            filtered_data = self._construct_filtered_data(self.test_base_samples_collection.squeeze(), sample)
            sample_ids = self._create_sample_ids(filtered_data)
            measures.append(1 - beta_diversity(metric='unweighted_unifrac', counts=filtered_data,
                                               ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree,
                                               validate=False)[1, 0])
        # store the results for the test subject
        results[self.test_key] = measures
        for key in self.base_samples_collections.keys():
            measures_surrogate = []
            for sample in self.test_samples:
                filtered_data = self._construct_filtered_data(self.base_samples_collections[key].squeeze(), sample)
                sample_ids = self._create_sample_ids(filtered_data)
                measures_surrogate.append(1 - beta_diversity(metric='unweighted_unifrac', counts=filtered_data,
                                               ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree,
                                               validate=False)[1, 0])
                results[key] = measures_surrogate
        return results