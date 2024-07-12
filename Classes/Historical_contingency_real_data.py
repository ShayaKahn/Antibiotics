from Classes.Historical_contingency import HC
import numpy as np
from Classes.overlap import Overlap
from skbio.diversity import beta_diversity, alpha_diversity
from sklearn.metrics import pairwise_distances

class HcRealData():
    """
    This class is responsible for the simulation of the historical contingency model in real data.
    """
    def __init__(self, ABX_cohort, post_ABX_cohort, rep, tree=None, otu_ids=None, measure="Jaccard"):#, baseline_cohort=None):
        self.ABX_cohort = ABX_cohort
        self.post_ABX_cohort = post_ABX_cohort
        self.rep = rep
        self.tree = tree
        self.otu_ids = otu_ids
        self.shuffled_state = self._generate_shuffled_state()
        self.measure = measure
        if self.measure == "Jaccard":
            (self.similarity_perturbed, self.similarity_post_perturbed,
             self.similarity_shuffled_matrices) = self._calculate_Jaccard_similarity()
        elif self.measure == "UnWeightedUnifrac":
            (self.similarity_perturbed, self.similarity_post_perturbed,
             self.similarity_shuffled_matrices) = self._calculate_uwunifrac_similarity()
        #self.baseline_state = baseline_cohort
        #self.jaccard_baseline = self._calculate_Jaccard_base()

    #def calculate_ss_similarity(self):
    #    if self.baseline_state is not None:
    #        ss_vals = []
    #        for i in range(self.perturbed_state.shape[0] - 1):
    #            for j in range(i + 1, self.perturbed_state.shape[0]):
    #                # Find the indexes that are not zero in both samples.
    #                nonzero_index_first = np.nonzero(self.perturbed_state[i, :])
    #                nonzero_index_second = np.nonzero(self.perturbed_state[j, :])
    #                intersect = np.intersect1d(nonzero_index_first, nonzero_index_second)
    #                # Find the complement of the intersection.
    #                full_array = np.arange(0, self.perturbed_state.shape[1])
    #                diff = np.setdiff1d(full_array, full_array[intersect])
    #                # Calculate zymkiewicz Simpson similarity.
    #                overlap_object = Overlap(self.filtered_post_perturbed_state[i, :][diff],
    #                                         self.baseline_state[i, :][diff],
    #                                         overlap_type="Szymkiewicz Simpson")
    #                ss_vals.append(overlap_object.calculate_overlap())
    #                overlap_object = Overlap(self.filtered_post_perturbed_state[j, :][diff],
    #                                         self.baseline_state[j, :][diff],
    #                                         overlap_type="Szymkiewicz Simpson")
    #                ss_vals.append(overlap_object.calculate_overlap())
    #        return ss_vals
    #    else:
    #        return None

    #def _calculate_Jaccard_base(self):
    #    if self.baseline_state is not None:
    #        jaccard_baseline = []
    #        for i in range(self.baseline_state.shape[0] - 1):
    #            for j in range(i + 1, self.baseline_state.shape[0]):
    #                overlap_object_perturbed = Overlap(self.perturbed_state[i, :],
    #                                                   self.perturbed_state[j, :], overlap_type="Jaccard")
    #                overlap_object_perturbed.calculate_overlap()
    #                intrsection_ind = overlap_object_perturbed.intersection_bool
    #                mask = np.ones(len(intrsection_ind), dtype=bool)
    #                mask[intrsection_ind] = False
    #                overlap_object_baseline = Overlap(self.baseline_state[i, :][mask],
    #                                                  self.baseline_state[j, :][mask], overlap_type="Jaccard")
    #                jaccard_baseline.append(overlap_object_baseline.calculate_overlap())
    #        return jaccard_baseline
    #    else:
    #        return None

    def _calculate_Jaccard_similarity(self):
        jaccard_perturbed = []
        jaccard_post_perturbed = []
        jaccard_shuffled_matrices = []
        for i in range(self.ABX_cohort.shape[0] - 1):
            for j in range(i + 1, self.ABX_cohort.shape[0]):
                overlap_object_perturbed = Overlap(self.ABX_cohort[i, :],
                                                   self.ABX_cohort[j, :],
                                                   overlap_type="Jaccard")

                jaccard_perturbed.append(overlap_object_perturbed.calculate_overlap())
                intrsection_ind = overlap_object_perturbed.intersection_bool

                mask = np.ones(len(intrsection_ind), dtype=bool)
                mask[intrsection_ind] = False
                overlap_object_post_perturbed = Overlap(self.post_ABX_cohort[i, :][mask],
                                                        self.post_ABX_cohort[j, :][mask],
                                                        overlap_type="Jaccard")
                jaccard_post_perturbed.append(overlap_object_post_perturbed.calculate_overlap())
                jaccard_shuffled_matrices.append(self._calculate_jaccard_similarity_matrix(self.shuffled_state, mask))
        return jaccard_perturbed, jaccard_post_perturbed, jaccard_shuffled_matrices

    def _generate_shuffled_state(self):
        shuffled_state = self.post_ABX_cohort.copy()
        extended_shuffled_state = np.tile(shuffled_state, (self.rep, 1))
        np.apply_along_axis(np.random.shuffle, 0, extended_shuffled_state)
        return extended_shuffled_state

    def _calculate_uwunifrac_similarity(self):
        unifrac_perturbed = []
        unifrac_post_perturbed = []
        unifrac_shuffled_matrices = []
        for i in range(self.ABX_cohort.shape[0] - 1):
            print(i)
            for j in range(i + 1, self.ABX_cohort.shape[0]):
                data = np.vstack([self.ABX_cohort[i, :], self.ABX_cohort[j, :]])
                sample_ids = np.arange(0, data.shape[0], 1).tolist()
                unj = beta_diversity(metric='unweighted_unifrac', counts=data,
                                     ids=sample_ids, taxa=self.otu_ids, tree=self.tree, validate=False)[1, 0]
                ind = np.where((self.ABX_cohort[i, :] == 0) | (self.ABX_cohort[j, :] == 0))
                unifrac_perturbed.append(1 - unj)
                data_post = np.vstack([self.post_ABX_cohort[i, :],
                                       self.post_ABX_cohort[j, :]])
                filtered_data_post = data_post[:, ind].squeeze()
                filtered_otu_ids = [self.otu_ids[i] for i in range(len(self.otu_ids)) if i in ind[0]]
                remove_otu_ids = [self.otu_ids[i] for i in range(len(self.otu_ids)) if i not in ind[0]]
                pruned_tree = self.tree.copy()
                for node in remove_otu_ids:
                    to_delete = pruned_tree.find(node)
                    pruned_tree.remove_deleted(lambda x: x == to_delete)
                    pruned_tree.prune()
                unj_post = beta_diversity(metric='unweighted_unifrac', counts=filtered_data_post,
                                          ids=sample_ids, taxa=filtered_otu_ids, tree=pruned_tree, validate=False)[1, 0]
                unifrac_post_perturbed.append(1 - unj_post)
                shuffled_unifrac_matrix = beta_diversity(metric='unweighted_unifrac',
                                                         counts=self.shuffled_state[:, ind].squeeze(),
                                                         ids=np.arange(0, self.shuffled_state.shape[0], 1).tolist(),
                                                         taxa=filtered_otu_ids,
                                                         tree=pruned_tree,
                                                         validate=True)
                shuffled_unifrac_matrix = shuffled_unifrac_matrix.to_data_frame()
                shuffled_unifrac_matrix = shuffled_unifrac_matrix.values
                unifrac_shuffled_matrices.append(1 - shuffled_unifrac_matrix)
        return unifrac_perturbed, unifrac_post_perturbed, unifrac_shuffled_matrices

    @staticmethod
    def _calculate_jaccard_similarity_matrix(cohort, mask):
        binary_matrix = cohort[:, mask].copy()
        ind = np.where(binary_matrix != 0)
        binary_matrix[ind] = 1
        jaccard_similarity_mat = 1 - pairwise_distances(binary_matrix.astype(bool), metric='jaccard')
        return jaccard_similarity_mat

    @staticmethod
    def _normalize_cohort(cohort):
        # normalization function
        if cohort.ndim == 1:
            cohort_normalized = cohort / cohort.sum()
        else:
            cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
        return cohort_normalized

    def get_results(self):
        results = {"similarity perturbed": self.similarity_perturbed,
                   "similarity post perturbed": self.similarity_post_perturbed,
                   #"jacard baseline": self.jaccard_baseline,
                   "similarity shuffled matrices": self.similarity_shuffled_matrices}
        return results