from scipy.spatial.distance import braycurtis
from scipy.spatial.distance import pdist, squareform
import numpy as np
class OptimalCohort:
    """
    This class calculates the IDOA values vector for a cohort or the IDOA value for a sample
    with respect to a reference cohort.
    """
    def __init__(self, samples_dict, criterion='lower'):
        """
        :param samples_dict: dictionary that contains n subjects as a keys and a matrix of baseline samples, the rows of
        the matrices represent the samples and the columns represent the species.
        """
        self.samples_dict = samples_dict
        self.dissimilarity_matrix_dict = {}
        self.criterion = criterion

        self._normalize_data()

        for key in self.samples_dict:
            y = pdist(self.samples_dict[key], metric='braycurtis')
            dissimilarity_matrix = squareform(y)
            self.dissimilarity_matrix_dict[key] = dissimilarity_matrix

    def _get_optimal_index(self):

        ind_container = {}

        if self.criterion == 'lower':

            for key in self.dissimilarity_matrix_dict:
                # Get the indices of the lower triangle of the matrix
                ind = np.tril_indices(self.dissimilarity_matrix_dict[key].shape[0], k=-1, m=None)
                # Find the index of the smallest value in the lower triangle
                index = np.argmin(self.dissimilarity_matrix_dict[key][ind])
                # Convert the flattened index to row and column indices of the original matrix
                row_index, col_index = ind[0][index], ind[1][index]
                indices = (row_index, col_index)
                ind_container[key] = indices
            return ind_container

        elif self.criterion == 'mean':

            for key in self.dissimilarity_matrix_dict:
                # calculate the mean distance for each sample
                mean_dist = np.mean(self.dissimilarity_matrix_dict[key], axis=1)
                # get the index of the sample with the minimum mean distance
                min_idx = np.argmin(mean_dist)
                ind_container[key] = min_idx
            return ind_container

    def get_optimal_samples(self):

        ind_container = self._get_optimal_index()
        optimal_samples = []
        chosen_indices = {}

        if self.criterion == 'lower':

            for key in self.samples_dict:
                # Calculate the mean dissimilarity of each row to the other rows
                mean_dissimilarity_row0 = np.mean(self.dissimilarity_matrix_dict[key][ind_container[key][0], :])
                mean_dissimilarity_row1 = np.mean(self.dissimilarity_matrix_dict[key][ind_container[key][1], :])
                if mean_dissimilarity_row0 <= mean_dissimilarity_row1:
                    optimal_samples.append(self.samples_dict[key][ind_container[key][0], :])
                    chosen_indices[key] = ind_container[key][0]
                else:
                    optimal_samples.append(self.samples_dict[key][ind_container[key][1], :])
                    chosen_indices[key] = ind_container[key][1]
            optimal_samples = np.vstack(optimal_samples)
            return optimal_samples, chosen_indices

        elif self.criterion == 'mean':

            for key in self.samples_dict:
                optimal_samples.append(self.samples_dict[key][ind_container[key], :])
                chosen_indices[key] = ind_container[key]
            optimal_samples = np.vstack(optimal_samples)
            return optimal_samples, chosen_indices

    def _normalize_data(self):
        for key in self.samples_dict:
            # Normalize the samples to sum to 1
            self.samples_dict[key] = self.samples_dict[key] / \
                                     self.samples_dict[key].sum(axis=1, keepdims=True)