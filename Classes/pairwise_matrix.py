from scipy.spatial.distance import pdist, squareform
from Classes.overlap import Overlap
import numpy as np

class Pairwise:
    def __init__(self, samples):
        self.samples = samples

    def create_pairwise_matrix(self, metric='braycurtis'):
        return squareform(pdist(self.samples, metric))

    def create_pairwise_matrix_similarity(self, metric='Overlap'):
        num_samples = len(self.samples)
        similarity_matrix = np.zeros((num_samples, num_samples))

        for i in range(num_samples):
            for j in range(i, num_samples):
                similarity = Overlap(self.samples[i], self.samples[j], metric).calculate_overlap()
                similarity_matrix[i, j] = similarity
                similarity_matrix[j, i] = similarity

        return similarity_matrix