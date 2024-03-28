import numpy as np
from Historical_contingency import HC


class ShufflePostABX(HC):
    """
    This class is responsible for shuffling the post-antibiotics cohort
     while keeping the surviving species present.
    """
    def __init__(self, ABX_cohort, post_ABX_cohort):
        self.ABX_cohort = ABX_cohort
        self.filtered_post_perturbed_state = post_ABX_cohort
        self.survived_matrix = self._find_survived_matrix()
        self.shuffled_state = self._generate_shuffled_state()

    def _find_survived_matrix(self):
        survived_matrix = []
        for i, (abx, post) in enumerate(zip(self.ABX_cohort,
                                            self.filtered_post_perturbed_state)):
            survived_matrix.append(np.where((abx != 0) & (post != 0))[0])
        return survived_matrix



