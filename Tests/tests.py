from unittest import TestCase
from optimal import OptimalCohort
from overlap import Overlap
from null_model import NullModel
from rarify import Rarify
import numpy as np
import pandas as pd
from Historical_contingency import HC
from shuffle_post_ABX import ShufflePostABX

class TestHC(TestCase):
    """
    This class tests the HC class.
    """
    def setUp(self) -> None:
        self.num_samples = 20
        self.pool_size = 20
        self.num_survived = 10
        self.mean = 0.0
        self.sigma = 0.1
        self.lower = 0.1
        self.higher = 0.7
        self.delta = 1e-2
        self.final_time = 1e5
        self.max_step = 0.1
        self.epsilon = 1e-3
        self.threshold = 1e-3
        self.min_skew = -15
        self.max_skew = 0
        self.interaction_strength = 2
        self.method = 'RK45'
        self.measure = 'Jaccard'
        self.HC_object = HC(self.num_samples, self.pool_size, self.num_survived, self.mean, self.sigma, self.lower,
                            self.higher, self.delta, self.final_time, self.max_step, self.epsilon, self.threshold,
                            self.min_skew, self.max_skew, self.interaction_strength, self.method, self.measure)
        self.results = self.HC_object.get_results()

    def test_set_initial_conditions(self):
        for smp in self.HC_object.y0:
            self.assertEqual(self.num_survived, np.size(np.nonzero(smp)))

    def test_insert_total_pool(self):
        for smp, ind in zip(self.HC_object.new_perturbed_state, self.HC_object.survived_matrix):
            mask = np.ones(smp.shape[0], dtype=bool)
            mask[ind] = False
            self.assertTrue(np.all(smp[mask] == smp[mask][0]))

class TestShufflePostABX(TestCase):
    """
    This class tests the ShufflePostABX class.
    """
    def setUp(self) -> None:
        self.ABX_cohort = np.array([[0, 0, 0, 1, 0, 0, 0, 0],
                                    [0.5, 0, 0, 0, 0, 0.5, 0, 0],
                                    [0.4, 0, 0.2, 0.4, 0, 0, 0, 0]])
        self.post_ABX_cohort = np.array([[0, 0, 0.3, 0.4, 0, 0, 0, 0.3],
                                        [0.2, 0, 0, 0.3, 0.2, 0, 0, 0.2],
                                        [0.1, 0, 0.2, 0.2, 0.1, 0, 0, 0.4]])
        self.shuffle = ShufflePostABX(self.ABX_cohort, self.post_ABX_cohort)

    def test_find_survived_matrix(self):
        self.assertEqual(len(self.shuffle.survived_matrix), self.ABX_cohort.shape[0])


class TestOptimalCohort(TestCase):
    """
    This class tests the OptimalCohort class.
    """
    def setUp(self) -> None:
        # Two default samples.
        self.samples_dict = {'a': np.array([[11, 0, 8], [3, 9, 2], [0, 1, 3]]),
                             'b': np.array([[7, 1, 2], [1, 6, 0], [2, 3, 8], [8, 2, 5], [0, 1, 0]]),
                             'c': np.array([[35, 0, 17], [3, 4, 3], [1, 0, 8]]),
                             'd': np.array([[12, 7, 4], [1, 0, 0], [7, 1, 0], [6, 6, 6]])}
        self.optimal = OptimalCohort(self.samples_dict)
        self.optimal_samples = self.optimal.get_optimal_samples()
    def test_get_optimal_samples(self):
        self.assertEqual(np.sum(self.optimal_samples[0], axis=1).tolist(),
                         np.ones(np.size(self.optimal_samples[0], axis=0)).tolist())

class Test_Rarify(TestCase):

    def setUp(self) -> None:
        self.df = pd.DataFrame({
                  'A': [1, 3, 0, 2],
                  'B': [2, 0, 2, 1],
                  'C': [0, 4, 0, 0],
                  'D': [50, 12, 0, 0],
                   })
        self.object_min = Rarify(self.df)
        self.object_depth = Rarify(self.df, depth=5)

    def test_rarify(self):
        rar_df_min = self.object_min.rarify()
        rar_df_depth = self.object_depth.rarify()
        print(rar_df_min)
        print(rar_df_depth)

class Test_Overlap(TestCase):
    def setUp(self) -> None:
        self.first_sample = np.array([0.1, 0, 0.2, 0.4, 0, 0, 0.1, 0.2])
        self.second_sample = np.array([0.1, 0, 0.2, 0.4, 0.1, 0, 0, 0.2])

    def test_calculate_overlap(self):
        overlap = Overlap(self.first_sample, self.second_sample,
                          overlap_type='Jaccard').calculate_overlap()
        print(overlap)

class Test_NullModel(TestCase):
    def setUp(self) -> None:
        self.baseline_sample = np.array([0.1, 0, 0, 0.2, 0.3, 0.1, 0, 0.1, 0.2, 0])
        self.test_sample = np.array([0.5, 0, 0.1, 0, 0.2, 0.1, 0, 0.1, 0, 0])
        self.ABX_sample = np.array([0, 0, 0.3, 0.1, 0, 0.1, 0, 0, 0, 0.5])
        self.baseline = np.array([[0.1, 0, 0, 0.2, 0.3, 0.1, 0, 0.1, 0.2, 0],
                                  [0.5, 0.1, 0, 0, 0.1, 0.1, 0, 0.1, 0.1, 0],
                                  [0.1, 0, 0.2, 0.1, 0, 0.1, 0, 0, 0.2, 0.3]])
        self.num_reals = 3
        self.null_model = NullModel(self.baseline_sample, self.ABX_sample, self.test_sample,
                                    self.baseline, self.num_reals)
        print(self.null_model.ARS)
        print(self.null_model.non_ARS)
        print(self.null_model.shuffled_samples)

    def test_distance(self):
        self.assertEqual(np.size(self.null_model.ARS) + np.size(self.null_model.non_ARS),
                         np.size(np.nonzero(np.sum(self.baseline, axis=0))))
        self.assertEqual(np.size(np.nonzero(self.null_model.shuffled_samples[0,
                         self.null_model.non_ARS])),
                         np.size(np.nonzero(self.test_sample[self.null_model.non_ARS])))