from unittest import TestCase
from optimal import OptimalCohort
from overlap import Overlap
from null_model import NullModel
from rarify import Rarify
import numpy as np
import pandas as pd

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
        self.null_model_include = NullModel(self.baseline_sample, self.ABX_sample, self.test_sample,
                                            self.baseline, self.num_reals, exclude_test=False)
        self.null_model_exclude = NullModel(self.baseline_sample, self.ABX_sample, self.test_sample,
                                            self.baseline, self.num_reals, exclude_test=True)

    def test_distance(self):
        self.assertEqual(np.size(self.null_model_include.ARS) + np.size(self.null_model_include.non_ARS),
                         np.size(np.nonzero(np.sum(self.baseline, axis=0))))
        self.assertEqual(np.size(self.null_model_exclude.ARS) + np.size(self.null_model_exclude.non_ARS),
                         np.size(np.nonzero(np.sum(self.baseline, axis=0))))
        self.assertEqual(np.size(np.nonzero(self.null_model_include.shuffled_samples[0,
                         self.null_model_include.non_ARS])),
                         np.size(np.nonzero(self.test_sample[self.null_model_exclude.non_ARS])))
        self.assertEqual(np.size(np.nonzero(self.null_model_exclude.shuffled_samples[0,
                         self.null_model_exclude.non_ARS])),
                         np.size(np.nonzero(self.test_sample[self.null_model_exclude.non_ARS])))