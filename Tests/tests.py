from unittest import TestCase
from Data_manipulation.optimal import OptimalCohort
from Classes.overlap import Overlap
from Classes.null_model import NullModel
from Data_manipulation.rarify import Rarify
import numpy as np
import pandas as pd
from Classes.Historical_contingency import HC
from Classes.shuffle_post_ABX import ShufflePostABX
from Classes.surrogate import Surrogate
from Classes.Stochastic_surrogate import StochasticSurrogate

class TestStochasticSurrogate(TestCase):
    """
    This class tests the StochasticSurrogate class.
    """
    def setUp(self) -> None:
        self.test_base_sample = {"Test subject": np.array([0.1, 0, 0.2, 0.4, 0, 0, 0.1, 0.2])}
        self.base_samples_collections = {"Subject A": np.array([0, 0.1, 0.3, 0.3, 0, 0.2, 0.1, 0]),
                                         "Subject B": np.array([0, 0.3, 0, 0.4, 0.1, 0, 0.1, 0.1])}
        self.test_ABX_sample = np.array([0.4, 0, 0, 0.6, 0, 0, 0, 0])
        self.test_post_ABX_sample = np.array([0.2, 0, 0.3, 0.4, 0, 0, 0.1, 0])
        self.reals = 3
        self.weighted = False
        self.stochastic = StochasticSurrogate(self.test_base_sample, self.base_samples_collections,
                                              self.test_post_ABX_sample, self.test_ABX_sample, self.reals,
                                              self.weighted)
        self.stochastic_samples = self.stochastic.stochastic_base_samples_collections
        self.results = self.stochastic.apply_surrogate_data_analysis(method="Szymkiewicz Simpson first")

    def test_find_ARS(self):
        self.assertListEqual(list(self.stochastic.non_ARS), [1, 2, 4, 5, 6, 7])
        self.assertListEqual(list(self.stochastic.ARS), [0, 3])

    def test_generate_stochastic_surrogate(self):
        for key in self.stochastic_samples.keys():
            nonzero_vals = [np.size(np.nonzero(smp)) for smp in self.stochastic_samples[key]]
            self.assertTrue(all(x == nonzero_vals[0] for x in nonzero_vals))
            self.assertEqual(self.stochastic_samples[key].shape[0], self.reals)
            self.assertEqual(self.stochastic_samples[key].shape[1], len(self.test_base_sample["Test subject"]))

class TestHC(TestCase):
    """
    This class tests the HC class.
    """
    def setUp(self) -> None:
        self.num_samples = 3
        self.pool_size = 40
        self.num_survived_min = 20
        self.num_survived_max = 35
        self.mean = 0
        self.sigma = 5
        self.c = 0.05
        self.delta = 1e-5
        self.final_time = 1000
        self.max_step = 0.05
        self.epsilon = 1e-3
        self.threshold = 1e-3
        self.min_growth = 0
        self.max_growth = 1
        self.rep = 10
        self.symmetric = True
        self.alpha = 1  # None
        self.method = 'RK45'
        self.multiprocess = False
        self.weighted = False

        self.HC_object = HC(self.num_samples, self.pool_size, self.num_survived_min, self.num_survived_max, self.mean,
                            self.sigma, self.c, self.delta, self.final_time, self.max_step, self.epsilon,
                            self.threshold, self.min_growth, self.max_growth, self.rep, self.symmetric, self.alpha, self.method,
                            self.multiprocess, self.weighted)
        self.results = self.HC_object.get_results()

    def test_set_initial_conditions(self):
        # Test if the number of survived species is correct.
        for index, smp in enumerate(self.results["perturbed state"]):
            self.assertEqual(len(self.results["survived matrix"][index]), np.size(np.nonzero(smp)))

    def test_insert_total_pool(self):
        # Test if the procedure of inserting the total pool is correct.
        for smp, ind in zip(self.results["new perturbed state"], self.results["survived matrix"]):
            mask = np.ones(smp.shape[0], dtype=bool)
            mask[ind] = False
            self.assertTrue(np.all(smp[mask] == smp[mask][0]))
            # check if all values are unique
            self.assertTrue(len(smp[~mask]) == len(set(smp[~mask])), "Not all values are unique")

    def test_generate_shuffled_state(self):
        self.assertTrue(np.all(np.isclose(self.results["shuffled state"].sum(axis=0),
                        np.multiply(self.results["filtered post perturbed state"].sum(axis=0),
                                    self.rep), atol=1e-6)))
        self.assertFalse(np.all(np.isclose(self.results["shuffled state"][0, :],
                                           self.results["filtered post perturbed state"][0, :],
                                           atol=1e-6)))

    def test_calculate_jaccard_similarity_matrix(self):
        cohort = np.array([[0.2, 0.4, 0, 0.3, 0.1], [0, 0.6, 0, 0.2, 0.2], [0.2, 0, 0.4, 0.3, 0.1]])
        mask = np.array([True, True, True, False, False])
        mat = self.HC_object._calculate_jaccard_similarity_matrix(cohort, mask)
        rows, cols = np.tril_indices(mat.shape[0], k=-1)
        jaccard_vactor_first_method = mat[rows, cols]
        jaccard_vactor_second_method = []
        for i in range(cohort.shape[0] - 1):
            for j in range(i + 1, cohort.shape[0]):
                overlap_object = Overlap(cohort[i, :][mask],
                                         cohort[j, :][mask],
                                         overlap_type="Jaccard")

                jaccard_vactor_second_method.append(overlap_object.calculate_overlap())
        jaccard_vactor_second_method = np.array(jaccard_vactor_second_method)
        self.assertTrue(np.all(np.isclose(jaccard_vactor_first_method, jaccard_vactor_second_method,
                                          atol=1e-6)))

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
        self.ABX_sample = np.array([0, 0, 0.3, 0.1, 0, 0.1, 0, 0, 0, 0.5])
        self.test_sample = np.array([0.5, 0, 0.1, 0, 0.2, 0.1, 0, 0.1, 0, 0])
        self.baseline = np.array([[0.1, 0, 0, 0.2, 0.3, 0.1, 0, 0.1, 0.2, 0],
                                  [0.5, 0.1, 0, 0, 0.1, 0.1, 0, 0.1, 0.1, 0],
                                  [0.1, 0, 0.2, 0.1, 0, 0.1, 0, 0, 0.2, 0.3]])
        self.num_reals = 3
        self.null_model = NullModel(self.baseline_sample, self.ABX_sample, self.test_sample,
                                    self.baseline, self.num_reals)
        print(f' ARS: {self.null_model.ARS}')
        print(f' non-ARS: {self.null_model.non_ARS}')
        print(self.null_model.shuffled_samples.sum(axis=1))

    def test_distance(self):
        self.assertEqual(np.size(self.null_model.ARS) + np.size(self.null_model.non_ARS),
                         np.size(np.nonzero(np.sum(self.baseline, axis=0))))
        self.assertEqual(np.size(np.nonzero(self.null_model.shuffled_samples[0,
                         self.null_model.non_ARS])),
                         np.size(np.nonzero(self.test_sample[self.null_model.non_ARS])))

class Test_Surrogate(TestCase):
    def setUp(self) -> None:
        self.test_samples = np.array([[0.1, 0, 0.2, 0.4, 0, 0, 0.1, 0.2],
                                      [0, 0, 0.5, 0.5, 0, 0, 0, 0],
                                      [0.1, 0, 0.4, 0.4, 0.1, 0, 0, 0],
                                      [0.1, 0, 0.3, 0.3, 0.1, 0, 0, 0.2]])
        self.test_base_samples_collection = {"Test subject": np.array([[0.1, 0, 0.2, 0.4, 0, 0, 0.1, 0.2]])}
        self.base_samples_collections = {"Subject A": np.array([[0, 0.1, 0.3, 0.3, 0, 0.2, 0.1, 0]]),
                                         "Subject B": np.array([[0.2, 0.1, 0, 0.5, 0.1, 0, 0, 0.1]])}
        self.upper = 2
        self.surrogate = Surrogate(self.test_base_samples_collection, self.base_samples_collections,
                                   self.test_samples, self.upper)
        self.results = self.surrogate.apply_surrogate_data_analysis(method="Szymkiewicz Simpson first")
        self.results_sec = self.surrogate.apply_surrogate_data_analysis(method="Szymkiewicz Simpson second")
        filtered_data = self.surrogate._construct_filtered_data(self.test_samples[0], self.test_samples[2])

    def test_find_non_ARS(self):
        self.assertListEqual(list(self.surrogate.non_ARS), [0, 1, 4, 5, 6, 7])
        self.assertListEqual(list(self.surrogate.ARS), [2, 3])

    def test_apply_surrogate_data_analysis(self):
        self.assertListEqual(self.results["Test subject"], [1, 0, 1/3, 2/3])
        self.assertListEqual(self.results_sec["Test subject"], [1, 0, 1 / 2, 2 / 3])

    def test_construct_filtered_data(self):
        filtered_data = self.surrogate._construct_filtered_data(self.test_samples[0], self.test_samples[2])
        print(filtered_data.shape)
        self.assertTrue(np.all(np.isclose(filtered_data, np.array([[0.25, 0, 0, 0, 0.25, 0.5], [0.5, 0, 0.5, 0, 0, 0]]),
                                          atol=1e-6)))

    def test_construct_pruned_tree(self):
        pass