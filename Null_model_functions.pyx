import numpy as np
cimport numpy as cnp
from libc.stdlib cimport rand
import random

def create_shuffled_samples(int num_reals, cnp.ndarray[cnp.float64_t, ndim=1] test_sample,
                            cnp.ndarray[cnp.int64_t, ndim=1]  non_ARS,
                            cnp.ndarray[cnp.float64_t, ndim=2] baseline_cohort,
                            cnp.ndarray[cnp.int64_t, ndim=1] ARS):
    # initialize shuffled samples
    cdef cnp.ndarray[cnp.float64_t, ndim=2] shuffled_samples = np.zeros((num_reals, np.shape(test_sample)[0]))

    # set the constraint on the number of species
    #cdef int stop = np.size(np.nonzero(test_sample[non_ARS]))

    cdef int stop = 0
    cdef int j

    for j in non_ARS:
        if test_sample[j] != 0:
            stop += 1

    cdef cnp.ndarray[cnp.float64_t, ndim=1] pool = baseline_cohort[baseline_cohort != 0]
    cdef cnp.ndarray[cnp.int64_t, ndim=1] pool_indices = np.where(baseline_cohort != 0)[1]
    cdef int i
    cdef int ind
    cdef list index_to_remove = []
    cdef set ARS_set = set(ARS)

    for i, ind in enumerate(pool_indices):
        if ind in ARS_set:
            index_to_remove.append(i)

    cdef cnp.ndarray[cnp.int64_t, ndim=1] temp_pool_indices = np.zeros(pool_indices.shape[0], dtype=np.int64)
    cdef Py_ssize_t n = 0

    for ind in pool_indices:
        if ind not in ARS_set:
            temp_pool_indices[n] = ind
            n += 1

    pool_indices = np.array(temp_pool_indices[:n])

    index_to_remove_set = set(index_to_remove)

    pool = np.array([pool[i] for i in range(len(pool)) if i not in index_to_remove_set])
    cdef cnp.ndarray[cnp.int64_t, ndim=1] index
    cdef int counter
    cdef cnp.ndarray[cnp.float64_t, ndim=1] pool_copy
    cdef cnp.ndarray[cnp.int64_t, ndim=1] pool_indices_copy
    cdef cnp.ndarray[cnp.int64_t, ndim=1] similar_specie
    for shuffled_sample in shuffled_samples:
        pool_copy = pool.copy()
        pool_indices_copy = pool_indices.copy()
        counter = 0
        while counter < stop:
            # Generate a random index in the range [0, len(pool_indices_copy) - 1]
            index = random.randint(0, pool_indices_copy.shape[0] - 1)
            shuffled_sample[pool_indices_copy[index]] = pool_copy[index]
            similar_specie = np.where(pool_indices_copy == pool_indices_copy[index])[0]
            pool_copy = np.delete(pool_copy, similar_specie)
            pool_indices_copy = np.delete(pool_indices_copy, similar_specie)
            counter += 1

    # normalize the shuffled_samples
    shuffled_samples = normalize_cohort(shuffled_samples)

    return shuffled_samples

def normalize_cohort(cohort):
    # normalization function
    if cohort.ndim == 1:
        cohort_normalized = cohort / cohort.sum()
    else:
        cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
    return cohort_normalized
