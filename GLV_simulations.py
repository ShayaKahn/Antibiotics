from Classes.Historical_contingency import HC
from scipy.stats import spearmanr, pearsonr
import numpy as np

num_samples = 100
pool_size = 50
num_survived = 20#40
mean = 0.0
sigma = 0.1
lower = 0.1
higher = 0.7
delta = 1e-2
final_time = 1e5
max_step = 0.1
epsilon = 1e-3
threshold = 1e-3
min_skew = -15
max_skew = 0
interaction_strength = 3.5#2
max_growth = 0.6
min_growth = 0
method = 'RK45'
measure = 'Jaccard'#'intersection'#

HC_object = HC(num_samples, pool_size, num_survived, mean, sigma, lower, higher, delta, final_time, max_step,
               epsilon, threshold, min_skew, max_skew, interaction_strength, max_growth, min_growth,
               method, measure)
results = HC_object.get_results()

ind = np.where(np.array(results["similarity_perturbed"]) > 0.25)

corr, p_value = spearmanr(np.array(results["similarity_perturbed"])[ind], np.array(results["similarity_real"])[ind])

print("Spearman's correlation coefficient real:", corr)
print("P-value real:", p_value)

corr, p_value = pearsonr(np.array(results["similarity_perturbed"])[ind], np.array(results["similarity_real"])[ind])

print("Pearsons's correlation coefficient real:", corr)
print("P-value real:", p_value)

corr, p_value = spearmanr(np.array(results["similarity_perturbed"])[ind], np.array(results["similarity_shuffled"])[ind])

print("Spearman's correlation coefficient shuffled:", corr)
print("P-value shuffled:", p_value)

corr, p_value = pearsonr(np.array(results["similarity_perturbed"])[ind], np.array(results["similarity_shuffled"])[ind])

print("Pearson's correlation coefficient shuffled:", corr)
print("P-value shuffled:", p_value)

import matplotlib.pyplot as plt
# Creating subplots
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Plotting on the first subplot
ax[0].scatter(np.array(results["similarity_perturbed"])[ind],
              np.array(results["similarity_real"])[ind], color='blue', label='Jaccard Real')
ax[0].set_title('Jaccard Perturbed vs. Real')
ax[0].set_xlabel('Jaccard Perturbed')
ax[0].set_ylabel('Jaccard Real')
ax[0].grid(True)
ax[0].legend()

# Plotting on the second subplot
ax[1].scatter(np.array(results["similarity_perturbed"])[ind],
              np.array(results["similarity_shuffled"])[ind],
              color='green', label='similarity Shuffled')
ax[1].set_title('Jaccard Perturbed vs. Shuffled')
ax[1].set_xlabel('Jaccard Perturbed')
ax[1].set_ylabel('Jaccard Shuffled')
ax[1].grid(True)
ax[1].legend()

plt.tight_layout()
plt.show()

print("No filtering")

corr, p_value = spearmanr(np.array(results["similarity_perturbed"]), np.array(results["similarity_real"]))

print("Spearman's correlation coefficient real:", corr)
print("P-value real:", p_value)

corr, p_value = pearsonr(np.array(results["similarity_perturbed"]), np.array(results["similarity_real"]))

print("Pearsons's correlation coefficient real:", corr)
print("P-value real:", p_value)

corr, p_value = spearmanr(np.array(results["similarity_perturbed"]), np.array(results["similarity_shuffled"]))

print("Spearman's correlation coefficient shuffled:", corr)
print("P-value shuffled:", p_value)

corr, p_value = pearsonr(np.array(results["similarity_perturbed"]), np.array(results["similarity_shuffled"]))

print("Pearson's correlation coefficient shuffled:", corr)
print("P-value shuffled:", p_value)

import matplotlib.pyplot as plt
# Creating subplots
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Plotting on the first subplot
ax[0].scatter(np.array(results["similarity_perturbed"]),
              np.array(results["similarity_real"]), color='blue', label='Jaccard Real')
ax[0].set_title('Jaccard Perturbed vs. Real')
ax[0].set_xlabel('Jaccard Perturbed')
ax[0].set_ylabel('Jaccard Real')
ax[0].grid(True)
ax[0].legend()

# Plotting on the second subplot
ax[1].scatter(np.array(results["similarity_perturbed"]),
              np.array(results["similarity_shuffled"]),
              color='green', label='similarity Shuffled')
ax[1].set_title('Jaccard Perturbed vs. Shuffled')
ax[1].set_xlabel('Jaccard Perturbed')
ax[1].set_ylabel('Jaccard Shuffled')
ax[1].grid(True)
ax[1].legend()

plt.tight_layout()
plt.show()

#import matplotlib.pyplot as plt
#import numpy as np

#ind = np.where(results["filtered_post_perturbed_state"] > 0.0)
#A = results["filtered_post_perturbed_state"].copy()
#A[ind] = 1

#plt.imshow(A, aspect='auto')

#plt.show()
a = 1