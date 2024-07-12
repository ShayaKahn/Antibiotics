from Classes.Historical_contingency_v2 import HC2
from scipy.stats import spearmanr, pearsonr, linregress
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

num_samples = 10
pool_size = 100
num_survived = 50
sigma = 5
mean = 0
c = 0.05
delta = 1e-5
final_time = 1000
max_step = 0.05
epsilon = 1e-3
threshold = 1e-3
min_growth = 0
max_growth = 1
alpha = 1#None
method = 'RK45'
measure = 'Jaccard'#'braycurtis'#
multiprocess = True

import time

init = time.time()

HC_object = HC2(num_samples, pool_size, num_survived, mean, sigma, c, delta, final_time, max_step, epsilon, threshold,
                max_growth, min_growth, alpha, method, measure, multiprocess)
results = HC_object.get_results()

final = time.time()

print("Total time in seconds: ", final - init)

### real ###

s_real, _ = spearmanr(np.array(results["similarity_perturbed"]),
                      np.array(results["real"]))
p_real, _ = pearsonr(np.array(results["similarity_perturbed"]),
                     np.array(results["real"]))
slope_real, intercept_real, _, _, _  = linregress(results["similarity_perturbed"],
                                                  results["real"])

### shuffled ###

s_shuffled, _ = spearmanr(np.array(results["similarity_perturbed"]),
                          np.array(results["shuffled"]))
p_shuffled, _ = pearsonr(np.array(results["similarity_perturbed"]),
                         np.array(results["shuffled"]))
slope_shuffled, intercept_shuffled, _, _, _ = linregress(results["similarity_perturbed"],
                                                         results["shuffled"])

x_values = np.linspace(min(results["similarity_perturbed"]),
                       max(results["similarity_perturbed"]), 100)
y_values_real = slope_real * x_values + intercept_real
y_values_shuffled = slope_shuffled * x_values + intercept_shuffled


### plots ###
fig = make_subplots(rows=1, cols=2, subplot_titles=('Real', 'Shuffled'))

fig.add_trace(go.Scatter(x=results["similarity_perturbed"], y=results["real"],
                         mode='markers', marker=dict(color='blue'), name='Real'),
                         row=1, col=1)
fig.add_trace(go.Scatter(x=x_values, y=y_values_real,
                         mode='lines', line=dict(color='red')),
                         row=1, col=1)
fig.add_trace(go.Scatter(x=results["similarity_perturbed"], y=results["shuffled"],
                         mode='markers', marker=dict(color='blue'), name='Shuffled'),
                         row=1, col=2)
fig.add_trace(go.Scatter(x=x_values, y=y_values_shuffled,
                         mode='lines', line=dict(color='red'),
                         showlegend=False),
                         row=1, col=2)

fig.update_layout(height=500, width=1000, plot_bgcolor='white', paper_bgcolor='white',
                  showlegend=False)

axis_updates = dict(showline=True, linecolor='black', showgrid=False,
                    zeroline=False)
fig.update_xaxes(**axis_updates, title_text='Jaccard ABX',  row=1, col=1)
fig.update_yaxes(**axis_updates, title_text='Jaccard post ABX', row=1, col=1)
fig.update_xaxes(**axis_updates, title_text='Jaccard ABX', row=1, col=2)
fig.update_yaxes(**axis_updates, title_text='Jaccard post ABX', row=1, col=2)

fig.show()



from Classes.Historical_contingency_v2 import HC2
from scipy.stats import spearmanr, pearsonr
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

num_samples = 25
pool_size = 100
num_survived = 50
sigma = 5
mean = 0
c = 0.05
delta = 1e-5
final_time = 1000
max_step = 0.05
epsilon = 1e-3
threshold = 1e-3
min_growth = 0
max_growth = 1
alpha = 1
method = 'RK45'
measure = 'Jaccard'#'braycurtis'#

iters = 5

results_dict = {}
import time
start_time = time.time()
HC_object = HC2(num_samples, pool_size, num_survived, mean, sigma, c, delta, final_time, max_step, epsilon, threshold,
                max_growth, min_growth, alpha, method, measure)
results = HC_object.get_results()

#for a in results["perturbed_state"]:
#    print(np.shape(np.where(a != 0)))

#for a in results["initial_conditions"]:
#    print(np.shape(np.where(a != 0)))

G = nx.from_numpy_array(results["interaction_matrix"])
G.edges(data=True)
nx.draw(G)
plt.show()

results = HC_object.get_results()
s, _ = spearmanr(np.array(results["similarity_perturbed"]),
                 np.array(results["real"]))

print(s)
from scipy.stats import linregress

# Linear regression to find the line of best fit
slope, intercept, r_value, p_value, std_err = linregress(results["similarity_perturbed"],
                                                         results["real"])

# Generate x values from the min to max of perturbed similarities for the line of best fit
x_values = np.linspace(min(results["similarity_perturbed"]), max(results["similarity_perturbed"]), 100)

# Calculate the corresponding y values from the linear regression equation
y_values = slope * x_values + intercept

print(slope)

# Creating the plot
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Scatter plot of the real vs. perturbed similarities
ax.scatter(results["similarity_perturbed"], results["real"], color='blue', label='Real')

# Line of best fit
ax.plot(x_values, y_values, color='red', label=f'Linear Fit: y={slope:.2f}x+{intercept:.2f}')

ax.set_xlabel('Jaccard Perturbed')
ax.set_ylabel('Real')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()

########### shuffled ###########
# Linear regression to find the line of best fit
slope, intercept, r_value, p_value, std_err = linregress(results["similarity_perturbed"],
                                                         results["shuffled"])

# Generate x values from the min to max of perturbed similarities for the line of best fit
x_values = np.linspace(min(results["similarity_perturbed"]), max(results["similarity_perturbed"]), 100)

# Calculate the corresponding y values from the linear regression equation
y_values = slope * x_values + intercept

print(slope)

# Creating the plot
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Scatter plot of the real vs. perturbed similarities
ax.scatter(results["similarity_perturbed"], results["shuffled"], color='blue', label='Shuffled')


# Line of best fit
ax.plot(x_values, y_values, color='red', label=f'Linear Fit: y={slope:.2f}x+{intercept:.2f}')

ax.set_xlabel('Jaccard Perturbed')
ax.set_ylabel('Shuffled')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()
###############################

ind = np.where(HC_object.filtered_post_perturbed_state != 0)
filtered_post_perturbed_state_copy = HC_object.filtered_post_perturbed_state.copy()
filtered_post_perturbed_state_copy[ind] = 1
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.imshow(filtered_post_perturbed_state_copy, aspect='auto')
plt.show()

"""
interaction_strength_vector = np.arange(0.1, 3.3, 0.1)
reals = 20

spearman_corrs = []

for interaction_strength in interaction_strength_vector:
    print(interaction_strength)
    s_val = 0
    for _ in range(reals):
        HC_object = HC(num_samples, pool_size, num_survived_min, num_survived_max, mean, sigma,
                       lower, higher, delta, final_time, max_step, epsilon, threshold, min_skew,
                       max_skew, interaction_strength, max_growth, min_growth, method, measure)
        results = HC_object.get_results()
        s, _ = spearmanr(np.array(results["similarity_perturbed"]),
                         np.array(results["similarity_real"]))
        s_val += s
    spearman_corrs.append(s_val / reals)

spearman_corrs = np.array(spearman_corrs)

plt.scatter(interaction_strength_vector, spearman_corrs, label='Spearman')

plt.show()
"""
#HC_object = HC(num_samples, pool_size, num_survived_min, num_survived_max, mean, sigma,
#                       lower, higher, delta, final_time, max_step, epsilon, threshold, min_skew,
#                       max_skew, interaction_strength, max_growth, min_growth, method, measure)
#results = HC_object.get_results()
#s, _ = spearmanr(np.array(results["similarity_perturbed"]), np.array(results["similarity_real"]))

#corr, p_value = spearmanr(np.array(results["similarity_perturbed"]), np.array(results["similarity_real"]))

#print("Spearman's correlation coefficient real:", corr)
#print("P-value real:", p_value)

#corr, p_value = pearsonr(np.array(results["similarity_perturbed"]), np.array(results["similarity_real"]))

#print("Pearsons's correlation coefficient real:", corr)
#print("P-value real:", p_value)

#corr, p_value = spearmanr(np.array(results["similarity_perturbed"]), np.array(results["similarity_shuffled"]))

#print("Spearman's correlation coefficient shuffled:", corr)
#print("P-value shuffled:", p_value)

#corr, p_value = pearsonr(np.array(results["similarity_perturbed"]), np.array(results["similarity_shuffled"]))

#print("Pearson's correlation coefficient shuffled:", corr)
#print("P-value shuffled:", p_value)

#import matplotlib.pyplot as plt
# Creating subplots
#fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Plotting on the first subplot
#ax[0].scatter(np.array(results["similarity_perturbed"]),
#              np.array(results["similarity_real"]), color='blue', label='Jaccard Real')
#ax[0].set_title('Jaccard Perturbed vs. Real')
#ax[0].set_xlabel('Jaccard Perturbed')
#ax[0].set_ylabel('Jaccard Real')
#ax[0].grid(True)
#ax[0].legend()

# Plotting on the second subplot
#ax[1].scatter(np.array(results["similarity_perturbed"]),
#              np.array(results["similarity_shuffled"]),
#              color='green', label='similarity Shuffled')
#ax[1].set_title('Jaccard Perturbed vs. Shuffled')
#ax[1].set_xlabel('Jaccard Perturbed')
#ax[1].set_ylabel('Jaccard Shuffled')
#ax[1].grid(True)
#ax[1].legend()

#plt.tight_layout()
#plt.show()

#import matplotlib.pyplot as plt
#import numpy as np

#ind = np.where(results["filtered_post_perturbed_state"] > 0.0)
#A = results["filtered_post_perturbed_state"].copy()
#A[ind] = 1

#plt.imshow(A, aspect='auto')

#plt.show()
a = 1





