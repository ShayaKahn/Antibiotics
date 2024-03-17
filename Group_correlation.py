import numpy as np
from scipy.spatial.distance import pdist, squareform
from data import *
from Gut import *

# Assuming matrices A and B are given, where rows represent samples for each subject in states A and B, respectively.
# For demonstration, let's create dummy data for A and B:
#A = baseline_cohort[14:, :]
#A = np.vstack([spo_1_D_after_ABX[-1 ,:], spo_2_D_after_ABX[-1 ,:], spo_3_D_after_ABX[-1 ,:],
#               spo_4_D_after_ABX[-1 ,:], spo_5_D_after_ABX[-1 ,:], spo_6_D_after_ABX[-1 ,:],
#               spo_7_D_after_ABX[-1 ,:]])
B = placebo_data_v5#rel_abund_rarefied_4#ABX_cohort[14:, :]
A = placebo_data_v2#placebo_data_v5 #baseline_rel_abund_rarefied_appear_4#rel_abund_rarefied_180_appear_4#baseline_cohort[14:, :]

def count_shared_nonzeros(vec1, vec2):
    # Identify non-zero elements in each vector
    nonzeros1 = vec1 != 0
    nonzeros2 = vec2 != 0

    # Find shared non-zero positions
    shared_nonzeros = nonzeros1 & nonzeros2

    # Count shared non-zero elements
    return np.sum(shared_nonzeros)


# Calculate shared non-zero counts for all pairs of rows in B
n = B.shape[0]  # Number of rows in B
shared_nonzeros_count = np.zeros((n, n))  # Initialize matrix to store counts

# Populate the shared_nonzeros_count matrix
for i in range(n):
    for j in range(i + 1, n):  # Only need to calculate for one triangle due to symmetry
        count = count_shared_nonzeros(B[i], B[j])
        shared_nonzeros_count[i, j] = count
        shared_nonzeros_count[j, i] = count

# Calculate pairwise Euclidean distances within each state
distances_A = squareform(pdist(A, metric='braycurtis'))
distances_B = shared_nonzeros_count#squareform(pdist(B, metric='braycurtis'))#

# Extracting the upper triangular parts without the diagonal (which are zeros) to avoid duplicate pairs and self-comparison
# This simplification assumes distances_A and distances_B are symmetric and have zeros on the diagonal
unique_distances_A = distances_A[np.triu_indices(len(A), k=1)]
unique_distances_B = distances_B[np.triu_indices(len(B), k=1)]

import plotly.express as px

# Using the previously calculated unique_distances_A and unique_distances_B
# Creating a DataFrame for Plotly
import pandas as pd
df = pd.DataFrame({'Distance in State A': unique_distances_A, 'Distance in State B': unique_distances_B})

# Recreating the scatter plot with the filtered data
fig_filtered = px.scatter(df, x='Distance in State A', y='Distance in State B',
                          title='Microbiome Sample Distances: State A vs. State B (Filtered)',
                          labels={'Distance in State A': 'Pairwise Distances in State A',
                                  'Distance in State B': 'Pairwise Distances in State B'})

# Update the layout to set figure size and remove grid lines for the filtered plot
fig_filtered.update_layout(width=600, height=600,
                           xaxis_showgrid=False, yaxis_showgrid=False,
                           xaxis_zeroline=False, yaxis_zeroline=False)

# Showing the filtered plot
fig_filtered.show()