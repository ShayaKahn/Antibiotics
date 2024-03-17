from null_model import NullModel
import plotly.graph_objects as go
from data import *

import time

start_time = time.time()

# null model for spontaneous 1
spo_base = baseline_cohort[20, :]
spo_base_sum = np.sum(spo_7_baseline, axis=0)
spo_abx = ABX_cohort[20, :]
spo_post = spo_7_samples[-1, :]

spo_obj = NullModel(spo_base, spo_abx, spo_post, baseline_cohort, num_reals=5000,
                    baseline_sample_sum=spo_base_sum, ABX_sample_sum=None, exclude_test=False)
bc_real_spo, bc_shuffled_spo = spo_obj.distance(method='Bray Curtis')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time} seconds")

p_value = np.sum(bc_shuffled_spo <= bc_real_spo) / len(bc_shuffled_spo)

fig = go.Figure()
fig.add_trace(go.Histogram(x=bc_shuffled_spo, marker_color='blue'))

# Add a red dot for the single numerical value
fig.add_trace(go.Scatter(x=[bc_real_spo], y=[0], mode='markers', marker_color='red',
                         marker_size=20))

fig.update_layout(
    xaxis_title='Bray Curtis Distance',
    yaxis_title='Counts',
    bargap=0.2,
    width=600,
    height=600,
    plot_bgcolor='white',
    xaxis=dict(
        showgrid=False,
        title_font=dict(size=30),
        tickfont=dict(size=30),
        #zeroline=True,
        #zerolinecolor='black',
        #zerolinewidth=2,
        showline=True,
        linecolor='black',
        linewidth=2
    ),
    yaxis=dict(
        showgrid=False,
        title_font=dict(size=30),
        tickfont=dict(size=30),
        #zeroline=True,
        #zerolinecolor='black',
        #zerolinewidth=2,
        showline=True,
        linecolor='black',
        linewidth=2
    ),
    showlegend=False,
    font=dict(
        family="Times New Roman, Times, serif",
        size=12,
        color="black"
    )
)

fig.show()

"""
# null model for subject 1

baseline = np.vstack((baseline_rel_abund_rarefied[0, :],
                      baseline_rel_abund_rarefied[1, :],
                      baseline_rel_abund_rarefied[2, :],
                      baseline_rel_abund_rarefied[3, :],
                      baseline_rel_abund_rarefied[4, :],
                      baseline_rel_abund_rarefied[5, :],
                      baseline_rel_abund_rarefied[6, :],
                      baseline_rel_abund_rarefied[7, :],
                      baseline_rel_abund_rarefied[8, :],
                      baseline_rel_abund_rarefied[9, :],
                      baseline_rel_abund_rarefied[10, :],
                      baseline_rel_abund_rarefied[11, :]))

base = baseline_rel_abund_rarefied_appear_4[0, :]
abx = rel_abund_rarefied_4[0, :]
post = rel_abund_rarefied_180_appear_4[0, :]

obj = NullModel(base, abx, post, baseline, num_reals=1000, ABX_sample_sum=None, exclude_test=False)
bc_real, bc_shuffled = obj.distance(method='Bray Curtis')

p_value_recovery = np.sum(bc_shuffled <= bc_real) / len(bc_shuffled)

fig = go.Figure()
fig.add_trace(go.Histogram(x=bc_shuffled, marker_color='blue'))

# Add a red dot for the single numerical value
fig.add_trace(go.Scatter(x=[bc_real], y=[0], mode='markers', marker_color='red',
                         marker_size=20))

fig.update_layout(
    xaxis_title='Bray Curtis Distance',
    yaxis_title='Counts',
    bargap=0.2,
    width=600,
    height=600,
    plot_bgcolor='white',
    xaxis=dict(
        showgrid=False,
        title_font=dict(size=30),
        tickfont=dict(size=30),
        #zeroline=True,
        #zerolinecolor='black',
        #zerolinewidth=2,
        showline=True,
        linecolor='black',
        linewidth=2
    ),
    yaxis=dict(
        showgrid=False,
        title_font=dict(size=30),
        tickfont=dict(size=30),
        #zeroline=True,
        #zerolinecolor='black',
        #zerolinewidth=2,
        showline=True,
        linecolor='black',
        linewidth=2
    ),
    showlegend=False,
    font=dict(
        family="Times New Roman, Times, serif",
        size=12,
        color="black"
    )
)

fig.show()
"""