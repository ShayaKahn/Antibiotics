from Classes.overlap import Overlap
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

n_min = 2
n_max = 100
n_vals = np.arange(n_min, n_max + 1, 1)

p_min = 0.1
p_max = 0.9
gup = 0.1
p_vals = np.arange(p_min, p_max + gup, gup)

reals = 20

fig = make_subplots(rows=3, cols=3, subplot_titles=[f"p = {p}" for p in p_vals],
                    vertical_spacing=0.1, horizontal_spacing=0.1)

for k, p in enumerate(p_vals):
    print(p)
    jaccard_means = []
    jaccard_stds = []

    for n in n_vals:
        mat = np.random.choice([0, 1], size=(reals, n), p=[1 - p, p])
        jaccard = []
        for i in range(mat.shape[0] - 1):
            for j in range(i + 1, mat.shape[0]):
                overlap_object = Overlap(mat[i, :], mat[j, :], overlap_type="Jaccard")
                jaccard.append(overlap_object.calculate_overlap())

        jaccard_means.append(np.mean(jaccard))
        jaccard_stds.append(np.std(jaccard))
    fig.add_trace(go.Scatter(x=n_vals, y=jaccard_means, error_y=dict(type='data', array=jaccard_stds, visible=True),
                  mode='markers+lines', marker=dict(color='RoyalBlue'), line=dict(color='RoyalBlue')),
                  row=k//3 + 1, col=k % 3 + 1)

fig.update_layout(
    template='plotly_white',
    width=1500,
    height=1500,
    xaxis_title='n',
    yaxis_title='Jaccard similarity',
    xaxis=dict(showline=True, linecolor='black', mirror=True),
    yaxis=dict(showline=True, linecolor='black', mirror=True),
    font=dict(size=18)
    )

fig.show()
"""
# probabilities
p_min = 0.1
p_max = 0.9
gup = 0.01
p_vals = np.arange(p_min, p_max + gup, gup)

n = 100

reals = 50

jaccard_means = []
jaccard_stds = []

for p in p_vals:
    mat = np.random.choice([0, 1], size=(reals, n), p=[1 - p, p])
    jaccard = []
    for i in range(mat.shape[0] - 1):
        for j in range(i + 1, mat.shape[0]):
            overlap_object = Overlap(mat[i, :], mat[j, :], overlap_type="Jaccard")
            jaccard.append(overlap_object.calculate_overlap())

    jaccard_means.append(np.mean(jaccard))
    jaccard_stds.append(np.std(jaccard))

fig = go.Figure(data=go.Scatter(
    x=p_vals,
    y=jaccard_means,
    error_y=dict(type='data', array=jaccard_stds, visible=True),
    mode='markers+lines',
    marker=dict(color='LightSkyBlue'),
    line=dict(color='RoyalBlue')  #
))

fig.update_layout(
    template='plotly_white',
    width=1500,
    height=1500,
    xaxis_title='p',
    yaxis_title='Jaccard similarity',
    xaxis=dict(showline=True, linecolor='black', mirror=True),
    yaxis=dict(showline=True, linecolor='black', mirror=True),
    font=dict(size=18)
)

fig.show()
"""

a = 1