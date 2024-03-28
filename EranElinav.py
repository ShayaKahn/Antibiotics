import os
os.chdir(r'C:\Users\USER\Desktop\Antibiotics')
import plotly.graph_objects as go
from data import *
from plotly.subplots import make_subplots
import numpy as np
from overlap import Overlap
from scipy.spatial.distance import braycurtis, jensenshannon
from pairwaise_matrix import Pairwise

# Define a function that plots the Bray Curtis, rJSD, Jaccard, and Overlap measures
# of a time series data from a subject to the baseline of other issues and himself.

def create_recovery_plot(rJSD_vals_list, BC_vals_list, jaccard_vals_list,
                         overlap_vals_list, labels, title, index):

    # Create subplots with shared x-axis
    fig = make_subplots(rows=4, cols=1, shared_xaxes=True, shared_yaxes=False,
                                 vertical_spacing=0, row_heights=[0.2, 0.2, 0.2, 0.2])

    for i, (rJSD_vals, BC_vals, jaccard_vals, overlap_vals) in enumerate(zip(rJSD_vals_list,
                                                                             BC_vals_list,
                                                                             jaccard_vals_list,
                                                                             overlap_vals_list)):
        if i == index:
            color = 'red'
        else:
            color = 'gray'

        fig.add_trace(go.Scatter(x=labels, y=jaccard_vals, mode='markers+lines',
                                          marker=dict(color=color, size=4, symbol="circle"),
                                          line=dict(color=color, width=2)), row=1, col=1)

        fig.add_trace(go.Scatter(x=labels, y=overlap_vals, mode='markers+lines',
                                          marker=dict(color=color, size=4, symbol="circle"),
                                          line=dict(color=color, width=2)), row=2, col=1)

        fig.add_trace(go.Scatter(x=labels, y=BC_vals, mode='markers+lines',
                                          marker=dict(color=color, size=4, symbol="circle"),
                                          line=dict(color=color, width=2)), row=3, col=1)

        fig.add_trace(go.Scatter(x=labels, y=rJSD_vals, mode='markers+lines',
                                          marker=dict(color=color, size=4, symbol="circle"),
                                          line=dict(color=color, width=2)), row=4, col=1)

    # Update subplot layout
    fig.update_layout(height=1000, width=800, title_text=title,
                               title_font=dict(size=25), template="plotly_dark")

    # Set y-axis labels
    latex_font = {"family": "serif", "size": 40, "color": "white"}
    fig.update_yaxes(title_text='Jaccard', row=1, col=1, title_font=latex_font,
                     tickfont=dict(family='Times New Roman', size=25))
    fig.update_yaxes(title_text='Overlap', row=2, col=1, title_font=latex_font,
                     tickfont=dict(family='Times New Roman', size=25))
    fig.update_yaxes(title_text='Bray Curtis', row=3, col=1, title_font=latex_font,
                     tickfont=dict(family='Times New Roman', size=25))
    fig.update_yaxes(title_text='rJSD', row=4, col=1, title_font=latex_font,
                     tickfont=dict(family='Times New Roman', size=25))
    fig.update_xaxes(title_text='Day', row=4, col=1, title_font=latex_font,
                     tickfont=dict(family='Times New Roman', size=25))

    # Set x and y axes lines to black
    fig.update_xaxes(showline=True, linewidth=1, row=1, col=1, linecolor='white',
                              mirror=False, tickvals=[])
    fig.update_xaxes(showline=True, linewidth=1, row=2, col=1, linecolor='white',
                              mirror=False, tickvals=[])
    fig.update_xaxes(showline=True, linewidth=1, row=3, col=1, linecolor='white',
                              mirror=False, tickvals=[])
    fig.update_xaxes(showline=True, linewidth=1, row=4, col=1, linecolor='white',
                              mirror=False)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='white', mirror=False)

    # Remove legend
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, showgrid=False)
    fig.update_yaxes(showline=True, zeroline=False, showgrid=False)

    # Show the plot
    return fig

def recovery(base_samples_list, samples, labels, title, index):
    """
    base_samples_list: list of baseline numpy matices of shape(# samples, # species)
    samples: numpy matrix of the samples of test subject of shape(# samples, # species),
             the samples are time ordered
    labels: the labels represent the time point corresponds to the samples
    title: the title of the plot
    index: the index of the test subject in the base_samples_list
    """
    rJSD_vals_list = []
    BC_vals_list = []
    jaccard_vals_list = []
    overlap_vals_list = []
    for base_samples in base_samples_list:
        # Calculate rJSD values
        rJSD_vals = [np.mean([jensenshannon(base_sample, sample
                                            ) for base_sample in base_samples if not np.array_equal(
            base_sample, sample)]) for sample in samples]
        # Calculate Bray Curtis values
        BC_vals = [np.mean([braycurtis(base_sample, sample
                                       ) for base_sample in base_samples if not np.array_equal(
            base_sample, sample)]) for sample in samples]
        # Calculate Jaccard values
        jaccard_vals = [np.mean([Overlap(base_sample, sample, overlap_type="Jaccard"
                                         ).calculate_overlap(
        ) for base_sample in base_samples if not np.array_equal(
            base_sample, sample)]) for sample in samples]
        # Calculate Overlap values
        overlap_vals = [np.mean([Overlap(base_sample, sample).calculate_overlap()
                                 for base_sample in base_samples if not np.array_equal(base_sample,
                                                                                       sample)]) for sample in samples]
        rJSD_vals_list.append(rJSD_vals)
        BC_vals_list.append(BC_vals)
        jaccard_vals_list.append(jaccard_vals)
        overlap_vals_list.append(overlap_vals)

    # Create the plot
    figure = create_recovery_plot(rJSD_vals_list, BC_vals_list, jaccard_vals_list,
                                                    overlap_vals_list, labels, title, index)

    return figure

# Labels of all subjects
spo_1_labels = ['-13', '-12', '-11', '-10', '-7', '-5', '-4', '-3', '-2', '-1', '0', '3', '4',
                '5', '7', '14', '21', '28', '42', '56', '146', '176', '206', '236']
spo_2_labels = ['-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2',
                '-1', '0', '1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56',
                '146', '176', '206', '236']
spo_3_labels = ['-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1',
                '0', '1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56', '176', '206',
                '236']
spo_4_labels = ['-13', '-12', '-11', '-10', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '0',
                '1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56']
spo_5_labels = ['-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1',
                '0', '1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56', '176', '206',
                '236']
spo_6_labels = ['-12', '-11', '-10', '-8', '-7', '-5', '-4', '-3', '-2', '-1', '0', '1', '2',
                '3', '4', '6', '7', '14', '21', '28', '42', '56']
spo_7_labels = ['-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1',
                '0', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56']

# Create the plots
spo_1_recovery_comp = recovery(baseline_samples_list, spo_1_samples,
                               labels=spo_1_labels, title='Spontaneous 1 recovery', index=14)
spo_2_recovery_comp = recovery(baseline_samples_list, spo_2_samples,
                               labels=spo_2_labels, title=' ',#title='Spontaneous 2 recovery',
                               index=15)
spo_3_recovery_comp = recovery(baseline_samples_list, spo_3_samples,
                               labels=spo_3_labels, title='Spontaneous 3 recovery', index=16)
spo_4_recovery_comp = recovery(baseline_samples_list, spo_4_samples,
                               labels=spo_4_labels, title='Spontaneous 4 recovery', index=17)
spo_5_recovery_comp = recovery(baseline_samples_list, spo_5_samples,
                               labels=spo_5_labels, title=' ',#title='Spontaneous 5 recovery',
                               index=18)
spo_6_recovery_comp = recovery(baseline_samples_list, spo_6_samples,
                               labels=spo_6_labels, title='Spontaneous 6 recovery', index=19)
spo_7_recovery_comp = recovery(baseline_samples_list, spo_7_samples,
                               labels=spo_7_labels, title='Spontaneous 7 recovery', index=20)

# Define function that creates heatmaps
def heatmap_from_distance_matrix(matrix, title=None, labels=None, similarity=False):
    # Set the values of the upper triangle to be 1
    if similarity:
        matrix[np.triu_indices_from(matrix)] = 0
    else:
        matrix[np.triu_indices_from(matrix)] = 1

    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=labels,
        y=labels,
        colorscale=[
            [0, 'darkgreen'],
            [0.25, 'green'],
            [0.5, 'lightgreen'],
            [0.75, 'palegreen'],
            [1, 'white']
        ],
        colorbar=dict(
            tickfont=dict(size=25)
        ),
        showscale=True
    ))

    # Update the layout
    fig.update_layout(
        xaxis=dict(ticks='', side='top', tickfont=dict(family='Times New Roman', size=30)),
        yaxis=dict(ticks='', side='left', autorange='reversed',
                   tickfont=dict(family='Times New Roman', size=30)),
        width=800,
        height=800,
        template="plotly_dark",
        title=dict(text=title, font=dict(color='white'))
    )

    return fig

# Create the heatmaps of Bray-Curtis distances
spo_1_pairwise_matrix_bc = Pairwise(spo_1_samples).create_pairwise_matrix()
spo_1_heatmap_time_bc = heatmap_from_distance_matrix(spo_1_pairwise_matrix_bc,
                                                     title='Spontaneous 1 Bray Curtis',
                                                     labels=spo_1_labels)
spo_2_pairwise_matrix_bc = Pairwise(spo_2_samples).create_pairwise_matrix()
spo_2_heatmap_time_bc = heatmap_from_distance_matrix(spo_2_pairwise_matrix_bc,
                                                     title=' ',#title='Spontaneous 2 Bray Curtis',
                                                     labels=spo_2_labels)
spo_3_pairwise_matrix_bc = Pairwise(spo_3_samples).create_pairwise_matrix()
spo_3_heatmap_time_bc = heatmap_from_distance_matrix(spo_3_pairwise_matrix_bc,
                                                     title='Spontaneous 3 Bray Curtis',
                                                     labels=spo_3_labels)
spo_4_pairwise_matrix_bc = Pairwise(spo_4_samples).create_pairwise_matrix()
spo_4_heatmap_time_bc = heatmap_from_distance_matrix(spo_4_pairwise_matrix_bc,
                                                     title='Spontaneous 4 Bray Curtis',
                                                     labels=spo_4_labels)
spo_5_pairwise_matrix_bc = Pairwise(spo_5_samples).create_pairwise_matrix()
spo_5_heatmap_time_bc = heatmap_from_distance_matrix(spo_5_pairwise_matrix_bc,
                                                     title=' ',#title='Spontaneous 5 Bray Curtis',
                                                     labels=spo_5_labels)
spo_6_pairwise_matrix_bc = Pairwise(spo_6_samples).create_pairwise_matrix()
spo_6_heatmap_time_bc = heatmap_from_distance_matrix(spo_6_pairwise_matrix_bc,
                                                     title='Spontaneous 6 Bray Curtis',
                                                     labels=spo_6_labels)
spo_7_pairwise_matrix_bc = Pairwise(spo_7_samples).create_pairwise_matrix()
spo_7_heatmap_time_bc = heatmap_from_distance_matrix(spo_7_pairwise_matrix_bc,
                                                     title='Spontaneous 7 Bray Curtis',
                                                     labels=spo_7_labels)

# Survived species function
def num_survived_species(samples, threshold=0):
    # Stack the samples into a matrix using numpy's vstack
    sample_matrix = np.vstack(samples)

    # Check if all elements in each column are non-zero and above the threshold
    bool_mat = np.all((sample_matrix != 0) & (sample_matrix > threshold), axis=0)

    # Get the indices where all elements are non-zero
    survived_species = np.where(bool_mat)[0]

    return len(survived_species)

# Species richness function
def species_richness(sample):
    """
    Count the number of non-zero values in a numpy array.

    Parameters:
    - sample (numpy array): A numpy array of species abundances.

    Returns:
    - int: Number of non-zero values (species).
    """
    return np.count_nonzero(sample)

# Normalization function
def normalize_cohort(cohort):
    if cohort.ndim == 1:
        cohort_normalized = cohort / cohort.sum()
    else:
        cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
    return cohort_normalized

# Mean baseline distance function
def mean_lower_triangle(matrix):
    """
    Calculate the mean value of the lower triangle of a square matrix.

    Parameters:
    - matrix (numpy array): A square numpy matrix.

    Returns:
    - float: Mean value of the lower triangle.
    """
    # Extract the lower triangle, excluding the diagonal
    lower_triangle = np.tril(matrix, k=-1)

    # Calculate the mean of the non-zero values in the lower triangle
    mean_value = lower_triangle[lower_triangle != 0].mean()

    return mean_value

# Bray-Curtis distance function
def calc_distance(sample_first, sample_sec):
    return braycurtis(sample_first, sample_sec)

# CM function
def mean_sample(baseline_cohort):
    return np.mean(baseline_cohort, axis=0)

# Correlation plot function
def correlation_plot(baseline_samples_list, ABX_list, follow_up_list, baseline_cohort,
                     samples, metric='braycurtis', method='optimal', delta=True, measure='fraction'):
    """
    baseline_samples_list: list of numpy matrices, each matrix contains in rows the baseline samples
                           of the subject
    ABX_list: list that contain the perturbated sample (last day of ABX) of each subject
    follow_up_list: list of numpy matrices, each matrix contains in rows the follow-up samples
                      (ordered in ascending order in time) of the subject
    baseline_cohort: cohort that contain the "optimal" baseline samples of the subjects
    samples: list of numpy matrices that contain in rows all samples of all the subjects
    metric: the metric method to calculate the Delta_D values
    delta: if True, calculate the change in distance
    method: the options are:
            - optimal: calculate the Delta_D values using the optimal baseline samples
            - cm: calculate the Delta_D values using the CM method
            - mean: calculate the Delta_D values using the mean of the baseline samples
    change: if True, calculate the change in species richness at ABX state
    """
    Delta_D_vals = []
    SR = []
    if method == 'optimal':
        optimal = baseline_cohort
    elif method == 'cm':
        optimal = np.vstack([mean_sample(smp) for smp in baseline_samples_list])
    elif method == 'mean':
        optimal = baseline_samples_list
    for (base, abx, follow, opt, sample) in zip(baseline_samples_list, ABX_list, follow_up_list,
                                                optimal, samples):

        # Calculate the mean distance between the baseline samples
        char_baseline_distance = mean_lower_triangle(Pairwise(base).create_pairwise_matrix(
            metric=metric))

        if not method == 'mean':
            follow_to_baseline_distance = np.mean([calc_distance(follow_sample,
                                                                 opt) for follow_sample in follow])
        else:
            follow_to_baseline_distance = np.mean([np.mean([calc_distance(follow_sample,
                                                                          base) for base in opt]
                                                           ) for follow_sample in follow])

        ########################
        if delta:
            Delta_D_vals.append(follow_to_baseline_distance - char_baseline_distance)
        else:
            Delta_D_vals.append(follow_to_baseline_distance)

        if measure == 'change':
            SR.append(species_richness(abx) - np.mean(
                [species_richness(base_smp) for base_smp in base]))
        elif measure == 'fraction':
            SR.append(np.divide(species_richness(abx), np.mean([species_richness(
                base_smp) for base_smp in base])))
        elif measure == 'richness':
            SR.append(species_richness(abx))
        #######################

    latex_font = {"family": "serif", "size": 30, "color": "black"}
    fig = go.Figure(data=go.Scatter(x=SR, y=Delta_D_vals, mode='markers', marker=dict(size=20)))
    fig.update_layout(height=600, width=600)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    if delta:
        fig.update_yaxes(title_text='Δ Bray Curtis', title_font=latex_font)
    else:
        fig.update_yaxes(title_text='Bray Curtis', title_font=latex_font)
    if measure == 'change':
        fig.update_xaxes(title_text='Δ Species richness', title_font=latex_font)
    elif measure == 'fraction':
        fig.update_xaxes(title_text='Species richness fraction at ABX state', title_font=latex_font)
    elif measure == 'richness':
        fig.update_xaxes(title_text='Species richness at ABX state', title_font=latex_font)

    return fig, SR, Delta_D_vals

# Create correlation plots

last = 3
spo_days_after_ABX_list_last = [mat[-last:] for mat in spo_days_after_ABX_list]

fig = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=False,
                         method='mean', measure='fraction')
fig_2 = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=False,
                         method='optimal', measure='fraction')
fig_3 = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=False,
                         method='cm', measure='fraction')
fig_4 = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=True,
                         method='mean', measure='fraction')
fig_5 = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=True,
                         method='optimal', measure='fraction')
fig_6, SR, Delta_D_vals = correlation_plot(baseline_samples_list=base_sam_spo_list,
                         ABX_list=list(ABX_cohort[14:, :]),
                         follow_up_list=spo_days_after_ABX_list_last,
                         baseline_cohort=baseline_cohort[14:, :],
                         samples=total_spo_samples,
                         metric='braycurtis', delta=True,
                         method='cm', measure='fraction')

# Correlation plot function for Survived species
def correlation_plot_Survived_species(baseline_samples_list, ABX_list, follow_up_list, baseline_cohort,
                         samples, metric='braycurtis', delta=True, method='optimal', include_post=True,
                         last=3):
    """
    baseline_samples_list: list of numpy matrices, each matrix contains in rows the baseline samples
                           of the subject
    ABX_list: list that contain the perturbated sample (last day of ABX) of each subject
    follow_up_list: list of numpy matrices, each matrix contains in rows the follow-up samples
                      (ordered in ascending order in time) of the subject
    baseline_cohort: cohort that contain the "optimal" baseline samples of the subjects
    samples: list of numpy matrices that contain in rows all samples of all the subjects
    metric: the metric method to calculate the Delta_D values
    delta: if True, calculate the change in distance
    method: the options are:
            - optimal: calculate the Delta_D values using the optimal baseline samples
            - cm: calculate the Delta_D values using the CM method
            - mean: calculate the Delta_D values using the mean of the baseline samples
    include_post: if True, the post ABX samples are included in the calculation of the Survived species
    last: integer, the number of days from the end of the experiment to include in the calculation of the
          Survived species
    """
    Delta_D_vals = []
    SR = []
    if method == 'optimal':
        optimal = baseline_cohort
    elif method == 'cm':
        optimal = np.vstack([mean_sample(smp) for smp in baseline_samples_list])
    elif method == 'mean':
        optimal = baseline_samples_list
    for (base, abx, follow, opt, sample) in zip(baseline_samples_list, ABX_list, follow_up_list,
                                                optimal, samples):

        char_baseline_distance = mean_lower_triangle(Pairwise(base).create_pairwise_matrix(
            metric=metric))

        if not method == 'mean':
            follow_to_baseline_distance = np.mean([calc_distance(follow_sample,
                                                                 opt) for follow_sample in follow])
        else:
            follow_to_baseline_distance = np.mean([np.mean([calc_distance(follow_sample,
                                                                          base) for base in opt]
                                                           ) for follow_sample in follow])
        if delta:
            Delta_D_vals.append(follow_to_baseline_distance - char_baseline_distance)
        else:
            Delta_D_vals.append(follow_to_baseline_distance)
        if include_post:
            SR.append(num_survived_species([abx, base, follow[-last, :]]))
        else:
            SR.append(num_survived_species([abx, base]))

    latex_font = {"family": "serif", "size": 30, "color": "black"}
    fig = go.Figure(data=go.Scatter(x=SR, y=Delta_D_vals, mode='markers', marker=dict(size=20)))
    fig.update_layout(height=600, width=600)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    if delta:
        fig.update_yaxes(title_text='Δ Bray Curtis', title_font=latex_font)
    else:
        fig.update_yaxes(title_text='Bray Curtis', title_font=latex_font)
    fig.update_xaxes(title_text='Survived species', title_font=latex_font)
    return fig

# Create the correlation plots for the Survived species
fig_10 = correlation_plot_Survived_species(baseline_samples_list=base_sam_spo_list,
                                           ABX_list=list(ABX_cohort[14:, :]),
                                           follow_up_list=spo_days_after_ABX_list,
                                           baseline_cohort=baseline_cohort[14:, :],
                                           samples=total_spo_samples,
                                           metric='braycurtis', delta=True,
                                           method='mean')
fig_11 = correlation_plot_Survived_species(baseline_samples_list=base_sam_spo_list,
                                           ABX_list=list(ABX_cohort[14:, :]),
                                           follow_up_list=spo_days_after_ABX_list,
                                           baseline_cohort=baseline_cohort[14:, :],
                                           samples=total_spo_samples,
                                           metric='braycurtis', delta=True,
                                           method='optimal')
fig_12 = correlation_plot_Survived_species(baseline_samples_list=base_sam_spo_list,
                                           ABX_list=list(ABX_cohort[14:, :]),
                                           follow_up_list=spo_days_after_ABX_list,
                                           baseline_cohort=baseline_cohort[14:, :],
                                           samples=total_spo_samples,
                                           metric='braycurtis', delta=True,
                                           method='cm')

# Calculate lists of the species richness at baseline, follow-up and ABX
last = 3
follow_up_list_last = [mat[-last:] for mat in spo_days_after_ABX_list]
num_species_base = [np.mean([species_richness(smp) for smp in baseline_samples]
                            ) for baseline_samples in base_sam_spo_list]
num_species_follow = [np.mean([species_richness(smp) for smp in follow_up]
                      ) for follow_up in follow_up_list_last]
num_species_ABX = [species_richness(smp) for smp in list(ABX_cohort[14:, :])]

# Save the lists
import os
os.chdir(r'C:\Users\USER\Desktop\Antibiotics\Eran Elinav\Data')

df_num_species_base = pd.DataFrame({'Species richness': num_species_base})
df_num_species_follow = pd.DataFrame({'Species richness': num_species_follow})
df_num_species_ABX = pd.DataFrame({'Species richness': num_species_ABX})

df_num_species_base.to_excel('num_species_base.xlsx', index=False)
df_num_species_follow.to_excel('num_species_follow.xlsx', index=False)
df_num_species_ABX.to_excel('num_species_ABX.xlsx', index=False)
