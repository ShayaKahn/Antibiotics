import numpy as np
from scipy.stats import pearsonr
from EranElinav import normalize_cohort, calc_distance, species_richness
from data import *

def num_survived_species_version_2(samples, threshold=0):
    # Stack the samples into a matrix using numpy's vstack
    sample_matrix = np.vstack(samples)

    # Check if all elements in each column are non-zero and above the threshold
    bool_mat = np.all((sample_matrix != 0) & (sample_matrix > threshold), axis=0)

    # Get the indices where all elements are non-zero
    survived_species = np.where(bool_mat)[0]

    return survived_species

def corr_graphs(baseline_cohort, ABX_cohort, follow_up_list, method='change', t=0.7):
    """
    baseline_cohort: numpy matrix, each row is the baseline sample of a subject
    ABX_cohort: numpy matrix, each row is the ABX sample of a subject
    follow_up_list: list of numpy matrices, each matrix contains in rows the follow-up samples
                      (ordered in accending order in time) of the subject
    method: optional values are: 'change', 'richness', 'survived', 'fraction'.
    """
    D_vals = []
    SR = []
    # condition
    non_zero_counts_base = np.count_nonzero(placebo_data_v2, axis=1)
    non_zero_counts_abx = np.count_nonzero(placebo_data_v3, axis=1)
    non_zero_counts_follow = np.count_nonzero(placebo_data_v5, axis=1)
    chosen_indices = np.where(((non_zero_counts_follow > non_zero_counts_base) & (np.divide(
        non_zero_counts_base, non_zero_counts_follow) > t) | (
            non_zero_counts_follow <= non_zero_counts_base) & (
            np.divide(non_zero_counts_follow, non_zero_counts_base) > t)) &
                              ((non_zero_counts_base > non_zero_counts_abx) &
                              (non_zero_counts_follow > non_zero_counts_abx)))[0]
    baseline_cohort = baseline_cohort[chosen_indices, :]
    ABX_cohort = ABX_cohort[chosen_indices, :]
    follow_up_list = follow_up_list[chosen_indices, :]
    for (base, abx, follow) in zip(baseline_cohort, ABX_cohort, follow_up_list):
        if method == 'change':
            SR.append(species_richness(abx) - species_richness(base))
        elif method == 'fraction':
            SR.append(np.divide(species_richness(abx), species_richness(base)))
        elif method == 'survived':
            ss = num_survived_species_version_2([base, abx, follow])
            SR.append(len(ss))
        elif method == 'richness':
            SR.append(species_richness(abx))
        if method == 'survived':
            follow_to_baseline_distance = calc_distance(normalize_cohort(np.delete(follow, ss)),
                                                        normalize_cohort(np.delete(base, ss)))

            D_vals.append(follow_to_baseline_distance)
        else:
            follow_to_baseline_distance = calc_distance(follow, base)

            D_vals.append(follow_to_baseline_distance)

    latex_font = {"family": "serif", "size": 40, "color": "black"}
    fig = go.Figure(data=go.Scatter(x=SR, y=D_vals, mode='markers', marker=dict(color='black', size=15)))

    fig.update_layout(height=600, width=600)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(title_text='Bray Curtis', title_font=latex_font)
    if method == 'change':
        correlation_coefficient, p_value = pearsonr(SR, D_vals)
        print(correlation_coefficient)
        print(p_value)
        fig.update_xaxes(title_text='Δ Species richness', title_font=latex_font)
    elif method == 'richness':
        correlation_coefficient, p_value = pearsonr(SR, D_vals)
        print(correlation_coefficient)
        print(p_value)
        fig.update_xaxes(title_text='Species richness ABX', title_font=latex_font)
    elif method == 'fraction':
        correlation_coefficient, p_value = pearsonr(SR, D_vals)
        print(correlation_coefficient)
        print(p_value)
        fig.update_xaxes(title_text='Species richness ABX fraction', title_font=latex_font)
    else:
        fig.update_xaxes(title_text='Survived species', title_font=latex_font)
        correlation_coefficient, p_value = pearsonr(SR, D_vals)
        print(correlation_coefficient)
        print(p_value)
    return fig, SR, D_vals

fig_1 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='change')
fig_2 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='richness')
fig_3 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='survived')
fig_4, SR, D_vals = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='fraction')