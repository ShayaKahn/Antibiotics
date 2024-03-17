import pandas as pd
import numpy as np
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\Gut\Data')
data = pd.read_excel('Gut_Bacterial_Microbiota_OTU_table.xlsx')
import plotly.graph_objects as go
from EranElinav import normalize_cohort, calc_distance, num_survived_species, species_richness
import openpyxl
from scipy.stats import pearsonr

def split_dataframe_by_column_name(df):
    # Get columns containing the letter 'A'
    columns_A = [col for col in df.columns if 'A' in col]
    df_A = df[columns_A]

    # Get columns containing the letter 'B'
    columns_B = [col for col in df.columns if 'B' in col]
    df_B = df[columns_B]

    return df_A, df_B

def plot_column_sums_histogram(df, title):
    # Calculate the sums of the columns
    column_sums = df.sum(axis=0)

    # Create the histogram
    fig = go.Figure(data=[go.Histogram(x=column_sums)])

    # Update layout for better visualization
    fig.update_layout(
        title=title,
        xaxis_title='Frequency per sample',
        yaxis_title='Number of samples',
        template='plotly_dark',
        width=600,
        height=600,
        font=dict(family='Computer Modern', size=14, color='white'),
        xaxis_showgrid=False,
        yaxis_showgrid=False
    )

    return fig

placebo_data, probiotics_data = split_dataframe_by_column_name(data)

col_sums_plot = plot_column_sums_histogram(placebo_data, title='Placebo')
#col_sums_plot.show()

from rarify import Rarify
placebo_data = Rarify(placebo_data).rarify()

def custom_sort(col_name):
    # Remove the '-' and the following letter at the end
    col_name = col_name.split('-')[0]

    # Split the cleaned column name by 'v' and convert the parts to integers
    num1, num2 = map(int, col_name.split('v'))
    return (num1, num2)

def order_dataframe_columns(df):
    # Sort the columns using the custom_sort function
    sorted_columns = sorted(df.columns, key=custom_sort)
    return df[sorted_columns]

placebo_data = order_dataframe_columns(placebo_data)

def split_to_visit(df):
    df_v2 = pd.DataFrame()
    df_v3 = pd.DataFrame()
    df_v4 = pd.DataFrame()
    df_v5 = pd.DataFrame()
    for col_name in df.columns:
        if 'v2' in col_name:
            df_v2[col_name] = df[col_name]
        elif 'v3' in col_name:
            df_v3[col_name] = df[col_name]
        elif 'v4' in col_name:
            df_v4[col_name] = df[col_name]
        else:
            df_v5[col_name] = df[col_name]
    return df_v2, df_v3, df_v4, df_v5

placebo_data_v2, placebo_data_v3, placebo_data_v4, placebo_data_v5 = split_to_visit(placebo_data)

#placebo_data_v2.to_excel('placebo_data_v2.xlsx', index=False)
#placebo_data_v3.to_excel('placebo_data_v3.xlsx', index=False)
#placebo_data_v4.to_excel('placebo_data_v4.xlsx', index=False)
#placebo_data_v5.to_excel('placebo_data_v5.xlsx', index=False)

placebo_data_v2 = placebo_data_v2.values.T
placebo_data_v3 = placebo_data_v3.values.T
placebo_data_v4 = placebo_data_v4.values.T
placebo_data_v5 = placebo_data_v5.values.T

placebo_data_v2 = normalize_cohort(placebo_data_v2)
placebo_data_v3 = normalize_cohort(placebo_data_v3)
placebo_data_v4 = normalize_cohort(placebo_data_v4)
placebo_data_v5 = normalize_cohort(placebo_data_v5)

#follow_up_list = [placebo_data_v5[0, :], placebo_data_v5[1, :], placebo_data_v5[2, :],
#                  placebo_data_v5[3, :], placebo_data_v5[4, :], placebo_data_v5[5, :],
#                  placebo_data_v5[6, :], placebo_data_v5[7, :], placebo_data_v5[8, :],
#                  placebo_data_v5[9, :], placebo_data_v5[10, :], placebo_data_v5[11, :],
#                  placebo_data_v5[12, :], placebo_data_v5[13, :], placebo_data_v5[14, :],
#                  placebo_data_v5[15, :], placebo_data_v5[16, :], placebo_data_v5[17, :],
#                  placebo_data_v5[18, :], placebo_data_v5[19, :], placebo_data_v5[20, :],
#                  placebo_data_v5[21, :], placebo_data_v5[22, :], placebo_data_v5[23, :],
#                  placebo_data_v5[24, :], placebo_data_v5[25, :], placebo_data_v5[26, :],
#                  placebo_data_v5[27, :], placebo_data_v5[28, :], placebo_data_v5[29, :],
#                  placebo_data_v5[30, :], placebo_data_v5[31, :], placebo_data_v5[32, :],
#                  placebo_data_v5[33, :], placebo_data_v5[34, :]]

follow_up_list = placebo_data_v5

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