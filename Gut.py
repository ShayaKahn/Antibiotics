import pandas as pd
import numpy as np
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\Gut\Data')
data = pd.read_excel('Gut_Bacterial_Microbiota_OTU_table.xlsx')
import plotly.graph_objects as go
from EranElinav import normalize_cohort, calc_distance, num_survived_species, species_richness
import openpyxl

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

placebo_data_v2.to_excel('placebo_data_v2.xlsx', index=False)
placebo_data_v3.to_excel('placebo_data_v3.xlsx', index=False)
placebo_data_v4.to_excel('placebo_data_v4.xlsx', index=False)
placebo_data_v5.to_excel('placebo_data_v5.xlsx', index=False)

placebo_data_v2 = placebo_data_v2.values.T
placebo_data_v3 = placebo_data_v3.values.T
placebo_data_v4 = placebo_data_v4.values.T
placebo_data_v5 = placebo_data_v5.values.T

placebo_data_v2 = normalize_cohort(placebo_data_v2)
placebo_data_v3 = normalize_cohort(placebo_data_v3)
placebo_data_v4 = normalize_cohort(placebo_data_v4)
placebo_data_v5 = normalize_cohort(placebo_data_v5)

follow_up_list = [placebo_data_v5[0, :], placebo_data_v5[1, :], placebo_data_v5[2, :],
                  placebo_data_v5[3, :], placebo_data_v5[4, :], placebo_data_v5[5, :],
                  placebo_data_v5[6, :], placebo_data_v5[7, :], placebo_data_v5[8, :],
                  placebo_data_v5[9, :], placebo_data_v5[10, :], placebo_data_v5[11, :],
                  placebo_data_v5[12, :], placebo_data_v5[13, :], placebo_data_v5[14, :],
                  placebo_data_v5[15, :], placebo_data_v5[16, :], placebo_data_v5[17, :],
                  placebo_data_v5[18, :], placebo_data_v5[19, :], placebo_data_v5[20, :],
                  placebo_data_v5[21, :], placebo_data_v5[22, :], placebo_data_v5[23, :],
                  placebo_data_v5[24, :], placebo_data_v5[25, :], placebo_data_v5[26, :],
                  placebo_data_v5[27, :], placebo_data_v5[28, :], placebo_data_v5[29, :],
                  placebo_data_v5[30, :], placebo_data_v5[31, :], placebo_data_v5[32, :],
                  placebo_data_v5[33, :], placebo_data_v5[34, :]]

from scipy.spatial.distance import braycurtis, jensenshannon
from scipy.stats import pearsonr


def permutation_test(x, y, n_permutations=10000):
    observed_corr, _ = pearsonr(x, y)
    permuted_corrs = np.zeros(n_permutations)

    for i in range(n_permutations):
        y_permuted = np.random.permutation(y)
        permuted_corrs[i], _ = pearsonr(x, y_permuted)

    p_value = np.sum(np.abs(permuted_corrs) >= np.abs(observed_corr)) / n_permutations
    return p_value


def partial_correlation(x, y, z):
    # Step 1: Fit a linear regression model for x and z, and find the residuals
    coef_xz = np.polyfit(z, x, 1)
    residuals_x = x - np.polyval(coef_xz, z)

    # Step 2: Fit a linear regression model for y and z, and find the residuals
    coef_yz = np.polyfit(z, y, 1)
    residuals_y = y - np.polyval(coef_yz, z)

    # Step 3: Calculate the correlation between the residuals
    partial_corr, _ = pearsonr(residuals_x, residuals_y)

    return partial_corr


def calculate_correlation(list1, list2):
    # Calculate Pearson correlation
    corr_coefficient, p_value = pearsonr(list1, list2)

    return corr_coefficient, p_value


def corr_graphs(baseline_cohort, ABX_cohort, follow_up_list, method='change'):
    """
    baseline_cohort: numpy matrix, each row is the baseline sample of a subject
    ABX_cohort: numpy matrix, each row is the ABX sample of a subject
    follow_up_list: list of numpy matrices, each matrix contains in rows the follow-up samples
                      (ordered in accending order in time) of the subject
    method: optional values are: 'change', 'richness', 'survived'.
    """
    D_vals = []
    SR = []
    for (base, abx, follow) in zip(baseline_cohort, ABX_cohort, follow_up_list):

        follow_to_baseline_distance = calc_distance(follow, base)

        D_vals.append(follow_to_baseline_distance)

        if method == 'change':
            SR.append(species_richness(abx) - species_richness(base))
        elif method == 'survived':
            SR.append(num_survived_species([abx, base, follow]))
        elif method == 'richness':
            SR.append(species_richness(abx))

    latex_font = {"family": "serif", "size": 40, "color": "black"}
    fig = go.Figure(data=go.Scatter(x=SR, y=D_vals, mode='markers', marker=dict(size=15)))

    fig.update_layout(height=600, width=600)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(title_text='Bray Curtis', title_font=latex_font)
    if method == 'change':
        fig.update_xaxes(title_text='Δ Species richness', title_font=latex_font)
    elif method == 'richness':
        fig.update_xaxes(title_text='Species richness', title_font=latex_font)
    else:
        fig.update_xaxes(title_text='Survived species', title_font=latex_font)
    return fig#, D_vals, SR

fig_1 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='change')
fig_2 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='richness')
fig_3 = corr_graphs(placebo_data_v2, placebo_data_v3, follow_up_list, method='survived')