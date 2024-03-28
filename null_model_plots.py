from Classes.null_model import NullModel
import plotly.graph_objects as go
from data import *
"""
# null model for spontaneous 1
spo_base = baseline_cohort[15, :]
spo_base_sum = np.sum(spo_2_baseline, axis=0)
spo_abx = ABX_cohort[15, :]
spo_post = spo_2_samples[-1, :]

spo_obj = NullModel(spo_base, spo_abx, spo_post, baseline_cohort, num_reals=1000,
                    baseline_sample_sum=spo_base_sum, ABX_sample_sum=None, exclude_test=False)
bc_real_spo, bc_shuffled_spo = spo_obj.distance(method='Bray Curtis')

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

import time

start_time = time.time()
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

base = baseline_rel_abund_rarefied_appear_4[3, :]
abx = rel_abund_rarefied_4[3, :]
post = rel_abund_rarefied_180_appear_4[3, :]

obj = NullModel(base, abx, post, baseline, num_reals=1000)
bc_real, bc_shuffled = obj.distance(method='Bray Curtis')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time} seconds")

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



#import os
#os.chdir(r'C:\Users\USER\Desktop\Antibiotics\Gut\Data')
"""
data = pd.read_excel('Gut_Bacterial_Microbiota_OTU_table.xlsx')
import plotly.graph_objects as go
from EranElinav import normalize_cohort


def split_dataframe_by_column_name(df):
    # Get columns containing the letter 'A'
    columns_A = [col for col in df.columns if 'A' in col]
    df_A = df[columns_A]

    # Get columns containing the letter 'B'
    columns_B = [col for col in df.columns if 'B' in col]
    df_B = df[columns_B]

    return df_A, df_B

placebo_data, probiotics_data = split_dataframe_by_column_name(data)

#from rarify import Rarify
#placebo_data = Rarify(placebo_data).rarify()

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

placebo_data_v2 = placebo_data_v2.values.T
placebo_data_v3 = placebo_data_v3.values.T
placebo_data_v4 = placebo_data_v4.values.T
placebo_data_v5 = placebo_data_v5.values.T

placebo_data_v2 = normalize_cohort(placebo_data_v2)
placebo_data_v3 = normalize_cohort(placebo_data_v3)
placebo_data_v4 = normalize_cohort(placebo_data_v4)
placebo_data_v5 = normalize_cohort(placebo_data_v5)

base = placebo_data_v2[5, :]
abx = placebo_data_v3[5, :]
post = placebo_data_v5[5, :]

import time
start_time = time.time()

obj = NullModel(base, abx, post, placebo_data_v2, num_reals=1000, ABX_sample_sum=None, exclude_test=False)
bc_real, bc_shuffled = obj.distance(method='Bray Curtis')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time} seconds")

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