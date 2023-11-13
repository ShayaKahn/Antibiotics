from data import *
from EranElinav import recovery, create_recovery_plot, calc_distance, species_richness,\
    num_survived_species
import plotly.graph_objects as go

recovery_baseline_samples_list = [[baseline_rel_abund_rarefied[0, :]],
                                  [baseline_rel_abund_rarefied[1, :]],
                                  [baseline_rel_abund_rarefied[2, :]],
                                  [baseline_rel_abund_rarefied[3, :]],
                                  [baseline_rel_abund_rarefied[4, :]],
                                  [baseline_rel_abund_rarefied[5, :]],
                                  [baseline_rel_abund_rarefied[6, :]],
                                  [baseline_rel_abund_rarefied[7, :]],
                                  [baseline_rel_abund_rarefied[8, :]],
                                  [baseline_rel_abund_rarefied[9, :]],
                                  [baseline_rel_abund_rarefied[10, :]],
                                  [baseline_rel_abund_rarefied[11, :]]]

labels_list = ['D4', 'D8', 'D42', 'D180']

# Create the plots
fig_1 = recovery(recovery_baseline_samples_list, subject_2_samples, labels=labels_list,
                 title='Subject 1 recovery', index=1)
fig_2 = recovery(recovery_baseline_samples_list, subject_3_samples, labels=labels_list,
                 title='Subject 2 recovery', index=2)
fig_3 = recovery(recovery_baseline_samples_list, subject_4_samples, labels=labels_list,
                 title='Subject 3 recovery', index=3)
fig_4 = recovery(recovery_baseline_samples_list, subject_5_samples, labels=labels_list,
                 title='Subject 4 recovery', index=4)
fig_5 = recovery(recovery_baseline_samples_list, subject_6_samples,
                 labels=labels_list, title='Subject 5 recovery', index=5)
fig_6 = recovery(recovery_baseline_samples_list, subject_7_samples, labels=labels_list,
                 title='Subject 6 recovery', index=6)
fig_7 = recovery(recovery_baseline_samples_list, subject_9_samples, labels=labels_list,
                 title='Subject 7 recovery', index=8)
fig_8 = recovery(recovery_baseline_samples_list, subject_11_samples, labels=labels_list,
                 title='Subject 8 recovery', index=10)
fig_9 = recovery(recovery_baseline_samples_list, subject_12_samples, labels=labels_list,
                 title='Subject 9 recovery', index=11)

# Correlation plot function
def correlation_plot_recovery(baseline_cohort, ABX_cohort, follow_up_list, change=False, survived=False):
    """
    baseline_cohort: numpy matrix, each row is the baseline sample of a subject
    ABX_cohort: numpy matrix, each row is the ABX sample of a subject
    follow_up_list: list of numpy matrices, each matrix contains in rows the follow-up samples
                      (ordered in accending order in time) of the subject
    change: if True, calculate the change in species richness at ABX state
    survived: if True, calculate the number of survived species
    """
    D_vals = []
    SR = []
    for (base, abx, follow) in zip(baseline_cohort, ABX_cohort, follow_up_list):

        if follow_up_list[0].shape[0] == 2:
            follow_to_baseline_distance = np.mean([calc_distance(follow_sample,
                                                                 base) for follow_sample in follow])
        else:
            follow_to_baseline_distance = calc_distance(follow, base)

        D_vals.append(follow_to_baseline_distance)

        if survived:
            SR.append(num_survived_species([base, abx]))
        else:
            if change:
                SR.append(species_richness(abx) - species_richness(base))
            else:
                SR.append(species_richness(abx))

    latex_font = {"family": "serif", "size": 30, "color": "black"}
    fig = go.Figure(data=go.Scatter(x=SR, y=D_vals, mode='markers', marker=dict(size=20)))
    fig.update_layout(height=600, width=600)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(showline=True, zeroline=False, linecolor='black', showgrid=False,
                     tickfont=dict(size=20))
    fig.update_yaxes(title_text='Bray Curtis', title_font=latex_font)
    if survived:
        fig.update_xaxes(title_text='Survived species', title_font=latex_font)
    else:
        if change:
            fig.update_xaxes(title_text='Δ Species richness', title_font=latex_font)
        else:
            fig.update_xaxes(title_text='Species richness', title_font=latex_font)
    return fig

follow_up_list_last = [np.array(rel_abund_rarefied_180_appear_4[0, :]),
                       np.array(rel_abund_rarefied_180_appear_4[1, :]),
                       np.array(rel_abund_rarefied_180_appear_4[2, :]),
                       np.array(rel_abund_rarefied_180_appear_4[3, :]),
                       np.array(rel_abund_rarefied_180_appear_4[4, :]),
                       np.array(rel_abund_rarefied_180_appear_4[5, :]),
                       np.array(rel_abund_rarefied_180_appear_4[6, :]),
                       np.array(rel_abund_rarefied_180_appear_4[7, :]),
                       np.array(rel_abund_rarefied_180_appear_4[8, :])]

# Create the correlation plots
fig_10 = correlation_plot_recovery(baseline_rel_abund_rarefied_appear_4, rel_abund_rarefied_4,
                                   follow_up_list_last, change=False)
fig_11 = correlation_plot_recovery(baseline_rel_abund_rarefied_appear_4, rel_abund_rarefied_4,
                                   follow_up_list_last, change=True)
fig_12 = correlation_plot_recovery(baseline_rel_abund_rarefied_appear_4, rel_abund_rarefied_4,
                                   follow_up_list_last, change=True, survived=True)

fig_10.show()
fig_11.show()
fig_12.show()

# Calculate lists of the species richness at baseline, follow-up and ABX
num_species_base_recovery = [species_richness(smp) for smp in baseline_rel_abund_rarefied_appear_4]
num_species_follow_recovery = [species_richness(smp) for smp in rel_abund_rarefied_180_appear_4]
num_species_ABX_recovery = [species_richness(smp) for smp in rel_abund_rarefied_4]

# Save the lists
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\Recovery\Data')

df_num_species_base_recovery = pd.DataFrame({'Species richness': num_species_base_recovery})
df_num_species_follow_recovery = pd.DataFrame({'Species richness': num_species_follow_recovery})
df_num_species_ABX_recovery = pd.DataFrame({'Species richness': num_species_ABX_recovery})

df_num_species_base_recovery.to_excel('num_species_base_recovery.xlsx', index=False)
df_num_species_follow_recovery.to_excel('num_species_follow_recovery.xlsx', index=False)
df_num_species_ABX_recovery.to_excel('num_species_ABX_recovery.xlsx', index=False)