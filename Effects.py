import pandas as pd
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\Effects of different amoxicillin\Data')
data = pd.read_excel('ASV_table.xlsx')
from EranElinav import normalize_cohort, recovery, species_richness
from Recovery import correlation_plot_recovery

# Separate the data into groups
A_columns = [col for col in data.columns if 'A' in col]
B_columns = [col for col in data.columns if 'B' in col]
C_columns = [col for col in data.columns if 'C' in col]
D_columns = [col for col in data.columns if 'D' in col]

group_A_cohort = data[A_columns]
group_B_cohort = data[B_columns]
group_C_cohort = data[C_columns]
group_D_cohort = data[D_columns]

# Separete each group into individuals
A1_columns = [col for col in group_A_cohort.columns if 'A1' in col]
A2_columns = [col for col in group_A_cohort.columns if 'A2' in col]
A3_columns = [col for col in group_A_cohort.columns if 'A3' in col]
A4_columns = [col for col in group_A_cohort.columns if 'A4' in col]
A5_columns = [col for col in group_A_cohort.columns if 'A5' in col]

A1_cohort = group_A_cohort[A1_columns]
A2_cohort = group_A_cohort[A2_columns]
A3_cohort = group_A_cohort[A3_columns]
A4_cohort = group_A_cohort[A4_columns]
A5_cohort = group_A_cohort[A5_columns]

B1_columns = [col for col in group_B_cohort.columns if 'B1' in col]
B2_columns = [col for col in group_B_cohort.columns if 'B2' in col]
B3_columns = [col for col in group_B_cohort.columns if 'B3' in col]
B4_columns = [col for col in group_B_cohort.columns if 'B4' in col]
B5_columns = [col for col in group_B_cohort.columns if 'B5' in col]

B1_cohort = group_B_cohort[B1_columns]
B2_cohort = group_B_cohort[B2_columns]
B3_cohort = group_B_cohort[B3_columns]
B4_cohort = group_B_cohort[B4_columns]
B5_cohort = group_B_cohort[B5_columns]

C1_columns = [col for col in group_C_cohort.columns if 'C1' in col]
C2_columns = [col for col in group_C_cohort.columns if 'C2' in col]
C3_columns = [col for col in group_C_cohort.columns if 'C3' in col]
C4_columns = [col for col in group_C_cohort.columns if 'C4' in col]
C5_columns = [col for col in group_C_cohort.columns if 'C5' in col]

C1_cohort = group_C_cohort[C1_columns]
C2_cohort = group_C_cohort[C2_columns]
C3_cohort = group_C_cohort[C3_columns]
C4_cohort = group_C_cohort[C4_columns]
C5_cohort = group_C_cohort[C5_columns]

D1_columns = [col for col in group_D_cohort.columns if 'D1' in col]
D2_columns = [col for col in group_D_cohort.columns if 'D2' in col]
D3_columns = [col for col in group_D_cohort.columns if 'D3' in col]
D4_columns = [col for col in group_D_cohort.columns if 'D4' in col]
D5_columns = [col for col in group_D_cohort.columns if 'D5' in col]

D1_cohort = group_D_cohort[D1_columns]
D2_cohort = group_D_cohort[D2_columns]
D3_cohort = group_D_cohort[D3_columns]
D4_cohort = group_D_cohort[D4_columns]
D5_cohort = group_D_cohort[D5_columns]

# Order the columns by the time of sampling

def order_columns(df):
    # Extract the suffix numbers from column names and convert to integers
    suffix_numbers = [int(col.split('S')[-1]) for col in df.columns]

    # Sort columns based on suffix numbers
    sorted_columns = [col for _, col in sorted(zip(suffix_numbers, df.columns))]

    # Create a new DataFrame with sorted columns
    ordered_df = df[sorted_columns]

    return ordered_df

A1_cohort = order_columns(A1_cohort)
A2_cohort = order_columns(A2_cohort)
A3_cohort = order_columns(A3_cohort)
A4_cohort = order_columns(A4_cohort)
A5_cohort = order_columns(A5_cohort)

B1_cohort = order_columns(B1_cohort)
B2_cohort = order_columns(B2_cohort)
B3_cohort = order_columns(B3_cohort)
B4_cohort = order_columns(B4_cohort)
B5_cohort = order_columns(B5_cohort)

C1_cohort = order_columns(C1_cohort)
C2_cohort = order_columns(C2_cohort)
C3_cohort = order_columns(C3_cohort)
C4_cohort = order_columns(C4_cohort)
C5_cohort = order_columns(C5_cohort)

D1_cohort = order_columns(D1_cohort)
D2_cohort = order_columns(D2_cohort)
D3_cohort = order_columns(D3_cohort)
D4_cohort = order_columns(D4_cohort)
D5_cohort = order_columns(D5_cohort)

days_A = [1, 5, 9, 13, 17, 21, 25, 29, 37]
days_B = [1, 5, 9, 13, 17, 21, 25]
days_C = [1, 5, 9, 13, 17, 21, 25, 29]
days_D = [1, 5, 9, 13, 17, 21, 25, 29, 37]

A1_cohort = normalize_cohort(A1_cohort.values.T)
A2_cohort = normalize_cohort(A2_cohort.values.T)
A3_cohort = normalize_cohort(A3_cohort.values.T)
A4_cohort = normalize_cohort(A4_cohort.values.T)
A5_cohort = normalize_cohort(A5_cohort.values.T)

B1_cohort = normalize_cohort(B1_cohort.values.T)
B2_cohort = normalize_cohort(B2_cohort.values.T)
B3_cohort = normalize_cohort(B3_cohort.values.T)
B4_cohort = normalize_cohort(B4_cohort.values.T)
B5_cohort = normalize_cohort(B5_cohort.values.T)

C1_cohort = normalize_cohort(C1_cohort.values.T)
C2_cohort = normalize_cohort(C2_cohort.values.T)
C3_cohort = normalize_cohort(C3_cohort.values.T)
C4_cohort = normalize_cohort(C4_cohort.values.T)
C5_cohort = normalize_cohort(C5_cohort.values.T)

D1_cohort = normalize_cohort(D1_cohort.values.T)
D2_cohort = normalize_cohort(D2_cohort.values.T)
D3_cohort = normalize_cohort(D3_cohort.values.T)
D4_cohort = normalize_cohort(D4_cohort.values.T)
D5_cohort = normalize_cohort(D5_cohort.values.T)

Effects_baseline_samples_list = [[A1_cohort[0, :]], [A2_cohort[0, :]], [A3_cohort[0, :]],
                                 [A4_cohort[0, :]], [A5_cohort[0, :]], [B1_cohort[0, :]],
                                 [B2_cohort[0, :]], [B3_cohort[0, :]], [B4_cohort[0, :]],
                                 [B5_cohort[0, :]], [C1_cohort[0, :]], [C2_cohort[0, :]],
                                 [C3_cohort[0, :]], [C4_cohort[0, :]], [C5_cohort[0, :]],
                                 [D1_cohort[0, :]], [D2_cohort[0, :]], [D3_cohort[0, :]],
                                 [D4_cohort[0, :]], [D5_cohort[0, :]]]

fig_1 = recovery(Effects_baseline_samples_list, B1_cohort[1:, :], labels=days_B[1:],
                 title='B1 recovery', index=5)
fig_2 = recovery(Effects_baseline_samples_list, B2_cohort[1:, :], labels=days_B[1:],
                 title='B2 recovery', index=6)
fig_3 = recovery(Effects_baseline_samples_list, B3_cohort[1:, :], labels=days_B[1:],
                 title='B3 recovery', index=7)
fig_4 = recovery(Effects_baseline_samples_list, B4_cohort[1:, :], labels=days_B[1:],
                 title='B4 recovery', index=8)
fig_5 = recovery(Effects_baseline_samples_list, B5_cohort[1:, :], labels=days_B[1:],
                 title='B5 recovery', index=9)
fig_6 = recovery(Effects_baseline_samples_list, C1_cohort[1:, :], labels=days_C[1:],
                 title='C1 recovery', index=10)
fig_7 = recovery(Effects_baseline_samples_list, C2_cohort[1:, :], labels=days_C[1:],
                 title='C2 recovery', index=11)
fig_8 = recovery(Effects_baseline_samples_list, C3_cohort[1:, :], labels=days_C[1:],
                 title='C3 recovery', index=12)
fig_9 = recovery(Effects_baseline_samples_list, C4_cohort[1:, :], labels=days_C[1:],
                 title='C4 recovery', index=13)
fig_10 = recovery(Effects_baseline_samples_list, C5_cohort[1:, :], labels=days_C[1:],
                 title='C5 recovery', index=14)
fig_11 = recovery(Effects_baseline_samples_list, D1_cohort[1:, :], labels=days_D[1:],
                 title='D1 recovery', index=15)
fig_12 = recovery(Effects_baseline_samples_list, D2_cohort[1:, :], labels=days_D[1:],
                 title='D2 recovery', index=16)
fig_13 = recovery(Effects_baseline_samples_list, D3_cohort[1:, :], labels=days_D[1:],
                 title='D3 recovery', index=17)
fig_14 = recovery(Effects_baseline_samples_list, D4_cohort[1:, :], labels=days_D[1:],
                 title='D4 recovery', index=18)
fig_15 = recovery(Effects_baseline_samples_list, D5_cohort[1:, :], labels=days_D[1:],
                 title='D5 recovery', index=19)


baseline = [B1_cohort[0, :], B2_cohort[0, :], B3_cohort[0, :], B4_cohort[0, :],
            B5_cohort[0, :], C1_cohort[0, :], C2_cohort[0, :], C3_cohort[0, :],
            C4_cohort[0, :], C5_cohort[0, :], D1_cohort[0, :], D2_cohort[0, :],
            D3_cohort[0, :], D4_cohort[0, :], D5_cohort[0, :]]

ABX = [B1_cohort[1, :], B2_cohort[1, :], B3_cohort[1, :], B4_cohort[1, :],
       B5_cohort[1, :], C1_cohort[2, :], C2_cohort[2, :], C3_cohort[2, :],
       C4_cohort[2, :], C5_cohort[2, :], D1_cohort[3, :], D2_cohort[3, :],
       D3_cohort[3, :], D4_cohort[3, :], D5_cohort[3, :]]

follow_up = [B1_cohort[-1, :], B2_cohort[-1, :], B3_cohort[-1, :], B4_cohort[-1, :],
             B5_cohort[-1, :], C1_cohort[-1, :], C2_cohort[-1, :], C3_cohort[-1, :],
             C4_cohort[-1, :], C5_cohort[-1, :], D1_cohort[-1, :], D2_cohort[-1, :],
             D3_cohort[-1, :], D4_cohort[-1, :], D5_cohort[-1, :]]

# Create the correlation plots
fig_16 = correlation_plot_recovery(baseline, ABX, follow_up, change=False)
fig_17 = correlation_plot_recovery(baseline, ABX, follow_up, change=True)
fig_18 = correlation_plot_recovery(baseline, ABX, follow_up, change=False, survived=True)

# Calculate lists of the species richness at baseline, follow-up and ABX
num_species_base_effects = [species_richness(smp) for smp in baseline]
num_species_follow_effects = [species_richness(smp) for smp in follow_up]
num_species_ABX_effects = [species_richness(smp) for smp in ABX]

# Save the lists
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\Effects of different amoxicillin\Data')

df_num_species_base_effects = pd.DataFrame({'Species richness': num_species_base_effects})
df_num_species_follow_effects = pd.DataFrame({'Species richness': num_species_follow_effects})
df_num_species_ABX_effects = pd.DataFrame({'Species richness': num_species_ABX_effects})

df_num_species_base_effects.to_excel('num_species_base_effects.xlsx', index=False)
df_num_species_follow_effects.to_excel('num_species_follow_effects.xlsx', index=False)
df_num_species_ABX_effects.to_excel('num_species_ABX_effects.xlsx', index=False)