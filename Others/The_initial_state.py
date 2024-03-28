import pandas as pd
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\The initial state\Data')
data = pd.read_excel('Species_table.xlsx')
from EranElinav import normalize_cohort, species_richness
from Recovery import correlation_plot_recovery

total_cohort = data.values

total_cohort = total_cohort[:, 1:]
total_cohort = normalize_cohort(total_cohort)

baseline = [total_cohort[0, :], total_cohort[3, :], total_cohort[6, :], total_cohort[9, :],
            total_cohort[12, :], total_cohort[15, :], total_cohort[18, :], total_cohort[21, :],
            total_cohort[24, :], total_cohort[27, :], total_cohort[30, :], total_cohort[33, :],
            total_cohort[36, :], total_cohort[39, :], total_cohort[42, :], total_cohort[45, :],
            total_cohort[48, :], total_cohort[51, :], total_cohort[54, :], total_cohort[57, :],
            total_cohort[60, :], total_cohort[63, :], total_cohort[66, :], total_cohort[69, :]]

ABX = [total_cohort[1, :], total_cohort[4, :], total_cohort[7, :], total_cohort[10, :],
       total_cohort[13, :], total_cohort[16, :], total_cohort[19, :], total_cohort[22, :],
       total_cohort[25, :], total_cohort[28, :], total_cohort[31, :], total_cohort[34, :],
       total_cohort[37, :], total_cohort[40, :], total_cohort[43, :], total_cohort[46, :],
       total_cohort[49, :], total_cohort[52, :], total_cohort[55, :], total_cohort[58, :],
       total_cohort[61, :], total_cohort[64, :], total_cohort[67, :], total_cohort[70, :]]

follow_up = [total_cohort[2, :], total_cohort[5, :], total_cohort[8, :], total_cohort[11, :],
             total_cohort[14, :], total_cohort[17, :], total_cohort[20, :], total_cohort[23, :],
             total_cohort[26, :], total_cohort[29, :], total_cohort[32, :], total_cohort[35, :],
             total_cohort[38, :], total_cohort[41, :], total_cohort[44, :], total_cohort[47, :],
             total_cohort[50, :], total_cohort[53, :], total_cohort[56, :], total_cohort[59, :],
             total_cohort[62, :], total_cohort[65, :], total_cohort[68, :], total_cohort[71, :]]

# Create the correlation plots
fig_1 = correlation_plot_recovery(baseline, ABX, follow_up, change=False)
fig_2 = correlation_plot_recovery(baseline, ABX, follow_up, change=True)
fig_3 = correlation_plot_recovery(baseline, ABX, follow_up, change=False, survived=True)

# Calculate lists of the species richness at baseline, follow-up and ABX
num_species_base_init = [species_richness(smp) for smp in baseline]
num_species_follow_init = [species_richness(smp) for smp in follow_up]
num_species_ABX_init = [species_richness(smp) for smp in ABX]

# Save the lists
import os
os.chdir(r'C:\Users\shaya\OneDrive\Desktop\Antibiotics project\The initial state\Data')

df_num_species_base_init = pd.DataFrame({'Species richness': num_species_base_init})
df_num_species_follow_init = pd.DataFrame({'Species richness': num_species_follow_init})
df_num_species_ABX_init = pd.DataFrame({'Species richness': num_species_ABX_init})

df_num_species_base_init.to_excel('num_species_base_init.xlsx', index=False)
df_num_species_follow_init.to_excel('num_species_follow_init.xlsx', index=False)
df_num_species_ABX_init.to_excel('num_species_ABX_init.xlsx', index=False)
