import os
import pandas as pd
import numpy as np
from Data_manipulation.subject import Subject
import warnings
from Data_manipulation.optimal import OptimalCohort
warnings.filterwarnings("ignore", category=UserWarning)

def normalize_cohort(cohort):
    # normalization function
    if cohort.ndim == 1:
        cohort_normalized = cohort / cohort.sum()
    else:
        cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
    return cohort_normalized

def filter_data(df):
    # filter the data
    columns_to_check = df.columns[1:]
    mean_values = df[columns_to_check].mean(axis=1)
    condition_mask = mean_values >= 0.0001
    df = df[condition_mask]
    return df

### Post-Antibiotic Gut Mucosal Microbiome Reconstitution Is Impaired by Probiotics and Improved by Autologous FMT ###

os.chdir(r"C:\Users\USER\OneDrive\Desktop\Antibiotics\Eran Elinav\Data")
df = pd.read_excel('Metaphlan_stool.xlsx')
Species_column = df[df.columns[0]]
Species = np.array(list(Species_column))

# aFMT data
aFMT_baseline_1 = df.iloc[:, 1:8].values  # 1,2,3,4,5,6,7 of baseline.
aFMT_antibiotics_1 = df.iloc[:, 8].values   # day 2 of antibiotics.
aFMT_intervention_1 = df.iloc[:, 9:19].values  # days 1,2,3,4,5,14,21,28,42,56 of intervention.
aFMT_month_1 = df.iloc[:, 19:24].values  # months after the end of intervention 1,2,3,4,6.

aFMT_baseline_2 = df.iloc[:, 24:30].values  # 1,2,3,5,6,7 of baseline.
aFMT_antibiotics_2 = df.iloc[:, 30:35].values   # days 2,4,5,6,7  of antibiotics.
aFMT_intervention_2 = df.iloc[:, 35:47].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
aFMT_month_2 = df.iloc[:, 47:52].values  # months after the end of intervention 2,3,4,5,6.

aFMT_baseline_3 = df.iloc[:, 52:58].values  # 1,2,3,4,5,7 of baseline.
aFMT_antibiotics_3 = df.iloc[:, 58:64].values   # days 2,3,4,5,6,7  of antibiotics.
aFMT_intervention_3 = df.iloc[:, 64:76].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
aFMT_month_3 = df.iloc[:, 76:80].values  # months after the end of intervention 2,3,4,5.

aFMT_baseline_4 = df.iloc[:, 80:86].values  # 1,3,4,5,6,7 of baseline.
aFMT_antibiotics_4 = df.iloc[:, 86:91].values   # days 1,2,5,6,7  of antibiotics.
aFMT_intervention_4 = df.iloc[:, 91:102].values  # days 3,4,5,6,7,14,21,24,28,42,56 of intervention.
aFMT_intervention_4 = np.delete(aFMT_intervention_4, 7, axis=1)
aFMT_month_4 = df.iloc[:, 102:106].values  # months after the end of intervention 3,4,5,6.

aFMT_baseline_5 = df.iloc[:, 106:113].values  # 1,2,3,4,5,6,7 of baseline.
aFMT_antibiotics_5 = df.iloc[:, 113:119].values   # days 1,2,3,4,5,6  of antibiotics.
aFMT_intervention_5 = df.iloc[:, 119:131].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
aFMT_month_5 = df.iloc[:, 131:134].values  # months after the end of intervention 4,5,6.

aFMT_baseline_6 = df.iloc[:, 134:141].values  # 1,2,3,4,5,6,7 of baseline.
aFMT_antibiotics_6 = df.iloc[:, 141:147].values   # days 2,3,4,5,6,7  of antibiotics.
aFMT_intervention_6 = df.iloc[:, 147:159].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.

# probiotics data

pro_baseline_1 = df.iloc[:, 159:166].values  # 1,2,3,4,5,6,7 of baseline.
pro_antibiotics_1 = df.iloc[:, 166:173].values   # days 1,2,3,4,5,6,7  of antibiotics.
pro_intervention_1 = df.iloc[:, 173:185].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
pro_month_1 = df.iloc[:, 185:190].values  # months after the end of intervention 2,3,4,5,6

pro_baseline_2 = df.iloc[:, 190:195].values  # 1,2,3,4,5 of baseline.
pro_antibiotics_2 = df.iloc[:, 195:201].values   # days 2,3,4,5,6,7  of antibiotics.
pro_intervention_2 = df.iloc[:, 201:208].values  # days 1,2,4,5,6,14,21 of intervention.
pro_month_2 = df.iloc[:, 208:212].values  # months after the end of intervention 2,3,4,5.

pro_baseline_3 = df.iloc[:, 212:218].values  # 1,2,3,4,5,7 of baseline.
pro_antibiotics_3 = df.iloc[:, 218:223].values   # days 1,2,3,4,6  of antibiotics.
pro_intervention_3 = df.iloc[:, 223:234].values  # days 1,2,3,4,5,6,14,21,28,42,56 of intervention.
pro_month_3 = df.iloc[:, 234:237].values  # months after the end of intervention 3,5,6.

pro_baseline_4 = df.iloc[:, 237:244].values  # 1,2,3,4,5,6,7 of baseline.
pro_antibiotics_4 = df.iloc[:, 244:251].values   # days 1,2,3,4,5,6,7  of antibiotics.
pro_intervention_4 = df.iloc[:, 251:261].values  # days 1,2,4,5,6,14,21,28,42,56 of intervention.
pro_month_4 = df.iloc[:, 261:265].values  # months after the end of intervention 3,4,5,6.

pro_baseline_5 = df.iloc[:, 265:272].values  # 1,2,3,4,5,6,7 of baseline.
pro_antibiotics_5 = df.iloc[:, 272:274].values   # days 1,2 of antibiotics.
pro_intervention_5 = df.iloc[:, 274:286].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
pro_month_5 = df.iloc[:, 286:289].values  # months after the end of intervention 4,5,6.

pro_baseline_6 = df.iloc[:, 289:296].values  # 1,2,3,4,5,6,7 of baseline.
pro_antibiotics_6 = df.iloc[:, 296:302].values   # days 1,2,4,5,6,7 of antibiotics.
pro_intervention_6 = df.iloc[:, 302:314].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.

pro_baseline_7 = df.iloc[:, 314:321].values  # 1,2,3,4,5,6,7 of baseline.
pro_antibiotics_7 = df.iloc[:, 321:328].values   # days 1,2,3,4,5,6,7 of antibiotics.
pro_intervention_7 = df.iloc[:, 328:340].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.

pro_baseline_8 = df.iloc[:, 340:343].values  # 1,2,3 of baseline.
pro_antibiotics_8 = df.iloc[:, 343:347].values   # days 2,4,6,7 of antibiotics.
pro_intervention_8 = df.iloc[:, 347:356].values  # days 3,5,6,7,14,21,28,42,56 of intervention.

# spontaneous data

spo_baseline_1 = df.iloc[:, 356:361].values  # 1,2,3,4,7 of baseline.
spo_antibiotics_1 = df.iloc[:, 361:367].values   # days 2,3,4,5,6,7 of antibiotics.
spo_intervention_1 = df.iloc[:, 367:376].values  # days 3,4,5,7,14,21,28,42,56 of intervention.
spo_month_1 = df.iloc[:, 376:380].values  # months after the end of intervention 3,4,5,6.

spo_baseline_2 = df.iloc[:, 380:387].values  # 1,2,3,4,5,6,7 of baseline.
spo_antibiotics_2 = df.iloc[:, 387:394].values   # days 1,2,3,4,5,6,7 of antibiotics.
spo_intervention_2 = df.iloc[:, 394:406].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
spo_month_2 = df.iloc[:, 406:410].values  # months after the end of intervention 3,4,5,6.

spo_baseline_3 = df.iloc[:, 410:417].values  # 1,2,3,4,5,6,7 of baseline.
spo_antibiotics_3 = df.iloc[:, 417:424].values   # days 1,2,3,4,5,6,7 of antibiotics.
spo_intervention_3 = df.iloc[:, 424:436].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
spo_month_3 = df.iloc[:, 436:439].values  # months after the end of intervention 4,5,6.

spo_baseline_4 = df.iloc[:, 439:445].values  # 1,2,3,4,6,7 of baseline.
spo_antibiotics_4 = df.iloc[:, 445:452].values   # days 1,2,3,4,5,6,7 of antibiotics.
spo_intervention_4 = df.iloc[:, 452:464].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.

spo_baseline_5 = df.iloc[:, 464:471].values  # 1,2,3,4,5,6,7 of baseline.
spo_antibiotics_5 = df.iloc[:, 471:478].values   # days 1,2,3,4,5,6,7 of antibiotics.
spo_intervention_5 = df.iloc[:, 478:490].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.
spo_month_5 = df.iloc[:, 490:493].values  # months after the end of intervention 4,5,6.

spo_baseline_6 = df.iloc[:, 493:498].values  # 2,3,4,6,7 of baseline.
spo_antibiotics_6 = df.iloc[:, 498:504].values   # days 2,3,4,5,6,7 of antibiotics.
spo_intervention_6 = df.iloc[:, 504:515].values  # days 1,2,3,4,6,7,14,21,28,42,56 of intervention.

spo_baseline_7 = df.iloc[:, 515:522].values  # 1,2,3,4,5,6,7 of baseline.
spo_antibiotics_7 = df.iloc[:, 522:529].values   # days 1,2,3,4,5,6,7 of antibiotics.
spo_intervention_7 = df.iloc[:, 529:541].values  # days 1,2,3,4,5,6,7,14,21,28,42,56 of intervention.

aFMT_subject_1 = Subject(base_array=aFMT_baseline_1, ant_array=aFMT_antibiotics_1, int_array=aFMT_intervention_1,
                         month_array=aFMT_month_1, base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['2'],
                         int_days=['1', '2', '3', '4', '5', '14', '21', '28', '42', '56'],
                         month_number=['1', '2', '3', '4', '6'])
aFMT_subject_1_base_dict = aFMT_subject_1.create_base_dict()
aFMT_subject_1_ant_dict = aFMT_subject_1.create_ant_dict()
aFMT_subject_1_int_dict = aFMT_subject_1.create_int_dict()
aFMT_subject_1_month_dict = aFMT_subject_1.create_month_dict()

aFMT_subject_2 = Subject(base_array=aFMT_baseline_2, ant_array=aFMT_antibiotics_2, int_array=aFMT_intervention_2,
                         month_array=aFMT_month_2, base_days=['1', '2', '3', '5', '6', '7'],
                         ant_days=['2', '4', '5', '6', '7'],
                         int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                         month_number=['2', '3', '4', '5', '6'])
aFMT_subject_2_base_dict = aFMT_subject_2.create_base_dict()
aFMT_subject_2_ant_dict = aFMT_subject_2.create_ant_dict()
aFMT_subject_2_int_dict = aFMT_subject_2.create_int_dict()
aFMT_subject_2_month_dict = aFMT_subject_2.create_month_dict()

aFMT_subject_3 = Subject(base_array=aFMT_baseline_3, ant_array=aFMT_antibiotics_3, int_array=aFMT_intervention_3,
                         month_array=aFMT_month_3, base_days=['1', '2', '3', '4', '5', '7'],
                         ant_days=['2', '3', '4', '5', '6', '7'],
                         int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                         month_number=['2', '3', '4', '5'])
aFMT_subject_3_base_dict = aFMT_subject_3.create_base_dict()
aFMT_subject_3_ant_dict = aFMT_subject_3.create_ant_dict()
aFMT_subject_3_int_dict = aFMT_subject_3.create_int_dict()
aFMT_subject_3_month_dict = aFMT_subject_3.create_month_dict()

aFMT_subject_4 = Subject(base_array=aFMT_baseline_4, ant_array=aFMT_antibiotics_4, int_array=aFMT_intervention_4,
                         month_array=aFMT_month_4, base_days=['1', '3', '4', '5', '6', '7'],
                         ant_days=['1', '2', '5', '6', '7'],
                         int_days=['3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                         month_number=['3', '4', '5', '6'])
aFMT_subject_4_base_dict = aFMT_subject_4.create_base_dict()
aFMT_subject_4_ant_dict = aFMT_subject_4.create_ant_dict()
aFMT_subject_4_int_dict = aFMT_subject_4.create_int_dict()
aFMT_subject_4_month_dict = aFMT_subject_4.create_month_dict()

aFMT_subject_5 = Subject(base_array=aFMT_baseline_5, ant_array=aFMT_antibiotics_5, int_array=aFMT_intervention_5,
                         month_array=aFMT_month_5, base_days=['1', '2', '3', '4', '5', '6', '7'],
                         ant_days=['1', '2', '3', '4', '5', '6'],
                         int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                         month_number=['4', '5', '6'])
aFMT_subject_5_base_dict = aFMT_subject_5.create_base_dict()
aFMT_subject_5_ant_dict = aFMT_subject_5.create_ant_dict()
aFMT_subject_5_int_dict = aFMT_subject_5.create_int_dict()
aFMT_subject_5_month_dict = aFMT_subject_5.create_month_dict()

aFMT_subject_6 = Subject(base_array=aFMT_baseline_6, ant_array=aFMT_antibiotics_6, int_array=aFMT_intervention_6,
                         base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['2', '3', '4', '5', '6', '7'],
                         int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'])
aFMT_subject_6_base_dict = aFMT_subject_6.create_base_dict()
aFMT_subject_6_ant_dict = aFMT_subject_6.create_ant_dict()
aFMT_subject_6_int_dict = aFMT_subject_6.create_int_dict()

pro_subject_1 = Subject(base_array=pro_baseline_1, ant_array=pro_antibiotics_1, int_array=pro_intervention_1,
                        month_array=pro_month_1,
                        base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        month_number=['2', '3', '4', '5', '6'])
pro_subject_1_base_dict = pro_subject_1.create_base_dict()
pro_subject_1_ant_dict = pro_subject_1.create_ant_dict()
pro_subject_1_int_dict = pro_subject_1.create_int_dict()
pro_subject_1_month_dict = pro_subject_1.create_month_dict()

pro_subject_2 = Subject(base_array=pro_baseline_2, ant_array=pro_antibiotics_2, int_array=pro_intervention_2,
                        month_array=pro_month_2,
                        base_days=['1', '2', '3', '4', '5'], ant_days=['2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '4', '5', '6', '14', '21'], month_number=['2', '3', '4', '5'])
pro_subject_2_base_dict = pro_subject_2.create_base_dict()
pro_subject_2_ant_dict = pro_subject_2.create_ant_dict()
pro_subject_2_int_dict = pro_subject_2.create_int_dict()
pro_subject_2_month_dict = pro_subject_2.create_month_dict()

pro_subject_3 = Subject(base_array=pro_baseline_3, ant_array=pro_antibiotics_3, int_array=pro_intervention_3,
                        month_array=pro_month_3,
                        base_days=['1', '2', '3', '4', '5', '7'], ant_days=['1', '2', '3', '4', '6'],
                        int_days=['1', '2', '3', '4', '5', '6', '14', '21', '28', '42', '56'],
                        month_number=['3', '5', '6'])
pro_subject_3_base_dict = pro_subject_3.create_base_dict()
pro_subject_3_ant_dict = pro_subject_3.create_ant_dict()
pro_subject_3_int_dict = pro_subject_3.create_int_dict()
pro_subject_3_month_dict = pro_subject_3.create_month_dict()

pro_subject_4 = Subject(base_array=pro_baseline_4, ant_array=pro_antibiotics_4, int_array=pro_intervention_4,
                        month_array=pro_month_4,
                        base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '4', '5', '6', '14', '21', '28', '42', '56'],
                        month_number=['3', '4', '5', '6'])
pro_subject_4_base_dict = pro_subject_4.create_base_dict()
pro_subject_4_ant_dict = pro_subject_4.create_ant_dict()
pro_subject_4_int_dict = pro_subject_4.create_int_dict()
pro_subject_4_month_dict = pro_subject_4.create_month_dict()

pro_subject_5 = Subject(base_array=pro_baseline_5, ant_array=pro_antibiotics_5, int_array=pro_intervention_5,
                        month_array=pro_month_5,
                        base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['1', '2'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        month_number=['4', '5', '6'])
pro_subject_5_base_dict = pro_subject_5.create_base_dict()
pro_subject_5_ant_dict = pro_subject_5.create_ant_dict()
pro_subject_5_int_dict = pro_subject_5.create_int_dict()
pro_subject_5_month_dict = pro_subject_5.create_month_dict()

pro_subject_6 = Subject(base_array=pro_baseline_6, ant_array=pro_antibiotics_6, int_array=pro_intervention_6,
                        base_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        ant_days=['1', '2', '4', '5', '6', '7'])
pro_subject_6_int_dict = pro_subject_6.create_int_dict()
pro_subject_6_base_dict = pro_subject_6.create_base_dict()
pro_subject_6_ant_dict = pro_subject_6.create_ant_dict()

pro_subject_7 = Subject(base_array=pro_baseline_7, ant_array=pro_antibiotics_7, int_array=pro_intervention_7,
                        base_days=['1', '2', '3', '4', '5', '6', '7'], ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'])
pro_subject_7_base_dict = pro_subject_7.create_base_dict()
pro_subject_7_ant_dict = pro_subject_7.create_ant_dict()
pro_subject_7_int_dict = pro_subject_7.create_int_dict()

pro_subject_8 = Subject(base_array=pro_baseline_8, int_array=pro_intervention_8, ant_array=pro_antibiotics_8,
                        base_days=['1', '2', '3'], int_days=['3', '5', '6', '7', '14', '21', '28', '42', '56'],
                        ant_days=['2', '4', '6', '7'])
pro_subject_8_int_dict = pro_subject_8.create_int_dict()
pro_subject_8_base_dict = pro_subject_8.create_base_dict()
pro_subject_8_ant_dict = pro_subject_8.create_ant_dict()

spo_subject_1 = Subject(base_array=spo_baseline_1, ant_array=spo_antibiotics_1, int_array=spo_intervention_1,
                        month_array=spo_month_1, base_days=['1', '2', '3', '4', '7'],
                        ant_days=['2', '3', '4', '5', '6', '7'],
                        int_days=['3', '4', '5', '7', '14', '21', '28', '42', '56'], month_number=['3', '4', '5', '6'])
spo_subject_1_base_dict = spo_subject_1.create_base_dict()
spo_subject_1_ant_dict = spo_subject_1.create_ant_dict()
spo_subject_1_int_dict = spo_subject_1.create_int_dict()
spo_subject_1_month_dict = spo_subject_1.create_month_dict()

spo_subject_2 = Subject(base_array=spo_baseline_2, ant_array=spo_antibiotics_2, int_array=spo_intervention_2,
                        month_array=spo_month_2, base_days=['1', '2', '3', '4', '5', '6', '7'],
                        ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        month_number=['3', '4', '5', '6'])
spo_subject_2_base_dict = spo_subject_2.create_base_dict()
spo_subject_2_ant_dict = spo_subject_2.create_ant_dict()
spo_subject_2_int_dict = spo_subject_2.create_int_dict()
spo_subject_2_month_dict = spo_subject_2.create_month_dict()

spo_subject_3 = Subject(base_array=spo_baseline_3, ant_array=spo_antibiotics_3, int_array=spo_intervention_3,
                        month_array=spo_month_3, base_days=['1', '2', '3', '4', '5', '6', '7'],
                        ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        month_number=['4', '5', '6'])
spo_subject_3_base_dict = spo_subject_3.create_base_dict()
spo_subject_3_ant_dict = spo_subject_3.create_ant_dict()
spo_subject_3_int_dict = spo_subject_3.create_int_dict()
spo_subject_3_month_dict = spo_subject_3.create_month_dict()

spo_subject_4 = Subject(base_array=spo_baseline_4, ant_array=spo_antibiotics_4, int_array=spo_intervention_4,
                        base_days=['1', '2', '3', '4', '6', '7'], ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'])
spo_subject_4_base_dict = spo_subject_4.create_base_dict()
spo_subject_4_ant_dict = spo_subject_4.create_ant_dict()
spo_subject_4_int_dict = spo_subject_4.create_int_dict()

spo_subject_5 = Subject(base_array=spo_baseline_5, ant_array=spo_antibiotics_5, int_array=spo_intervention_5,
                        month_array=spo_month_5, base_days=['1', '2', '3', '4', '5', '6', '7'],
                        ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'],
                        month_number=['4', '5', '6'])
spo_subject_5_base_dict = spo_subject_5.create_base_dict()
spo_subject_5_ant_dict = spo_subject_5.create_ant_dict()
spo_subject_5_int_dict = spo_subject_5.create_int_dict()
spo_subject_5_month_dict = spo_subject_5.create_month_dict()

spo_subject_6 = Subject(base_array=spo_baseline_6, ant_array=spo_antibiotics_6, int_array=spo_intervention_6,
                        base_days=['2', '3', '4', '6', '7'],
                        ant_days=['2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '6', '7', '14', '21', '28', '42', '56'])
spo_subject_6_base_dict = spo_subject_6.create_base_dict()
spo_subject_6_ant_dict = spo_subject_6.create_ant_dict()
spo_subject_6_int_dict = spo_subject_6.create_int_dict()

spo_subject_7 = Subject(base_array=spo_baseline_7, ant_array=spo_antibiotics_7, int_array=spo_intervention_7,
                        base_days=['1', '2', '3', '4', '5', '6', '7'],
                        ant_days=['1', '2', '3', '4', '5', '6', '7'],
                        int_days=['1', '2', '3', '4', '5', '6', '7', '14', '21', '28', '42', '56'])
spo_subject_7_base_dict = spo_subject_7.create_base_dict()
spo_subject_7_ant_dict = spo_subject_7.create_ant_dict()
spo_subject_7_int_dict = spo_subject_7.create_int_dict()

# Create cohorts
baseline_dict = {'aFMT_1': aFMT_baseline_1.T, 'aFMT_2': aFMT_baseline_2.T, 'aFMT_3': aFMT_baseline_3.T,
                 'aFMT_4': aFMT_baseline_4.T, 'aFMT_5': aFMT_baseline_5.T, 'aFMT_6': aFMT_baseline_6[:, 2:].T,
                 'pro_1': pro_baseline_1.T, 'pro_2': pro_baseline_2.T, 'pro_3': pro_baseline_3.T,
                 'pro_4': pro_baseline_4.T, 'pro_5': pro_baseline_5.T, 'pro_6': pro_baseline_6.T,
                 'pro_7': pro_baseline_7.T, 'pro_8': pro_baseline_8.T, 'spo_1': spo_baseline_1.T,
                 'spo_2': spo_baseline_2.T, 'spo_3': spo_baseline_3.T, 'spo_4': spo_baseline_4.T,
                 'spo_5': spo_baseline_5.T, 'spo_6': spo_baseline_6.T, 'spo_7': spo_baseline_7.T}

baseline_list = list(baseline_dict.values())
baseline_total_matrix = np.concatenate(baseline_list, axis=0)
baseline_total_matrix = normalize_cohort(baseline_total_matrix)

baseline_cohort_opt = OptimalCohort(baseline_dict)
baseline_cohort, chosen_indices = baseline_cohort_opt.get_optimal_samples()

baseline_cohort = normalize_cohort(baseline_cohort)

spo_1_samples = np.vstack([spo_subject_1_base_dict['1'], spo_subject_1_base_dict['2'],
                           spo_subject_1_base_dict['3'], spo_subject_1_base_dict['4'],
                           spo_subject_1_base_dict['7'], spo_subject_1_ant_dict['2'],
                           spo_subject_1_ant_dict['3'], spo_subject_1_ant_dict['4'],
                           spo_subject_1_ant_dict['5'], spo_subject_1_ant_dict['6'],
                           spo_subject_1_ant_dict['7'], spo_subject_1_int_dict['3'],
                           spo_subject_1_int_dict['4'], spo_subject_1_int_dict['5'],
                           spo_subject_1_int_dict['7'], spo_subject_1_int_dict['14'],
                           spo_subject_1_int_dict['21'], spo_subject_1_int_dict['28'],
                           spo_subject_1_int_dict['42'], spo_subject_1_int_dict['56'],
                           spo_subject_1_month_dict['3'], spo_subject_1_month_dict['4'],
                           spo_subject_1_month_dict['5'], spo_subject_1_month_dict['6']])

spo_1_samples = normalize_cohort(spo_1_samples)

spo_2_samples = np.vstack([spo_subject_2_base_dict['1'], spo_subject_2_base_dict['2'],
                           spo_subject_2_base_dict['3'], spo_subject_2_base_dict['4'],
                           spo_subject_2_base_dict['5'], spo_subject_2_base_dict['6'],
                           spo_subject_2_base_dict['7'], spo_subject_2_ant_dict['1'],
                           spo_subject_2_ant_dict['2'], spo_subject_2_ant_dict['3'],
                           spo_subject_2_ant_dict['4'], spo_subject_2_ant_dict['5'],
                           spo_subject_2_ant_dict['6'], spo_subject_2_ant_dict['7'],
                           spo_subject_2_int_dict['1'], spo_subject_2_int_dict['2'],
                           spo_subject_2_int_dict['3'], spo_subject_2_int_dict['4'],
                           spo_subject_2_int_dict['5'], spo_subject_2_int_dict['6'],
                           spo_subject_2_int_dict['7'], spo_subject_2_int_dict['14'],
                           spo_subject_2_int_dict['21'], spo_subject_2_int_dict['28'],
                           spo_subject_2_int_dict['42'], spo_subject_2_int_dict['56'],
                           spo_subject_2_month_dict['3'], spo_subject_2_month_dict['4'],
                           spo_subject_2_month_dict['5'], spo_subject_2_month_dict['6']])

spo_2_samples = normalize_cohort(spo_2_samples)

spo_3_samples = np.vstack([spo_subject_3_base_dict['1'], spo_subject_3_base_dict['2'],
                           spo_subject_3_base_dict['3'], spo_subject_3_base_dict['4'],
                           spo_subject_3_base_dict['5'], spo_subject_3_base_dict['6'],
                           spo_subject_3_base_dict['7'], spo_subject_3_ant_dict['1'],
                           spo_subject_3_ant_dict['2'], spo_subject_3_ant_dict['3'],
                           spo_subject_3_ant_dict['4'], spo_subject_3_ant_dict['5'],
                           spo_subject_3_ant_dict['6'], spo_subject_3_ant_dict['7'],
                           spo_subject_3_int_dict['1'], spo_subject_3_int_dict['2'],
                           spo_subject_3_int_dict['3'], spo_subject_3_int_dict['4'],
                           spo_subject_3_int_dict['5'], spo_subject_3_int_dict['6'],
                           spo_subject_3_int_dict['7'], spo_subject_3_int_dict['14'],
                           spo_subject_3_int_dict['21'], spo_subject_3_int_dict['28'],
                           spo_subject_3_int_dict['42'], spo_subject_3_int_dict['56'],
                           spo_subject_3_month_dict['4'], spo_subject_3_month_dict['5'],
                           spo_subject_3_month_dict['6']])

spo_3_samples = normalize_cohort(spo_3_samples)

spo_4_samples = np.vstack([spo_subject_4_base_dict['1'], spo_subject_4_base_dict['2'],
                           spo_subject_4_base_dict['3'], spo_subject_4_base_dict['4'],
                           spo_subject_4_base_dict['6'],
                           spo_subject_4_base_dict['7'], spo_subject_4_ant_dict['1'],
                           spo_subject_4_ant_dict['2'], spo_subject_4_ant_dict['3'],
                           spo_subject_4_ant_dict['4'], spo_subject_4_ant_dict['5'],
                           spo_subject_4_ant_dict['6'], spo_subject_4_ant_dict['7'],
                           spo_subject_4_int_dict['1'], spo_subject_4_int_dict['2'],
                           spo_subject_4_int_dict['3'], spo_subject_4_int_dict['4'],
                           spo_subject_4_int_dict['5'], spo_subject_4_int_dict['6'],
                           spo_subject_4_int_dict['7'], spo_subject_4_int_dict['14'],
                           spo_subject_4_int_dict['21'], spo_subject_4_int_dict['28'],
                           spo_subject_4_int_dict['42'], spo_subject_4_int_dict['56']])

spo_4_samples = normalize_cohort(spo_4_samples)

spo_5_samples = np.vstack([spo_subject_5_base_dict['1'], spo_subject_5_base_dict['2'],
                           spo_subject_5_base_dict['3'], spo_subject_5_base_dict['4'],
                           spo_subject_5_base_dict['5'], spo_subject_5_base_dict['6'],
                           spo_subject_5_base_dict['7'], spo_subject_5_ant_dict['1'],
                           spo_subject_5_ant_dict['2'], spo_subject_5_ant_dict['3'],
                           spo_subject_5_ant_dict['4'], spo_subject_5_ant_dict['5'],
                           spo_subject_5_ant_dict['6'], spo_subject_5_ant_dict['7'],
                           spo_subject_5_int_dict['1'], spo_subject_5_int_dict['2'],
                           spo_subject_5_int_dict['3'], spo_subject_5_int_dict['4'],
                           spo_subject_5_int_dict['5'], spo_subject_5_int_dict['6'],
                           spo_subject_5_int_dict['7'], spo_subject_5_int_dict['14'],
                           spo_subject_5_int_dict['21'], spo_subject_5_int_dict['28'],
                           spo_subject_5_int_dict['42'], spo_subject_5_int_dict['56'],
                           spo_subject_5_month_dict['4'], spo_subject_5_month_dict['5'],
                           spo_subject_5_month_dict['6']])

spo_5_samples = normalize_cohort(spo_5_samples)

spo_6_samples = np.vstack([spo_subject_6_base_dict['2'], spo_subject_6_base_dict['3'],
                           spo_subject_6_base_dict['4'], spo_subject_6_base_dict['6'],
                           spo_subject_6_base_dict['7'], spo_subject_6_ant_dict['2'],
                           spo_subject_6_ant_dict['3'], spo_subject_6_ant_dict['4'],
                           spo_subject_6_ant_dict['5'], spo_subject_6_ant_dict['6'],
                           spo_subject_6_ant_dict['7'], spo_subject_6_int_dict['1'],
                           spo_subject_6_int_dict['2'], spo_subject_6_int_dict['3'],
                           spo_subject_6_int_dict['4'], spo_subject_6_int_dict['6'],
                           spo_subject_6_int_dict['7'], spo_subject_6_int_dict['14'],
                           spo_subject_6_int_dict['21'], spo_subject_6_int_dict['28'],
                           spo_subject_6_int_dict['42'], spo_subject_6_int_dict['56']])

spo_6_samples = normalize_cohort(spo_6_samples)

spo_7_samples = np.vstack([spo_subject_7_base_dict['1'], spo_subject_7_base_dict['2'],
                           spo_subject_7_base_dict['3'], spo_subject_7_base_dict['4'],
                           spo_subject_7_base_dict['5'], spo_subject_7_base_dict['6'],
                           spo_subject_7_base_dict['7'], spo_subject_7_ant_dict['1'],
                           spo_subject_7_ant_dict['2'], spo_subject_7_ant_dict['3'],
                           spo_subject_7_ant_dict['4'], spo_subject_7_ant_dict['5'],
                           spo_subject_7_ant_dict['6'], spo_subject_7_ant_dict['7'],
                           spo_subject_7_int_dict['2'], spo_subject_7_int_dict['3'],
                           spo_subject_7_int_dict['4'], spo_subject_7_int_dict['5'],
                           spo_subject_7_int_dict['6'], spo_subject_7_int_dict['7'],
                           spo_subject_7_int_dict['14'], spo_subject_7_int_dict['21'],
                           spo_subject_7_int_dict['28'], spo_subject_7_int_dict['42'],
                           spo_subject_7_int_dict['56']])

spo_7_samples = normalize_cohort(spo_7_samples)

spo_samples_list = [spo_1_samples, spo_2_samples, spo_3_samples, spo_4_samples,
                    spo_5_samples, spo_6_samples, spo_7_samples]

aFMT_1_baseline = np.vstack([aFMT_subject_1_base_dict['1'], aFMT_subject_1_base_dict['2'],
                             aFMT_subject_1_base_dict['3'], aFMT_subject_1_base_dict['4'],
                             aFMT_subject_1_base_dict['5'], aFMT_subject_1_base_dict['6'],
                             aFMT_subject_1_base_dict['7']])

aFMT_1_baseline = normalize_cohort(aFMT_1_baseline)

aFMT_1_samples = np.vstack([aFMT_subject_1_base_dict['1'], aFMT_subject_1_base_dict['2'],
                            aFMT_subject_1_base_dict['3'], aFMT_subject_1_base_dict['4'],
                            aFMT_subject_1_base_dict['5'], aFMT_subject_1_base_dict['6'],
                            aFMT_subject_1_base_dict['7'], aFMT_subject_1_ant_dict['2'],
                            aFMT_subject_1_int_dict['1'],  aFMT_subject_1_int_dict['2'],
                            aFMT_subject_1_int_dict['3'],  aFMT_subject_1_int_dict['4'],
                            aFMT_subject_1_int_dict['5'],  aFMT_subject_1_int_dict['14'],
                            aFMT_subject_1_int_dict['21'], aFMT_subject_1_int_dict['28'],
                            aFMT_subject_1_int_dict['42'], aFMT_subject_1_int_dict['56'],
                            aFMT_subject_1_month_dict['1'], aFMT_subject_1_month_dict['2'],
                            aFMT_subject_1_month_dict['3'], aFMT_subject_1_month_dict['4'],
                            aFMT_subject_1_month_dict['6']])

aFMT_1_samples = normalize_cohort(aFMT_1_samples)

aFMT_2_baseline = np.vstack([aFMT_subject_2_base_dict['1'], aFMT_subject_2_base_dict['2'],
                             aFMT_subject_2_base_dict['3'], aFMT_subject_2_base_dict['5'],
                             aFMT_subject_2_base_dict['6'], aFMT_subject_2_base_dict['7']])

aFMT_2_baseline = normalize_cohort(aFMT_2_baseline)

aFMT_2_samples = np.vstack([aFMT_subject_2_base_dict['1'], aFMT_subject_2_base_dict['2'],
                            aFMT_subject_2_base_dict['3'], aFMT_subject_2_base_dict['5'],
                            aFMT_subject_2_base_dict['6'], aFMT_subject_2_base_dict['7'],
                            aFMT_subject_2_ant_dict['2'], aFMT_subject_2_ant_dict['4'],
                            aFMT_subject_2_ant_dict['5'], aFMT_subject_2_ant_dict['6'],
                            aFMT_subject_2_ant_dict['7'], aFMT_subject_2_int_dict['1'],
                            aFMT_subject_2_int_dict['2'], aFMT_subject_2_int_dict['3'],
                            aFMT_subject_2_int_dict['4'], aFMT_subject_2_int_dict['5'],
                            aFMT_subject_2_int_dict['6'], aFMT_subject_2_int_dict['7'],
                            aFMT_subject_2_int_dict['14'], aFMT_subject_2_int_dict['21'],
                            aFMT_subject_2_int_dict['28'], aFMT_subject_2_int_dict['42'],
                            aFMT_subject_2_int_dict['56'], aFMT_subject_2_month_dict['2'],
                            aFMT_subject_2_month_dict['3'], aFMT_subject_2_month_dict['4'],
                            aFMT_subject_2_month_dict['5'], aFMT_subject_2_month_dict['6']])

aFMT_2_samples = normalize_cohort(aFMT_2_samples)

aFMT_3_baseline = np.vstack([aFMT_subject_3_base_dict['1'], aFMT_subject_3_base_dict['2'],
                             aFMT_subject_3_base_dict['3'], aFMT_subject_3_base_dict['4'],
                             aFMT_subject_3_base_dict['5'], aFMT_subject_3_base_dict['7']])

aFMT_3_baseline = normalize_cohort(aFMT_3_baseline)

aFMT_3_samples = np.vstack([aFMT_subject_3_base_dict['1'], aFMT_subject_3_base_dict['2'],
                            aFMT_subject_3_base_dict['3'], aFMT_subject_3_base_dict['4'],
                            aFMT_subject_3_base_dict['5'], aFMT_subject_3_base_dict['7'],
                            aFMT_subject_3_ant_dict['2'], aFMT_subject_3_ant_dict['3'],
                            aFMT_subject_3_ant_dict['4'], aFMT_subject_3_ant_dict['5'],
                            aFMT_subject_3_ant_dict['6'], aFMT_subject_3_ant_dict['7'],
                            aFMT_subject_3_int_dict['1'], aFMT_subject_3_int_dict['2'],
                            aFMT_subject_3_int_dict['3'], aFMT_subject_3_int_dict['4'],
                            aFMT_subject_3_int_dict['5'], aFMT_subject_3_int_dict['6'],
                            aFMT_subject_3_int_dict['7'], aFMT_subject_3_int_dict['14'],
                            aFMT_subject_3_int_dict['21'], aFMT_subject_3_int_dict['28'],
                            aFMT_subject_3_int_dict['42'], aFMT_subject_3_int_dict['56'],
                            aFMT_subject_3_month_dict['2'], aFMT_subject_3_month_dict['3'],
                            aFMT_subject_3_month_dict['4'], aFMT_subject_3_month_dict['5']])

aFMT_3_samples = normalize_cohort(aFMT_3_samples)

aFMT_4_baseline = np.vstack([aFMT_subject_4_base_dict['1'], aFMT_subject_4_base_dict['3'],
                             aFMT_subject_4_base_dict['4'], aFMT_subject_4_base_dict['5'],
                             aFMT_subject_4_base_dict['6'], aFMT_subject_4_base_dict['7']])

aFMT_4_baseline = normalize_cohort(aFMT_4_baseline)

aFMT_4_samples = np.vstack([aFMT_subject_4_base_dict['1'], aFMT_subject_4_base_dict['3'],
                            aFMT_subject_4_base_dict['4'], aFMT_subject_4_base_dict['5'],
                            aFMT_subject_4_base_dict['6'], aFMT_subject_4_base_dict['7'],
                            aFMT_subject_4_ant_dict['1'], aFMT_subject_4_ant_dict['2'],
                            aFMT_subject_4_ant_dict['5'], aFMT_subject_4_ant_dict['6'],
                            aFMT_subject_4_ant_dict['7'], aFMT_subject_4_int_dict['3'],
                            aFMT_subject_4_int_dict['4'], aFMT_subject_4_int_dict['5'],
                            aFMT_subject_4_int_dict['6'], aFMT_subject_4_int_dict['7'],
                            aFMT_subject_4_int_dict['14'], aFMT_subject_4_int_dict['21'],
                            aFMT_subject_4_int_dict['28'], aFMT_subject_4_int_dict['42'],
                            aFMT_subject_4_int_dict['56'], aFMT_subject_4_month_dict['3'],
                            aFMT_subject_4_month_dict['4'], aFMT_subject_4_month_dict['5'],
                            aFMT_subject_4_month_dict['6']])

aFMT_4_samples = normalize_cohort(aFMT_4_samples)

aFMT_5_baseline = np.vstack([aFMT_subject_5_base_dict['1'], aFMT_subject_5_base_dict['2'],
                             aFMT_subject_5_base_dict['3'], aFMT_subject_5_base_dict['4'],
                             aFMT_subject_5_base_dict['5'], aFMT_subject_5_base_dict['6'],
                             aFMT_subject_5_base_dict['7']])

aFMT_5_baseline = normalize_cohort(aFMT_5_baseline)

aFMT_5_samples = np.vstack([aFMT_subject_5_base_dict['1'], aFMT_subject_5_base_dict['2'],
                            aFMT_subject_5_base_dict['3'], aFMT_subject_5_base_dict['4'],
                            aFMT_subject_5_base_dict['5'], aFMT_subject_5_base_dict['6'],
                            aFMT_subject_5_base_dict['7'], aFMT_subject_5_ant_dict['1'],
                            aFMT_subject_5_ant_dict['2'], aFMT_subject_5_ant_dict['3'],
                            aFMT_subject_5_ant_dict['4'], aFMT_subject_5_ant_dict['5'],
                            aFMT_subject_5_ant_dict['6'], aFMT_subject_5_int_dict['1'],
                            aFMT_subject_5_int_dict['2'], aFMT_subject_5_int_dict['3'],
                            aFMT_subject_5_int_dict['4'], aFMT_subject_5_int_dict['5'],
                            aFMT_subject_5_int_dict['6'], aFMT_subject_5_int_dict['7'],
                            aFMT_subject_5_int_dict['14'], aFMT_subject_5_int_dict['21'],
                            aFMT_subject_5_int_dict['28'], aFMT_subject_5_int_dict['42'],
                            aFMT_subject_5_int_dict['56'], aFMT_subject_5_month_dict['4'],
                            aFMT_subject_5_month_dict['5'], aFMT_subject_5_month_dict['6']])

aFMT_5_samples = normalize_cohort(aFMT_5_samples)

aFMT_6_baseline = np.vstack([aFMT_subject_6_base_dict['3'], aFMT_subject_6_base_dict['4'],
                             aFMT_subject_6_base_dict['5'], aFMT_subject_6_base_dict['6'],
                             aFMT_subject_6_base_dict['7']])

aFMT_6_baseline = normalize_cohort(aFMT_6_baseline)

aFMT_6_samples = np.vstack([aFMT_subject_6_base_dict['3'], aFMT_subject_6_base_dict['4'],
                            aFMT_subject_6_base_dict['5'], aFMT_subject_6_base_dict['6'],
                            aFMT_subject_6_base_dict['7'],
                            aFMT_subject_6_ant_dict['2'], aFMT_subject_6_ant_dict['3'],
                            aFMT_subject_6_ant_dict['4'], aFMT_subject_6_ant_dict['5'],
                            aFMT_subject_6_ant_dict['6'], aFMT_subject_6_int_dict['1'],
                            aFMT_subject_6_int_dict['2'], aFMT_subject_6_int_dict['3'],
                            aFMT_subject_6_int_dict['4'], aFMT_subject_6_int_dict['5'],
                            aFMT_subject_6_int_dict['6'], aFMT_subject_6_int_dict['7'],
                            aFMT_subject_6_int_dict['14'], aFMT_subject_6_int_dict['21'],
                            aFMT_subject_6_int_dict['28'], aFMT_subject_6_int_dict['42'],
                            aFMT_subject_6_int_dict['56']])

aFMT_6_samples = normalize_cohort(aFMT_6_samples)

pro_1_baseline = np.vstack([pro_subject_1_base_dict['1'], pro_subject_1_base_dict['2'],
                            pro_subject_1_base_dict['3'], pro_subject_1_base_dict['4'],
                            pro_subject_1_base_dict['5'], pro_subject_1_base_dict['6'],
                            pro_subject_1_base_dict['7']])

pro_1_baseline = normalize_cohort(pro_1_baseline)

pro_1_samples = np.vstack([pro_subject_1_base_dict['1'], pro_subject_1_base_dict['2'], pro_subject_1_base_dict['3'],
                           pro_subject_1_base_dict['4'], pro_subject_1_base_dict['5'], pro_subject_1_base_dict['6'],
                           pro_subject_1_base_dict['7'], pro_subject_1_ant_dict['1'], pro_subject_1_ant_dict['2'],
                           pro_subject_1_ant_dict['3'], pro_subject_1_ant_dict['4'], pro_subject_1_ant_dict['5'],
                           pro_subject_1_ant_dict['6'], pro_subject_1_ant_dict['7'], pro_subject_1_int_dict['1'],
                           pro_subject_1_int_dict['2'], pro_subject_1_int_dict['3'], pro_subject_1_int_dict['4'],
                           pro_subject_1_int_dict['5'], pro_subject_1_int_dict['6'], pro_subject_1_int_dict['7'],
                           pro_subject_1_int_dict['14'], pro_subject_1_int_dict['21'], pro_subject_1_int_dict['28'],
                           pro_subject_1_int_dict['42'], pro_subject_1_int_dict['56'], pro_subject_1_month_dict['2'],
                           pro_subject_1_month_dict['3'], pro_subject_1_month_dict['4'], pro_subject_1_month_dict['5'],
                           pro_subject_1_month_dict['6']])

pro_1_samples = normalize_cohort(pro_1_samples)

pro_2_baseline = np.vstack([pro_subject_2_base_dict['1'], pro_subject_2_base_dict['2'],
                            pro_subject_2_base_dict['3'], pro_subject_2_base_dict['4'],
                            pro_subject_2_base_dict['5']])

pro_2_baseline = normalize_cohort(pro_2_baseline)

pro_2_samples = np.vstack([pro_subject_2_base_dict['1'], pro_subject_2_base_dict['2'], pro_subject_2_base_dict['3'],
                           pro_subject_2_base_dict['4'], pro_subject_2_base_dict['5'], pro_subject_2_ant_dict['2'],
                           pro_subject_2_ant_dict['3'], pro_subject_2_ant_dict['4'], pro_subject_2_ant_dict['5'],
                           pro_subject_2_ant_dict['6'], pro_subject_2_ant_dict['7'], pro_subject_2_int_dict['1'],
                           pro_subject_2_int_dict['2'], pro_subject_2_int_dict['4'], pro_subject_2_int_dict['5'],
                           pro_subject_2_int_dict['6'], pro_subject_2_int_dict['14'], pro_subject_2_int_dict['21'],
                           pro_subject_2_month_dict['2'], pro_subject_2_month_dict['3'], pro_subject_2_month_dict['4'],
                           pro_subject_2_month_dict['5']])

pro_2_samples = normalize_cohort(pro_2_samples)

pro_3_baseline = np.vstack([pro_subject_3_base_dict['1'], pro_subject_3_base_dict['2'],
                            pro_subject_3_base_dict['3'], pro_subject_3_base_dict['4'],
                            pro_subject_3_base_dict['5'], pro_subject_3_base_dict['7']])

pro_3_baseline = normalize_cohort(pro_3_baseline)

pro_3_samples = np.vstack([pro_subject_3_base_dict['1'], pro_subject_3_base_dict['2'], pro_subject_3_base_dict['3'],
                           pro_subject_3_base_dict['4'], pro_subject_3_base_dict['5'], pro_subject_3_base_dict['7'],
                           pro_subject_3_ant_dict['1'], pro_subject_3_ant_dict['2'], pro_subject_3_ant_dict['3'],
                           pro_subject_3_ant_dict['4'], pro_subject_3_ant_dict['6'], pro_subject_3_int_dict['1'],
                           pro_subject_3_int_dict['2'], pro_subject_3_int_dict['3'], pro_subject_3_int_dict['4'],
                           pro_subject_3_int_dict['5'], pro_subject_3_int_dict['6'], pro_subject_3_int_dict['14'],
                           pro_subject_3_int_dict['21'], pro_subject_3_int_dict['28'], pro_subject_3_int_dict['42'],
                           pro_subject_3_int_dict['56'], pro_subject_3_month_dict['3'], pro_subject_3_month_dict['5'],
                           pro_subject_3_month_dict['6']])

pro_3_samples = normalize_cohort(pro_3_samples)

pro_4_baseline = np.vstack([pro_subject_4_base_dict['1'], pro_subject_4_base_dict['2'],
                            pro_subject_4_base_dict['3'], pro_subject_4_base_dict['4'],
                            pro_subject_4_base_dict['5'], pro_subject_4_base_dict['6'],
                            pro_subject_4_base_dict['7']])

pro_4_baseline = normalize_cohort(pro_4_baseline)

pro_4_samples = np.vstack([pro_subject_4_base_dict['1'], pro_subject_4_base_dict['2'], pro_subject_4_base_dict['3'],
                           pro_subject_4_base_dict['4'], pro_subject_4_base_dict['5'], pro_subject_4_base_dict['6'],
                           pro_subject_4_base_dict['7'], pro_subject_4_ant_dict['1'], pro_subject_4_ant_dict['2'],
                           pro_subject_4_ant_dict['3'], pro_subject_4_ant_dict['4'], pro_subject_4_ant_dict['5'],
                           pro_subject_4_ant_dict['6'], pro_subject_4_ant_dict['7'], pro_subject_4_int_dict['1'],
                           pro_subject_4_int_dict['2'], pro_subject_4_int_dict['4'], pro_subject_4_int_dict['5'],
                           pro_subject_4_int_dict['6'], pro_subject_4_int_dict['14'], pro_subject_4_int_dict['21'],
                           pro_subject_4_int_dict['28'], pro_subject_4_int_dict['42'], pro_subject_4_int_dict['56'],
                           pro_subject_4_month_dict['3'], pro_subject_4_month_dict['4'], pro_subject_4_month_dict['5'],
                           pro_subject_4_month_dict['6']])

pro_4_samples = normalize_cohort(pro_4_samples)

pro_5_baseline = np.vstack([pro_subject_5_base_dict['1'], pro_subject_5_base_dict['2'],
                            pro_subject_5_base_dict['3'], pro_subject_5_base_dict['4'],
                            pro_subject_5_base_dict['5'], pro_subject_5_base_dict['6'],
                            pro_subject_5_base_dict['7']])

pro_5_baseline = normalize_cohort(pro_5_baseline)

pro_5_samples = np.vstack([pro_subject_5_base_dict['1'], pro_subject_5_base_dict['2'], pro_subject_5_base_dict['3'],
                           pro_subject_5_base_dict['4'], pro_subject_5_base_dict['5'], pro_subject_5_base_dict['6'],
                           pro_subject_5_base_dict['7'], pro_subject_5_ant_dict['1'], pro_subject_5_ant_dict['2'],
                           pro_subject_5_int_dict['1'], pro_subject_5_int_dict['2'], pro_subject_5_int_dict['3'],
                           pro_subject_5_int_dict['4'], pro_subject_5_int_dict['5'], pro_subject_5_int_dict['6'],
                           pro_subject_5_int_dict['7'], pro_subject_5_int_dict['14'], pro_subject_5_int_dict['21'],
                           pro_subject_5_int_dict['28'], pro_subject_5_int_dict['42'], pro_subject_5_int_dict['56'],
                           pro_subject_5_month_dict['4'], pro_subject_5_month_dict['5'], pro_subject_5_month_dict['6']])

pro_5_samples = normalize_cohort(pro_5_samples)

pro_6_baseline = np.vstack([pro_subject_6_base_dict['1'], pro_subject_6_base_dict['2'],
                            pro_subject_6_base_dict['3'], pro_subject_6_base_dict['4'],
                            pro_subject_6_base_dict['5'], pro_subject_6_base_dict['6'],
                            pro_subject_6_base_dict['7']])

pro_6_baseline = normalize_cohort(pro_6_baseline)

pro_6_samples = np.vstack([pro_subject_6_base_dict['1'], pro_subject_6_base_dict['2'], pro_subject_6_base_dict['3'],
                           pro_subject_6_base_dict['4'], pro_subject_6_base_dict['5'], pro_subject_6_base_dict['6'],
                           pro_subject_6_base_dict['7'], pro_subject_6_ant_dict['1'], pro_subject_6_ant_dict['2'],
                           pro_subject_6_ant_dict['4'], pro_subject_6_ant_dict['5'], pro_subject_6_ant_dict['6'],
                           pro_subject_6_ant_dict['7'], pro_subject_6_int_dict['1'], pro_subject_6_int_dict['2'],
                           pro_subject_6_int_dict['3'], pro_subject_6_int_dict['4'], pro_subject_6_int_dict['5'],
                           pro_subject_6_int_dict['6'], pro_subject_6_int_dict['7'], pro_subject_6_int_dict['14'],
                           pro_subject_6_int_dict['21'], pro_subject_6_int_dict['28'], pro_subject_6_int_dict['42'],
                           pro_subject_6_int_dict['56']])

pro_6_samples = normalize_cohort(pro_6_samples)

pro_7_baseline = np.vstack([pro_subject_7_base_dict['1'], pro_subject_7_base_dict['2'],
                            pro_subject_7_base_dict['3'], pro_subject_7_base_dict['4'],
                            pro_subject_7_base_dict['5'], pro_subject_7_base_dict['6'],
                            pro_subject_7_base_dict['7']])

pro_7_baseline = normalize_cohort(pro_7_baseline)

pro_7_samples = np.vstack([pro_subject_7_base_dict['1'], pro_subject_7_base_dict['2'], pro_subject_7_base_dict['3'],
                           pro_subject_7_base_dict['4'], pro_subject_7_base_dict['5'], pro_subject_7_base_dict['6'],
                           pro_subject_7_base_dict['7'], pro_subject_7_ant_dict['1'], pro_subject_7_ant_dict['2'],
                           pro_subject_7_ant_dict['3'], pro_subject_7_ant_dict['4'], pro_subject_7_ant_dict['5'],
                           pro_subject_7_ant_dict['6'], pro_subject_7_ant_dict['7'], pro_subject_7_int_dict['1'],
                           pro_subject_7_int_dict['2'], pro_subject_7_int_dict['3'], pro_subject_7_int_dict['4'],
                           pro_subject_7_int_dict['5'], pro_subject_7_int_dict['6'], pro_subject_7_int_dict['7'],
                           pro_subject_7_int_dict['14'], pro_subject_7_int_dict['21'], pro_subject_7_int_dict['28'],
                           pro_subject_7_int_dict['42'], pro_subject_7_int_dict['56']])

pro_7_samples = normalize_cohort(pro_7_samples)

pro_8_baseline = np.vstack([pro_subject_8_base_dict['1'], pro_subject_8_base_dict['2'],
                            pro_subject_8_base_dict['3']])

pro_8_baseline = normalize_cohort(pro_8_baseline)

pro_8_samples = np.vstack([pro_subject_8_base_dict['1'], pro_subject_8_base_dict['2'], pro_subject_8_base_dict['3'],
                           pro_subject_8_ant_dict['2'], pro_subject_8_ant_dict['4'], pro_subject_8_ant_dict['6'],
                           pro_subject_8_ant_dict['7'], pro_subject_8_int_dict['3'], pro_subject_8_int_dict['5'],
                           pro_subject_8_int_dict['6'], pro_subject_8_int_dict['7'], pro_subject_8_int_dict['14'],
                           pro_subject_8_int_dict['21'], pro_subject_8_int_dict['28'], pro_subject_8_int_dict['42'],
                           pro_subject_8_int_dict['56']])

pro_8_samples = normalize_cohort(pro_8_samples)

ABX_cohort = np.vstack([aFMT_subject_1_ant_dict['2'], aFMT_subject_2_ant_dict['7'], aFMT_subject_3_ant_dict['7'],
                        aFMT_subject_4_ant_dict['7'], aFMT_subject_5_ant_dict['6'], aFMT_subject_6_ant_dict['7'],
                        pro_subject_1_ant_dict['7'], pro_subject_2_ant_dict['7'], pro_subject_3_ant_dict['6'],
                        pro_subject_4_ant_dict['7'], pro_subject_5_ant_dict['2'], pro_subject_6_ant_dict['7'],
                        pro_subject_7_ant_dict['7'], pro_subject_8_ant_dict['7'], spo_subject_1_ant_dict['7'],
                        spo_subject_2_ant_dict['7'], spo_subject_3_ant_dict['7'], spo_subject_4_ant_dict['7'],
                        spo_subject_5_ant_dict['7'], spo_subject_6_ant_dict['7'], spo_subject_7_ant_dict['7']])

ABX_cohort = normalize_cohort(ABX_cohort)

#ABX_dict = {'aFMT_1': aFMT_antibiotics_1.reshape(-1, 1).T, 'aFMT_2': aFMT_antibiotics_2.T, 'aFMT_3': aFMT_antibiotics_3.T,
#            'aFMT_4': aFMT_antibiotics_4.T, 'aFMT_5': aFMT_antibiotics_5.T, 'aFMT_6': aFMT_antibiotics_6.T,
#            'pro_1': pro_antibiotics_1.T, 'pro_2': pro_antibiotics_2.T, 'pro_3': pro_antibiotics_3.T,
#            'pro_4': pro_antibiotics_4.T, 'pro_5': pro_antibiotics_5.T, 'pro_6': pro_antibiotics_6.T,
#            'pro_7': pro_antibiotics_7.T, 'pro_8': pro_antibiotics_8.T, 'spo_1': spo_antibiotics_1.T,
#            'spo_2': spo_antibiotics_2.T, 'spo_3': spo_antibiotics_3.T, 'spo_4': spo_antibiotics_4.T,
#            'spo_5': spo_antibiotics_5.T, 'spo_6': spo_antibiotics_6.T, 'spo_7': spo_antibiotics_7.T}

spo_1_baseline = np.vstack([spo_subject_1_base_dict['1'], spo_subject_1_base_dict['2'],
                            spo_subject_1_base_dict['3'], spo_subject_1_base_dict['4'],
                            spo_subject_1_base_dict['7']])

spo_2_baseline = np.vstack([spo_subject_2_base_dict['1'], spo_subject_2_base_dict['2'],
                           spo_subject_2_base_dict['3'], spo_subject_2_base_dict['4'],
                           spo_subject_2_base_dict['5'], spo_subject_2_base_dict['6'],
                           spo_subject_2_base_dict['7']])

spo_3_baseline = np.vstack([spo_subject_3_base_dict['1'], spo_subject_3_base_dict['2'],
                            spo_subject_3_base_dict['3'], spo_subject_3_base_dict['4'],
                            spo_subject_3_base_dict['5'], spo_subject_3_base_dict['6'],
                            spo_subject_3_base_dict['7']])

spo_4_baseline = np.vstack([spo_subject_4_base_dict['1'], spo_subject_4_base_dict['2'],
                            spo_subject_4_base_dict['3'], spo_subject_4_base_dict['4'],
                            spo_subject_4_base_dict['6'],
                            spo_subject_4_base_dict['7']])

spo_5_baseline = np.vstack([spo_subject_5_base_dict['1'], spo_subject_5_base_dict['2'],
                           spo_subject_5_base_dict['3'], spo_subject_5_base_dict['4'],
                           spo_subject_5_base_dict['5'], spo_subject_5_base_dict['6'],
                           spo_subject_5_base_dict['7']])

spo_6_baseline = np.vstack([spo_subject_6_base_dict['2'], spo_subject_6_base_dict['3'],
                            spo_subject_6_base_dict['4'], spo_subject_6_base_dict['6'],
                            spo_subject_6_base_dict['7']])

spo_7_baseline = np.vstack([spo_subject_7_base_dict['1'], spo_subject_7_base_dict['2'],
                            spo_subject_7_base_dict['3'], spo_subject_7_base_dict['4'],
                            spo_subject_7_base_dict['5'], spo_subject_7_base_dict['6'],
                            spo_subject_7_base_dict['7']])

spo_1_baseline = normalize_cohort(spo_1_baseline)
spo_2_baseline = normalize_cohort(spo_2_baseline)
spo_3_baseline = normalize_cohort(spo_3_baseline)
spo_4_baseline = normalize_cohort(spo_4_baseline)
spo_5_baseline = normalize_cohort(spo_5_baseline)
spo_6_baseline = normalize_cohort(spo_6_baseline)
spo_7_baseline = normalize_cohort(spo_7_baseline)

baseline_samples_list = [aFMT_1_baseline, aFMT_2_baseline, aFMT_3_baseline, aFMT_4_baseline, aFMT_5_baseline,
                         aFMT_6_baseline, pro_1_baseline, pro_2_baseline, pro_3_baseline, pro_4_baseline,
                         pro_5_baseline, pro_6_baseline, pro_7_baseline, pro_8_baseline, spo_1_baseline,
                         spo_2_baseline, spo_3_baseline, spo_4_baseline, spo_5_baseline, spo_6_baseline,
                         spo_7_baseline]

base_sam_spo_list = [spo_1_baseline, spo_2_baseline, spo_3_baseline, spo_4_baseline,
                     spo_5_baseline, spo_6_baseline, spo_7_baseline]

total_spo_samples = [spo_1_samples, spo_2_samples, spo_3_samples, spo_4_samples, spo_5_samples,
                     spo_6_samples, spo_7_samples]

spo_1_D_after_ABX = np.vstack([spo_subject_1_int_dict['42'], spo_subject_1_int_dict['56'],
                               spo_subject_1_month_dict['3'], spo_subject_1_month_dict['4'],
                               spo_subject_1_month_dict['5'], spo_subject_1_month_dict['6']])

spo_1_D_after_ABX = normalize_cohort(spo_1_D_after_ABX)

spo_2_D_after_ABX = np.vstack([spo_subject_2_int_dict['14'], spo_subject_2_int_dict['21'],
                               spo_subject_2_int_dict['28'], spo_subject_2_int_dict['42'],
                               spo_subject_2_int_dict['56'], spo_subject_2_month_dict['3'],
                               spo_subject_2_month_dict['4'], spo_subject_2_month_dict['5'],
                               spo_subject_2_month_dict['6']])

spo_2_D_after_ABX = normalize_cohort(spo_2_D_after_ABX)

spo_3_D_after_ABX = np.vstack([spo_subject_3_int_dict['14'], spo_subject_3_int_dict['21'],
                               spo_subject_3_int_dict['28'], spo_subject_3_int_dict['42'],
                               spo_subject_3_int_dict['56'], spo_subject_3_month_dict['4'],
                               spo_subject_3_month_dict['5'], spo_subject_3_month_dict['6']])

spo_3_D_after_ABX = normalize_cohort(spo_3_D_after_ABX)

spo_4_D_after_ABX = np.vstack([spo_subject_4_int_dict['14'], spo_subject_4_int_dict['21'],
                               spo_subject_4_int_dict['28'], spo_subject_4_int_dict['42'],
                               spo_subject_4_int_dict['56']])

spo_4_D_after_ABX = normalize_cohort(spo_4_D_after_ABX)

spo_5_D_after_ABX = np.vstack([spo_subject_5_int_dict['14'], spo_subject_5_int_dict['21'],
                               spo_subject_5_int_dict['28'], spo_subject_5_int_dict['42'],
                               spo_subject_5_int_dict['56'], spo_subject_5_month_dict['4'],
                               spo_subject_5_month_dict['5'], spo_subject_5_month_dict['6']])

spo_5_D_after_ABX = normalize_cohort(spo_5_D_after_ABX)

spo_6_D_after_ABX = np.vstack([spo_subject_6_int_dict['14'], spo_subject_6_int_dict['21'],
                               spo_subject_6_int_dict['28'], spo_subject_6_int_dict['42'],
                               spo_subject_6_int_dict['56']])

spo_6_D_after_ABX = normalize_cohort(spo_6_D_after_ABX)

spo_7_D_after_ABX = np.vstack([spo_subject_7_int_dict['14'], spo_subject_7_int_dict['21'],
                               spo_subject_7_int_dict['28'], spo_subject_7_int_dict['42'],
                               spo_subject_7_int_dict['56']])

spo_7_D_after_ABX = normalize_cohort(spo_7_D_after_ABX)

spo_days_after_ABX_list = [spo_1_D_after_ABX, spo_2_D_after_ABX, spo_3_D_after_ABX, spo_4_D_after_ABX,
                           spo_5_D_after_ABX, spo_6_D_after_ABX, spo_7_D_after_ABX]


### Recovery of gut microbiota of healthy adults following antibiotic exposure ###
os.chdir(r"C:\Users\USER\OneDrive\Desktop\Antibiotics\Recovery\Data")

#rel_abund_rarefied = pd.read_csv('annotated.mOTU.rel_abund.rarefied.tsv', sep='\t')
rel_abund_rarefied = pd.read_csv('annotated.mOTU.rel_abund.tsv', sep='\t')

rel_abund_rarefied = filter_data(rel_abund_rarefied)

baseline_columns = ['ERAS1_Dag0', 'ERAS2_Dag0', 'ERAS3_Dag0', 'ERAS4_Dag0', 'ERAS5_Dag0',
                    'ERAS6_Dag0', 'ERAS7_Dag0', 'ERAS8_Dag0', 'ERAS9_Dag0', 'ERAS10_Dag0',
                    'ERAS11_Dag0', 'ERAS12_Dag0']

baseline_columns_appear_4 = ['ERAS2_Dag0', 'ERAS3_Dag0', 'ERAS4_Dag0', 'ERAS5_Dag0',
                             'ERAS6_Dag0', 'ERAS7_Dag0', 'ERAS9_Dag0', 'ERAS11_Dag0', 'ERAS12_Dag0']

columns_4 = ['ERAS2_Dag4opt', 'ERAS3_Dag4', 'ERAS4_Dag4opt', 'ERAS5_Dag4', 'ERAS6_Dag4opt',
             'ERAS7_Dag4opt', 'ERAS9_Dag4', 'ERAS11_Dag4opt', 'ERAS12_Dag4opt']

columns_8 = ['ERAS1_Dag8',  'ERAS2_Dag8', 'ERAS3_Dag8', 'ERAS4_Dag8opt', 'ERAS5_Dag8', 'ERAS6_Dag8opt', 'ERAS7_Dag8',
             'ERAS8_Dag8', 'ERAS9_Dag8', 'ERAS10_Dag8', 'ERAS11_Dag8', 'ERAS12_Dag8']

columns_8_appear_4 = ['ERAS2_Dag8', 'ERAS3_Dag8', 'ERAS4_Dag8opt', 'ERAS5_Dag8', 'ERAS6_Dag8opt',
                        'ERAS7_Dag8', 'ERAS9_Dag8', 'ERAS11_Dag8', 'ERAS12_Dag8']

columns_42 = ['ERAS1_Dag42', 'ERAS2_Dag42', 'ERAS3_Dag42', 'ERAS4_Dag42', 'ERAS5_Dag42', 'ERAS6_Dag42',
              'ERAS7_Dag42',  'ERAS8_Dag42', 'ERAS9_Dag42', 'ERAS10_Dag42', 'ERAS11_Dag42', 'ERAS12_Dag42']

columns_42_appear_4 = ['ERAS2_Dag42', 'ERAS3_Dag42', 'ERAS4_Dag42', 'ERAS5_Dag42', 'ERAS6_Dag42',
                       'ERAS7_Dag42', 'ERAS9_Dag42', 'ERAS11_Dag42', 'ERAS12_Dag42']

columns_180 = ['ERAS1_Dag180', 'ERAS2_Dag180', 'ERAS3_Dag180', 'ERAS4_Dag180', 'ERAS5_Dag180', 'ERAS6_Dag180',
                'ERAS7_Dag180', 'ERAS8_Dag180', 'ERAS9_Dag180', 'ERAS10_Dag180', 'ERAS11_Dag180', 'ERAS12_Dag180']

columns_180_appear_4 = ['ERAS2_Dag180', 'ERAS3_Dag180', 'ERAS4_Dag180', 'ERAS5_Dag180', 'ERAS6_Dag180',
                        'ERAS7_Dag180', 'ERAS9_Dag180', 'ERAS11_Dag180', 'ERAS12_Dag180']

baseline_rel_abund_rarefied = rel_abund_rarefied[baseline_columns].values
baseline_rel_abund_rarefied = baseline_rel_abund_rarefied.T
baseline_rel_abund_rarefied = normalize_cohort(baseline_rel_abund_rarefied)

baseline_rel_abund_rarefied_appear_4 = rel_abund_rarefied[baseline_columns_appear_4].values
baseline_rel_abund_rarefied_appear_4 = baseline_rel_abund_rarefied_appear_4.T
baseline_rel_abund_rarefied_appear_4 = normalize_cohort(baseline_rel_abund_rarefied_appear_4)

rel_abund_rarefied_4 = rel_abund_rarefied[columns_4].values
rel_abund_rarefied_4 = rel_abund_rarefied_4.T
rel_abund_rarefied_4 = normalize_cohort(rel_abund_rarefied_4)

rel_abund_rarefied_8 = rel_abund_rarefied[columns_8].values
rel_abund_rarefied_8 = rel_abund_rarefied_8.T
rel_abund_rarefied_8 = normalize_cohort(rel_abund_rarefied_8)

rel_abund_rarefied_8_appear_4 = rel_abund_rarefied[columns_8_appear_4].values
rel_abund_rarefied_8_appear_4 = rel_abund_rarefied_8_appear_4.T
rel_abund_rarefied_8_appear_4 = normalize_cohort(rel_abund_rarefied_8_appear_4)

rel_abund_rarefied_42 = rel_abund_rarefied[columns_42].values
rel_abund_rarefied_42 = rel_abund_rarefied_42.T
rel_abund_rarefied_42 = normalize_cohort(rel_abund_rarefied_42)

rel_abund_rarefied_42_appear_4 = rel_abund_rarefied[columns_42_appear_4].values
rel_abund_rarefied_42_appear_4 = rel_abund_rarefied_42_appear_4.T
rel_abund_rarefied_42_appear_4 = normalize_cohort(rel_abund_rarefied_42_appear_4)

rel_abund_rarefied_180 = rel_abund_rarefied[columns_180].values
rel_abund_rarefied_180 = rel_abund_rarefied_180.T
rel_abund_rarefied_180 = normalize_cohort(rel_abund_rarefied_180)

rel_abund_rarefied_180_appear_4 = rel_abund_rarefied[columns_180_appear_4].values
rel_abund_rarefied_180_appear_4 = rel_abund_rarefied_180_appear_4.T
rel_abund_rarefied_180_appear_4 = normalize_cohort(rel_abund_rarefied_180_appear_4)

subject_2_samples = np.vstack([rel_abund_rarefied_4[0, :], rel_abund_rarefied_8_appear_4[0, :],
                               rel_abund_rarefied_42_appear_4[0, :], rel_abund_rarefied_180_appear_4[0, :]])

subject_2_samples = normalize_cohort(subject_2_samples)

subject_3_samples = np.vstack([rel_abund_rarefied_4[1, :], rel_abund_rarefied_8_appear_4[1, :],
                               rel_abund_rarefied_42_appear_4[1, :], rel_abund_rarefied_180_appear_4[1, :]])

subject_3_samples = normalize_cohort(subject_3_samples)

subject_4_samples = np.vstack([rel_abund_rarefied_4[2, :], rel_abund_rarefied_8_appear_4[2, :],
                               rel_abund_rarefied_42_appear_4[2, :], rel_abund_rarefied_180_appear_4[2, :]])

subject_4_samples = normalize_cohort(subject_4_samples)

subject_5_samples = np.vstack([rel_abund_rarefied_4[3, :], rel_abund_rarefied_8_appear_4[3, :],
                               rel_abund_rarefied_42_appear_4[3, :], rel_abund_rarefied_180_appear_4[3, :]])

subject_5_samples = normalize_cohort(subject_5_samples)

subject_6_samples = np.vstack([rel_abund_rarefied_4[4, :], rel_abund_rarefied_8_appear_4[4, :],
                               rel_abund_rarefied_42_appear_4[4, :], rel_abund_rarefied_180_appear_4[4, :]])

subject_6_samples = normalize_cohort(subject_6_samples)

subject_7_samples = np.vstack([rel_abund_rarefied_4[5, :], rel_abund_rarefied_8_appear_4[5, :],
                               rel_abund_rarefied_42_appear_4[5, :], rel_abund_rarefied_180_appear_4[5, :]])

subject_7_samples = normalize_cohort(subject_7_samples)

subject_9_samples = np.vstack([rel_abund_rarefied_4[6, :], rel_abund_rarefied_8_appear_4[6, :],
                               rel_abund_rarefied_42_appear_4[6, :], rel_abund_rarefied_180_appear_4[6, :]])

subject_9_samples = normalize_cohort(subject_9_samples)

subject_11_samples = np.vstack([rel_abund_rarefied_4[7, :], rel_abund_rarefied_8_appear_4[7, :],
                                rel_abund_rarefied_42_appear_4[7, :], rel_abund_rarefied_180_appear_4[7, :]])

subject_11_samples = normalize_cohort(subject_11_samples)

subject_12_samples = np.vstack([rel_abund_rarefied_4[8, :], rel_abund_rarefied_8_appear_4[8, :],
                                rel_abund_rarefied_42_appear_4[8, :], rel_abund_rarefied_180_appear_4[8, :]])

subject_12_samples = normalize_cohort(subject_12_samples)

Spo_ABX_cohort = ABX_cohort[14:, :]

Spo_baseline_cohort = baseline_cohort[14:, :]

Spo_post_ABX_cohort = np.vstack([spo_1_D_after_ABX[-1, :], spo_2_D_after_ABX[-1, :], spo_3_D_after_ABX[-1, :],
                                 spo_4_D_after_ABX[-1, :], spo_5_D_after_ABX[-1, :], spo_6_D_after_ABX[-1, :],
                                 spo_7_D_after_ABX[-1, :]])

aFMT_post_ABX_cohort = np.vstack([aFMT_1_samples[-1, :], aFMT_2_samples[-1, :], aFMT_3_samples[-1, :],
                                  aFMT_4_samples[-1, :], aFMT_5_samples[-1, :], aFMT_6_samples[-1, :]])

#### Gut Bacterial Microbiota and its Resistome Rapidly Recover to Basal State Levels after Short-term ###

os.chdir(r'C:\Users\USER\OneDrive\Desktop\Antibiotics\Gut\Data')
data = pd.read_excel('Gut_Bacterial_Microbiota_OTU_table.xlsx')
import plotly.graph_objects as go

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

from Data_manipulation.rarify import Rarify
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