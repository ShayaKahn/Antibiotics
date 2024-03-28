import numpy as np
from scipy.spatial.distance import braycurtis
import plotly.graph_objects as go
from data import *

def normalize_cohort(cohort):
    # normalization function
    if cohort.ndim == 1:
        cohort_normalized = cohort / cohort.sum()
    else:
        cohort_normalized = cohort / np.linalg.norm(cohort, ord=1, axis=1, keepdims=True)
    return cohort_normalized

def create_cohort(pool_size, num_samples):
    prob_vector = np.random.uniform(0.1, 0.9, pool_size)
    cohort = np.zeros((num_samples, pool_size))
    for i in range(0, num_samples):
        for j in range(0, pool_size):
            if np.random.uniform(0, 1) < prob_vector[j]:
                cohort[i, j] = np.random.normal(0.5, 0.15)#np.random.uniform(0, 1)
            else:
                cohort[i, j] = 0
    cohort = normalize_cohort(cohort)
    return cohort

def create_ABX_sample(baseline_sample, prob):
    ABX_sample = baseline_sample.copy()
    mask = np.random.rand(ABX_sample.shape[0]) < prob
    ABX_sample[mask] = 0
    ABX_sample = normalize_cohort(ABX_sample)
    return ABX_sample

def create_test_sample(baseline_sample, ABX_sample, cohort):
    test_sample = np.zeros(baseline_sample.shape[0])
    indexes = np.where((baseline_sample != 0) & (ABX_sample != 0))[0]
    test_sample[indexes] = baseline_sample[indexes]
    random_indices = np.random.choice(cohort.shape[0], size=cohort.shape[1])
    for i in range(0, test_sample.shape[0]):
        if i not in indexes:
            test_sample[i] = cohort[random_indices[i], i]
    test_sample = normalize_cohort(test_sample)
    return test_sample

def create_noisy_baseline(baseline_sample, sigma, k):
    noisy_baseline = np.tile(baseline_sample, (k, 1))
    noise = np.random.normal(0, sigma, noisy_baseline.shape)
    zero_columns = np.where(np.sum(noisy_baseline, axis=0) == 0)[0]
    noise[:, zero_columns] = 0
    noisy_baseline += noise
    noisy_baseline = normalize_cohort(noisy_baseline)
    return noisy_baseline

prob = 0.75
pool_size = 100
num_samples = 1000
cohort = create_cohort(pool_size, num_samples)
baseline_bc = np.zeros(num_samples)
Base_bc = []
ABX_bc = []
test_bc = []
n = 3
sigma = 0.01
k = 10

for base in cohort[0:n+1, :]:
    baseline_sample = base
    noisy_baseline = create_noisy_baseline(baseline_sample, sigma, k)
    ABX_sample = create_ABX_sample(baseline_sample, prob)
    test_sample = create_test_sample(baseline_sample, ABX_sample, cohort)
    Base_bc.append(np.mean([braycurtis(baseline_sample, base) for base in noisy_baseline]))
    ABX_bc.append(braycurtis(baseline_sample, ABX_sample))
    test_bc.append(braycurtis(baseline_sample, test_sample))

# Simulations plot
categories = ['Baseline', 'ABX', 'Post-ABX', ' ']
bc_vals_1 = [Base_bc[0], ABX_bc[0], test_bc[0], test_bc[0]]
bc_vals_2 = [Base_bc[1], ABX_bc[1], test_bc[1], test_bc[1]]
bc_vals_3 = [Base_bc[2], ABX_bc[2], test_bc[2], test_bc[2]]
bc_vals_4 = [Base_bc[3], ABX_bc[3], test_bc[3], test_bc[3]]

step_categories = [categories[0]] + [cat for cat in categories for _ in (0, 1)][1:-1] + [categories[-1]]
step_values1 = [bc_vals_1[0]] + [val for val in bc_vals_1 for _ in (0, 1)][1:] + [bc_vals_1[-1]]
step_values2 = [bc_vals_2[0]] + [val for val in bc_vals_2 for _ in (0, 1)][1:] + [bc_vals_2[-1]]
step_values3 = [bc_vals_3[0]] + [val for val in bc_vals_3 for _ in (0, 1)][1:] + [bc_vals_3[-1]]
step_values4 = [bc_vals_4[0]] + [val for val in bc_vals_4 for _ in (0, 1)][1:] + [bc_vals_4[-1]]

fig = go.Figure()

start_category = 'ABX'
end_category = 'Post-ABX'

# Add a grey background from the second to the third category
fig.add_vrect(x0=start_category, x1=end_category, fillcolor="grey",
              opacity=0.25, layer="below", line_width=0)

fig.add_trace(go.Scatter(x=step_categories, y=step_values1,
                         mode='lines',
                         line=dict(color='red', width=2, shape='hv')))

fig.add_trace(go.Scatter(x=step_categories, y=step_values2,
                         mode='lines',
                         line=dict(color='blue', width=2, shape='hv')))

fig.add_trace(go.Scatter(x=step_categories, y=step_values3,
                         mode='lines',
                         line=dict(color='green', width=2, shape='hv')))

fig.add_trace(go.Scatter(x=step_categories, y=step_values4,
                         mode='lines',
                         line=dict(color='black', width=2, shape='hv')))

fig.update_layout(
    xaxis=dict(title='', showline=True, linewidth=2, linecolor='#000',
               tickfont=dict(family='Times New Roman', size=30), showticklabels=False, showgrid=False,
               zeroline=False),
    yaxis=dict(title='Distance from baseline (BC)', showline=True, linewidth=2, linecolor='#000',
               tickfont=dict(family='Times New Roman', size=30),
               title_font=dict(size=30, family='Times New Roman'),
               showgrid=False, zeroline=False),
    plot_bgcolor='white',
    width=600,
    height=600,
    showlegend=False,
    title=dict(text="Simplified model of 'general recovery'", x=0.5, xanchor='center',
               font=dict(size=35, family='Times New Roman'))
)

fig.show()

# Real data plot
BC_spo_1 = np.array([0.4381722927692182, 0.7161840953683432, 0.3792010968054209, 0.4524277036555165,
                    0.5146196321785041, 0.37883315921869537, 0.4559889709821168,
                    0.5599109954338625, 0.746749436711404, 0.8397540224349662, 0.907936388123862,
                    0.8038333885489596, 0.78709808786586, 0.7001143489949009, 0.6642356862262513,
                    0.6871088684425041, 0.7126399642722918, 0.5913921407724099, 0.5810222752342167,
                    0.6066422286138706, 0.8076927654471137, 0.7751710519455377, 0.6210252987773932,
                    0.5844826427115435])
BC_spo_2 = np.array([0.2695585896232256, 0.23573894515108354, 0.21401330429775187, 0.30784762705994445,
                    0.28061000155283256, 0.2453866546723531, 0.28407023697809525, 0.3331659570104275,
                    0.28406442296785267, 0.914101267587213, 0.9475552231136027, 0.9536485189411132,
                    0.9372188324341344, 0.9171700828946213, 0.901095317650553, 0.766618703754275,
                    0.5703212747129778, 0.5870670152366436, 0.5337525512265516, 0.4611571956959236,
                    0.4386018220216054, 0.4341354198290753, 0.4648400100457194, 0.4100675934841149,
                    0.47745870339935453, 0.4374416795036567, 0.4536064283014144, 0.40500286291299703,
                    0.454054009814847, 0.43651884817571773])
BC_spo_2_ind = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
BC_spo_3 = np.array([0.2544935748412424, 0.3457213801108197, 0.2854532316922424, 0.23295113048848912,
                    0.2305854073155147, 0.28555875408839543, 0.3191896982777785, 0.25949868252342884,
                    0.30960780145224887, 0.7161812263103172, 0.8823288649085731, 0.6433633409272544,
                    0.6817433926473738, 0.6443052854108211, 0.37071296411122784, 0.38013681451048686,
                    0.39119283352956247, 0.35861877789546603, 0.43850253319504723, 0.38619045058313445,
                    0.34416218806631627, 0.4669538536531782, 0.6014778864188182, 0.4442602216402718,
                    0.43595370130765715, 0.4347473537583602, 0.3832292274130206, 0.42402904374816913,
                    0.3918920611677202])
BC_spo_3_ind = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 23, 24, 25, 26]
BC_spo_4 = np.array([0.2198020064766865, 0.2220938847742024, 0.21109552899263032, 0.24602139817953309,
                     0.3824297740490442, 0.317433345717522, 0.3297913339477537, 0.28626963808907996,
                     0.8636188875278822, 0.7357837889689969, 0.704659238842735, 0.6070909370352453,
                     0.6122355507524918, 0.6327597284703913, 0.4591975800528892, 0.45461443063984314,
                     0.45972366763246425, 0.4598041142190723, 0.4490123884516213, 0.4478024056820249,
                     0.36490959889221103, 0.4164688845588933, 0.41933070445346626, 0.3894480351838785,
                     0.47145571608313586])
BC_spo_4_ind = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
BC_spo_5 = np.array([0.28016879878611733, 0.20964441040701734, 0.2771782686252533, 0.22174030631091482,
                     0.24808435200677115, 0.2401857913125556, 0.22930773981923294, 0.28154165286891925,
                     0.5957443913195188, 0.9304188647532489, 0.9986870572482843, 0.99853358582907,
                     0.9830774994270525, 0.9969253432906956, 0.9729506585152751, 0.9223251017370349,
                     0.8517655534372228, 0.813522817830474, 0.815300668837513, 0.8637479502358099,
                     0.8357236704680806, 0.6873710118681492, 0.6614347497002554, 0.6173471767685059,
                     0.5891686241770256, 0.5982603388312585, 0.6269986783714004, 0.5733436131506647,
                     0.6315723362699316])
BC_spo_5_ind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
BC_spo_6 = np.array([0.22067683286878362, 0.20476990684065988, 0.20918474258385983, 0.23062453475725458,
                     0.21276183127899218, 0.21442947438099194, 0.4371922937668229, 0.9688405443656471,
                     0.9852895647092531, 0.9800283848181888, 0.9813951629047933, 0.9515419927220791,
                     0.9258624555647526, 0.8821244850672543, 0.8628593488227395, 0.8337045571505126,
                     0.8396912996071265, 0.7350137832648637, 0.6574506533545611, 0.6489985912275359,
                     0.512129655279107, 0.5920484725467415])
BC_spo_7 = np.array([0.11718714486663985, 0.11885613575538012, 0.1218454345608244, 0.16680608125763638,
                     0.135663542061151, 0.12469498746656044, 0.15281563137822332, 0.16420197451257845,
                     0.16520109870287575, 0.423224719854402, 0.642218443702192, 0.719443836869992,
                     0.5870356721837834, 0.6018119242376917, 0.5569496347974869, 0.544466171292894,
                     0.5357381972808305, 0.5564484466483453, 0.5471413360930925, 0.5261411886724925,
                     0.45572332270967936, 0.49034676965113266, 0.436718470855983, 0.5211917681951722,
                     0.4752700040231054])

x = np.arange(0, 22, 1)
print(np.size(BC_spo_2))
print(np.size(BC_spo_3))
print(np.size(BC_spo_4))
print(np.size(BC_spo_5))

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=BC_spo_2[BC_spo_2_ind], mode='lines', line=dict(color='red')))
fig.add_trace(go.Scatter(x=x, y=BC_spo_3[BC_spo_3_ind], mode='lines', line=dict(color='red')))
fig.add_trace(go.Scatter(x=x, y=BC_spo_4[BC_spo_4_ind], mode='lines', line=dict(color='red')))
fig.add_trace(go.Scatter(x=x, y=BC_spo_5[BC_spo_5_ind], mode='lines', line=dict(color='red')))

fig.update_layout(
    width=600,
    height=600,
    xaxis=dict(
        showticklabels=True,
        showgrid=False,
        zeroline=True,
        showline=True,
        linecolor='black',
        linewidth=2,
        title_font=dict(size=30),
        tickfont=dict(color='white', size=1)
    ),
    yaxis=dict(
        title='Distance from baseline (BC)',
        showgrid=False,
        zeroline=True,
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(family='Times New Roman', size=30),
        title_font=dict(size=30, family='Times New Roman')
    ),
    plot_bgcolor='white',
    showlegend=False,
    title=dict(text='Post-ABX recovery', x=0.5, xanchor='center',
               font=dict(size=35, family='Times New Roman'))
)

start_highlight = x[6]
end_highlight = x[11]

fig.add_vrect(
    x0=start_highlight,
    x1=end_highlight,
    fillcolor="grey",
    opacity=0.25,
    layer="below",
    line_width=0,
)

fig.show()
