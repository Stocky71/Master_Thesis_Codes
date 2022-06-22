import sys
import os
import os.path

base_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, base_path + "/resources/external/")
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style

matplotlib.style.use("./resources/mpl/paper.mplstyle")
from matplotlib.font_manager import FontProperties
import scipy.stats

import pickle as pkl

import data_loader
import weighter
import binning
import fc

# Load data/MC. By default load_mc loads events at energies >60 TeV, but we want to plot all events.
mc_filenames = [
    "./resources/data/HESE_mc_observable.json",
    "./resources/data/HESE_mc_flux.json",
    "./resources/data/HESE_mc_truth.json",
]
mc = data_loader.load_mc(mc_filenames, emin=10.0e3)
mc_above1 = mc_above2 = mc
ang_100 = math.pi * 10.0 / 18.0 #100 is 80 from Earth vertical
ang_110 = math.pi * 11.0 / 18.0 #110 is 70 from Earth vertical
mc_above1 = mc_above1[np.where(mc_above1['recoZenith'] >= ang_100)]
mc_above2 = mc_above2[np.where(mc_above2['recoZenith'] >= ang_110)]

data = data_loader.load_data("./resources/data/HESE_data.json", emin=10.0e3)
data_above1 = data_above2 = data
data_above1 = data_above1[np.where(data_above1['recoZenith'] >= ang_100)]
data_above2 = data_above2[np.where(data_above2['recoZenith'] >= ang_110)]

e_edges, _, _ = binning.get_bins(emin=10.0e3)
bin_centers = 10.0 ** (0.5 * (np.log10(e_edges[:-1]) + np.log10(e_edges[1:])))

n_events_above1, _ = np.histogram(data_above1["recoDepositedEnergy"], bins=e_edges)
n_events_above2, _ = np.histogram(data_above2["recoDepositedEnergy"], bins=e_edges)

###############################################################################
parameter_names = [
    "cr_delta_gamma",
    "nunubar_ratio",
    "anisotropy_scale",
    "astro_gamma",
    "astro_norm",
    "conv_norm",
    "epsilon_dom",
    "epsilon_head_on",
    "muon_norm",
    "kpi_ratio",
    "prompt_norm",
]

# We plot using the best fit parameters found by HESE_fit.py
params = np.array(
    [
        -0.05309302,
        0.99815326,
        1.000683,
        2.87375956,
        6.36488608,
        1.00621679,
        0.95192328,
        -0.0548763,
        1.18706341,
        1.00013744,
        0.0,
    ]
)

livtm = 227708167.68

params_dict = dict(zip(parameter_names, params))

params_zeroed = params_dict.copy()

component_order = [
    ("muon_norm", "Atmo. Muons"),
    ("conv_norm", "Atmo. Conv."),
    ("prompt_norm", "Atmo. Prompt"),
    ("astro_norm", "Astro."),
]

for (zeroed_norm, _) in component_order:
    params_zeroed[zeroed_norm] = 0.0

# We want to separate the histogram by components, so we separately get the weights
# where the all normalization parameters but one are set to zero
colors_above = []
labels_above = []
cm = plt.get_cmap("inferno")
color_scale = [cm(x) for x in [0.2, 0.55, 0.75, 0.9]]
for i, (zeroed_norm, zeroed_label) in enumerate(component_order):
    if params_dict[zeroed_norm] == 0.0:
        continue
    colors_above.append(color_scale[i])
    labels_above.append(zeroed_label)

xerr = [bin_centers - e_edges[:-1], e_edges[1:] - bin_centers]
# The error bars are obtained using fc.py, a function writen by Austin Schneider
# and found at https://github.com/austinschneider/feldman_cousins/blob/master/fc.py
one_sigma_proportion = scipy.special.erf(1.0 / np.sqrt(2.0))
yerr_above1 = np.array(
    [
        (lambda x: [k - x[0], x[1] - k])(
            fc.poisson_interval(k, alpha=one_sigma_proportion)
        )
        for k in n_events_above1
    ]
).T
yerr_above2 = np.array(
    [
        (lambda x: [k - x[0], x[1] - k])(
            fc.poisson_interval(k, alpha=one_sigma_proportion)
        )
        for k in n_events_above2
    ]
).T

###############################################################################
###############################################################################

mass_ran = list(np.linspace(500, 1500, 11))
coupling_ran = [0.5, 1.0, 1.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]

        wghtm_above1 = weighter.Weighter(mc_above1) 
        wghtm_above2 = weighter.Weighter(mc_above2) 

        wghtm_above1.set_lpmass_coupling(mass_lq, coupling, 'A', ang=ang_100)
        wghtm_above2.set_lpmass_coupling(mass_lq, coupling, 'A', ang=ang_110)

# We want to separate the histogram by components, so we separately get the weights
# where the all normalization parameters but one are set to zero
        weights_above1 = []
        weights_above2 = []
        for k, (zeroed_norm, zeroed_label) in enumerate(component_order):
            if params_dict[zeroed_norm] == 0.0:
                continue
            p_cpy = params_zeroed.copy()
            p_cpy[zeroed_norm] = params_dict[zeroed_norm]
            weights_above1.append(
                wghtm_above1.get_weights(livtm, p_cpy.keys(), p_cpy.values())[0]
            )
            weights_above2.append(
                wghtm_above2.get_weights(livtm, p_cpy.keys(), p_cpy.values())[0]
            )

###############################################################################
#Save diag data into file
        above_error_bar1 = []
        above_error_bar1.extend(xerr)
        above_error_bar1.extend(yerr_above1)
        above_error_bar1.extend(bin_centers.tolist())
        above_error_bar1.extend(n_events_above1.tolist())

        file_name_1 = ('datos_diff_en/ang_100/lq_' + str(int(mass_lq)) + 
                       '/above_diag_errorbar-' + str(coupling) + '.pkl')
        with open(file_name_1, 'wb') as f:
            pkl.dump(above_error_bar1, f)
            
        above_error_bar2 = []
        above_error_bar2.extend(xerr)
        above_error_bar2.extend(yerr_above2)
        above_error_bar2.extend(bin_centers.tolist())
        above_error_bar2.extend(n_events_above2.tolist())

        file_name_2 = ('datos_diff_en/ang_110/lq_' + str(int(mass_lq)) + 
                       '/above_diag_errorbar-' + str(coupling) + '.pkl')
        with open(file_name_2, 'wb') as f:
            pkl.dump(above_error_bar2, f)
            
        above_data1 = []
        above_data1.extend(mc_above1["recoDepositedEnergy"]) 
        above_data1.extend(weights_above1) #.tolist())
        above_data1.extend(e_edges.tolist())

        file_name_3 = ('datos_diff_en/ang_100/lq_' + str(int(mass_lq)) + 
                       '/above_diag_data-' + str(coupling) + '.pkl')
        with open(file_name_3, 'wb') as f:
            pkl.dump(above_data1, f)
            
        above_data2 = []
        above_data2.extend(mc_above2["recoDepositedEnergy"]) 
        above_data2.extend(weights_above2) #.tolist())
        above_data2.extend(e_edges.tolist())

        file_name_4 = ('datos_diff_en/ang_110/lq_' + str(int(mass_lq)) + 
                       '/above_diag_data-' + str(coupling) + '.pkl')
        with open(file_name_4, 'wb') as f:
            pkl.dump(above_data2, f)
            

#eof