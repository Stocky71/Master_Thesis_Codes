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
mc_above = mc_below = mc
mc_above = mc_above[np.where(mc_above['recoZenith'] >= math.pi/2)]
mc_below = mc_below[np.where(mc_below['recoZenith'] <  math.pi/2)]

data = data_loader.load_data("./resources/data/HESE_data.json", emin=10.0e3)
data_above = data_below = data
data_above = data_above[np.where(data_above['recoZenith'] >= math.pi/2)]
data_below = data_below[np.where(data_below['recoZenith'] <  math.pi/2)]

e_edges, _, _ = binning.get_bins(emin=10.0e3)
bin_centers = 10.0 ** (0.5 * (np.log10(e_edges[:-1]) + np.log10(e_edges[1:])))

n_events_above, _ = np.histogram(data_above["recoDepositedEnergy"], bins=e_edges)
n_events_below, _ = np.histogram(data_below["recoDepositedEnergy"], bins=e_edges)

weight_maker_above = weighter.Weighter(mc_above)
weight_maker_below = weighter.Weighter(mc_below)

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

livetime = 227708167.68

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
weights_above = []
weights_below = []
colors_above = []
colors_below = []
labels_above = []
labels_below = []
cm = plt.get_cmap("inferno")
color_scale = [cm(x) for x in [0.2, 0.55, 0.75, 0.9]]
for i, (zeroed_norm, zeroed_label) in enumerate(component_order):
    if params_dict[zeroed_norm] == 0.0:
        continue
    p_copy = params_zeroed.copy()
    p_copy[zeroed_norm] = params_dict[zeroed_norm]
    weights_above.append(
        weight_maker_above.get_weights(livetime, p_copy.keys(), p_copy.values())[0]
    )
    colors_above.append(color_scale[i])
    labels_above.append(zeroed_label)

    weights_below.append(
        weight_maker_below.get_weights(livetime, p_copy.keys(), p_copy.values())[0]
    )
    colors_below.append(color_scale[i])
    labels_below.append(zeroed_label)

xerr = [bin_centers - e_edges[:-1], e_edges[1:] - bin_centers]
# The error bars are obtained using fc.py, a function writen by Austin Schneider
# and found at https://github.com/austinschneider/feldman_cousins/blob/master/fc.py
one_sigma_proportion = scipy.special.erf(1.0 / np.sqrt(2.0))
yerr_above = np.array(
    [
        (lambda x: [k - x[0], x[1] - k])(
            fc.poisson_interval(k, alpha=one_sigma_proportion)
        )
        for k in n_events_above
    ]
).T
yerr_below = np.array(
    [
        (lambda x: [k - x[0], x[1] - k])(
            fc.poisson_interval(k, alpha=one_sigma_proportion)
        )
        for k in n_events_below
    ]
).T

#Save plot data into files
"""
above_error_bar = []
above_error_bar.extend(xerr)
above_error_bar.extend(yerr_above)
above_error_bar.extend(bin_centers.tolist())
above_error_bar.extend(n_events_above.tolist())

file_name_1 = 'diff_data/SM/above_diag_errorbar.pkl'
with open(file_name_1, 'wb') as f:
    pkl.dump(above_error_bar, f)
    
below_error_bar = []
below_error_bar.extend(xerr)
below_error_bar.extend(yerr_below)
below_error_bar.extend(bin_centers.tolist())
below_error_bar.extend(n_events_below.tolist())

file_name_2 = 'diff_data/SM/below_diag_errorbar.pkl'
with open(file_name_2, 'wb') as f:
    pkl.dump(below_error_bar, f)

above_data = []
above_data.extend(mc_above["recoDepositedEnergy"]) 
above_data.extend(weights_above) #.tolist())
above_data.extend(e_edges.tolist())

file_name_3 = 'diff_data/SM/above_diag_data.pkl'
with open(file_name_3, 'wb') as f:
    pkl.dump(above_data, f)
    
below_data = []
below_data.extend(mc_below["recoDepositedEnergy"])
below_data.extend(weights_below) #.tolist())
below_data.extend(e_edges.tolist())

file_name_4 = 'diff_data/SM/below_diag_data.pkl'
with open(file_name_4, 'wb') as f:
    pkl.dump(below_data, f)
"""

# Plot above histogram
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1.0e-2, 4.0e1)
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Events per 2635 days")

plt.errorbar(
    bin_centers,
    n_events_above,
    xerr=xerr,
    yerr=yerr_above,
    color="black",
    marker=None,
    label="Data",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_a,
)

labels_above[2] = 'Astro. Above'

plt.hist(
    len(weights_above) * [mc_above["recoDepositedEnergy"]],
    weights=weights_above,
    bins=e_edges,
    histtype="bar",
    stacked=True,
    label=labels_above,
    color=colors_above,
#	ax=ax_a,
#	alpha=0.4,
)

# This applies the visual filter that covers events that are not used
# in the analysis
mask_color = "#7ab9f3"
ax_a.fill_between(
    [10e3, 60e3],
    [1e3, 1e3],
    [0, 0],
    edgecolor="none",
    linewidth=0.0,
    facecolor=mask_color,
    zorder=3,
    alpha=0.7,
)
plt.axvline(x=60e3, linestyle="dashed")

font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

# Simply reverses the order of the legend labels.
handles_above, labels_above = ax_a.get_legend_handles_labels()
ax_a.legend(handles_above[::-1], labels_above[::-1])

plt.tight_layout()
fig_name_above = 'above_diag.png'
plt.savefig(fig_name_above)
#plt.show()

# Plot below histogram
fit_b, ax_b = plt.subplots(figsize=(7, 5))
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1.0e-2, 4.0e1)
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Events per 2635 days")

plt.errorbar(
    bin_centers,
    n_events_below,
    xerr=xerr,
    yerr=yerr_below,
    color="black",
    marker=None,
    label="Data",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_b,
)

labels_below[2] = 'Astro. Below'

plt.hist(
    len(weights_below) * [mc_below["recoDepositedEnergy"]],
    weights=weights_below,
    bins=e_edges,
    histtype="bar",
    stacked=True,
    label=labels_below,
    color=colors_below,
#	ax=ax_b,
#	alpha=0.4,
)

# This applies the visual filter that covers events that are not used
# in the analysis
mask_color = "#7ab9f3"
ax_b.fill_between(
    [10e3, 60e3],
    [1e3, 1e3],
    [0, 0],
    edgecolor="none",
    linewidth=0.0,
    facecolor=mask_color,
    zorder=3,
    alpha=0.2,
)
plt.axvline(x=60e3, linestyle="dashed")

font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

# Simply reverses the order of the legend labels.
handles_below, labels_below = ax_b.get_legend_handles_labels()
ax_b.legend(handles_below[::-1], labels_below[::-1])

plt.tight_layout()
#fig_name_below = 'below_diag.png'
#plt.savefig(fig_name_below)
plt.show()

#eof