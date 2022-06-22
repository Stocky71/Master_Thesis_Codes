import numpy  as np
import pickle as pkl
import matplotlib.pyplot as plt
#import matplotlib.style
#import matplotlib
#matplotlib.style.use("./resources/mpl/paper.mplstyle")
from matplotlib.font_manager import FontProperties

mass_lq  = 500.0
coupling = 0.5

# Load data/MC. By default load_mc loads events at energies >60 TeV, but we want to plot all events.
###############################################################################
# LOAD SM DATA FOR PLOT
#file_name_1 = 'datos_diff_en/SM/above_diag_errorbar.pkl'
file_name_1 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
               '/above_diag_errorbar-' + str(coupling) + '.pkl')
with open(file_name_1, 'rb') as f:
    arr_er_a = pkl._Unpickler(f)
    arr_er_a.encoding = 'latin1'
    arr_er_a = arr_er_a.load()
    
xerr_a        = [arr_er_a[0], arr_er_a[1]]   
yerr_a        = [arr_er_a[2], arr_er_a[3]]  
bin_centers_a = arr_er_a[4:33]
n_events_er_a = arr_er_a[33:]

#file_name_2 = 'datos_diff_en/SM/below_diag_errorbar.pkl'
file_name_2 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
               '/below_diag_errorbar-' + str(coupling) + '.pkl')
with open(file_name_2, 'rb') as f:
    arr_er_b = pkl._Unpickler(f)
    arr_er_b.encoding = 'latin1'
    arr_er_b = arr_er_b.load()
    
xerr_b        = [arr_er_b[0], arr_er_b[1]]   
yerr_b        = [arr_er_b[2], arr_er_b[3]]  
bin_centers_b = arr_er_b[4:33]
n_events_er_b = arr_er_b[33:]

#file_name_3 = 'datos_diff_en/SM/above_diag_data.pkl'
file_name_3 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
               '/above_diag_data-' + str(coupling) + '.pkl')
with open(file_name_3, 'rb') as f:
    arr_dt_a = pkl._Unpickler(f)
    arr_dt_a.encoding = 'latin1'
    arr_dt_a = arr_dt_a.load()
    
mc_a       = arr_dt_a[:236365]   
weights_a  = [arr_dt_a[236365], arr_dt_a[236366], arr_dt_a[236367]] 
n_events_a = arr_dt_a[236368:]

#file_name_4 = 'datos_diff_en/SM/below_diag_data.pkl'
file_name_4 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
               '/below_diag_data-' + str(coupling) + '.pkl')
with open(file_name_4, 'rb') as f:
    arr_dt_b = pkl._Unpickler(f)
    arr_dt_b.encoding = 'latin1'
    arr_dt_b = arr_dt_b.load()
    
mc_b       = arr_dt_b[:583005]    
weights_b  = [arr_dt_b[583005], arr_dt_b[583006], arr_dt_b[583007]] 
n_events_b = arr_dt_b[583008:]

###############################################################################
component_order = [("muon_norm", "Atmo. Muons"),
                   ("conv_norm", "Atmo. Conv."),
                   #("prompt_norm", "Atmo. Prompt"),
                   ("astro_norm", "Astro."), ]

colors_above = []
colors_below = []
labels_above = []
labels_below = []
cm = plt.get_cmap("inferno")
color_scale = [cm(x) for x in [0.2, 0.55, 0.75, 0.9]]
for i, (zeroed_norm, zeroed_label) in enumerate(component_order):
    colors_above.append(color_scale[i])
    labels_above.append(zeroed_label)

    colors_below.append(color_scale[i])
    labels_below.append(zeroed_label)

###############################################################################

###############################################################################
# Plot above histogram
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1.0e-2, 4.0e1)
plt.xlabel('Deposited Energy [GeV]')
plt.ylabel('Events per 2635 days')
fig_title = ('Deposited energy for leptoquark - Mass:' + str(int(mass_lq)) +
             ' - Coupling:' + str(coupling))
plt.title(fig_title)

plt.errorbar(
    bin_centers_a,
    n_events_er_a,
    xerr=xerr_a,
    yerr=yerr_a,
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
    len(weights_a) * [mc_a],
    weights=weights_a,
    bins=n_events_a,
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
#fig_name_above = 'above_diag.png'
#fig_name_above = 'above_diag_' + str(mass_lq) + '_' + str(coupling) + '.png'
#plt.savefig(fig_name_above, dpi=200)
plt.show()

###############################################################################

###############################################################################
# Plot below histogram
fit_b, ax_b = plt.subplots(figsize=(7, 5))
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1.0e-2, 4.0e1)
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Events per 2635 days")
fig_title = ('Deposited energy for leptoquark - Mass:' + str(int(mass_lq)) +
             ' - Coupling:' + str(coupling))
plt.title(fig_title)

plt.errorbar(
    bin_centers_b,
    n_events_er_b,
    xerr=xerr_b,
    yerr=yerr_b,
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
    len(weights_b) * [mc_b],
    weights=weights_b,
    bins=n_events_b,
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
#fig_name_below = 'below_diag_' + str(mass_lq) + '_' + str(coupling) + '.png'
#plt.savefig(fig_name_below, dpi=200)
plt.show()

#eof