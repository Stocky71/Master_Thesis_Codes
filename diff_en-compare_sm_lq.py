import numpy  as np
import pickle as pkl
import matplotlib.pyplot as plt
#import matplotlib.style
#import matplotlib
#matplotlib.style.use("./resources/mpl/paper.mplstyle")
from matplotlib.font_manager import FontProperties

# Load data/MC. By default load_mc loads events at energies >60 TeV, but we want to plot all events.
###############################################################################
# LOAD SM DATA FOR PLOT
file_name_1 = 'datos_diff_en/SM/above_diag_errorbar.pkl'
with open(file_name_1, 'rb') as f:
    arr_er_a = pkl._Unpickler(f)
    arr_er_a.encoding = 'latin1'
    arr_er_a = arr_er_a.load()
    
xerr_a_sm        = [arr_er_a[0], arr_er_a[1]]   
yerr_a_sm        = [arr_er_a[2], arr_er_a[3]]  
bin_centers_a_sm = arr_er_a[4:33]
n_events_er_a_sm = arr_er_a[33:]

file_name_2 = 'datos_diff_en/SM/below_diag_errorbar.pkl'
with open(file_name_2, 'rb') as f:
    arr_er_b = pkl._Unpickler(f)
    arr_er_b.encoding = 'latin1'
    arr_er_b = arr_er_b.load()
    
xerr_b_sm        = [arr_er_b[0], arr_er_b[1]]   
yerr_b_sm        = [arr_er_b[2], arr_er_b[3]]  
bin_centers_b_sm = arr_er_b[4:33]
n_events_er_b_sm = arr_er_b[33:]

file_name_3 = 'datos_diff_en/SM/above_diag_data.pkl'
with open(file_name_3, 'rb') as f:
    arr_dt_a = pkl._Unpickler(f)
    arr_dt_a.encoding = 'latin1'
    arr_dt_a = arr_dt_a.load()
    
mc_a_sm       = arr_dt_a[:236365]   
weights_a_sm  = [arr_dt_a[236365], arr_dt_a[236366], arr_dt_a[236367]] 
n_events_a_sm = arr_dt_a[236368:]

file_name_4 = 'datos_diff_en/SM/below_diag_data.pkl'
with open(file_name_4, 'rb') as f:
    arr_dt_b = pkl._Unpickler(f)
    arr_dt_b.encoding = 'latin1'
    arr_dt_b = arr_dt_b.load()
    
mc_b_sm       = arr_dt_b[:583005]    
weights_b_sm  = [arr_dt_b[583005], arr_dt_b[583006], arr_dt_b[583007]] 
n_events_b_sm = arr_dt_b[583008:]

counts_a_sm = np.histogram(3*[mc_a_sm], weights=weights_a_sm,
                           bins=n_events_a_sm)[0]
counts_b_sm = np.histogram(3*[mc_b_sm], weights=weights_b_sm,
                           bins=n_events_b_sm)[0]

###############################################################################
# LOAD LQ DATA FOR PLOT
norm_count_a_lq = []
norm_count_b_lq = []

mass_ran = list(np.linspace(500, 1500, 3))
coupling_ran = [0.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]
        
        norm_count_a_lq_comp = []
        norm_count_b_lq_comp = []
        
        arr_er_a = []
        arr_er_b = []
        arr_dt_a = []
        arr_dt_b = []
        
        file_name_1 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
                       '/above_diag_errorbar-' + str(coupling) + '.pkl')
        with open(file_name_1, 'rb') as f:
            arr_er_a = pkl._Unpickler(f)
            arr_er_a.encoding = 'latin1'
            arr_er_a = arr_er_a.load()
            
        xerr_a_lq        = [arr_er_a[0], arr_er_a[1]]   
        yerr_a_lq        = [arr_er_a[2], arr_er_a[3]]  
        bin_centers_a_lq = arr_er_a[4:33]
        n_events_er_a_lq = arr_er_a[33:]

        file_name_2 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
                       '/below_diag_errorbar-' + str(coupling) + '.pkl')
        'datos_diff_en/SM/below_diag_errorbar.pkl'
        with open(file_name_2, 'rb') as f:
            arr_er_b = pkl._Unpickler(f)
            arr_er_b.encoding = 'latin1'
            arr_er_b = arr_er_b.load()
            
        xerr_b_lq        = [arr_er_b[0], arr_er_b[1]]   
        yerr_b_lq        = [arr_er_b[2], arr_er_b[3]]  
        bin_centers_b_lq = arr_er_b[4:33]
        n_events_er_b_lq = arr_er_b[33:]
        
        file_name_3 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
                       '/above_diag_data-' + str(coupling) + '.pkl')
        with open(file_name_3, 'rb') as f:
            arr_dt_a = pkl._Unpickler(f)
            arr_dt_a.encoding = 'latin1'
            arr_dt_a = arr_dt_a.load()
            
        mc_a_lq       = arr_dt_a[:236365]   
        weights_a_lq  = [arr_dt_a[236365], arr_dt_a[236366], arr_dt_a[236367]] 
        n_events_a_lq = arr_dt_a[236368:]

        file_name_4 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
                       '/below_diag_data-' + str(coupling) + '.pkl')
        with open(file_name_4, 'rb') as f:
            arr_dt_b = pkl._Unpickler(f)
            arr_dt_b.encoding = 'latin1'
            arr_dt_b = arr_dt_b.load()
            
        mc_b_lq       = arr_dt_b[:583005]    
        weights_b_lq  = [arr_dt_b[583005], arr_dt_b[583006], arr_dt_b[583007]] 
        n_events_b_lq = arr_dt_b[583008:]
        
    
        counts_a_lq = np.histogram(3*[mc_a_lq], weights=weights_a_lq,
                                   bins=n_events_a_lq)[0]
        counts_b_lq = np.histogram(3*[mc_b_lq], weights=weights_b_lq,
                                   bins=n_events_b_lq)[0]

        norm_count_a_lq_comp = []
        norm_count_b_lq_comp = []
        #error propagation of x/a is x*x_err/a**2
        for k in range(0, len(bin_centers_a_lq)):
            if counts_a_sm[k] == 0.0:
                norm_count_a_lq_comp.append( counts_a_lq[k] / 1e-6 )
            else:
                norm_count_a_lq_comp.append( counts_a_lq[k] / counts_a_sm[k] )
        for k in range(0, len(bin_centers_b_lq)):
            if counts_b_sm[k] == 0.0:
                norm_count_b_lq_comp.append( counts_b_lq[k] / 1e-6 )
            else:
                norm_count_b_lq_comp.append( counts_b_lq[k] / counts_b_sm[k] )

        norm_count_a_lq.append((mass_lq, coupling, norm_count_a_lq_comp, 
                                bin_centers_a_lq))
        norm_count_b_lq.append((mass_lq, coupling, norm_count_b_lq_comp,
                                bin_centers_b_lq))

###############################################################################
#ABOVE DIAGRAM

y_events_er_a_sm = []
y_errorbar_a_sm  = []
y_errorbar_a_sm1 = []
y_errorbar_a_sm2 = []
#error propagation of x/a is x*x_err/a**2
for i in range(0, len(bin_centers_a_sm)):
    if counts_a_sm[i] == 0.0:
        count_val_a = 1e-5
    else:
        count_val_a = counts_a_sm[i]
    y_events_er_a_sm.append( n_events_er_a_sm[i] / count_val_a )
    y_errorbar_a_sm1.append( yerr_a_sm[0][i] / count_val_a )
    y_errorbar_a_sm2.append( yerr_a_sm[1][i] / count_val_a )
y_errorbar_a_sm = [np.array(y_errorbar_a_sm1), np.array(y_errorbar_a_sm2)]   


# Plot above histogram
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.title("Astro. Above")
plt.grid(True)
plt.xscale('log')
plt.xlim(10.0e3, 10.0e6)
#plt.xlim(10.0e3, 2.5e6)
plt.yscale('log') #'linear')
plt.ylim(1e-2, 1.0e2) 
#plt.ylim(-0.5e0, 4.0e0)
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Normalized events per 2635 days")
plt.axhline(y=1, linestyle="dashed", label='Theoretical SM')

plt.errorbar(
    bin_centers_a_sm,
    y_events_er_a_sm, #n_events_er_a_sm,
    #xerr=xerr_a_sm,
    yerr=y_errorbar_a_sm, #yerr_a_sm,
    color="black",
    marker=".",#None,
    label="Data",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_a,
)

color_lq = ['violet', 'fuchsia', 'purple']
style_lq = ['-', '-.', ':']
for i in range(0, len(norm_count_a_lq)):
    label_a_lq = ('LQ - Mass: ' + str(int(norm_count_a_lq[i][0])) + 
                  ' - Coupling: ' + str(norm_count_a_lq[i][1]))
    plt.plot(norm_count_a_lq[i][3], norm_count_a_lq[i][2], marker='.', 
             ls=style_lq[i], label=label_a_lq, color=color_lq[i])

font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

plt.tight_layout()
plt.legend(fontsize=8, loc=9)
start_name = 'diagrams/above-below_'
#fig_name_a = start_name + 'coupling2/diff_en_above_' + str(coupling) + '.png'
#fig_name_a = start_name + 'coupling2/diff_en_above_' + str(coupling) + '_1.png'
#fig_name_a = start_name + 'mass2/diff_en_above_' + str(int(mass_lq)) + '.png'
#fig_name_a = start_name + 'mass2/diff_en_above_' + str(int(mass_lq)) + '_1.png'
#plt.savefig(fig_name_a, dpi=200)
plt.show()

###############################################################################
#BELOW DIAGRAM

y_events_er_b_sm = []
y_errorbar_b_sm  = []
y_errorbar_b_sm1 = []
y_errorbar_b_sm2 = []
#error propagation of x/a is x*x_err/a**2
for i in range(0, len(bin_centers_b_sm)):
    if counts_b_sm[i] == 0.0:
        count_val_b = 1e-5
    else:
        count_val_b = counts_b_sm[i]
    y_events_er_b_sm.append( n_events_er_b_sm[i] / count_val_b )
    y_errorbar_b_sm1.append( yerr_b_sm[0][i] / count_val_b )
    y_errorbar_b_sm2.append( yerr_b_sm[1][i] / count_val_b )
y_errorbar_b_sm = [np.array(y_errorbar_b_sm1), np.array(y_errorbar_b_sm2)]   


# Plot below histogram
fit_b, ax_b = plt.subplots(figsize=(7, 5))
plt.title("Astro. Below")
plt.grid(True)
plt.xscale('log')
plt.xlim(10.0e3, 10.0e6)
#plt.xlim(10.0e3, 2.5e6)
plt.yscale('log') #'linear')
plt.ylim(1e-2, 1.0e2) 
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Normalized events per 2635 days")
plt.axhline(y=1, linestyle="dashed", label='Theoretical SM')

plt.errorbar(
    bin_centers_b_sm,
    y_events_er_b_sm, #n_events_er_b_sm,
    #xerr=xerr_b_sm,
    yerr=y_errorbar_b_sm, #yerr_b_sm,
    color="black",
    marker=".",#None,
    label="Data",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_b,
)

color_lq = ['violet', 'fuchsia', 'purple']
style_lq = ['-', '-.', ':']
for i in range(0, len(norm_count_b_lq)):
    label_b_lq = ('LQ - Mass: ' + str(int(norm_count_b_lq[i][0])) + 
                  ' - Coupling: ' + str(norm_count_b_lq[i][1]))
    plt.plot(norm_count_b_lq[i][3], norm_count_b_lq[i][2], marker='.', 
             ls=style_lq[i], label=label_b_lq, color=color_lq[i])

font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

plt.tight_layout()
plt.legend(fontsize=8)
start_name = 'diagrams/above-below_'
#fig_name_b = start_name + 'coupling2/diff_en_below_' + str(coupling) + '.png'
#fig_name_b = start_name + 'coupling2/diff_en_below_' + str(coupling) + '_1.png'
#fig_name_b = start_name + 'mass2/diff_en_below_' + str(int(mass_lq)) + '.png'
#fig_name_b = start_name + 'mass2/diff_en_below_' + str(int(mass_lq)) + '_1.png'
#plt.savefig(fig_name_b, dpi=200)
plt.show()


#eof