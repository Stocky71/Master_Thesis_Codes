import numpy  as np
import pickle as pkl
import matplotlib.pyplot as plt
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

file_name_3 = 'datos_diff_en/SM/above_diag_data.pkl'
with open(file_name_3, 'rb') as f:
    arr_dt_a = pkl._Unpickler(f)
    arr_dt_a.encoding = 'latin1'
    arr_dt_a = arr_dt_a.load()
    
mc_a_sm       = arr_dt_a[:236365]   
weights_a_sm  = [arr_dt_a[236365], arr_dt_a[236366], arr_dt_a[236367]] 
n_events_a_sm = arr_dt_a[236368:]

counts_a_sm = np.histogram(len(weights_a_sm) * [mc_a_sm], weights=weights_a_sm,
                           bins=n_events_a_sm)[0]

#counts_a_sm = n_a[0][2]
counts_a_py = counts_a_sm / 7.5 #counts_a_sm / 7.5

counts_a_10  = []
counts_a_20  = []
counts_a_50  = []
counts_a_100 = []

error_a_10  = []
error_a_20  = []
error_a_50  = []
error_a_100 = []

for i in range(0, len(counts_a_py)):
    counts_a_10.append(np.random.poisson(lam=counts_a_py[i]*10, size=1)[0])
    counts_a_20.append(np.random.poisson(lam=counts_a_py[i]*20, size=1)[0])
    counts_a_50.append(np.random.poisson(lam=counts_a_py[i]*50, size=1)[0])
    counts_a_100.append(np.random.poisson(lam=counts_a_py[i]*100, size=1)[0])
    
    error_a_10.append(np.sqrt(counts_a_10[i]))
    error_a_20.append(np.sqrt(counts_a_20[i]))
    error_a_50.append(np.sqrt(counts_a_50[i]))
    error_a_100.append(np.sqrt(counts_a_100[i]))

###############################################################################
# LOAD LQ DATA FOR PLOT
norm_count_a_lq = []

mass_ran = list(np.linspace(500, 1500, 3))
coupling_ran = [1.0] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]
        
        norm_count_a_lq_comp = []
        
        arr_er_a = []
        arr_dt_a = []
        
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
        
        file_name_3 = ('datos_diff_en/lq_' + str(int(mass_lq)) + 
                       '/above_diag_data-' + str(coupling) + '.pkl')
        with open(file_name_3, 'rb') as f:
            arr_dt_a = pkl._Unpickler(f)
            arr_dt_a.encoding = 'latin1'
            arr_dt_a = arr_dt_a.load()
            
        mc_a_lq       = arr_dt_a[:236365]   
        weights_a_lq  = [arr_dt_a[236365], arr_dt_a[236366], arr_dt_a[236367]] 
        n_events_a_lq = arr_dt_a[236368:]

        
        counts_a_lq = np.histogram(3*[mc_a_lq], weights=weights_a_lq,
                                   bins=n_events_a_lq)[0]

        norm_count_a_lq_comp = []
        #error propagation of x/a is x*x_err/a**2
        for k in range(0, len(bin_centers_a_lq)):
            if counts_a_sm[k] == 0.0:
                norm_count_a_lq_comp.append( counts_a_lq[k] / 1e-5 )
            else:
                norm_count_a_lq_comp.append( counts_a_lq[k] / counts_a_sm[k] )
            
                
        norm_count_a_lq.append((mass_lq, coupling, norm_count_a_lq_comp, 
                                bin_centers_a_lq))

###############################################################################
###PLOT ABOVE DIAGRAM - 10 YEARS
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.title('Simulated above events for coming years')
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1e-2, 1.0e2)
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Normalized events per 2635 days")
plt.grid(axis='both', alpha=0.6)
plt.axhline(y=1, linestyle="dashed", label='Theoretical SM')
counts_a = counts_a_10
counts_a_err = error_a_10


sm_data_norm  = []
sm_dter1_norm = []
for i in range(0, len(bin_centers_a_sm)):
    if counts_a_py[i]==0:
        count_val = 1e-5
    else:
        count_val = counts_a_py[i]
    count_val = count_val * 10
    sm_data_norm.append( counts_a[i] / count_val )
    sm_dter1_norm.append( counts_a_err[i] / count_val )
sm_dter_norm = [np.array(sm_dter1_norm), np.array(sm_dter1_norm)] 


plt.errorbar(
    bin_centers_a_sm,
    sm_data_norm, #n_events_er_a_sm,
    #xerr=xerr_a_sm,
    yerr=sm_dter_norm, #yerr_a_sm,
    color="black",
    marker=".",#None,
    label="Forecast Data - 10 years (1 year Gen2)",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_a,
)

color_lq = ['purple', 'mediumvioletred', 'fuchsia', 'violet']
style_lq = ['-', '--', '-.', ':']
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
plt.legend(fontsize=8, loc=2)
fig_name = 'diagrams/above_diag_future1-10.png'
plt.savefig(fig_name, dpi=200)
plt.show()

###############################################################################
###PLOT ABOVE DIAGRAM - 100 YEARS
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.title('Simulated above events for coming years')
plt.loglog()
plt.xlim(10.0e3, 10.0e6)
plt.ylim(1e-2, 1.0e2) 
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Normalized events per 2635 days")
plt.grid(axis='both', alpha=0.6)
plt.axhline(y=1, linestyle="dashed", label='Theoretical SM')
counts_a1 = counts_a_100
counts_a1_err = error_a_100


sm_data_norm1  = []
sm_dter1_norm1 = []
for i in range(0, len(bin_centers_a_sm)):
    if counts_a_py[i]==0:
        count_val = 1e-5
    else:
        count_val = counts_a_py[i]
    count_val = count_val * 100
    sm_data_norm1.append( counts_a1[i] / count_val )
    sm_dter1_norm1.append( counts_a1_err[i] / count_val )
sm_dter_norm = [np.array(sm_dter1_norm1), np.array(sm_dter1_norm1)] 


plt.errorbar(
    bin_centers_a_sm,
    sm_data_norm1, #n_events_er_a_sm,
    #xerr=xerr_a_sm,
    yerr=sm_dter1_norm1, #yerr_a_sm,
    color="black",
    marker=".",#None,
    label="Forecast Data - 100 years (10 years Gen2)",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_a,
)

color_lq = ['purple', 'mediumvioletred', 'fuchsia', 'violet']
style_lq = ['-', '--', '-.', ':']
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
plt.legend(fontsize=8, loc=2)
fig_name = 'diagrams/above_diag_future1-100.png'
plt.savefig(fig_name, dpi=200)
plt.show()


#eof    