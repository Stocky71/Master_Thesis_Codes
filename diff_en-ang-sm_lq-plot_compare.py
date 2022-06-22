import numpy  as np
import pickle as pkl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

#ANGLE
ang = 100 

# Load data/MC. By default load_mc loads events at energies >60 TeV, but we want to plot all events.
###############################################################################
# LOAD SM DATA FOR PLOT
#### ANG-100 -> 80
file_name_1 = 'datos_diff_en/SM/above_diag_errorbar-ang_' + str(ang) + '.pkl'
with open(file_name_1, 'rb') as f:
    arr_er_sm = pkl._Unpickler(f)
    arr_er_sm.encoding = 'latin1'
    arr_er_sm = arr_er_sm.load()
    
xerr_sm        = [arr_er_sm[0], arr_er_sm[1]]   
yerr_sm        = [arr_er_sm[2], arr_er_sm[3]]  
bin_centers_sm = arr_er_sm[4:33]
n_events_er_sm = arr_er_sm[33:]

file_name_2 = 'datos_diff_en/SM/above_diag_data-ang_' + str(ang) + '.pkl'
with open(file_name_2, 'rb') as f:
    arr_dt_sm = pkl._Unpickler(f)
    arr_dt_sm .encoding = 'latin1'
    arr_dt_sm = arr_dt_sm .load()
    
len_dt = len(arr_dt_sm) #100º->183902; 110º->

mc_sm       = arr_dt_sm [:len_dt-33]   
weights_sm  = [arr_dt_sm [len_dt-33], arr_dt_sm [len_dt-32], arr_dt_sm [len_dt-31]] 
n_events_sm = arr_dt_sm [len_dt-30:]
    
#counts_sm     = np.histogram(mc_sm, weights=weights_sm[2],
#                             bins=n_events_sm)[0]
counts_sm     = np.histogram(3*[mc_sm], weights=weights_sm,
                             bins=n_events_sm)[0]


###############################################################################
# LOAD LQ DATA FOR PLOT
norm_count_lq = []

mass_ran = list(np.linspace(500, 500, 1))
coupling_ran = [0.5, 1.0, 1.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]
        
        norm_count_lq_comp = []
        arr_er = []
        arr_dt = []
        
        file_name_1 = ('datos_diff_en/ang_' + str(ang) + '/lq_' + 
                       str(int(mass_lq)) + '/above_diag_errorbar-' + 
                       str(coupling) + '.pkl')
        with open(file_name_1, 'rb') as f:
            arr_er = pkl._Unpickler(f)
            arr_er.encoding = 'latin1'
            arr_er = arr_er.load()
            
        xerr_lq        = [arr_er[0], arr_er[1]]   
        yerr_lq        = [arr_er[2], arr_er[3]]  
        bin_centers_lq = arr_er[4:33]
        n_events_er_lq = arr_er[33:]

        
        file_name_2 = ('datos_diff_en/ang_' + str(ang) + '/lq_' + 
                       str(int(mass_lq)) + '/above_diag_data-' + 
                       str(coupling) + '.pkl')
        with open(file_name_2, 'rb') as f:
            arr_dt = pkl._Unpickler(f)
            arr_dt.encoding = 'latin1'
            arr_dt = arr_dt.load()
            
        len_dt_lq = len(arr_dt) #100º->183902; 110º->
        
        mc_lq       = arr_dt[:len_dt-33]   
        weights_lq  = [arr_dt[len_dt-33], arr_dt[len_dt-32], arr_dt[len_dt-31]] 
        n_events_lq = arr_dt[len_dt-30:]
        
        
        #counts_lq = np.histogram(mc_lq, weights=weights_lq[2],
        #                         bins=n_events_lq)[0]

        counts_lq = np.histogram(3*[mc_lq], weights=weights_lq,
                                 bins=n_events_lq)[0]
        
        norm_count_lq_comp = []
        #error propagation of x/a is x*x_err/a**2
        for k in range(0, len(bin_centers_lq)):
            if counts_sm[k] == 0.0:
                norm_count_lq_comp.append( counts_lq[k] / 1e-6 )
            else:
                norm_count_lq_comp.append( counts_lq[k] / counts_sm[k] )
        

        norm_count_lq.append((mass_lq, coupling, norm_count_lq_comp, 
                              bin_centers_lq))

###############################################################################
#ABOVE DIAGRAM

y_events_er_sm = []
y_errorbar_sm  = []
y_errorbar_sm1 = []
y_errorbar_sm2 = []
#error propagation of x/a is x*x_err/a**2
for i in range(0, len(bin_centers_sm)):
    if counts_sm[i] == 0.0:
        count_val = 1e-5
    else:
        count_val = counts_sm[i]
    y_events_er_sm.append( n_events_er_sm[i] / count_val )
    y_errorbar_sm1.append( yerr_sm[0][i] / count_val )
    y_errorbar_sm2.append( yerr_sm[1][i] / count_val )
y_errorbar_sm = [np.array(y_errorbar_sm1), np.array(y_errorbar_sm2)]   


y_sm_data = []
y_sm_mc   = []
for i in range(0, len(bin_centers_sm)):
    if counts_sm[i] == 0.0:
        count_val = 1e-5
    else:
        count_val = counts_sm[i]
    y_sm_data.append( n_events_er_sm[i] / count_val )
    y_sm_mc.append( counts_sm[i] / count_val )
    

# Plot above histogram
fit_a, ax_a = plt.subplots(figsize=(7, 5))
plt.title("Astro. Above -- Zenith >= " + str(ang) + "º")
plt.grid(True)
plt.xscale('log')
plt.xlim(10.0e3, 10.0e6)
#plt.xlim(10.0e3, 2.5e6)
plt.yscale('log')
#plt.ylim(3e-1, 0.5e1) 
plt.ylim(1e-2, 1.0e2) 
plt.xlabel("Deposited Energy [GeV]")
plt.ylabel("Normalized events per 2635 days")
plt.axhline(y=1, linestyle="dashed", label='Theoretical SM')

plt.errorbar(
    bin_centers_sm,
    y_events_er_sm, #n_events_er_a_sm,
    #xerr=xerr_a_sm,
    yerr=y_errorbar_sm, #yerr_a_sm,
    color="black",
    marker=".",#None,
    label="Data",
    linestyle="None",
    capsize=3,
    elinewidth=1,
#    ax=ax_a,
)

color_lq  = ['violet', 'fuchsia', 'purple']
style_lq  = ['-', '-.', ':']
marker_lq = ['+', 'x'] #['^', 's']
for i in range(0, len(norm_count_lq)):
    label_lq = ('LQ - Mass: ' + str(int(norm_count_lq[i][0])) + 
                  ' - Coupling: ' + str(norm_count_lq[i][1]))
    plt.plot(norm_count_lq[i][3], norm_count_lq[i][2], marker='.', 
             ls=style_lq[i], label=label_lq, color=color_lq[i])

font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")


plt.tight_layout()

plt.legend(fontsize=8, loc=9)
start_name = 'diagrams/above-below_'
mid_name   = '/diff_en_above_'
end_name   = '-ang_' + str(ang) + '.png'
end_name1  = '_1-ang_' + str(ang) + '.png'
#fig_name = start_name + 'coupling2' + mid_name + str(coupling) + end_name
#fig_name = start_name + 'coupling2' + mid_name + str(coupling) + end_name1
#fig_name = start_name + 'mass2' + mid_name + str(int(mass_lq)) + end_name
#fig_name = start_name + 'mass2' + mid_name + str(int(mass_lq)) + end_name1
#plt.savefig(fig_name, dpi=200)
plt.show()

#eof