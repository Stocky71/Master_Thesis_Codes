import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

###############################################################################
###############################################################################
#                          PLOT CS_LQ VS CS_LQ_CHECK
###############################################################################
###############################################################################
plt.grid(True)
plt.xlabel('Energy (GeV)')
plt.ylabel('Cross section (cm^2)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1.0e5, 10.0e6)
plt.ylim(1.0e-39, 1.0e-32)

color1 = ['violet', 'fuchsia', 'purple']
style1 = ['-', '-.', ':']

#Change mass range and coupling for different plots 
#Usually, 1 mass with 3 couplings or viceversa
mass_ran = list(np.linspace(500, 1500, 3)) #1500, 21))
coupling_ran = [1.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]
        
        file_name_lq = ('datos_cs/cross_section_lq-' + str(mass_lq) + '_' + 
                        str(coupling) + '.pkl')
        with open(file_name_lq, 'rb') as f:
            arr_lq = pkl._Unpickler(f)
            arr_lq.encoding = 'latin1'
            arr_lq = arr_lq.load()

        energy_lq = []
        cs_lq     = []
        for k in range(0, len(arr_lq)):
            energy_lq.append(arr_lq[k][0])
            cs_lq.append(arr_lq[k][1])
            
        
        file_name_lq_c = ('datos_cs_check/cross_section_lq-' + str(mass_lq) + 
                          '_' + str(coupling) + '.pkl')
        with open(file_name_lq_c, 'rb') as f:
            arr_lq_c = pkl._Unpickler(f)
            arr_lq_c.encoding = 'latin1'
            arr_lq_c = arr_lq_c.load()

        energy_lq_c = []
        cs_lq_c     = []
        for k in range(0, len(arr_lq_c)):
            energy_lq_c.append(arr_lq_c[k][0])
            cs_lq_c.append(arr_lq_c[k][1])
        
        label_lq   = 'CS_LQ - ' + str(mass_lq) + ', ' + str(coupling)
        plt.plot(energy_lq, cs_lq, label=label_lq, ls=style1[0], lw=1.0, 
                 marker='', mew=0.001, c=color1[i])#'g', mfc='r')
        
        label_lq_c = 'CS_LQ_C - ' + str(mass_lq) + ', ' + str(coupling)
        plt.plot(energy_lq_c, cs_lq_c, label=label_lq_c, ls=style1[1], lw=1.0, 
                 marker='', mew=0.001, c=color1[i]) #'y', mfc='b')
    
font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

plt.tight_layout()
plt.legend(fontsize=9)
#fig_name = 'diagrams/check/CS_check-' + str(int(mass_lq)) + '.png'
#fig_name = 'diagrams/check/CS_check-' + str(coupling) + '.png'
#plt.savefig(fig_name, dpi=200)
plt.show()

#eof