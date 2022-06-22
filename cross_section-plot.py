import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

#Only have one part of the code in use -> CS_SM VS CS_LQ or (CS_SM + CS_LQ) / CS_SM
#Also, inside them, there is the option to plot only one CS_LQ (usually commented) 
#or multiple

mass_lq_ini = 1000.0
mass_lq_fin = 1000.0
step        = 1
"""
###############################################################################
###############################################################################
#                          PLOT CS_SM VS CS_LQ
###############################################################################
###############################################################################
file_name_sm = 'datos_cs/cross_section_sm.pkl'
with open(file_name_sm, 'rb') as f:
    arr_sm = pkl._Unpickler(f)
    arr_sm.encoding = 'latin1'
    arr_sm = arr_sm.load()

energy_sm = []
cs_sm     = []
for i in range(0, len(arr_sm)):
    energy_sm.append(arr_sm[i][0])
    cs_sm.append(arr_sm[i][1])

color1 = ['violet', 'fuchsia', 'purple']
style1 = ['-', '-.', ':']
plt.grid(True)
plt.xlabel('Energy (GeV)')
plt.ylabel('Cross section (cm^2)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10.0e3, 10.0e7)
plt.plot(energy_sm, cs_sm, label='CS_SM', ls='--', lw=1.0, marker='', 
         mew=0.001, c='k') #, mfc='b')
"""
"""
####JUST ONE CS_LQ#############################################################
mass_lq = 1500.0
coupling = 10.0

file_name_lq = ('datos_cs/cross_section_lq-' + str(mass_lq) + '_' + 
                str(coupling) + '.pkl')
with open(file_name_lq, 'rb') as f:
    arr_lq = pkl._Unpickler(f)
    arr_lq.encoding = 'latin1'
    arr_lq = arr_lq.load()

energy_lq = []
cs_lq     = []
for i in range(0, len(arr_lq)):
    energy_lq.append(arr_lq[i][0])
    cs_lq.append(arr_lq[i][1])

plt.plot(energy_lq, cs_lq, label='CS_LQ', lw=1.5, marker='', mew=0.001, 
         c=color1[0]) #'g', mfc='r')
"""
"""
####EVERY CS_LQ################################################################
#Change mass range and coupling for different plots 
#Usually, 1 mass with 3 couplings or viceversa
mass_ran = list(np.linspace(mass_lq_ini, mass_lq_fin, step))
coupling_ran = [0.5, 1.0, 1.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
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
        
        label_lq = 'CS_LQ - ' + str(mass_lq) + ', ' + str(coupling)
        plt.plot(energy_lq, cs_lq, label=label_lq, lw=1.0, ls=style1[j],
                 marker='', mew=0.001, c=color1[j])#'g', mfc='r')


"""
###############################################################################
###############################################################################
#                      PLOT (CS_SM + CS_LQ) / CS_SM
###############################################################################
###############################################################################
file_name_sm = 'datos_cs/cross_section_sm.pkl'
with open(file_name_sm, 'rb') as f:
    arr_sm = pkl._Unpickler(f)
    arr_sm.encoding = 'latin1'
    arr_sm = arr_sm.load()

energy_sm = []
cs_sm     = []
for i in range(0, len(arr_sm)):
    energy_sm.append(arr_sm[i][0])
    cs_sm.append(arr_sm[i][1])

color1 = ['violet', 'fuchsia', 'purple']
style1 = ['-', '-.', ':']
plt.grid(True)
plt.xlabel('Energy (GeV)')
plt.ylabel('CS_TOT/CS_SM')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(10.0e3, 10.0e7)
"""
####JUST ONE CS_LQ#############################################################
mass_lq = 500.0
coupling = 10

file_name_lq = ('datos_cs/cross_section_lq-' + str(mass_lq) + '_' + 
                str(coupling) + '.pkl')
with open(file_name_lq, 'rb') as f:
    arr_lq = pkl._Unpickler(f)
    arr_lq.encoding = 'latin1'
    arr_lq = arr_lq.load()

energy_lq = []
cs_lq     = []
for i in range(0, len(arr_lq)):
    energy_lq.append(arr_lq[i][0])
    cs_lq.append(arr_lq[i][1])

x = []
y = []
for k in range(0, len(arr_lq)):
    x.append(arr_lq[k][0])
    cs_oper = (arr_lq[k][1] + arr_sm[k][1]) / (arr_sm[k][1])
    y.append(cs_oper)

label_lq = 'CS_TOT/CS_SM - ' + str(mass_lq) + ', ' + str(coupling)
plt.plot(x, y, label=label_lq, lw=1.5, marker=' ', mew=0.001, c=color1[0]) #'g', mfc='r')
"""
####EVERY CS_LQ################################################################
mass_lq_ini = 1000.0
mass_lq_fin = 1000.0
step        = 1

#Change mass range and coupling for different plots 
#Usually, 1 mass with 3 couplings or viceversa
mass_ran = list(np.linspace(mass_lq_ini, mass_lq_fin, step))
coupling_ran = [0.5, 1.0, 1.5] #[0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
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
        
        x = []
        y = []
        for k in range(0, len(arr_lq)):
            x.append(arr_lq[k][0])
            cs_oper = (arr_lq[k][1] + arr_sm[k][1]) / (arr_sm[k][1])
            y.append(cs_oper)

        label_lq = 'CS_TOT/CS_SM - ' + str(mass_lq) + ', ' + str(coupling)
        plt.plot(x, y, label=label_lq, lw=1.0, ls=style1[j], marker=' ', 
                 mew=0.001, c=color1[j]) #'g', mfc='r')


############################################################################### 
font = FontProperties()
font.set_size("medium")
font.set_family("sans-serif")
font.set_weight("bold")

plt.tight_layout()
plt.legend(fontsize=9)
#fig_name = 'diagrams/cs-sm_lq/cs-sm_lq-' + str(int(mass_lq)) + '.png'
#fig_name = 'diagrams/cs-sm_lq/cs-sm_lq-' + str(coupling) + '.png'
#fig_name = 'diagrams/cs-ratio/cs-ratio-' + str(int(mass_lq)) + '.png'
#fig_name = 'diagrams/cs-ratio/cs-ratio-' + str(coupling) + '.png'
#plt.savefig(fig_name, dpi=200)
plt.show()