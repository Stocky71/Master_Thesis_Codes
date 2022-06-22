import numpy as np
import math
import pickle as pkl
import lhapdf as l

###############################################
#CONSTANTS
###############################################
mass_pr = 0.9382720813 #GeV
mass_ne = 0.9395654133 #GeV
mass_n  = mass_ne      #CHOOSE PROTON OR NEUTRON MASS!!!!

#Conversion GeV^-2 to cm^2
conv_gev_cm = 0.389379e-27

#Coupling x between quarks and neutrinos
coupling = 1.0 #Default coupling, it can be changed

#Leptoquark mass
mass_lq = 1100 #GeV -- Default lq mass, it's the lower limit

pdf = l.mkPDF("CT10nlo", 0)

###############################################
#FUNCTIONS
###############################################
def change_lq_mass(mass_leptoquark):
    global mass_lq 
    mass_lq = mass_leptoquark
    
def change_coupling(coupling_val):
    global coupling
    coupling = coupling_val

def get_x(energy):
    return mass_lq**2 / (2 * mass_n * energy)

###############################################
#START --- OBTAIN CROSS SECTION
###############################################
#quark: 1 to 5 -> (d,u,s,c,b)
#neutrino: 0 to cs neutrino contribution, 1 to antineutrino 
def cs_only_lq(energy, neutrino):
    print('CROSS SECTION LEPTOQUARK', neutrino)
    cs_lq = math.pi * coupling**2 / (2 * 2 * mass_n * energy)
    cs = 0
    x = get_x(energy)
    nu = mass_lq**2
    
    if neutrino == 0:
        cs_lq   = cs_lq * pdf.xfxQ2( 1, x, nu) / x
        cs_lq_a = cs_lq * pdf.xfxQ2(-1, x, nu) / x
    elif neutrino == 1:
        cs_lq   = cs_lq * pdf.xfxQ2(-1, x, nu) / x
        cs_lq_a = cs_lq * pdf.xfxQ2( 1, x, nu) / x

    cs = cs_lq + cs_lq_a

    #CS in GeV^-2 units, need conversion to cm^2 -> 1GeV^-2 = 0.389e-27 cm^2
    cs = cs * conv_gev_cm 

    return cs

##############################################################################
energy_ran = np.linspace(math.log10(100000.0), math.log10(100000000.0), 30)

print("###########CROSS SECTION LQ############")
mass_ran = list(np.linspace(500, 1500, 11))
coupling_ran = [0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    change_lq_mass(mass_ran[i])
    
    #THIS PART ONLY IF COUPLING WANTING TO CHANGE
    
    for j in range(0,len(coupling_ran)):
        change_coupling(coupling_ran[j])
        print("#### CS_LQ - m_lq: ", mass_lq, " - coupling: ", coupling)
        
        cross_sect_lq = [[0.0, 0.0]] #Energy, CS            
        for k in range(0, len(energy_ran)):
            energy_val = 10**energy_ran[k]
            print('Energia: ', k, '/', len(energy_ran), ' - ', energy_val, 'GeV')
            if get_x(energy_val) <= 1 and get_x(energy_val) >= 0:
                print('Calcula CS')
                cs_lq = cs_only_lq(energy_val, 0) + cs_only_lq(energy_val, 1)
                print('Datos calculo LQ (energia, cs_lq): ', [energy_val, cs_lq])
                cross_sect_aux = [[energy_val, cs_lq]]
                cross_sect_lq = np.vstack([cross_sect_lq, cross_sect_aux])
            else:
                print('No calcula CS')

        print("Acabado proceso energias LQ")
        cross_sect_lq = cross_sect_lq[1:, 0:]
        file_name_lq = ('datos_cs_check/cross_section_lq-' + str(mass_lq) + 
                        '_' + str(coupling) + '.pkl')
        with open(file_name_lq, 'wb') as f:
            pkl.dump(cross_sect_lq, f)
        print("Archivo datos LQ guardado")
        

    #THIS PART ONLY IF COUPLING IS CONSTANT    
"""
    print("#### CS_LQ - m_lq: ", mass_lq, " - coupling: ", coupling)
    
    cross_sect_lq = [[0.0, 0.0]] #Energy, CS            
    for k in range(0, len(energy_ran)):
        energy_val = 10**energy_ran[k]
        print('Energia: ', k, '/', len(energy_ran), ' - ', energy_val, 'GeV')
        cs_lq = cs_only_lq(energy_val, 0) + cs_only_lq(energy_val, 1)
        print('Datos calculo LQ (energia, cs_lq): ', [energy_val, cs_lq])
        cross_sect_aux = [[energy_val, cs_lq]]
        cross_sect_lq = np.vstack([cross_sect_lq, cross_sect_aux])

    print("Acabado proceso energias LQ")
    cross_sect_lq = cross_sect_lq[1:, 0:]
    file_name_lq = 'datos_cs/cross_section_lq-' + str(mass_lq) + '.pkl'
    with open(file_name_lq, 'wb') as f:
        pkl.dump(cross_sect_lq, f)
    print("Archivo datos LQ guardado")
"""    
    
#eof