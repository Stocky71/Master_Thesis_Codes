import numpy as np
import math
import scipy.integrate as ig
import pickle as pkl
import lhapdf as l

###############################################
#CONSTANTS
###############################################
r_earth = 6367.0e+5    #cm

g_fermi = 1.1663787e-5 #GeV^-2
mass_w  = 80.379       #GeV
mass_z  = 91.1876      #GeV
mass_pr = 0.9382720813 #GeV
mass_ne = 0.9395654133 #GeV
mass_n  = mass_ne      #CHOOSE PROTON OR NEUTRON MASS!!!!

mass_pr_g = 1.67262192369e-24 #mass in g
mass_ne_g = 1.67492749804e-24 #mass in g
mass_n_g  = mass_ne_g         #CHOOSE PROTON OR NEUTRON MASS!!!!

#Conversion GeV^-2 to cm^2
conv_gev_cm = 0.389379e-27
    
#Weinberg's angle -> sin_w_angle = sin^2 W angle = 1-(mass_w/mass_z)^2
sin_w_angle = 1 - (mass_w / mass_z)**2

#Coefficient constant g^2/(1-sin^2 W angle)
const1 = (8 * g_fermi * mass_w**2) / (np.sqrt(2) * (1 - sin_w_angle))

#Coupling x between quarks and neutrinos
coupling = 1.0 #Default coupling, it can be changed

#Leptoquark mass
mass_lq = 1100 #GeV -- Default lq mass, it's the lower limit

#Other constants
L_f_u = -(2/3) * sin_w_angle + 0.5
L_f_d =  (1/3) * sin_w_angle - 0.5
R_f_u = -(2/3) * sin_w_angle
R_f_d =  (1/3) * sin_w_angle

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

def obtain_s(energy):
    return 2 * mass_n * energy

def obtain_y_factor(x, Q2, energy):#Obtain the (1-y)**2 factor
    fact = Q2 / (2 * x * mass_n * energy)
    z = 1 - fact
    return z**2

###############################################
#START --- OBTAIN COEFFICIENTS
###############################################
##Structure of functions: c_obtain_x_quark or c_obtain_x_quark_anti 
##(+lq for leptoquark part) with x=(a or b), quark=(u,d,c,s,b)

#Obtain a for up, charm and anti's
def c_obtain_a_up(x, Q2):
    a = (const1 * L_f_u) / (Q2 + mass_z**2)
    return a

def c_obtain_a_up_anti(x, Q2):
    return c_obtain_a_up(x, Q2)

def c_obtain_a_charm(x, Q2):
    return c_obtain_a_up(x, Q2)

def c_obtain_a_charm_anti(x, Q2):
    return c_obtain_a_up_anti(x, Q2)

#Obtain b for up, charm and anti's
def c_obtain_b_up(x, Q2):
    b = (const1 * R_f_u) / (Q2 + mass_z**2)
    return b

def c_obtain_b_up_anti(x, Q2):
    return c_obtain_b_up(x, Q2)

def c_obtain_b_charm(x, Q2):
    return c_obtain_b_up(x, Q2)

def c_obtain_b_charm_anti(x, Q2):
    return c_obtain_b_up_anti(x, Q2)

#Obtain a for down, strange and bottom and anti's
def c_obtain_a_down(x, Q2):
    a = (const1 * L_f_d) / (Q2 + mass_z**2)
    return a

def c_obtain_a_down_anti(x, Q2):
    return c_obtain_a_down(x, Q2)

def c_obtain_a_strange(x, Q2):
    return c_obtain_a_down(x, Q2)

def c_obtain_a_strange_anti(x, Q2):
    return c_obtain_a_down_anti(x, Q2)

def c_obtain_a_bottom(x, Q2):
    return c_obtain_a_down(x, Q2)

def c_obtain_a_bottom_anti(x, Q2):
    return c_obtain_a_down_anti(x, Q2)

#Obtain b for down, strange and bottom and anti's
def c_obtain_b_down(x, Q2):
    b = (const1 * R_f_d) / (Q2 + mass_z**2)
    return b

def c_obtain_b_down_anti(x, Q2):
    return c_obtain_b_down(x, Q2)

def c_obtain_b_down_lq2(x, energy):
    a = (x * obtain_s(energy) - mass_lq**2)**2 
    a = a + (mass_lq**2 * coupling**2 / (16.0 * math.pi))**2
    b = coupling**2 / a   
    return b

def c_obtain_b_down_anti_lq2(x, Q2, energy):
    a = (Q2 - x * obtain_s(energy) - mass_lq**2)**2
    #a = a + (mass_lq**2 * coupling**2 / (16.0 * math.pi))**2
    b = coupling**2 / a
    return b

def c_obtain_b_strange(x, Q2):
    return c_obtain_b_down(x, Q2)

def c_obtain_b_strange_anti(x, Q2):
    return c_obtain_b_down_anti(x, Q2)

def c_obtain_b_strange_lq2(x, energy):
    return c_obtain_b_down_lq2(x, energy)

def c_obtain_b_strange_anti_lq2(x, Q2, energy):
    return c_obtain_b_down_anti_lq2(x, Q2, energy)

def c_obtain_b_bottom(x, Q2):
    return c_obtain_b_down(x, Q2)

def c_obtain_b_bottom_anti(x, Q2):
    return c_obtain_b_down_anti(x, Q2)

def c_obtain_b_bottom_lq2(x, energy):
    return c_obtain_b_down_lq2(x, energy)

def c_obtain_b_bottom_anti_lq2(x, Q2, energy):
    return c_obtain_b_down_anti_lq2(x, Q2, energy)
###############################################
#END --- OBTAIN COEFFICIENTS
###############################################

###############################################
#START --- OBTAIN CROSS SECTION
###############################################
#quark: 1 to 5 -> (d,u,s,c,b)
#neutrino: 0 to cs neutrino contribution, 1 to antineutrino 
def quark_cont_cs_sm(x, Q2, energy, quark, neutrino):
    cte1 = obtain_y_factor(x, Q2, energy) #Obtain factor (1-y)**2
    
    quark    = int(quark)
    neutrino = int(neutrino)
    
    if quark == 1: #1->down
        cs_sm   = c_obtain_a_down(x, Q2)**2 
        cs_sm   = cs_sm + c_obtain_b_down(x, Q2)**2 * cte1
        cs_sm_a = c_obtain_b_down_anti(x, Q2)**2 
        cs_sm_a = cs_sm_a + c_obtain_a_down_anti(x, Q2)**2 * cte1
    elif quark == 2: #2->up
        cs_sm   = c_obtain_a_up(x, Q2)**2 
        cs_sm   = cs_sm + c_obtain_b_up(x, Q2)**2 * cte1
        cs_sm_a = c_obtain_b_up_anti(x, Q2)**2 
        cs_sm_a = cs_sm_a + c_obtain_a_up_anti(x, Q2)**2 * cte1
    elif quark == 3: #3->strange
        cs_sm   = c_obtain_a_strange(x, Q2)**2 
        cs_sm   = cs_sm + c_obtain_b_strange(x, Q2)**2 * cte1
        cs_sm_a = c_obtain_b_strange_anti(x, Q2)**2 
        cs_sm_a = cs_sm_a + c_obtain_a_strange_anti(x, Q2)**2 * cte1
    elif quark == 4: #4->charm
        cs_sm   = c_obtain_a_charm(x, Q2)**2 
        cs_sm   = cs_sm + c_obtain_b_charm(x, Q2)**2 * cte1
        cs_sm_a = c_obtain_b_charm_anti(x, Q2)**2 
        cs_sm_a = cs_sm_a + c_obtain_a_charm_anti(x, Q2)**2 * cte1
    elif quark == 5: #5->bottom
        cs_sm   = c_obtain_a_bottom(x, Q2)**2 
        cs_sm   = cs_sm + c_obtain_b_bottom(x, Q2)**2 * cte1
        cs_sm_a = c_obtain_b_bottom_anti(x, Q2)**2 
        cs_sm_a = cs_sm_a + c_obtain_a_bottom_anti(x, Q2)**2 * cte1
    
 
    if neutrino == 0:
        cs_sm   = cs_sm   * pdf.xfxQ2( quark, x, Q2) / x
        cs_sm_a = cs_sm_a * pdf.xfxQ2(-quark, x, Q2) / x
    elif neutrino == 1:
        cs_sm   = cs_sm   * pdf.xfxQ2(-quark, x, Q2) / x
        cs_sm_a = cs_sm_a * pdf.xfxQ2( quark, x, Q2) / x
                
    return (cs_sm + cs_sm_a)

def quark_cont_cs_lq(x, Q2, energy, quark, neutrino):
    cte1 = obtain_y_factor(x, Q2, energy) #Obtain factor (1-y)**2
    
    quark    = int(quark)
    neutrino = int(neutrino)
    
    cs_lq = cs_lq_a = 0
    if quark == 1: #1->down
        cs_lq   = c_obtain_b_down_lq2(x, energy) * cte1
        cs_lq_a = c_obtain_b_down_anti_lq2(x, Q2, energy) 
    elif quark == 2: #2->up
        pass
    elif quark == 3: #3->strange
        cs_lq   = c_obtain_b_strange_lq2(x, energy) * cte1
        cs_lq_a = c_obtain_b_strange_anti_lq2(x, Q2, energy) 
    elif quark == 4: #4->charm
        pass
    elif quark == 5: #5->bottom
        cs_lq   = c_obtain_b_bottom_lq2(x, energy) * cte1
        cs_lq_a = c_obtain_b_bottom_anti_lq2(x, Q2, energy) 

    
    if neutrino == 0:
        cs_lq   = cs_lq   * pdf.xfxQ2( quark, x, Q2) / x
        cs_lq_a = cs_lq_a * pdf.xfxQ2(-quark, x, Q2) / x
    elif neutrino == 1:
        cs_lq   = cs_lq   * pdf.xfxQ2(-quark, x, Q2) / x
        cs_lq_a = cs_lq_a * pdf.xfxQ2( quark, x, Q2) / x
        
    return (cs_lq + cs_lq_a)
 
#Part. ID in https://particle.wiki/wiki/PDG_particle_numbering_scheme
#neutrino: 0 to cs neutrino contribution, 1 to antineutrino 
def cross_section_sm(Q2, x, energy, neutrino):
    neutrino = int(neutrino)
    
    cs_sm =         quark_cont_cs_sm(x, Q2, energy, 1, neutrino) #q down
    cs_sm = cs_sm + quark_cont_cs_sm(x, Q2, energy, 2, neutrino) #q up
    cs_sm = cs_sm + quark_cont_cs_sm(x, Q2, energy, 3, neutrino) #q strange
    cs_sm = cs_sm + quark_cont_cs_sm(x, Q2, energy, 4, neutrino) #q charm
    cs_sm = cs_sm + quark_cont_cs_sm(x, Q2, energy, 5, neutrino) #q bottom

    return cs_sm

def cross_section_lq(Q2, x, energy, neutrino):
    neutrino = int(neutrino)
    
    cs_lq =         quark_cont_cs_lq(x, Q2, energy, 1, neutrino) #q down
    cs_lq = cs_lq + quark_cont_cs_lq(x, Q2, energy, 2, neutrino) #q up
    cs_lq = cs_lq + quark_cont_cs_lq(x, Q2, energy, 3, neutrino) #q strange
    cs_lq = cs_lq + quark_cont_cs_lq(x, Q2, energy, 4, neutrino) #q charm
    cs_lq = cs_lq + quark_cont_cs_lq(x, Q2, energy, 5, neutrino) #q bottom
    
    return cs_lq

def pre_cross_section_sm(x, energy, neutrino):
    lim_up = float(2 * mass_n * energy) #Upper limit for Q^2
    cs_sm = ig.quad(cross_section_sm, 0, x*lim_up, args=(x, energy, neutrino), 
                    limit=50, points=None)

    return cs_sm[0]

def pre_cross_section_lq(x, energy, neutrino):
    lim_up = float(2 * mass_n * energy) #Upper limit for Q^2
    cs_lq = ig.quad(cross_section_lq, 0, x*lim_up, args=(x, energy, neutrino), 
                    limit=50, points=None)
 
    return cs_lq[0]

def cs_only_sm(energy, neutrino):
    print('CROSS SECTION STANDARD MODEL', neutrino)
    cs_sm = ig.quad(pre_cross_section_sm, 0, 1, args=(energy, neutrino),
                     limit=50, points=None)
    #CS in GeV^-2 units, need conversion to cm^2 -> 1GeV^-2 = 0.389e-27 cm^2
    cs_sm = cs_sm[0] * conv_gev_cm / (32 * math.pi)

    return cs_sm

def cs_only_lq(energy, neutrino):
    print('CROSS SECTION LEPTOQUARK', neutrino)
    cs_lq = ig.quad(pre_cross_section_lq, 0, 1, args=(energy, neutrino),
                    limit=50, points=None)
 
    #CS in GeV^-2 units, need conversion to cm^2 -> 1GeV^-2 = 0.389e-27 cm^2
    cs_lq = cs_lq[0] * conv_gev_cm / (32 * math.pi)

    return cs_lq

#neutrino: 0 to cs neutrino contribution, 1 to antineutrino 
def cross_section(energy, neutrino):
    cs_sm = cs_only_sm(energy, neutrino)
    cs_lq = cs_only_lq(energy, neutrino)
    return [cs_sm, cs_lq]


##############################################################################
#GET CROSS SECTION SM
energy_ran = np.linspace(math.log10(10000.0), math.log10(100000000.0), 30)

print("###########CROSS SECTION SM############")
cross_sect_sm = [[0.0, 0.0]] #Energy, CS
for i in range(0, len(energy_ran)):
    energy_val = 10**energy_ran[i]
    print('Energia: ', i, '/', len(energy_ran), ' - ', energy_val, 'GeV')
    cs_sm = cs_only_sm(energy_val, 0) + cs_only_sm(energy_val, 1)
    print('Datos calculo SM (energia, cs_sm): ', [energy_val, cs_sm])
    cross_sect_aux = [[energy_val, cs_sm]]
    cross_sect_sm = np.vstack([cross_sect_sm, cross_sect_aux])

print("Acabado proceso energias SM")
cross_sect_sm = cross_sect_sm[1:, 0:]
file_name_sm = 'datos_cs/cross_section_sm.pkl'
with open(file_name_sm, 'wb') as f:
    pkl.dump(cross_sect_sm, f)
print("Archivo datos SM guardado")


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
            cs_lq = cs_only_lq(energy_val, 0) + cs_only_lq(energy_val, 1)
            print('Datos calculo LQ (energia, cs_lq): ', [energy_val, cs_lq])
            cross_sect_aux = [[energy_val, cs_lq]]
            cross_sect_lq = np.vstack([cross_sect_lq, cross_sect_aux])

        print("Acabado proceso energias LQ")
        cross_sect_lq = cross_sect_lq[1:, 0:]
        file_name_lq = ('datos_cs/cross_section_lq-' + str(mass_lq) + '_' +
                        str(coupling) + '.pkl')
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