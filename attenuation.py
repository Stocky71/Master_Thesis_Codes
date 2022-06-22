import numpy  as np
import math
import scipy.interpolate as ipol
import scipy.integrate   as ig
import pickle as pkl


class Attenuation:
    def __init__(self, mass_lq, coupling):
        #####CONSTANTS#####
        self.r_earth   = 6367.0e+5         #cm
        self.mass_pr_g = 1.67262192369e-24 #mass in g
        self.mass_ne_g = 1.67492749804e-24 #mass in g
        self.mass_n_g  = self.mass_ne_g    #CHOOSE PROTON OR NEUTRON MASS!!!!
        #####
        
        self.mass_lq   = round(float(mass_lq), 1)
        self.coupling  = coupling
        
        file_name = ('/datos_cs/cross_section_lq-' + str(mass_lq) + '_' + 
                     str(coupling) + '.pkl')
        base_path = ('c:/users/anton/desktop/máster/tfm/' + 
                     'hese-7-year-data-release-main/hese-7-year-data-release')
        
        full_filepath = base_path + file_name
        with open(full_filepath, 'rb') as f:
            arr = pkl._Unpickler(f)
            arr.encoding = 'latin1'
            arr = arr.load()
        
        self.arr = arr
    
    def set_interpolation(self):
        energy = []
        cs_lq  = []
        
        for i in range(0, len(self.arr)):
            energy.append(self.arr[i][0])
            cs_lq.append(self.arr[i][1])
            
        f_interp_lq = ipol.interp1d(energy, cs_lq, kind='quadratic')
        
        self.energy = energy
        self.cs_lq  = cs_lq
        self.f_interp_lq = f_interp_lq
        
    def get_interpolation_value(self, energy):
        return self.f_interp_lq(energy)
    
    def get_cross_section_value(self, energy):
        return self.get_interpolation_value(energy)
    
    def density_prof(self, s, zenith):
        x1 = 1 + (s/self.r_earth)**2 + 2 * (s/self.r_earth) * math.cos(zenith)
        x = math.sqrt(x1)
        
        if x >= 0 and x <= 0.191:
            density = 13.0885 - 8.8381 * x**2
        elif x > 0.191 and x <= 0.546:
            density = 12.5815 - 1.2638 * x - 3.6426 * x**2 - 5.5281 * x**3
        elif x > 0.546 and x <= 0.895:
            density =  7.9565 - 6.4761 * x + 5.5283 * x**2 - 3.0807 * x**3
        elif x > 0.895 and x <= 0.905:
            density =  5.3197 - 1.4836 * x
        elif x > 0.905 and x <= 0.937:
            density = 11.2494 - 8.0298 * x
        elif x > 0.937 and x <= 0.965:
            density =  7.1089 - 3.8045 * x
        elif x > 0.965 and x <= 0.996:
            density =  2.6910 + 0.6924 * x
        elif x > 0.996 and x <= 0.997:
            density =  2.900
        elif x > 0.997 and x <= 0.999:
            density =  2.600
        elif x > 0.999 and x <= 1.000:
            density =  1.020
        
        return density

    def obtain_trav_wght(self, zenith):
        l_max = -2 * self.r_earth * math.cos(zenith) #Units are in cm^-2   
        x = ig.quad(self.density_prof, 0, l_max, args=(zenith))[0]
        
        return (x / self.mass_n_g)
    
    def get_attenuation(self, energy, zenith):
        att = np.array([[0.0, 0.0, 0.0]])
        
        for i in range(0, len(energy)):    
            if zenith[i] < (math.pi/2):
                att = np.append(att, np.array([[energy[i], zenith[i], 1.0]]), 
                                axis=0)
            elif zenith[i] >= (math.pi/2):
                cross_sect = self.get_cross_section_value(energy[i])
                exp_value  = self.obtain_trav_wght(zenith[i]) * cross_sect
                att_val    = math.exp(- exp_value)
                att = np.append(att, np.array([[energy[i], zenith[i], att_val]]), 
                                axis=0)  
    	    #Units in exp -> cm^-2 * cm^2 = 1
        
        return att[1:, 0:]
                
#eof   