import sys
import os
import os.path

base_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, base_path + "/resources/external/")
import numpy as np

import pickle as pkl
import data_loader
import attenuation


# Load data/MC. By default load_mc loads events at energies >60 TeV, but we want to plot all events.
mc_filenames = [
    "./resources/data/HESE_mc_observable.json",
    "./resources/data/HESE_mc_flux.json",
    "./resources/data/HESE_mc_truth.json",
]
mc = data_loader.load_mc(mc_filenames, emin=10.0e3)
data = data_loader.load_data("./resources/data/HESE_data.json", emin=10.0e3)


mass_ran = list(np.linspace(500, 1500, 11))
coupling_ran = [0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0]
for i in range(0, len(mass_ran)):
    mass_lq = mass_ran[i]
    for j in range(0, len(coupling_ran)):
        coupling = coupling_ran[j]
        
        print("PROCESO M_LQ:", mass_lq, ' - COUPLING:', coupling)
        att = attenuation.Attenuation(mass_lq, coupling)
        att.set_interpolation()

        att_val = att.get_attenuation(mc['primaryEnergy'], 
                                      mc['recoZenith'])

        file_name = ('datos_diff_en/attenuation/att' + str(mass_lq) + '_' + 
                     str(coupling) + '.pkl') 
        with open(file_name, 'wb') as f:
            pkl.dump(att_val, f)
            print("Archivo guardado")
            
#eof