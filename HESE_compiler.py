# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 21:32:18 2022

@author: anton
"""

import sys
import argparse
import numpy as np
import HESE_fit_comp

def obtain_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cr_delta_gamma", default=-0.05, type=float, help="set initial cosmic ray slope"
    )
    parser.add_argument(
        "--nunubar_ratio",
        default=1.0,
        type=float,
        help="set initial neutrino/antineutrino ratio",
    )
    parser.add_argument(
        "--anisotropy_scale",
        default=1.0,
        type=float,
        help="set initial ice anisotropy scale",
    )
    parser.add_argument(
        "--astro_gamma",
        default=2.5,
        type=float,
        help="set initial astrophysical spectral index",
    )
    parser.add_argument(
        "--astro_norm",
        default=6.0,
        type=float,
        help="set initial astrophysical six-neutrino flux normalization",
    )
    parser.add_argument(
        "--conv_norm",
        default=1.0,
        type=float,
        help="set initial atmospheric conventional neutrino flux normalization",
    )
    parser.add_argument(
        "--epsilon_dom",
        default=0.99,
        type=float,
        help="set initial DOM absolute energy scale",
    )
    parser.add_argument(
        "--epsilon_head_on",
        default=0.0,
        type=float,
        help="set initial DOM angular response",
    )
    parser.add_argument(
        "--muon_norm",
        default=1.0,
        type=float,
        help="set initial atmospheric muon flux normalization",
    )
    parser.add_argument(
        "--kpi_ratio",
        default=1.0,
        type=float,
        help="set initial kaon/pion ratio correction",
    )
    parser.add_argument(
        "--prompt_norm",
        default=1.0,
        type=float,
        help="set initial atmospheric prompt neutrino flux normalization",
    )

    parser.add_argument(
        "--fix_cr_delta_gamma", action="store_true", help="fix cosmic ray slope in fit"
    )
    parser.add_argument(
        "--fix_nunubar_ratio",
        action="store_true",
        help="fix neutrino/antineutrino ratio in fit",
    )
    parser.add_argument(
        "--fix_anisotropy_scale",
        action="store_true",
        help="fix ice anisotropy scale in fit",
    )
    parser.add_argument(
        "--fix_astro_gamma",
        action="store_true",
        help="fix astrophysical spectral index in fit",
    )
    parser.add_argument(
        "--fix_astro_norm",
        action="store_true",
        help="fix astrophysical six-neutrino flux normalization in fit",
    )
    parser.add_argument(
        "--fix_conv_norm",
        action="store_true",
        help="fix atmospheric conventional neutrino flux normalization in fit",
    )
    parser.add_argument(
        "--fix_epsilon_dom",
        action="store_true",
        help="fix DOM absolute energy scale in fit",
    )
    parser.add_argument(
        "--fix_epsilon_head_on", action="store_true", help="fix DOM angular response in fit"
    )
    parser.add_argument(
        "--fix_muon_norm",
        action="store_true",
        help="fix atmospheric muon flux normalization in fit",
    )
    parser.add_argument(
        "--fix_kpi_ratio", action="store_true", help="fix kaon/pion ratio correction in fit"
    )
    parser.add_argument(
        "--fix_prompt_norm",
        action="store_true",
        help="fix atmospheric prompt neutrino flux normalization in fit",
    )

    args = parser.parse_args()  
    
    return args

###############################################################################
###############################################################################

args = obtain_parser()

astro_gamma_min = 2.0
astro_gamma_int = 0.04
astro_norm_min = 0
astro_norm_int = 0.3

for i in range(0, 51):
    astro_gamma = astro_gamma_min + i*astro_gamma_int
    for j in range(0, 41):
        astro_norm = astro_norm_min + j*astro_norm_int
        
        args.astro_gamma = astro_gamma
        args.astro_norm  = astro_norm
        args.fix_astro_gamma = True
        args.fix_astro_norm  = True
        
        HESE_fit_comp.HESE_fit_comp(args)
    
    
#EOF