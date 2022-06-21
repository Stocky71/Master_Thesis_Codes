import numpy as np
import math

import autodiff as ad
import det_sys_weights
#import attenuation
import pickle as pkl


class Weighter:
    def __init__(self, mc):

        self.mc = mc
        # The SysWeighter class handles the weight corrections arising from varying detector
        # systematic parameters.
        self.sys_weighter = det_sys_weights.SysWeighter(self.mc)

#####NEW ADDED
    def set_lpmass_coupling(self, mass_leptoquark, coupling, type_z, ang=math.pi/2 ):
        self.mass_leptoquark = mass_leptoquark
        self.coupling        = coupling
        
        file_name = ('/datos_diff_en/attenuation/att' + str(mass_leptoquark) + 
                     '_' + str(coupling) + '.pkl')
        base_path = ('c:/users/anton/desktop/máster/tfm/' + 
                     'hese-7-year-data-release-main/hese-7-year-data-release')
        #base_path = ('/home/antonio/TFM/HESE-7-year-data-release/' + 
        #             'HESE-7-year-data-release')
        full_filepath = base_path + file_name
        with open(full_filepath, 'rb') as f:
            arr = pkl._Unpickler(f)
            arr.encoding = 'latin1'
            arr = arr.load()
        if ang != math.pi/2:
            arr3 = arr[np.where(arr[:,1] >= ang)]
        else:
            if type_z == 'A': #Above
                arr3 = arr[np.where(arr[:,1] >= math.pi/2)]
            elif type_z == 'B': #Below
                arr3 = arr[np.where(arr[:,1] <  math.pi/2)]
        
        self.attenuation = arr3[:,2:]
        #self.attenuation = attenuation.Attenuation(self.mass_leptoquark, 
        #                                           self.coupling)
        #self.attenuation.set_interpolation()
#####NEW ADDED    

    def flux_power_law(self, energy, zenith, norm, gamma, pivot): #ADDED ZENITH
        #NEW ADDED
        #attenuation = self.attenuation.get_attenuation(energy, zenith)
        #NEW ADDED
        e_scale = energy / pivot
        spectrum = ad.pow_r(e_scale, ad.mul(gamma, -1.0))
        flux = ad.mul_grad(norm, spectrum)
        #NEW ADDED
        #flux = ad.mul_r(flux, attenuation)
        flux_new = flux
        for i in range(0, len(self.attenuation)):
            flux_new[0][i] = flux[0][i] * self.attenuation[i][0]
            flux_new[1][i] = flux[1][i] * self.attenuation[i][0]
        flux = flux_new
        #NEW ADDED
        return flux

    def flux_spl(self, mc, astro_norm, astro_gamma, pivot_point=1e5):
        energy = mc["primaryEnergy"]
        flux = self.flux_power_law(energy, mc["recoZenith"], astro_norm, astro_gamma, pivot_point)

        # astro_norm is the 6 neutrino normalization so we need to convert it to the flux for 1 neutrino
        astro_flux = 1e-18 / 6.0
        flux = ad.mul(flux, astro_flux)

        return flux

    def weight_nunubar_ratio(self, mc, nunubar_ratio):
        p_id = mc["primaryType"]

        weights = np.empty(len(mc))
        gradients = np.empty((len(mc), len(nunubar_ratio[1])))

        weights[p_id > 0], gradients[p_id > 0] = nunubar_ratio
        weights[p_id < 0], gradients[p_id < 0] = ad.minus_r(2.0, nunubar_ratio)

        return weights, gradients

    def flux_conv(
        self,
        mc,
        conv_norm,
        kpi_ratio,
        cr_delta_gamma,
        nunubar_ratio,
        pivot_point=2020.0,
    ):
        energy = mc["primaryEnergy"]
        tilt_flux = self.flux_power_law(energy, mc["recoZenith"], conv_norm, cr_delta_gamma, pivot_point)

        pion_flux = mc["pionFlux"]
        kaon_flux = mc["kaonFlux"]

        total_flux = ad.plus_r(pion_flux, ad.mul_r(kaon_flux, kpi_ratio))

        flux = ad.mul_grad(total_flux, tilt_flux)
        flux = ad.mul_grad(flux, self.weight_nunubar_ratio(mc, nunubar_ratio))

        return flux

    def flux_prompt(
        self,
        mc,
        prompt_norm,
        cr_delta_gamma,
        nunubar_ratio,
        pivot_point=7887.0,
    ):
        energy = mc["primaryEnergy"]
        flux = self.flux_power_law(energy, mc["recoZenith"], prompt_norm, cr_delta_gamma, pivot_point)

        prompt_flux = mc["promptFlux"]

        flux = ad.mul(flux, prompt_flux)
        flux = ad.mul_grad(flux, self.weight_nunubar_ratio(mc, nunubar_ratio))

        return flux

    def weight_muon(self, mc, muon_norm):

        return ad.mul(muon_norm, mc["muonWeightOverLivetime"])

    def astro_detector_correction(self, epsilon_dom, epsilon_head_on, anisotropy_scale):

        hole_ice_weight = self.sys_weighter.get_hole_ice_weights(
            "Astro", epsilon_head_on
        )
        dom_eff_weight = self.sys_weighter.get_dom_eff_weights("Astro", epsilon_dom)
        anisotropy_weight = self.sys_weighter.get_anisotropy_weights(
            "Astro", anisotropy_scale
        )

        return ad.mul_grad(
            ad.mul_grad(hole_ice_weight, dom_eff_weight), anisotropy_weight
        )

    def conv_detector_correction(self, epsilon_dom, epsilon_head_on, anisotropy_scale):

        hole_ice_weight = self.sys_weighter.get_hole_ice_weights(
            "Conv", epsilon_head_on
        )
        dom_eff_weight = self.sys_weighter.get_dom_eff_weights("Conv", epsilon_dom)
        anisotropy_weight = self.sys_weighter.get_anisotropy_weights(
            "Conv", anisotropy_scale
        )

        return ad.mul_grad(
            ad.mul_grad(hole_ice_weight, dom_eff_weight), anisotropy_weight
        )

    def prompt_detector_correction(
        self, epsilon_dom, epsilon_head_on, anisotropy_scale
    ):

        hole_ice_weight = self.sys_weighter.get_hole_ice_weights(
            "Prompt", epsilon_head_on
        )
        dom_eff_weight = self.sys_weighter.get_dom_eff_weights("Prompt", epsilon_dom)
        anisotropy_weight = self.sys_weighter.get_anisotropy_weights(
            "Prompt", anisotropy_scale
        )

        return ad.mul_grad(
            ad.mul_grad(hole_ice_weight, dom_eff_weight), anisotropy_weight
        )

    def get_weights(self, livetime, parameter_names, params):
        """
        Return
        ---------
        tuple
            Zeroth element contains the list of weights
            First element contains the list of the gradients for each weight
        """

        n_params = len(params)

        p = dict()

        # Initialize parameter vector with gradient
        for i, (name, param) in enumerate(zip(parameter_names, params)):
            p_grad = np.zeros(shape=n_params).astype(float)
            p_grad[i] = 1.0
            p[name] = [param, p_grad]

        # Each element is a tuple. The zeroth element is the value, and the
        # first element is the corresponding gradient
        astro_norm = p["astro_norm"]
        astro_gamma = p["astro_gamma"]
        conv_norm = p["conv_norm"]
        prompt_norm = p["prompt_norm"]
        kpi_ratio = p["kpi_ratio"]
        cr_delta_gamma = p["cr_delta_gamma"]
        epsilon_dom = p["epsilon_dom"]
        epsilon_head_on = p["epsilon_head_on"]
        anisotropy_scale = p["anisotropy_scale"]
        nunubar_ratio = p["nunubar_ratio"]
        muon_norm = p["muon_norm"]

        # Calculate the expected neutrino flux from each component
        astro_fluxes = self.flux_spl(
            self.mc, astro_norm=astro_norm, astro_gamma=astro_gamma
        )
        conv_fluxes = self.flux_conv(
            self.mc,
            conv_norm=conv_norm,
            kpi_ratio=kpi_ratio,
            cr_delta_gamma=cr_delta_gamma,
            nunubar_ratio=nunubar_ratio,
        )
        prompt_fluxes = self.flux_prompt(
            self.mc,
            prompt_norm=prompt_norm,
            cr_delta_gamma=cr_delta_gamma,
            nunubar_ratio=nunubar_ratio,
        )

        # Calculate the muon weights
        muon_weights = self.weight_muon(self.mc, muon_norm=muon_norm)

        # Correct atmospheric weights with self-veto
        conv_fluxes = ad.mul(conv_fluxes, self.mc["conventionalSelfVetoCorrection"])
        prompt_fluxes = ad.mul(prompt_fluxes, self.mc["promptSelfVetoCorrection"])

        # Weight modifications due to detector systematics
        astro_weights = ad.mul_grad(
            astro_fluxes,
            self.astro_detector_correction(
                epsilon_dom, epsilon_head_on, anisotropy_scale
            ),
        )
        conv_weights = ad.mul_grad(
            conv_fluxes,
            self.conv_detector_correction(
                epsilon_dom, epsilon_head_on, anisotropy_scale
            ),
        )
        prompt_weights = ad.mul_grad(
            prompt_fluxes,
            self.prompt_detector_correction(
                epsilon_dom, epsilon_head_on, anisotropy_scale
            ),
        )

        neutrino_flux = ad.plus_grad(
            astro_weights, ad.plus_grad(conv_weights, prompt_weights)
        )

        neutrino_weights = ad.mul(neutrino_flux, self.mc["weightOverFluxOverLivetime"])

        weights = ad.plus_grad(neutrino_weights, muon_weights)

        weights = ad.mul(weights, livetime)

        return weights
