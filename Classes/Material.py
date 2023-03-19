import numpy as np
from typing import List

from Tools.csv_tool import read_csv
from matplotlib import colors


class PhaseChangeMaterial:
    def __init__(self, name: str, fusion_temp: float, rho_l: float, lambda_l: float, c_l: float, latent_heat: float,
                 rho_s: float,
                 lambda_s: float, c_s: float, solid_color: np.array, liquid_color: np.array):
        self.name = name
        # Material properties
        self.fusion_temp = fusion_temp
        self.rho_l = rho_l  # Volumetric mass density: kg/m^3
        self.lambda_l = lambda_l  # Thermal conductivity: W/(m.K)
        self.c_l = c_l  # Specific heat capacity: J/(kg.K)
        self.latent_heat = latent_heat  # Specific latent heat: J/kg
        self.rho_s = rho_s
        self.lambda_s = lambda_s
        self.c_s = c_s

        self.alpha_l = self.lambda_l / (self.rho_l * self.c_l)  # Thermal diffusivity: m^2/s
        self.alpha_s = self.lambda_s / (self.rho_s * self.c_s)
        self.n = self.rho_l * self.latent_heat  # Used in calculations

        self.solid_color = colors.to_rgba_array(solid_color)
        self.liquid_color = colors.to_rgba_array(liquid_color)


def create_pcm_list_from_csv(csv_path: str) -> List[PhaseChangeMaterial]:
    pcm_list = []
    pcm_settings = read_csv(csv_path)
    for pcm in pcm_settings:
        pcm_list.append(PhaseChangeMaterial(name=pcm[0],
                                            fusion_temp=float(pcm[1]) + 273.15,
                                            rho_l=float(pcm[2]),
                                            lambda_l=float(pcm[3]),
                                            c_l=float(pcm[4]),
                                            latent_heat=float(pcm[5]),
                                            rho_s=float(pcm[6]),
                                            lambda_s=float(pcm[7]),
                                            c_s=float(pcm[8]),
                                            solid_color=str(pcm[9]),
                                            liquid_color=str(pcm[10])))
    return pcm_list
