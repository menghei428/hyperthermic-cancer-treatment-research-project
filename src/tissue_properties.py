"""
All values sourced from:
Haemmerich et al. — https://pmc.ncbi.nlm.nih.gov/articles/PMC12835268/table/Tab2/
"""

import numpy as np

# Layer boundary depths (m)
SKIN_END   = 0.002
FAT_END    = 0.015

#   k       thermal conductivity        [W / (m·K)]
#   rho     density                     [kg / m³]
#   cp      specific heat capacity      [J / (kg·K)]
#   omega   blood perfusion rate        [1 / s]
#   qm      metabolic heat generation   [W / m³]

LAYERS = {
    "skin": {
        "k":     0.42,
        "rho":   1109.0,
        "cp":    3500.0,
        "omega": 0.0022,
        "qm":    1620.0,
    },
    "fat": {
        "k":     0.25,
        "rho":   911.0,
        "cp":    2500.0,
        "omega": 0.00045,
        "qm":    300.0,
    },
    "muscle": {
        "k":     0.50,
        "rho":   1090.4,
        "cp":    3600.0,
        "omega": 0.001,
        "qm":    480.0,
    },
}


def get_tissue_properties(z_arr: np.ndarray) -> tuple:

    Nz = len(z_arr)

    k_arr     = np.empty(Nz)
    rho_arr   = np.empty(Nz)
    cp_arr    = np.empty(Nz)
    omega_arr = np.empty(Nz)
    qm_arr    = np.empty(Nz)

    for i, z in enumerate(z_arr):
        if z < SKIN_END:
            layer = LAYERS["skin"]
        elif z < FAT_END:
            layer = LAYERS["fat"]
        else:
            layer = LAYERS["muscle"]

        k_arr[i] = layer["k"]
        rho_arr[i] = layer["rho"]
        cp_arr[i] = layer["cp"]
        omega_arr[i] = layer["omega"]
        qm_arr[i] = layer["qm"]

    return k_arr, rho_arr, cp_arr, omega_arr, qm_arr