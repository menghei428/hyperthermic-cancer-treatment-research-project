# 1D Pennes bioheat model for superficial hyperthermia
# FTCS in time + flux-conservative diffusion in space

import numpy as np
from tissue_properties import get_tissue_properties


def run_simulation(T_bolus=40, SAR_multiplier=1.0,delta=0.01):
    # Global settings

    D = 0.04  # Total tissues depth in meters
    treatment_time = 3600  # second
    Nz = 201 # Number of spatial grid points
    z = np.linspace(0, D, Nz)
    dz = z[1] - z[0] # Spatial step size

    dt = 0.02 # Time step (seconds)
    Nt = int(treatment_time / dt)

    # Parameters

    rho_b = 1060 # blood density
    c_b = 3770 # blood specific heat
    T_body = 37

    Qapp_arr = np.zeros(Nz)      # applicator heating term [W/m^3]

    # Layer

    skin_end = 0.002    # 2 mm
    fat_end = 0.015     # 15 mm
    # muscle: remainder


    # Water bolus cooling

    h = 100  # W/(m^2 K),  heat transfer coefficient between bolus and skin. The literature suggested 70-152

    # Tissue properties

    k_arr, rho_arr, cp_arr, omega_arr, qm_arr = get_tissue_properties(z)

    rhoCp = rho_arr * cp_arr
    beta = (omega_arr * rho_b * c_b) / rhoCp
    metabolic = qm_arr / rhoCp


    Q0       = 2e4
    Qapp_arr = Q0 * SAR_multiplier * np.exp(-z / delta)
    Qapp_norm = Qapp_arr / rhoCp


    #Stability check

    alpha = k_arr / rhoCp
    stability = np.max(alpha) * dt / dz**2

    if stability > 0.5:
        raise ValueError(f"Stability condition violated: alpha*dt/dz^2 = {stability:.3f} > 0.5. Reduce dt or increase dz.")


    # Initial condition

    T = np.ones(Nz) * T_body # Initial temperature everywhere = body temperature
    T_new = T.copy()

    k_face = 2.0 * k_arr[:-1] * k_arr[1:] / (k_arr[:-1] + k_arr[1:])


    # Time stepping

    for _ in range(Nt):

        # Interior nodes: flux-conservative diffusion

        flux_r = k_face[1:] * (T[2:] - T[1:-1]) / dz
        flux_l = k_face[:-1] * (T[1:-1] - T[:-2]) / dz

        diffusion = dt * (flux_r - flux_l) / (rhoCp[1:-1] * dz)
        perfusion = -beta[1:-1] * dt * (T[1:-1] - T_body)

        T_new[1:-1] = T[1:-1] + diffusion + perfusion + metabolic[1:-1] * dt + Qapp_norm[1:-1] * dt

        # Surface boundary at z = 0:
        # -k dT/dz = h (T_surface - T_bolus)
        # Use ghost node with central difference
        k0 = k_arr[0]
        ghost = T[1] - 2.0 * dz * (h / k0) * (T[0] - T_bolus)

        diff0 = k0 * dt * (T[1] - 2.0 * T[0] + ghost) / (rhoCp[0] * dz**2)
        perf0 = -beta[0] * dt * (T[0] - T_body)

        T_new[0] = T[0] + diff0 + perf0 + metabolic[0] * dt + Qapp_norm[0] * dt
        T_new[-1] = T_body

        T[:] = T_new[:]


    # Post-processing

    z_mm = z * 1000

    # Therapeutic window mask (40–43 °C)
    mask = (T >= 40.0) & (T <= 43.0)

    if np.any(mask):
        z_min = z_mm[mask][0]
        z_max = z_mm[mask][-1]
        treatable_depth = z_max - z_min
    else:
        treatable_depth = 0.0

    T_max = np.max(T)

    return treatable_depth, T_max

if __name__ == "__main__":
    depth, T_max = run_simulation(T_bolus=40, SAR_multiplier=1.0,delta=0.01)
    print(f"Treatable depth: {depth:.2f} mm")
    print(f"Maximum temperature: {T_max:.2f} °C")