import numpy as np
import matplotlib.pyplot as plt
from sar_bolus_cooling import run_simulation

T_bolus = 41.0

penetration_depths = np.linspace(0.003, 0.02, 10)  # 3mm–20mm
sar_multipliers = np.linspace(0.5, 2.0, 10) # 10 steps between 50% and 200% power

thickness_results = np.zeros((len(penetration_depths), len(sar_multipliers)))
Tmax_results = np.zeros((len(penetration_depths), len(sar_multipliers)))

for i, delta in enumerate(penetration_depths):
    for j, sar in enumerate(sar_multipliers):

        thickness, Tmax = run_simulation(
            T_bolus=T_bolus,
            SAR_multiplier=sar,
            delta=delta
        )

        thickness_results[i, j] = thickness
        Tmax_results[i, j] = Tmax

        print(f"delta={delta*1000:.1f} mm, SAR={sar:.2f} -> Thickness: {thickness:.2f} mm, Tmax: {Tmax:.2f}")


plt.figure(figsize=(8, 6))

contour = plt.contourf(
    sar_multipliers,
    penetration_depths * 1000,
    thickness_results,
    cmap='inferno',
    levels=20
)

overheat_mask = Tmax_results > 43.0

plt.contourf(
    sar_multipliers,
    penetration_depths * 1000,
    overheat_mask,
    levels=[0.5, 1],
    colors='cyan',
    alpha=0.3
)


plt.colorbar(contour, label='Treatable Depth (mm)')

plt.xlabel('SAR Multiplier')
plt.ylabel('SAR Penetration Depth (mm)')
plt.title(f'Treatable Zone Depth (Bolus Temperature={T_bolus}°C)')

plt.text(
    0.55, 
    0.95, 
    "Cyan region: Tmax > 43°C (Overheating)",
    transform=plt.gca().transAxes,
    fontsize=9,
    verticalalignment='top',
    bbox=dict(facecolor='white', alpha=0.7)
)

plt.tight_layout()
plt.savefig("treatable_depth_map.png", dpi=300)
plt.show()