from srxraylib.plot.gol import set_qt
set_qt()


import numpy as np
import matplotlib.pyplot as plt


def get_phi_exact(beta, alpha, thetaD):
    cos_phi = np.cos(alpha) * np.cos(thetaD) + np.cos(beta) * np.sin(alpha) * np.sin(thetaD)
    # Clip to avoid numerical errors with arccos
    cos_phi = np.clip(cos_phi, -1, 1)
    return np.arccos(cos_phi)


def get_s0_s1(alpha, thetaD):
    cos_phi0 = np.cos(alpha) * np.cos(thetaD) + np.sin(alpha) * np.sin(thetaD)
    # s0 = arccos(cos(alpha-thetaD))^2 = (alpha-thetaD)^2
    s0 = (alpha - thetaD) ** 2

    # Calculation of s1 using the formula in the 2013 report (Eq 9)
    # s1 = sin(alpha) sin(thetaD) arccos(cos_phi0) / sqrt(1 - cos_phi0^2)
    sin_alpha = np.sin(alpha)
    sin_thetaD = np.sin(thetaD)

    if np.abs(alpha - thetaD) < 1e-10:
        s1 = sin_alpha * sin_thetaD  # Limit as u -> 1
    else:
        u = np.cos(alpha) * np.cos(thetaD) + np.sin(alpha) * np.sin(thetaD)
        s1 = (sin_alpha * sin_thetaD * np.arccos(u)) / np.sqrt(1 - u ** 2)
    return s0, s1


# Parameters: Graphite 002 @ 8 keV
thetaB = np.deg2rad(13.38257103)
# Symmetrical ray
alpha_sym = np.pi / 2 - thetaB
thetaD_sym = alpha_sym

# Slightly off-Bragg ray (e.g. 0.1 degree off)
alpha_off = np.pi / 2 - (thetaB + np.deg2rad(0.1))
thetaD_off = np.pi / 2 - thetaB  # thetaD is fixed by lambda and d-spacing

# Beta range (azimuthal angle)
# Typically tau is 0.1 deg, so we look at beta around tau' (approx 0.1 deg)
beta = np.linspace(-np.deg2rad(1.0), np.deg2rad(1.0), 1000)

# Symmetrical Case
s0_s, s1_s = get_s0_s1(alpha_sym, thetaD_sym)
phi_exact_s = get_phi_exact(beta, alpha_sym, thetaD_sym)
phi_app_s = np.sqrt(s0_s + s1_s * beta ** 2)

# Off-Bragg Case
s0_o, s1_o = get_s0_s1(alpha_off, thetaD_off)
phi_exact_o = get_phi_exact(beta, alpha_off, thetaD_off)
phi_app_o = np.sqrt(s0_o + s1_o * beta ** 2)

# Error calculation
error_s = phi_exact_s - phi_app_s
error_o = phi_exact_o - phi_app_o

# Plotting
plt.figure(figsize=(12, 10))

# Subplot 1: Exact vs Approx
plt.subplot(2, 1, 1)
plt.plot(np.rad2deg(beta), np.rad2deg(phi_exact_s), 'b-', label='Exact (Symmetrical)')
plt.plot(np.rad2deg(beta), np.rad2deg(phi_app_s), 'b--', label='Expansion (Symmetrical)')
plt.plot(np.rad2deg(beta), np.rad2deg(phi_exact_o), 'r-', label='Exact (0.1° Off)')
plt.plot(np.rad2deg(beta), np.rad2deg(phi_app_o), 'r--', label='Expansion (0.1° Off)')
plt.title('Comparison of Exact $\phi(\\beta)$ and Series Expansion')
plt.xlabel('Azimuthal Angle $\\beta$ (deg)')
plt.ylabel('Mosaic Angle $\phi$ (deg)')
plt.legend()
plt.grid(True)

# Subplot 2: Error
plt.subplot(2, 1, 2)
plt.plot(np.rad2deg(beta), np.rad2deg(error_s), 'b-', label='Error (Symmetrical)')
plt.plot(np.rad2deg(beta), np.rad2deg(error_o), 'r-', label='Error (0.1° Off)')
plt.title('Residual Error ($\phi_{exact} - \phi_{expansion}$)')
plt.xlabel('Azimuthal Angle $\\beta$ (deg)')
plt.ylabel('Error (deg)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('phi_error_check.png')
plt.show()

print(f"Symmetrical Case: s0={s0_s:.6e}, s1={s1_s:.6f}")
print(f"Off-Bragg Case: s0={s0_o:.6e}, s1={s1_o:.6f}")