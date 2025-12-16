from srxraylib.plot.gol import set_qt
set_qt()

import numpy as np
import matplotlib.pyplot as plt

# Graphite002 at 8 keV
eta_deg = 0.1 # mosaicity, degrees
eta_rad = np.deg2rad(eta_deg)
Q = 0.16693157449312368
mu = 9.6590899717931258  # cm^-1
t0 = 0.1 # cm
thetaB = np.deg2rad(13.38257103)
gamma0 = np.sin(thetaB)


def get_W(delta, eta, dist_type='gaussian'):
    """Calculates the mosaic distribution W(delta)."""
    if dist_type == 'gaussian':
        return (1 / (eta * np.sqrt(2 * np.pi))) * np.exp(-(delta ** 2) / (2 * eta ** 2))
    elif dist_type == 'lorentzian':
        return (1 / np.pi) * (eta / (delta ** 2 + eta ** 2))


def reflectivity_calc(delta, method='reflexion', dist_type='gaussian'):
    """Implementation of Equation 17 and 18 from Bacon & Lowde (1948)."""
    W = get_W(delta, eta_rad, dist_type)
    a = (Q / mu) * W # [cite: 331]
    A = mu * t0 / gamma0 # [cite: 330]

    if method == 'reflexion':
        # Precise integrand from Equation 17 [cite: 325]
        sqrt_term = np.sqrt(1 + 2 * a)
        # Using 1/tanh for coth
        denominator = (1 + a) + sqrt_term * (1.0 / np.tanh(A * sqrt_term))
        return a / denominator
    elif method == 'transmission':
        # Precise integrand from Equation 18 [cite: 329]
        return np.sinh(A * a) * np.exp(-A * (1 + a))


# X-axis range in degrees
delta_deg = np.linspace(-0.5, 0.5, 1000)
delta_rad = np.deg2rad(delta_deg)

# Calculation
R_refl_gauss = reflectivity_calc(delta_rad, 'reflexion', 'gaussian')
R_trans_gauss = reflectivity_calc(delta_rad, 'transmission', 'gaussian')
R_refl_lorentz = reflectivity_calc(delta_rad, 'reflexion', 'lorentzian')

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(delta_deg, R_refl_gauss, label='Reflexion (Gaussian)', color='blue', lw=2)
plt.plot(delta_deg, R_trans_gauss, label='Transmission (Gaussian)', color='red', linestyle='--')
plt.plot(delta_deg, R_refl_lorentz, label='Reflexion (Lorentzian)', color='green', alpha=0.7)

plt.title(f'Graphite 002 Reflectivity Profile\n($Q={Q:.3f}, \mu={mu:.3f}, \eta={eta_deg}^\circ$)')
plt.xlabel('$\Delta\Theta$ (Degrees)')
plt.ylabel('Reflectivity')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
#
#
# def get_W(delta, eta, dist_type='gaussian'):
#     """Calculates the mosaic distribution function W(delta)."""
#     if dist_type == 'gaussian':
#         # Eq. 144: Standard deviation eta
#         return (1 / (eta * np.sqrt(2 * np.pi))) * np.exp(-(delta ** 2) / (2 * eta ** 2))
#     elif dist_type == 'lorentzian':
#         # Lorentzian/Cauchy: eta as HWHM
#         return (1 / np.pi) * (eta / (delta ** 2 + eta ** 2))
#
#
# def calculate_reflectivity(delta, Q, mu, t0, gamma0, eta, method='reflexion', dist_type='gaussian'):
#     """
#     Calculates reflectivity based on Equations 17 (Reflexion) and 18 (Transmission).
#     """
#     W = get_W(delta, eta, dist_type)
#     a = (Q / mu) * W  # Eq. 331
#     A = mu * t0 / gamma0  # Eq. 330
#
#     if method == 'reflexion':
#         # Eq. 17 Integrand: Symmetrical Reflexion
#         sqrt_term = np.sqrt(1 + 2 * a)
#         denominator = (1 + a) + sqrt_term * (1.0 / np.tanh(A * sqrt_term))
#         return a / denominator
#
#     elif method == 'transmission':
#         # Eq. 18 Integrand: Symmetrical Transmission
#         return np.sinh(A * a) * np.exp(-A * (1 + a))
#
#
# # Setup Parameters
# eta_rad = np.radians(5.0 / 60.0)  # 5 arcminutes [cite: 246]
# Q = 0.02  # cm^-1 [cite: 246]
# mu = 0.1  # cm^-1 (Typical for neutrons) [cite: 311]
# t0 = 0.5  # 0.5 cm thickness
# gamma0 = 0.4  # cos(theta) [cite: 246]
#
# # eta_rad = np.deg2rad(0.1)  #
# # Q = 0.16693157449312368 #        cm^-1
# # mu = 9515.6342193855689 #     cm^-1
# # t0 = 0.1   # cm
# # thetaB = np.deg2rad(13.38257103)
# # gamma0 = np.sin(thetaB)
#
# delta_range = np.linspace(-10 * eta_rad, 10 * eta_rad, 1000)
# delta_arcmin = np.degrees(delta_range) * 60
#
# # Plotting
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
#
# # Plot 1: Reflexion vs Transmission (Gaussian)
# ax1.plot(delta_arcmin, calculate_reflectivity(delta_range, Q, mu, t0, gamma0, eta_rad, 'reflexion'),
#          label='Reflexion (Eq. 17)', color='blue')
# ax1.plot(delta_arcmin, calculate_reflectivity(delta_range, Q, mu, t0, gamma0, eta_rad, 'transmission'),
#          label='Transmission (Eq. 18)', color='red', linestyle='--')
# ax1.set_title('Reflexion vs Transmission Methods (Gaussian)')
# ax1.set_xlabel('$\Delta$ (arcminutes)')
# ax1.set_ylabel('Reflectivity')
# ax1.legend()
# ax1.grid(True, alpha=0.3)
#
# # Plot 2: Gaussian vs Lorentzian Distribution (Reflexion)
# ax2.plot(delta_arcmin, calculate_reflectivity(delta_range, Q, mu, t0, gamma0, eta_rad, 'reflexion', 'gaussian'),
#          label='Gaussian Mosaic', color='blue')
# ax2.plot(delta_arcmin, calculate_reflectivity(delta_range, Q, mu, t0, gamma0, eta_rad, 'reflexion', 'lorentzian'),
#          label='Lorentzian Mosaic', color='green')
# ax2.set_title('Effect of Mosaic Distribution Shape (Reflexion)')
# ax2.set_xlabel('$\Delta$ (arcminutes)')
# ax2.set_ylabel('Reflectivity')
# ax2.legend()
# ax2.grid(True, alpha=0.3)
#
# plt.tight_layout()
# plt.show()