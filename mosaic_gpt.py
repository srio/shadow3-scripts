from srxraylib.plot.gol import set_qt
set_qt()

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Mosaic distribution W(Δ)
# -----------------------------
def W_mosaic(Delta, eta, kind="gaussian"):
    """
    Normalized mosaic block angular distribution W(Δ), so ∫ W dΔ = 1.

    Parameters
    ----------
    Delta : array [rad]
    eta   : spread parameter [rad]
        - gaussian: eta is the standard deviation (σ)
        - lorentzian: eta is the half-width at half-maximum (HWHM), i.e. gamma
    kind  : "gaussian" or "lorentzian"
    """
    Delta = np.asarray(Delta, dtype=float)
    kind = kind.lower()

    if eta <= 0:
        raise ValueError("eta must be > 0.")

    if kind == "gaussian":
        return (1.0 / (eta * np.sqrt(2.0 * np.pi))) * np.exp(-Delta**2 / (2.0 * eta**2))

    if kind == "lorentzian":
        # Cauchy/Lorentz distribution with HWHM = eta
        return (1.0 / np.pi) * (eta / (Delta**2 + eta**2))

    raise ValueError("kind must be 'gaussian' or 'lorentzian'.")

def coth_safe(x):
    """
    Numerically stable coth(x) for arrays.
    For small |x|, uses coth(x) ~ 1/x + x/3.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)

    ax = np.abs(x)
    small = ax < 1e-6

    out[small] = (1.0 / x[small]) + (x[small] / 3.0)
    out[~small] = np.cosh(x[~small]) / np.sinh(x[~small])
    return out

# -----------------------------
# Eq. (17) integrand (reflection method)
# r_refl(Δ) = a / [ (1+a) + sqrt(1+2a) * coth(A * sqrt(1+2a)) ]
# A = μ t0 / γ0,   a = (Q/μ) W(Δ)
# -----------------------------
def r_reflection_eq17(Delta, Q, mu, t0, gamma0, eta, dist="gaussian"):
    if mu <= 0:
        raise ValueError("Eq. (17) uses a=(Q/mu)W, so mu must be > 0.")

    A = (mu * t0) / gamma0
    a = (Q / mu) * W_mosaic(Delta, eta, kind=dist)

    s = np.sqrt(1.0 + 2.0 * a)
    denom = (1.0 + a) + s * coth_safe(A * s)
    return a / denom

# -----------------------------
# Eq. (18) integrand (transmission method)
# r_trans(Δ) = sinh(A a) * exp(-A(1+a))
# A = μ t0 / γ0,   a = (Q/μ) W(Δ)
# -----------------------------
def r_transmission_eq18(Delta, Q, mu, t0, gamma0, eta, dist="gaussian"):
    if mu <= 0:
        raise ValueError("Eq. (18) uses a=(Q/mu)W, so mu must be > 0.")

    A = (mu * t0) / gamma0
    a = (Q / mu) * W_mosaic(Delta, eta, kind=dist)
    return np.sinh(A * a) * np.exp(-A * (1.0 + a))

# -----------------------------
# Demo plot (edit parameters)
# -----------------------------
thetaB = np.deg2rad(30.0)
eta = np.deg2rad(5.0 / 60.0)  # 5 arcmin in radians
Q = 0.02   # 1/cm
mu = 0.3   # 1/cm
t0 = 0.1   # cm

thetaB = np.deg2rad(13.38257103)
eta = np.deg2rad(0.1)  # mosaicity, degrees
Q = 0.16693157449312368 #        cm^-1
mu = 9.6590899717931258 #     cm^-1
t0 = 0.1   # cm

# Symmetric geometries: paper uses γ0=sinθB for reflection, γ0=cosθB for transmission in captions.
gamma0_refl = np.sin(thetaB)
gamma0_trans = np.cos(thetaB)

Delta = np.linspace(-6*eta, 6*eta, 4001)

plt.figure(figsize=(7.4, 4.8))

for dist in ["gaussian", "lorentzian"]:
    r_refl = r_reflection_eq17(Delta, Q, mu, t0, gamma0_refl, eta, dist=dist)
    r_tran = r_transmission_eq18(Delta, Q, mu, t0, gamma0_trans, eta, dist=dist)

    plt.plot(Delta/eta, r_refl, lw=2, label=f"Reflection (Eq.17), {dist}")
    plt.plot(Delta/eta, r_tran, lw=2, ls="--", label=f"Transmission (Eq.18), {dist}")

plt.xlabel(r"Angular deviation $\Delta/\eta$")
plt.ylabel("Reflectivity integrand")
plt.title("Mosaic crystal reflectivity: Eq. (17) vs Eq. (18)\nGaussian vs Lorentzian mosaic spread")
plt.grid(True)
plt.legend(fontsize=9)
plt.tight_layout()
plt.show()
