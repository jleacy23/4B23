import numpy as np


def non_linear_coefficient(A_eff, n_2, wavelength):
    """Returns the non-linear coefficient gamma
    Parameters
    A_eff: effective area of the fiber [m^2]
    n_2: non_linear refractive index [m^2/W]
    wavelength: wavelength of the signal [m]
    """
    c = 3e8
    return (2 * np.pi * n_2) / (wavelength * c * A_eff)


def effective_length(alpha, length):
    """Returns the effective length of the fiber
    Parameters
    alpha: attenuation of the fiber [dB/km]
    length: length of the fiber [km]
    """
    alpha_nepers = 0.23 * alpha
    return (1 - np.exp(-alpha_nepers * length)) / alpha_nepers


def get_C_nli(alpha, dispersion, length, bandwidth, A_eff, n_2, wavelength):
    """Returns C_nli
    Parameters
    alpha: attenuation of the fiber [dB/km]
    dispersion: dispersion of the fiber [ps/nm/km]
    length: length of the fiber [km]
    bandwidth: bandwidth of the channel [Hz]
    A_eff: effective area of the fiber [m^2]
    n_2: non_linear refractive index [m^2/W]
    wavelength: wavelength of the signal [m]
    """
    alpha_nepers = 0.23 * alpha
    beta_2 = 1 / 0.78 * dispersion * 1e-27  # ps^2/km -> s^2/m
    L_eff = effective_length(alpha, length) * 1e3  # km -> m
    gamma = non_linear_coefficient(A_eff, n_2, wavelength)
    return (
        (8 * gamma**2 * L_eff**2 * alpha_nepers)
        / (27 * np.pi * beta_2)
        * np.log(beta_2 / alpha_nepers * np.pi**2 * bandwidth**2)
    )


def amplifier_noise(gain, n_sp):
    """Returns the power spectral density of an amplifier
    Parameters
    gain: gain of the amplifier [dB]
    n_sp: spontaneous emission factor
    """
    gain_linear = 10 ** (gain / 10)
    hv = 1.3 * 10 ** (-19)  # J
    return 2 * n_sp * hv * (gain_linear - 1)  # W/Hz


def optimal_launch_psd(
    gain, n_sp, alpha, dispersion, length, bandwidth, A_eff, n_2, wavelength
):
    """Returns the optimal launch power spectral density
    Parameters
    gain: gain of the amplifier [dB]
    n_sp: spontaneous emission factor
    alpha: attenuation of the fiber [dB/km]
    dispersion: dispersion of the fiber [ps/nm/km]
    length: length of the fiber [km]
    bandwidth: bandwidth of the channel [Hz]
    A_eff: effective area of the fiber [m^2]
    n_2: non_linear refractive index [m^2/W]
    wavelength: wavelength of the signal [m]
    """
    N_ase = amplifier_noise(gain, n_sp)
    C_nli = get_C_nli(alpha, dispersion, length, bandwidth, A_eff, n_2, wavelength)
    return (N_ase / (2 * C_nli)) ** (1 / 3)  # W/Hz


def get_span_length(alpha, gain):
    """Returns the span length for a given gain and attenuation
    Parameters
    alpha: attenuation of the fiber [dB/km]
    gain: gain of the amplifier [dB]
    """
    return gain / alpha  # km
