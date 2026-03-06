import numpy as np


def get_gamma(n2, wavelength, effective_area):
    """
    n2: m^2/W
    wavelength: m
    effective_area: m^2
    """
    k = 2 * np.pi / wavelength  # m^-1
    return (k * n2 / effective_area) * 1e3  # W^-1 km^-1


def get_effective_length(length, attenuation):
    """
    length: km
    attenuation: dB/km
    """
    attenuation_nepers = 0.23 * attenuation
    l_eff = (1 - np.exp(-attenuation_nepers * length)) / attenuation_nepers  # km
    return l_eff


def get_C_nli(gamma, attenuation, dispersion, effective_length, bandwidth):
    """
    gamma: W^-1 m^-1
    attenuation: dB/km
    dispersion: ps/nm/km
    effective_length: km
    bandwidth: GHz
    """
    attenuation_nepers = 0.23 * attenuation  # km^-1
    beta_2 = 1 / 0.78 * dispersion  # ps^2/km

    C_nli = (
        (8 * gamma**2 * effective_length**2 * attenuation_nepers)
        / (27 * np.pi * beta_2)
        * np.log(beta_2 * np.pi**2 * (bandwidth * 1e-3) ** 2 / attenuation_nepers)
    )  # pJ^-2

    return C_nli


def get_amplifier_noise(n_sp, hv, gain):
    """
    n_sp: population inversion factor
    hv: photon energy /J
    gain: amplifier gain /dB
    """
    gain_linear = 10 ** (gain / 10)
    return 2 * n_sp * hv * (gain_linear - 1) * 1e12  # pJ


def link_nsr_amps(amplifier_gain, psd_tx, num_spans, n_sp, symbol_rate, hv):
    """
    amplifier_gain: dB
    psd_tx: pJ
    num_spans: number of full spans in the link
    n_sp: population inversion factor
    symbol_rate: GBd
    hv: photon energy /J
    """
    noise_factor = 10 * np.log10(2 * n_sp)
    pow_tx = 10 * np.log10(psd_tx * symbol_rate)  # dBm
    symbol_rate_bd = symbol_rate * 1e9

    nsr = (
        noise_factor
        + 10 * np.log10(num_spans)
        + amplifier_gain
        + 10 * np.log10(symbol_rate_bd * hv * 1e3)
        - pow_tx
    )  # dB
    nsr_linear = 10 ** (nsr / 10)
    return nsr_linear


def link_nsr_nli(C_nli, psd_tx, num_spans):
    """
    C_nli: non-linear coefficient of a span pJ^-2
    psd_tx: pJ
    num_spans: number of spans on a link
    """
    nsr_span = C_nli * psd_tx**2
    nsr_total = num_spans * nsr_span
    return nsr_total  # linear


def get_amp_nsr(
    amp_noise_psd,
    psd_rx,
    gain,
):
    """
    amp_noise_psd : amplifier noise PSD /pJ
    psd_rx: input power to amplifier /mW
    gain: amplifier gain /dB
    """
    gain_linear = 10 ** (gain / 10)
    nsr = amp_noise_psd / (gain_linear * psd_rx)
    return nsr


def get_true_gain(unsaturated_gain, output_power, saturation_power):
    """
    unsaturated_gain: dB
    output_power: mW
    saturation_power: mW
    """
    return unsaturated_gain - 10 * np.log10(1 + output_power / saturation_power)


def get_pow_rx(attenuation, span_length, gain, num_spans, pow_tx):
    """
    attenuation: dB/km
    span_length: km
    gain: true amplfier gain /dB
    num_spans: number of spans on link
    pow_tx: transmitted power /dBm
    """

    gain_per_span = gain - attenuation * span_length  # dB
    gain = gain_per_span + 10 * np.log10(num_spans)
    pow_rx = pow_tx + gain
    return pow_rx


def get_opt_ptx(C_nli, N_ase, bandwidth):
    """
    C_nli: non-linear coefficient of span /pj^-2
    N_ase: amplifier noise PSD /pJ
    bandwidth: GHz
    """

    psd_opt = (N_ase / (2 * C_nli)) ** (1 / 3)  # pJ
    ptx_opt_linear = psd_opt * bandwidth  # mw
    ptx_opt = 10 * np.log10(ptx_opt_linear)

    return psd_opt, ptx_opt_linear, ptx_opt


def get_ptx_from_gain(gain_req, unsaturated_gain, p_sat, bandwidth, num_channels):
    """
    gain_req: required gain /dB
    unsaturated_gain: dB
    p_sat: saturation power /mW
    bandwidth: GHz
    num_channels: number of channels
    """
    if gain_req > unsaturated_gain:
        print("Warning: required gain is above the unsaturated gain")
        ratio = 2
    else:
        ratio = 10 ** ((unsaturated_gain - gain_req) / 10)
    ptx_total_linear = (ratio - 1) * p_sat  # mW
    ptx_linear = ptx_total_linear / num_channels
    psd_linear = ptx_linear / bandwidth
    return psd_linear, ptx_linear, 10 * np.log10(ptx_linear)
