import numpy as np


def get_span_length(attenuation, gain):
    """
    span length to allow EDFA to compensate for the loss

    attenuation: dB/km
    gain: dB
    """
    return gain / attenuation  # km


def get_num_spans(link_length, span_length):
    """
    link_length: km
    span_length: km
    """
    spans = link_length / span_length
    int_spans = np.floor(spans)
    excess = (spans - int_spans) * span_length
    return int_spans, excess  # -, km


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


def get_optimal_launch_power(amplifier_noise, C_nli, bandwidth):
    """
    amplifier_noise: psd of amplifier noise /pJ
    C_nli: pJ^-2
    bandwidth: GHz
    """
    psd_opt = (amplifier_noise / (2 * C_nli)) ** (1 / 3)  # pJ
    power_opt = psd_opt * bandwidth  # mW
    power_opt_dbm = 10 * np.log10(power_opt)  # dBm
    return power_opt_dbm, psd_opt


def integer_span_nsr(
    amplifier_gain, opt_launch_power, num_spans, n_sp, symbol_rate, hv
):
    """
    amplifier_gain: dB
    opt_launch_power: dBm
    num_spans: number of full spans in the link
    n_sp: population inversion factor
    symbol_rate: GBd
    hv: photon energy /J
    """
    noise_factor = 10 * np.log10(2 * n_sp)
    symbol_rate_bd = symbol_rate * 1e9
    nsr = (
        noise_factor
        + 10 * np.log10(num_spans)
        + amplifier_gain
        + 10 * np.log10(symbol_rate_bd * hv * 1e3)
        - opt_launch_power
    )  # dB
    nsr_linear = 10 ** (nsr / 10)
    return 10 * np.log10(nsr_linear * 1.5)  # account for N_nli = N_ase / 2


def get_excess_nsr(
    psd_tx,
    C_nli,
):
    """
    psd_opt: optimum launch psd /pJ
    C_nli: C_nli of the excess length /pJ^-2
    amplifier_noise: psd of amplifier noise /pJ
    amplifier_gain: dB
    add_amplifier: bool to determine if amplifier added at end of span
    excess_length: km
    attenuation: dB/km
    """
    excess_nsr = C_nli * psd_tx**2  # assume noise applied at input
    return 10 * np.log10(excess_nsr)


def get_amp_nsr(
    amp_noise_psd,
    bandwidth,
    pow_rx,
    gain,
):
    """
    amp_noise_psd : amplifier noise PSD /pJ
    bandwidth: GHz
    pow_rx: input power to amplifier /mW
    gain: amplifier gain /dB
    """
    noise = amp_noise_psd * bandwidth  # mW
    gain_linear = 10 ** (gain / 10)
    nsr = noise / (gain_linear * pow_rx)
    return nsr


def get_true_gain(unsaturated_gain, output_power, saturation_power):
    """
    unsaturated_gain: dB
    output_power: mW
    saturation_power: mW
    """
    return unsaturated_gain - 10 * np.log10(1 + output_power / saturation_power)


def get_rx_power(power_tx, attenuation, excess_length, extra_amplfier, amplifier_gain):
    """
    power_tx: transmitted power / dBm
    attenuation: dB/km
    excess_length: km
    extra_amplifier: bool
    amplifier_gain: dB
    """
    power_rx = power_tx - attenuation * excess_length
    if extra_amplfier:
        power_rx += amplifier_gain
    return power_rx
