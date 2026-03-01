from src.noise import *
from src.params import *

for network in NETWORKS:

    print(f"Testing {network['name']}")

    for transceiver in TRANSCEIVERS:

        for fiber in FIBERS:

            attenuation = fiber["attenuation"]
            effective_area = fiber["effective_area"]
            dispersion = fiber["dispersion"]

            data_rate = transceiver["data_rate"]
            max_nsr = transceiver["max_nsr"]
            min_power = transceiver["min_power"]
            transceiver_nsr = transceiver["transceiver_nsr"]

            print(
                f"Using fiber: {fiber['type']}, transceiver type {transceiver['type']}"
            )

            for distance, channels in network["links"]:
                span_length = get_span_length(attenuation, UNSATURATED_GAIN)
                spans, excess = get_num_spans(distance, span_length)
                gamma = get_gamma(N_2, WAVELENGTH, effective_area)
                span_eff_length = get_effective_length(span_length, attenuation)
                span_Cnli = get_C_nli(
                    gamma, attenuation, dispersion, span_eff_length, CHANNEL_BANDWIDTH
                )
                amplifier_noise = get_amplifier_noise(N_SP, HV, UNSATURATED_GAIN)
                opt_ptx = get_optimal_launch_power(
                    amplifier_noise, span_Cnli, CHANNEL_BANDWIDTH
                )
                span_nsr = integer_span_nsr(
                    UNSATURATED_GAIN, opt_ptx, spans, N_SP, SYMBOL_RATE, HV
                )
