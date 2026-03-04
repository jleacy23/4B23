import numpy as np
from src.noise import *
from src.params import *

results = {}

for network in NETWORKS:

    print(f"Testing {network['name']}")

    for transceiver in TRANSCEIVERS:
        transceiver_valid = True

        for fiber in FIBERS:
            valid = True

            total_channels = network["total_channels"]

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

            nsr_table = []
            amplifiers_used = 0
            total_length = 0

            for distance, channels in network["links"]:

                extra_amplifier = False
                while True:

                    print(f"distance: {distance}, channels: {channels}")

                    gamma = get_gamma(N_2, WAVELENGTH, effective_area)

                    gain_prev = UNSATURATED_GAIN
                    while True:
                        span_length = get_span_length(attenuation, gain_prev)
                        span_eff_length = get_effective_length(span_length, attenuation)
                        span_Cnli = get_C_nli(
                            gamma,
                            attenuation,
                            dispersion,
                            span_eff_length,
                            CHANNEL_BANDWIDTH,
                        )
                        amplifier_noise = get_amplifier_noise(N_SP, HV, gain_prev)
                        opt_ptx, opt_psd = get_optimal_launch_power(
                            amplifier_noise, span_Cnli, CHANNEL_BANDWIDTH
                        )
                        opt_ptx_linear = 10 ** (opt_ptx / 10)
                        total_ptx_linear = channels * opt_ptx_linear
                        gain = get_true_gain(
                            UNSATURATED_GAIN, total_ptx_linear, P_SAT_LINEAR
                        )
                        if abs(gain - gain_prev) < 0.01:
                            break
                        gain_prev = gain
                    print(f"    span_length [km]: {span_length}")
                    print(f"    true amplifier gain [dB]: {gain_prev}")
                    print(f"    opt_ptx [dBm]: {opt_ptx}")
                    spans, excess = get_num_spans(distance, span_length)
                    span_nsr = integer_span_nsr(
                        UNSATURATED_GAIN, opt_ptx, spans, N_SP, SYMBOL_RATE, HV
                    )
                    excess_eff_length = get_effective_length(excess, attenuation)
                    excess_Cnli = get_C_nli(
                        gamma,
                        attenuation,
                        dispersion,
                        excess_eff_length,
                        CHANNEL_BANDWIDTH,
                    )
                    excess_nsr = get_excess_nsr(
                        opt_psd,
                        excess_Cnli,
                        amplifier_noise,
                        gain_prev,
                        extra_amplifier,
                        excess,
                        attenuation,
                    )
                    total_nsr_linear = (
                        10 ** (span_nsr / 10)
                        + 10 ** (excess_nsr / 10)
                        + 10 ** (transceiver_nsr / 10)
                    )
                    total_nsr = 10 * np.log10(total_nsr_linear)
                    print(
                        f"     span NSR: {span_nsr}, excess NSR: {excess_nsr}, total NSR in link:  {total_nsr} \n"
                    )
                    rx_power = get_rx_power(
                        opt_ptx, attenuation, excess, extra_amplifier, gain_prev
                    )
                    if rx_power > min_power:
                        break
                    if rx_power < min_power and extra_amplifier:
                        print("    Warning: power still too low with extra amplifier")
                        break
                    extra_amplifier = True

                if total_nsr > max_nsr:
                    valid = False
                nsr_table.append(total_nsr)
                total_length += distance
                amplifiers_used += spans + (1 if extra_amplifier else 0)

            capacity = total_channels * data_rate
            results[(network["name"], fiber["type"], transceiver["type"])] = (
                valid,
                capacity,
                nsr_table,
                total_length,
                amplifiers_used,
            )

# choose result with highest capacity that is valid for each network, amplifiers used as tie breaker.
final_results = {}
for network in NETWORKS:
    best_capacity = 0
    best_amplifiers_used = float("inf")
    best_result = None
    for fiber in FIBERS:
        for transceiver in TRANSCEIVERS:
            result = results[(network["name"], fiber["type"], transceiver["type"])]
            if (
                result[0]
                and result[1] > best_capacity
                and result[4] < best_amplifiers_used
            ):
                best_capacity = result[1]
                best_amplifiers_used = result[4]
                best_result = (
                    fiber["type"],
                    transceiver["type"],
                    best_capacity,
                    result[2],
                    result[3],
                    result[4],
                )
    final_results[network["name"]] = best_result

for network_name, (
    fiber_type,
    transceiver_type,
    capacity,
    nsr_table,
    total_length,
    amplifiers_used,
) in final_results.items():
    print(f"Best result for {network_name}:")
    print(f"    Fiber type: {fiber_type}")
    print(f"    Transceiver type: {transceiver_type}")
    print(f"    Capacity: {capacity} Tbps")
    print(f"    NSR table: {nsr_table}")
    print(f"    Total length: {total_length} km")
    print(f"    Amplifiers used: {amplifiers_used}\n")
