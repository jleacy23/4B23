import numpy as np
from src.noise import *
from src.params import *
from src.routing import generate_network_params

# Generate network parameters from routing optimization
NETWORKS = []
for config in NETWORK_CONFIGS:
    network = generate_network_params(
        name=config["name"],
        edges=config["edges"],
        demands=DEMANDS,
        max_channels=N_CHANNELS,
        max_paths=config["max_paths"],
        min_split_fraction=config["min_split_fraction"],
    )
    NETWORKS.append(network)

# ===========================================================================
# Phase 1: Pre-compute per-link properties for each fiber type
#   For each (network, fiber): opt_ptx, opt_psd, rx_power (no final amp),
#   link NSR (no final amp), number of in-link amplifiers (spans)
# ===========================================================================

# link_info[(network_name, fiber_type)][link_idx] = dict of link properties
link_info = {}

for network in NETWORKS:
    for fiber in FIBERS:
        attenuation = fiber["attenuation"]
        effective_area = fiber["effective_area"]
        dispersion = fiber["dispersion"]
        gamma = get_gamma(N_2, WAVELENGTH, effective_area)

        info_list = []
        for link_idx, (distance, channels) in enumerate(network["links"]):
            # Iteratively find true gain and optimum launch power
            gain_prev = UNSATURATED_GAIN
            while True:
                span_length = get_span_length(attenuation, gain_prev)
                span_eff_length = get_effective_length(span_length, attenuation)
                span_Cnli = get_C_nli(
                    gamma, attenuation, dispersion, span_eff_length, CHANNEL_BANDWIDTH
                )
                amplifier_noise = get_amplifier_noise(N_SP, HV, gain_prev)
                opt_ptx, opt_psd = get_optimal_launch_power(
                    amplifier_noise, span_Cnli, CHANNEL_BANDWIDTH
                )
                opt_ptx_linear = 10 ** (opt_ptx / 10)
                total_ptx_linear = channels * opt_ptx_linear
                gain = get_true_gain(UNSATURATED_GAIN, total_ptx_linear, P_SAT_LINEAR)
                if abs(gain - gain_prev) < 0.01:
                    break
                gain_prev = gain

            spans, excess = get_num_spans(distance, span_length)

            # NSR from integer spans (no extra amplifier)
            span_nsr = integer_span_nsr(
                UNSATURATED_GAIN, opt_ptx, spans, N_SP, SYMBOL_RATE, HV
            )

            # NSR from excess fiber (no extra amplifier)
            excess_eff_length = get_effective_length(excess, attenuation)
            excess_Cnli = get_C_nli(
                gamma, attenuation, dispersion, excess_eff_length, CHANNEL_BANDWIDTH
            )
            excess_nsr = get_excess_nsr(opt_psd, excess_Cnli)

            # Total link NSR in linear (no inter-link or final amplifiers)
            link_nsr_linear = 10 ** (span_nsr / 10) + 10 ** (excess_nsr / 10)

            # Received power at end of link with no extra amplifier
            rx_power = opt_ptx - attenuation * excess

            info_list.append(
                {
                    "distance": distance,
                    "channels": channels,
                    "opt_ptx": opt_ptx,  # dBm
                    "opt_ptx_linear": opt_ptx_linear,  # mW
                    "opt_psd": opt_psd,
                    "total_ptx_linear": total_ptx_linear,  # mW
                    "rx_power": rx_power,  # dBm (no final amp)
                    "link_nsr_linear": link_nsr_linear,
                    "spans": int(spans),
                    "gain": gain,  # true in-link amplifier gain dB
                }
            )

        link_info[(network["name"], fiber["type"])] = info_list

# ===========================================================================
# Phase 2: For each network, find the best fiber type
#   For each fiber: iterate routes, compute inter-link amplifier NSRs,
#   find the best valid transceiver per route, sum route capacities.
# ===========================================================================

final_results = {}

for network in NETWORKS:
    best_fiber_capacity = 0
    best_fiber_result = None

    all_fiber_results = []
    for fiber in FIBERS:
        info = link_info[(network["name"], fiber["type"])]
        fiber_capacity = 0
        route_details = []  # (route_idx, transceiver_type, capacity, total_nsr_dB)
        fiber_valid = True

        for route_idx, (link_indices, route_channels) in enumerate(network["routes"]):
            best_route_capacity = 0
            best_route_detail = None

            for transceiver in TRANSCEIVERS:
                max_nsr = transceiver["max_nsr"]
                min_power = transceiver["min_power"]
                transceiver_nsr = transceiver["transceiver_nsr"]
                data_rate = transceiver["data_rate"]

                # Accumulate NSR in linear domain
                total_nsr_linear = 10 ** (transceiver_nsr / 10)
                valid = True
                amplifiers_on_route = 0

                for pos, li in enumerate(link_indices):
                    link = info[li]

                    # Add link NSR (spans + excess, no inter-link amp)
                    total_nsr_linear += link["link_nsr_linear"]
                    amplifiers_on_route += link["spans"]

                    # Inter-link amplifier: boost rx power of this link
                    # to the opt_ptx of the next link (not needed for last link)
                    if pos < len(link_indices) - 1:
                        next_link = info[link_indices[pos + 1]]
                        # Gain needed: next link's opt_ptx - this link's rx_power
                        needed_gain = next_link["opt_ptx"] - link["rx_power"]  # dB
                        if needed_gain > 0:
                            # Compute amplifier NSR
                            rx_linear = 10 ** (link["rx_power"] / 10)  # mW
                            amp_noise_psd = get_amplifier_noise(N_SP, HV, needed_gain)
                            amp_nsr = get_amp_nsr(
                                amp_noise_psd,
                                CHANNEL_BANDWIDTH,
                                rx_linear,
                                needed_gain,
                            )
                            total_nsr_linear += amp_nsr
                            amplifiers_on_route += 1

                # Check final link received power
                final_link = info[link_indices[-1]]
                rx_power_final = final_link["rx_power"]

                if rx_power_final < min_power:
                    # Add a final amplifier with minimum gain to reach min_power
                    needed_gain = min_power - rx_power_final  # dB
                    rx_linear = 10 ** (rx_power_final / 10)  # mW
                    amp_noise_psd = get_amplifier_noise(N_SP, HV, needed_gain)
                    amp_nsr = get_amp_nsr(
                        amp_noise_psd, CHANNEL_BANDWIDTH, rx_linear, needed_gain
                    )
                    total_nsr_linear += amp_nsr
                    amplifiers_on_route += 1

                total_nsr_dB = 10 * np.log10(total_nsr_linear)

                if total_nsr_dB > max_nsr:
                    valid = False

                if valid:
                    route_capacity = route_channels * data_rate
                    if route_capacity > best_route_capacity:
                        best_route_capacity = route_capacity
                        best_route_detail = (
                            route_idx,
                            transceiver["type"],
                            route_capacity,
                            total_nsr_dB,
                            amplifiers_on_route,
                        )

            if best_route_detail is None:
                fiber_valid = False
                break

            fiber_capacity += best_route_capacity
            route_details.append(best_route_detail)

        if fiber_valid:
            if fiber_capacity > best_fiber_capacity:
                best_fiber_capacity = fiber_capacity
            all_fiber_results.append(
                {
                    "fiber": fiber["type"],
                    "capacity": fiber_capacity,
                    "routes": route_details,
                }
            )

    final_results[network["name"]] = {
        "best_capacity": best_fiber_capacity,
        "fibers": all_fiber_results,
    }

# ===========================================================================
# Print results
# ===========================================================================

for network_name, result in final_results.items():
    if not result["fibers"]:
        print(f"\nNo valid configuration found for {network_name}")
        continue
    print(f"\nResults for {network_name}:")
    for fiber_result in result["fibers"]:
        best_marker = (
            " <-- best" if fiber_result["capacity"] == result["best_capacity"] else ""
        )
        print(
            f"  Fiber type: {fiber_result['fiber']}  |  Total capacity: {fiber_result['capacity']:.1f} Tb/s{best_marker}"
        )
        for route_idx, tx_type, cap, nsr, amps in fiber_result["routes"]:
            print(
                f"      Route {route_idx}: Transceiver {tx_type}, "
                f"{cap:.1f} Tb/s, NSR={nsr:.2f} dB, {amps} amplifiers"
            )
