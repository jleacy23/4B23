import numpy as np
from src.noise import *
from src.params import *
from src.routing import generate_network_params


def build_networks(configs, excluded_edge=None):
    """
    Generate network parameter dicts from NETWORK_CONFIGS.
    If excluded_edge is (config_idx, edge_idx), that undirected edge is omitted.
    Returns None for a network if routing fails (e.g. disconnected graph).
    """
    networks = []
    for ci, config in enumerate(configs):
        edges = config["edges"]
        if excluded_edge is not None and excluded_edge[0] == ci:
            edges = [e for i, e in enumerate(edges) if i != excluded_edge[1]]
        try:
            network = generate_network_params(
                name=config["name"],
                edges=edges,
                demands=DEMANDS,
                max_channels=N_CHANNELS,
                max_paths=config["max_paths"],
                min_split_fraction=config["min_split_fraction"],
            )
            networks.append(network)
        except Exception:
            networks.append(None)
    return networks


def get_results(networks):
    """
    Run Phase 1 and Phase 2 for a list of networks.
    Returns a list of per-network result dicts (None if the network is None).
    Each result dict has keys: name, fibers.
    Each fiber entry has keys: fiber, total_capacity, total_amps, route_nsrs.
    """
    # ------------------------------------------------------------------
    # Phase 1: Compute NSR and transmission power (per channel) of each link
    # ------------------------------------------------------------------
    link_results = {}

    for network in networks:
        if network is None:
            continue

        for fiber in FIBERS:
            link_results[(network["name"], fiber["type"])] = []

            for distance, channels in network["links"]:
                attenuation = fiber["attenuation"]
                dispersion = fiber["dispersion"]
                effective_area = fiber["effective_area"]

                max_span_length = UNSATURATED_GAIN / attenuation  # km
                min_spans = np.ceil(distance / max_span_length)
                span_length = distance / min_spans

                gamma = get_gamma(N_2, WAVELENGTH, effective_area)
                L_eff = get_effective_length(span_length, attenuation)
                C_nli = get_C_nli(
                    gamma, attenuation, dispersion, L_eff, CHANNEL_BANDWIDTH
                )

                # Method 1: Sweep per-channel launch power, pick lowest total NSR
                # Attenuator after each amp ensures every span sees the same launch power,
                # so only the first amp gain (from saturation) needs to be computed.
                ptx_sweep_db = np.linspace(-10, 10, 500)  # dBm per channel
                best_nsr_1 = np.inf
                ptx_opt_db = ptx_sweep_db[0]
                for ptx_db in ptx_sweep_db:
                    ptx_linear = 10 ** (ptx_db / 10)  # mW per channel
                    psd = ptx_linear / CHANNEL_BANDWIDTH  # pJ (PSD per channel)
                    ptx_total = ptx_linear * channels  # mW total across all channels
                    gain = get_true_gain(
                        UNSATURATED_GAIN, ptx_total, P_SAT_LINEAR
                    )  # dB
                    nsr_amps = link_nsr_amps(
                        gain, psd, min_spans, N_SP, SYMBOL_RATE, HV
                    )
                    nsr_nli = link_nsr_nli(C_nli, psd, min_spans)
                    nsr_total = nsr_amps + nsr_nli
                    if nsr_total < best_nsr_1:
                        best_nsr_1 = nsr_total
                        ptx_opt_db = ptx_db
                nsr_1 = best_nsr_1

                # Method 2: Reduce gain to avoid attenuators
                gain_req = span_length * attenuation
                psd_tx, ptx_linear, ptx_db = get_ptx_from_gain(
                    gain_req,
                    UNSATURATED_GAIN,
                    P_SAT_LINEAR,
                    CHANNEL_BANDWIDTH,
                    channels,
                )
                nsr_amp_2 = link_nsr_amps(
                    gain_req, psd_tx, min_spans, N_SP, SYMBOL_RATE, HV
                )
                nsr_nli_2 = link_nsr_nli(C_nli, psd_tx, min_spans)
                nsr_2 = nsr_amp_2 + nsr_nli_2

                if nsr_1 < nsr_2:
                    nsr, ptx_final_db = nsr_1, ptx_opt_db
                else:
                    nsr, ptx_final_db = nsr_2, ptx_db

                link_results[(network["name"], fiber["type"])].append(
                    (nsr, ptx_final_db, min_spans)
                )

    # ------------------------------------------------------------------
    # Phase 2: Compute NSR of each route and assign transceivers
    # ------------------------------------------------------------------
    all_results = []

    for network in networks:
        if network is None:
            all_results.append(None)
            continue

        fiber_results = []

        for fiber in FIBERS:
            results = link_results[(network["name"], fiber["type"])]
            total_capacity = 0
            route_nsrs = []

            # Per-link amplifier counts (in-link span amps) + inter-link amp flags
            link_amp_counts = [
                int(results[li][2]) for li in range(len(network["links"]))
            ]
            link_amp_flags = [False] * len(network["links"])
            transceiver_amps = 0

            for route_idx, (link_indices, route_channels) in enumerate(
                network["routes"]
            ):
                nsr_linear = sum(results[li][0] for li in link_indices)

                for pos in range(len(link_indices) - 1):
                    ptx_cur = results[link_indices[pos]][1]
                    ptx_next = results[link_indices[pos + 1]][1]
                    if ptx_next > ptx_cur:
                        next_li = link_indices[pos + 1]
                        next_link_channels = network["links"][next_li][1]
                        ptx_next_total = 10 ** (ptx_next / 10) * next_link_channels
                        amp_gain = get_true_gain(
                            UNSATURATED_GAIN, ptx_next_total, P_SAT_LINEAR
                        )
                        amp_noise_psd = get_amplifier_noise(N_SP, HV, amp_gain)
                        psd_rx = 10 ** (ptx_cur / 10) / CHANNEL_BANDWIDTH
                        nsr_linear += get_amp_nsr(amp_noise_psd, psd_rx, amp_gain)
                        if not link_amp_flags[next_li]:
                            link_amp_counts[next_li] += 1
                            link_amp_flags[next_li] = True

                ptx_last = results[link_indices[-1]][1]
                best_data_rate = 0
                best_route_nsr_db = float("nan")
                best_needs_final_amp = False

                for transceiver in TRANSCEIVERS:
                    route_nsr = nsr_linear

                    needs_final_amp = ptx_last < transceiver["min_power"]
                    if needs_final_amp:
                        p_out_total = (
                            10 ** (transceiver["min_power"] / 10) * route_channels
                        )
                        amp_gain = get_true_gain(
                            UNSATURATED_GAIN, p_out_total, P_SAT_LINEAR
                        )
                        amp_noise_psd = get_amplifier_noise(N_SP, HV, amp_gain)
                        psd_rx = 10 ** (ptx_last / 10) / CHANNEL_BANDWIDTH
                        route_nsr += get_amp_nsr(amp_noise_psd, psd_rx, amp_gain)

                    route_nsr += 10 ** (transceiver["transceiver_nsr"] / 10)
                    route_nsr_db = 10 * np.log10(route_nsr)

                    if route_nsr_db <= transceiver["max_nsr"]:
                        if transceiver["data_rate"] > best_data_rate:
                            best_data_rate = transceiver["data_rate"]
                            best_route_nsr_db = route_nsr_db
                            best_needs_final_amp = needs_final_amp

                if best_data_rate > 0 and best_needs_final_amp:
                    transceiver_amps += 1

                total_capacity += best_data_rate * route_channels
                route_nsrs.append(best_route_nsr_db)

            total_amps = sum(link_amp_counts) + transceiver_amps
            fiber_results.append(
                {
                    "fiber": fiber["type"],
                    "total_capacity": total_capacity,
                    "total_amps": total_amps,
                    "route_nsrs": route_nsrs,
                }
            )

        all_results.append({"name": network["name"], "fibers": fiber_results})

    return all_results


def print_results(results, label=""):
    """Print results in the standard format."""
    header = f"Results{' — ' + label if label else ''}"
    print(f"\n{'='*60}")
    print(f"  {header}")
    print(f"{'='*60}")
    for net_result in results:
        if net_result is None:
            print("  (network infeasible — disconnected graph)")
            continue
        print(f"\n  Network: {net_result['name']}")
        for fr in net_result["fibers"]:
            print(
                f"    Fiber: {fr['fiber']}  |  "
                f"Total capacity: {fr['total_capacity']:.2f} Tb/s  |  "
                f"Amplifiers: {fr['total_amps']}"
            )
            for i, nsr in enumerate(fr["route_nsrs"]):
                nsr_str = f"{nsr:.2f} dB" if not np.isnan(nsr) else "invalid"
                print(f"      Route {i}: NSR = {nsr_str}")


def main():
    # ------------------------------------------------------------------
    # Baseline results
    # ------------------------------------------------------------------
    networks = build_networks(NETWORK_CONFIGS)
    baseline = get_results(networks)
    print_results(baseline, label="Baseline")

    # ------------------------------------------------------------------
    # Resilience: repeat with each undirected edge removed in turn
    # For each network, track the edge deletion that gives the worst
    # (lowest) total capacity across all fiber types.
    # ------------------------------------------------------------------
    print(f"\n\n{'#'*60}")
    print("  Resilience Analysis — worst case per edge deletion")
    print(f"{'#'*60}")

    for ci, config in enumerate(NETWORK_CONFIGS):
        worst_capacity = float("inf")
        worst_label = ""
        worst_result = None

        for ei, edge in enumerate(config["edges"]):
            edge_label = f"{edge[0]}–{edge[1]}"
            reduced_networks = build_networks(NETWORK_CONFIGS, excluded_edge=(ci, ei))
            result = get_results(reduced_networks)

            # Best total capacity for this network with this edge removed
            net_result = result[ci]
            if net_result is None:
                cap = 0.0
            else:
                cap = max(fr["total_capacity"] for fr in net_result["fibers"])

            if cap < worst_capacity:
                worst_capacity = cap
                worst_label = f"{config['name']} without {edge_label}"
                worst_result = result

        if worst_result is not None:
            print_results(
                worst_result,
                label=f"Worst case for {config['name']} (edge removed: {worst_label})",
            )


if __name__ == "__main__":
    main()
