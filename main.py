import numpy as np
from src.noise import *
from src.params import *
from src.routing import generate_network_params

# Single-letter abbreviations for each city (printed as a key in results)
CITY_KEYS = {
    "London": "A",
    "Birmingham": "B",
    "Manchester": "C",
    "Leeds": "D",
    "Glasgow": "E",
}


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
            route_capacities = []
            node_capacity = {}

            for route_idx, (link_indices, route_channels, src, dst, *_) in enumerate(
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

                route_cap = best_data_rate * route_channels
                total_capacity += route_cap
                route_capacities.append(route_cap)
                node_capacity[dst] = node_capacity.get(dst, 0) + route_cap
                route_nsrs.append(best_route_nsr_db)

            total_amps = sum(link_amp_counts) + transceiver_amps
            fiber_results.append(
                {
                    "fiber": fiber["type"],
                    "total_capacity": total_capacity,
                    "total_amps": total_amps,
                    "route_nsrs": route_nsrs,
                    "route_capacities": route_capacities,
                    "node_capacity": node_capacity,
                }
            )

        all_results.append(
            {
                "name": network["name"],
                "fibers": fiber_results,
                "channels": [rt[1] for rt in network["routes"]],
                "route_nodes": [rt[4] for rt in network["routes"]],
            }
        )

    return all_results


def get_resilience_results(networks, baseline_results, edge_a, edge_b):
    """
    Compute resilience results without re-routing.
    Routes passing through the removed edge (edge_a <-> edge_b) have capacity
    set to zero; all other routes keep their baseline capacity.
    Returns result dicts with node_capacity but no route_nsrs / total_amps.
    """
    all_results = []
    for network, base_net_result in zip(networks, baseline_results):
        if network is None or base_net_result is None:
            all_results.append(None)
            continue

        # Find link indices for the removed undirected edge (both directions)
        removed_links = {
            li
            for li, (a, b) in enumerate(network["link_nodes"])
            if (a == edge_a and b == edge_b) or (a == edge_b and b == edge_a)
        }

        fiber_results = []
        for base_fr in base_net_result["fibers"]:
            total_capacity = 0
            node_capacity = {}

            for route_idx, (link_indices, route_channels, src, dst, *_) in enumerate(
                network["routes"]
            ):
                if removed_links and any(li in removed_links for li in link_indices):
                    cap = 0.0
                else:
                    cap = base_fr["route_capacities"][route_idx]

                total_capacity += cap
                node_capacity[dst] = node_capacity.get(dst, 0) + cap

            fiber_results.append(
                {
                    "fiber": base_fr["fiber"],
                    "total_capacity": total_capacity,
                    "node_capacity": node_capacity,
                }
            )

        all_results.append(
            {
                "name": network["name"],
                "fibers": fiber_results,
                "channels": base_net_result["channels"],
                "route_nodes": base_net_result["route_nodes"],
            }
        )

    return all_results


def print_results(results, label="", show_nsr=True, file=None):
    """Write results in the standard format to file (default: stdout)."""
    header = f"Results{' — ' + label if label else ''}"
    print(f"\n{'='*60}", file=file)
    print(f"  {header}", file=file)
    print(f"{'='*60}", file=file)
    for net_result in results:
        if net_result is None:
            print("  (network infeasible — disconnected graph)", file=file)
            continue
        print(f"\n  Network: {net_result['name']}", file=file)
        # Print city abbreviation key
        key_str = ", ".join(
            f"{v} = {k}" for k, v in sorted(CITY_KEYS.items(), key=lambda x: x[1])
        )
        print(f"  Key: {key_str}", file=file)
        route_nodes = net_result.get("route_nodes", [])
        if show_nsr:
            for i, ch in enumerate(net_result["channels"]):
                if i < len(route_nodes):
                    path_str = (
                        " ("
                        + "\u2192".join(CITY_KEYS.get(n, n) for n in route_nodes[i])
                        + ")"
                    )
                else:
                    path_str = ""
                print(f"    Route {i}{path_str}: Channels = {ch}", file=file)
        for fr in net_result["fibers"]:
            amps_str = (
                f"  |  Amplifiers: {fr['total_amps']}" if "total_amps" in fr else ""
            )
            print(
                f"    Fiber: {fr['fiber']}  |  "
                f"Total capacity: {fr['total_capacity']:.2f} Tb/s" + amps_str,
                file=file,
            )
            for node, cap in sorted(fr["node_capacity"].items()):
                node_k = CITY_KEYS.get(node, node)
                print(f"      Node {node} ({node_k}): {cap:.2f} Tb/s", file=file)
            if show_nsr and "route_nsrs" in fr:
                for i, nsr in enumerate(fr["route_nsrs"]):
                    nsr_str = f"{nsr:.2f} dB" if not np.isnan(nsr) else "invalid"
                    if i < len(route_nodes):
                        path_str = (
                            " ("
                            + "\u2192".join(CITY_KEYS.get(n, n) for n in route_nodes[i])
                            + ")"
                        )
                    else:
                        path_str = ""
                    print(f"      Route {i}{path_str}: NSR = {nsr_str}", file=file)


def main():
    output_path = "report/results.txt"
    with open(output_path, "w", encoding="utf-8") as f:
        # ------------------------------------------------------------------
        # Baseline results
        # ------------------------------------------------------------------
        networks = build_networks(NETWORK_CONFIGS)
        baseline = get_results(networks)
        print_results(baseline, label="Baseline", file=f)

        # ------------------------------------------------------------------
        # Resilience: repeat with each undirected edge removed in turn
        # For each network, track the edge deletion that gives the worst
        # (lowest) total capacity across all fiber types.
        # ------------------------------------------------------------------
        print(f"\n\n{'#'*60}", file=f)
        print("  Resilience Analysis — worst case per edge deletion", file=f)
        print(f"{'#'*60}", file=f)

        for ci, config in enumerate(NETWORK_CONFIGS):
            worst_capacity = float("inf")
            worst_label = ""
            worst_result = None
            worst_ei = None

            for ei, edge in enumerate(config["edges"]):
                edge_a, edge_b = edge[0], edge[1]
                edge_label = f"{edge_a}\u2013{edge_b}"
                result = get_resilience_results(networks, baseline, edge_a, edge_b)

                net_result = result[ci]
                if net_result is None:
                    cap = 0.0
                else:
                    cap = max(fr["total_capacity"] for fr in net_result["fibers"])

                if cap < worst_capacity:
                    worst_capacity = cap
                    worst_label = f"{config['name']} without {edge_label}"
                    worst_result = result
                    worst_ei = ei

            if worst_result is not None:
                print_results(
                    worst_result,
                    label=f"Worst case for {config['name']} — same channel assignment (edge removed: {worst_label})",
                    show_nsr=False,
                    file=f,
                )

                # Re-run routing and capacity calculation with the worst edge removed
                rerouted_networks = build_networks(
                    NETWORK_CONFIGS, excluded_edge=(ci, worst_ei)
                )
                rerouted_results = get_results(rerouted_networks)
                print_results(
                    rerouted_results,
                    label=f"Worst case for {config['name']} — re-routed (edge removed: {worst_label})",
                    show_nsr=False,
                    file=f,
                )

    print(f"Results written to {output_path}")


if __name__ == "__main__":
    main()
