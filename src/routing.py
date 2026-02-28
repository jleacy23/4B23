"""
Optical Network Traffic Routing Optimizer
==========================================
Solves the multi-commodity flow problem to find optimal traffic splitting
ratios across routes, minimising maximum link utilisation.

Dependencies:
    pip install networkx cvxpy numpy

Usage:
    See the example at the bottom of this file, or import and call
    optimize_traffic_routing() directly.
"""

import itertools
import networkx as nx
import cvxpy as cp
import numpy as np
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class RoutingResult:
    """Results from the optimisation."""

    status: str  # 'optimal', 'infeasible', 'unbounded', etc.
    max_link_utilisation: float  # The minimised MLU value (0–1+)
    routing: dict  # {(src, dst): [(path, fraction), ...]}
    link_utilisation: dict  # {(u, v): utilisation_fraction}

    def print_summary(self):
        print(f"\n{'='*60}")
        print(f"  Optimisation Status : {self.status}")
        if self.status != "optimal":
            print("  No feasible solution found.")
            return
        print(f"  Max Link Utilisation: {self.max_link_utilisation:.1%}")
        print(f"{'='*60}")

        print("\n--- Routing Splits ---")
        for (src, dst), splits in sorted(self.routing.items()):
            print(f"\n  {src} → {dst}:")
            for path, frac in splits:
                path_str = " → ".join(str(n) for n in path)
                print(f"    {frac:6.1%}  via  {path_str}")

        print("\n--- Link Utilisations ---")
        for (u, v), util in sorted(self.link_utilisation.items(), key=lambda x: -x[1]):
            bar = "█" * int(util * 20)
            print(f"  {u} → {v}:  {util:6.1%}  {bar}")
        print()


# ---------------------------------------------------------------------------
# Core optimiser
# ---------------------------------------------------------------------------


def optimize_traffic_routing(
    graph: nx.DiGraph,
    demands: dict,
    capacity_attr: str = "capacity",
    max_paths: int = 8,
    min_split_fraction: float = 0.0,
) -> RoutingResult:
    """
    Find optimal traffic splitting ratios to minimise maximum link utilisation.

    Parameters
    ----------
    graph : nx.DiGraph
        Directed graph. Each edge must have a numeric attribute for capacity
        (name given by capacity_attr). For undirected optical links, add both
        directions as separate edges.

    demands : dict
        Traffic demands as fractions of total capacity, e.g.:
            {('A', 'B'): 0.30, ('C', 'D'): 0.20}
        Values should be in [0, 1] — they represent the proportion of a
        reference unit of traffic that must be routed.

    capacity_attr : str
        Name of the edge attribute storing link capacity.

    max_paths : int
        Maximum number of simple paths to consider per demand pair.
        Paths are found in order of increasing hop count.

    min_split_fraction : float
        If > 0, any path that is used must carry at least this fraction
        of the demand (avoids tiny impractical splits). Introduces binary
        variables, turning the problem into a MIP — still fast at small scale.

    Returns
    -------
    RoutingResult
    """

    # ------------------------------------------------------------------
    # 1. Enumerate candidate paths for each demand pair
    # ------------------------------------------------------------------
    demand_paths = {}  # {(s, d): [list_of_node_tuples]}
    for src, dst in demands:
        if src not in graph or dst not in graph:
            raise ValueError(f"Node '{src}' or '{dst}' not found in graph.")
        paths = list(
            itertools.islice(nx.shortest_simple_paths(graph, src, dst), max_paths)
        )
        if not paths:
            raise ValueError(f"No path exists between {src} and {dst}.")
        demand_paths[(src, dst)] = [tuple(p) for p in paths]

    # ------------------------------------------------------------------
    # 2. Build edge → path incidence structure
    # ------------------------------------------------------------------
    # Flat list of (demand_key, path_index, path_tuple)
    all_paths = []
    for sd, paths in demand_paths.items():
        for i, p in enumerate(paths):
            all_paths.append((sd, i, p))

    num_vars = len(all_paths)

    # edges present in graph
    edges = list(graph.edges())
    edge_index = {e: i for i, e in enumerate(edges)}

    # incidence[edge_idx] = list of (var_idx, demand_key)
    incidence = {i: [] for i in range(len(edges))}
    for var_idx, (sd, _, path) in enumerate(all_paths):
        for u, v in zip(path[:-1], path[1:]):
            if (u, v) in edge_index:
                incidence[edge_index[(u, v)]].append((var_idx, sd))

    # ------------------------------------------------------------------
    # 3. Build and solve LP
    # ------------------------------------------------------------------
    x = cp.Variable(num_vars, nonneg=True)  # flow fractions
    u = cp.Variable()  # max utilisation (to minimise)

    constraints = []

    # (a) Flow conservation: fractions for each demand must sum to 1
    ptr = 0
    for sd, paths in demand_paths.items():
        n = len(paths)
        constraints.append(cp.sum(x[ptr : ptr + n]) == 1)
        ptr += n

    # (b) Capacity constraints: load on each edge ≤ u * capacity
    for edge_idx, (u_node, v_node) in enumerate(edges):
        cap = graph[u_node][v_node].get(capacity_attr)
        if cap is None:
            raise ValueError(
                f"Edge ({u_node}, {v_node}) is missing the '{capacity_attr}' attribute."
            )
        contributors = incidence[edge_idx]
        if not contributors:
            continue
        load = cp.sum([demands[sd] * x[var_idx] for var_idx, sd in contributors])
        constraints.append(load <= u * cap)

    # (c) Optional: minimum split fraction (turns LP into MIP)
    if min_split_fraction > 0:
        y = cp.Variable(num_vars, boolean=True)  # 1 if path is used
        constraints += [
            x >= min_split_fraction * y,
            x <= y,
        ]

    objective = cp.Minimize(u)
    problem = cp.Problem(objective, constraints)
    problem.solve(solver=cp.CLARABEL)

    # ------------------------------------------------------------------
    # 4. Package results
    # ------------------------------------------------------------------
    if problem.status not in ("optimal", "optimal_inaccurate"):
        return RoutingResult(
            status=problem.status,
            max_link_utilisation=float("inf"),
            routing={},
            link_utilisation={},
        )

    x_val = np.array(x.value)
    x_val = np.clip(x_val, 0, 1)  # numerical safety

    # Build routing splits, dropping negligible paths (< 0.1%)
    routing = {}
    ptr = 0
    for sd, paths in demand_paths.items():
        n = len(paths)
        fracs = x_val[ptr : ptr + n]
        splits = [
            (path, float(frac)) for path, frac in zip(paths, fracs) if frac > 1e-3
        ]
        # Renormalise to exactly 1.0
        total = sum(f for _, f in splits)
        splits = [(p, f / total) for p, f in splits]
        routing[sd] = splits
        ptr += n

    # Compute per-link utilisation
    link_util = {}
    for edge_idx, (u_node, v_node) in enumerate(edges):
        cap = graph[u_node][v_node][capacity_attr]
        contributors = incidence[edge_idx]
        if not contributors:
            link_util[(u_node, v_node)] = 0.0
            continue
        load = sum(demands[sd] * float(x_val[var_idx]) for var_idx, sd in contributors)
        link_util[(u_node, v_node)] = load / cap

    return RoutingResult(
        status=problem.status,
        max_link_utilisation=float(u.value),
        routing=routing,
        link_utilisation=link_util,
    )


# ---------------------------------------------------------------------------
# Channel assignment
# ---------------------------------------------------------------------------


def assign_channels(result: RoutingResult, demands: dict, max_channels: int) -> dict:
    """
    Assign wavelength channels to each transmission route.

    The link with the highest utilisation is fully loaded with max_channels.
    Routes passing through that link are assigned channels proportionally to
    their share of traffic on it.  All remaining routes are scaled from any
    one of those reference assignments using:

        channels = reference_channels * (route_demand / reference_demand)

    Parameters
    ----------
    result : RoutingResult
        Output of optimize_traffic_routing().
    demands : dict
        Original demand dict {(src, dst): volume} passed to the optimiser.
    max_channels : int
        Total channels available on every link (all links have the same limit).

    Returns
    -------
    dict  {(src, dst, path): channel_count}
        Channel count for every active route (path is a tuple of nodes).
    """

    # ------------------------------------------------------------------
    # 1. Compute per-route effective demand
    #    route_demand = total_demand(s,d) * split_fraction
    # ------------------------------------------------------------------
    route_demands = {}  # {(src, dst, path): effective_demand}
    for (src, dst), splits in result.routing.items():
        for path, fraction in splits:
            route_demands[(src, dst, path)] = demands[(src, dst)] * fraction

    # ------------------------------------------------------------------
    # 2. Find the most utilised link
    # ------------------------------------------------------------------
    busiest_link = max(result.link_utilisation, key=result.link_utilisation.get)
    print(
        f"  Most utilised link: {busiest_link[0]} → {busiest_link[1]} "
        f"({result.link_utilisation[busiest_link]:.1%})"
    )

    # ------------------------------------------------------------------
    # 3. Identify routes that pass through the busiest link
    # ------------------------------------------------------------------
    def uses_link(path, link):
        u, v = link
        return any(a == u and b == v for a, b in zip(path[:-1], path[1:]))

    # Total traffic flowing through the busiest link
    total_load_on_busiest = sum(
        rd
        for (s, d, path), rd in route_demands.items()
        if uses_link(path, busiest_link)
    )

    # Assign channels to routes through the busiest link
    reference_channels = {}  # will hold one entry to use as reference for others
    channel_assignments = {}

    for (s, d, path), rd in route_demands.items():
        if uses_link(path, busiest_link):
            channels = round(max_channels * (rd / total_load_on_busiest))
            channel_assignments[(s, d, path)] = channels
            reference_channels[(s, d, path)] = (channels, rd)  # store all as candidates

    # ------------------------------------------------------------------
    # 4. Assign channels to remaining routes using any reference route
    # ------------------------------------------------------------------
    # Pick the reference route with the largest demand for numerical stability
    ref_key = max(reference_channels, key=lambda k: reference_channels[k][1])
    ref_ch, ref_demand = reference_channels[ref_key]

    for (s, d, path), rd in route_demands.items():
        if (s, d, path) not in channel_assignments:
            channels = round(ref_ch * (rd / ref_demand))
            channel_assignments[(s, d, path)] = channels

    # ------------------------------------------------------------------
    # 5. Correct rounding overflows
    #    Rounding independently can push multiple links over max_channels.
    #    The loop processes ALL overloaded links: each iteration picks the
    #    most overloaded link and removes one channel from the route through
    #    it with the most channels.  Because the assignment is stored per
    #    route, that decrement automatically applies to every other link the
    #    same route passes through, potentially resolving multiple overloads
    #    in one step.  The loop recomputes all link totals after every
    #    decrement, so it always has a complete picture of remaining
    #    overloads.  Termination is guaranteed: each iteration reduces the
    #    total channel count by 1, and all counts are bounded below by 0.
    # ------------------------------------------------------------------
    def compute_link_totals(assignments):
        totals = {}
        for (s, d, path), ch in assignments.items():
            for u, v in zip(path[:-1], path[1:]):
                totals[(u, v)] = totals.get((u, v), 0) + ch
        return totals

    link_totals = compute_link_totals(channel_assignments)
    overloaded = {lnk: tot for lnk, tot in link_totals.items() if tot > max_channels}

    while overloaded:
        # Most overloaded link first
        worst_link = max(overloaded, key=overloaded.get)

        # Route through that link with the most channels
        candidate = max(
            (
                key
                for key in channel_assignments
                if uses_link(key[2], worst_link) and channel_assignments[key] > 0
            ),
            key=lambda k: channel_assignments[k],
        )
        channel_assignments[candidate] -= 1
        print(
            f"  [rounding fix] -{1} channel from "
            f"{candidate[0]} → {candidate[1]} "
            f"via {" → ".join(str(n) for n in candidate[2])} "
            f"(link {worst_link[0]} → {worst_link[1]} "
            f"was {overloaded[worst_link]}/{max_channels})"
        )

        # Recompute ALL link totals and ALL overloaded links
        link_totals = compute_link_totals(channel_assignments)
        overloaded = {
            lnk: tot for lnk, tot in link_totals.items() if tot > max_channels
        }

    return channel_assignments


def print_channel_summary(channel_assignments: dict, max_channels: int):
    print(f"\n--- Channel Assignments (max per link: {max_channels}) ---")
    for (src, dst, path), ch in sorted(
        channel_assignments.items(), key=lambda x: -x[1]
    ):
        path_str = " → ".join(str(n) for n in path)
        print(f"  {src} → {dst}  via  {path_str}:  {ch} channels")


def check_link_channel_usage(channel_assignments: dict, max_channels: int) -> dict:
    """
    Sum the channels used on every link across all routes that traverse it.

    The busiest link should sum to exactly max_channels (within floating point
    precision). All other links will be <= max_channels.

    Parameters
    ----------
    channel_assignments : dict
        Output of assign_channels(): {(src, dst, path): channel_count}.
    max_channels : int
        The per-link channel limit, used to compute headroom.

    Returns
    -------
    dict  {(u, v): total_channels_used}
    """
    link_channels = {}

    for (src, dst, path), ch in channel_assignments.items():
        for u, v in zip(path[:-1], path[1:]):
            link_channels[(u, v)] = link_channels.get((u, v), 0) + ch

    print(f"\n--- Channel Usage per Link (max: {max_channels}) ---")
    for (u, v), total in sorted(link_channels.items(), key=lambda x: -x[1]):
        headroom = max_channels - total
        bar = "█" * int((total / max_channels) * 20)
        if total > max_channels:
            status = " ⚠ OVER LIMIT"
        elif headroom == 0:
            status = " ◀ busiest"
        else:
            status = ""
        print(
            f"  {u} → {v}:  {total:3d} / {max_channels} channels  "
            f"({headroom:+d} spare)  {bar}{status}"
        )

    return link_channels


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # -------------------------------------------------------------------
    # Define your network here
    # -------------------------------------------------------------------

    # Directed graph: add an edge in each direction for bidirectional links.
    # Each edge needs a 'capacity' value (use any consistent unit, e.g. Gbps,
    # or just a normalised value between 0 and 1).
    G = nx.DiGraph()

    # Edges: (from, to, capacity)
    edges = [
        ("A", "B", 1),
        ("A", "E", 1),
        ("B", "A", 1),
        ("B", "C", 1),
        ("B", "E", 1),
        ("C", "B", 1),
        ("C", "D", 1),
        ("C", "E", 1),
        ("D", "C", 1),
        ("D", "E", 1),
        ("E", "A", 1),
        ("E", "B", 1),
        ("E", "C", 1),
        ("E", "D", 1),
    ]
    for u, v, cap in edges:
        G.add_edge(u, v, capacity=cap)

    # Demands: {(source, destination): traffic_volume}
    # Use the same units as capacity (e.g. Gbps), or normalised fractions.
    demands = {
        ("A", "B"): 0.12,
        ("A", "C"): 0.09,
        ("A", "D"): 0.06,
        ("B", "A"): 0.12,
        ("B", "C"): 0.04,
        ("B", "D"): 0.03,
        ("B", "E"): 0.03,
        ("C", "A"): 0.09,
        ("C", "B"): 0.04,
        ("C", "D"): 0.02,
        ("C", "E"): 0.02,
        ("D", "A"): 0.06,
        ("D", "B"): 0.03,
        ("D", "C"): 0.02,
        ("D", "E"): 0.02,
        ("E", "A"): 0.07,
        ("E", "B"): 0.03,
        ("E", "C"): 0.02,
        ("E", "D"): 0.02,
    }

    # -------------------------------------------------------------------
    # Run optimisation
    # -------------------------------------------------------------------
    result = optimize_traffic_routing(
        graph=G,
        demands=demands,
        capacity_attr="capacity",
        max_paths=6,
        min_split_fraction=0.0,  # Set e.g. 0.05 to avoid splits < 5%
    )

    result.print_summary()

    # -------------------------------------------------------------------
    # Assign channels
    # -------------------------------------------------------------------
    MAX_CHANNELS = 25

    channels = assign_channels(result, demands, max_channels=MAX_CHANNELS)
    print_channel_summary(channels, max_channels=MAX_CHANNELS)
    check_link_channel_usage(channels, max_channels=MAX_CHANNELS)

    CHANNEL_CAPACITY = 1.6  # Tbps per channel

    total_capacity = sum(channels.values()) * CHANNEL_CAPACITY
    print(f"\nTotal assigned capacity: {total_capacity:.1f} Tbps")
