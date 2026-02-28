N_2 = 2.6 * 10 ** (-20)  # m^2/W
WAVELENGTH = 1550 * 10 ** (-9)  # m
SYMBOL_RATE = 200 * 10**9  # Baud
CHANNEL_BANDWIDTH = 100 * 10**9  # Hz
N_CHANNELS = 25
CHANNEL_SPACING = 200 * 10**9  # Hz

# Amplifier parameters
UNSATURATED_GAIN = 25  # dB
SATURATED_POWER = 23  # dBm
N_SP = 2


# Fiber parameters
ATTENUATIONS = {"NZ-DSF": 0.19, "PCSF": 0.15, "NDSF": 0.18}  # dB/km
DISPERSIONS = {"NZ-DSF": 3, "PCSF": 22, "NDSF": 17}  # ps/nm/km
EFFECTIVE_AREAS = {
    "NZ-DSF": 72 * 10 ** (-12),
    "PCSF": 125 * 10 ** (-12),
    "NDSF": 85 * 10 ** (-12),
}  # m^2

# Transceiver parameters
NET_DATA_RATES = {
    "A": 1600 * 10**9,
    "B": 1200 * 10**9,
    "C": 800 * 10**9,
    "D": 400 * 10**9,
}  # bps
MAX_NSR = {"A": -12.8, "B": -9.5, "C": -5.8, "D": -1}  # dB
MIN_POWER = {"A": -22, "B": -23, "C": -25, "D": -28}  # dBm
TRANSCEIVER_NSR = {"A": -20, "B": -20, "C": -20, "D": -20}  # dB

# Network parameters
CITIES = {
    "London": 1,
    "Manchester": 2,
    "Birmingham": 3,
    "Leeds": 4,
    "Glasgow": 5,
}  # cities and identifying indices for distances, traffic and network graph.
DISTANCES = [
    [0, 262, 162, 272, 555],
    [262, 0, 113, 59, 295],
    [162, 113, 0, 148, 407],
    [272, 59, 148, 0, 288],
    [555, 295, 407, 288, 0],
]  # km
TRAFFIC = [
    [0, 9, 12, 7, 6],
    [9, 0, 4, 2, 2],
    [12, 4, 0, 3, 3],
    [7, 2, 3, 0, 2],
    [6, 2, 3, 2, 0],
]  # %
NETWORKS = [
    [[3, 4], [3, 4, 5], [1, 2, 4], [1, 2, 3, 5], [2, 4]],
    [[2, 3, 4], [1, 3, 4, 5], [1, 2], [1, 2, 5], [2, 4]],
]  # undirected graph of the network, where each sublist contains the indices of the neighboring cities.
