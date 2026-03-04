import numpy as np

N_2 = 2.6 * 10 ** (-20)  # m^2/W
WAVELENGTH = 1550 * 10 ** (-9)  # m
SYMBOL_RATE = 200  # GBaud
CHANNEL_BANDWIDTH = 200  # GHz
N_CHANNELS = 25
CHANNEL_SPACING = 200  # GHz
SPEED_OF_LIGHT = 3 * 10**8  # m/s
HV = 1.3 * 1e-19

# Amplifier parameters
UNSATURATED_GAIN = 25  # dB
SATURATED_POWER = 23  # dBm
P_SAT = SATURATED_POWER - 10 * np.log10(0.69)
P_SAT_LINEAR = 10 ** (P_SAT / 10)
N_SP = 2


# Fiber parameters
FIBERS = [
    {
        "type": "NZ-DSF",
        "attenuation": 0.19,  # dB/km
        "dispersion": 3,  # ps/nm/km
        "effective_area": 72 * 1e-12,  # m^2
    },
    {
        "type": "PCSF",
        "attenuation": 0.15,  # dB/km
        "dispersion": 22,  # ps/nm/km
        "effective_area": 125 * 1e-12,
    },
    {
        "type": "NDSF",
        "attenuation": 0.18,  # dB/km
        "dispersion": 17,  # ps/nm/km
        "effective_area": 85 * 1e-12,
    },
]

# Transceiver parameters
TRANSCEIVERS = [
    {
        "type": "A",
        "data_rate": 1.6,  # Tb/s
        "max_nsr": -12.8,  # dB
        "min_power": -22,
        "transceiver_nsr": -20,
    },
    {
        "type": "B",
        "data_rate": 1.2,
        "max_nsr": -9.5,  # dB
        "min_power": -23,
        "transceiver_nsr": -20,
    },
    {
        "type": "C",
        "data_rate": 0.8,
        "max_nsr": -5.8,  # dB
        "min_power": -25,
        "transceiver_nsr": -20,
    },
    {
        "type": "D",
        "data_rate": 0.4,
        "max_nsr": -1,  # dB
        "min_power": -28,
        "transceiver_nsr": -20,
    },
]

# Network graph definitions
NETWORK_CONFIGS = [
    {
        "name": "Network 1",
        "edges": [
            ("London", "Birmingham", 162),
            ("London", "Leeds", 272),
            ("Birmingham", "Leeds", 113),
            ("Birmingham", "Manchester", 148),
            ("Manchester", "Leeds", 59),
            ("Manchester", "Glasgow", 295),
            ("Glasgow", "Leeds", 288),
        ],
        "max_paths": 6,
        "min_split_fraction": 0.05,
    },
    {
        "name": "Network 2",
        "edges": [
            ("London", "Birmingham", 162),
            ("London", "Leeds", 272),
            ("London", "Manchester", 262),
            ("Birmingham", "Manchester", 148),
            ("Manchester", "Glasgow", 295),
            ("Glasgow", "Leeds", 288),
        ],
        "max_paths": 6,
        "min_split_fraction": 0.05,
    },
]

# Traffic demands between city pairs
DEMANDS = {
    ("London", "Birmingham"): 0.12,
    ("London", "Manchester"): 0.09,
    ("London", "Glasgow"): 0.06,
    ("Birmingham", "London"): 0.12,
    ("Birmingham", "Manchester"): 0.04,
    ("Birmingham", "Glasgow"): 0.03,
    ("Birmingham", "Leeds"): 0.03,
    ("Manchester", "London"): 0.09,
    ("Manchester", "Birmingham"): 0.04,
    ("Manchester", "Glasgow"): 0.02,
    ("Manchester", "Leeds"): 0.02,
    ("Glasgow", "London"): 0.06,
    ("Glasgow", "Birmingham"): 0.03,
    ("Glasgow", "Manchester"): 0.02,
    ("Glasgow", "Leeds"): 0.02,
    ("Leeds", "London"): 0.07,
    ("Leeds", "Birmingham"): 0.03,
    ("Leeds", "Manchester"): 0.02,
    ("Leeds", "Glasgow"): 0.02,
}
