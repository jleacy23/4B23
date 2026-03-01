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

# Network parameters
NETWORKS = [{"name": "Network 1"}, {"name": "Network 2"}]
