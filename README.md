4B23 – Optical Fibre Communication  
Coursework – Optical Network Design – marks 15/60 
Learning objectives: 
• To be able to calculate the throughput of an optical network 
• To understand impact of topology on network throughput 
• To understand the design decisions and trade-offs that occur in optical network design 
 
The coursework exercise is to design an optical network to link the cities of London, Manchester, 
Birmingham, Leeds and Glasgow. Students are required to write an individual report detailing their 
proposed design and expected performance. 
 
The physical distance (in km) between each node pairs is as follows (with the percentage normalized 
traffic flows in parentheses)

London Manchester Birmingham Leeds Glasgow 
London / - / 262 (9) / 162 (12) / 272 (7) / 555 (6) 
Manchester / 262 (9) / - / 113 (4) / 59 (2) / 295 (2) 
Birmingham / 162 (12) / 113 (4) / - / 148 (3) / 407 (3) 
Leeds  / 272 (7) / 59 (2) / 148 (3) / - / 288 (2) 
Glasgow  / 555 (6) / 295 (2) / 407 (3) / 288 (2) / - 

Design an optical network to connect these cities and calculate the maximum network throughput for 
the given traffic matrix. Optical routing is assumed to minimize the number of transceivers 
utilized (the insertion loss of any wavelength selective switch may be neglected). Within the report 
two different possible topologies should be compared both in terms of network throughput, resilience 
and total length of fibre deployed. For each topology you should include an NSR table giving the 
NSR (in dB) including the transceiver NSR. Any assumptions made should be stated. 
 
The report should be no more than 10 sides of A4 with minimum font size of 11. Detailed calculations 
regarding design choices such as amplifier spacing, launch power etc. may be included in a technical 
appendix that is not subject to page limits. Deadline for submission via moodle is 4pm on 25/03/26.  
 
Specifications for photonic subsystems 
Optical amplifier 
C-band optical amplifier (suitable for amplifying up to 25 channels spaced at 200 GHz), with 
unsaturated gain of 25 dB, saturated output power of 23 dBm, with nsp=2 (noise figure of 6 dB). Input 
and output variable attenuators included. 
  
Optical Fibres 
NZ-DSF (LEAF) / PSCF (EX2500) / NDSF (SMF-28) 
Attenuation (dB/km) / 0.19 / 0.15 / 0.18 
Dispersion ps/nm/km / 3 / 22 / 17 
Effective area µm2 / 72 / 125 / 85 
 
Adaptive transceiver modes (all operate at 200 GBd with 200 GHz channels spacing)  
Mode / A / Mode B / Mode C / Mode D 
Net data rate (Gb/s) / 1600 / 1200 / 800 / 400 
Maximum NSR (dB) / -12.8 / -9.5 / -5.8 / -1 
Minimum required optical power (dBm) / -22 / -23 / -25 / -28 
Transceiver NSR (dB) / -20.0 / -20.0 / -20.0 / -20.0

OPTIMIZATION GOAL:
Maximize network throughput subject to the constraints on NSR and optical power of the transceivers.

KEY INFORMATION:
- at optimal launch psd the non-linear noise is half that of amplifier noise
- NSR is additive in the linear domain (not decibels), will only need to consider non-linear NSR, amplifier NSR, Transceiver NSR.
 