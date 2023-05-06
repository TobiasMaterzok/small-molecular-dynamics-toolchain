"""
Author: Tobias Materzok https://github.com/TobiasMaterzok
Date: 2022

This script computes the work performed during a pull-off simulation.
The simulation consists of an adhesive material being pulled from a surface,
where the force linearly increases with time and then suddenly drops to below
zero as the adhesive detaches. The work is calculated by integrating the force
up to the point of maximum force (just before detachment).

The 'gd' data contains Gaussian-denoised force values, where the force is
represented as a function of the virtual cantilever position. The virtual
cantilever position is calculated by multiplying the pulling velocity by time.
The data has been denoised using a Gaussian kernel to remove thermal noise and improve
the quality of the data.

Usage: python py_compute_work.py <gd_data_file> <f_data_file> <x_data_file>
"""

import numpy as np
from scipy import integrate
import sys
from scipy.signal import find_peaks

# Load input data files from command line arguments
gaussian_denoised_data = np.loadtxt(sys.argv[1], comments=["@", "#"])
#@    xaxis  label "Position (nm)"
#@    yaxis  label "Force (kJ/mol/nm)"

force_data = np.loadtxt(sys.argv[2], comments=["@", "#"])
#@    xaxis  label "Time (ps)"
#@    yaxis  label "Force (kJ/mol/nm)"

position_data = np.loadtxt(sys.argv[3], comments=["@", "#"])
#@    xaxis  label "Time (ps)"
#@    yaxis  label "Position (nm)"

# Identify the index of the maximum force before detachment (prominence=30)
peak_indices = find_peaks(gaussian_denoised_data[:, 1], prominence=30)[0]
max_force_index = peak_indices[0]

# Extract force and position data up to the max_force_index
force_up_to_max = force_data[:np.argmax(force_data), 1]
gd_force_up_to_max = gaussian_denoised_data[:max_force_index, 1]
gd_virtual_positions = gaussian_denoised_data[:max_force_index, 0]
sim_positions = position_data[:max_force_index, 1]

# Calculate work performed using the trapezoidal rule for numerical integration
work_at_sim_positions = np.trapz(force_up_to_max, sim_positions)
work_at_gd_positions = np.trapz(gd_force_up_to_max, gd_virtual_positions)

# Output work performed
print(work_at_sim_positions)

# Write work performed to output files
with open("work_performed.xvg", "w") as output_file:
    print(work_at_sim_positions, file=output_file)

with open("work_performed_gd_positions.xvg", "w") as output_file:
    print(work_at_gd_positions, file=output_file)
