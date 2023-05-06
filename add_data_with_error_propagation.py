"""
Author: Tobias Materzok https://github.com/TobiasMaterzok/
Date: 2021

This script combines two data files with the same number of rows and a specific structure.
For example the output of https://github.com/TobiasMaterzok/Time-Position-Correlation-Calculator
Each file should have columns of numeric values, where the second and third columns represent
values and their uncertainties, respectively.

The first column is left unchanged, while the second columns of both files are added together,
and the third column represents the combined uncertainties of the added values.

The script uses the root-sum-square method (also known as sum in quadrature) to combine the
uncertainties. This method is commonly used in error propagation when the uncertainties are
considered to be independent and random. The formula for combining uncertainties in quadrature is:

    combined_uncertainty = sqrt((uncertainty1)^2 + (uncertainty2)^2)

Usage:
    python combine_data_files.py <input_file1> <input_file2> <output_file>
"""

import numpy as np
import sys

# Read command line arguments for input and output filenames
file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = sys.argv[3]

# Load data from the files, ignoring lines starting with '@' or '#'
data1 = np.loadtxt(file1, comments=['@', '#'])
data2 = np.loadtxt(file2, comments=['@', '#'])

# Add the second columns of both data files
data1[:, 1] += data2[:, 1]

# Combine uncertainties using the root-sum-square method
data1[:, 2] = np.sqrt(data1[:, 2]**2 + data2[:, 2]**2)

# Save the combined data to the output file
np.savetxt(output_file, data1, fmt='%16.8f')
