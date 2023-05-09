"""
Author: Tobias Materzok https://github.com/TobiasMaterzok/
Date: 2022

Description:
This script analyzes the number of close van der Waals contacts between two groups of atoms within a 
specified distance threshold for molecular dynamics simulations. It takes the atom selections, distance 
threshold (Angstrom), topology file, and trajectory file as command-line arguments and iterates through the time steps 
of the trajectory, printing the number of close contacts between the specified atom groups at each time step.

Usage:
python script_name.py atom_group1 atom_group2 distance_threshold topology_file trajectory_file

Example:
python script_name.py "protein" "resname SOL" 15 run_pulling.tpr run_pulling.xtc
"""

import numpy as np
import sys
import MDAnalysis
from MDAnalysis.analysis import distances

def analyze_vdw_contacts(atom_group1, atom_group2, distance_threshold, topology_file, trajectory_file):
    """
    Analyze the number of close van der Waals contacts between two groups of atoms.

    This function calculates the number of close van der Waals contacts between
    two groups of atoms within a specified distance threshold (Angstrom) for molecular dynamics
    simulations. It iterates through the time steps in the trajectory, printing the
    time and number of close contacts between the specified atom groups at each time step.

    Parameters
    ----------
    atom_group1 : str
        Selection string for the first group of atoms.
    atom_group2 : str
        Selection string for the second group of atoms.
    distance_threshold : float
        Distance threshold for determining close contacts.
    topology_file : str
        Path to the topology file for the molecular dynamics simulation.
    trajectory_file : str
        Path to the trajectory file for the molecular dynamics simulation.

    Returns
    -------
    None
        Prints the time and number of close contacts between the two atom groups
        at each time step in the trajectory.
    """

    # Create a Universe object using the provided topology and trajectory files
    univ = MDAnalysis.Universe(topology_file, trajectory_file)
    
    # Select the specified atom groups from the universe
    group1 = univ.select_atoms(atom_group1)
    group2 = univ.select_atoms(atom_group2)

    results = []

    for ts in univ.trajectory:
        count_group = 0
        
        # Calculate the number of atoms in group1 within the specified distance from group2
        num_group = len(np.where(distances.distance_array(group1.positions, group2.positions, box=univ.dimensions, backend='OpenMP') < distance_threshold)[0])
        
        # Update count if the calculated value is greater than 0
        if num_group > 0:
            count_group = num_group
        
        # Append the time and count for each time step to the results list
        results.append((univ.trajectory.time, count_group))

    return results


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script_name.py atom_group1 atom_group2 distance_threshold topology_file trajectory_file")
        sys.exit(1)

    atom_group1 = sys.argv[1]
    atom_group2 = sys.argv[2]
    distance_threshold = float(sys.argv[3])
    topology_file = sys.argv[4]
    trajectory_file = sys.argv[5]

    results = analyze_vdw_contacts(atom_group1, atom_group2, distance_threshold, topology_file, trajectory_file)
    
    for time, count_group in results:
        print("%8f " % time, count_group)
