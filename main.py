import dope
import pdb_analysis
import os
import argparse
import math
import numpy as np
import random

# Constante de Boltzmann en Joules par Kelvin
BOLTZMANN_CONSTANT = 1.380649e-23

def get_index_distance(distance):
    """Return the index of the list corresponding to the 
    value of the distance.
    """
    possible_distances = [0.25 + 0.5 * i for i in range(30)]
    closest_distance = min(possible_distances, key=lambda x: abs(x - distance))
    return possible_distances.index(closest_distance)

def get_dope_value_for_distance(interaction, distance):
    """Return the DOPE value according to the distance."""
    index_distance = get_index_distance(distance)
    return interaction[index_distance]

def get_list_distance(pdb_file):
    """Return the list of distances between every atom from the Molecule 
    class created from the pdb file.
    """
    mol = pdb_analysis.Molecule(os.path.splitext(os.path.basename(pdb_file))[0])
    mol.build_mlc_from_pdb(pdb_file)
    return mol.calc_distance_list()

def calculate_free_energy_boltzmann(pseudo_energy_sum, temperature):
    """Calculate the free energy using Boltzmann's constant and temperature."""
    pseudo_energy_joules = pseudo_energy_sum * 4184  # 1 kcal/mol = 4184 J/mol
    free_energy = pseudo_energy_joules / temperature
    return free_energy

def shuffle_residues(atom_pairs):
    """Shuffle the residues in the atom pairs while keeping atom positions the same."""
    residues = [(a1, r1, a2, r2) for a1, r1, a2, r2, dist in atom_pairs]
    random.shuffle(residues)
    shuffled_pairs = [(a1, r1, a2, r2, dist) for (a1, r1, a2, r2), (_, _, _, _, dist) in zip(residues, atom_pairs)]
    return shuffled_pairs

def calculate_pseudo_energy(atom_pairs, dope_file):
    """Calculate the sum of pseudo-energy for a given list of atom pairs."""
    sum_pseudo_energy = 0
    for atom1, residue1, atom2, residue2, distance in atom_pairs:
        interaction = dope_file.get_interaction(residue1, atom1, residue2, atom2)
        if interaction:
            dope_value = get_dope_value_for_distance(interaction, distance)
            if dope_value is not None:
                sum_pseudo_energy += dope_value
    return sum_pseudo_energy

def calculate_z_score(real_energy, random_energies):
    """Calculate the z-score given the real energy and a list of random energies."""
    mean_random = np.mean(random_energies)
    std_random = np.std(random_energies)
    z_score = (real_energy - mean_random) / std_random
    return z_score

def process_pdb_file(pdb_file, dope_file, temperature, num_random):
    """Process a single PDB file to calculate pseudo-energy, z-score, and free energy."""
    pdb_distance = get_list_distance(pdb_file)
    real_pseudo_energy = calculate_pseudo_energy(pdb_distance, dope_file)
    print(f"Fichier: {pdb_file}")
    print(f"somme des pseudo_energies (réel): {real_pseudo_energy:.3f}")

    # Generate random shuffles and calculate pseudo-energies
    random_energies = []
    for _ in range(num_random):
        shuffled_pairs = shuffle_residues(pdb_distance)
        random_energy = calculate_pseudo_energy(shuffled_pairs, dope_file)
        random_energies.append(random_energy)

    # Calculate the z-score
    z_score = calculate_z_score(real_pseudo_energy, random_energies)
    print(f"z-score : {z_score:.3f}")

    # Optionally calculate free energy with Boltzmann's constant
    free_energy_boltzmann = calculate_free_energy_boltzmann(real_pseudo_energy, temperature)
    print(f"énergie libre (Boltzmann approximation) : {free_energy_boltzmann:.3e} J")
    print("-" * 40)
    
    return z_score

if __name__ == '__main__':
    """Main function"""

    parser = argparse.ArgumentParser(description="Calculate DOPE pseudo-energy, free energy, and z-score for multiple PDB files in a folder.")
    parser.add_argument("pdb_folder", help="Path to the folder containing PDB files")
    parser.add_argument("-T", "--temperature", type=float, default=298.15, help="Temperature in Kelvin (default: 298.15 K)")
    parser.add_argument("-n", "--num_random", type=int, default=100, help="Number of random shuffles for z-score (default: 100)")
    args = parser.parse_args()

    dope_file = dope.DopeFile()
    dope_file.parse_file("donnees/dope.par")
    
    # List to store all z-scores
    z_scores = []

    # Process all PDB files in the given folder
    for pdb_file in os.listdir(args.pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(args.pdb_folder, pdb_file)
            z_score = process_pdb_file(pdb_path, dope_file, args.temperature, args.num_random)
            z_scores.append(z_score)

    # Calculate and display the mean and standard deviation of the z-scores
    if z_scores:
        mean_z_score = np.mean(z_scores)
        stddev_z_score = np.std(z_scores)
        print(f"Z-scores - Moyenne: {mean_z_score:.3f}, Écart-type: {stddev_z_score:.3f}")
    else:
        print("Aucun fichier PDB trouvé ou aucune donnée valide.")