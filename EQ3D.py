import dope
import pdb_analysis
import os
import argparse

def get_index_distance(distance):
    """Return the index of the list corresponding to the 
    value of the distance.
    """

    # The DOPE file contain pseudoenergy for distance between
    # 0.25 et 15 Å with a step of 0.5.
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

if __name__ == '__main__':
    """This is the main function"""

    parser = argparse.ArgumentParser(description="Calculate DOPE pseudo-energy from a PDB file.")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    args = parser.parse_args()

    dope_file = dope.DopeFile()
    dope_file.parse_file("donnees/dope.par")
    
    pdb_distance = get_list_distance(args.pdb_file)

    sum_pseudo_energy = 0
    for atom1, residue1, atom2, residue2, distance in pdb_distance:
        interaction = dope_file.get_interaction(residue1, atom1, residue2, atom2)
        if interaction:
            dope_value = get_dope_value_for_distance(interaction, distance)
            if dope_value is not None:
                sum_pseudo_energy += dope_value

    print(f"somme des pseudo_energies : {sum_pseudo_energy:.3f}")