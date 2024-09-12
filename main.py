import dope
import pdb_analysis
import os


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

def get_list_distance():
    """Return the list of ditance between every atoms from the Molecule 
    class create from the pdb file.
    """
    file_path = input("Saisir le chemin complet du fichier PDB: ")
    mol = pdb_analysis.Molecule(os.path.splitext(os.path.basename(file_path))[0])
    mol.build_mlc_from_pdb(file_path)
    return mol.calc_distance_list()

if __name__ == '__main__':
    """This is the main function"""

    dope_file = dope.DopeFile()
    dope_file.parse_file("donnees/dope.par")
    
    pdb_distance = get_list_distance()

    sum_pseudo_energy = 0
    for atom1, residue1, atom2, residue2, distance in pdb_distance:
        interaction = dope_file.get_interaction(residue1, atom1, residue2, atom2)
        if interaction:
            dope_value = get_dope_value_for_distance(interaction, distance)
            if dope_value is not None:
                sum_pseudo_energy += dope_value

    print(f"somme des pseudo_energies : {sum_pseudo_energy:.3f}")