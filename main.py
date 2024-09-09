import dope
import pdb_analysis
import os

def get_closest_distance(distance):
    """Arrondit une distance à la plus proche valeur avec un pas de 0.5, allant de 0.25 à 15 Å."""
    possible_distances = [0.25 + 0.5 * i for i in range(30)]  # Génère les distances entre 0.25 et 15 Å
    closest_distance = min(possible_distances, key=lambda x: abs(x - distance))
    index = possible_distances.index(closest_distance)
    return index

def get_dope_value_for_distance(interaction, distance):
    """Retourne la valeur DOPE correspondant à la distance arrondie."""
    closest_distance = get_closest_distance(distance)
    return interaction.dope_values[closest_distance]

# Initialisation des fichiers PDB et DOPE
dope_file = dope.DopeFile()
dope_file.parse_file("donnees/dope.par")

# Saisie du fichier PDB par l'utilisateur
file_path = input("Saisir le chemin complet du fichier PDB: ")
file_name = os.path.basename(file_path) 
molecule_name = os.path.splitext(os.path.basename(file_name))[0]

# Création de la molécule à partir du fichier PDB
mol = pdb_analysis.Molecule(molecule_name)
mol.build_mlc_from_pdb(file_path)

# Calcul des distances entre les atomes de la molécule
list_distance = mol.calc_distance_list()

sum_pseudo_energy = 0

for atom1, residue1, atom2, residue2, distance in list_distance:
    # Recherche de l'interaction DOPE correspondante
    interaction = dope_file.get_interaction(residue1, atom1, residue2, atom2)  # Utilise les résidus et noms d'atomes comme critère
    if interaction:
        # Récupère la valeur DOPE correspondant à la distance arrondie
        dope_value = get_dope_value_for_distance(interaction, distance)
        if dope_value is not None:
            sum_pseudo_energy += dope_value

print(f"somme des pseudo_energies : {sum_pseudo_energy:.3f}")

#if __name__ == '__main__':
