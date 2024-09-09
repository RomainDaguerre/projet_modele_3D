import os
import math

class Atom:
    """This is the class Atom."""

    def __init__(self, name, residue_name, x=0.0, y=0.0, z=0.0):
        """Initialisation atom object.

        Coor x, y, z (floats) are 0.0 by default.
        """
        self.name = name
        self.residue_name = residue_name
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """Display atom object information."""
        chaine = f"atom {self.name}, res {self.residue_name}, " \
                 f"coor({self.x:7.3f}, {self.y:7.3f}, {self.z:7.3f})\n"
        return chaine

    def calc_dist(self, atom):
        """Distance between two atoms.

        First atom is self (the object itself).
        Second atom is passed as an argument to the method.
        """
        distance2 = (self.x - atom.x)**2 + (self.y - atom.y)**2 + \
                    (self.z - atom.z)**2
        return math.sqrt(distance2)

class Molecule:
    """This is the classe Molecule"""

    def __init__(self, name=None):
        """Initialisation molecule object.

        self.list_atoms will be filled with elements from the Atom class.
        self.list_distance will be filled by all the distance between every two atoms.
        """
        self.name = name
        self.list_atoms = []
        self.list_distance = []

    def __str__(self):
        """Display molecule object information."""
        message = f"Molecule {self.name}\n"
        if  self.list_atoms:
            for atom in self.list_atoms:
                message += atom.__str__()
        else:
            message += "No atom for the moment\n"
        if self.list_distance:
            for at1, at2, dist in self.list_distance:
                message += f"{at1} connected to {at2} : distance : {dist:.3f}\n"
        else:
            message += "No connectivity for the moment\n"
        return message
        
    def build_mlc_from_pdb(self, pdbfilename):
        """Recover information from a pdb file.

        The number, the name, the residue name and the coordinates for each atoms.
        """
        with open(pdbfilename, "r") as pdbfile:
            for line in pdbfile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atome_num = int(line[6:11])
                    atome_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    self.list_atoms.append(Atom(atome_name, residue_name, x, y, z))

    def calc_distance_list(self):
        """Distance between two atoms for every atoms in the list.

        Every atoms from the self.list_atoms.
        Return the list with every pair of atoms with the name, the residue name for 
        both atoms and the distance between the two atoms.
        """
        if not self.list_atoms:
            return None
        nb_atoms = len(self.list_atoms)
        for i in range(nb_atoms-1):
            for j in range(i+1, nb_atoms):
                distance = self.list_atoms[i].calc_dist(self.list_atoms[j])
                if (distance <= 15.5): # We don't need pair of atoms with distance over 15.5 angstrom (from 1027461 to 532414)
                    bond = (self.list_atoms[i].name, self.list_atoms[i].residue_name, 
                    self.list_atoms[j].name, self.list_atoms[j].residue_name , distance)
                    self.list_distance.append(bond)
        return self.list_distance