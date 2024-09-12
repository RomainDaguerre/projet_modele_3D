from concurrent.futures import ThreadPoolExecutor

class Dope:
    """This is the Dope class.
    
    Pseudoenergy for a pair of atoms depending on the distance.
    """
    def __init__(self, residue1, atom1, residue2, atom2, dope_values):
        """Initialisation dope object.

        residue name and atom element for pair of atoms.
        """
        self.residue1 = residue1
        self.atom1 = atom1
        self.residue2 = residue2
        self.atom2 = atom2
        self.dope_values = dope_values


class DopeFile:
    """This is the DopeFile classs"""    

    def __init__(self):
        """Initialisation dopefile object.

        Dictionnary of interactions between 2 atoms.
        """
        self.interactions = {}

    def parse_file(self, filename):
        """Parses the DOPE file and stores the interactions in a dictionary."""
        with open(filename, 'r') as file:
            for line in file:
                elements = line.split()
                residue1 = elements[0]
                atom1 = elements[1]
                residue2 = elements[2]
                atom2 = elements[3]
                dope_values = [float(value) for value in elements[4:]]
                key = (residue1, atom1, residue2, atom2)
                self.interactions[key] = dope_values

    def get_interaction(self, residue1, atom1, residue2, atom2):
        """Return the interaction between 2 atoms if it exists in the DOPE file."""
        key = (residue1, atom1, residue2, atom2)
        return self.interactions.get(key, None)