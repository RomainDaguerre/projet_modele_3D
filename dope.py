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
        

    def __str__(self):
        """Display dope object information."""
        return f"{self.residue1} {self.atom1} - {self.residue2} {self.atom2} : {self.dope_values}"


class DopeFile:
    """This is the DopeFile classs"""    

    def __init__(self):
        """Initialisation dopefile object.

        List of Dope class.
        """
        self.interactions = []

    def parse_file(self, filename):
        """Parses the DOPE file and stores the interactions."""
        with open(filename, 'r') as file:
            for line in file:
                elements = line.split()
                residue1 = elements[0]
                atom1 = elements[1]
                residue2 = elements[2]
                atom2 = elements[3]
                dope_values = [float(value) for value in elements[4:]]
                interaction = Dope(residue1, atom1, residue2, atom2, dope_values)
                self.interactions.append(interaction)

    def get_interaction(self, residue1, atom1, residue2, atom2):
        """Return the interaction between 2 atoms if it exists in the DOFE file"""
        for interaction in self.interactions:
            if (interaction.residue1 == residue1 and interaction.atom1 == atom1 and
                interaction.residue2 == residue2 and interaction.atom2 == atom2):
                return interaction
        return None

    def __str__(self):
        """Display dopefile object information."""
        message = ""
        for interaction in self.interactions:
            message += str(interaction) + "\n"
        return message