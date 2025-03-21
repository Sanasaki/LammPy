import re
from io import StringIO

from Python.ChemPy.AtomicSystems import Atom, MolecularSystem, Molecule
from Python.FileTypes.File import File


class XSDFile(File):
    @staticmethod
    def get_property_value(line: str, propertyName: str) -> str:
        """
        Retourne la valeur de la propriété propertyName dans la ligne donnée.
        Prévu pour être utilisé avec les fichiers XSD.
        Exemple: get_property_value('<Atom3d ID="1" Name="H" XYZ="0.0,0.0,0.0" />', 'XYZ') retourne '0.0,0.0,0.0'
        """
        match = re.search(rf'{propertyName}="([^"]+)"', line)
        if match:
            return match.group(1)
        return " "

    @staticmethod
    def find_molecules(chemicalBonds: list[tuple[Atom, Atom]]) -> list[Molecule]:
        """
        Regroupe les atomes en molécules basées sur les liaisons.
        bonds : liste de tuples (atom1, atom2) représentant les liaisons.
        Retourne une liste de molécules (listes d'atomes).
        """

        linkedAtoms: dict[Atom, list[Atom]] = {}
        # for atom1, atom2 in chemicalBonds:
        #     atom1.bond_to(atom2)

        # Recherche des molécules en parcourant le graphe
        visited_atoms: set[Atom] = set()

        # Fonction récursive pour explorer une molécule
        def appendAtomChain(atom: Atom, molecule: Molecule):
            """
            Add one atom (and recursively its neighbors) to the molecule.
            """
            if atom in visited_atoms:
                return
            visited_atoms.add(atom)
            molecule.addAtom(atom)
            for neighbor in linkedAtoms[atom]:
                appendAtomChain(neighbor, molecule)

        molecules: list[Molecule] = []

        for atom in linkedAtoms:
            if atom not in visited_atoms:
                molecule: Molecule = Molecule(atomList=[])
                appendAtomChain(atom, molecule)
                molecules.append(molecule)

        return molecules

    def get_bonds(self) -> list[tuple[Atom, Atom]]:
        self.atoms: dict[int, Atom] = self.get_atoms()
        bonds: list[tuple[Atom, Atom]] = []
        with open(self.filePath, "r") as xsdFile:
            # atomDict: dict[str, Atom] = {}
            for line in xsdFile:
                if ("<Bond" and "Connects") in line:
                    ## Hotfix to exclude HBond
                    if "HBond" in line:
                        continue
                    xsdIdAtom1, xsdIdAtom2 = self.get_property_value(line, "Connects").split(",")
                    # atomPairs += [(Atom(xsdIdAtom1), Atom(xsdIdAtom2))]
                    print(xsdIdAtom1, xsdIdAtom2)
                    # print(self.atoms[])
                    bonds.append((self.atoms[str(xsdIdAtom1)], self.atoms[str(xsdIdAtom2)]))
        return bonds

    def get_cell_parameters(self) -> list[list[str]]:
        with open(self.filePath, "r") as xsdFile:
            aVector: list[str] = []
            bVector: list[str] = []
            cVector: list[str] = []
            for line in xsdFile:
                if "SpaceGroup" in line:
                    aVector: list[str] = self.get_property_value(line, "AVector").split(",")
                    bVector: list[str] = self.get_property_value(line, "BVector").split(",")
                    cVector: list[str] = self.get_property_value(line, "CVector").split(",")
            return [aVector, bVector, cVector]

    def get_atoms(self) -> dict[int, Atom]:
        with open(self.filePath, "r") as xsdFile:
            self.atoms: dict[int, Atom] = {}
            for line in xsdFile:
                if ("<Atom3d" and "UserID") in line:
                    newAtom = Atom()
                    newAtom.id = self.get_property_value(line, "ID")
                    newAtom.label = self.get_property_value(line, "Name")
                    newAtom.x, newAtom.y, newAtom.z = self.get_property_value(line, "XYZ").split(",")
                    self.atoms[newAtom.id] = newAtom
                    print(newAtom, newAtom.id, newAtom.label, newAtom.x, newAtom.y, newAtom.z)
        return self.atoms

    def getCrystal(self) -> str:
        textBuffer = StringIO()
        atoms = self.get_atoms()
        bonds = self.get_bonds()
        aVector, bVector, cVector = self.get_cell_parameters()

        molecularSystem = MolecularSystem()
        self.molecules = molecularSystem.find_molecules(bonds)

        # print(atomPairs)
        # self.molecules = self.find_molecules(atomPairs)
        crystalMatrix: list[float] = [float(aVector[0]), float(bVector[1]), float(cVector[2])]

        textBuffer.write("\n")
        lammpsIdCorrespondingTo: dict[int, int] = {}
        i = 1
        for xsdId, atom in atoms.items():
            createAtom: str = f"create_atoms {atom.label} single {float(atom.x) * crystalMatrix[0]} {float(atom.y) * crystalMatrix[1]} {float(atom.z) * crystalMatrix[2]} remap yes #{xsdId}\n"
            lammpsIdCorrespondingTo[xsdId] = i
            i += 1
            textBuffer.write(createAtom)

        print(self.molecules)
        for i, molecule in enumerate(self.molecules):
            textBuffer.write("\n")
            for atom in molecule:
                atomBindsToMolecule: str = f"set atom {lammpsIdCorrespondingTo.get(atom.id)} mol {i + 1}\n"
                textBuffer.write(atomBindsToMolecule)

        # for molecule, atoms in correctedMolDict.items():
        #     for atomWithXsdID in atoms:

        # s'asssurer que ces groups soient à jour, ou au moins dynamiques.
        atomGroups: str = """
    group NitricHydrogenAtoms type 1
    group NitricNitrogenAtoms type 2
    group NitricOxygen1Atoms type 3
    group NitricOxygen2Atoms type 4
    group WaterOxygenAtoms type 5
    group WaterHydrogenAtoms type 6
    group NitrateNitrogenAtoms type 7
    group NitrateOxygenAtoms type 8
    group HydroniumHydrogenAtoms type 9
    group HydroniumOxygenAtoms type 10
    """
        textBuffer.write(atomGroups)

        setCharges: str = """
    set group NitricHydrogenAtoms charge 0.497
    set group NitricNitrogenAtoms charge 0.964
    set group NitricOxygen1Atoms charge -0.445
    set group NitricOxygen2Atoms charge -0.571
    set group NitrateNitrogenAtoms charge 0.65
    set group NitrateOxygenAtoms charge -0.55
    set group HydroniumHydrogenAtoms charge 0.578
    set group HydroniumOxygenAtoms charge -0.734
    """
        textBuffer.write(setCharges)

        createBonds: str = """
    create_bonds many NitricOxygen2Atoms NitricHydrogenAtoms 1 0.95 1.0
    create_bonds many NitricNitrogenAtoms NitricOxygen1Atoms 2 1.19 1.21
    create_bonds many NitricNitrogenAtoms NitricOxygen2Atoms 2 1.19 1.21
    create_bonds many WaterOxygenAtoms WaterHydrogenAtoms 3 0.98 1.1
    create_bonds many NitrateNitrogenAtoms NitrateOxygenAtoms 4 1.25 1.3
    create_bonds many HydroniumOxygenAtoms HydroniumHydrogenAtoms 5 0.98 1.1
    """
        textBuffer.write(createBonds)

        changeBox: str = f"\nchange_box all x final 0.0 {crystalMatrix[0]} y final 0.0 {crystalMatrix[1]} z final 0.0 {crystalMatrix[2]}\n"
        textBuffer.write(changeBox)

        return textBuffer.getvalue()


if __name__ == "__main__":
    xsdExample: str = r"C:\Archives\Thesis\Lab\VisualStudioCode\Python\LammPy\NAMCrystalEditedForLAMMPS.xsd"
    xsdFile = XSDFile(xsdExample)
    print(xsdFile.getCrystal())
