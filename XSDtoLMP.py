import re
from io import StringIO

xsdExample: str = r"C:\Users\JL252842\Documents\Thesis\Lab\MaterialsStudio\CreatingLammpsCrystals_Files\Documents\NAMCrystalEditedForLAMMPS.xsd"


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


def find_molecules(bonds: list[tuple[str, str]]) -> list[list[str]]:
    """
    Regroupe les atomes en molécules basées sur les liaisons.
    bonds : liste de tuples (atom1, atom2) représentant les liaisons.
    Retourne une liste de molécules (listes d'atomes).
    """
    # Construire un graphe sous forme de dictionnaire d'adjacence
    graph: dict[str, str] = {}
    for a1, a2 in bonds:
        graph[a1].add(a2)
        graph[a2].add(a1)

    # Recherche des molécules en parcourant le graphe
    visited_atoms: set[str] = set()

    # Fonction récursive pour explorer une molécule
    def merge_with_neighbors(atom: str, molecule: list[str]):
        """
        Modifie la liste renseignée en argument pour qu'elle inclue tous les atomes voisins (et leurs propres voisins) de l'atome donné.
        Input: un atome, molécule mère
        Output: molécule mère modifiée
        """
        if atom in visited_atoms:
            return
        visited_atoms.add(atom)
        molecule.append(atom)
        for neighbor in graph[atom]:
            merge_with_neighbors(neighbor, molecule)

    molecules: list[list[str]] = []

    for atom in graph:
        if atom not in visited_atoms:
            molecule: list[str] = []
            merge_with_neighbors(atom, molecule)
            molecules.append(molecule)

    return molecules


def getCrystal(xsdFilePath: str) -> str:
    textBuffer = StringIO()
    aVector: list[str] = []
    bVector: list[str] = []
    cVector: list[str] = []
    with open(xsdFilePath, "r") as xsdFile:
        atomDict: dict[str, tuple[str, str, str, str]] = {}
        lmpIDxsdIDdict: dict[str, str] = {}
        i = 1
        atomPairs: list[tuple[str, str]] = []
        for line in xsdFile:
            # Retriving atoms
            if ("<Atom3d" and "UserID") in line:
                atomID: str = get_property_value(line, "ID")
                lammpsType: str = get_property_value(line, "Name")
                x, y, z = get_property_value(line, "XYZ").split(",")
                atomDict[atomID] = (lammpsType, x, y, z)
                lmpIDxsdIDdict[atomID] = str(i)
                i += 1

            # Retriving molecules
            if ("<Bond" and "Connects") in line:
                ## Hotfix to exclude HBond
                if "HBond" in line:
                    continue
                atom1, atom2 = get_property_value(line, "Connects").split(",")
                atomPairs += [(atom1, atom2)]

            if "SpaceGroup" in line:
                aVector: list[str] = get_property_value(line, "AVector").split(",")
                bVector: list[str] = get_property_value(line, "BVector").split(",")
                cVector: list[str] = get_property_value(line, "CVector").split(",")

    print(atomPairs)
    molecules = find_molecules(atomPairs)
    crystalMatrix: list[float] = [float(aVector[0]), float(bVector[1]), float(cVector[2])]

    textBuffer.write("\n")
    for atomWithXsdID in atomDict.values():
        createAtom: str = f"create_atoms {atomWithXsdID[0]} single {float(atomWithXsdID[1]) * crystalMatrix[0]} {float(atomWithXsdID[2]) * crystalMatrix[1]} {float(atomWithXsdID[3]) * crystalMatrix[2]} remap yes\n"
        textBuffer.write(createAtom)

    print(molecules)
    for i, molecule in enumerate(molecules):
        textBuffer.write("\n")
        for atomWithXsdID in molecule:
            atomBindsToMolecule: str = f"set atom {lmpIDxsdIDdict.get(atomWithXsdID)} mol {i + 1} #{atomDict[atomWithXsdID][0]}\n"
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
    print(getCrystal(xsdExample))
