import re
from collections import defaultdict
from io import StringIO

xsdExample: str = r"C:\Users\JL252842\Documents\Thesis\Lab\MaterialsStudio\CreatingLammpsCrystals_Files\Documents\NAMCrystalEditedForLAMMPS.xsd"


def getPropertyValue(line: str, propertyName: str) -> str:
    match = re.search(rf'{propertyName}="([^"]+)"', line)
    if match:
        return match.group(1)
    return " "


def find_molecules(bonds):
    """
    Regroupe les atomes en molécules basées sur les liaisons.
    bonds : liste de tuples (atom1, atom2) représentant les liaisons.
    Retourne une liste de molécules (listes d'atomes).
    """
    # Construire un graphe sous forme de dictionnaire d'adjacence
    graph = defaultdict(set)
    for a1, a2 in bonds:
        graph[a1].add(a2)
        graph[a2].add(a1)

    # Fonction récursive pour explorer une molécule
    def dfs(atom, molecule):
        if atom in visited:
            return
        visited.add(atom)
        molecule.append(atom)
        for neighbor in graph[atom]:
            dfs(neighbor, molecule)

    # Recherche des molécules en parcourant le graphe
    visited = set()
    molecules = []

    for atom in graph:
        if atom not in visited:
            molecule = []
            dfs(atom, molecule)
            molecules.append(molecule)

    return molecules


def getCrystal(xsdFilePath: str) -> str:
    textBuffer = StringIO()
    with open(xsdFilePath, "r") as xsdFile:
        atomDict: dict[str, tuple[str, str, str, str]] = {}
        lmpIDxsdIDdict: dict[str, str] = {}
        moleculeDict: dict[str, list[str]] = {}
        i = 1
        bondList = []
        for line in xsdFile:
            # Retriving atoms
            if ("<Atom3d" and "UserID") in line:
                atomID: str = getPropertyValue(line, "ID")
                lammpsType: str = getPropertyValue(line, "Name")
                x, y, z = getPropertyValue(line, "XYZ").split(",")
                atomDict[atomID] = (lammpsType, x, y, z)
                lmpIDxsdIDdict[atomID] = str(i)
                i += 1

            # Retriving molecules
            if ("<Bond" and "Connects") in line:
                if "HBond" in line:
                    continue
                atom1, atom2 = getPropertyValue(line, "Connects").split(",")
                bondList += [(atom1, atom2)]
                # found = False
                # for molecule in moleculeDict.values():
                #     if atom1 in molecule or atom2 in molecule:
                #         molecule.extend([atom for atom in [atom1, atom2] if atom not in molecule])
                #         found = True
                #         break
                # if not found:
                #     moleculeDict[f"{len(moleculeDict) + 1}"] = [atom1, atom2]

            if "SpaceGroup" in line:
                aVector: list[str] = getPropertyValue(line, "AVector").split(",")
                bVector: list[str] = getPropertyValue(line, "BVector").split(",")
                cVector: list[str] = getPropertyValue(line, "CVector").split(",")

    print(bondList)
    molecules = find_molecules(bondList)

    crystalMatrix = [float(aVector[0]), float(bVector[1]), float(cVector[2])]
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
