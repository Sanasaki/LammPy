from io import StringIO
from random import randint


class LammpsScriptFactory:
    labelAtoms = {
        1: "H[Nitric]",
        2: "N[Nitric]",
        3: "O1[Nitric]",
        4: "O2[Nitric]",
        5: "O[Water]",
        6: "H[Water]",
        7: "N[Nitrate]",
        8: "O[Nitrate]",
        9: "H[Hydronium]",
        10: "O[Hydronium]",
    }
    labelBonds = {
        1: "OH[Nitric]",
        2: "NO[Nitric]",
        3: "OH[Water]",
        4: "NO[Nitrate]",
        5: "OH[Hydronium]",
    }
    labelAngles = {
        1: "NOH[Nitric]",
        2: "ONO_1[Nitric]",
        3: "ONO_2[Nitric]",
        4: "ONO_3[Nitric]",
        5: "HOH[Water]",
        6: "ONO[Nitrate]",
        7: "HOH[Hydronium]",
    }

    @property
    def rng(self):
        return randint(0, 100_000)

    def __init__(self):
        self.units: str = "real"
        self.boundary: str = "p p p"
        self.timestep: float = 0.5
        self.atomStyle: str = "full"
        self.pairStyle: str = "lj/cut/coul/wolf 0.2 10"
        self.bondStyle: str = "zero"
        self.angleStyle: str = "zero"
        self.dihedralStyle: str = "none"
        self.improperStyle: str = "none"
        self.regionName: str = "system"
        self.regionType: str = "block"
        self.xlo: float = -20.0
        self.xhi: float = 20.0
        self.ylo: float = -20.0
        self.yhi: float = 20.0
        self.zlo: float = -20.0
        self.zhi: float = 20.0
        self.system: str = "#No system loaded"
        self.replicates: list[str] = []
        self.fixes: list[str] = []
        self.atomTypes: int = 10
        self.bondTypes: int = 5
        self.angleTypes: int = 7
        self.dihedralTypes: int = 1
        self.improperTypes: int = 1
        self.extraBondPerAtom: int = 3
        self.extraAnglePerAtom: int = 3
        self.extraSpecialPerAtom: int = 2
        self.extraDihedralPerAtom: int = 1
        self.extraImproperPerAtom: int = 1

    def buildJobAtPath(self, finalScriptPath: str) -> None:
        with open(finalScriptPath, "w") as file:
            file.write(self._getScript(name=finalScriptPath.split("/")[-1]))

    def _getScript(self, name: str) -> str:
        self._script = StringIO()

        self._script.write(f"units {self.units}\n")
        self._script.write(f"boundary {self.boundary}\n")
        self._script.write(f"timestep {self.timestep}\n")

        self._script.write(f"atom_style {self.atomStyle}\n")
        self._script.write(f"pair_style {self.pairStyle}\n")
        self._script.write(f"bond_style {self.bondStyle}\n")
        self._script.write(f"angle_style {self.angleStyle}\n")
        self._script.write(f"dihedral_style {self.dihedralStyle}\n")
        self._script.write(f"improper_style {self.improperStyle}\n")

        self._script.write(f"region {self.regionName} block {self.xlo} {self.xhi} {self.ylo} {self.yhi} {self.zlo} {self.zhi}\n")

        self._script.write(f"create_box {self.atomTypes} {self.regionName} &\n")
        self._script.write(f"bond/types {self.bondTypes} &\n")
        self._script.write(f"angle/types {self.angleTypes} &\n")
        self._script.write(f"dihedral/types {self.dihedralTypes} &\n")
        self._script.write(f"improper/types {self.improperTypes} &\n")
        self._script.write(f"extra/bond/per/atom {self.extraBondPerAtom} &\n")
        self._script.write(f"extra/angle/per/atom {self.extraAnglePerAtom} &\n")
        self._script.write(f"extra/special/per/atom {self.extraSpecialPerAtom} &\n")
        self._script.write(f"extra/dihedral/per/atom {self.extraDihedralPerAtom} &\n")
        self._script.write(f"extra/improper/per/atom {self.extraImproperPerAtom}\n")

        self._script.write("labelmap atom")
        for key, value in self.labelAtoms.items():
            self._script.write(f" &\n\t{key} {value}")
        self._script.write("\n")

        self._script.write("labelmap bond")
        for key, value in self.labelBonds.items():
            self._script.write(f" &\n\t{key} {value}")
        self._script.write("\n")

        self._script.write("labelmap angle")
        for key, value in self.labelAngles.items():
            self._script.write(f" &\n\t{key} {value}")
        self._script.write("\n")

        self._script.write(MASSES)
        self._script.write(FORCEFIELD)

        self._script.write(self.system)
        for replicate in self.replicates:
            self._script.write(replicate)

        self._script.write(f"""
variable H equal 4.184*enthalpy
variable Ec equal 4.184*ke
variable Ep equal 4.184*pe
variable Et equal 4.184*etotal
variable P equal 1.01325*press
variable T equal temp
variable d equal density
variable Vol equal vol
variable cpu_time equal cpu
variable sim_time equal time/1000

shell mkdir output
if $(is_os(^Windows)) then &
    "shell copy {name} output" &
    else &
    "shell cp {name} output"
""")

        self._script.write("""fix dataOutput all ave/time 1 100 100 v_sim_time v_cpu_time v_T v_P v_d v_Vol v_H file ./output/FixDataGlobal.csv &
        title2 "TimeStep VirtualTime(s) CpuTime(s) T(K) P(bar) Density(-) Volume(A^3) H(kJ/mol.at)"
        """)

        self._script.write("""
dump trajectory all xyz 100 ./output/Trajectory.xyz
dump_modify trajectory element H N O O O H N O H O
""")
        self._script.write("""
thermo 100
thermo_style custom step v_sim_time cpu cpuremain temp press density econserve ke pe enthalpy
thermo_modify norm yes
thermo_modify &
colname 1 "Timestep" &
colname 2 "Time(ps)" &
colname 3 "CpuTime(s)" &
colname 4 "TimeToEnd(s)" &
colname 5 "T(K)" &
colname 6 "P(atm)" &
colname 7 "Density(-)" &
colname 8 "Econserve(kcal/mol.at)" &
colname 9 "Ec(kcal/mol.at)" &
colname 10 "Ep(kcal/mol.at)" &
colname 11 "H(kcal/mol.at)"
""")

        for fix in self.fixes:
            self._script.write(f"{fix}\n")

        self._script.write("""
if $(is_os(^Windows)) then &
"shell copy log.lammps output" &
else &
"shell cp log.lammps output"
""")

        return self._script.getvalue()

    def loadSystem(self, lammpsSystem: str) -> None:
        setGroupCharges: str = """
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

set group WaterOxygenAtoms charge -0.8476
set group WaterHydrogenAtoms charge 0.4238
set group NitricHydrogenAtoms charge 0.497
set group NitricNitrogenAtoms charge 0.964
set group NitricOxygen1Atoms charge -0.445
set group NitricOxygen2Atoms charge -0.571
set group NitrateNitrogenAtoms charge 0.65
set group NitrateOxygenAtoms charge -0.55
set group HydroniumHydrogenAtoms charge 0.578
set group HydroniumOxygenAtoms charge -0.734

create_bonds many NitricOxygen2Atoms NitricHydrogenAtoms 1 0.95 1.0
create_bonds many NitricNitrogenAtoms NitricOxygen1Atoms 2 1.19 1.21
create_bonds many NitricNitrogenAtoms NitricOxygen2Atoms 2 1.19 1.21
create_bonds many WaterOxygenAtoms WaterHydrogenAtoms 3 0.98 1.1
create_bonds many NitrateNitrogenAtoms NitrateOxygenAtoms 4 1.25 1.3
create_bonds many HydroniumOxygenAtoms HydroniumHydrogenAtoms 5 0.98 1.1
"""
        self.system = lammpsSystem + setGroupCharges

    def replicate(self, x: int, y: int, z: int) -> None:
        self.replicates.append(f"replicate {x} {y} {z}\n")

    def addNVE(
        self,
        fixDurationPs: int,
    ) -> None:
        self.fixes.append(f"""
fix DataNVE all ave/time 1 100 100 v_sim_time v_cpu_time v_T v_P v_d v_Vol v_H file ./output/NVE.csv &
title2 "TimeStep VirtualTime(s) CpuTime(s) T(K) P(bar) Density(-) Volume(A^3) H(kJ/mol.at)"
                          
fix NVE all rigid/nve/small molecule
run $(1000*{fixDurationPs}/dt) #NVE for {fixDurationPs}ps
unfix NVE

unfix DataNVE
""")

    def addNVT(
        self,
        Temp1K: float,
        Temp2K: float,
        fixDurationPs: int,
    ) -> None:
        self.fixes.append(f"""
fix DataNVT all ave/time 1 100 100 v_sim_time v_cpu_time v_T v_P v_d v_Vol v_H file ./output/NVT-{int(Temp1K)}K.csv &
title2 "TimeStep VirtualTime(s) CpuTime(s) T(K) P(bar) Density(-) Volume(A^3) H(kJ/mol.at)"
                          
fix NVT all rigid/nvt/small molecule temp {Temp1K} {Temp2K} $(100*dt)
run $(1000*{fixDurationPs}/dt) #NVT from {Temp1K}K to {Temp2K}K in {fixDurationPs}ps
unfix NVT

unfix DataNVT
""")

    def addNPT(
        self,
        Temp1K: float,
        Temp2K: float,
        PressureBar: float,
        fixDurationPs: int,
    ) -> None:
        self.fixes.append(f"""
fix DataNPT all ave/time 1 100 100 v_sim_time v_cpu_time v_T v_P v_d v_Vol v_H file ./output/NPT-{int(Temp1K)}K-{int(PressureBar)}bar.csv &
title2 "TimeStep VirtualTime(s) CpuTime(s) T(K) P(bar) Density(-) Volume(A^3) H(kJ/mol.at)"

fix NPT all rigid/npt/small molecule temp {Temp1K} {Temp2K} $(100*dt) iso {PressureBar * 0.987} {PressureBar * 0.987} $(1000*dt)
run $(1000*{fixDurationPs}/dt) #NPT from {Temp1K}K to {Temp2K}K at {PressureBar}bar in {fixDurationPs}ps
unfix NPT

unfix DataNPT
""")


WATER_CRYSTAL: str = """
    # Water Crystal Conventional Cell (Ice-11)
    # Geometry optimized with CASTEP/RSCAN
        create_atoms H[Water] single -0.8164863157501667 2.0590504918748076 1.1638139454665992 remap yes
        create_atoms H[Water] single 1.3726310775769335 5.878564458111827 1.1638139454665992 remap yes
        create_atoms H[Water] single 0.8164863157501667 5.579977440599233 4.728053468694649 remap yes
        create_atoms H[Water] single 3.005603709077267 1.7604634743622123 4.728053468694649 remap yes
        create_atoms H[Water] single -0.8164863157501667 5.579977440599233 4.728053468694649 remap yes
        create_atoms H[Water] single 1.3726310775769335 1.7604634743622123 4.728053468694649 remap yes
        create_atoms H[Water] single 0.8164863157501667 2.0590504918748076 1.1638139454665992 remap yes
        create_atoms H[Water] single 3.005603709077267 5.878564458111827 1.1638139454665992 remap yes
        create_atoms O[Water] single 2.1891173933271 1.2513417109430336 1.7046670751815933 remap yes
        create_atoms H[Water] single 2.1891173933271 1.251879965657427 2.704666930322516 remap yes
        create_atoms H[Water] single 2.1891173933271 0.30834629735991537 1.3718613511304354 remap yes
        create_atoms O[Water] single 1.139254000162102e-16 5.070855677180054 1.7046670751815933 remap yes
        create_atoms H[Water] single 1.139254000162102e-16 5.071393931894447 2.704666930322516 remap yes
        create_atoms H[Water] single 2.2310390836507833e-16 4.127860263596936 1.3718613511304354 remap yes
        create_atoms O[Water] single 2.1891173933271 6.3876862215310055 5.268906598409643 remap yes
        create_atoms H[Water] single 2.1891173933271 6.387147966816612 6.268906453550565 remap yes
        create_atoms H[Water] single 2.1891173933271 7.330681635114123 4.936100874358486 remap yes
        create_atoms O[Water] single 4.044351700575464e-16 2.568172255293986 5.268906598409643 remap yes
        create_atoms H[Water] single 4.044351700575464e-16 2.5676340005795923 6.268906453550565 remap yes
        create_atoms H[Water] single 2.952566617086783e-16 3.511167668877103 4.936100874358486 remap yes
        create_atoms O[Water] single -3.009529317094888e-16 2.5860338814019683 0.9279468369540018 remap yes
        create_atoms O[Water] single 2.1891173933271 6.405547847638988 0.9279468369540018 remap yes
        create_atoms O[Water] single 3.009529317094888e-16 5.052994051072072 4.492186360182052 remap yes
        create_atoms O[Water] single 2.1891173933271 1.2334800848350518 4.492186360182052 remap yes

        set atom 1 mol 1 #H[Water]
        set atom 21 mol 1 #O[Water]
        set atom 7 mol 1 #H[Water]

        set atom 2 mol 2 #H[Water]
        set atom 22 mol 2 #O[Water]
        set atom 8 mol 2 #H[Water]

        set atom 3 mol 3 #H[Water]
        set atom 23 mol 3 #O[Water]
        set atom 5 mol 3 #H[Water]

        set atom 4 mol 4 #H[Water]
        set atom 24 mol 4 #O[Water]
        set atom 6 mol 4 #H[Water]

        set atom 10 mol 5 #H[Water]
        set atom 9 mol 5 #O[Water]
        set atom 11 mol 5 #H[Water]

        set atom 13 mol 6 #H[Water]
        set atom 12 mol 6 #O[Water]
        set atom 14 mol 6 #H[Water]

        set atom 16 mol 7 #H[Water]
        set atom 15 mol 7 #O[Water]
        set atom 17 mol 7 #H[Water]

        set atom 19 mol 8 #H[Water]
        set atom 18 mol 8 #O[Water]
        set atom 20 mol 8 #H[Water]

        change_box all x final 0.0 4.3782347866542 y final 0.0 7.63902793247404 z final 0.0 7.1284790464561
"""

NITRIC_CRYSTAL: str = """
    # Nitric Crystal Cell 
    # Geometry optimized with CASTEP/RSCAN
    # Optimization initialized from experimental cell of Allan et al. (2010). The crystal structures of the low-temperature and high-pressure polymorphs of nitric acid. Dalton Trans. Volume Number 39(15), 3736-3743. DOI: https://doi.org/10.1039/B923975H
        create_atoms H[Nitric] single 3.029263056686041 2.9836161709699334 17.100567829294825 remap yes
        create_atoms H[Nitric] single 0.03804715133467395 1.7056730791425403 15.01854538108791 remap yes
        create_atoms H[Nitric] single 0.0636631874671588 2.535364722848075 10.639572013843967 remap yes
        create_atoms H[Nitric] single 2.9651005736290252 3.8970596076535395 12.728609113594409 remap yes
        create_atoms H[Nitric] single 3.0787951950782784 5.418171767274984 8.767378908803247 remap yes
        create_atoms H[Nitric] single 6.0948992326213265 6.7253291526893335 10.87502436815919 remap yes
        create_atoms H[Nitric] single 6.137068490198665 5.907983889190034 15.1580133421852 remap yes
        create_atoms H[Nitric] single 3.125567979937462 8.098870886034726 13.016345447079084 remap yes
        create_atoms N[Nitric] single 3.009571514361713 4.842528192835353 17.215402820811384 remap yes
        create_atoms N[Nitric] single 0.006856172145505774 3.5675418345220016 15.058566430128325 remap yes
        create_atoms N[Nitric] single 0.01298790944833865 0.6812070448843923 10.808848515133025 remap yes
        create_atoms N[Nitric] single 3.0448737705631927 2.044547437212355 12.904609856757746 remap yes
        create_atoms N[Nitric] single 3.11244983760357 3.560724636469783 8.633646641234844 remap yes
        create_atoms N[Nitric] single 6.119011905747985 4.863495133688807 10.828952739618552 remap yes
        create_atoms N[Nitric] single 6.13215531177211 7.769052327937299 15.083641439406255 remap yes
        create_atoms N[Nitric] single 3.0841726879945526 6.237100783407409 12.981266739207618 remap yes
        create_atoms O1[Nitric] single 3.7682605314696187 4.680532733868712 18.136433352933743 remap yes
        create_atoms O1[Nitric] single 2.562660347037633 5.857178468743443 16.74579497539405 remap yes
        create_atoms O2[Nitric] single 2.578464934345625 3.66098096904876 16.589045615969436 remap yes
        create_atoms O1[Nitric] single -0.745362724066744 3.4219910302536336 15.987006384337612 remap yes
        create_atoms O1[Nitric] single 0.413511884423978 4.574049679793646 14.537982113399874 remap yes
        create_atoms O2[Nitric] single 0.4928257106034023 2.373981867724898 14.499065733939414 remap yes
        create_atoms O1[Nitric] single -0.7392096026920181 0.8987866388511214 11.7231167150881 remap yes
        create_atoms O1[Nitric] single 0.43089944551983383 -0.36222463440125147 10.37749889499562 remap yes
        create_atoms O2[Nitric] single 0.4739629157552642 1.8260884569868125 10.13764940009218 remap yes
        create_atoms O1[Nitric] single 3.7851870466728337 2.2767153878855075 13.825571271964298 remap yes
        create_atoms O2[Nitric] single 2.578793224148784 3.179278122569871 12.21946561277676 remap yes
        create_atoms O1[Nitric] single 2.6403321489057605 0.9935888259849293 12.478488114416864 remap yes
        create_atoms O1[Nitric] single 2.346437196109478 3.7258663190225927 7.719259929559632 remap yes
        create_atoms O1[Nitric] single 3.572774977659724 2.545424151527364 9.088665030370633 remap yes
        create_atoms O2[Nitric] single 3.535164748437684 4.7390829020635055 9.27162272850982 remap yes
        create_atoms O1[Nitric] single 6.876485647257146 5.009159727692499 9.904813610257852 remap yes
        create_atoms O1[Nitric] single 5.706578021322649 3.8568845401602325 11.344770077758907 remap yes
        create_atoms O2[Nitric] single 5.634863404396803 6.05709347406878 11.389948532980938 remap yes
        create_atoms O1[Nitric] single 6.88457795145532 7.621451035706739 14.155690339965021 remap yes
        create_atoms O1[Nitric] single 5.681331470675029 8.775616733970486 15.56635901904073 remap yes
        create_atoms O2[Nitric] single 5.71344585633327 6.578353065569369 15.70066739458119 remap yes
        create_atoms O1[Nitric] single 2.3422512302568004 6.383721327301792 12.044143098349691 remap yes
        create_atoms O1[Nitric] single 3.4795935692029425 5.230103772415533 13.509825541242375 remap yes
        create_atoms O2[Nitric] single 3.5725365427807265 7.429494198126892 13.541563533053269 remap yes
        create_atoms H[Nitric] single 3.109998894739789 -1.2417199153485674 0.15409148269347506 remap yes
        create_atoms H[Nitric] single 6.101214800091155 5.93100916546104 2.2361139309003883 remap yes
        create_atoms H[Nitric] single 6.075598763958672 6.7607008091665755 6.615087298144331 remap yes
        create_atoms H[Nitric] single 3.1741613777968047 8.12239569397204 4.526050198393888 remap yes
        create_atoms H[Nitric] single 3.0604667563475516 9.643507853593452 8.48728040318505 remap yes
        create_atoms H[Nitric] single 0.04436271880450237 2.4999930663708327 6.3796349438291085 remap yes
        create_atoms H[Nitric] single 0.0021934612271659946 1.6826478028715341 2.0966459698031 remap yes
        create_atoms H[Nitric] single 3.013693971488368 3.873534799716225 4.238313864909214 remap yes
        create_atoms N[Nitric] single 3.129690437064117 0.6171921065168533 0.03925649117690974 remap yes
        create_atoms N[Nitric] single 6.132405779280325 7.7928779208405015 2.196092881859973 remap yes
        create_atoms N[Nitric] single 6.126274041977489 4.906543131202896 6.445810796855273 remap yes
        create_atoms N[Nitric] single 3.0943881808626372 6.269883523530855 4.350049455230553 remap yes
        create_atoms N[Nitric] single 3.0268121138222597 7.7860607227882825 8.621012670753455 remap yes
        create_atoms N[Nitric] single 0.020250045677845985 0.6381590473703066 6.425706572369746 remap yes
        create_atoms N[Nitric] single 0.007106639653719005 3.5437162416187995 2.171017872582044 remap yes
        create_atoms N[Nitric] single 3.055089263431278 2.0117646970889087 4.273392572780682 remap yes
        create_atoms O1[Nitric] single 2.3710014199562113 0.4551966475502114 -0.8817740409454434 remap yes
        create_atoms O1[Nitric] single 3.5766016043881965 1.6318423824249424 0.5088643365942561 remap yes
        create_atoms O2[Nitric] single 3.5607970170802052 -0.5643551172697411 0.6656136960188622 remap yes
        create_atoms O1[Nitric] single 6.884624675492543 7.647327116572133 1.2676529276506878 remap yes
        create_atoms O1[Nitric] single 5.725750067001854 8.799385766112104 2.7166771985884264 remap yes
        create_atoms O2[Nitric] single 5.6464362408224265 6.599317954043398 2.755593578048884 remap yes
        create_atoms O1[Nitric] single 6.878471554117847 5.124122725169621 5.5315425969001994 remap yes
        create_atoms O1[Nitric] single 5.708362505905996 3.8631114519172463 6.8771604169926785 remap yes
        create_atoms O2[Nitric] single 5.665299035670568 6.051424543305313 7.117009911896119 remap yes
        create_atoms O1[Nitric] single 2.3540749047529963 6.502051474204007 3.4290880400240007 remap yes
        create_atoms O2[Nitric] single 3.5604687272770454 7.404614208888371 5.035193699211539 remap yes
        create_atoms O1[Nitric] single 3.49892980252007 5.21892491230343 4.776171197571435 remap yes
        create_atoms O1[Nitric] single 3.792824755316352 7.9512024053410935 9.535399382428666 remap yes
        create_atoms O1[Nitric] single 2.566486973766106 6.770760237845864 8.165994281617664 remap yes
        create_atoms O2[Nitric] single 2.604097202988146 8.964418988381981 7.983036583478478 remap yes
        create_atoms O1[Nitric] single -0.7372236958313159 0.7838236413739995 7.349845701730446 remap yes
        create_atoms O1[Nitric] single 0.4326839301031805 -0.36845154615826714 5.90988923422939 remap yes
        create_atoms O2[Nitric] single 0.5043985470290253 1.8317573877502793 5.864710779007359 remap yes
        create_atoms O1[Nitric] single -0.7453160000294912 3.396114949388239 3.0989689720232776 remap yes
        create_atoms O1[Nitric] single 0.4579304807508018 4.550280647651984 1.6883002929475648 remap yes
        create_atoms O2[Nitric] single 0.4258160950925607 2.3530169792508695 1.5539919174071066 remap yes
        create_atoms O1[Nitric] single 3.79701072116903 2.158385240983292 5.210516213638607 remap yes
        create_atoms O1[Nitric] single 2.6596683822228875 1.0047676860970323 3.7448337707459247 remap yes
        create_atoms O2[Nitric] single 2.5667254086451035 3.204158111808391 3.7130957789350285 remap yes

        set atom 19 mol 1 #O2[Nitric]
        set atom 1 mol 1 #H[Nitric]
        set atom 9 mol 1 #N[Nitric]
        set atom 17 mol 1 #O1[Nitric]
        set atom 18 mol 1 #O1[Nitric]

        set atom 22 mol 2 #O2[Nitric]
        set atom 2 mol 2 #H[Nitric]
        set atom 10 mol 2 #N[Nitric]
        set atom 20 mol 2 #O1[Nitric]
        set atom 21 mol 2 #O1[Nitric]

        set atom 25 mol 3 #O2[Nitric]
        set atom 3 mol 3 #H[Nitric]
        set atom 11 mol 3 #N[Nitric]
        set atom 23 mol 3 #O1[Nitric]
        set atom 24 mol 3 #O1[Nitric]

        set atom 27 mol 4 #O2[Nitric]
        set atom 4 mol 4 #H[Nitric]
        set atom 12 mol 4 #N[Nitric]
        set atom 26 mol 4 #O1[Nitric]
        set atom 28 mol 4 #O1[Nitric]

        set atom 31 mol 5 #O2[Nitric]
        set atom 5 mol 5 #H[Nitric]
        set atom 13 mol 5 #N[Nitric]
        set atom 30 mol 5 #O1[Nitric]
        set atom 29 mol 5 #O1[Nitric]

        set atom 34 mol 6 #O2[Nitric]
        set atom 6 mol 6 #H[Nitric]
        set atom 14 mol 6 #N[Nitric]
        set atom 33 mol 6 #O1[Nitric]
        set atom 32 mol 6 #O1[Nitric]

        set atom 37 mol 7 #O2[Nitric]
        set atom 15 mol 7 #N[Nitric]
        set atom 36 mol 7 #O1[Nitric]
        set atom 35 mol 7 #O1[Nitric]
        set atom 7 mol 7 #H[Nitric]

        set atom 40 mol 8 #O2[Nitric]
        set atom 16 mol 8 #N[Nitric]
        set atom 38 mol 8 #O1[Nitric]
        set atom 39 mol 8 #O1[Nitric]
        set atom 8 mol 8 #H[Nitric]

        set atom 59 mol 9 #O2[Nitric]
        set atom 49 mol 9 #N[Nitric]
        set atom 58 mol 9 #O1[Nitric]
        set atom 57 mol 9 #O1[Nitric]
        set atom 41 mol 9 #H[Nitric]

        set atom 62 mol 10 #O2[Nitric]
        set atom 42 mol 10 #H[Nitric]
        set atom 50 mol 10 #N[Nitric]
        set atom 60 mol 10 #O1[Nitric]
        set atom 61 mol 10 #O1[Nitric]

        set atom 65 mol 11 #O2[Nitric]
        set atom 51 mol 11 #N[Nitric]
        set atom 64 mol 11 #O1[Nitric]
        set atom 63 mol 11 #O1[Nitric]
        set atom 43 mol 11 #H[Nitric]

        set atom 67 mol 12 #O2[Nitric]
        set atom 52 mol 12 #N[Nitric]
        set atom 68 mol 12 #O1[Nitric]
        set atom 66 mol 12 #O1[Nitric]
        set atom 44 mol 12 #H[Nitric]

        set atom 71 mol 13 #O2[Nitric]
        set atom 45 mol 13 #H[Nitric]
        set atom 53 mol 13 #N[Nitric]
        set atom 70 mol 13 #O1[Nitric]
        set atom 69 mol 13 #O1[Nitric]

        set atom 74 mol 14 #O2[Nitric]
        set atom 54 mol 14 #N[Nitric]
        set atom 72 mol 14 #O1[Nitric]
        set atom 73 mol 14 #O1[Nitric]
        set atom 46 mol 14 #H[Nitric]

        set atom 77 mol 15 #O2[Nitric]
        set atom 55 mol 15 #N[Nitric]
        set atom 75 mol 15 #O1[Nitric]
        set atom 76 mol 15 #O1[Nitric]
        set atom 47 mol 15 #H[Nitric]

        set atom 80 mol 16 #O2[Nitric]
        set atom 56 mol 16 #N[Nitric]
        set atom 78 mol 16 #O1[Nitric]
        set atom 79 mol 16 #O1[Nitric]
        set atom 48 mol 16 #H[Nitric]

        change_box all x final 0.0 6.13926195142583 y final 0.0 8.450672172637 z final 0.0 17.2546593119883
"""

NAM_CRYSTAL: str = """
    # Nitric Acid Monohydrate (NAM) Crystal Cell
    # Geometry optimized with CASTEP/RSCAN
    # Optimization initialized from experimental cell of Lebrun et al. (2001). Kinetic behaviour investigations and crystal structure of nitric acid dihydrate. Acta Cryst B. Volume Number 57(1). 27-35. DOI: https://doi.org/10.1107/S0108768100014506
        create_atoms H[Hydronium] single 3.7672952833150797 2.165297947956077 2.5230642741943523 remap yes
        create_atoms H[Hydronium] single 2.210946378257988 1.6256462369147173 2.2449698970246432 remap yes
        create_atoms H[Hydronium] single 2.5289216559868324 3.251724598985764 2.327000234935788 remap yes
        create_atoms O[Hydronium] single 2.772683683438196 2.3370792271013494 2.737840892108796 remap yes
        create_atoms O[Nitrate] single 1.2373556436491262 0.4682053237594338 1.8289338500595516 remap yes
        create_atoms O[Nitrate] single -0.2719562171773933 1.9609015450450593 2.392228573831209 remap yes
        create_atoms O[Nitrate] single -0.8889343306425594 0.0814328733148276 1.4380757980717076 remap yes
        create_atoms N[Nitrate] single 0.025412139066260747 0.8364755602798499 1.8871916786545018 remap yes
        create_atoms H[Hydronium] single 0.9848772793908797 6.796589449000454 3.8953415611223976 remap yes
        create_atoms H[Hydronium] single -0.571471625666212 7.336241160041813 4.173435938292107 remap yes
        create_atoms H[Hydronium] single -0.2534963479373675 5.710162797970766 4.091405600380962 remap yes
        create_atoms O[Hydronium] single -0.00973432048600388 6.624808169855181 3.680564943207954 remap yes
        create_atoms O[Nitrate] single 4.009212238606976 8.453816866921978 4.571465820499312 remap yes
        create_atoms O[Nitrate] single 2.4999003777804556 6.961120645636354 4.008171096727655 remap yes
        create_atoms O[Nitrate] single 1.8829222643152848 8.840589317366584 4.962323872487156 remap yes
        create_atoms N[Nitrate] single 2.797268734024108 8.085546630401563 4.513207991904362 remap yes
        create_atoms H[Hydronium] single 3.7672952833150797 2.3156457505221884 5.732267191852728 remap yes
        create_atoms H[Hydronium] single 2.210946378257988 2.855297461563548 5.454172814683018 remap yes
        create_atoms H[Hydronium] single 2.5289216559868324 1.2292190994925012 5.5362031525941635 remap yes
        create_atoms O[Hydronium] single 2.772683683438196 2.143864471376916 5.947043809767171 remap yes
        create_atoms O[Nitrate] single 1.2373556436491262 4.012738374718831 5.038136767717926 remap yes
        create_atoms O[Nitrate] single -0.27195621717739166 2.5200421534331974 5.601431491489584 remap yes
        create_atoms O[Nitrate] single -0.8889343306425594 4.399510825163436 4.647278715730082 remap yes
        create_atoms N[Nitrate] single 0.025412139066261302 3.644468138198415 5.0963945963128765 remap yes
        create_atoms H[Hydronium] single 0.9848772793908797 6.646241646434342 0.6861386434640226 remap yes
        create_atoms H[Hydronium] single -0.571471625666212 6.106589935392983 0.964233020633732 remap yes
        create_atoms H[Hydronium] single -0.2534963479373675 7.73266829746403 0.8822026827225865 remap yes
        create_atoms O[Hydronium] single -0.00973432048600388 6.818022925579615 0.47136202554957785 remap yes
        create_atoms O[Nitrate] single 4.019773647573326 4.949149022237709 1.3802690675988236 remap yes
        create_atoms O[Nitrate] single 2.5104617867468058 6.441845243523334 0.8169743438271657 remap yes
        create_atoms O[Nitrate] single 1.8934836732816347 4.562376571793094 1.7711271195866671 remap yes
        create_atoms N[Nitrate] single 2.8078301429904635 5.3174192587581155 1.3220112390038732 remap yes

        set atom 4 mol 1 #O[Hydronium]
        set atom 1 mol 1 #H[Hydronium]
        set atom 2 mol 1 #H[Hydronium]
        set atom 3 mol 1 #H[Hydronium]

        set atom 8 mol 2 #N[Nitrate]
        set atom 5 mol 2 #O[Nitrate]
        set atom 6 mol 2 #O[Nitrate]
        set atom 7 mol 2 #O[Nitrate]

        set atom 12 mol 3 #O[Hydronium]
        set atom 9 mol 3 #H[Hydronium]
        set atom 10 mol 3 #H[Hydronium]
        set atom 11 mol 3 #H[Hydronium]

        set atom 16 mol 4 #N[Nitrate]
        set atom 13 mol 4 #O[Nitrate]
        set atom 14 mol 4 #O[Nitrate]
        set atom 15 mol 4 #O[Nitrate]

        set atom 20 mol 5 #O[Hydronium]
        set atom 17 mol 5 #H[Hydronium]
        set atom 18 mol 5 #H[Hydronium]
        set atom 19 mol 5 #H[Hydronium]

        set atom 24 mol 6 #N[Nitrate]
        set atom 21 mol 6 #O[Nitrate]
        set atom 22 mol 6 #O[Nitrate]
        set atom 23 mol 6 #O[Nitrate]

        set atom 28 mol 7 #O[Hydronium]
        set atom 25 mol 7 #H[Hydronium]
        set atom 26 mol 7 #H[Hydronium]
        set atom 27 mol 7 #H[Hydronium]

        set atom 32 mol 8 #N[Nitrate]
        set atom 29 mol 8 #O[Nitrate]
        set atom 30 mol 8 #O[Nitrate]
        set atom 31 mol 8 #O[Nitrate]

        change_box all x final 0.0 5.5648360078484 y final 0.0 8.96188739695653 z final 0.0 6.41840583531675
"""

MASSES: str = """
    # Hydrogens
    mass H[Water] 1.008
    mass H[Nitric] 1.008
    mass H[Hydronium] 1.008

    # Nitrogens
    mass N[Nitric] 14.0067
    mass N[Nitrate] 14.0067

    # Oxygens
    mass O[Water] 15.9994
    mass O1[Nitric] 15.9994
    mass O2[Nitric] 15.9994
    mass O[Hydronium] 15.9994
    mass O[Nitrate] 15.9994
"""

FORCEFIELD: str = """
    # Forcefield parameters from Cordeiro et al. (2020). Parametrization and Molecular Dynamics Simulations of Nitrogen Oxyanions and Oxyacids for Applications in Atmospheric and Biomolecular Sciences. J. Phys. Chem. B. Volume Number 124(6), 1082-1089. DOI: https://doi.org/10.1021/acs.jpcb.9b08172
        # H2O
        pair_coeff H[Water]     *               0    0
        pair_coeff O[Water]     O[Water]        0.16275 3.16435   # OW water Oxygen type

        # HNO3
        pair_coeff H[Nitric]    *               0      0
        pair_coeff N[Nitric]    N[Nitric]       0.1407 2.8569   # NQ reactive type
        pair_coeff O1[Nitric]   O1[Nitric]      0.0979 3.0352   # OQ ozone Oxygen type
        pair_coeff O2[Nitric]   O2[Nitric]      0.1720 3.0384   # OP peroxyde oxygen type

        # H3O+
        pair_coeff H[Hydronium] *               0      0
        pair_coeff O[Hydronium] O[Hydronium]    0.1554 3.1655   # OW water Oxygen type

        # NO3-
        pair_coeff N[Nitrate]   N[Nitrate]      0.1047 3.3411   # NR aromatic type
        pair_coeff O[Nitrate]   O[Nitrate]      0.0903 3.3828   # OM carboxyl Oxygen type

        ### Bonds ###        
        bond_coeff OH[Water]    0

        bond_coeff OH[Nitric]   0
        bond_coeff NO[Nitric]   0

        bond_coeff OH[Hydronium]   0

        bond_coeff NO[Nitrate]   0

        ### Angles ##
        angle_coeff HOH[Water]      109.47

        angle_coeff ONO_1[Nitric]   130.27
        angle_coeff ONO_2[Nitric]   113.85
        angle_coeff ONO_3[Nitric]   115.88
        angle_coeff NOH[Nitric]     102.15

        angle_coeff HOH[Hydronium]  113.0

        angle_coeff ONO[Nitrate]    120.0
"""


if __name__ == "__main__":
    script = LammpsScriptFactory()
    script.loadSystem(WATER_CRYSTAL)
    script.replicate(5, 3, 3)
    script.replicate(2, 2, 2)
    # script.addNVE(fixDurationPs=50)
    script.addNVT(Temp1K=10, Temp2K=10, fixDurationPs=10)
    Temperatures = [
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        90,
        100,
        110,
        120,
        130,
        140,
        150,
        160,
        170,
        180,
        190,
        200,
        210,
        220,
        230,
        240,
        250,
        260,
        270,
        280,
        290,
        300,
    ]
    for temp in Temperatures:
        script.addNPT(Temp1K=temp, Temp2K=temp, PressureBar=1.0, fixDurationPs=10)
    script.buildJobAtPath("C:/Archives/Thesis/Lab/VisualStudioCode/Packages/LammPy/Water-NPT.lammps")
