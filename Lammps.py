import os

from lammps import PyLammps
from LammPy.XSDtoLMP import getCrystal


def NVTtransition(
    nvtFixName: str,
    startTempKelvin: float,
    finalTempKelvin: float,
    fixDurationPicoseconds: int,
) -> str:
    return f"""
    fix {nvtFixName} all rigid/nvt/small molecule temp {startTempKelvin} {finalTempKelvin} $(100*dt)
    run $({int(fixDurationPicoseconds * 1000)}/dt) #{fixDurationPicoseconds} ps from {startTempKelvin}K to {finalTempKelvin}K
    unfix {nvtFixName}
    """


def NptFix(
    nptFixName: str,
    startTempKelvin: float,
    finalTempKelvin: float,
    startPressureAtm: float,
    finalPressureAtm: float,
    fixDurationPicoseconds: int,
) -> str:
    return f"""
    fix {nptFixName} all rigid/npt/small molecule temp {startTempKelvin} {finalTempKelvin} $(100*dt) iso {startPressureAtm} {finalPressureAtm} $(1000*dt)
    run $({int(fixDurationPicoseconds * 1000)}/dt) #{fixDurationPicoseconds} ps from {startTempKelvin}K to {finalTempKelvin}K
    unfix {nptFixName}
    """


Masses: str = """
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

ForceField: str = """
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

NitricAcidMonohydrateAtoms: str = """
    variable x equal 5.5648
    variable y equal 8.9619
    variable z equal 6.4184

    create_atoms H[Hydronium] single $(0.67698*v_x) $(0.24161*v_y) $(0.39310*v_z) remap yes #H1
    create_atoms H[Hydronium] single $(0.39731*v_x) $(0.18140*v_y) $(0.34977*v_z) remap yes #H2
    create_atoms H[Hydronium] single $(0.45445*v_x) $(0.36284*v_y) $(0.36255*v_z) remap yes #H3
    create_atoms O[Hydronium] single $(0.49825*v_x) $(0.26078*v_y) $(0.42656*v_z) remap yes #O4

    set atom 1 mol 1
    set atom 2 mol 1
    set atom 3 mol 1
    set atom 4 mol 1

    create_atoms O[Nitrate] single $(0.23008*v_x) $(0.05078*v_y) $(0.28464*v_z) remap yes #O5
    create_atoms O[Nitrate] single $(-0.05077*v_x) $(0.22325*v_y) $(0.37552*v_z) remap yes #O6
    create_atoms O[Nitrate] single $(-0.16557*v_x) $(0.00609*v_y) $(0.22158*v_z) remap yes #O7
    create_atoms N[Nitrate] single $(0.00457*v_x) $(0.09333*v_y) $(0.29404*v_z) remap yes #N8

    set atom 5 mol 2
    set atom 6 mol 2
    set atom 7 mol 2
    set atom 8 mol 2

    create_atoms H[Hydronium] single $(0.17698*v_x) $(0.75839*v_y) $(0.60690*v_z) remap yes #H9
    create_atoms H[Hydronium] single $(-0.10269*v_x) $(0.81860*v_y) $(0.65023*v_z) remap yes #H10
    create_atoms H[Hydronium] single $(-0.04555*v_x) $(0.63716*v_y) $(0.63745*v_z) remap yes #H11
    create_atoms O[Hydronium] single $(-0.00175*v_x) $(0.73922*v_y) $(0.57344*v_z) remap yes #O12

    set atom 9 mol 3
    set atom 10 mol 3
    set atom 11 mol 3
    set atom 12 mol 3

    create_atoms O[Nitrate] single $(0.73008*v_x) $(0.94922*v_y) $(0.71536*v_z) remap yes #O13
    create_atoms O[Nitrate] single $(0.44923*v_x) $(0.77675*v_y) $(0.62448*v_z) remap yes #O14
    create_atoms O[Nitrate] single $(0.33443*v_x) $(0.99391*v_y) $(0.77842*v_z) remap yes #O15
    create_atoms N[Nitrate] single $(0.50457*v_x) $(0.90667*v_y) $(0.70596*v_z) remap yes #N16

    set atom 13 mol 4
    set atom 14 mol 4
    set atom 15 mol 4
    set atom 16 mol 4

    create_atoms H[Hydronium] single $(0.67698*v_x) $(0.25839*v_y) $(0.89310*v_z) remap yes #H17
    create_atoms H[Hydronium] single $(0.39731*v_x) $(0.31860*v_y) $(0.84977*v_z) remap yes #H18
    create_atoms H[Hydronium] single $(0.45445*v_x) $(0.13716*v_y) $(0.86255*v_z) remap yes #H19
    create_atoms O[Hydronium] single $(0.49825*v_x) $(0.23922*v_y) $(0.92656*v_z) remap yes #O20

    set atom 17 mol 5
    set atom 18 mol 5
    set atom 19 mol 5
    set atom 20 mol 5

    create_atoms O[Nitrate] single $(0.23008*v_x) $(0.44922*v_y) $(0.78464*v_z) remap yes #O21
    create_atoms O[Nitrate] single $(-0.05077*v_x) $(0.27675*v_y) $(0.87552*v_z) remap yes #O22
    create_atoms O[Nitrate] single $(-0.16557*v_x) $(0.49391*v_y) $(0.72158*v_z) remap yes #O23
    create_atoms N[Nitrate] single $(0.00457*v_x) $(0.40667*v_y) $(0.79404*v_z) remap yes #N24

    set atom 21 mol 6
    set atom 22 mol 6
    set atom 23 mol 6
    set atom 24 mol 6

    create_atoms H[Hydronium] single $(0.17698*v_x) $(0.74161*v_y) $(0.10690*v_z) remap yes #H25
    create_atoms H[Hydronium] single $(-0.10269*v_x) $(0.68140*v_y) $(0.15023*v_z) remap yes #H26
    create_atoms H[Hydronium] single $(-0.04555*v_x) $(0.86284*v_y) $(0.13745*v_z) remap yes #H27
    create_atoms O[Hydronium] single $(-0.00175*v_x) $(0.76078*v_y) $(0.07344*v_z) remap yes #O28

    set atom 25 mol 7
    set atom 26 mol 7
    set atom 27 mol 7
    set atom 28 mol 7

    create_atoms O[Nitrate] single $(0.73008*v_x) $(0.55078*v_y) $(0.21536*v_z) remap yes #O29
    create_atoms O[Nitrate] single $(0.44923*v_x) $(0.72325*v_y) $(0.12448*v_z) remap yes #O30
    create_atoms O[Nitrate] single $(0.33443*v_x) $(0.50609*v_y) $(0.27842*v_z) remap yes #O31
    create_atoms N[Nitrate] single $(0.50457*v_x) $(0.59333*v_y) $(0.20596*v_z) remap yes #N32

    set atom 29 mol 8
    set atom 30 mol 8
    set atom 31 mol 8
    set atom 32 mol 8

    group NitrateNitrogenAtoms type 7
    group NitrateOxygenAtoms type 8
    group HydroniumHydrogenAtoms type 9
    group HydroniumOxygenAtoms type 10

    set group NitrateNitrogenAtoms charge 0.65
    set group NitrateOxygenAtoms charge -0.55
    set group HydroniumHydrogenAtoms charge 0.578
    set group HydroniumOxygenAtoms charge -0.734

    create_bonds many NitrateNitrogenAtoms NitrateOxygenAtoms 4 1.3 1.4
    create_bonds many HydroniumOxygenAtoms HydroniumHydrogenAtoms 5 0.98 1.1

    change_box all x final 0 $x y final 0 $y z final 0 $z

    replicate 10 6 5

    group HNO3 type 1 2 3 4
    group H2O type 5 6
    group NO3 type 7 8
    group H3O type 9 10
"""

ThermoVariables: str = """
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
"""


def RelaxSystem(lmpScript: PyLammps, temperatureKelvin: float, durationPicosecond: int) -> None:
    lmpScript.append_cmd_history(
        NVTtransition(
            nvtFixName="NvtRelaxation",
            startTempKelvin=temperatureKelvin,
            finalTempKelvin=temperatureKelvin,
            fixDurationPicoseconds=durationPicosecond,
        )
    )
    lmpScript.append_cmd_history(
        NptFix(
            nptFixName="NptRelaxation",
            startTempKelvin=temperatureKelvin,
            finalTempKelvin=temperatureKelvin,
            fixDurationPicoseconds=durationPicosecond,
            startPressureAtm=1.0,
            finalPressureAtm=1.0,
        )
    )
    return


def MeasureCp(
    lmpScript: PyLammps,
    temperatureKelvin: float,
    diffTempKelvin: int = 5,
    rampTimePicosecond: int = 100,
    pressureBar: float = 1.0,
) -> None:
    lowTemp: float = temperatureKelvin - diffTempKelvin
    highTemp: float = temperatureKelvin + diffTempKelvin

    if lowTemp <= 0 or (lowTemp + diffTempKelvin <= 0):
        lowTemp = 1.0 + diffTempKelvin
        highTemp = 1.0 + 2 * diffTempKelvin
    pressureAtm: float = pressureBar / 1.01325

    # Setting up the system at temperature T
    # lmpScript.append_cmd_history(NVTtransition("NvtEquilbriumT", temperatureKelvin, lowTemp, 5))
    DataFixName: str = f"DataWritingT{int(temperatureKelvin)}"
    lmpScript.append_cmd_history(DataFix(DataFixName))
    lmpScript.append_cmd_history(
        NptFix(
            nptFixName=f"NptEquilibrationT{int(lowTemp)}",
            startTempKelvin=lowTemp,
            finalTempKelvin=lowTemp,
            startPressureAtm=pressureAtm,
            finalPressureAtm=pressureAtm,
            fixDurationPicoseconds=10,
        )
    )
    # lmpScript.append_cmd_history(
    #     NptFix(
    #         "NptHeating",
    #         lowTemp,
    #         highTemp,
    #         pressureAtm,
    #         pressureAtm,
    #         5,
    #     )
    # )
    lmpScript.append_cmd_history(
        NptFix(
            nptFixName=f"NptEquilibrationT{int(highTemp)}",
            startTempKelvin=highTemp,
            finalTempKelvin=highTemp,
            startPressureAtm=pressureAtm,
            finalPressureAtm=pressureAtm,
            fixDurationPicoseconds=10,
        )
    )
    lmpScript.append_cmd_history(f"unfix {DataFixName}")
    return


def DataFix(
    dataFixName: str,
) -> str:
    return f"""
fix {dataFixName} all ave/time 1 100 100 v_sim_time v_cpu_time v_T v_P v_d v_Vol v_H file ./output/{dataFixName}.csv &
title2 "TimeStep VirtualTime(s) CpuTime(s) T(K) P(bar) Density(-) Volume(A^3) H(kJ/mol(atom))"
"""


def main() -> None:
    # NOTE: argv[0] is set by the lammps class constructor
    args = ["-log", "none"]
    timestep = 0.5

    # create LAMMPS instance
    lammpsSimulation = PyLammps(cmdargs=args)
    lammpsSimulation.append_cmd_history("units real")
    lammpsSimulation.append_cmd_history("boundary p p p")
    lammpsSimulation.append_cmd_history(f"timestep {timestep}")
    lammpsSimulation.append_cmd_history("atom_style full")
    lammpsSimulation.append_cmd_history("pair_style lj/cut/coul/wolf 0.2 10")
    lammpsSimulation.append_cmd_history("bond_style zero")
    lammpsSimulation.append_cmd_history("angle_style zero")
    lammpsSimulation.append_cmd_history("dihedral_style none")
    lammpsSimulation.append_cmd_history("improper_style none")

    # Region
    boxEdgeSize: float = 20.0
    xlo: float = -boxEdgeSize
    xhi: float = boxEdgeSize
    ylo: float = -boxEdgeSize
    yhi: float = boxEdgeSize
    zlo: float = -boxEdgeSize
    zhi: float = boxEdgeSize
    lammpsSimulation.append_cmd_history(f"region box block {xlo} {xhi} {ylo} {yhi} {zlo} {zhi}")

    atomTypes: int = 10
    bondTypes: int = 5
    angleTypes: int = 7
    dihedralTypes: int = 1
    improperTypes: int = 1
    extraBondPerAtom: int = 3
    extraAnglePerAtom: int = 3
    extraSpecialPerAtom: int = 0
    extraDihedralPerAtom: int = 1
    extraImproperPerAtom: int = 1

    lammpsSimulation.append_cmd_history(f"create_box {atomTypes} box &")
    lammpsSimulation.append_cmd_history(f"bond/types {bondTypes} &")
    lammpsSimulation.append_cmd_history(f"angle/types {angleTypes} &")
    lammpsSimulation.append_cmd_history(f"dihedral/types {dihedralTypes} &")
    lammpsSimulation.append_cmd_history(f"improper/types {improperTypes} &")
    lammpsSimulation.append_cmd_history(f"extra/bond/per/atom {extraBondPerAtom} &")
    lammpsSimulation.append_cmd_history(f"extra/angle/per/atom {extraAnglePerAtom} &")
    lammpsSimulation.append_cmd_history(f"extra/special/per/atom {extraSpecialPerAtom} &")
    lammpsSimulation.append_cmd_history(f"extra/dihedral/per/atom {extraDihedralPerAtom} &")
    lammpsSimulation.append_cmd_history(f"extra/improper/per/atom {extraImproperPerAtom}")

    labelMaps: str = """
labelmap atom &
    1 H[Nitric] &
    2 N[Nitric] &
    3 O1[Nitric] &
    4 O2[Nitric] &
    5 O[Water] &
    6 H[Water] &
    7 N[Nitrate] &
    8 O[Nitrate] &
    9 H[Hydronium] &
    10 O[Hydronium]
labelmap bond &
    1 OH[Nitric] &
    2 NO[Nitric] &
    3 OH[Water] &
    4 NO[Nitrate] &
    5 OH[Hydronium]
labelmap angle &
    1 NOH[Nitric] &
    2 ONO_1[Nitric] &
    3 ONO_2[Nitric] &
    4 ONO_3[Nitric] &
    5 HOH[Water] &
    6 ONO[Nitrate] &
    7 HOH[Hydronium]
    """

    lammpsSimulation.append_cmd_history(labelMaps)
    lammpsSimulation.append_cmd_history(Masses)
    lammpsSimulation.append_cmd_history(ForceField)

    # cellParameterA: float = 5.5648
    # cellParameterB: float = 8.9619
    # cellParameterC: float = 6.4184

    # for atom in atoms:
    #     lineToWrite: str = f"create_atoms "
    #     lammpsSimulation.append_cmd_history(lineToWrite)
    crystalFile: str = r"C:\Archives\Thesis\Lab\ReferenceFiles\NitricAcidCrystalEditedForLAMMPS-P1.xsd"
    outputDir: str = crystalFile.split("\\")[-1].split(".")[0]

    lammpsSimulation.append_cmd_history(getCrystal(crystalFile))
    lammpsSimulation.append_cmd_history("replicate 6 6 6")
    lammpsSimulation.append_cmd_history(ThermoVariables)

    MakeOutputDir: str = f"""
shell mkdir output
if $(is_os(^Windows)) then &
    "shell copy {outputDir}.lammps output" &
    else &
    "shell cp {outputDir}.lammps output"
"""
    lammpsSimulation.append_cmd_history(MakeOutputDir)

    lammpsSimulation.append_cmd_history(DataFix("FixDataGlobal"))

    DumpTrajectory: str = """
dump trajectory all xyz 200 ./output/Trajectory.xyz
dump_modify trajectory element H N O O O H N O H O
"""

    lammpsSimulation.append_cmd_history(DumpTrajectory)

    ThermoStyle: str = """
thermo 100
thermo_style custom step v_sim_time cpu cpuremain temp press density vol ke pe enthalpy
thermo_modify norm yes
thermo_modify &
colname 1 "Timestep" &
colname 2 "Time(ps)" &
colname 3 "CpuTime(s)" &
colname 4 "TimeToEnd(s)" &
colname 5 "T(K)" &
colname 6 "P(atm)" &
colname 7 "Density(-)" &
colname 8 "Volume(A^3)" &
colname 9 "Ec(kcal/mol.at)" &
colname 10 "Ep(kcal/mol.at)" &
colname 11 "H(kcal/mol.at)"
"""
    lammpsSimulation.append_cmd_history(ThermoStyle)

    WriteData: str = """
write_data ./output/restart.data
write_dump all xyz ./output/FinalState.xyz modify element H N O O O H N O H O
"""

    RelaxSystem(lmpScript=lammpsSimulation, temperatureKelvin=5, durationPicosecond=5)
    coldTemperatures: list[int] = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    hotTemperatures: list[int] = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300]
    for T in coldTemperatures:
        MeasureCp(
            lmpScript=lammpsSimulation,
            temperatureKelvin=T,
            diffTempKelvin=2,
            pressureBar=1.0,
        )
    for T in hotTemperatures:
        MeasureCp(
            lmpScript=lammpsSimulation,
            temperatureKelvin=T,
            diffTempKelvin=5,
            pressureBar=1.0,
        )
    lammpsSimulation.append_cmd_history(WriteData)

    Repetability: str = """
if $(is_os(^Windows)) then &
"shell copy log.lammps output" &
else &
"shell cp log.lammps output"
"""

    lammpsSimulation.append_cmd_history(Repetability)

    # explicitly close and delete LAMMPS instance (optional)
    os.system(f"mkdir {outputDir}")
    lammpsSimulation.write_script(f"{outputDir}/{outputDir}.lammps")
    lammpsSimulation.close()


if __name__ == "__main__":
    main()
