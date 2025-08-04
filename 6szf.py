from openmm.app import PDBFile, Modeller, ForceField
from openmm.openmm import Platform


pdb = PDBFile("6szf.pdb")
model = Modeller(pdb.topology, pdb.positions)
# model.addHydrogens()
platform = Platform.getPlatformByName("HIP")
forcefield = ForceField("amber14-all.xml")

system = forcefield.createSystem(model.topology)

# Save the system to a PDB file (this won't save forces, just atom positions)
with open("system.pdb", "w") as f:
    PDBFile.writeFile(model.topology, model.positions, f)

print(system)
