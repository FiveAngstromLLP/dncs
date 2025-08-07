import sys
from openmm.app import PDBFile, Modeller, ForceField
from openmm import Platform, VerletIntegrator, Context
from openmm.unit import elementary_charge, kilojoule_per_mole

pdb = PDBFile(sys.argv[1])
forcefield = ForceField("amber99sb.xml")
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(modeller.topology)
for i, force in enumerate(system.getForces()): force.setForceGroup(i)
platform = Platform.getPlatformByName("CPU")
integrator = VerletIntegrator(0.001)
context = Context(system, integrator, platform)
context.setPositions(modeller.positions)

print("Individual energy terms (kJ/mol)\n", "-" * 50)
force_energies = {}
for i, force in enumerate(system.getForces()):
    force_name = force.__class__.__name__
    state = context.getState(getEnergy=True, groups=2**i)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    force_energies[force_name] = energy
    print(f"{force_name:<25} {energy:12.3f} kJ/mol")

state = context.getState(getEnergy=True)
openmm_total = state.getPotentialEnergy()
sum_individual = sum(force_energies.values())

print(f"{'-'*50}\n{'Individual energy sum':<25} {sum_individual:12.3f} kJ/mol")
print(f"{'Total system energy':<25} {openmm_total}\n{'-'*50}")
