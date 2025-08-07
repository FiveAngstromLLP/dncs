#!/usr/bin/env python3
"""
Debug OpenMM forcefield application - FIXED VERSION for unit conflicts
"""
import sys
# import numpy as np  # Not needed anymore
from openmm.app import PDBFile, Modeller, ForceField
from openmm import *
from openmm.unit import *

# Use Python's built-in sum to avoid OpenMM's unit system conflicts
import builtins

def count_residue_atoms(residue):
    """Count atoms in a residue by iterating over them."""
    return builtins.sum(1 for atom in residue.atoms())

pdb_file = sys.argv[1] if len(sys.argv) > 1 else "yggfm.pdb"
pdb = PDBFile(pdb_file)

print("=== STRUCTURE DEBUG ===")
print(f"Original atoms: {pdb.topology.getNumAtoms()}")
print("Residues:")
for i, residue in enumerate(pdb.topology.residues()):
    atom_count = count_residue_atoms(residue)
    print(f"  {i}: {residue.name} - {atom_count} atoms")

# Check what's in your PDB - show ALL atoms to see the structure
print(f"\nAll atoms in structure:")
for i, atom in enumerate(pdb.topology.atoms()):
    print(f"  {i:2d}: {atom.name:>4s} ({atom.element.symbol:>2s}) in {atom.residue.name}")

print(f"\nBonds in original PDB:")
bonds = list(pdb.topology.bonds())
print(f"  Total bonds: {len(bonds)}")
for i, (atom1, atom2) in enumerate(bonds[:10]):  # Show first 10 bonds
    print(f"    {i}: {atom1.name}-{atom2.name}")
if len(bonds) > 10:
    print(f"    ... and {len(bonds)-10} more")

print("\n" + "="*50)
print("=== FORCEFIELD APPLICATION ===")

try:
    forcefield = ForceField("amber99sb.xml")
    print("✓ Forcefield loaded successfully")
except Exception as e:
    print(f"✗ Forcefield loading failed: {e}")
    sys.exit(1)

modeller = Modeller(pdb.topology, pdb.positions)

try:
    print("Adding hydrogens...")
    modeller.addHydrogens(forcefield)
    print(f"✓ After adding H: {modeller.topology.getNumAtoms()} atoms")
except Exception as e:
    print(f"✗ Adding hydrogens failed: {e}")
    print("This usually means residue names don't match the forcefield")
    
    # Show residue names that might be problematic
    print("\nResidue names in your PDB:")
    for residue in pdb.topology.residues():
        print(f"  - {residue.name}")
    
    # Show standard amino acid codes for comparison
    print("\nExpected AMBER residue names:")
    standard_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                        "THR", "TRP", "TYR", "VAL", "ACE", "NME"]
    print(f"  Standard: {', '.join(standard_residues[:10])}...")
    sys.exit(1)

# Try to create the system
try:
    print("Creating system...")
    system = forcefield.createSystem(modeller.topology)
    print("✓ System created successfully")
    
    # Assign each force to a different group for energy decomposition
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
except Exception as e:
    print(f"✗ System creation failed: {e}")
    print("This means some atoms don't have forcefield parameters")
    
    # Try to identify which atoms are problematic
    print("\nAtom types that might be missing:")
    atom_types = set()
    for atom in modeller.topology.atoms():
        atom_types.add(f"{atom.name}({atom.element.symbol}) in {atom.residue.name}")
    for atom_type in sorted(atom_types)[:10]:
        print(f"  - {atom_type}")
    sys.exit(1)

print("\n=== FORCE ANALYSIS ===")
total_forces = 0
for i, force in enumerate(system.getForces()):
    force_name = force.__class__.__name__
    print(f"\nForce {i}: {force_name}")
    
    # Check if force has any terms
    if hasattr(force, 'getNumBonds') and force.getNumBonds() > 0:
        print(f"  ✓ Bonds: {force.getNumBonds()}")
        # Show first bond as example
        bond = force.getBondParameters(0)
        print(f"    Example: atoms {bond[0]}-{bond[1]}, k={bond[2]}, r0={bond[3]}")
        total_forces += 1
    elif hasattr(force, 'getNumAngles') and force.getNumAngles() > 0:
        print(f"  ✓ Angles: {force.getNumAngles()}")
        angle = force.getAngleParameters(0)
        print(f"    Example: atoms {angle[0]}-{angle[1]}-{angle[2]}, k={angle[3]}, theta0={angle[4]}")
        total_forces += 1
    elif hasattr(force, 'getNumTorsions') and force.getNumTorsions() > 0:
        print(f"  ✓ Torsions: {force.getNumTorsions()}")
        torsion = force.getTorsionParameters(0)
        print(f"    Example: atoms {torsion[0]}-{torsion[1]}-{torsion[2]}-{torsion[3]}")
        total_forces += 1
    elif hasattr(force, 'getNumParticles') and force.getNumParticles() > 0:
        print(f"  ✓ Particles: {force.getNumParticles()}")
        # Check if particles have non-zero parameters
        particle = force.getParticleParameters(0)
        if len(particle) >= 3:  # charge, sigma, epsilon
            charge, sigma, epsilon = particle[:3]
            print(f"    Example: q={charge}, σ={sigma}, ε={epsilon}")
            if abs(charge.value_in_unit(elementary_charge)) > 1e-6 or abs(epsilon.value_in_unit(kilojoule_per_mole)) > 1e-6:
                total_forces += 1
                print("    ✓ Non-zero parameters found")
            else:
                print("    ⚠ All parameters are zero!")
    else:
        print(f"  ⚠ No terms found")

print(f"\nSummary: {total_forces} forces have active terms")

if total_forces == 0:
    print("❌ CRITICAL: No forces have active terms!")
    print("This explains why most energies are zero.")
    sys.exit(1)

# Now compute energies with detailed analysis
print("\n=== ENERGY COMPUTATION ===")
platform = Platform.getPlatformByName("CPU")
integrator = VerletIntegrator(0.001)
context = Context(system, integrator, platform)
context.setPositions(modeller.positions)



print("Individual energy terms (kJ/mol)")
print("-" * 50)

# Use force groups to get individual energies  
force_energies = {}
for i, force in enumerate(system.getForces()):
    force_name = force.__class__.__name__
    
    # Skip CMMotionRemover - it never contributes energy
    if force_name == "CMMotionRemover":
        continue
        
    # Use 2**i as the group mask to select only this force group
    state = context.getState(getEnergy=True, groups=2**i)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    force_energies[force_name] = energy
    
    # Add diagnostic info
    status = ""
    if abs(energy) < 1e-6:
        status = " ⚠ ZERO"
    elif energy > 1000:
        status = " ⚠ HIGH"
    elif energy < -1000:
        status = " ⚠ LOW"
    
    print(f"{force_name:<25} {energy:12.3f}{status}")

# Compare with total system energy
state = context.getState(getEnergy=True)
openmm_total = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
sum_individual = builtins.sum(force_energies.values())

print("-" * 50)
print(f"{'Sum of individual':<25} {sum_individual:12.3f}")
print(f"{'Total system energy':<25} {openmm_total:12.3f}")
print("-" * 50)
print(f"\nIndividual force energies sum: {sum_individual:.3f} kJ/mol")
print(f"Total system energy: {openmm_total:.3f} kJ/mol")
if abs(sum_individual - openmm_total) > 0.001:
    print(f"⚠ Difference: {abs(sum_individual - openmm_total):.3f} kJ/mol")
else:
    print("✓ Energies match - OpenMM working correctly!")

# Final diagnosis
print("\n=== DIAGNOSIS ===")
if abs(openmm_total) < 1e-3:
    print("❌ Problem: Total energy is essentially zero!")
    print("Possible causes:")
    print("- Atoms are too far apart (no interactions)")
    print("- Forcefield parameters missing or zero")
    print("- Structure has topology issues")
    
    # Additional checks
    print("\nAdditional checks:")
    coords = context.getState(getPositions=True).getPositions(asNumpy=True)
    distances = []
    for i in range(min(10, len(coords))):
        for j in range(i+1, min(10, len(coords))):
            dist = ((coords[i] - coords[j])**2).sum()**0.5
            distances.append(float(dist))
    
    if distances:
        min_dist = min(distances)
        max_dist = max(distances)
        print(f"  Distance range: {min_dist:.3f} - {max_dist:.3f} nm")
        if min_dist > 2.0:
            print("  ⚠ Atoms are very far apart - no interactions expected")
        elif min_dist < 0.05:
            print("  ⚠ Atoms are very close - might cause issues")
    
elif openmm_total > 10000:
    print("⚠ Problem: Very high energy - structure may have clashes")
elif openmm_total < -10000:
    print("⚠ Problem: Very negative energy - check parameters")
else:
    print("✅ Energy seems reasonable")

# Show how to fix common issues
print("\n=== SUGGESTIONS ===")
print("To debug further:")
print("1. Check your PDB residue names match AMBER14 expectations")
print("2. Verify all atoms have proper connectivity")  
print("3. Try a simpler forcefield like 'amber99sb.xml' first")
print("4. Use a known working structure to test your setup")