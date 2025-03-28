# Digital Nets Conformational Sampling (DNCS)
# Copyright [2024] [Abraham Rebairo J., Satheeshkumar S, Sam Paul D., Stephen A. and FiveAngstromLLP]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
import re
import dncs
import math
import logging
import datetime
import multiprocessing
from openmm.app import ForceField, PDBFile, Simulation, Modeller, StateDataReporter
from openmm.openmm import Platform, LangevinMiddleIntegrator
from openmm.unit import kelvin, nano, pico

class OptimizedMolecularDynamics:
    def __init__(self, config):
        self.config = config
        self.forcefield = ForceField(*config.forcefield)
        self.inpfolder = f"{config.folder}/{config.moleculename}/sample"
        self.outfolder = f"{config.folder}/{config.moleculename}"
        os.makedirs(self.outfolder, exist_ok=True)
        os.makedirs(f"{self.outfolder}/Langevin", exist_ok=True)
        os.makedirs(f"{self.outfolder}/Minimized", exist_ok=True)

        # Get PDB files
        allpdb = sorted([f for f in os.listdir(self.inpfolder) if f.endswith('.pdb')])
        self.pdbs = allpdb[:self.config.md_simulation]

        # Setup logging
        self.log = self._setup_logger()

        # Results storage
        self.results = []

    def _setup_logger(self):
        logging.basicConfig(
            filename=f"{self.outfolder}/dncs.log",
            filemode="w",
            level=logging.INFO,
            format='%(message)s'
        )
        return logging.getLogger(f"dncs_{self.config.moleculename}")

    def _log_parameters(self):
        date = datetime.datetime.now()
        message = (f"{date}\nNo of Samples: {len(self.pdbs)}\nForceField: {self.config.forcefield}\n"
                   f"Integrator Parameters:\nSteps: {self.config.steps}\ntemperature: {self.config.temp} Kelvin\n"
                   f"dt: {self.config.dt} picoseconds\nNo of Solvent: {self.config.solvent}\n"
                   f"FrictionalCoefficient: {self.config.gamma} picosecond^(-1)")
        self.log.info(message)

    def process_pdb(self, pdb_path, index):
        """
        Process a single PDB file for molecular dynamics simulation

        Args:
            pdb_path (str): Path to the PDB file
            index (int): Model index

        Returns:
            tuple: Simulation results and metadata
        """
        try:
            # Load PDB
            model = PDBFile(pdb_path)
            modeller = Modeller(model.topology, model.positions)

            # Platform selection
            platform = Platform.getPlatformByName(self.config.device)

            # Prepare system
            modeller.addHydrogens()
            modeller.addSolvent(self.forcefield, padding=1.0 * nano.factor)
            modeller.addSolvent(self.forcefield, numAdded=self.config.solvent)

            # Create system and integrator
            system = self.forcefield.createSystem(modeller.topology, ignoreExternalBonds=True)
            integrator = LangevinMiddleIntegrator(
                self.config.temp * kelvin,
                self.config.gamma / pico.factor,
                self.config.dt * pico.factor
            )

            # Run simulation
            simulation = Simulation(modeller.topology, system, integrator, platform)
            simulation.context.setPositions(modeller.getPositions())

            # Get initial state
            initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
            initial_energy = float(str(initial_state.getPotentialEnergy()).split(" ")[0])

            # Minimize energy
            simulation.minimizeEnergy()
            minimized_state = simulation.context.getState(getEnergy=True, getPositions=True)
            min_energy = float(str(minimized_state.getPotentialEnergy()).split(" ")[0])

            # Run equilibration
            simulation.step(self.config.steps)
            equilibrated_state = simulation.context.getState(getEnergy=True, getPositions=True)
            equil_energy = float(str(equilibrated_state.getPotentialEnergy()).split(" ")[0])

            # Calculate weight
            weight = math.exp(-equil_energy / (self.config.temp * 1.380649e-23 * 6.02214076e23))

            # Paths for saving
            equil_pdb_path = f"{self.outfolder}/Langevin/Equilibrated_{index:04}.pdb"
            min_pdb_path = f"{self.outfolder}/Minimized/Minimized_{index:04}.pdb"

            # Save equilibrated PDB
            PDBFile.writeFile(
                simulation.topology,
                equilibrated_state.getPositions(),
                open(equil_pdb_path, "w")
            )

            # Save minimized PDB if weight condition met
            min_pdb_path = None
            if weight > 1.0:
                min_pdb_path = f"{self.outfolder}/Minimized/Minimized_{index:04}.pdb"
                PDBFile.writeFile(
                    simulation.topology,
                    minimized_state.getPositions(),
                    open(min_pdb_path, "w")
                )

            # Calculate angles
            angles_str = dncs.pdb_to_angle(equil_pdb_path)

            # Clean up large objects
            del simulation, system, integrator

            return (
                index,
                initial_energy,
                min_energy,
                equil_energy,
                min_pdb_path,
                equil_pdb_path,
                angles_str
            )

        except Exception as e:
            print(f"Error processing {pdb_path}: {e}")
            return None

    def run_parallel_simulation(self):
        """
        Run molecular dynamics simulations in parallel
        """
        # Log parameters
        self._log_parameters()

        # Prepare input paths
        input_paths = [os.path.join(self.inpfolder, pdb) for pdb in self.pdbs]

        # Use multiprocessing to run simulations
        with multiprocessing.Pool(processes=max(1, multiprocessing.cpu_count() - 1)) as pool:
            # Use starmap to pass both path and index
            results = pool.starmap(
                self.process_pdb,
                [(path, i+1) for i, path in enumerate(input_paths)]
            )

        # Filter out None results (failed runs)
        self.results = [r for r in results if r is not None]

        # Sort results by equilibration energy
        self.results.sort(key=lambda x: x[3])

        # Write output files
        self._write_output_files()

    def _write_output_files(self):
        """
        Write combined equilibrated and minimized PDB files
        """
        # Write equilibrated PDB
        with open(f"{self.outfolder}/equilibrated.pdb", "w") as equil_file:
            for i, result in enumerate(self.results, 1):
                if result[5]:  # equilibrated PDB path
                    pdb = PDBFile(result[5])
                    PDBFile.writeModel(pdb.topology, pdb.positions, equil_file, modelIndex=i)

        # Write minimized PDB
        with open(f"{self.outfolder}/minimized.pdb", "w") as min_file:
            for i, result in enumerate(self.results, 1):
                if result[4]:  # minimized PDB path
                    pdb = PDBFile(result[4])
                    PDBFile.writeModel(pdb.topology, pdb.positions, min_file, modelIndex=i)

        # Write energy and angle data
        with open(f"{self.outfolder}/equilibrated.out", "w") as equil_out, \
             open(f"{self.outfolder}/minimized.out", "w") as min_out:
            for result in self.results:
                index, _, min_energy, equil_energy, _, _, angles_str = result
                equil_out.write(f"{index}, {equil_energy}, {angles_str}\n")
                if result[4]:  # Only write minimized if path exists
                    min_out.write(f"{index}, {min_energy}, {angles_str}\n")

class MDSimulation:
    def __init__(self, config):
        self.config = config
        self.forcefield = ForceField(*config.forcefield)
        self.inpfolder = f"{config.folder}/{self.config.moleculename}/Langevin"
        self.outfolder = f"{config.folder}/{self.config.moleculename}/MDSimulation"
        os.makedirs(self.outfolder, exist_ok=True)
        self.pdbs = sorted([f for f in os.listdir(self.inpfolder) if f.endswith('.pdb')])

    def run_simulation(self):
        platform = Platform.getPlatformByName(self.config.device)
        steps_per_segment = int(self.config.md_steps / len(self.pdbs))

        for i, pdb in enumerate(self.pdbs):
            print(f"Processing {i}th structure..")
            pdbdata = PDBFile(f"{self.inpfolder}/{pdb}")

            # Create system with error handling
            try:
                system = self.forcefield.createSystem(pdbdata.topology)
                integrator = LangevinMiddleIntegrator(
                    self.config.temp * kelvin,
                    self.config.gamma/pico.factor,
                    self.config.dt*pico.factor
                )

                simulation = Simulation(pdbdata.topology, system, integrator, platform)
                simulation.context.setPositions(pdbdata.positions)

                # Optional: Add reporters if detailed tracking is needed
                simulation.reporters.append(StateDataReporter(
                    sys.stdout, 100, step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    temperature=True
                ))

                simulation.step(steps_per_segment)

                # Save final positions
                position = simulation.context.getState(getPositions=True).getPositions()
                self._save_pdb(f"{self.outfolder}/simulated_{i}.pdb", simulation.topology, position)

                # Clean up resources
                del simulation, system, integrator

            except Exception as e:
                print(f"Error in MD simulation for {pdb}: {e}")

    def _save_pdb(self, filename: str, topology, positions):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        PDBFile.writeFile(topology, positions, open(filename, "w"))
