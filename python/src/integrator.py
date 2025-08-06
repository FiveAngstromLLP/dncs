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
import logging
import datetime
import concurrent.futures
from openmm.app import ForceField, PDBFile,Simulation,Modeller, StateDataReporter
from openmm.openmm import Platform, LangevinMiddleIntegrator
from openmm.unit import kelvin, nano, pico


class DncsIntegrator:
    def __init__(self, config):
        self.config = config
        self.forcefield = ForceField(*config.forcefield)
        self.inpfolder = f"{config.folder}/{config.moleculename}/sample"
        self.outfolder = f"{config.folder}/{config.moleculename}"
        os.makedirs(self.outfolder, exist_ok=True)
        allpdb = sorted([f for f in os.listdir(self.inpfolder) if f.endswith('.pdb')])
        self.pdbs = allpdb[0:self.config.md_simulation]
        self.log = self.setup_logger()

    def setup_logger(self):
        logging.basicConfig(
            filename=f"{self.outfolder}/dncs.log",
            filemode="w",
            level=logging.INFO,
            format='%(message)s'
        )
        return logging.getLogger(f"dncs_{self.config.moleculename}")

    def run_integrator(self):
        self.log_parameters()
        # Run sampling and minimization in parallel as before
        with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            futures = []
            for i, pdb in enumerate(self.pdbs):
                model = PDBFile(os.path.join(self.inpfolder, pdb))
                modeller = Modeller(model.topology, model.positions)
                futures.append(executor.submit(self.run_minimization, modeller, i+1))
            concurrent.futures.wait(futures)
        for future in futures:
            if future.exception():
                print(f"Error occurred: {future.exception()}")

        # Run equilibration with single context
        self.run_equilibration()

    def log_parameters(self):
        date = datetime.datetime.now()
        message = (f"{date}\nNo of Samples: {len(self.pdbs)}\nForceField: {self.config.forcefield}\n"
                   f"Integrator Parameters:\nSteps: {self.config.steps}\ntemperature: {self.config.temp} Kelvin\n"
                   f"dt: {self.config.dt} picoseconds\nNo of Solvent: {self.config.solvent}\n"
                   f"FrictionalCoefficient: {self.config.gamma} picosecond^(-1)")
        self.log.info(message)

    def run_minimization(self, modeller: Modeller, i: int):
        try:
            platform = Platform.getPlatformByName(self.config.device)

            modeller.addHydrogens()
            modeller.addSolvent(self.forcefield, numAdded=self.config.solvent)

            system = self.forcefield.createSystem(modeller.topology, ignoreExternalBonds=True)

            integrator = LangevinMiddleIntegrator(self.config.temp * kelvin, self.config.gamma / pico.factor, self.config.dt * pico.factor)

            simulation = Simulation(modeller.topology, system, integrator, platform)
            simulation.context.setPositions(modeller.getPositions())

            self.run_and_save_minimization(simulation, i)

        except Exception as e:
            print(f"Error in run_simulation for model {i}: {e}")

    def run_and_save_minimization(self, simulation: Simulation, i: int):
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        self.log.info(f"ENERGY FOR MODEL {i} = {state.getPotentialEnergy()}")
        print(f"Initial energy for model {i}: {state.getPotentialEnergy()}")

        os.makedirs(f"{self.outfolder}/Minimized", exist_ok=True)
        simulation.minimizeEnergy()
        minimized_state = simulation.context.getState(getEnergy=True, getPositions=True)
        self.log.info(f"MINIMIZED ENERGY FOR MODEL {i} = {minimized_state.getPotentialEnergy()}")
        print(f"Minimized energy for model {i}: {minimized_state.getPotentialEnergy()}")
        # Save minimized structure
        self.save_pdb(
            f"{self.outfolder}/Minimized/Minimized_{i:04}.pdb",
            simulation.topology,
            minimized_state.getPositions()
        )

    def run_equilibration(self):
        """Run equilibration with a single OpenMM context, cycling through minimized structures."""
        # Get all minimized structures
        minimized_dir = f"{self.outfolder}/Minimized"
        minimized_files = sorted([f for f in os.listdir(minimized_dir) if f.endswith('.pdb')])

        if not minimized_files:
            print("No minimized structures found. Skipping equilibration.")
            return

        print(f"Running equilibration with {len(minimized_files)} structures, {self.config.steps} steps each, reusing a single OpenMM context.")

        # Create Langevin directory
        os.makedirs(f"{self.outfolder}/Langevin", exist_ok=True)

        platform = Platform.getPlatformByName(self.config.device)

        # Load the topology and initial positions from the first PDB to set up the single OpenMM context.
        # This assumes all subsequent PDBs in minimized_files have an identical topology.
        first_pdb_path = os.path.join(minimized_dir, minimized_files[0])
        pdbdata_initial = PDBFile(first_pdb_path)

        system = self.forcefield.createSystem(pdbdata_initial.topology, ignoreExternalBonds=True)
        integrator = LangevinMiddleIntegrator(
            self.config.temp * kelvin,
            self.config.gamma / pico.factor,
            self.config.dt * pico.factor
        )

        simulation = Simulation(pdbdata_initial.topology, system, integrator, platform)

        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True,
                                            potentialEnergy=True,
                                            kineticEnergy=True,
                                            temperature=True))

        for i, pdb_file in enumerate(minimized_files):
            model_num = i + 1
            print(f"Equilibrating model {model_num}/{len(minimized_files)}")

            # Load positions for the current minimized structure.
            # The topology of this PDB is expected to be identical to the one used to initialize 'simulation'.
            current_pdb_data = PDBFile(os.path.join(minimized_dir, pdb_file))

            # Set positions of the existing OpenMM context to the current structure's positions
            simulation.context.setPositions(current_pdb_data.positions)

            # Re-initialize velocities to the target temperature for each new equilibration run segment
            temperature = self.config.temp * kelvin
            simulation.context.setVelocitiesToTemperature(temperature)

            # Run equilibration steps
            simulation.step(self.config.steps)

            # Get equilibrated state and save
            equilibrated_state = simulation.context.getState(getEnergy=True, getPositions=True)
            self.log.info(f"EQUILIBRATED ENERGY AFTER {self.config.steps} STEPS FOR MODEL {model_num} = {equilibrated_state.getPotentialEnergy()}")
            print(f"Equilibrated energy for model {model_num}: {equilibrated_state.getPotentialEnergy()}")

            # Save equilibrated structure. The topology used for saving is from the 'simulation' object,
            # which was set up with the first PDB's topology. This is consistent if all topologies are the same.
            self.save_pdb(
                f"{self.outfolder}/Langevin/Equilibrated_{model_num:04}.pdb",
                simulation.topology,
                equilibrated_state.getPositions()
            )

    @staticmethod
    def save_pdb(filename: str, topology, positions):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        PDBFile.writeFile(topology, positions, open(filename, "w"))


class CleanUp:
    def __init__(self, config):
        self.config = config
        self.inpfolder = f"{config.folder}/{self.config.moleculename}"
        self.process_files()
        self.write_equilibrated()
        self.write_minimized()

    def process_files(self):
        self.process_directory(f"{self.inpfolder}/Langevin", "equilibrated.out")
        self.process_directory(f"{self.inpfolder}/Minimized", "minimized.out")


    def process_directory(self, directory: str, output_file: str):
        pdb_files = sorted([f for f in os.listdir(directory) if f.endswith('.pdb')])
        with open(os.path.join(directory, output_file), 'w') as src:
            for file in pdb_files:
                i = int(file.split(".")[0].split("_")[1])
                pdb = PDBFile(os.path.join(directory, file))
                model = Modeller(pdb.topology, pdb.positions)
                self.write_model_info(src, i, model)
                del pdb

    def write_model_info(self, file, model_num: int, model: Modeller):
        pattern = self.get_energy_pattern(file.name)
        energy = self.find_energy(model_num, pattern)
        fname = ""
        if "equilibrated.out" in file.name:
            fname = f"Langevin/Equilibrated_{model_num:04}.pdb"
        elif "minimized.out" in file.name:
            fname = f"Minimized/Minimized_{model_num:04}.pdb"
        pdbfile = f"{self.inpfolder}/{fname}"
        angles_str = dncs.pdb_to_angle(pdbfile)
        file.write(f"{model_num}, {energy}, {angles_str}\n")

    def get_energy_pattern(self, filename: str) -> str:
        if "equilibrated.out" in filename:
            return r"EQUILIBRATED ENERGY AFTER \d+ STEPS FOR MODEL {model_num} = (.*)"
        elif "minimized.out" in filename:
            return r"MINIMIZED ENERGY FOR MODEL {model_num} = (.*)"
        else:
            raise ValueError(f"Unknown file type: {filename}")

    def find_energy(self, model_num: int, pattern: str) -> str:
        with open(f"{self.inpfolder}/dncs.log", "r") as src:
            data = src.read()
        match = re.search(pattern.format(model_num=model_num), data)
        return f"{match.group().split('=')[-1].strip()}" if match else "N/A"

    def write_equilibrated(self):
        with open(f"{self.inpfolder}/Langevin/equilibrated.out", "r") as f:
            minimized_lines = f.readlines()
        weng = []
        for line in minimized_lines:
            data = line.split(",")
            weng.append((data[0], data[1]))

        with open(f"{self.inpfolder}/equilibrated.pdb", "a") as file:
            for i,(m,e) in enumerate(sorted(weng, key=lambda x: x[1])):
                f = f"{self.inpfolder}/Langevin/Equilibrated_{int(m):04}.pdb"
                pdb = PDBFile(f)
                PDBFile.writeModel(pdb.topology, pdb.positions , file, modelIndex=i+1)


    def write_minimized(self):
        with open(f"{self.inpfolder}/Minimized/minimized.out", "r") as f:
            minimized_lines = f.readlines()
        weng = []
        for line in minimized_lines:
            data = line.split(",")
            weng.append((data[0], data[1]))

        with open(f"{self.inpfolder}/minimized.pdb", "a") as file:
            for i,(m,e) in enumerate(sorted(weng, key=lambda x: x[1])):
                f = f"{self.inpfolder}/Minimized/Minimized_{int(m):04}.pdb"
                pdb = PDBFile(f)
                PDBFile.writeModel(pdb.topology, pdb.positions , file, modelIndex=i+1)


class MDSimulation:
    def __init__(self, config):
        self.config = config
        self.forcefield = ForceField(*config.forcefield)
        self.inpfolder = f"{config.folder}/{self.config.moleculename}/Langevin"
        self.outfolder = f"{config.folder}/{self.config.moleculename}/MDSimulation"
        os.makedirs(self.outfolder, exist_ok=True)
        self.pdbs = sorted([f for f in os.listdir(self.inpfolder) if f.endswith('.pdb')])


    def run_simulation(self):
        """Run production MD with a single OpenMM context, cycling through equilibrated structures."""
        if not self.pdbs:
            print("No equilibrated structures found. Skipping production MD.")
            return

        print(f"Running production MD with {len(self.pdbs)} structures, {self.config.md_steps} steps each, reusing a single OpenMM context.")

        platform = Platform.getPlatformByName(self.config.device)

        # Load the topology and initial positions from the first PDB to set up the single OpenMM context.
        # This assumes all subsequent PDBs in self.pdbs have an identical topology.
        first_pdb_path = os.path.join(self.inpfolder, self.pdbs[0])
        pdbdata_initial = PDBFile(first_pdb_path)

        system = self.forcefield.createSystem(pdbdata_initial.topology)
        integrator = LangevinMiddleIntegrator(
            self.config.temp * kelvin,
            self.config.gamma / pico.factor,
            self.config.dt * pico.factor
        )

        simulation = Simulation(pdbdata_initial.topology, system, integrator, platform)

        # Add reporters once to the single simulation object
        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True,
                                            potentialEnergy=True,
                                            kineticEnergy=True,
                                            temperature=True))

        for i, pdb_file in enumerate(self.pdbs):
            print(f"Production MDSimulation {i+1}/{len(self.pdbs)}")

            # Load positions for the current equilibrated structure.
            # The topology of this PDB is expected to be identical to the one used to initialize 'simulation'.
            current_pdb_data = PDBFile(os.path.join(self.inpfolder, pdb_file))

            # Set positions of the existing OpenMM context to the current structure's positions
            simulation.context.setPositions(current_pdb_data.positions)

            # Re-initialize velocities to the target temperature for each new production run segment
            temperature = self.config.temp * kelvin
            simulation.context.setVelocitiesToTemperature(temperature)

            # Run the production MD steps for this structure
            simulation.step(self.config.md_steps)

            # Get the final state (positions) from the simulation context
            final_state = simulation.context.getState(getPositions=True)

            # Save the final structure. The topology used for saving is from the 'simulation' object,
            # which was set up with the first PDB's topology. This is consistent if all topologies are the same.
            save_pdb(f"{self.outfolder}/md_{i+1:04}.pdb", simulation.topology, final_state.getPositions())

def save_pdb(filename: str, topology, positions):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    PDBFile.writeFile(topology, positions, open(filename, "w"))
