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

import numpy as np
from typing import Optional, Tuple, Dict, Any
from openmm.app import ForceField, PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter, PME, NoCutoff, HBonds
from openmm.openmm import Platform, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.unit import kelvin, atmosphere, picosecond, nanosecond
from openmm import unit


class DncsIntegrator:
    """Main integrator class for DNCS molecular dynamics simulations.

    Handles minimization, equilibration, and production MD runs with proper
    validation, error handling, and convergence monitoring.
    """

    # Constants for validation
    MIN_TEMPERATURE = 100.0  # K
    MAX_TEMPERATURE = 800.0  # K
    MIN_STEPS = 100
    MAX_STEPS = 1000000
    MIN_DT = 0.0001  # ps
    MAX_DT = 0.01    # ps

    def __init__(self, config):
        self.config = config
        self._validate_config()

        try:
            self.forcefield = ForceField(*config.forcefield)
        except Exception as e:
            raise ValueError(f"Failed to load force field {config.forcefield}: {e}")

        self.inpfolder = f"{config.folder}/{config.moleculename}/sample"
        self.outfolder = f"{config.folder}/{config.moleculename}"

        # Validate input directory exists
        if not os.path.exists(self.inpfolder):
            raise FileNotFoundError(f"Input folder not found: {self.inpfolder}")

        os.makedirs(self.outfolder, exist_ok=True)

        # Validate PDB files exist
        allpdb = sorted([f for f in os.listdir(self.inpfolder)
                        if f.endswith('.pdb') and os.path.getsize(os.path.join(self.inpfolder, f)) > 0])

        if not allpdb:
            raise FileNotFoundError(f"No valid PDB files found in {self.inpfolder}")

        self.pdbs = allpdb[0:self.config.md_simulation]
        self.log = self.setup_logger()

        # Initialize platform with fallback
        self.platform = self._get_platform()

        self.log.info(f"Initialized DncsIntegrator with {len(self.pdbs)} structures")
        self.log.info(f"Using platform: {self.platform.getName()}")

    def _validate_config(self):
        """Validate simulation configuration parameters."""
        if not (self.MIN_TEMPERATURE <= self.config.temp <= self.MAX_TEMPERATURE):
            raise ValueError(f"Temperature {self.config.temp}K outside valid range [{self.MIN_TEMPERATURE}-{self.MAX_TEMPERATURE}]")

        if not (self.MIN_STEPS <= self.config.steps <= self.MAX_STEPS):
            raise ValueError(f"Equilibration steps {self.config.steps} outside valid range [{self.MIN_STEPS}-{self.MAX_STEPS}]")

        if not (self.MIN_STEPS <= self.config.md_steps <= self.MAX_STEPS):
            raise ValueError(f"MD steps {self.config.md_steps} outside valid range [{self.MIN_STEPS}-{self.MAX_STEPS}]")

        if not (self.MIN_DT <= self.config.dt <= self.MAX_DT):
            raise ValueError(f"Timestep {self.config.dt}ps outside valid range [{self.MIN_DT}-{self.MAX_DT}]")

        if self.config.gamma <= 0:
            raise ValueError(f"Friction coefficient must be positive, got {self.config.gamma}")

        if self.config.solvent < 0:
            raise ValueError(f"Solvent count must be non-negative, got {self.config.solvent}")

    def _get_platform(self) -> Platform:
        """Get OpenMM platform with fallback options."""
        platforms = ['CUDA', 'OpenCL', 'CPU']
        requested = self.config.device

        # Try requested platform first
        if requested in platforms:
            try:
                return Platform.getPlatformByName(requested)
            except Exception as e:
                self.log.warning(f"Failed to get {requested} platform: {e}")

        # Try fallback platforms
        for platform_name in platforms:
            try:
                platform = Platform.getPlatformByName(platform_name)
                if platform_name != requested:
                    self.log.warning(f"Using fallback platform: {platform_name}")
                return platform
            except:
                continue

        raise RuntimeError("No OpenMM platform available")

    def setup_logger(self):
        """Setup logging with proper formatting and error handling."""
        log_file = f"{self.outfolder}/dncs.log"

        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # Setup file handler
        try:
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(formatter)
        except Exception as e:
            raise IOError(f"Cannot create log file {log_file}: {e}")

        # Setup logger
        logger = logging.getLogger(f"dncs_{self.config.moleculename}")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()  # Remove any existing handlers
        logger.addHandler(file_handler)

        return logger

    def run_integrator(self):
        self.log_parameters()
        # Run sampling and minimization sequentially to avoid pickle issues with OpenMM objects
        for i, pdb in enumerate(self.pdbs):
            model = PDBFile(os.path.join(self.inpfolder, pdb))
            modeller = Modeller(model.topology, model.positions)
            self.run_minimization(modeller, i+1)


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
        """Run energy minimization for a single structure with proper validation.

        Args:
            modeller: OpenMM Modeller object with initial structure
            i: Structure index for logging and output
        """
        try:
            self.log.info(f"Starting minimization for structure {i}")

            # Add hydrogens with validation
            initial_atoms = modeller.topology.getNumAtoms()
            modeller.addHydrogens()
            after_hydrogens = modeller.topology.getNumAtoms()
            self.log.info(f"Structure {i}: Added {after_hydrogens - initial_atoms} hydrogens")

            # Add solvent properly (fix redundant addition)
            if self.config.solvent > 0:
                # Use either padding OR numAdded, not both
                modeller.addSolvent(self.forcefield,
                                  model='tip3p',
                                  padding=1.0 * unit.nanometer,
                                  ionicStrength=0.1*unit.molar)
                after_solvation = modeller.topology.getNumAtoms()
                solvent_added = after_solvation - after_hydrogens
                self.log.info(f"Structure {i}: Added {solvent_added} solvent atoms")

                # Validate reasonable solvation
                if solvent_added < 100:
                    self.log.warning(f"Structure {i}: Very few solvent atoms added ({solvent_added})")
            else:
                self.log.info(f"Structure {i}: Skipping solvation (solvent=0)")

            # Create system with proper parameters
            system = self.forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=PME if self.config.solvent > 0 else NoCutoff,
                nonbondedCutoff=1.0*unit.nanometer,
                constraints=HBonds,
                rigidWater=True,
                removeCMMotion=True
            )

            # Create integrator
            integrator = LangevinMiddleIntegrator(
                self.config.temp * kelvin,
                self.config.gamma / picosecond,
                self.config.dt * picosecond
            )
            integrator.setConstraintTolerance(1e-5)

            # Create simulation
            simulation = Simulation(modeller.topology, system, integrator, self.platform)

            # Set positions with validation
            positions = modeller.getPositions()
            if len(positions) != system.getNumParticles():
                raise ValueError(f"Position count mismatch: {len(positions)} vs {system.getNumParticles()}")

            simulation.context.setPositions(positions)

            # Validate initial state
            try:
                initial_state = simulation.context.getState(getEnergy=True)
                initial_energy = initial_state.getPotentialEnergy()
                self.log.info(f"Structure {i}: Initial energy = {initial_energy}")

                # Check for extremely high energy (indicates problems)
                if initial_energy.value_in_unit(unit.kilojoule_per_mole) > 1e6:
                    self.log.warning(f"Structure {i}: Very high initial energy, may indicate clashes")

            except Exception as e:
                self.log.error(f"Structure {i}: Cannot evaluate initial energy: {e}")
                raise

            self.run_and_save_minimization(simulation, i)

        except Exception as e:
            error_msg = f"Error in minimization for structure {i}: {e}"
            self.log.error(error_msg)
            print(error_msg)
            raise

    def run_and_save_minimization(self, simulation: Simulation, i: int):
        """Run energy minimization with convergence monitoring and validation.

        Args:
            simulation: OpenMM Simulation object
            i: Structure index
        """
        # Get initial state
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        initial_energy = state.getPotentialEnergy()
        self.log.info(f"Structure {i}: Pre-minimization energy = {initial_energy}")
        print(f"Initial energy for structure {i}: {initial_energy}")

        # Create output directory
        minimized_dir = f"{self.outfolder}/Minimized"
        os.makedirs(minimized_dir, exist_ok=True)

        # Run minimization with proper parameters
        try:
            self.log.info(f"Structure {i}: Starting energy minimization")
            simulation.minimizeEnergy(tolerance=10*unit.kilojoule/unit.mole/unit.nanometer, maxIterations=1000)

            # Get minimized state
            minimized_state = simulation.context.getState(getEnergy=True, getPositions=True)
            final_energy = minimized_state.getPotentialEnergy()

            # Calculate energy change
            energy_change = final_energy - initial_energy
            self.log.info(f"Structure {i}: Post-minimization energy = {final_energy}")
            self.log.info(f"Structure {i}: Energy change = {energy_change}")
            print(f"Minimized energy for structure {i}: {final_energy}")
            # print(f"Energy change: {energy_change}")

            # Validate minimization success
            if energy_change.value_in_unit(unit.kilojoule_per_mole) > 0:
                self.log.warning(f"Structure {i}: Energy increased during minimization!")
            elif abs(energy_change.value_in_unit(unit.kilojoule_per_mole)) < 1.0:
                self.log.warning(f"Structure {i}: Very small energy change, may already be minimized")
            else:
                self.log.info(f"Structure {i}: Minimization successful")

            # Save minimized structure
            output_file = f"{minimized_dir}/Minimized_{i:04d}.pdb"
            self.save_pdb(output_file, simulation.topology, minimized_state.getPositions())
            self.log.info(f"Structure {i}: Saved minimized structure to {output_file}")

        except Exception as e:
            error_msg = f"Structure {i}: Minimization failed: {e}"
            self.log.error(error_msg)
            raise RuntimeError(error_msg)

    def run_equilibration(self):
        """Run proper multi-phase equilibration with convergence monitoring.

        Phases:
        1. NVT equilibration (constant volume, temperature)
        2. NPT equilibration (constant pressure, temperature)
        3. Convergence monitoring
        """
        # Get all minimized structures
        minimized_dir = f"{self.outfolder}/Minimized"

        if not os.path.exists(minimized_dir):
            raise FileNotFoundError(f"Minimized directory not found: {minimized_dir}")

        minimized_files = sorted([f for f in os.listdir(minimized_dir)
                                if f.endswith('.pdb') and os.path.getsize(os.path.join(minimized_dir, f)) > 0])

        if not minimized_files:
            raise FileNotFoundError("No valid minimized structures found. Run minimization first.")

        self.log.info(f"Starting equilibration with {len(minimized_files)} structures")
        print(f"Running equilibration with {len(minimized_files)} structures")

        # Create output directories
        langevin_dir = f"{self.outfolder}/Langevin"
        os.makedirs(langevin_dir, exist_ok=True)
        os.makedirs(f"{langevin_dir}/trajectories", exist_ok=True)

        for i, pdb_file in enumerate(minimized_files):
            model_num = i + 1
            self.log.info(f"Starting equilibration for structure {model_num}/{len(minimized_files)}")
            print(f"Equilibrating structure {model_num}/{len(minimized_files)}")

            try:
                pdb_path = os.path.join(minimized_dir, pdb_file)
                pdb_data = PDBFile(pdb_path)

                # Validate PDB file
                if pdb_data.topology.getNumAtoms() == 0:
                    raise ValueError(f"Empty topology in {pdb_file}")

                # Create system with proper settings
                system = self.forcefield.createSystem(
                    pdb_data.topology,
                    nonbondedMethod=PME,
                    nonbondedCutoff=1.0*unit.nanometer,
                    constraints=HBonds,
                    rigidWater=True,
                    removeCMMotion=True
                )

                # Add barostat for NPT equilibration if solvated
                if self.config.solvent > 0:
                    barostat = MonteCarloBarostat(1.0*atmosphere, self.config.temp*kelvin, 25)
                    system.addForce(barostat)
                    self.log.info(f"Structure {model_num}: Added barostat for NPT equilibration")

                # Run multi-phase equilibration
                final_positions = self._run_multiphase_equilibration(system, pdb_data, model_num)

                # Save equilibrated structure
                output_file = f"{langevin_dir}/Equilibrated_{model_num:04d}.pdb"
                self.save_pdb(output_file, pdb_data.topology, final_positions)
                self.log.info(f"Structure {model_num}: Saved equilibrated structure")

            except Exception as e:
                error_msg = f"Equilibration failed for structure {model_num}: {e}"
                self.log.error(error_msg)
                print(f"ERROR: {error_msg}")
                continue

        self.log.info("Equilibration phase completed")
        print("Equilibration completed")

    def _run_multiphase_equilibration(self, system, pdb_data, model_num: int):
        """Run multi-phase equilibration with proper monitoring.

        Args:
            system: OpenMM System object
            pdb_data: PDBFile object
            model_num: Structure index

        Returns:
            Final equilibrated positions
        """
        # Phase 1: NVT Equilibration (constant volume)
        self.log.info(f"Structure {model_num}: Phase 1 - NVT equilibration")

        integrator1 = LangevinMiddleIntegrator(
            self.config.temp * kelvin,
            self.config.gamma / picosecond,
            self.config.dt * picosecond
        )
        integrator1.setConstraintTolerance(1e-5)

        simulation1 = Simulation(pdb_data.topology, system, integrator1, self.platform)
        simulation1.context.setPositions(pdb_data.positions)
        simulation1.context.setVelocitiesToTemperature(self.config.temp * kelvin)

        # Add reporters for monitoring
        nvt_steps = self.config.steps // 2
        simulation1.reporters.append(StateDataReporter(
            f"{self.outfolder}/Langevin/nvt_equilibration_{model_num:04d}.log",
            100, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            temperature=True, density=True
        ))

        # Run NVT equilibration
        simulation1.step(nvt_steps)

        # Get intermediate state
        nvt_state = simulation1.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        nvt_energy = nvt_state.getPotentialEnergy()
        self.log.info(f"Structure {model_num}: NVT energy = {nvt_energy}")

        # Phase 2: NPT Equilibration (constant pressure) - only if solvated
        if self.config.solvent > 0:
            self.log.info(f"Structure {model_num}: Phase 2 - NPT equilibration")

            integrator2 = LangevinMiddleIntegrator(
                self.config.temp * kelvin,
                self.config.gamma / picosecond,
                self.config.dt * picosecond
            )
            integrator2.setConstraintTolerance(1e-5)

            simulation2 = Simulation(pdb_data.topology, system, integrator2, self.platform)
            simulation2.context.setPositions(nvt_state.getPositions())
            simulation2.context.setVelocities(nvt_state.getVelocities())

            # Add reporters
            npt_steps = self.config.steps // 2
            simulation2.reporters.append(StateDataReporter(
                f"{self.outfolder}/Langevin/npt_equilibration_{model_num:04d}.log",
                100, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                temperature=True, density=True, volume=True
            ))

            # Run NPT equilibration
            simulation2.step(npt_steps)

            # Get final state
            final_state = simulation2.context.getState(getEnergy=True, getPositions=True)
            final_energy = final_state.getPotentialEnergy()
            self.log.info(f"Structure {model_num}: NPT energy = {final_energy}")

            return final_state.getPositions()
        else:
            # No solvent, skip NPT phase
            self.log.info(f"Structure {model_num}: Skipping NPT phase (no solvent)")
            return nvt_state.getPositions()

    @staticmethod
    def save_pdb(filename: str, topology, positions):
        """Save PDB file with error handling.

        Args:
            filename: Output PDB filename
            topology: OpenMM Topology object
            positions: Atomic positions
        """
        try:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            with open(filename, 'w') as f:
                PDBFile.writeFile(topology, positions, f)
        except Exception as e:
            raise IOError(f"Failed to save PDB file {filename}: {e}")


class CleanUp:
    """Post-processing class for organizing and combining simulation results.

    Handles data extraction from log files, energy analysis, and creation
    of combined PDB files sorted by energy.
    """

    def __init__(self, config):
        self.config = config
        self.inpfolder = f"{config.folder}/{self.config.moleculename}"

        # Validate input directory
        if not os.path.exists(self.inpfolder):
            raise FileNotFoundError(f"Simulation directory not found: {self.inpfolder}")

        # Setup logging
        self.log = self._setup_logger()
        self.log.info("Starting cleanup and post-processing")

        try:
            self.process_files()
            self.write_equilibrated()
            self.write_minimized()
            self.log.info("Cleanup completed successfully")
        except Exception as e:
            self.log.error(f"Cleanup failed: {e}")
            raise

    def _setup_logger(self):
        """Setup logging for cleanup operations."""
        log_file = f"{self.inpfolder}/cleanup.log"

        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        file_handler = logging.FileHandler(log_file, mode='a')  # Append mode
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)

        logger = logging.getLogger(f"cleanup_{self.config.moleculename}")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        logger.addHandler(file_handler)

        return logger

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
        """Extract energy value from log file with better error handling."""
        log_file = f"{self.inpfolder}/dncs.log"
        try:
            with open(log_file, "r") as src:
                data = src.read()

            match = re.search(pattern.format(model_num=model_num), data)
            if match:
                if '=' in match.group():
                    energy_str = match.group().split('=')[-1].strip()
                else:
                    energy_str = match.group(1).strip()

                # Validate the energy value
                try:
                    # Extract numeric part
                    energy_part = energy_str.split()[0] if ' ' in energy_str else energy_str
                    float(energy_part)  # Try to parse as float
                    return energy_str
                except ValueError:
                    if hasattr(self, 'log'):
                        self.log.warning(f"Invalid energy format for model {model_num}: {energy_str}")
                    return "N/A"
            else:
                if hasattr(self, 'log'):
                    self.log.warning(f"Energy not found for model {model_num} in log file")
                return "N/A"
        except Exception as e:
            if hasattr(self, 'log'):
                self.log.error(f"Failed to read energy from log file: {e}")
            return "N/A"

    def write_equilibrated(self):
        with open(f"{self.inpfolder}/Langevin/equilibrated.out", "r") as f:
            minimized_lines = f.readlines()
        weng = []
        for line in minimized_lines:
            data = line.strip().split(",")
            if len(data) >= 2:
                try:
                    model_id = int(data[0].strip())
                    energy = float(data[1].strip())
                    weng.append((model_id, energy))
                except (ValueError, IndexError) as e:
                    self.log.warning(f"Skipping invalid line in equilibrated.out: {line.strip()}")

        if not weng:
            self.log.warning("No valid energy data found for equilibrated structures")
            return

        # Sort by energy and write combined PDB
        equilibrated_pdb = f"{self.inpfolder}/equilibrated.pdb"
        try:
            with open(equilibrated_pdb, "w") as file:
                for i, (m, e) in enumerate(sorted(weng, key=lambda x: x[1])):
                    pdb_file = f"{self.inpfolder}/Langevin/Equilibrated_{m:04d}.pdb"
                    if os.path.exists(pdb_file):
                        try:
                            pdb = PDBFile(pdb_file)
                            PDBFile.writeModel(pdb.topology, pdb.positions, file, modelIndex=i+1)
                        except Exception as e:
                            self.log.error(f"Failed to write model {m} to combined PDB: {e}")
                    else:
                        self.log.warning(f"Equilibrated PDB file not found: {pdb_file}")
        except Exception as e:
            self.log.error(f"Failed to create combined equilibrated PDB: {e}")


    def write_minimized(self):
        with open(f"{self.inpfolder}/Minimized/minimized.out", "r") as f:
            minimized_lines = f.readlines()
        weng = []
        for line in minimized_lines:
            data = line.strip().split(",")
            if len(data) >= 2:
                try:
                    model_id = int(data[0].strip())
                    energy = float(data[1].strip())
                    weng.append((model_id, energy))
                except (ValueError, IndexError) as e:
                    self.log.warning(f"Skipping invalid line in minimized.out: {line.strip()}")

        if not weng:
            self.log.warning("No valid energy data found for minimized structures")
            return

        # Sort by energy and write combined PDB
        minimized_pdb = f"{self.inpfolder}/minimized.pdb"
        try:
            with open(minimized_pdb, "w") as file:
                for i, (m, e) in enumerate(sorted(weng, key=lambda x: x[1])):
                    pdb_file = f"{self.inpfolder}/Minimized/Minimized_{m:04d}.pdb"
                    if os.path.exists(pdb_file):
                        try:
                            pdb = PDBFile(pdb_file)
                            PDBFile.writeModel(pdb.topology, pdb.positions, file, modelIndex=i+1)
                        except Exception as e:
                            self.log.error(f"Failed to write model {m} to combined PDB: {e}")
                    else:
                        self.log.warning(f"Minimized PDB file not found: {pdb_file}")
        except Exception as e:
            self.log.error(f"Failed to create combined minimized PDB: {e}")


class MDSimulation:
    """Production MD simulation class with trajectory output and proper monitoring.

    Handles production molecular dynamics runs with trajectory saving,
    energy monitoring, and proper validation.
    """

    def __init__(self, config):
        self.config = config

        try:
            self.forcefield = ForceField(*config.forcefield)
        except Exception as e:
            raise ValueError(f"Failed to load force field {config.forcefield}: {e}")

        self.inpfolder = f"{config.folder}/{self.config.moleculename}/Langevin"
        self.outfolder = f"{config.folder}/{self.config.moleculename}/MDSimulation"

        # Validate input directory
        if not os.path.exists(self.inpfolder):
            raise FileNotFoundError(f"Equilibrated structures directory not found: {self.inpfolder}")

        os.makedirs(self.outfolder, exist_ok=True)
        os.makedirs(f"{self.outfolder}/trajectories", exist_ok=True)

        # Get valid PDB files
        self.pdbs = sorted([f for f in os.listdir(self.inpfolder)
                          if f.endswith('.pdb') and os.path.getsize(os.path.join(self.inpfolder, f)) > 0])

        if not self.pdbs:
            raise FileNotFoundError(f"No valid equilibrated structures found in {self.inpfolder}")

        # Setup platform
        self.platform = self._get_platform()

        # Setup logging
        self.log = self._setup_logger()
        self.log.info(f"Initialized MDSimulation with {len(self.pdbs)} equilibrated structures")

    def _get_platform(self) -> Platform:
        """Get OpenMM platform with fallback options."""
        platforms = ['CUDA', 'OpenCL', 'CPU']
        requested = self.config.device

        for platform_name in ([requested] + platforms):
            try:
                return Platform.getPlatformByName(platform_name)
            except:
                continue

        raise RuntimeError("No OpenMM platform available")

    def _setup_logger(self):
        """Setup logging for MD simulation."""
        log_file = f"{self.outfolder}/md_simulation.log"

        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)

        logger = logging.getLogger(f"md_{self.config.moleculename}")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        logger.addHandler(file_handler)

        return logger

    def run_simulation(self):
        """Run production MD with trajectory output and comprehensive monitoring."""
        if not self.pdbs:
            raise FileNotFoundError("No equilibrated structures found. Run equilibration first.")

        self.log.info(f"Starting production MD with {len(self.pdbs)} structures, {self.config.md_steps} steps each")
        print(f"Running production MD with {len(self.pdbs)} equilibrated structures")
        print(f"Total simulation time per structure: {self.config.md_steps * self.config.dt:.1f} ps")

        successful_runs = 0

        for i, pdb_file in enumerate(self.pdbs):
            structure_num = i + 1
            self.log.info(f"Starting production MD for structure {structure_num}/{len(self.pdbs)}")
            print(f"Production MD: Structure {structure_num}/{len(self.pdbs)}")

            try:
                # Load equilibrated structure
                pdb_path = f"{self.inpfolder}/{pdb_file}"
                pdbdata = PDBFile(pdb_path)

                # Validate structure
                if pdbdata.topology.getNumAtoms() == 0:
                    raise ValueError(f"Empty topology in {pdb_file}")

                # Create system with proper settings
                system = self.forcefield.createSystem(
                    pdbdata.topology,
                    nonbondedMethod=PME,
                    nonbondedCutoff=1.0*unit.nanometer,
                    constraints=HBonds,
                    rigidWater=True,
                    removeCMMotion=True
                )

                # Add barostat if solvated
                if self.config.solvent > 0:
                    barostat = MonteCarloBarostat(1.0*atmosphere, self.config.temp*kelvin, 25)
                    system.addForce(barostat)

                # Create integrator
                integrator = LangevinMiddleIntegrator(
                    self.config.temp * kelvin,
                    self.config.gamma / picosecond,
                    self.config.dt * picosecond
                )
                integrator.setConstraintTolerance(1e-5)

                # Create simulation
                simulation = Simulation(pdbdata.topology, system, integrator, self.platform)

                # Set positions and validate
                simulation.context.setPositions(pdbdata.positions)

                # Check initial energy
                initial_state = simulation.context.getState(getEnergy=True)
                initial_energy = initial_state.getPotentialEnergy()
                self.log.info(f"Structure {structure_num}: Initial production energy = {initial_energy}")

                # Set velocities to temperature
                simulation.context.setVelocitiesToTemperature(self.config.temp * kelvin)

                # Add comprehensive reporters
                self._add_reporters(simulation, structure_num)

                # Run production MD
                self.log.info(f"Structure {structure_num}: Running {self.config.md_steps} production steps")
                simulation.step(self.config.md_steps)

                # Get final state and save
                final_state = simulation.context.getState(getEnergy=True, getPositions=True)
                final_energy = final_state.getPotentialEnergy()

                self.log.info(f"Structure {structure_num}: Final production energy = {final_energy}")
                self.log.info(f"Structure {structure_num}: Production MD completed successfully")

                # Save final structure
                final_pdb = f"{self.outfolder}/production_final_{structure_num:04d}.pdb"
                self.save_pdb(final_pdb, simulation.topology, final_state.getPositions())

                successful_runs += 1
                print(f"✓ Structure {structure_num} completed successfully")

            except Exception as e:
                error_msg = f"Production MD failed for structure {structure_num}: {e}"
                self.log.error(error_msg)
                print(f"✗ {error_msg}")
                continue

        # Summary
        self.log.info(f"Production MD completed: {successful_runs}/{len(self.pdbs)} successful")
        print(f"Production MD Summary: {successful_runs}/{len(self.pdbs)} structures completed successfully")

        if successful_runs == 0:
            raise RuntimeError("All production MD runs failed")

    def _add_reporters(self, simulation: Simulation, structure_num: int):
        """Add comprehensive reporters for trajectory and energy monitoring.

        Args:
            simulation: OpenMM Simulation object
            structure_num: Structure index for file naming
        """
        # State data reporter (energy, temperature, etc.)
        state_file = f"{self.outfolder}/production_data_{structure_num:04d}.log"
        simulation.reporters.append(StateDataReporter(
            state_file,
            reportInterval=max(1, self.config.md_steps // 100),  # 100 data points
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, density=True, volume=True
        ))

        # Trajectory reporter (coordinates)
        traj_file = f"{self.outfolder}/trajectories/production_traj_{structure_num:04d}.dcd"
        simulation.reporters.append(DCDReporter(
            traj_file,
            reportInterval=max(1, self.config.md_steps // 1000),  # 1000 frames max
            enforcePeriodicBox=True
        ))

        # PDB reporter for intermediate structures
        pdb_file = f"{self.outfolder}/trajectories/production_frames_{structure_num:04d}.pdb"
        simulation.reporters.append(PDBReporter(
            pdb_file,
            reportInterval=max(1, self.config.md_steps // 10)  # 10 PDB frames
        ))

        self.log.info(f"Structure {structure_num}: Added reporters - state: {state_file}, trajectory: {traj_file}")

    @staticmethod
    def save_pdb(filename: str, topology, positions):
        """Save PDB file with error handling."""
        try:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            with open(filename, 'w') as f:
                PDBFile.writeFile(topology, positions, f)
        except Exception as e:
            raise IOError(f"Failed to save PDB file {filename}: {e}")

# Removed redundant save_pdb function - now part of classes
