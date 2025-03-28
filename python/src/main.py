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
import gc
import toml
import time
import dncs
import logging
from typing import List
from dataclasses import dataclass
from integrator import OptimizedMolecularDynamics, MDSimulation

@dataclass
class SimulationConfig:
    sequence: str
    moleculename: str
    folder: str
    interface: str
    n_samples: int
    md_simulation: int
    temp: float
    forcefield: List[str]
    device: str
    solvent: int
    steps: int
    gamma: float
    dt: float
    md_steps: int
    method: str

class GenerateSamples:
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.generate_samples()

    def generate_samples(self):
        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        logger = logging.getLogger(__name__)

        try:
            # Ensure output directory exists
            output_dir = f"{self.config.folder}/{self.config.moleculename}"
            os.makedirs(output_dir, exist_ok=True)

            # Generate polymer samples
            sample = dncs.Polymer(self.config.sequence, "amberfb15.xml")
            dncs.SobolSampler(
                sample,
                self.config.n_samples,
                self.config.method,
                self.config.temp,
                output_dir
            )
            logger.info("Sample generation completed successfully")

        except Exception as e:
            logger.error(f"Error in sample generation: {e}")
            raise

def run_openmm(config: SimulationConfig):
    """
    Run OpenMM molecular dynamics simulation

    Args:
        config (SimulationConfig): Simulation configuration
    """
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    try:
        # Run molecular dynamics integrator
        logger.info("Starting molecular dynamics simulation")
        dynamics = OptimizedMolecularDynamics(config)
        dynamics.run_parallel_simulation()

        # Run MD simulation
        logger.info("Starting long MD simulation")
        md_sim = MDSimulation(config)
        md_sim.run_simulation()

        logger.info("Molecular dynamics simulation completed")

    except Exception as e:
        logger.error(f"Error in OpenMM simulation: {e}")
        raise

def main():
    """
    Main entry point for DNCS simulation
    """
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('simulation.log'),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)

    try:
        # Load configuration from TOML file
        config_path = 'dncs.toml'
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file {config_path} not found")

        with open(config_path, 'r') as f:
            toml_config = toml.load(f)

        # Create configuration object
        config = SimulationConfig(**toml_config['simulation'])

        # Start timing
        start_time = time.time()

        # Generate conformational samples
        logger.info("Starting sample generation")
        GenerateSamples(config)

        # Run simulation if interface is OpenMM
        if config.interface == "openmm":
            run_openmm(config)

        # Calculate and log total simulation time
        end_time = time.time()
        simulation_time = end_time - start_time

        if simulation_time > 60:
            minutes = simulation_time / 60
            logger.info(f"Total simulation time: {minutes:.5f} minutes")
        else:
            logger.info(f"Total simulation time: {simulation_time:.5f} seconds")

    except Exception as e:
        logger.error(f"Simulation failed: {e}")

    finally:
        # Final cleanup
        gc.collect()

if __name__ == "__main__":
    main()
