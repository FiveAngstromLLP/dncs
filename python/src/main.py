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


import toml
import time
import dncs
from typing import List
from dataclasses import dataclass
from integrator import DncsIntegrator, CleanUp, MDSimulation

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
        sample = dncs.Polymer(self.config.sequence, "amberfb15.xml")
        dncs.SobolSampler(
            sample,
            self.config.n_samples,
            self.config.method,
            self.config.temp,
            f"{self.config.folder}/{self.config.moleculename}"
        )

def runopenmm(config):
    # Run Integrator
    integrator = DncsIntegrator(config)
    integrator.run_integrator()
    CleanUp(config)

    md = MDSimulation(config)
    md.run_simulation()



if __name__ == "__main__":
    # Load configuration from dncs.toml file
    with open('dncs.toml', 'r') as f:
        toml_config = toml.load(f)

    config = SimulationConfig(**toml_config['simulation'])

    start_time = time.time()


    #Generate conformational samples
    GenerateSamples(config)
    print("Completed Sampling")

    if config.interface == "openmm":
        runopenmm(config)

    end_time = time.time()
    simulation_time = end_time - start_time

    if simulation_time > 60:
        minutes = simulation_time / 60
        print(f"Total simulation time: {minutes:.5f} minutes")
    else:
        print(f"Total simulation time: {simulation_time:.5f} seconds")

    # # Final cleanup
    import gc
    gc.collect()
