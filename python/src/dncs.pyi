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

"""Python interface for dncs."""

def getPDB(seq: str, filename: str) -> None:
    """Generate PDB file from amino acid sequence.

    Args:
        seq: Amino acid sequence string
        filename: Output PDB filename
    """
    ...

def pdb_to_angle(filename: str) -> str:
    """Convert PDB file to dihedral angles.

    Args:
        filename: Input PDB filename

    Returns:
        str: Comma-separated string of angles
    """
    ...

class Polymer:
    """Represents a polymer system with force field parameters."""

    def __init__(self, seq: str, forcefield: str) -> None:
        """Initialize polymer from amino acid sequence.

        Args:
            seq: Amino acid sequence string
            forcefield: Force field name, must be one of:
                       - amber03.xml
                       - amber10.xml
                       - amber96.xml
                       - amber99sb.xml
                       - amberfb15.xml
        """
        ...

    def getEnergy(self) -> float:
        """Calculate energy of the polymer system.

        Returns:
            float: Energy value
        """
        ...

    def toPDB(self, filename: str) -> None:
        """Save polymer structure to PDB file.

        Args:
            filename: Output PDB filename
        """
        ...

    def dihedral(self, foldername: str) -> None:
        """Log dihedral angles to a folder.

        Args:
            foldername: Output folder path
        """
        ...

class SobolSampler:
    """Sobol sequence sampler for polymer conformations."""

    def __init__(self, system: Polymer, no_of_samples: int, grid: int, temp: float) -> None:
        """Initialize Sobol sampler.

        Args:
            system: Polymer system to sample
            no_of_samples: Number of conformations to generate
        """
        ...

    def write_angles(self, filename: str) -> None:
        """Write sampled angles to file.

        Args:
            filename: Output filename for angles
        """
        ...

    def toPDB(self, filename: str) -> None:
        """Save sampled conformations to PDB file.

        Args:
            filename: Output PDB filename
        """
        ...

    def toPDBFiles(self, prefix: str) -> None:
        """Save each sampled conformation to separate PDB files.

        Args:
            prefix: Prefix for output PDB filenames
        """
        ...
