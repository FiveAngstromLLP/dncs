

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

    def __init__(self, system: Polymer, no_of_samples: int) -> None:
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
