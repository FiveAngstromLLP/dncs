# Digital Nets Conformational Sample (DNCS) CLI

A command-line tool for molecular conformational sampling using Digital Nets.

## Installation

```bash
cargo install --path .
```

## Usage

```bash
dncs [OPTIONS]
```

## Options

- `-c, --config`: Generates a dncs.json configuration file
- `-N, --molecule <NAME>`: Molecule name
- `-s, --sequence <SEQ>`: Amino acid sequence
- `-n, --samples <NUM>`: Number of samples to generate
- `-f, --forcefield <FF>`: Force field selection:
  - amber03.xml
  - amber10.xml
  - amber96.xml
  - amber99sb.xml
  - amberfb15.xml
- `-m, --minimize`: Enable energy minimization
- `-g, --grid <NUM>`: Number of grids to divide the sample

## Example

Run sampling for a molecule:

```bash
dncs -N test -s YGGFM -n 10 -f amber03.xml -m -g 4
```

Or use a configuration file:

```bash
# Generate sample config
dncs -c

# Edit dncs.json as needed
dncs # Run with config file
```

## Config File Format

The `dncs.json` configuration file has the following format:

```json
{
  "Generate": {
    "molecule": "Sample",
    "sequence": "YGGFM",
    "n_samples": 10,
    "forcefield": "amberfb15.xml",
    "minimize": true,
    "grid": 4
  }
}
```
