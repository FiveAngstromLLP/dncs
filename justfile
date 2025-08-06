# set shell:=["poweshell.exe","-c"]
@run:
    python python/src/main.py

@install:
    pip install -r python/requirements.txt
    maturin develop --release -m python/Cargo.toml

@install-dncs:
    maturin develop --release -m python/Cargo.toml

# View top N samples in PyMOL
@view FOLDERNAME TOP_N:
    python scripts/pymol_script.py "{{FOLDERNAME}}" "{{TOP_N}}"

# List available result folders and files
@list *FOLDERNAME:
    python scripts/list_results.py {{FOLDERNAME}}

# For Windows, use this command to install the dependencies:
#just --shell powershell.exe --shell-arg -c install
# cargo install maturin
# maturin develop --release -m python/Cargo.toml
