# set shell:=["poweshell.exe","-c"]
@run:
    python python/src/main.py

@install:
    pip install -r python/requirements.txt
    maturin develop --release -m python/Cargo.toml

@install-dncs:
    maturin develop --release -m python/Cargo.toml
# For Windows, use this command to install the dependencies:
#just --shell powershell.exe --shell-arg -c install
# cargo install maturin
# maturin develop --release -m python/Cargo.toml
