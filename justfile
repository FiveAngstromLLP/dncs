
@run:
    python python/src/main.py

@install:
    pip install -r python/requirements.txt
    maturin develop --release -m python/Cargo.toml

@install-dncs:
    maturin develop --release -m python/Cargo.toml
