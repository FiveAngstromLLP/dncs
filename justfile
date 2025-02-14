
@run:
    python python/src/main.py

@install:
    pip install -r python/requirements.txt
    pip install maturin
    maturin develop --release -m python/Cargo.toml
