
@run:
    python python/src/main.py

@install:
    pip3 install -r python/requirements.txt
    pip3 install maturin
    maturin develop --release -m python/Cargo.toml
