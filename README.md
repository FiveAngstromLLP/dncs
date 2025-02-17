
## Compile dncs LIBRARY
```bash
cargo build -p cabi --release
```

## Compile CPP LIBRARY
```bash
g++ main.cpp -o main -ldncs -L ./target/release -Wl,-rpath=$PWD/target/release
```
