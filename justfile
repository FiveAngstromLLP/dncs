@run *ARGS:
    cargo run -p dncs-cli {{ARGS}}

@lib ARGS:
    cd dncs-lib && cargo run --example {{ARGS}}

@docs:
    cd docs && mdbook serve
