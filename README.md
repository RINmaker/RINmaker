# RINmaker

configure:

```bash
cmake -S . -B build
```

make_instance:

```bash
cmake --build build --target RINmaker
```

make_instance tests:

```bash
cmake --build build --target RINmaker_test
```

executable should be in `make_instance/sources/Debug/`.
