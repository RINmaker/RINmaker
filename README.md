# RINmaker

configure:

```bash
cmake -S . -B build
```

build:

```bash
cmake --build build --target RINmaker
```

build tests:

```bash
cmake --build build --target RINmaker_tests
```

executable should be in `build/sources/Debug/`.
