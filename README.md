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
cmake --build build --target RINmaker_test
```

executable should be in `build/sources/Debug/`.
