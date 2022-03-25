# RINmaker

configure:

```bash
cmake -S . -B build
```

configure tests:

```bash
cmake -S . -B build -DUSE_TESTS=ON
```

build:

```bash
cmake --build build
```

executable should be in `build/sources/Debug/`.
