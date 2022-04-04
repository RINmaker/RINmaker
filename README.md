# RINmaker

### How to build

Configure project:

```bash
cmake -S . -B cmake-build-debug
```

#### App

```bash
cmake --build cmake-build-debug --target RINmaker
```

#### Test

```bash
cmake --build cmake-build-debug --target RINmaker_test
```

Executable will be in `cmake-build-debug/app/RINmaker`.
