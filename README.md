# RINmaker

configure:

```bash
cmake -S . -B make_instance
```

make_instance:

```bash
cmake --build make_instance --target RINmaker
```

make_instance tests:

```bash
cmake --build make_instance --target RINmaker_test
```

executable should be in `make_instance/sources/Debug/`.
