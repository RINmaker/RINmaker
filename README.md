# RINmaker

## Clone this repo

```bash
git clone -b dev https://github.com/RINmaker/RINmaker.git
git submodules init
git submodules update
```

## Configure project:

```bash
cmake -S . -B build
```

## App

#### Build

```bash
cmake --build build --target RINmaker
```

#### Run

```bash
./build/app/RINmaker -h
```

## Test

#### Build

```bash
cmake --build build --target RINmaker_test
```

#### Run

```bash
./build/test/RINmaker_test
```
