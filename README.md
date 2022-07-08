# RINmaker

## Clone this repo

```bash
git clone -b dev https://github.com/RINmaker/RINmaker.git
```

```bash
cd RINmaker
git submodule init
git submodule update
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
