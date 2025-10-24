# Local Testing Guide - Week 3

This guide contains instructions for running the full test suite locally after completing the exercises on Godbolt.

## Prerequisites

- [Pixi](https://pixi.sh) package manager

Pixi manages all dependencies including:
- C++ compiler (GCC 14)
- CMake 3.25+
- Ninja build system
- GoogleTest 1.14

## Quick Start

```bash
# Run all: configure, build, and test
pixi run all
```

## Available Commands

### Configure the build
```bash
pixi run configure
```
This runs: `cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release`

### Build the project
```bash
pixi run build
```
This runs: `cmake --build build --parallel`

### Run tests
```bash
pixi run test
```
This runs: `ctest --test-dir build --output-on-failure`

### Clean build artifacts
```bash
pixi run clean
```

## Manual Build (without Pixi tasks)

If you prefer to run commands manually:

```bash
# Enter the pixi environment
pixi shell

# Configure
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build --parallel

# Test
ctest --test-dir build --output-on-failure
```

## Project Structure

```
week3/
├── pixi.toml                       # Pixi project configuration
├── pixi.lock                       # Locked dependency versions
├── CMakeLists.txt                  # CMake build configuration
├── test_transform3d.cpp            # Exercise 1 tests
├── test_sensor_array.cpp           # Exercise 2 tests
├── test_trajectory.cpp             # Exercise 3 tests
├── test_multi_robot.cpp            # Exercise 4 tests
├── README.md                       # Exercise descriptions
├── week3_solutions.md              # Solutions with explanations
└── TESTING.md                      # This file
```

## Compiler Requirements

- **C++20 support required**
- Tested with GCC 14.3.0 (provided by Pixi)
- Also works with Clang 14+ and MSVC 19.29+

## Troubleshooting

### CMake cache issues
If you see errors about CMakeCache.txt directory mismatches:
```bash
pixi run clean
pixi run configure
pixi run build
```

### Pixi environment issues
```bash
# Remove the Pixi environment and reinstall
rm -rf .pixi
pixi install
```

### Test failures
Run tests with verbose output:
```bash
pixi shell
ctest --test-dir build --output-on-failure --verbose
```

## CI/CD

This project is tested automatically in GitHub Actions. See `.github/workflows/ci.yml` in the repository root.
