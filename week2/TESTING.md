# Local Testing Guide - Week 2

This guide contains instructions for running the full test suite locally after completing the exercises on Godbolt.

## Prerequisites

- [Pixi](https://pixi.sh) package manager

Pixi manages all dependencies including:
- C++ compiler (GCC 14)
- CMake 3.25+
- Ninja build system
- GoogleTest 1.17

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
week2/
├── pixi.toml                       # Pixi project configuration
├── pixi.lock                       # Locked dependency versions
├── CMakeLists.txt                  # CMake build configuration
├── test_normalize_whitespace.cpp   # Exercise 1 tests
├── test_split.cpp                  # Exercise 2 tests
├── test_serial_protocol.cpp        # Exercise 3 tests
├── test_command_joint_state.cpp    # Exercise 4 tests
├── week2_exercises.md              # Exercise descriptions
├── week2_solutions.md              # Solutions with explanations
├── README.md                       # Exercise descriptions
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
