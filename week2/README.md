# Week 2: Serial Protocol Parser - Testing Project

This project contains tests for the modern C++20 serial protocol parser implementation from Week 2 exercises.

## Prerequisites

- [Pixi](https://pixi.sh) package manager

## Quick Start

All dependencies are managed by Pixi. Simply run:

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
.
├── pixi.toml                    # Pixi project configuration
├── CMakeLists.txt               # CMake build configuration
├── test_serial_protocol.cpp     # GTest test suite
├── week2_exercises.md           # Exercise descriptions
├── week2_solutions.md           # Solutions
└── README.md                    # This file
```

## What's Tested

The test suite validates the modern C++20 implementation of `read_joint_state()`:

- ✅ Valid message parsing with correct checksum
- ✅ Valid message without newline terminator
- ✅ Invalid header detection
- ✅ Invalid checksum detection
- ✅ Malformed data handling
- ✅ Whitespace trimming

## Implementation Details

The implementation uses modern C++20 features:
- `std::ranges` and views for lazy evaluation
- `std::from_chars` for fast, exception-free parsing
- `std::optional` for error handling
- `constexpr` for compile-time constants
- Alternative operators (`and`, `or`, `not`)
