# C++ STL Learning Exercises

[![CI](https://github.com/griswald/cpp-exercises/workflows/CI/badge.svg)](https://github.com/griswald/cpp-exercises/actions)

A structured curriculum for learning the C++ Standard Template Library using C++20, designed for junior developers.

## Project Structure

```
cpp-exercises/
├── week2/           # std::string and std::string_view
├── week3/           # std::array and std::span (coming soon)
├── week4/           # std::deque and std::list (coming soon)
└── ...
```

Each week is a self-contained project with:
- Exercise descriptions (`weekN_exercises.md`)
- Solution implementations (`weekN_solutions.md`)
- Test suite (GoogleTest)
- Build configuration (CMake + Pixi)

## Prerequisites

- [Pixi](https://pixi.sh) package manager

Pixi manages all dependencies including:
- C++ compiler (GCC 14)
- CMake 3.25+
- Ninja build system
- GoogleTest

## Quick Start

Navigate to any week and run:

```bash
cd week2
pixi run all
```

## Available Commands (per week)

```bash
# Run all: configure, build, and test
pixi run all

# Individual steps
pixi run configure   # Configure CMake
pixi run build       # Build tests
pixi run test        # Run tests
pixi run clean       # Clean build artifacts
```

## Curriculum Overview

See [cpp_tutoring_plan.md](cpp_tutoring_plan.md) for the complete 7-week curriculum covering:

- **Week 2**: `std::string` and `std::string_view`
- **Week 3**: `std::array` and `std::span`
- **Week 4**: `std::deque` and `std::list`
- **Week 5**: `std::map` and `std::unordered_map`
- **Week 6**: `std::set` and `std::unordered_set`
- **Week 7**: Algorithm Basics (`<algorithm>`)
- **Week 8**: More Algorithms and C++20 Ranges

## CI/CD

All weeks are automatically tested in GitHub Actions CI. The workflow:
1. Checks out code
2. Sets up Pixi environment
3. Configures, builds, and tests each week independently

## Contributing

When adding a new week:

1. Create `weekN/` directory with the standard structure
2. Add `weekN` to the matrix in `.github/workflows/ci.yml`
3. Ensure `pixi run all` works locally before pushing

## License

Educational use only.
