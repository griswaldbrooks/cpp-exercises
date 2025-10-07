# Week 2 Serial Protocol Parser - Development Context

## Overview

This is a C++20 serial protocol parser implementation for a robotics joint state communication protocol. The project demonstrates modern C++ features including ranges, views, `std::optional`, `std::from_chars`, and const-correctness with east const style.

## Project Structure

```
week2/
├── pixi.toml                    # Pixi package manager config
├── CMakeLists.txt               # CMake build configuration
├── test_serial_protocol.cpp     # Implementation + GTest test suite
├── week2_exercises.md           # Student exercises
├── week2_solutions.md           # Reference solutions
├── README.md                    # Build/run instructions
└── claude.md                    # This file
```

## Protocol Specification

**Format:** `JS:<angle1>,<angle2>:<checksum>\n`

- `JS` - Header (2 bytes)
- `:` - Section delimiter (separates major protocol parts)
- `<angle1>,<angle2>` - Two joint angles in radians (floating point)
- `,` - Angle delimiter (separates angle values)
- `:` - Section delimiter before checksum
- `<checksum>` - 8-bit checksum as hex (0xXX format)
- `\n` - Optional newline terminator

**Checksum:** Sum all ASCII byte values from start of message (including "JS:") up to but not including the last `:`, then modulo 256.

**Example:** `JS:1.57,-0.785:0xFD\n`

## Code Style Guidelines

### Naming Conventions
- **Types:** `snake_case` (e.g., `joint_state`, `protocol_config`)
- **Variables:** `snake_case` with east const (e.g., `std::string_view const message`)
- **Constants:** `constexpr` with descriptive names (e.g., `constexpr char section_delimiter`)

### Modern C++ Patterns Used
- **East const:** `Type const` instead of `const Type`
- **Alternative operators:** `and`, `or`, `not` instead of `&&`, `||`, `!`
- **Ranges & views:** Lazy evaluation with `std::views::split`, `std::views::transform`
- **`std::from_chars`:** Fast, exception-free parsing (instead of `std::stod`)
- **`std::optional`:** Error handling without exceptions
- **`std::string_view`:** Zero-copy string operations
- **Structured bindings:** `auto const [ptr, ec] = std::from_chars(...)`

### Key Implementation Details

1. **Views are NOT const:** View variables must be mutable (`auto view = ...` not `auto const view = ...`) because `std::ranges::begin/end` requires non-const access.

2. **Protocol is configurable:** All protocol parameters are in `protocol_config` struct:
   ```cpp
   struct protocol_config {
       char section_delimiter = ':';
       char angle_delimiter = ',';
       std::string_view header = "JS";
       size_t expected_parts = 3;
       size_t expected_angles = 2;
   };
   ```

3. **Helper functions:**
   - `trim()`: Remove leading/trailing whitespace (constexpr)
   - `split_to_vector()`: Split string_view into vector using ranges
   - `calculate_checksum()`: Sum ASCII bytes modulo 256

4. **Validation order:**
   1. Trim whitespace
   2. Validate header
   3. Split by section delimiter
   4. Check expected parts count
   5. Split angles by angle delimiter
   6. Check expected angles count
   7. Parse angles with `std::from_chars`
   8. Validate checksum

## Building and Testing

```bash
# All in one (configure + build + test)
pixi run all

# Individual steps
pixi run configure    # cmake -S . -B build -GNinja
pixi run build        # cmake --build build --parallel
pixi run test         # ctest --test-dir build --output-on-failure
pixi run clean        # rm -rf build
```

## Current Test Coverage

10 tests, all passing:
- ✅ ValidMessage
- ✅ ValidMessageNoNewline
- ✅ InvalidHeader
- ✅ InvalidChecksum
- ✅ MalformedData
- ✅ TrimWhitespace
- ✅ CustomDelimiters (using `|` and `;`)
- ✅ WrongExpectedParts (misconfigured)
- ✅ WrongExpectedAngles (misconfigured)
- ✅ WrongHeader (misconfigured)

## Common Gotchas

1. **Checksums must be calculated correctly:**
   - Include "JS:" header
   - Include all data up to (but not including) the last `:`
   - Example: For `JS:1.57,-0.785:0xFD`, checksum `"JS:1.57,-0.785"`

2. **View lifetime issues:**
   - Views created from `std::views::split` are lazy and don't allocate
   - Must materialize into container (e.g., `std::vector`) to use later
   - Views can't be `const` when passed to `std::ranges::begin/end`

3. **`std::from_chars` for floating point:**
   - Available in C++17 but implementations vary
   - Check `std::errc{}` for errors, not exceptions
   - Works directly with `std::string_view::data()` and size

## Next Steps / Potential Improvements

1. **Add `command_joint_state()` function** - Inverse operation to write messages
2. **Support variable number of joints** - Make angles a vector instead of fixed pair
3. **Better error reporting** - Return error details instead of just `std::nullopt`
4. **Benchmark performance** - Compare ranges vs manual parsing
5. **Add checksum_type enum** - Support different checksum algorithms
6. **Header validation** - Support multiple header types
7. **Add to solutions.md** - Update with final implementation

## Dependencies

- CMake >= 3.25
- C++20 compiler (GCC 13+, Clang 16+, or equivalent)
- Ninja build system
- GTest >= 1.14
- Pixi for environment management

## Educational Context

This is week 2 of a C++ tutoring curriculum focused on:
- `std::string` and `std::string_view`
- String parsing and manipulation
- Modern C++20 features (ranges, views)
- Real-world protocol design (checksums, delimiters, validation)

Week 1 covered `std::vector`. Future weeks will cover containers, algorithms, and ranges in more depth.
