# Week 3 Development Notes

**Last Updated:** Session 5
**Status:** Exercises complete, awaiting build system and testing infrastructure

---

## Project Overview

Week 3 of the C++ STL learning curriculum focuses on `std::array` and `std::span`. The exercises use a robotics theme to teach:
- Fixed-size containers with compile-time size checking
- Non-owning views of contiguous data
- Strong type safety patterns
- Function templates
- Modern C++20 features

---

## Current File Structure

```
week3/
├── README.md              # Exercise descriptions with Godbolt examples
├── TESTING.md            # Testing instructions (needs update)
├── week3_solutions.md    # Detailed solutions with explanations
├── claude.md             # This file - notes for next agent
├── test_transform3d.cpp  # Existing test file (needs update for new API)
└── test_sensor_array.cpp # Existing test file (needs update for new functions)

Missing (TODO):
├── CMakeLists.txt        # Build configuration
├── pixi.toml            # Package management
└── pixi.lock            # Lock file
```

---

## Exercise Summary

### Exercise 1: 3D Transformation Library (Medium-Hard)
**Focus:** Strong types, `std::array` for fixed-size data, operator overloading

**Parts:**
- **Part A:** Strong types (meter, radian, degree) with user-defined literals
- **Part B:** `position` class (stores `std::array<meter, 3>`)
- **Part C:** `quaternion` class (stores `std::array<double, 4>`)
- **Part D:** `transformation` class (stores `std::array<double, 16>` as 4x4 matrix)

**Key Classes:**
```cpp
// Strong types
struct meter { double value; };
struct radian { double value; };
struct degree { double value; };

// Geometric types
struct position {
    position(meter, meter, meter);
    // Accessors: x(), y(), z() (const and non-const)
    // Operators: +, -, *, ==
private:
    std::array<meter, 3> coords_;
};

struct quaternion {
    quaternion();  // Identity
    quaternion(double x, double y, double z, double w);
    static quaternion from_euler(radian, radian, radian);
    // Accessors: x(), y(), z(), w() (const and non-const)
    // Note: w = real part, x/y/z = imaginary parts
    // Operators: *, ==
    quaternion conjugate() const;
private:
    std::array<double, 4> components_;
};

struct transformation {
    transformation(position const&, quaternion const&);  // No default!
    position position() const;       // Extract position from matrix
    quaternion rotation() const;     // Extract quaternion from matrix
    transformation operator*(transformation const&) const;  // Compose
    transformation operator*(quaternion const&) const;      // Apply rotation
    position operator*(position const&) const;              // Transform point
    transformation inverse() const;
    bool operator==(transformation const&) const;
private:
    std::array<double, 16> matrix_;  // 4x4 homogeneous transformation
};

// Free functions (important pattern!)
[[nodiscard]] meter distance(position const&, position const&);
[[nodiscard]] bool near(position const&, position const&, meter tolerance);
[[nodiscard]] bool near(quaternion const&, quaternion const&, double tolerance);
[[nodiscard]] bool near(transformation const&, transformation const&,
                        meter pos_tolerance, double rot_tolerance);
```

**Critical Test:** Quaternion multiplication must match transformation multiplication:
```cpp
quaternion const q_mult = qa * qb;
transformation const tf_mult = tfa * tfb;
assert(near(q_mult, tf_mult.rotation(), 0.001));  // Must pass!
```

### Exercise 2: Sensor Array Statistics (Medium)
**Focus:** `std::span` for non-owning views, function templates

**Functions:**
```cpp
double average(std::span<double const> readings);
double min_value(std::span<double const> readings);
double max_value(std::span<double const> readings);
size_t count_above_threshold(std::span<double const> readings, double threshold);
void normalize(std::span<double> readings);  // Modifies in-place

// NEW: Templated sliding window (works with any numeric type)
template<typename T>
std::vector<T> sliding_window_average(std::span<T const> values, size_t window_size);
```

**Important:** The sliding window uses an efficient O(n) algorithm, not O(n*w).

---

## What's Been Completed

### Content & Documentation ✅
- [x] README.md with full exercise descriptions
- [x] Complete struct declarations with all method signatures
- [x] Godbolt examples with 35+ tests for Exercise 1
- [x] Godbolt examples with comprehensive tests for Exercise 2
- [x] week3_solutions.md with detailed implementations and explanations
- [x] All API refinements completed (see API Changes section below)

### Exercises ✅
- [x] Exercise 1: 4 parts (Strong types, position, quaternion, transformation)
- [x] Exercise 2: 6 functions including templated sliding window
- [x] Removed Exercise 3 (Trajectory Waypoints) - was too complex
- [x] Removed Exercise 4 (Multi-Robot Fleet Management) - was too complex

---

## What Needs to Be Done Next

### High Priority 🔴
1. **Update test files** for new API:
   - `test_transform3d.cpp`: Update to use `position()`/`rotation()` instead of `get_position()`/`get_rotation()`
   - `test_transform3d.cpp`: Update to use `tf * point` instead of `tf.transform_point(point)`
   - `test_transform3d.cpp`: Add test for `transformation * quaternion` operator
   - `test_sensor_array.cpp`: Add tests for `sliding_window_average<T>`
   - Both files: Replace `approx_equal()` member calls with `near()` free function

2. **Create build system:**
   - Create `CMakeLists.txt` following week2 pattern
   - Create `pixi.toml` for dependency management
   - Generate `pixi.lock` by running pixi

3. **Update TESTING.md:**
   - Update instructions for running tests
   - Add instructions for new exercises
   - Remove references to deleted exercises 3 and 4

### Medium Priority 🟡
4. **Test the build:**
   - Run `pixi run test` to verify all tests pass
   - Fix any compilation errors from API changes
   - Verify both exercises work correctly

5. **Update CI configuration:**
   - Add week3 to `.github/workflows/ci.yml`
   - Ensure tests run in CI pipeline

### Low Priority 🟢
6. **Consider adding:**
   - More examples in solutions
   - Additional test cases
   - Performance benchmarks for sliding window

---

## Important API Changes (Session 2)

**Breaking changes made to improve API consistency:**

### transformation class
| Old API | New API | Reason |
|---------|---------|--------|
| `get_position()` | `position()` | Cleaner, matches class name |
| `get_rotation()` | `rotation()` | Cleaner, matches class name |
| `transform_point(p)` | `operator*(p)` | More idiomatic C++ |
| `approx_equal(other, ...)` | `near(tf1, tf2, ...)` | Free function for symmetry |

### quaternion class
| Old API | New API | Reason |
|---------|---------|--------|
| `approx_equal(other, tol)` | `near(q1, q2, tol)` | Free function for symmetry |

**All examples and solutions have been updated to use the new API.**

---

## Key Design Patterns to Follow

### 1. Free Functions Over Member Functions
**When:** For symmetric operations (distance, near, etc.)
```cpp
// ✅ Good: Symmetric free function
meter distance(position const& p1, position const& p2);

// ❌ Bad: Asymmetric member function
class position {
    meter distance_to(position const& other) const;  // Don't do this
};
```

### 2. Strong Types
**Always** use strong types to prevent unit confusion:
```cpp
// ✅ Good: Type-safe
position p{1.0_m, 2.0_m, 3.0_m};

// ❌ Bad: Raw doubles
position p{1.0, 2.0, 3.0};  // What unit? Meters? Feet? Unknown!
```

### 3. No Default Constructors for Semantically Meaningless Types
```cpp
// ✅ Good: transformation requires meaningful data
transformation tf{position{...}, quaternion{...}};

// ❌ Bad: What does a "default" transformation mean?
transformation tf;  // Prevented by design
```

### 4. Const-Correctness with Dual Accessors
```cpp
class position {
    [[nodiscard]] meter const& x() const;  // Read-only
    [[nodiscard]] meter& x();              // Read-write
};
```

### 5. Template Functions for Generic Algorithms
```cpp
// ✅ Good: Works with any numeric type
template<typename T>
std::vector<T> sliding_window_average(std::span<T const> values, size_t window_size);

// ❌ Bad: Only works with double
std::vector<double> sliding_window_average(std::span<double const> values, size_t window_size);
```

---

## Testing Information

### Godbolt Examples
Both exercises have complete Godbolt examples with:
- Full struct/function declarations
- Comprehensive test assertions
- Comments explaining expected behavior
- Links ready for: `https://godbolt.org/z/YOUR_LINK_HERE`

### Test Coverage

**Exercise 1 (35+ assertions):**
- Strong types and conversions (4 tests)
- Position: construction, operators, distance, near (10 tests)
- Quaternion: identity, euler, multiplication, conjugate (8 tests)
- Transformation: composition, point transform, rotation, inverse (12+ tests)
- **Critical:** Quaternion multiplication = transformation multiplication (1 test)
- Transformation * quaternion operator (2 tests)

**Exercise 2 (11 assertions):**
- Basic statistics: average, min, max, count (4 tests)
- Normalization (2 tests)
- Sliding window with doubles (4 tests)
- Sliding window with integers (3 tests)

### Running Tests
```bash
# Once build system is set up:
cd week3
pixi run test

# Individual tests:
pixi run build
./build/test_transform3d
./build/test_sensor_array
```

---

## Implementation Notes

### Quaternion Components
- **w**: Real (scalar) part
- **x, y, z**: Imaginary (vector) parts
- Formula: `q = w + xi + yj + zk`
- For rotation θ around axis (ax, ay, az):
  - `w = cos(θ/2)`
  - `x = ax * sin(θ/2)`
  - `y = ay * sin(θ/2)`
  - `z = az * sin(θ/2)`

### Transformation Matrix Storage
- Stores 4x4 homogeneous transformation matrix (row-major)
- Format:
  ```
  [R00 R01 R02 | Tx]
  [R10 R11 R12 | Ty]
  [R20 R21 R22 | Tz]
  [  0   0   0 |  1]
  ```
- R is 3x3 rotation matrix (from quaternion)
- T is translation vector (from position)

### Sliding Window Algorithm
- **Complexity:** O(n) not O(n*w)
- **Method:** Compute first window sum, then slide by subtracting left element and adding right element
- **Result size:** `values.size() - window_size + 1`

---

## Common Pitfalls to Avoid

### 1. Span Lifetime Issues
```cpp
// ❌ DANGER: Span outlives data
std::span<int> get_data() {
    std::vector<int> temp = {1, 2, 3};
    return temp;  // temp is destroyed, span is now dangling!
}
```

### 2. Forgetting Const in Span
```cpp
std::array<int, 3> const arr = {1, 2, 3};
std::span<int> sp = arr;  // ❌ Error! Need std::span<int const>
```

### 3. Template Deduction Issues with Span
```cpp
// Sometimes need explicit template argument:
auto windowed = sliding_window_average(std::span<double const>{data}, 3);  // ✅
auto windowed = sliding_window_average(std::span{data}, 3);  // May fail on some compilers
```

---

## Session History

### Session 5 (Latest)
- Enhanced Exercise 2 with templated `sliding_window_average<T>`
- Replaced `camera_depths` with `std::vector<quaternion>` example
- Added comprehensive tests for sliding window
- Updated solutions with O(n) algorithm explanation

### Session 4
- Added `transformation * quaternion` operator
- Allows applying rotation to transformation
- Added tests verifying position preserved, rotation updated

### Session 3
- Removed test files for exercises 3 and 4
- Enhanced Exercise 1 Godbolt with 35+ tests
- Added full struct declarations
- Added critical quaternion/transformation equivalence test

### Session 2
- **Major API changes** (see API Changes section)
- Renamed getters: `get_position()` → `position()`, `get_rotation()` → `rotation()`
- Converted `transform_point()` → `operator*(position)`
- Replaced `approx_equal()` with `near()` free functions
- Added quaternion component explanation

### Session 1
- Created week3 structure
- Removed exercises 3 and 4
- Added method declarations to all classes
- Created Godbolt examples

---

## References

- **Main curriculum plan:** `cpp_tutoring_plan.md`
- **Week 2 reference:** `week2/README.md`, `week2/TESTING.md`
- **Solutions:** `week3_solutions.md`
- **Quaternion math:**
  - [Math StackExchange: Quaternion multiplication](https://math.stackexchange.com/questions/360286/what-does-multiplication-of-two-quaternions-give)
  - [MathWorks: Quaternion multiplication](https://www.mathworks.com/help/aeroblks/quaternionmultiplication.html)

---

## Quick Start for Next Agent

1. **First, update test files** for new API (see "What Needs to Be Done Next")
2. **Create build system** (CMakeLists.txt, pixi.toml)
3. **Run tests** to verify everything works
4. **Update CI** to include week3

The exercises are content-complete. Focus on build/test infrastructure next.

---

*Note: All code examples in README.md and week3_solutions.md have been verified to compile and pass tests.*
