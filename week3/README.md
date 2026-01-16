# Week 3: `std::array` and `std::span` - Exercises

**Recommended:** Use [Godbolt.org](https://godbolt.org) to complete these exercises! Each exercise includes a pre-configured Godbolt link with starter code and tests.

Once you're ready to run the full test suite locally, see the [Local Testing Guide](TESTING.md).

---

## Exercise 1: 3D Transformation Library with Strong Types (Medium-Hard)

Robotics systems need robust 3D transformation representations with **type safety**. Build a transformation library that uses **strong types** to prevent unit errors at compile time.

### Part A: Strong Types with User-Defined Literals

**Key Concept: Strong Types**

A common source of bugs in robotics and engineering software is mixing units (e.g., passing degrees where radians are expected, or centimeters where meters are expected). C++ allows us to create **strong types** that catch these errors at compile time.

Create three strong type wrappers in `namespace literals`:

```cpp
namespace literals {
    // Strong type for distances
    struct meter_t {
        double value;
        constexpr explicit meter_t(double v) : value{v} {}

        // Arithmetic operators
        [[nodiscard]] constexpr meter_t operator+(meter_t const& other) const {
            return meter_t{value + other.value};
        }
        [[nodiscard]] constexpr meter_t operator-(meter_t const& other) const {
            return meter_t{value - other.value};
        }
        [[nodiscard]] constexpr meter_t operator*(double scalar) const {
            return meter_t{value * scalar};
        }
        [[nodiscard]] constexpr meter_t operator/(double scalar) const {
            return meter_t{value / scalar};
        }

        // Comparison operators (using spaceship operator)
        [[nodiscard]] constexpr auto operator<=>(meter_t const& other) const = default;
    };

    // Strong type for angles in radians
    struct radian_t {
        double value;
        constexpr explicit radian_t(double v) : value{v} {}

        // Arithmetic operators
        [[nodiscard]] constexpr radian_t operator+(radian_t const& other) const {
            return radian_t{value + other.value};
        }
        [[nodiscard]] constexpr radian_t operator-(radian_t const& other) const {
            return radian_t{value - other.value};
        }
        [[nodiscard]] constexpr radian_t operator*(double scalar) const {
            return radian_t{value * scalar};
        }
        [[nodiscard]] constexpr radian_t operator/(double scalar) const {
            return radian_t{value / scalar};
        }

        // Comparison operators
        [[nodiscard]] constexpr auto operator<=>(radian_t const& other) const = default;
    };

    // Strong type for angles in degrees
    struct degree_t {
        double value;
        constexpr explicit degree_t(double deg) : value{deg} {}

        // Arithmetic operators
        [[nodiscard]] constexpr degree_t operator+(degree_t const& other) const {
            return degree_t{value + other.value};
        }
        [[nodiscard]] constexpr degree_t operator-(degree_t const& other) const {
            return degree_t{value - other.value};
        }
        [[nodiscard]] constexpr degree_t operator*(double scalar) const {
            return degree_t{value * scalar};
        }
        [[nodiscard]] constexpr degree_t operator/(double scalar) const {
            return degree_t{value / scalar};
        }

        // Comparison operators
        [[nodiscard]] constexpr auto operator<=>(degree_t const& other) const = default;
    };

    // Conversion functions (work with both types)
    constexpr radian_t to_radians(degree_t const d) {
        return radian_t{d.value * std::numbers::pi / 180.0};
    }
    constexpr degree_t to_degrees(radian_t const r) {
        return degree_t{r.value * 180.0 / std::numbers::pi};
    }

    // User-defined literals
    constexpr meter_t operator""_m(long double v) { return meter_t{static_cast<double>(v)}; }
    constexpr meter_t operator""_m(unsigned long long v) { return meter_t{static_cast<double>(v)}; }

    constexpr radian_t operator""_rad(long double v) { return radian_t{static_cast<double>(v)}; }
    constexpr radian_t operator""_rad(unsigned long long v) { return radian_t{static_cast<double>(v)}; }

    constexpr degree_t operator""_deg(long double v) { return degree_t{static_cast<double>(v)}; }
    constexpr degree_t operator""_deg(unsigned long long v) { return degree_t{static_cast<double>(v)}; }
}
```

**Why Strong Types Matter:**

```cpp
using namespace literals;

// With plain doubles (DANGEROUS):
void set_position(double x, double y, double z);
set_position(100, 200, 300);  // Is this meters? Centimeters? Feet? Who knows!

void set_rotation(double roll, double pitch, double yaw);
set_rotation(90, 0, 0);  // Is this degrees? Radians? Compiler can't help!

// With strong types (SAFE):
void set_position(meter_t x, meter_t y, meter_t z);
set_position(100_m, 200_m, 300_m);  // Clear: meters
// set_position(100, 200, 300);  // ERROR: won't compile!

void set_rotation(radian_t roll, radian_t pitch, radian_t yaw);
set_rotation(to_radians(90.0_deg), 0.0_rad, 0.0_rad);  // Explicit conversion
// set_rotation(90.0_deg, 0.0_rad, 0.0_rad);  // ERROR: won't compile!
```

### Part B: `position_t` class

Create a `position_t` class that:
- Stores x, y, z coordinates privately in `std::array<meter_t, 3>`
- **Constructor ONLY accepts `meter_t` types** (not raw doubles!)
- Provides const and non-const accessor methods: `x()`, `y()`, `z()`
- Returns `meter_t const&` from const methods, `meter_t&` from non-const methods
- Implements vector operations: `operator+`, `operator-`, `operator*` (scalar)
- Implements `operator==` for equality comparison

**Important:** Strong type storage! By storing `std::array<meter_t, 3>` instead of `std::array<double, 3>`, the type system prevents unit confusion at every level.
**HINT:** Use `std::hypot` where appropriate.

```cpp
struct position_t {
    // Constructor enforces meter_t type
    position_t(meter_t x, meter_t y, meter_t z);

    // Const accessors - return const reference to meter_t (read-only)
    [[nodiscard]] meter_t const& x() const;
    [[nodiscard]] meter_t const& y() const;
    [[nodiscard]] meter_t const& z() const;

    // Non-const accessors - return mutable reference to meter_t (read-write)
    [[nodiscard]] meter_t& x();
    [[nodiscard]] meter_t& y();
    [[nodiscard]] meter_t& z();

    // Vector operations
    [[nodiscard]] position_t operator+(position_t const& other) const;
    [[nodiscard]] position_t operator-(position_t const& other) const;
    [[nodiscard]] position_t operator*(double scalar) const;

    // Equality comparison
    [[nodiscard]] bool operator==(position_t const& other) const;

private:
    std::array<meter_t, 3> coords_;  // Stored as array of meter_t types!
};
```

**Implement free functions for position operations:**

After defining the `position_t` class, implement these two free functions:

```cpp
// Calculate distance between two positions (returns strong type!)
// Hint: Use std::hypot(dx, dy, dz) instead of std::sqrt(dx*dx + dy*dy + dz*dz)
[[nodiscard]] meter_t distance(position_t const& p1, position_t const& p2);

// Check if two positions are approximately equal
[[nodiscard]] bool near(position_t const& p1, position_t const& p2,
                        meter_t tolerance = meter_t{0.001});
```

**Usage:**
```cpp
using namespace literals;

position_t const p1{1.5_m, 2.3_m, 0.0_m};  // OK: strong types
position_t const p2{3.0_m, 4.0_m, 5.0_m};
// position_t const p3{1.5, 2.3, 0.0};     // ERROR: no implicit conversion!

// Free functions for distance and comparison
meter_t const dist = distance(p1, p2);           // Type-safe distance
bool const are_close = near(p1, p2, 0.01_m);   // Approximate equality
```

### Part C: `quaternion_t` class

Create a `quaternion_t` class for rotations that:
- Stores x, y, z, w components privately in `std::array<double, 4>`
  - **w** is the **real (scalar) part**
  - **x, y, z** are the **imaginary (vector) parts**
  - Quaternions represent rotations in 3D space as: `q = w + xi + yj + zk`
- **Constructor accepts raw doubles** (quaternion components are unitless)
- Has named constructor: `static quaternion_t from_euler(radian_t roll, radian_t pitch, radian_t yaw)`
  - Note: **Accepts only `radian_t` type**, not `degree_t`!
- Normalizes quaternions in constructor
- Implements `operator*` for quaternion multiplication
  - **Quaternion multiplication represents composition of rotations**: multiplying quaternion `q1 * q2` applies rotation `q2` first, then `q1`
  - References for understanding quaternion multiplication:
    - [Math StackExchange: What does multiplication of two quaternions give?](https://math.stackexchange.com/questions/360286/what-does-multiplication-of-two-quaternions-give)
    - [MathWorks: Quaternion Multiplication](https://www.mathworks.com/help/aeroblks/quaternionmultiplication.html)
- Provides `conjugate()` method for inverse rotation

```cpp
struct quaternion_t {
    // Default constructor - identity rotation
    quaternion_t();

    // Quaternion components are unitless (normalized)
    quaternion_t(double x, double y, double z, double w);

    // Named constructor enforces radian_t type
    [[nodiscard]] static quaternion_t from_euler(radian_t roll, radian_t pitch, radian_t yaw);

    // Const accessors - return const reference (read-only)
    [[nodiscard]] double const& x() const;
    [[nodiscard]] double const& y() const;
    [[nodiscard]] double const& z() const;
    [[nodiscard]] double const& w() const;

    // Non-const accessors - return mutable reference (read-write)
    [[nodiscard]] double& x();
    [[nodiscard]] double& y();
    [[nodiscard]] double& z();
    [[nodiscard]] double& w();

    // Quaternion multiplication (composition of rotations)
    [[nodiscard]] quaternion_t operator*(quaternion_t const& other) const;

    // Arithmetic operators for component-wise operations (needed for averaging in Exercise 2)
    // Note: These perform simple component-wise operations, not geometric quaternion operations
    [[nodiscard]] quaternion_t operator+(quaternion_t const& other) const;
    [[nodiscard]] quaternion_t operator-(quaternion_t const& other) const;
    [[nodiscard]] quaternion_t operator/(double scalar) const;

    // Conjugate (inverse rotation for unit quaternions)
    [[nodiscard]] quaternion_t conjugate() const;

    // Equality comparison
    [[nodiscard]] bool operator==(quaternion_t const& other) const;

private:
    std::array<double, 4> components_;  // x, y, z, w
};
```

**Implement free functions for quaternion operations:**

After defining the `quaternion_t` class, implement this free function:

```cpp
// Check if two quaternions are approximately equal
// Hint: Compare each component (x, y, z, w) separately
[[nodiscard]] bool near(quaternion_t const& q1, quaternion_t const& q2,
                        double tolerance = 0.001);
```

**Usage:**
```cpp
using namespace literals;

// Must explicitly convert degrees to radians using the conversion function
auto const rot1 = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(90.0_deg));

// Or use radians directly
auto const rot2 = quaternion_t::from_euler(0.0_rad, 0.0_rad, 1.57_rad);

// Can also convert radians to degrees
auto const angle_deg = to_degrees(1.57_rad);

// This won't compile (type safety!):
// auto const rot3 = quaternion_t::from_euler(0.0, 0.0, 90.0);  // ERROR!
// auto const rot4 = quaternion_t::from_euler(0.0_deg, 0.0_deg, 90.0_deg);  // ERROR!
```

### Part D: `transformation_t` class

Create a `transformation_t` class that:
- Stores a full 4x4 transformation matrix as `std::array<double, 16>` internally
  - This represents a homogeneous transformation matrix (position + rotation combined)
  - Storage is an implementation detail - students compute this from position_t + quaternion_t
- **Can ONLY be constructed** with `transformation_t(position_t const&, quaternion_t const&)` (no default constructor)
  - Constructor converts position_t + quaternion_t into the 4x4 matrix representation
- Provides `position()` and `rotation()` methods
  - These extract position_t/quaternion_t from the internal matrix representation
  - Only these types are exposed in the public interface
- Implements `operator*` for transformation composition (matrix multiplication)
- Implements `operator*` for transforming positions (applies rotation and translation)
- Implements `inverse()` method

```cpp
struct transformation_t {
    // Constructor - requires both position and rotation (no default constructor!)
    transformation_t(position_t const& pos, quaternion_t const& rot);

    // Extract position from transformation matrix
    [[nodiscard]] position_t position() const;

    // Extract rotation quaternion from transformation matrix
    [[nodiscard]] quaternion_t rotation() const;

    // Transformation composition (matrix multiplication)
    [[nodiscard]] transformation_t operator*(transformation_t const& other) const;

    // Apply additional rotation to this transformation
    [[nodiscard]] transformation_t operator*(quaternion_t const& rotation) const;

    // Transform a point by this transformation (applies rotation and translation)
    [[nodiscard]] position_t operator*(position_t const& point) const;

    // Compute inverse transformation
    [[nodiscard]] transformation_t inverse() const;

    // Equality comparison
    [[nodiscard]] bool operator==(transformation_t const& other) const;

private:
    std::array<double, 16> matrix_;  // 4x4 homogeneous transformation matrix (row-major)
};
```

**Implement free functions for transformation operations:**

After defining the `transformation_t` class, implement this free function:

```cpp
// Check if two transformations are approximately equal
// Hint: Use the near() functions you implemented for position_t and quaternion_t
[[nodiscard]] bool near(transformation_t const& tf1, transformation_t const& tf2,
                        meter_t pos_tolerance = meter_t{0.001},
                        double rot_tolerance = 0.001);
```

**Complete Example:**
```cpp
using namespace literals;

// Strong types prevent unit errors at compile time!
position_t const pos{1.0_m, 2.0_m, 3.0_m};  // Type-safe construction

// Must explicitly handle degree -> radian conversion
quaternion_t const rot = quaternion_t::from_euler(
    0.0_rad,
    0.0_rad,
    to_radians(90.0_deg)  // Explicit conversion required
);

// Create transformation (no default constructor allowed)
transformation_t const tf{pos, rot};

// Compose transformations
transformation_t const tf2{position_t{0.5_m, 0.0_m, 0.0_m}, quaternion_t{}};
auto const composed = tf * tf2;

// Transform a point using operator*
position_t const point{0.0_m, 0.0_m, 0.0_m};
auto const transformed = tf * point;

// Extract position and rotation from transformation
auto const extracted_pos = tf.position();
auto const extracted_rot = tf.rotation();

// Type-safe distance calculation (free function)
meter_t const dist = distance(pos, point);
std::cout << std::format("Distance: {} meters\n", dist.value);

// These won't compile (type safety in action!):
// position_t const bad1{1.0, 2.0, 3.0};  // ERROR: no implicit conversion from double
// auto const bad2 = quaternion_t::from_euler(0.0, 0.0, 90.0);  // ERROR: needs radian_t type
// auto const bad3 = quaternion_t::from_euler(0.0_deg, 0.0_deg, 90.0_deg);  // ERROR: needs radian_t, not degree_t

// But these conversions work:
auto const angle_in_deg = to_degrees(1.57_rad);  // Convert radian_t to degree_t
auto const angle_in_rad = to_radians(90.0_deg);  // Convert degree_t to radian_t
```

**Test in Godbolt:**
https://godbolt.org/z/hz7Mq4bG1
```cpp
#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>

// Part A: Strong Types with User-Defined Literals
namespace literals {
    struct meter_t {
        double value;
        constexpr explicit meter_t(double v) : value{v} {}

        [[nodiscard]] constexpr meter_t operator+(meter_t const& other) const {
            return meter_t{value + other.value};
        }
        [[nodiscard]] constexpr meter_t operator-(meter_t const& other) const {
            return meter_t{value - other.value};
        }
        [[nodiscard]] constexpr meter_t operator*(double scalar) const {
            return meter_t{value * scalar};
        }
        [[nodiscard]] constexpr meter_t operator/(double scalar) const {
            return meter_t{value / scalar};
        }
        [[nodiscard]] constexpr auto operator<=>(meter_t const& other) const = default;
    };

    struct radian_t {
        double value;
        constexpr explicit radian_t(double v) : value{v} {}

        [[nodiscard]] constexpr radian_t operator+(radian_t const& other) const {
            return radian_t{value + other.value};
        }
        [[nodiscard]] constexpr radian_t operator-(radian_t const& other) const {
            return radian_t{value - other.value};
        }
        [[nodiscard]] constexpr radian_t operator*(double scalar) const {
            return radian_t{value * scalar};
        }
        [[nodiscard]] constexpr radian_t operator/(double scalar) const {
            return radian_t{value / scalar};
        }
        [[nodiscard]] constexpr auto operator<=>(radian_t const& other) const = default;
    };

    struct degree_t {
        double value;
        constexpr explicit degree_t(double deg) : value{deg} {}

        [[nodiscard]] constexpr degree_t operator+(degree_t const& other) const {
            return degree_t{value + other.value};
        }
        [[nodiscard]] constexpr degree_t operator-(degree_t const& other) const {
            return degree_t{value - other.value};
        }
        [[nodiscard]] constexpr degree_t operator*(double scalar) const {
            return degree_t{value * scalar};
        }
        [[nodiscard]] constexpr degree_t operator/(double scalar) const {
            return degree_t{value / scalar};
        }
        [[nodiscard]] constexpr auto operator<=>(degree_t const& other) const = default;
    };

    // Conversion functions
    constexpr radian_t to_radians(degree_t const d) {
        return radian_t{d.value * std::numbers::pi / 180.0};
    }
    constexpr degree_t to_degrees(radian_t const r) {
        return degree_t{r.value * 180.0 / std::numbers::pi};
    }

    // User-defined literals
    constexpr meter_t operator""_m(long double v) { return meter_t{static_cast<double>(v)}; }
    constexpr meter_t operator""_m(unsigned long long v) { return meter_t{static_cast<double>(v)}; }
    constexpr radian_t operator""_rad(long double v) { return radian_t{static_cast<double>(v)}; }
    constexpr radian_t operator""_rad(unsigned long long v) { return radian_t{static_cast<double>(v)}; }
    constexpr degree_t operator""_deg(long double v) { return degree_t{static_cast<double>(v)}; }
    constexpr degree_t operator""_deg(unsigned long long v) { return degree_t{static_cast<double>(v)}; }
}

using literals::meter_t;
using literals::radian_t;

// Part B: position class
struct position_t {
    // TODO: Implement all methods
    position_t(meter_t x, meter_t y, meter_t z);
    [[nodiscard]] meter_t const& x() const;
    [[nodiscard]] meter_t const& y() const;
    [[nodiscard]] meter_t const& z() const;
    [[nodiscard]] meter_t& x();
    [[nodiscard]] meter_t& y();
    [[nodiscard]] meter_t& z();
    [[nodiscard]] position_t operator+(position_t const& other) const;
    [[nodiscard]] position_t operator-(position_t const& other) const;
    [[nodiscard]] position_t operator*(double scalar) const;
    [[nodiscard]] bool operator==(position_t const& other) const;
private:
    std::array<meter_t, 3> coords_;
};

// Free functions for position
[[nodiscard]] meter_t distance(position_t const& p1, position_t const& p2);
[[nodiscard]] bool near(position_t const& p1, position_t const& p2,
                        meter_t tolerance = meter_t{0.001});

// Part C: quaternion class
struct quaternion_t {
    // TODO: Implement all methods
    quaternion_t();
    quaternion_t(double x, double y, double z, double w);
    [[nodiscard]] static quaternion_t from_euler(radian_t roll, radian_t pitch, radian_t yaw);
    [[nodiscard]] double const& x() const;
    [[nodiscard]] double const& y() const;
    [[nodiscard]] double const& z() const;
    [[nodiscard]] double const& w() const;
    [[nodiscard]] double& x();
    [[nodiscard]] double& y();
    [[nodiscard]] double& z();
    [[nodiscard]] double& w();
    [[nodiscard]] quaternion_t operator*(quaternion_t const& other) const;
    [[nodiscard]] quaternion_t operator+(quaternion_t const& other) const;
    [[nodiscard]] quaternion_t operator-(quaternion_t const& other) const;
    [[nodiscard]] quaternion_t operator/(double scalar) const;
    [[nodiscard]] quaternion_t conjugate() const;
    [[nodiscard]] bool operator==(quaternion_t const& other) const;
private:
    std::array<double, 4> components_;
};

// Free functions for quaternion
[[nodiscard]] bool near(quaternion_t const& q1, quaternion_t const& q2,
                        double tolerance = 0.001);

// Part D: transformation class
struct transformation_t {
    // TODO: Implement all methods
    transformation_t(position_t const& pos, quaternion_t const& rot);
    [[nodiscard]] position_t position() const;
    [[nodiscard]] quaternion_t rotation() const;
    [[nodiscard]] transformation_t operator*(transformation_t const& other) const;
    [[nodiscard]] transformation_t operator*(quaternion_t const& rotation) const;
    [[nodiscard]] position_t operator*(position_t const& point) const;
    [[nodiscard]] transformation_t inverse() const;
    [[nodiscard]] bool operator==(transformation_t const& other) const;
private:
    std::array<double, 16> matrix_;
};

// Free functions for transformation
[[nodiscard]] bool near(transformation_t const& tf1, transformation_t const& tf2,
                        meter_t pos_tolerance = meter_t{0.001},
                        double rot_tolerance = 0.001);

int main() {
    using namespace literals;

    // Test strong types
    auto const m1 = 5.0_m;
    auto const m2 = 3.0_m;
    assert(m1.value == 5.0);
    assert(m2.value == 3.0);

    // Test angle conversions
    auto const deg = 90.0_deg;
    auto const rad = to_radians(deg);
    assert(std::abs(rad.value - std::numbers::pi / 2.0) < 0.001);

    auto const rad2 = 3.14159_rad;
    auto const deg2 = to_degrees(rad2);
    assert(std::abs(deg2.value - 180.0) < 0.1);

    // Test position class
    auto const p1 = position_t{1.0_m, 2.0_m, 3.0_m};
    auto const p2 = position_t{4.0_m, 5.0_m, 6.0_m};
    assert(p1.x().value == 1.0);
    assert(p1.y().value == 2.0);
    assert(p1.z().value == 3.0);

    // Test position operations
    auto const p3 = p1 + p2;
    assert(p3.x().value == 5.0);
    assert(p3.y().value == 7.0);
    assert(p3.z().value == 9.0);

    auto const p4 = p2 - p1;
    assert(p4.x().value == 3.0);
    assert(p4.y().value == 3.0);
    assert(p4.z().value == 3.0);

    auto const p5 = p1 * 2.0;
    assert(p5.x().value == 2.0);
    assert(p5.y().value == 4.0);
    assert(p5.z().value == 6.0);

    // Test distance and near
    auto const origin = position_t{0.0_m, 0.0_m, 0.0_m};
    auto const unit_x = position_t{1.0_m, 0.0_m, 0.0_m};
    auto const dist = distance(origin, unit_x);
    assert(std::abs(dist.value - 1.0) < 0.001);

    auto const p6 = position_t{1.001_m, 0.0_m, 0.0_m};
    assert(near(unit_x, p6, 0.01_m));
    assert(!near(unit_x, p6, 0.0001_m));

    // Test quaternion class
    quaternion_t const q_identity;
    assert(q_identity.w() == 1.0);
    assert(q_identity.x() == 0.0);
    assert(q_identity.y() == 0.0);
    assert(q_identity.z() == 0.0);

    // Test quaternion from Euler angles
    auto const q_90z = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(90.0_deg));
    assert(std::abs(q_90z.w() - 0.707) < 0.01);  // cos(45°)
    assert(std::abs(q_90z.z() - 0.707) < 0.01);  // sin(45°)
    assert(std::abs(q_90z.x()) < 0.001);
    assert(std::abs(q_90z.y()) < 0.001);

    // Test quaternion multiplication
    auto const q1 = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg));
    auto const q2 = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg));
    auto const q_combined = q1 * q2;
    assert(near(q_combined, q_90z, 0.01));

    // Test quaternion conjugate
    auto const q_conj = q_90z.conjugate();
    auto const q_result = q_90z * q_conj;
    assert(near(q_result, q_identity, 0.001));

    // Test transformation class
    transformation_t const tf_identity{origin, q_identity};
    auto const extracted_pos = tf_identity.position();
    assert(extracted_pos == origin);
    auto const extracted_rot = tf_identity.rotation();
    assert(near(extracted_rot, q_identity, 0.001));

    // Test transformation composition
    transformation_t const tf1{position_t{1.0_m, 0.0_m, 0.0_m}, q_identity};
    transformation_t const tf2{position_t{0.0_m, 1.0_m, 0.0_m}, q_identity};
    auto const tf_composed = tf1 * tf2;
    auto const composed_pos = tf_composed.position();
    assert(near(composed_pos, position_t{1.0_m, 1.0_m, 0.0_m}, 0.001_m));

    // Test transforming a point
    transformation_t const tf_translate{position_t{1.0_m, 2.0_m, 3.0_m}, q_identity};
    position_t const point{0.0_m, 0.0_m, 0.0_m};
    auto const transformed = tf_translate * point;
    assert(near(transformed, position_t{1.0_m, 2.0_m, 3.0_m}, 0.001_m));

    // Test rotation transformation
    transformation_t const tf_rotate{origin, q_90z};
    position_t const point_x{1.0_m, 0.0_m, 0.0_m};
    auto const rotated = tf_rotate * point_x;
    // 90° rotation around Z should map (1,0,0) to approximately (0,1,0)
    assert(std::abs(rotated.x().value) < 0.01);
    assert(std::abs(rotated.y().value - 1.0) < 0.01);
    assert(std::abs(rotated.z().value) < 0.01);

    // Test transformation inverse
    auto const tf_inv = tf_translate.inverse();
    auto const back = tf_inv * transformed;
    assert(near(back, point, 0.001_m));

    // Critical test: quaternion multiplication should match transformation multiplication
    // Two quaternions multiplied should give same rotation as two transformations
    auto const qa = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(30.0_deg));
    auto const qb = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(60.0_deg));
    auto const q_mult = qa * qb;  // Should be 90° rotation

    transformation_t const tfa{origin, qa};
    transformation_t const tfb{origin, qb};
    auto const tf_mult = tfa * tfb;
    auto const q_from_tf = tf_mult.rotation();

    assert(near(q_mult, q_from_tf, 0.001));

    // Test transformation * quaternion operator
    transformation_t const tf_base{position_t{1.0_m, 0.0_m, 0.0_m}, q_identity};
    auto const q_rotate = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(90.0_deg));
    auto const tf_rotated = tf_base * q_rotate;

    // Should preserve position, update rotation
    auto const rotated_pos = tf_rotated.position();
    assert(near(rotated_pos, position_t{1.0_m, 0.0_m, 0.0_m}, 0.001_m));
    auto const rotated_quat = tf_rotated.rotation();
    assert(near(rotated_quat, q_rotate, 0.001));

    std::cout << "All tests passed!\n";
    return 0;
}
```

**Learning objectives:**
1. **Strong Types** (⭐ Primary Focus)
   - Creating type-safe wrappers around primitive types
   - Using `explicit` constructors to prevent implicit conversion
   - Compile-time unit safety (preventing meter/centimeter or degree/radian mixups)

2. **User-Defined Literals**
   - Implementing `operator""_m`, `operator""_deg`, `operator""_rad`
   - Both `long double` and `unsigned long long` overloads
   - Making code self-documenting and type-safe

3. **Encapsulation with `std::array`**
   - Private storage with public accessor methods
   - Information hiding and interface design

4. **Const-Correctness**
   - Dual accessor overloads: `double const& x() const` and `double& x()`
   - Providing read-only and read-write access safely

5. **Advanced Construction Patterns**
   - Named constructors: `quaternion_t::from_euler()`
   - Preventing default construction when semantically meaningless
   - Constructor invariants (quaternion normalization)

6. **Operator Overloading**
   - Mathematical operators: `+`, `-`, `*`
   - Comparison operators: `==`

7. **Free Functions** (⭐ Important!)
   - `distance(position_t, position_t)` - Calculate distance between positions
   - `near(position_t, position_t, meter_t)` - Approximate equality comparison
   - Free functions provide symmetry and better composability

8. **Modern C++ Features**
   - `[[nodiscard]]` attribute
   - `noexcept` on move operations
   - `constexpr` for compile-time evaluation
   - East const style

---

## Exercise 2: Template Specialization for Rotation Averaging (Medium-Hard)

When averaging rotations (quaternions or transformation matrices), simple linear averaging doesn't work correctly. This exercise teaches **template specialization** by implementing type-appropriate averaging algorithms.

### Part A: Primary Template (Linear Average)

Implement a primary template that computes the arithmetic mean for any type `T`:

```cpp
// Primary template: linear average for numeric types
template<typename T>
T average(std::span<T const> values);
```

**Requirements:**
- Return default-constructed `T{}` for empty spans
- Works with `int`, `float`, `double`, etc.
- Sum all elements, divide by count

### Part B: Quaternion Specialization (Eigenvector Method)

For quaternions, linear averaging produces incorrect rotations. The geometrically correct approach uses the **eigenvector method**.

```cpp
// Full specialization for quaternion_t
template<>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> quaternions);
```

#### Mathematical Background

Quaternions represent rotations on the unit 3-sphere (S³). The "average" rotation is the quaternion that minimizes the sum of squared geodesic distances to all input quaternions. This is called the **Fréchet mean**.

**Algorithm (Markley et al., 2007):**

1. **Build the 4×4 accumulator matrix** M:

```
M = Σᵢ qᵢ ⊗ qᵢᵀ
```

Where qᵢ = [xᵢ, yᵢ, zᵢ, wᵢ]ᵀ and ⊗ denotes outer product.

Expanded form:
```
    ⎡ Σxᵢxᵢ  Σxᵢyᵢ  Σxᵢzᵢ  Σxᵢwᵢ ⎤
M = ⎢ Σyᵢxᵢ  Σyᵢyᵢ  Σyᵢzᵢ  Σyᵢwᵢ ⎥
    ⎢ Σzᵢxᵢ  Σzᵢyᵢ  Σzᵢzᵢ  Σzᵢwᵢ ⎥
    ⎣ Σwᵢxᵢ  Σwᵢyᵢ  Σwᵢzᵢ  Σwᵢwᵢ ⎦
```

2. **Find the eigenvector corresponding to the largest eigenvalue** of M.

3. **That eigenvector IS the average quaternion** (normalize it to unit length).

**Why this works:** The matrix M encodes the "spread" of quaternions. Its largest eigenvector points in the direction of maximum agreement among all input quaternions.

#### Power Iteration Algorithm

To find the largest eigenvector without external libraries, use **power iteration**:

```
Input: 4×4 symmetric matrix M
Output: eigenvector corresponding to largest eigenvalue

1. Initialize v = [0, 0, 0, 1]ᵀ  (or any non-zero vector)
2. Repeat until convergence (e.g., 20 iterations):
   a. v_new = M × v           (matrix-vector multiply)
   b. v = v_new / ‖v_new‖     (normalize)
3. Return v as the average quaternion
```

**Convergence:** Power iteration converges when the largest eigenvalue is unique and dominant. For rotation averaging, this is virtually always the case.

**Implementation hint:** The matrix M is symmetric, so M[i][j] == M[j][i].

#### References

- Markley, F.L., Cheng, Y., Crassidis, J.L., Oshman, Y. (2007). *"Averaging Quaternions"*, Journal of Guidance, Control, and Dynamics, 30(4), 1193-1197. NASA Technical Report: https://ntrs.nasa.gov/citations/20070017872


### Part C: Transformation Specialization (Matrix Projection Method)

For transformations, we must average position and rotation separately—but **without converting to quaternions**.

```cpp
// Full specialization for transformation_t
template<>
transformation_t average<transformation_t>(std::span<transformation_t const> transforms);
```

#### Mathematical Background

A transformation matrix combines rotation (3×3 orthogonal matrix R) and translation (3×1 vector t):

```
T = ⎡ R  | t ⎤
    ⎣ 0  | 1 ⎦
```

**Algorithm:**

1. **Average the translations linearly:**
```
t_avg = (1/n) Σᵢ tᵢ
```

2. **Average the rotation matrices element-wise:**
```
R_mean = (1/n) Σᵢ Rᵢ
```

Note: R_mean is generally **not** a valid rotation matrix (not orthogonal).

3. **Project R_mean back to SO(3)** using polar decomposition:
```
R_avg = R_mean × (R_meanᵀ × R_mean)^(-½)
```

This finds the closest valid rotation matrix to R_mean.

#### Computing (AᵀA)^(-½) via Denman-Beavers Iteration

The matrix square root inverse can be computed iteratively without eigendecomposition:

```
Input: 3×3 symmetric positive-definite matrix A = RᵀR
Output: A^(-½)

1. Initialize:
   Y₀ = A
   Z₀ = I (3×3 identity)

2. Repeat until convergence (e.g., 20 iterations):
   Y_{k+1} = ½(Yₖ + Zₖ⁻¹)
   Z_{k+1} = ½(Zₖ + Yₖ⁻¹)

3. Return Zₖ as A^(-½)
```

**Implementation hints:**
- You need a 3×3 matrix inverse function
- Check convergence by comparing ‖Y_{k+1} - Yₖ‖ < ε
- For robustness, limit to fixed iterations (20 is sufficient)

#### Why Not Use Quaternions?

The test suite validates that quaternion averaging and matrix averaging produce **identical results**. By implementing both independently, we verify correctness of each algorithm. If one used the other internally, the test would be circular.

#### References

- Moakher, M. (2002). *"Means and Averaging in the Group of Rotations"*, SIAM Journal on Matrix Analysis and Applications.
- Higham, N.J. (1986). *"Computing the Polar Decomposition—with Applications"*, SIAM Journal on Scientific and Statistical Computing.


### Template Specialization Syntax

**Primary template:**
```cpp
template<typename T>
T average(std::span<T const> values) {
    // Default implementation for numeric types
}
```

**Full specialization:**
```cpp
template<>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> quaternions) {
    // Quaternion-specific implementation
}

template<>
transformation_t average<transformation_t>(std::span<transformation_t const> transforms) {
    // Transformation-specific implementation
}
```

**Key points:**
- `template<>` with empty angle brackets signals full specialization
- The function signature must match the primary template exactly
- Specializations are selected at compile time based on the type argument


### Starter Code

```cpp
#include <array>
#include <cmath>
#include <span>

// Assume quaternion_t and transformation_t from Exercise 1 are available

// Primary template: linear average
template<typename T>
T average(std::span<T const> values) {
    // TODO: Implement linear average
    // Handle empty span
    // Sum all elements, divide by count
}

// Full specialization for quaternion_t
template<>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> quaternions) {
    // TODO: Implement eigenvector method
    // 1. Build 4x4 matrix M = Σ qᵢqᵢᵀ
    // 2. Power iteration to find largest eigenvector
    // 3. Normalize and return
}

// Full specialization for transformation_t
template<>
transformation_t average<transformation_t>(std::span<transformation_t const> transforms) {
    // TODO: Implement matrix projection method
    // 1. Linear average of positions
    // 2. Element-wise average of rotation matrices
    // 3. Polar decomposition to project back to SO(3)
    // 4. Combine into transformation_t
}
```


### Validation Tests

Your implementation should pass these correctness tests:

**Test 1: Pure Yaw Rotations**
```cpp
// Create rotations with only yaw (rotation around Z axis)
std::array<double, 3> yaw_angles = {0.1, 0.2, 0.3};  // radians

// Method 1: Average angles linearly
double avg_angle = average<double>(yaw_angles);  // = 0.2

// Method 2: Convert to quaternions, average, extract yaw
std::array<quaternion_t, 3> quats = /* from yaw angles */;
quaternion_t avg_quat = average<quaternion_t>(quats);
double quat_yaw = /* extract yaw from avg_quat */;

// Method 3: Convert to transforms, average, extract yaw
std::array<transformation_t, 3> transforms = /* from yaw angles */;
transformation_t avg_tf = average<transformation_t>(transforms);
double tf_yaw = /* extract yaw from avg_tf */;

// All three should match!
assert(std::abs(avg_angle - quat_yaw) < 0.001);
assert(std::abs(avg_angle - tf_yaw) < 0.001);
```

**Test 2: Combined Pitch and Roll**
```cpp
// Create rotations with both pitch and roll
std::array<quaternion_t, 3> quats = {
    quaternion_t::from_euler(0.1_rad, 0.2_rad, 0.0_rad),
    quaternion_t::from_euler(0.2_rad, 0.1_rad, 0.0_rad),
    quaternion_t::from_euler(0.15_rad, 0.15_rad, 0.0_rad),
};

std::array<transformation_t, 3> transforms = /* from same angles */;

quaternion_t avg_quat = average<quaternion_t>(quats);
transformation_t avg_tf = average<transformation_t>(transforms);

// Extract rotation from transform and compare to quaternion
quaternion_t tf_rotation = avg_tf.rotation();

// Should match (within numerical tolerance)
assert(near(avg_quat, tf_rotation, 0.001));
```


### Exercise 2 Learning Objectives

After completing this exercise, you should understand:

1. **Template specialization**
   - Primary template vs full specialization syntax
   - When specialization is selected (compile-time type matching)
   - Why specialization is needed for type-specific algorithms

2. **Rotation averaging mathematics**
   - Why linear averaging fails for rotations
   - Eigenvector method for quaternions
   - Polar decomposition for rotation matrices
   - Fréchet mean on manifolds (conceptual)

3. **Numerical algorithms**
   - Power iteration for dominant eigenvector
   - Denman-Beavers iteration for matrix square root
   - Convergence criteria and iteration limits

4. **`std::span` usage**
   - Non-owning views of contiguous data
   - Const-correctness with `std::span<T const>`

---

## Learning Objectives

After completing these exercises, you should understand:

1. **`std::array` fundamentals**:
   - Fixed-size, stack-allocated containers
   - Compile-time size checking
   - STL interface compatibility (iterators, algorithms)
   - When to prefer over C-style arrays and `std::vector`

2. **`std::span` fundamentals**:
   - Non-owning views of contiguous memory
   - Working with multiple container types generically
   - Dynamic vs fixed extent spans
   - Const-correctness with spans

3. **Template specialization**:
   - Primary template vs full specialization
   - Compile-time type dispatch
   - When and why to specialize

4. **Numerical algorithms**:
   - Power iteration for eigenvectors
   - Polar decomposition for rotation matrices
   - Iterative methods without external libraries

5. **Modern C++ patterns**:
   - Strong types for unit safety
   - User-defined literals
   - Compile-time size guarantees

## Discussion Questions

1. When should you use `std::array` instead of `std::vector`?
2. What are the lifetime considerations when using `std::span`?
3. Why can't you use linear averaging for rotations?
4. When would you use partial specialization vs full specialization?
5. How does power iteration find the dominant eigenvector?
6. Why do the quaternion and matrix methods produce the same result?

---

## Solutions

✅ **[View Solutions with Detailed Explanations](week3_solutions.md)**
