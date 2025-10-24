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
    };

    // Strong type for angles in radians
    struct radian_t {
        double value;
        constexpr explicit radian_t(double v) : value{v} {}
    };

    // Strong type for angles in degrees
    struct degree_t {
        double value;
        constexpr explicit degree_t(double deg) : value{deg} {}
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
        [[nodiscard]] constexpr bool operator==(meter_t const& other) const = default;
    };

    struct radian_t {
        double value;
        constexpr explicit radian_t(double v) : value{v} {}
        [[nodiscard]] constexpr bool operator==(radian_t const& other) const = default;
    };

    struct degree_t {
        double value;
        constexpr explicit degree_t(double deg) : value{deg} {}
        [[nodiscard]] constexpr bool operator==(degree_t const& other) const = default;
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

5. **Rule of 5**
   - Copy constructor, move constructor, copy/move assignment, destructor
   - When to use `= default` vs explicit implementation

6. **Advanced Construction Patterns**
   - Named constructors: `quaternion::from_euler()`
   - Preventing default construction when semantically meaningless
   - Constructor invariants (quaternion normalization)

7. **Operator Overloading**
   - Mathematical operators: `+`, `-`, `*`
   - Comparison operators: `==`

8. **Free Functions** (⭐ Important!)
   - `distance(position, position)` - Calculate distance between positions
   - `near(position, position, meter)` - Approximate equality comparison
   - Free functions provide symmetry and better composability

9. **Modern C++ Features**
   - `[[nodiscard]]` attribute
   - `noexcept` on move operations
   - `constexpr` for compile-time evaluation
   - East const style

---

## Exercise 2: Sensor Array Statistics (Medium)

Your robot has multiple sensors (LIDAR, cameras, IMU, etc.). You need generic functions that can compute statistics over any contiguous sequence of sensor readings.

Write functions using `std::span` to work with sensor data:

```cpp
// Calculate the average of sensor readings
double average(std::span<double const> readings);

// Find the minimum value in sensor readings
double min_value(std::span<double const> readings);

// Find the maximum value in sensor readings
double max_value(std::span<double const> readings);

// Count how many readings exceed a threshold
size_t count_above_threshold(std::span<double const> readings, double threshold);

// Normalize readings to [0, 1] range (modifies in-place)
void normalize(std::span<double> readings);

// Sliding window average - computes average of each window of size 'window_size'
// Returns a vector with (values.size() - window_size + 1) elements
// Template this function to work with any numeric type!
template<typename T>
std::vector<T> sliding_window_average(std::span<T const> values, size_t window_size);
```

**Example:**
```cpp
using namespace literals;

// Full 360-degree LIDAR scan (one reading per degree)
std::array<meter, 360> lidar_scan = { /* ... 360 distance readings ... */ };

// Use span to analyze just the front 90 degrees (indices 0-89)
auto const front_readings = std::span{lidar_scan}.subspan(0, 90);

// Extract values from strong types for averaging
std::array<double, 90> front_values;
for (size_t i = 0; i < front_readings.size(); ++i) {
    front_values[i] = front_readings[i].value;
}
double const front_avg = average(front_values);

// Sliding window example: smooth noisy IMU orientation data
// Note: quaternion_t needs operator+, operator-, and operator/ for averaging
std::array<quaternion_t, 20> imu_orientations = {
    quaternion_t::from_euler(0.01_rad, 0.02_rad, 0.10_rad),
    quaternion_t::from_euler(0.02_rad, 0.01_rad, 0.15_rad),
    quaternion_t::from_euler(0.01_rad, 0.03_rad, 0.12_rad),
    // ... 17 more noisy quaternion readings with roll, pitch, and yaw variations
};

// Apply 5-point moving average filter directly to quaternions!
auto const smoothed_orientations = sliding_window_average(
    std::span<quaternion_t const>{imu_orientations}, 5);
// smoothed_orientations.size() == 16 (20 - 5 + 1)
```

**Test in Godbolt:**
https://godbolt.org/z/YOUR_LINK_HERE
```cpp
#include <array>
#include <vector>
#include <span>
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numbers>

// Strong types and quaternion from Exercise 1 (simplified for this test)
namespace literals {
    struct radian_t {
        double value;
        constexpr explicit radian_t(double v) : value{v} {}
    };

    constexpr radian_t operator""_rad(long double v) {
        return radian_t{static_cast<double>(v)};
    }
    constexpr radian_t operator""_rad(unsigned long long v) {
        return radian_t{static_cast<double>(v)};
    }
}

struct quaternion_t {
    quaternion_t() : components_{0.0, 0.0, 0.0, 1.0} {}
    quaternion_t(double x, double y, double z, double w) : components_{x, y, z, w} {}

    static quaternion_t from_euler(literals::radian_t roll, literals::radian_t pitch, literals::radian_t yaw) {
        double const cy = std::cos(yaw.value * 0.5);
        double const sy = std::sin(yaw.value * 0.5);
        double const cp = std::cos(pitch.value * 0.5);
        double const sp = std::sin(pitch.value * 0.5);
        double const cr = std::cos(roll.value * 0.5);
        double const sr = std::sin(roll.value * 0.5);
        return quaternion_t{
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
            cr * cp * cy + sr * sp * sy
        };
    }

    double const& x() const { return components_[0]; }
    double const& y() const { return components_[1]; }
    double const& z() const { return components_[2]; }
    double const& w() const { return components_[3]; }

    // Arithmetic operators needed for averaging
    quaternion_t operator+(quaternion_t const& other) const {
        return quaternion_t{
            x() + other.x(),
            y() + other.y(),
            z() + other.z(),
            w() + other.w()
        };
    }

    quaternion_t operator-(quaternion_t const& other) const {
        return quaternion_t{
            x() - other.x(),
            y() - other.y(),
            z() - other.z(),
            w() - other.w()
        };
    }

    quaternion_t operator/(double scalar) const {
        return quaternion_t{
            x() / scalar,
            y() / scalar,
            z() / scalar,
            w() / scalar
        };
    }

private:
    std::array<double, 4> components_;
};

// Your solutions here
double average(std::span<double const> readings);
double min_value(std::span<double const> readings);
double max_value(std::span<double const> readings);
size_t count_above_threshold(std::span<double const> readings, double threshold);
void normalize(std::span<double> readings);

template<typename T>
std::vector<T> sliding_window_average(std::span<T const> values, size_t window_size);

int main() {
    // Test basic statistics
    std::array<double, 5> const data1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    assert(std::abs(average(data1) - 3.0) < 0.001);
    assert(std::abs(min_value(data1) - 1.0) < 0.001);
    assert(std::abs(max_value(data1) - 5.0) < 0.001);
    assert(count_above_threshold(data1, 3.0) == 2);

    // Test normalization
    std::vector<double> data2 = {2.0, 4.0, 6.0, 8.0};
    normalize(data2);
    assert(std::abs(data2[0] - 0.0) < 0.001);  // min -> 0
    assert(std::abs(data2[3] - 1.0) < 0.001);  // max -> 1

    // Test sliding window average
    std::array<double, 5> const data3 = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto const windowed = sliding_window_average(std::span<double const>{data3}, 3);
    assert(windowed.size() == 3);  // 5 - 3 + 1 = 3
    assert(std::abs(windowed[0] - 2.0) < 0.001);  // (1+2+3)/3
    assert(std::abs(windowed[1] - 3.0) < 0.001);  // (2+3+4)/3
    assert(std::abs(windowed[2] - 4.0) < 0.001);  // (3+4+5)/3

    // Test sliding window with integers
    std::vector<int> const int_data = {10, 20, 30, 40};
    auto const int_windowed = sliding_window_average(std::span<int const>{int_data}, 2);
    assert(int_windowed.size() == 3);
    assert(int_windowed[0] == 15);  // (10+20)/2
    assert(int_windowed[1] == 25);  // (20+30)/2
    assert(int_windowed[2] == 35);  // (30+40)/2

    // ROBOTICS EXAMPLE 1: LIDAR scan analysis with span subranges
    // Simulate a 360-degree LIDAR scan (simplified to just front 10 degrees for testing)
    std::array<double, 10> lidar_distances = {
        2.5, 2.6, 2.4, 2.7, 2.5,  // Front-left: ~2.5m
        2.5, 2.6, 2.5, 2.4, 2.6   // Front-right: ~2.5m
    };

    // Analyze just the front-left section using span (indices 0-4)
    auto const front_left = std::span{lidar_distances}.subspan(0, 5);
    double const front_left_avg = average(front_left);
    assert(std::abs(front_left_avg - 2.54) < 0.01);  // Average of front-left readings

    // Find minimum obstacle distance in front-right section (indices 5-9)
    auto const front_right = std::span{lidar_distances}.subspan(5, 5);
    double const closest_obstacle = min_value(front_right);
    assert(std::abs(closest_obstacle - 2.4) < 0.001);

    // ROBOTICS EXAMPLE 2: IMU orientation filtering with sliding window
    // Noisy IMU orientation readings (20 samples at 100Hz = 0.2 seconds of data)
    // Simulating a robot with noisy roll, pitch, and yaw measurements
    using namespace literals;
    std::array<quaternion_t, 20> const imu_orientations = {
        quaternion_t::from_euler(0.01_rad, 0.02_rad, 0.10_rad),
        quaternion_t::from_euler(0.02_rad, 0.01_rad, 0.15_rad),
        quaternion_t::from_euler(0.01_rad, 0.03_rad, 0.12_rad),
        quaternion_t::from_euler(0.03_rad, 0.02_rad, 0.18_rad),
        quaternion_t::from_euler(0.02_rad, 0.01_rad, 0.14_rad),
        quaternion_t::from_euler(0.01_rad, 0.04_rad, 0.20_rad),
        quaternion_t::from_euler(0.03_rad, 0.03_rad, 0.16_rad),
        quaternion_t::from_euler(0.02_rad, 0.02_rad, 0.22_rad),
        quaternion_t::from_euler(0.01_rad, 0.01_rad, 0.18_rad),
        quaternion_t::from_euler(0.04_rad, 0.03_rad, 0.24_rad),
        quaternion_t::from_euler(0.02_rad, 0.02_rad, 0.20_rad),
        quaternion_t::from_euler(0.03_rad, 0.01_rad, 0.26_rad),
        quaternion_t::from_euler(0.01_rad, 0.04_rad, 0.22_rad),
        quaternion_t::from_euler(0.02_rad, 0.03_rad, 0.28_rad),
        quaternion_t::from_euler(0.03_rad, 0.02_rad, 0.24_rad),
        quaternion_t::from_euler(0.01_rad, 0.01_rad, 0.30_rad),
        quaternion_t::from_euler(0.04_rad, 0.02_rad, 0.26_rad),
        quaternion_t::from_euler(0.02_rad, 0.03_rad, 0.32_rad),
        quaternion_t::from_euler(0.03_rad, 0.01_rad, 0.28_rad),
        quaternion_t::from_euler(0.01_rad, 0.02_rad, 0.34_rad)
    };

    // Apply 5-point moving average filter directly to quaternions!
    auto const smoothed_orientations = sliding_window_average(
        std::span<quaternion_t const>{imu_orientations}, 5);
    assert(smoothed_orientations.size() == 16);  // 20 - 5 + 1 = 16

    // Verify first smoothed quaternion is average of first 5
    quaternion_t const expected_first = (
        imu_orientations[0] + imu_orientations[1] + imu_orientations[2] +
        imu_orientations[3] + imu_orientations[4]
    ) / 5.0;

    // Check components are approximately equal (simple component-wise averaging)
    assert(std::abs(smoothed_orientations[0].x() - expected_first.x()) < 0.001);
    assert(std::abs(smoothed_orientations[0].y() - expected_first.y()) < 0.001);
    assert(std::abs(smoothed_orientations[0].z() - expected_first.z()) < 0.001);
    assert(std::abs(smoothed_orientations[0].w() - expected_first.w()) < 0.001);

    std::cout << "All tests passed!\n";
    return 0;
}
```

**Learning objectives:**
- Non-owning views with `std::span`
- Generic functions that work with multiple container types
- Const-correctness with `std::span<T const>` vs `std::span<T>`
- Performance benefits of avoiding copies
- Function templates for generic algorithms
- Sliding window algorithms for signal processing

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

3. **Performance characteristics**:
   - Stack vs heap allocation
   - Zero-copy operations with views
   - Cache locality benefits

4. **Modern C++ patterns**:
   - Structured bindings with arrays
   - Using spans for function parameters
   - Compile-time size guarantees

## Discussion Questions

1. When should you use `std::array` instead of `std::vector`?
2. What are the lifetime considerations when using `std::span`?
3. How does `std::span` improve on the old pattern of passing `T* ptr, size_t length`?
4. Why does `std::span` support both dynamic and fixed extents?
5. What are the trade-offs between `std::array<T, N>` and C-style `T[N]` arrays?
6. How can `std::span` help with interfacing to C libraries?

---

## Solutions

✅ **[View Solutions with Detailed Explanations](week3_solutions.md)**
