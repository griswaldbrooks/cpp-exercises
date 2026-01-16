# Week 3: `std::array` and `std::span` - Solutions

## Solution 1: 3D Transformation Library

### Overview

This exercise demonstrates several modern C++ concepts:
1. **Strong Types**: Type-safe wrappers (meter_t, radian_t, degree_t) preventing unit errors
2. **Free Functions**: Non-member functions for operations (`distance`, `near`)
3. **Conversion Functions**: Standalone functions for type-safe conversions
4. **Encapsulation**: Private `std::array` storage with public accessors
5. **Const-correctness**: Separate const/non-const overloads
6. **User-defined literals**: Custom suffixes for units
7. **Operator overloading**: Mathematical operations
8. **Named constructors**: Alternative construction patterns

### Part A: `position_t` Class

**Key Concept: Const-Correctness with Accessor Overloads**

The position_t class demonstrates a fundamental C++ pattern: providing both const and non-const accessors to the same data.

```cpp
class position_t {
private:
    std::array<meter_t, 3> coords_;  // Strong type storage!

public:
    // Default constructor
    position_t() : coords_{meter_t{0.0}, meter_t{0.0}, meter_t{0.0}} {}

    // Parameterized constructor (accepts ONLY meter_t types - strong type safety!)
    position_t(meter_t const x, meter_t const y, meter_t const z)
        : coords_{x, y, z} {}

    // CONST accessors: return const reference to meter_t (read-only)
    [[nodiscard]] meter_t const& x() const { return coords_[0]; }
    [[nodiscard]] meter_t const& y() const { return coords_[1]; }
    [[nodiscard]] meter_t const& z() const { return coords_[2]; }

    // NON-CONST accessors: return mutable reference to meter_t (read-write)
    [[nodiscard]] meter_t& x() { return coords_[0]; }
    [[nodiscard]] meter_t& y() { return coords_[1]; }
    [[nodiscard]] meter_t& z() { return coords_[2]; }

    // Mathematical operations (use meter_t's operators)
    [[nodiscard]] position_t operator+(position_t const& other) const {
        return position_t{
            x() + other.x(),
            y() + other.y(),
            z() + other.z()
        };
    }

    [[nodiscard]] position_t operator-(position_t const& other) const {
        return position_t{
            x() - other.x(),
            y() - other.y(),
            z() - other.z()
        };
    }

    [[nodiscard]] position_t operator*(double const scalar) const {
        return position_t{
            x() * scalar,
            y() * scalar,
            z() * scalar
        };
    }

    [[nodiscard]] bool operator==(position_t const& other) const {
        return coords_ == other.coords_;
    }
};

// Free functions for position_t operations
[[nodiscard]] meter_t distance(position_t const& p1, position_t const& p2) {
    double const dx = (p1.x() - p2.x()).value;
    double const dy = (p1.y() - p2.y()).value;
    double const dz = (p1.z() - p2.z()).value;
    return meter_t{std::hypot(dx, dy, dz)};
}

[[nodiscard]] bool near(position_t const& p1, position_t const& p2,
                        meter_t const tolerance = meter_t{0.001}) {
    return std::abs((p1.x() - p2.x()).value) <= tolerance.value and
           std::abs((p1.y() - p2.y()).value) <= tolerance.value and
           std::abs((p1.z() - p2.z()).value) <= tolerance.value;
}
```

**Why Use `std::hypot` Instead of `std::sqrt`?**

The 3-argument `std::hypot(x, y, z)` (C++17) is preferred over `std::sqrt(x*x + y*y + z*z)` because:
1. **Avoids overflow/underflow**: `std::hypot` uses a more robust algorithm that won't overflow with large values or underflow with small values
2. **Better accuracy**: Avoids numerical issues from intermediate squared values
3. **Clearer intent**: Explicitly computing Euclidean distance
4. **Standard library optimization**: May use hardware-specific optimizations

**Why Free Functions?**

Using free functions `distance` and `near` (instead of member functions) provides several advantages:
1. **Symmetry**: Both positions are treated equally as parameters
2. **Discoverability**: Free functions appear in namespace scope, easier to find
3. **Extensibility**: Can be used with generic algorithms and concepts
4. **Non-member non-friend idiom**: Reduces coupling with the class

**Usage Examples:**

```cpp
using namespace literals;

position_t const origin{0.0_m, 0.0_m, 0.0_m};
position_t const p1{3.0_m, 4.0_m, 0.0_m};
position_t const p2{3.001_m, 4.001_m, 0.001_m};

// Calculate distance using free function (returns strong type!)
meter_t const dist = distance(origin, p1);  // dist.value == 5.0 (3-4-5 triangle)

// Check approximate equality with strong type tolerance
bool const are_close = near(p1, p2, 0.01_m);     // Returns true
bool const are_equal = near(p1, p2, 0.0001_m);   // Returns false

// Compare with symmetry (works with strong types)
distance(p1, p2) == distance(p2, p1);      // Always true (symmetric)
```

**Why Store `std::array<meter_t, 3>` Instead of `std::array<double, 3>`?**

By storing strong types directly, we enforce type safety at every level:
1. **No accidental raw doubles**: Can't accidentally assign a raw number
2. **Type-safe comparisons**: `operator==` works correctly with meter_t types
3. **Clearer intent**: The storage itself documents what it contains

**Note:** For `std::array<meter_t, 3>` to support `operator==`, the `meter` type must define `operator==`. We use defaulted comparison:

```cpp
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
```

**Why Const-Correctness Matters:**

```cpp
void process_position(position_t const& p) {
    meter_t const& val = p.x();  // OK: calls const overload, returns meter_t const&
    // p.x() = 5.0_m;  // ERROR: cannot modify through const reference
}

void modify_position(position_t& p) {
    p.x() = 10.0_m;  // OK: calls non-const overload, assigns meter
    meter_t const& val = p.x();  // OK: calls non-const overload (converts to const&)
}
```

### Part B: `quaternion_t` Class

**Key Concept: Quaternion Representation**

**Quaternion Components:**
- **w** is the **real (scalar) part**
- **x, y, z** are the **imaginary (vector) parts**
- Quaternions represent rotations in 3D space as: `q = w + xi + yj + zk`
- For a unit quaternion_t representing a rotation of angle θ around axis (ax, ay, az):
  - `w = cos(θ/2)`
  - `x = ax * sin(θ/2)`
  - `y = ay * sin(θ/2)`
  - `z = az * sin(θ/2)`

**Note on Normalization:** For calculating the magnitude of a 4D quaternion_t, we use nested `std::hypot` calls. Since C++17 only provides up to 3-argument `hypot`, we compute `hypot(hypot(x, y, z), w)` to get the same robustness benefits for all four components.

```cpp
class quaternion_t {
private:
    std::array<double, 4> components_;  // x, y, z, w

    // Private helper for normalization
    void normalize() {
        // Use nested hypot for 4D magnitude calculation
        double const norm = std::hypot(
            std::hypot(components_[0], components_[1], components_[2]),
            components_[3]
        );
        if (norm > 1e-10) {
            for (auto& component : components_) {
                component /= norm;
            }
        }
    }

public:
    // Default constructor: identity rotation
    quaternion_t() : components_{0.0, 0.0, 0.0, 1.0} {}

    // Parameterized constructor
    quaternion_t(double const x, double const y, double const z, double const w)
        : components_{x, y, z, w} {
        normalize();  // Always maintain unit quaternions
    }

    // Named constructor: alternative construction method (accepts ONLY radian_t type!)
    [[nodiscard]] static quaternion_t from_euler(radian_t const roll,
                                               radian_t const pitch,
                                               radian_t const yaw) {
        // Euler to quaternion_t conversion (extract values from strong types)
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

    // Accessors (same pattern as position_t)
    [[nodiscard]] double const& x() const { return components_[0]; }
    [[nodiscard]] double const& y() const { return components_[1]; }
    [[nodiscard]] double const& z() const { return components_[2]; }
    [[nodiscard]] double const& w() const { return components_[3]; }

    [[nodiscard]] double& x() { return components_[0]; }
    [[nodiscard]] double& y() { return components_[1]; }
    [[nodiscard]] double& z() { return components_[2]; }
    [[nodiscard]] double& w() { return components_[3]; }

    // Quaternion multiplication (composition of rotations)
    // Note: q1 * q2 applies rotation q2 first, then q1
    // References:
    // - https://math.stackexchange.com/questions/360286/what-does-multiplication-of-two-quaternions-give
    // - https://www.mathworks.com/help/aeroblks/quaternionmultiplication.html
    [[nodiscard]] quaternion_t operator*(quaternion_t const& other) const {
        return quaternion_t{
            w() * other.x() + x() * other.w() + y() * other.z() - z() * other.y(),
            w() * other.y() - x() * other.z() + y() * other.w() + z() * other.x(),
            w() * other.z() + x() * other.y() - y() * other.x() + z() * other.w(),
            w() * other.w() - x() * other.x() - y() * other.y() - z() * other.z()
        };
    }

    // Arithmetic operators for component-wise operations (needed for averaging)
    // Note: These perform simple component-wise operations, not geometric quaternion_t operations
    // This is a simplified approach suitable for basic sensor filtering in Exercise 2
    // For production code, consider using SLERP or other quaternion_t interpolation methods
    [[nodiscard]] quaternion_t operator+(quaternion_t const& other) const {
        return quaternion_t{
            x() + other.x(),
            y() + other.y(),
            z() + other.z(),
            w() + other.w()
        };
    }

    [[nodiscard]] quaternion_t operator-(quaternion_t const& other) const {
        return quaternion_t{
            x() - other.x(),
            y() - other.y(),
            z() - other.z(),
            w() - other.w()
        };
    }

    [[nodiscard]] quaternion_t operator/(double const scalar) const {
        return quaternion_t{
            x() / scalar,
            y() / scalar,
            z() / scalar,
            w() / scalar
        };
    }

    [[nodiscard]] quaternion_t conjugate() const {
        return quaternion_t{-x(), -y(), -z(), w()};
    }

    [[nodiscard]] bool operator==(quaternion_t const& other) const {
        return components_ == other.components_;
    }
};

// Free function for approximate quaternion_t equality
[[nodiscard]] bool near(quaternion_t const& q1, quaternion_t const& q2,
                        double const tolerance = 0.001) {
    return std::abs(q1.x() - q2.x()) <= tolerance and
           std::abs(q1.y() - q2.y()) <= tolerance and
           std::abs(q1.z() - q2.z()) <= tolerance and
           std::abs(q1.w() - q2.w()) <= tolerance;
}
```

**Why Use `= default`?**

- **Explicitly documents intent**: "I want the compiler-generated version"
- **Maintains triviality**: Compiler-generated special members can be trivial
- **Future-proof**: If you add members later, compiler updates behavior automatically
- **Performance**: May enable optimizations

### Part C: Conversion Functions

**Key Concept: Type-Safe Conversions Between Strong Types**

Before we define the user-defined literals, let's add conversion functions to convert between `degree` and `radian` types safely:

```cpp
namespace literals {
    // Strong type for angles in radians
    struct radian_t {
        double value;  // Stored in radians
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
        double value;  // Stored in degrees!
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

    // Conversion functions
    [[nodiscard]] constexpr radian_t to_radians(degree_t const d) {
        return radian_t{d.value * std::numbers::pi / 180.0};  // Convert degrees to radians
    }

    [[nodiscard]] constexpr degree_t to_degrees(radian_t const r) {
        return degree_t{r.value * 180.0 / std::numbers::pi};  // Convert radians to degrees
    }
}
```

**Why Conversion Functions?**

These standalone functions provide several advantages:
1. **Type-safe conversion**: Explicit conversion between degree_t and radian_t
2. **Bidirectional**: Can convert both ways (`to_radians` and `to_degrees`)
3. **Discoverable**: Easier to find via IDE autocomplete than member functions
4. **Clear intent**: Makes conversion explicit in code

**Usage Examples:**

```cpp
using namespace literals;

// Converting degrees to radians
auto const deg = 90.0_deg;
auto const rad = to_radians(deg);
EXPECT_NEAR(rad.value, std::numbers::pi / 2.0, 0.001);

// Converting radians to degrees
auto const rad2 = 3.14159_rad;
auto const deg2 = to_degrees(rad2);
EXPECT_NEAR(deg2.value, 180.0, 0.1);  // Now stored as degrees!

// Use in quaternion_t construction
auto const q = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg));
```

### Part D: User-Defined Literals for Strong Types

**Key Concept: Custom Suffixes for Type Safety and Readability**

We define user-defined literals for our basic strong types (meter_t, radian, degree):

```cpp
namespace literals {
    // Meter literals (floating-point and integer overloads required)
    [[nodiscard]] constexpr meter_t operator""_m(long double value) {
        return meter_t{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr meter_t operator""_m(unsigned long long value) {
        return meter_t{static_cast<double>(value)};
    }

    // Radian literals (returns radian_t type)
    [[nodiscard]] constexpr radian_t operator""_rad(long double value) {
        return radian_t{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr radian_t operator""_rad(unsigned long long value) {
        return radian_t{static_cast<double>(value)};
    }

    // Degree literals (returns degree_t type - stores as degrees!)
    [[nodiscard]] constexpr degree_t operator""_deg(long double value) {
        return degree_t{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr degree_t operator""_deg(unsigned long long value) {
        return degree_t{static_cast<double>(value)};
    }
}

// Usage:
using namespace literals;

position_t const p1{1.5_m, 2.3_m, 0.0_m};
auto const angle = 90.0_deg;  // Creates degree_t type (stored as 90.0)
auto const rot = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(angle));
```

**Benefits:**
- **Type safety**: Units are self-documenting
- **Compile-time conversion**: No runtime cost
- **Readability**: `90.0_deg` is clearer than `1.5707963...`

**Combined Usage with Conversion Functions:**

```cpp
using namespace literals;

// Create angles with literals
auto const angle1 = 90.0_deg;
auto const angle2 = 1.57_rad;

// Convert as needed
auto const rot = quaternion_t::from_euler(
    0.0_rad,
    0.0_rad,
    to_radians(angle1)  // Explicit, type-safe conversion
);

// Can also convert back
auto const angle1_in_rad = to_radians(angle1);
auto const angle2_in_deg = to_degrees(angle2);
```

### Part E: `transformation_t` Class

**Key Concept: 4x4 Homogeneous Transformation Matrix**

The transformation_t class stores a full 4x4 homogeneous transformation_t matrix internally, but only exposes `position_t` and `quaternion_t` types in its public interface. This demonstrates:
- **Information hiding**: Internal representation (matrix) vs. external interface (position_t/quaternion_t)
- **Type safety**: Users only work with position_t and quaternion_t, never raw matrix data
- **Mathematical correctness**: Uses proper SE(3) transformation_t matrix representation

**4x4 Transformation Matrix Format:**
```
[R00 R01 R02 | Tx]
[R10 R11 R12 | Ty]
[R20 R21 R22 | Tz]
[  0   0   0 |  1]
```
Where R is the 3x3 rotation matrix (from quaternion_t) and T is the translation vector (from position_t).

```cpp
class transformation_t {
private:
    std::array<double, 16> matrix_;  // 4x4 matrix in row-major order

    // Private constructor: directly initialize from a matrix
    // Used internally by operator* to avoid creating temporary position/quaternion
    explicit transformation_t(std::array<double, 16> const& mat) : matrix_{mat} {}

    // Helper: Convert quaternion_t to rotation matrix and build 4x4 transform
    static std::array<double, 16> build_matrix(position_t const& pos, quaternion_t const& rot) {
        // Extract quaternion_t components
        double const x = rot.x();
        double const y = rot.y();
        double const z = rot.z();
        double const w = rot.w();

        // Convert quaternion_t to rotation matrix (standard formula)
        double const xx = x * x, yy = y * y, zz = z * z;
        double const xy = x * y, xz = x * z, yz = y * z;
        double const wx = w * x, wy = w * y, wz = w * z;

        // Build 4x4 matrix in row-major order
        return {
            1 - 2*(yy + zz), 2*(xy - wz),     2*(xz + wy),     pos.x().value,  // Row 0
            2*(xy + wz),     1 - 2*(xx + zz), 2*(yz - wx),     pos.y().value,  // Row 1
            2*(xz - wy),     2*(yz + wx),     1 - 2*(xx + yy), pos.z().value,  // Row 2
            0.0,             0.0,             0.0,             1.0             // Row 3
        };
    }

public:
    // ONLY constructor: requires both position_t and rotation
    transformation_t(position_t const& pos, quaternion_t const& rot)
        : matrix_{build_matrix(pos, rot)} {}

    // Rule of 5 (all defaulted)
    transformation_t(transformation_t const& other) = default;
    transformation_t(transformation_t&& other) noexcept = default;
    transformation_t& operator=(transformation_t const& other) = default;
    transformation_t& operator=(transformation_t&& other) noexcept = default;
    ~transformation_t() = default;

    // Extract position_t (from last column of matrix)
    [[nodiscard]] position_t position() const {
        return position_t{
            meter_t{matrix_[3]},   // m03
            meter_t{matrix_[7]},   // m13
            meter_t{matrix_[11]}   // m23
        };
    }

    // Extract quaternion_t (convert rotation matrix to quaternion_t)
    [[nodiscard]] quaternion_t rotation() const {
        // Extract 3x3 rotation matrix elements
        double const m00 = matrix_[0], m01 = matrix_[1], m02 = matrix_[2];
        double const m10 = matrix_[4], m11 = matrix_[5], m12 = matrix_[6];
        double const m20 = matrix_[8], m21 = matrix_[9], m22 = matrix_[10];

        // Convert rotation matrix to quaternion using trace-based method
        // The trace determines which quaternion component is largest and most numerically stable
        double const trace = m00 + m11 + m22;

        // Case 1: w is the largest component (trace > 0)
        // This is the most common case for small rotations
        // We compute w first, then derive x, y, z from off-diagonal elements
        if (trace > 0) {
            double const s = 0.5 / std::sqrt(trace + 1.0);
            return quaternion_t{
                (m21 - m12) * s,
                (m02 - m20) * s,
                (m10 - m01) * s,
                0.25 / s
            };
        }

        // Case 2: x is the largest component (m00 is the largest diagonal element)
        // This occurs when rotation is primarily around the X axis
        // We compute x first, then derive y, z, w
        if (m00 > m11 && m00 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m00 - m11 - m22);
            return quaternion_t{
                0.25 * s,
                (m01 + m10) / s,
                (m02 + m20) / s,
                (m21 - m12) / s
            };
        }

        // Case 3: y is the largest component (m11 is the largest diagonal element)
        // This occurs when rotation is primarily around the Y axis
        // We compute y first, then derive x, z, w
        if (m11 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m11 - m00 - m22);
            return quaternion_t{
                (m01 + m10) / s,
                0.25 * s,
                (m12 + m21) / s,
                (m02 - m20) / s
            };
        }

        // Case 4: z is the largest component (m22 is the largest diagonal element)
        // This occurs when rotation is primarily around the Z axis
        // We compute z first, then derive x, y, w
        double const s = 2.0 * std::sqrt(1.0 + m22 - m00 - m11);
        return quaternion_t{
            (m02 + m20) / s,
            (m12 + m21) / s,
            0.25 * s,
            (m10 - m01) / s
        };
    }

    // Transformation composition (4x4 matrix multiplication)
    [[nodiscard]] transformation_t operator*(transformation_t const& other) const {
        std::array<double, 16> result{};

        // Multiply: result = this->matrix_ * other.matrix_
        // Row-major order: result[row*4 + col]
        // Composition applies 'other' transformation first, then 'this' transformation
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += matrix_[row * 4 + k] * other.matrix_[k * 4 + col];
                }
                result[row * 4 + col] = sum;
            }
        }

        // Use private constructor to avoid temporary position/quaternion creation
        return transformation_t{result};
    }

    // Apply additional rotation to this transformation
    // Equivalent to: *this * transformation_t{position_t{}, rotation}
    [[nodiscard]] transformation_t operator*(quaternion_t const& rotation) const {
        // Create a rotation-only transformation_t (identity position_t)
        transformation_t const rot_tf{position_t{meter_t{0}, meter_t{0}, meter_t{0}}, rotation};
        // Compose: apply rotation after this transformation
        return *this * rot_tf;
    }

    // Transform a position_t (apply rotation and translation)
    [[nodiscard]] position_t operator*(position_t const& point) const {
        // Treat point as [x, y, z, 1] in homogeneous coordinates
        // This allows us to apply both rotation and translation in a single matrix multiplication
        double const x = point.x().value;
        double const y = point.y().value;
        double const z = point.z().value;

        // Matrix-vector multiplication: result = matrix_ * [x, y, z, 1]
        // The '1' in the homogeneous coordinate enables translation
        // Result: rotated_point + translation
        return position_t{
            meter_t{matrix_[0]*x + matrix_[1]*y + matrix_[2]*z + matrix_[3]},
            meter_t{matrix_[4]*x + matrix_[5]*y + matrix_[6]*z + matrix_[7]},
            meter_t{matrix_[8]*x + matrix_[9]*y + matrix_[10]*z + matrix_[11]}
        };
    }

    // Inverse transformation_t (invert 4x4 matrix)
    [[nodiscard]] transformation_t inverse() const {
        // For SE(3) transformations, the inverse is:
        // [R^T | -R^T*t]
        // [0   | 1     ]
        // where R is the rotation matrix and t is the translation vector
        // This works because R is orthogonal (R^T * R = I for rotation matrices)

        // Extract rotation matrix (upper-left 3x3)
        double const r00 = matrix_[0], r01 = matrix_[1], r02 = matrix_[2];
        double const r10 = matrix_[4], r11 = matrix_[5], r12 = matrix_[6];
        double const r20 = matrix_[8], r21 = matrix_[9], r22 = matrix_[10];

        // Extract translation vector (last column of upper 3x4 part)
        double const tx = matrix_[3];
        double const ty = matrix_[7];
        double const tz = matrix_[11];

        // Compute -R^T * t
        // This gives the new translation after rotating the inverse direction
        double const new_tx = -(r00*tx + r10*ty + r20*tz);
        double const new_ty = -(r01*tx + r11*ty + r21*tz);
        double const new_tz = -(r02*tx + r12*ty + r22*tz);

        // Build inverse matrix with transposed rotation and new translation
        std::array<double, 16> const inv_matrix = {
            r00, r10, r20, new_tx,   // Row 0: R^T[0] | -R^T*t[0]
            r01, r11, r21, new_ty,   // Row 1: R^T[1] | -R^T*t[1]
            r02, r12, r22, new_tz,   // Row 2: R^T[2] | -R^T*t[2]
            0.0, 0.0, 0.0, 1.0       // Row 3
        };

        // Use private constructor to avoid temporary position/quaternion creation
        return transformation_t{inv_matrix};
    }

    [[nodiscard]] bool operator==(transformation_t const& other) const {
        return matrix_ == other.matrix_;
    }
};

// Free function for approximate transformation_t equality
[[nodiscard]] bool near(transformation_t const& tf1, transformation_t const& tf2,
                        meter_t const pos_tolerance = meter_t{0.001},
                        double const rot_tolerance = 0.001) {
    return near(tf1.position(), tf2.position(), pos_tolerance) and
           near(tf1.rotation(), tf2.rotation(), rot_tolerance);
}
```

**Why Store as 4x4 Matrix?**

1. **Mathematically correct**: Proper homogeneous transformation_t representation used in robotics/graphics
2. **Efficient composition**: Matrix multiplication is the natural way to compose transformations
3. **Industry standard**: SE(3) transformations are universally represented this way
4. **Encapsulation**: Users work with position_t/quaternion_t (intuitive) while internally using matrices (efficient)

**Why Use a Private Constructor?**

The private constructor `transformation_t(std::array<double, 16> const& mat)` is an optimization technique:

```cpp
// Without private constructor - inefficient:
transformation_t operator*(transformation_t const& other) const {
    std::array<double, 16> result = /* matrix multiplication */;
    transformation_t temp{position_t{}, quaternion_t{}};  // Creates temp position/quaternion
    temp.matrix_ = result;  // Then overwrites the matrix
    return temp;
}

// With private constructor - efficient:
transformation_t operator*(transformation_t const& other) const {
    std::array<double, 16> result = /* matrix multiplication */;
    return transformation_t{result};  // Directly construct from matrix
}
```

Benefits:
- **No temporary objects**: Avoids creating unnecessary `position_t` and `quaternion_t` objects
- **Better performance**: Fewer allocations and constructions
- **Cleaner code**: Direct construction from computed matrix
- **Encapsulation maintained**: Still private, so users can't misuse it

**Why No Default Constructor?**

```cpp
using namespace literals;

// transformation_t tf;  // ERROR: no default constructor
// This is intentional! A transformation_t without both position
// and rotation is meaningless in our domain.

transformation_t const tf{position_t{1.0_m, 2.0_m, 3.0_m}, quaternion_t{}};  // OK: explicit
```

**Why Use Independent `if` Statements in `rotation()` Instead of `else if`?**

The `rotation()` method converts a rotation matrix back to a quaternion using a trace-based algorithm. Notice it uses independent `if` statements rather than `else if`:

```cpp
if (trace > 0) { /* Case 1 */ return ...; }
if (m00 > m11 && m00 > m22) { /* Case 2 */ return ...; }
if (m11 > m22) { /* Case 3 */ return ...; }
/* Case 4 */ return ...;
```

**Why not use `else if`?**

Using independent `if` statements makes the code more explicit and self-documenting:

1. **Clarity**: Each case clearly shows its condition, making it obvious what triggers each branch
2. **Early returns**: Each case returns immediately, so `else` is unnecessary
3. **Easier debugging**: You can add breakpoints on each condition independently
4. **Pattern visibility**: The cascading structure shows we're checking for the "largest" component

**How the algorithm works:**

The algorithm chooses which quaternion component to compute first based on numerical stability:
- **Case 1** (`trace > 0`): w is largest → compute w, derive x, y, z (most common for small rotations)
- **Case 2** (`m00 largest`): x is largest → compute x, derive others (rotation around X axis)
- **Case 3** (`m11 largest`): y is largest → compute y, derive others (rotation around Y axis)
- **Case 4** (otherwise): z is largest → compute z, derive others (rotation around Z axis)

Each case avoids division by small numbers, ensuring numerical stability.

### Complete Usage Example

```cpp
using namespace literals;

// Robot base at origin
transformation_t const base{
    position_t{0.0_m, 0.0_m, 0.0_m},
    quaternion_t{}
};

// First joint: moves up 1 meter
transformation_t const joint1{
    position_t{0.0_m, 0.0_m, 1.0_m},
    quaternion_t{}
};

// Second joint: extends forward 0.5m, rotates 45 degrees
transformation_t const joint2{
    position_t{0.5_m, 0.0_m, 0.0_m},
    quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg))
};

// Compose transformations: base -> joint1 -> joint2
auto const end_effector = base * joint1 * joint2;

// Query end effector position
auto const pos = end_effector.position();
std::cout << std::format("End effector at: ({}, {}, {})\n",
                         pos.x(), pos.y(), pos.z());

// Example using conversion functions
auto const angle_deg = 180.0_deg;
auto const angle_rad = to_radians(angle_deg);
auto const rotation = quaternion_t::from_euler(angle_rad, 0.0_rad, 0.0_rad);
```

---

## Solution 2: Template Specialization for Rotation Averaging

### Overview

This exercise demonstrates **template specialization** by implementing type-appropriate averaging algorithms:

1. **Primary template**: Linear average for numeric types (`int`, `float`, `double`)
2. **Quaternion specialization**: Eigenvector method (Markley et al., 2007)
3. **Transformation specialization**: Matrix projection via polar decomposition

### Part A: Primary Template (Linear Average)

The primary template implements simple arithmetic mean for numeric types:

```cpp
// Primary template: linear average for numeric types
template<typename T>
T average(std::span<T const> values) {
    if (values.empty()) {
        return T{};
    }

    T sum = values[0];
    for (size_t i = 1; i < values.size(); ++i) {
        sum = sum + values[i];
    }
    return sum / static_cast<double>(values.size());
}
```

**Why start from `values[0]` instead of `T{0}`?**
- Not all types have a meaningful zero value
- For numeric types this is equivalent
- Avoids assumptions about default construction

### Part B: Quaternion Specialization (Eigenvector Method)

For quaternions, linear averaging produces geometrically incorrect results. The correct approach finds the quaternion that minimizes the sum of squared geodesic distances to all inputs.

**Algorithm (Markley et al., 2007):**
1. Build 4×4 matrix M = Σ qᵢqᵢᵀ (sum of outer products)
2. Find eigenvector corresponding to largest eigenvalue using power iteration
3. Normalize to unit quaternion

```cpp
// Full specialization for quaternion_t
template<>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> quaternions) {
    if (quaternions.empty()) {
        return quaternion_t{};  // Identity rotation
    }

    if (quaternions.size() == 1) {
        return quaternions[0];
    }

    // Build 4x4 accumulator matrix M = Σ qᵢ ⊗ qᵢᵀ
    // M is symmetric, so we only compute upper triangle
    std::array<std::array<double, 4>, 4> M{};

    for (auto const& q : quaternions) {
        std::array<double, 4> const qv = {q.x(), q.y(), q.z(), q.w()};

        // Add outer product qᵢ ⊗ qᵢᵀ to M
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                M[i][j] += qv[i] * qv[j];
            }
        }
    }

    // Power iteration to find eigenvector with largest eigenvalue
    // Initialize with [0, 0, 0, 1] (identity quaternion direction)
    std::array<double, 4> v = {0.0, 0.0, 0.0, 1.0};

    constexpr int max_iterations = 20;
    for (int iter = 0; iter < max_iterations; ++iter) {
        // v_new = M × v
        std::array<double, 4> v_new{};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                v_new[i] += M[i][j] * v[j];
            }
        }

        // Normalize v_new
        double const norm = std::sqrt(
            v_new[0] * v_new[0] + v_new[1] * v_new[1] +
            v_new[2] * v_new[2] + v_new[3] * v_new[3]
        );

        if (norm < 1e-10) {
            break;  // Degenerate case
        }

        for (size_t i = 0; i < 4; ++i) {
            v[i] = v_new[i] / norm;
        }
    }

    // Return as normalized quaternion
    return quaternion_t{v[0], v[1], v[2], v[3]};
}
```

**Why power iteration works:**
- The matrix M is symmetric positive semi-definite
- Its largest eigenvector points in the direction of maximum "agreement" among quaternions
- Power iteration converges to this dominant eigenvector
- 20 iterations is sufficient for typical precision requirements

### Part C: Transformation Specialization (Matrix Projection Method)

For transformations, we average position linearly and rotation via polar decomposition—**without converting to quaternions**.

**Algorithm:**
1. Average translations linearly: t_avg = (1/n) Σ tᵢ
2. Average rotation matrices element-wise: R_mean = (1/n) Σ Rᵢ
3. Project R_mean back to SO(3) via polar decomposition: R_avg = R_mean × (R_meanᵀ × R_mean)^(-½)

```cpp
// Helper: 3x3 matrix type
using Mat3 = std::array<std::array<double, 3>, 3>;

// Helper: 3x3 matrix multiplication
Mat3 mat3_multiply(Mat3 const& A, Mat3 const& B) {
    Mat3 result{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Helper: 3x3 matrix transpose
Mat3 mat3_transpose(Mat3 const& A) {
    Mat3 result{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

// Helper: 3x3 matrix inverse (using cofactor method)
Mat3 mat3_inverse(Mat3 const& A) {
    // Compute determinant
    double const det =
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    if (std::abs(det) < 1e-10) {
        // Return identity for singular matrix
        return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    }

    double const inv_det = 1.0 / det;

    Mat3 result{};
    result[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * inv_det;
    result[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * inv_det;
    result[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_det;
    result[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * inv_det;
    result[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_det;
    result[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * inv_det;
    result[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * inv_det;
    result[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * inv_det;
    result[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * inv_det;

    return result;
}

// Helper: 3x3 identity matrix
Mat3 mat3_identity() {
    return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
}

// Helper: Matrix addition and scalar multiplication
Mat3 mat3_add(Mat3 const& A, Mat3 const& B) {
    Mat3 result{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

Mat3 mat3_scale(Mat3 const& A, double s) {
    Mat3 result{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result[i][j] = A[i][j] * s;
        }
    }
    return result;
}

// Denman-Beavers iteration for matrix square root inverse
// Computes A^(-1/2) for symmetric positive-definite A
Mat3 mat3_sqrt_inverse(Mat3 const& A) {
    Mat3 Y = A;
    Mat3 Z = mat3_identity();

    constexpr int max_iterations = 20;
    for (int iter = 0; iter < max_iterations; ++iter) {
        Mat3 const Y_inv = mat3_inverse(Y);
        Mat3 const Z_inv = mat3_inverse(Z);

        Mat3 const Y_new = mat3_scale(mat3_add(Y, Z_inv), 0.5);
        Mat3 const Z_new = mat3_scale(mat3_add(Z, Y_inv), 0.5);

        Y = Y_new;
        Z = Z_new;
    }

    return Z;  // Z converges to A^(-1/2)
}

// Extract 3x3 rotation matrix from transformation_t
// (Assumes transformation_t stores 4x4 matrix in row-major order)
Mat3 extract_rotation(transformation_t const& tf) {
    // Get rotation by extracting the quaternion and converting to matrix
    quaternion_t const q = tf.rotation();
    double const x = q.x(), y = q.y(), z = q.z(), w = q.w();

    double const xx = x * x, yy = y * y, zz = z * z;
    double const xy = x * y, xz = x * z, yz = y * z;
    double const wx = w * x, wy = w * y, wz = w * z;

    return {{
        {1 - 2*(yy + zz), 2*(xy - wz),     2*(xz + wy)    },
        {2*(xy + wz),     1 - 2*(xx + zz), 2*(yz - wx)    },
        {2*(xz - wy),     2*(yz + wx),     1 - 2*(xx + yy)}
    }};
}

// Convert rotation matrix to quaternion
quaternion_t mat3_to_quaternion(Mat3 const& R) {
    double const trace = R[0][0] + R[1][1] + R[2][2];

    if (trace > 0) {
        double const s = 0.5 / std::sqrt(trace + 1.0);
        return quaternion_t{
            (R[2][1] - R[1][2]) * s,
            (R[0][2] - R[2][0]) * s,
            (R[1][0] - R[0][1]) * s,
            0.25 / s
        };
    }

    if (R[0][0] > R[1][1] && R[0][0] > R[2][2]) {
        double const s = 2.0 * std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]);
        return quaternion_t{
            0.25 * s,
            (R[0][1] + R[1][0]) / s,
            (R[0][2] + R[2][0]) / s,
            (R[2][1] - R[1][2]) / s
        };
    }

    if (R[1][1] > R[2][2]) {
        double const s = 2.0 * std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]);
        return quaternion_t{
            (R[0][1] + R[1][0]) / s,
            0.25 * s,
            (R[1][2] + R[2][1]) / s,
            (R[0][2] - R[2][0]) / s
        };
    }

    double const s = 2.0 * std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]);
    return quaternion_t{
        (R[0][2] + R[2][0]) / s,
        (R[1][2] + R[2][1]) / s,
        0.25 * s,
        (R[1][0] - R[0][1]) / s
    };
}

// Full specialization for transformation_t
template<>
transformation_t average<transformation_t>(std::span<transformation_t const> transforms) {
    if (transforms.empty()) {
        return transformation_t{
            position_t{meter_t{0}, meter_t{0}, meter_t{0}},
            quaternion_t{}
        };
    }

    if (transforms.size() == 1) {
        return transforms[0];
    }

    double const n = static_cast<double>(transforms.size());

    // Step 1: Average positions linearly
    double sum_x = 0, sum_y = 0, sum_z = 0;
    for (auto const& tf : transforms) {
        position_t const pos = tf.position();
        sum_x += pos.x().value;
        sum_y += pos.y().value;
        sum_z += pos.z().value;
    }
    position_t const avg_pos{
        meter_t{sum_x / n},
        meter_t{sum_y / n},
        meter_t{sum_z / n}
    };

    // Step 2: Average rotation matrices element-wise
    Mat3 R_mean{};
    for (auto const& tf : transforms) {
        Mat3 const R = extract_rotation(tf);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                R_mean[i][j] += R[i][j] / n;
            }
        }
    }

    // Step 3: Project R_mean to SO(3) via polar decomposition
    // R_avg = R_mean × (R_meanᵀ × R_mean)^(-1/2)
    Mat3 const Rt = mat3_transpose(R_mean);
    Mat3 const RtR = mat3_multiply(Rt, R_mean);
    Mat3 const RtR_inv_sqrt = mat3_sqrt_inverse(RtR);
    Mat3 const R_avg = mat3_multiply(R_mean, RtR_inv_sqrt);

    // Convert rotation matrix to quaternion
    quaternion_t const avg_rot = mat3_to_quaternion(R_avg);

    return transformation_t{avg_pos, avg_rot};
}
```

**Why polar decomposition works:**
- Any matrix A can be decomposed as A = UP where U is orthogonal and P is positive semi-definite
- For our rotation average, R_mean = R_avg × P where R_avg ∈ SO(3)
- The formula R_avg = R_mean × (R_meanᵀ × R_mean)^(-½) extracts the orthogonal part
- This finds the rotation matrix closest to R_mean in the Frobenius norm

### Why Two Independent Methods?

The test suite validates that quaternion averaging and matrix averaging produce **identical results**:

```cpp
// Both methods should give the same answer!
quaternion_t avg_quat = average<quaternion_t>(quats);
transformation_t avg_tf = average<transformation_t>(transforms);
assert(near(avg_quat, avg_tf.rotation(), 0.001));
```

By implementing both independently:
- We verify mathematical correctness of each algorithm
- If one used the other internally, the test would be circular
- Students learn two different approaches to the same problem

### Usage Examples

```cpp
using namespace literals;

// Linear average for numeric types (uses primary template)
std::array<double, 3> values = {1.0, 2.0, 3.0};
double avg = average<double>(values);  // = 2.0

// Quaternion average (uses specialization)
std::array<quaternion_t, 3> quats = {
    quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.1_rad),
    quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.2_rad),
    quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.3_rad)
};
quaternion_t avg_quat = average<quaternion_t>(quats);

// Transformation average (uses specialization)
std::array<transformation_t, 3> transforms = {
    transformation_t{position_t{1.0_m, 0.0_m, 0.0_m},
                     quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.1_rad)},
    transformation_t{position_t{2.0_m, 0.0_m, 0.0_m},
                     quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.2_rad)},
    transformation_t{position_t{3.0_m, 0.0_m, 0.0_m},
                     quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.3_rad)}
};
transformation_t avg_tf = average<transformation_t>(transforms);
// avg_tf.position() = {2.0m, 0.0m, 0.0m}  (linear average)
// avg_tf.rotation() ≈ avg_quat            (spherical average)
```

### References

- Markley, F.L., Cheng, Y., Crassidis, J.L., Oshman, Y. (2007). *"Averaging Quaternions"*, Journal of Guidance, Control, and Dynamics, 30(4), 1193-1197.
- Moakher, M. (2002). *"Means and Averaging in the Group of Rotations"*, SIAM Journal on Matrix Analysis and Applications.
- Higham, N.J. (1986). *"Computing the Polar Decomposition—with Applications"*, SIAM Journal on Scientific and Statistical Computing.

---

## Key Takeaways

### When to Use Template Specialization

✅ **Use when:**
- Different types need fundamentally different algorithms (e.g., rotation averaging)
- The primary template makes no sense for certain types
- You need to optimize for specific types without changing the interface

❌ **Don't use when:**
- Simple overloading would suffice
- The algorithm is truly generic across all types
- You're just handling edge cases (use `if constexpr` instead)

### Rotation Averaging Methods

| Method | Used For | Algorithm | Complexity |
|--------|----------|-----------|------------|
| Linear | Scalars, vectors | (Σxᵢ) / n | O(n) |
| Eigenvector | Quaternions | Power iteration on M = Σqᵢqᵢᵀ | O(n + k) |
| Polar decomposition | Rotation matrices | Denman-Beavers iteration | O(n + k) |

Where k = number of iterations (typically 20).

### Why Linear Averaging Fails for Rotations

```cpp
// Two 180° rotations around different axes
quaternion_t q1 = from_euler(180_deg, 0, 0);   // [1, 0, 0, 0]
quaternion_t q2 = from_euler(0, 180_deg, 0);   // [0, 1, 0, 0]

// Linear average: [0.5, 0.5, 0, 0] - NOT a valid rotation!
// Eigenvector method: finds geometrically meaningful average
```

Linear averaging fails because:
- Quaternions lie on the unit 3-sphere (S³)
- Linear combination leaves the sphere
- The "average" must stay on the manifold

### Template Specialization Syntax

```cpp
// Primary template
template<typename T>
T average(std::span<T const> values);

// Full specialization - note empty angle brackets
template<>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> values);
```

Key points:
- `template<>` signals full specialization
- Function signature must match primary template exactly
- Compiler selects most specific match at compile time

---

## Discussion Questions & Answers

**Q1: When should you use `std::array` instead of `std::vector`?**

A: Use `std::array` when:
- Size is known at compile time (e.g., 3D coordinates, rotation matrices)
- Performance is critical (stack allocation is faster)
- You want compile-time size checking
- Working in embedded/constrained environments

**Q2: What are the lifetime considerations when using `std::span`?**

A: `std::span` is a non-owning view. You must ensure:
- The underlying data outlives the span
- The data isn't reallocated (e.g., `vector::push_back` invalidates)
- You don't return spans to local variables

**Q3: Why can't you use linear averaging for rotations?**

A: Rotations form a curved space (manifold), not a flat vector space:
- Quaternions lie on the unit 3-sphere (S³)
- Rotation matrices lie on SO(3)
- Linear combinations leave these surfaces
- The average must be "intrinsic" to the manifold (Fréchet mean)

**Q4: When would you use partial vs full specialization?**

A: Full specialization (what we used) is for a single specific type.
Partial specialization (not supported for functions) would be for a family of types, e.g., all pointers.

For functions, use overloading or `if constexpr` for partial-specialization-like behavior:
```cpp
template<typename T>
T average(std::span<T const> values) {
    if constexpr (std::is_same_v<T, quaternion_t>) {
        // Quaternion-specific code
    } else {
        // Generic code
    }
}
```

**Q5: How does power iteration find the dominant eigenvector?**

A: Power iteration repeatedly multiplies a vector by M and normalizes:
1. Start with any non-zero vector v
2. Compute v' = Mv / ‖Mv‖
3. Repeat until convergence

The vector converges to the eigenvector with largest eigenvalue because:
- Each multiplication amplifies components in the dominant eigenvector's direction
- Normalization prevents overflow
- Other components shrink relative to the dominant one

**Q6: Why do the quaternion and matrix methods produce the same result?**

A: Both methods are solving the same optimization problem:
- Find the rotation that minimizes the sum of squared geodesic distances to all inputs

The eigenvector method works in quaternion space (S³).
The polar decomposition works in rotation matrix space (SO(3)).
Since quaternions and rotation matrices are different representations of the same rotations, the solutions must agree.
