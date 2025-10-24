# Week 3: `std::array` and `std::span` - Solutions

## Solution 1: 3D Transformation Library

### Overview

This exercise demonstrates several modern C++ concepts:
1. **Strong Types**: Type-safe wrappers (meter, radian, degree) preventing unit errors
2. **Free Functions**: Non-member functions for operations (`distance_to`, `near`)
3. **Conversion Functions**: Standalone functions for type-safe conversions
4. **Encapsulation**: Private `std::array` storage with public accessors
5. **Const-correctness**: Separate const/non-const overloads
6. **Rule of 5**: Explicit special member functions
7. **User-defined literals**: Custom suffixes for units
8. **Operator overloading**: Mathematical operations
9. **Named constructors**: Alternative construction patterns

### Part A: `position` Class

**Key Concept: Const-Correctness with Accessor Overloads**

The position class demonstrates a fundamental C++ pattern: providing both const and non-const accessors to the same data.

```cpp
class position {
private:
    std::array<meter, 3> coords_;  // Strong type storage!

public:
    // Default constructor
    position() : coords_{meter{0.0}, meter{0.0}, meter{0.0}} {}

    // Parameterized constructor (accepts ONLY meter types - strong type safety!)
    position(meter const x, meter const y, meter const z)
        : coords_{x, y, z} {}

    // CONST accessors: return const reference to meter (read-only)
    [[nodiscard]] meter const& x() const { return coords_[0]; }
    [[nodiscard]] meter const& y() const { return coords_[1]; }
    [[nodiscard]] meter const& z() const { return coords_[2]; }

    // NON-CONST accessors: return mutable reference to meter (read-write)
    [[nodiscard]] meter& x() { return coords_[0]; }
    [[nodiscard]] meter& y() { return coords_[1]; }
    [[nodiscard]] meter& z() { return coords_[2]; }

    // Mathematical operations (work with meter's value member)
    [[nodiscard]] position operator+(position const& other) const {
        return position{
            meter{x().value + other.x().value},
            meter{y().value + other.y().value},
            meter{z().value + other.z().value}
        };
    }

    [[nodiscard]] position operator-(position const& other) const {
        return position{
            meter{x().value - other.x().value},
            meter{y().value - other.y().value},
            meter{z().value - other.z().value}
        };
    }

    [[nodiscard]] position operator*(double const scalar) const {
        return position{
            meter{x().value * scalar},
            meter{y().value * scalar},
            meter{z().value * scalar}
        };
    }

    [[nodiscard]] bool operator==(position const& other) const {
        return coords_ == other.coords_;
    }
};

// Free functions for position operations
[[nodiscard]] meter distance(position const& p1, position const& p2) {
    double const dx = p1.x().value - p2.x().value;
    double const dy = p1.y().value - p2.y().value;
    double const dz = p1.z().value - p2.z().value;
    return meter{std::hypot(dx, dy, dz)};
}

[[nodiscard]] bool near(position const& p1, position const& p2,
                        meter const tolerance = meter{0.001}) {
    return std::abs(p1.x().value - p2.x().value) <= tolerance.value and
           std::abs(p1.y().value - p2.y().value) <= tolerance.value and
           std::abs(p1.z().value - p2.z().value) <= tolerance.value;
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

position const origin{0.0_m, 0.0_m, 0.0_m};
position const p1{3.0_m, 4.0_m, 0.0_m};
position const p2{3.001_m, 4.001_m, 0.001_m};

// Calculate distance using free function (returns strong type!)
meter const dist = distance(origin, p1);  // dist.value == 5.0 (3-4-5 triangle)

// Check approximate equality with strong type tolerance
bool const are_close = near(p1, p2, 0.01_m);     // Returns true
bool const are_equal = near(p1, p2, 0.0001_m);   // Returns false

// Compare with symmetry (works with strong types)
distance(p1, p2) == distance(p2, p1);      // Always true (symmetric)
```

**Why Store `std::array<meter, 3>` Instead of `std::array<double, 3>`?**

By storing strong types directly, we enforce type safety at every level:
1. **No accidental raw doubles**: Can't accidentally assign a raw number
2. **Type-safe comparisons**: `operator==` works correctly with meter types
3. **Clearer intent**: The storage itself documents what it contains

**Note:** For `std::array<meter, 3>` to support `operator==`, the `meter` type must define `operator==`. We use defaulted comparison:

```cpp
struct meter {
    double value;
    constexpr explicit meter(double v) : value{v} {}

    [[nodiscard]] constexpr bool operator==(meter const& other) const = default;
};
```

**Why Const-Correctness Matters:**

```cpp
void process_position(position const& p) {
    meter const& val = p.x();  // OK: calls const overload, returns meter const&
    // p.x() = 5.0_m;  // ERROR: cannot modify through const reference
}

void modify_position(position& p) {
    p.x() = 10.0_m;  // OK: calls non-const overload, assigns meter
    meter const& val = p.x();  // OK: calls non-const overload (converts to const&)
}
```

### Part B: `quaternion` Class

**Key Concept: Rule of 5**

The quaternion class explicitly implements all five special member functions, demonstrating when and why to use `= default`.

**Quaternion Components:**
- **w** is the **real (scalar) part**
- **x, y, z** are the **imaginary (vector) parts**
- Quaternions represent rotations in 3D space as: `q = w + xi + yj + zk`
- For a unit quaternion representing a rotation of angle θ around axis (ax, ay, az):
  - `w = cos(θ/2)`
  - `x = ax * sin(θ/2)`
  - `y = ay * sin(θ/2)`
  - `z = az * sin(θ/2)`

**Note on Normalization:** For calculating the magnitude of a 4D quaternion, we use nested `std::hypot` calls. Since C++17 only provides up to 3-argument `hypot`, we compute `hypot(hypot(x, y, z), w)` to get the same robustness benefits for all four components.

```cpp
class quaternion {
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
    quaternion() : components_{0.0, 0.0, 0.0, 1.0} {}

    // Parameterized constructor
    quaternion(double const x, double const y, double const z, double const w)
        : components_{x, y, z, w} {
        normalize();  // Always maintain unit quaternions
    }

    // Named constructor: alternative construction method (accepts ONLY radian type!)
    [[nodiscard]] static quaternion from_euler(radian const roll,
                                               radian const pitch,
                                               radian const yaw) {
        // Euler to quaternion conversion (extract values from strong types)
        double const cy = std::cos(yaw.value * 0.5);
        double const sy = std::sin(yaw.value * 0.5);
        double const cp = std::cos(pitch.value * 0.5);
        double const sp = std::sin(pitch.value * 0.5);
        double const cr = std::cos(roll.value * 0.5);
        double const sr = std::sin(roll.value * 0.5);

        return quaternion{
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
            cr * cp * cy + sr * sp * sy
        };
    }

    // Rule of 5: Explicitly defaulted special member functions
    quaternion(quaternion const& other) = default;              // Copy constructor
    quaternion(quaternion&& other) noexcept = default;          // Move constructor
    quaternion& operator=(quaternion const& other) = default;    // Copy assignment
    quaternion& operator=(quaternion&& other) noexcept = default; // Move assignment
    ~quaternion() = default;                                     // Destructor

    // Accessors (same pattern as position)
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
    [[nodiscard]] quaternion operator*(quaternion const& other) const {
        return quaternion{
            w() * other.x() + x() * other.w() + y() * other.z() - z() * other.y(),
            w() * other.y() - x() * other.z() + y() * other.w() + z() * other.x(),
            w() * other.z() + x() * other.y() - y() * other.x() + z() * other.w(),
            w() * other.w() - x() * other.x() - y() * other.y() - z() * other.z()
        };
    }

    // Arithmetic operators for component-wise operations (needed for averaging)
    // Note: These perform simple component-wise operations, not geometric quaternion operations
    // This is a simplified approach suitable for basic sensor filtering in Exercise 2
    // For production code, consider using SLERP or other quaternion interpolation methods
    [[nodiscard]] quaternion operator+(quaternion const& other) const {
        return quaternion{
            x() + other.x(),
            y() + other.y(),
            z() + other.z(),
            w() + other.w()
        };
    }

    [[nodiscard]] quaternion operator-(quaternion const& other) const {
        return quaternion{
            x() - other.x(),
            y() - other.y(),
            z() - other.z(),
            w() - other.w()
        };
    }

    [[nodiscard]] quaternion operator/(double const scalar) const {
        return quaternion{
            x() / scalar,
            y() / scalar,
            z() / scalar,
            w() / scalar
        };
    }

    [[nodiscard]] quaternion conjugate() const {
        return quaternion{-x(), -y(), -z(), w()};
    }

    [[nodiscard]] bool operator==(quaternion const& other) const {
        return components_ == other.components_;
    }
};

// Free function for approximate quaternion equality
[[nodiscard]] bool near(quaternion const& q1, quaternion const& q2,
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
    struct radian {
        double value;  // Stored in radians
        constexpr explicit radian(double v) : value{v} {}
    };

    // Strong type for angles in degrees
    struct degree {
        double value;  // Stored in degrees!
        constexpr explicit degree(double deg) : value{deg} {}
    };

    // Conversion functions (work with both types)
    [[nodiscard]] constexpr radian to_radians(degree const d) {
        return radian{d.value * std::numbers::pi / 180.0};  // Convert degrees to radians
    }

    [[nodiscard]] constexpr radian to_radians(radian const r) {
        return r;  // No conversion needed (idempotent)
    }

    [[nodiscard]] constexpr degree to_degrees(radian const r) {
        return degree{r.value * 180.0 / std::numbers::pi};  // Convert radians to degrees
    }

    [[nodiscard]] constexpr degree to_degrees(degree const d) {
        return d;  // No conversion needed (idempotent)
    }
}
```

**Why Conversion Functions?**

These standalone functions provide several advantages:
1. **Consistent interface**: Same function name works for both types
2. **Idempotent**: `to_radians(radian)` returns the same radian (useful for generic code)
3. **Bidirectional**: Can convert both ways (`to_radians` and `to_degrees`)
4. **Discoverable**: Easier to find via IDE autocomplete than member functions

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

// Idempotent operations (useful for generic code)
auto const rad3 = 1.57_rad;
auto const rad4 = to_radians(rad3);  // Returns same value
EXPECT_DOUBLE_EQ(rad3.value, rad4.value);

// Use in quaternion construction
auto const q = quaternion::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg));
```

### Part D: User-Defined Literals for Strong Types

**Key Concept: Custom Suffixes for Type Safety and Readability**

We define user-defined literals for our basic strong types (meter, radian, degree):

```cpp
namespace literals {
    // Meter literals (floating-point and integer overloads required)
    [[nodiscard]] constexpr meter operator""_m(long double value) {
        return meter{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr meter operator""_m(unsigned long long value) {
        return meter{static_cast<double>(value)};
    }

    // Radian literals (returns radian type)
    [[nodiscard]] constexpr radian operator""_rad(long double value) {
        return radian{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr radian operator""_rad(unsigned long long value) {
        return radian{static_cast<double>(value)};
    }

    // Degree literals (returns degree type - stores as degrees!)
    [[nodiscard]] constexpr degree operator""_deg(long double value) {
        return degree{static_cast<double>(value)};
    }

    [[nodiscard]] constexpr degree operator""_deg(unsigned long long value) {
        return degree{static_cast<double>(value)};
    }
}

// Usage:
using namespace literals;

position const p1{1.5_m, 2.3_m, 0.0_m};
auto const angle = 90.0_deg;  // Creates degree type (stored as 90.0)
auto const rot = quaternion::from_euler(0.0_rad, 0.0_rad, to_radians(angle));
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
auto const rot = quaternion::from_euler(
    0.0_rad,
    0.0_rad,
    to_radians(angle1)  // Explicit, type-safe conversion
);

// Can also convert back
auto const angle1_in_rad = to_radians(angle1);
auto const angle2_in_deg = to_degrees(angle2);
```

### Part E: `transformation` Class

**Key Concept: 4x4 Homogeneous Transformation Matrix**

The transformation class stores a full 4x4 homogeneous transformation matrix internally, but only exposes `position` and `quaternion` types in its public interface. This demonstrates:
- **Information hiding**: Internal representation (matrix) vs. external interface (position/quaternion)
- **Type safety**: Users only work with position and quaternion, never raw matrix data
- **Mathematical correctness**: Uses proper SE(3) transformation matrix representation

**4x4 Transformation Matrix Format:**
```
[R00 R01 R02 | Tx]
[R10 R11 R12 | Ty]
[R20 R21 R22 | Tz]
[  0   0   0 |  1]
```
Where R is the 3x3 rotation matrix (from quaternion) and T is the translation vector (from position).

```cpp
class transformation {
private:
    std::array<double, 16> matrix_;  // 4x4 matrix in row-major order

    // Helper: Convert quaternion to rotation matrix and build 4x4 transform
    static std::array<double, 16> build_matrix(position const& pos, quaternion const& rot) {
        // Extract quaternion components
        double const x = rot.x();
        double const y = rot.y();
        double const z = rot.z();
        double const w = rot.w();

        // Convert quaternion to rotation matrix (standard formula)
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
    // ONLY constructor: requires both position and rotation
    transformation(position const& pos, quaternion const& rot)
        : matrix_{build_matrix(pos, rot)} {}

    // Rule of 5 (all defaulted)
    transformation(transformation const& other) = default;
    transformation(transformation&& other) noexcept = default;
    transformation& operator=(transformation const& other) = default;
    transformation& operator=(transformation&& other) noexcept = default;
    ~transformation() = default;

    // Extract position (from last column of matrix)
    [[nodiscard]] position position() const {
        return position{
            meter{matrix_[3]},   // m03
            meter{matrix_[7]},   // m13
            meter{matrix_[11]}   // m23
        };
    }

    // Extract quaternion (convert rotation matrix to quaternion)
    [[nodiscard]] quaternion rotation() const {
        // Extract 3x3 rotation matrix
        double const m00 = matrix_[0], m01 = matrix_[1], m02 = matrix_[2];
        double const m10 = matrix_[4], m11 = matrix_[5], m12 = matrix_[6];
        double const m20 = matrix_[8], m21 = matrix_[9], m22 = matrix_[10];

        // Convert rotation matrix to quaternion (standard algorithm)
        double const trace = m00 + m11 + m22;

        if (trace > 0) {
            double const s = 0.5 / std::sqrt(trace + 1.0);
            return quaternion{
                (m21 - m12) * s,
                (m02 - m20) * s,
                (m10 - m01) * s,
                0.25 / s
            };
        } else if (m00 > m11 && m00 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m00 - m11 - m22);
            return quaternion{
                0.25 * s,
                (m01 + m10) / s,
                (m02 + m20) / s,
                (m21 - m12) / s
            };
        } else if (m11 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m11 - m00 - m22);
            return quaternion{
                (m01 + m10) / s,
                0.25 * s,
                (m12 + m21) / s,
                (m02 - m20) / s
            };
        } else {
            double const s = 2.0 * std::sqrt(1.0 + m22 - m00 - m11);
            return quaternion{
                (m02 + m20) / s,
                (m12 + m21) / s,
                0.25 * s,
                (m10 - m01) / s
            };
        }
    }

    // Transformation composition (4x4 matrix multiplication)
    [[nodiscard]] transformation operator*(transformation const& other) const {
        std::array<double, 16> result{};

        // Multiply: result = this->matrix_ * other.matrix_
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += matrix_[row * 4 + k] * other.matrix_[k * 4 + col];
                }
                result[row * 4 + col] = sum;
            }
        }

        // Create transformation from resulting matrix
        transformation temp{position{}, quaternion{}};
        temp.matrix_ = result;
        return temp;
    }

    // Apply additional rotation to this transformation
    // Equivalent to: *this * transformation{position{}, rotation}
    [[nodiscard]] transformation operator*(quaternion const& rotation) const {
        // Create a rotation-only transformation (identity position)
        transformation const rot_tf{position{meter{0}, meter{0}, meter{0}}, rotation};
        // Compose: apply rotation after this transformation
        return *this * rot_tf;
    }

    // Transform a position (apply 4x4 matrix to homogeneous coordinate)
    [[nodiscard]] position operator*(position const& point) const {
        // Treat point as [x, y, z, 1] in homogeneous coordinates
        double const x = point.x().value;
        double const y = point.y().value;
        double const z = point.z().value;

        // Matrix-vector multiplication
        return position{
            meter{matrix_[0]*x + matrix_[1]*y + matrix_[2]*z + matrix_[3]},
            meter{matrix_[4]*x + matrix_[5]*y + matrix_[6]*z + matrix_[7]},
            meter{matrix_[8]*x + matrix_[9]*y + matrix_[10]*z + matrix_[11]}
        };
    }

    // Inverse transformation (for SE(3): [R^T | -R^T*t])
    [[nodiscard]] transformation inverse() const {
        // Extract rotation and translation
        double const r00 = matrix_[0], r01 = matrix_[1], r02 = matrix_[2];
        double const r10 = matrix_[4], r11 = matrix_[5], r12 = matrix_[6];
        double const r20 = matrix_[8], r21 = matrix_[9], r22 = matrix_[10];
        double const tx = matrix_[3], ty = matrix_[7], tz = matrix_[11];

        // Compute -R^T * t
        double const new_tx = -(r00*tx + r10*ty + r20*tz);
        double const new_ty = -(r01*tx + r11*ty + r21*tz);
        double const new_tz = -(r02*tx + r12*ty + r22*tz);

        std::array<double, 16> inv_matrix = {
            r00, r10, r20, new_tx,   // Row 0
            r01, r11, r21, new_ty,   // Row 1
            r02, r12, r22, new_tz,   // Row 2
            0.0, 0.0, 0.0, 1.0       // Row 3
        };

        transformation result{position{}, quaternion{}};
        result.matrix_ = inv_matrix;
        return result;
    }

    [[nodiscard]] bool operator==(transformation const& other) const {
        return matrix_ == other.matrix_;
    }
};

// Free function for approximate transformation equality
[[nodiscard]] bool near(transformation const& tf1, transformation const& tf2,
                        meter const pos_tolerance = meter{0.001},
                        double const rot_tolerance = 0.001) {
    return near(tf1.position(), tf2.position(), pos_tolerance) and
           near(tf1.rotation(), tf2.rotation(), rot_tolerance);
}
```

**Why Store as 4x4 Matrix?**

1. **Mathematically correct**: Proper homogeneous transformation representation used in robotics/graphics
2. **Efficient composition**: Matrix multiplication is the natural way to compose transformations
3. **Industry standard**: SE(3) transformations are universally represented this way
4. **Encapsulation**: Users work with position/quaternion (intuitive) while internally using matrices (efficient)

**Why No Default Constructor?**

```cpp
using namespace literals;

// transformation tf;  // ERROR: no default constructor
// This is intentional! A transformation without both position
// and rotation is meaningless in our domain.

transformation const tf{position{1.0_m, 2.0_m, 3.0_m}, quaternion{}};  // OK: explicit
```

### Complete Usage Example

```cpp
using namespace literals;

// Robot base at origin
transformation const base{
    position{0.0_m, 0.0_m, 0.0_m},
    quaternion{}
};

// First joint: moves up 1 meter
transformation const joint1{
    position{0.0_m, 0.0_m, 1.0_m},
    quaternion{}
};

// Second joint: extends forward 0.5m, rotates 45 degrees
transformation const joint2{
    position{0.5_m, 0.0_m, 0.0_m},
    quaternion::from_euler(0.0_rad, 0.0_rad, to_radians(45.0_deg))
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
auto const rotation = quaternion::from_euler(angle_rad, 0.0_rad, 0.0_rad);
```

---

## Solution 2: Sensor Array Statistics

### Key Concepts

**std::span benefits:**
- Zero-copy view of contiguous data
- Works with `std::array`, `std::vector`, C-arrays
- Const-correctness: `std::span<T const>` vs `std::span<T>`
- No allocations (just pointer + size)

```cpp
double average(std::span<double const> const readings) {
    if (readings.empty()) {
        return 0.0;
    }
    double const sum = std::accumulate(readings.begin(), readings.end(), 0.0);
    return sum / static_cast<double>(readings.size());
}

double min_value(std::span<double const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<double>::max();
    }
    return *std::min_element(readings.begin(), readings.end());
}

double max_value(std::span<double const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<double>::lowest();
    }
    return *std::max_element(readings.begin(), readings.end());
}

size_t count_above_threshold(std::span<double const> const readings,
                             double const threshold) {
    return std::count_if(readings.begin(), readings.end(),
                        [threshold](double const value) { return value > threshold; });
}

void normalize(std::span<double> const readings) {
    if (readings.empty()) {
        return;
    }

    double const min = min_value(readings);
    double const max = max_value(readings);
    double const range = max - min;

    if (range == 0.0) {
        std::fill(readings.begin(), readings.end(), 0.5);
        return;
    }

    for (auto& value : readings) {
        value = (value - min) / range;
    }
}

// Sliding window average - templated for any numeric type
template<typename T>
std::vector<T> sliding_window_average(std::span<T const> const values, size_t const window_size) {
    if (values.size() < window_size || window_size == 0) {
        return {};
    }

    std::vector<T> result;
    result.reserve(values.size() - window_size + 1);

    // Compute sum of first window
    T window_sum = T{};
    for (size_t i = 0; i < window_size; ++i) {
        window_sum += values[i];
    }
    result.push_back(window_sum / static_cast<T>(window_size));

    // Slide the window: remove leftmost, add new rightmost
    for (size_t i = window_size; i < values.size(); ++i) {
        window_sum = window_sum - values[i - window_size] + values[i];
        result.push_back(window_sum / static_cast<T>(window_size));
    }

    return result;
}
```

**Why Template the Sliding Window?**

By templating `sliding_window_average`, it works with any numeric type:
- `int`, `float`, `double` for sensor readings
- Custom numeric types (if they support `+`, `-`, `/`)
- Maintains type safety (no implicit conversions)

**Sliding Window Algorithm Explanation:**

Instead of recalculating the sum for each window (O(n*w) complexity), we use a sliding approach:
1. Calculate sum for first window
2. For each subsequent window, subtract the element leaving and add the element entering
3. This reduces complexity to O(n)

**Usage with different containers:**
```cpp
std::array<double, 5> arr = {1, 2, 3, 4, 5};
std::vector<double> vec = {1, 2, 3, 4, 5};
double c_arr[] = {1, 2, 3, 4, 5};

// All work with the same function!
auto avg1 = average(arr);
auto avg2 = average(vec);
auto avg3 = average(std::span{c_arr});

// Sliding window examples
std::array<double, 7> noisy = {1.0, 5.0, 2.0, 6.0, 3.0, 7.0, 4.0};
auto smoothed = sliding_window_average(std::span<double const>{noisy}, 3);
// smoothed = {2.67, 4.33, 3.67, 5.33, 4.67}

// Works with integers too!
std::vector<int> int_data = {10, 20, 30, 40, 50};
auto int_averaged = sliding_window_average(std::span<int const>{int_data}, 2);
// int_averaged = {15, 25, 35, 45}

// Subspans (slicing)
auto first_half = std::span{arr}.subspan(0, 3);  // {1, 2, 3}
auto avg_half = average(first_half);
```

**Modern C++20 ranges alternative:**
```cpp
#include <ranges>

double average(std::span<double const> const readings) {
    if (readings.empty()) return 0.0;

    auto const sum = std::ranges::fold_left(readings, 0.0, std::plus{});
    return sum / static_cast<double>(readings.size());
}

size_t count_above_threshold(std::span<double const> const readings,
                             double const threshold) {
    return std::ranges::count_if(readings,
        [threshold](double v) { return v > threshold; });
}
```

---

## Key Takeaways

### When to use `std::array`

✅ **Use when:**
- Size is known at compile time
- Want stack allocation (performance)
- Need compile-time size guarantees
- Replacing C-style arrays

❌ **Don't use when:**
- Size varies at runtime → use `std::vector`
- Need frequent insertions/deletions → use `std::deque` or `std::list`
- Size is very large → might overflow stack, use `std::vector`

### When to use `std::span`

✅ **Use when:**
- Function needs to work with multiple container types
- Want to avoid copying large containers
- Need a view/slice of contiguous data
- Interfacing with C APIs (replaces `T*, size_t` pairs)

❌ **Don't use when:**
- Need ownership of data → use `std::vector` or `std::array`
- Data isn't contiguous → use iterators or ranges
- Lifetime management is complex → be very careful with dangling

### Performance Characteristics

| Operation | `std::array` | `std::vector` | `std::span` |
|-----------|--------------|---------------|-------------|
| Construction | O(1) | O(n) | O(1) |
| Random access | O(1) | O(1) | O(1) |
| Size | Compile-time | Runtime | Runtime |
| Storage | Stack | Heap | No storage (view) |
| Reallocation | Never | Possible | N/A (view) |

### Memory Layout

```cpp
std::array<int, 3> arr = {1, 2, 3};
// Stack: [1][2][3]

std::vector<int> vec = {1, 2, 3};
// Stack: [ptr][size][capacity]
// Heap:  [1][2][3]

std::span<int> sp = arr;
// Stack: [ptr][size]
// Points to: arr's stack memory
```

### Common Pitfalls

**1. Dangling span:**
```cpp
std::span<int> get_data() {
    std::vector<int> temp = {1, 2, 3};
    return temp;  // DANGER: span outlives vector!
}
```

**2. Modifying const data:**
```cpp
std::array<int, 3> const arr = {1, 2, 3};
std::span<int> sp = arr;  // Compile error! Need std::span<int const>
```

**3. Fixed extent vs dynamic extent:**
```cpp
std::span<int, 3> fixed = arr;     // Size known at compile time
std::span<int> dynamic = arr;      // Size known at runtime
std::span<int, 3> wrong = vec;     // Error if vec.size() != 3
```

---

## Discussion Questions & Answers

**Q1: When should you use `std::array` instead of `std::vector`?**

A: Use `std::array` when:
- Size is known at compile time (e.g., 3D coordinates, RGB colors, fixed joints)
- Performance is critical (stack allocation is faster than heap)
- You want compile-time size checking
- Working in embedded/constrained environments

**Q2: What are the lifetime considerations when using `std::span`?**

A: `std::span` is a non-owning view. It's your responsibility to ensure:
- The underlying data outlives the span
- The data isn't reallocated (e.g., `vector::push_back` can invalidate)
- You don't return spans to local variables

**Q3: How does `std::span` improve on the old pattern of passing `T* ptr, size_t length`?**

A: `std::span` provides:
- Type safety (single object, no separate pointer and size)
- Iterator interface (works with algorithms)
- Bounds checking in debug mode
- Const-correctness (`std::span<T const>`)
- Subspan operations (slicing)
- Better compiler optimization opportunities

**Q4: Why does `std::span` support both dynamic and fixed extents?**

A: Fixed extent `std::span<T, N>` allows:
- Compile-time size checks
- Potential optimizations (no runtime size storage)
- API contracts (function requires exactly N elements)

Dynamic extent `std::span<T>` allows:
- Flexibility with variable-sized data
- Works like a "view" of any sized container

**Q5: What are the trade-offs between `std::array<T, N>` and C-style `T[N]` arrays?**

A: `std::array` advantages:
- No decay to pointer
- Size accessible via `.size()`
- Works with STL algorithms
- Copy/assignment support
- Structured bindings

C-array advantages:
- Slightly simpler syntax
- Guaranteed aggregate initialization
- (Very minor) Less verbose in some contexts

**Q6: How can `std::span` help with interfacing to C libraries?**

A: Many C APIs use `void* data, size_t size` or `T* array, int count` patterns:

```cpp
// C API
void process_data(double const* data, size_t count);

// Before: manual pointer + size
std::vector<double> vec = {1, 2, 3};
process_data(vec.data(), vec.size());

// With span: cleaner wrapper
void process_data_cpp(std::span<double const> data) {
    process_data(data.data(), data.size());
}

process_data_cpp(vec);  // Works with vector
process_data_cpp(arr);  // Works with array
```
