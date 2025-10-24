// Exercise 1: 3D Transformation Library with Strong Types
// This exercise demonstrates strong types, encapsulation, const-correctness,
// rule of 5, user-defined literals, and mathematical operator overloading

#include <array>
#include <cmath>
#include <numbers>
#include <gtest/gtest.h>

// ===== Strong Types and User-Defined Literals =====

namespace literals {
    // Strong type for distances (meters)
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

        // Comparison operators
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

    // Conversion functions (work with both types)
    [[nodiscard]] constexpr radian_t to_radians(degree_t const d) {
        return radian_t{d.value * std::numbers::pi / 180.0};
    }

    [[nodiscard]] constexpr radian_t to_radians(radian_t const r) {
        return r;  // No conversion needed
    }

    [[nodiscard]] constexpr degree_t to_degrees(radian_t const r) {
        return degree_t{r.value * 180.0 / std::numbers::pi};
    }

    [[nodiscard]] constexpr degree_t to_degrees(degree_t const d) {
        return d;  // No conversion needed
    }

    // User-defined literals
    [[nodiscard]] constexpr meter_t operator""_m(long double v) {
        return meter_t{static_cast<double>(v)};
    }

    [[nodiscard]] constexpr meter_t operator""_m(unsigned long long v) {
        return meter_t{static_cast<double>(v)};
    }

    [[nodiscard]] constexpr radian_t operator""_rad(long double v) {
        return radian_t{static_cast<double>(v)};
    }

    [[nodiscard]] constexpr radian_t operator""_rad(unsigned long long v) {
        return radian_t{static_cast<double>(v)};
    }

    [[nodiscard]] constexpr degree_t operator""_deg(long double v) {
        return degree_t{static_cast<double>(v)};
    }

    [[nodiscard]] constexpr degree_t operator""_deg(unsigned long long v) {
        return degree_t{static_cast<double>(v)};
    }
}

// ===== Position Class =====

/**
 * @brief Represents a 3D position_t in space.
 *
 * Stores coordinates privately as std::array and provides accessor methods
 * with const and non-const overloads to demonstrate const-correctness.
 */
class position_t {
private:
    std::array<literals::meter_t, 3> coords_;

public:
    // Default constructor: origin
    position_t() : coords_{literals::meter_t{0.0}, literals::meter_t{0.0}, literals::meter_t{0.0}} {}

    // Construct from strong types (meters only!)
    position_t(literals::meter_t const x, literals::meter_t const y, literals::meter_t const z)
        : coords_{x, y, z} {}

    // Accessors (const overloads return const reference)
    [[nodiscard]] literals::meter_t const& x() const { return coords_[0]; }
    [[nodiscard]] literals::meter_t const& y() const { return coords_[1]; }
    [[nodiscard]] literals::meter_t const& z() const { return coords_[2]; }

    // Mutators (non-const overloads return non-const reference)
    [[nodiscard]] literals::meter_t& x() { return coords_[0]; }
    [[nodiscard]] literals::meter_t& y() { return coords_[1]; }
    [[nodiscard]] literals::meter_t& z() { return coords_[2]; }

    // Equality comparison
    [[nodiscard]] bool operator==(position_t const& other) const {
        return coords_ == other.coords_;
    }

    // Vector addition
    [[nodiscard]] position_t operator+(position_t const& other) const {
        return position_t{
            literals::meter_t{x().value + other.x().value},
            literals::meter_t{y().value + other.y().value},
            literals::meter_t{z().value + other.z().value}
        };
    }

    // Vector subtraction
    [[nodiscard]] position_t operator-(position_t const& other) const {
        return position_t{
            literals::meter_t{x().value - other.x().value},
            literals::meter_t{y().value - other.y().value},
            literals::meter_t{z().value - other.z().value}
        };
    }

    // Scalar multiplication
    [[nodiscard]] position_t operator*(double const scalar) const {
        return position_t{
            literals::meter_t{x().value * scalar},
            literals::meter_t{y().value * scalar},
            literals::meter_t{z().value * scalar}
        };
    }
};

// ===== Position Free Functions =====

/**
 * @brief Calculate distance between two positions.
 *
 * @param p1 First position
 * @param p2 Second position
 * @return Distance in meters (strong type)
 */
[[nodiscard]] literals::meter_t distance(position_t const& p1, position_t const& p2) {
    double const dx = p1.x().value - p2.x().value;
    double const dy = p1.y().value - p2.y().value;
    double const dz = p1.z().value - p2.z().value;
    return literals::meter_t{std::hypot(dx, dy, dz)};
}

/**
 * @brief Check if two positions are approximately equal within a tolerance.
 *
 * @param p1 First position
 * @param p2 Second position
 * @param tolerance Maximum allowed difference per coordinate (default: 0.001m)
 * @return true if positions are within tolerance
 */
[[nodiscard]] bool near(position_t const& p1, position_t const& p2,
                        literals::meter_t const tolerance = literals::meter_t{0.001}) {
    return std::abs(p1.x().value - p2.x().value) <= tolerance.value and
           std::abs(p1.y().value - p2.y().value) <= tolerance.value and
           std::abs(p1.z().value - p2.z().value) <= tolerance.value;
}

// ===== Quaternion Class =====

/**
 * @brief Represents a rotation as a quaternion_t.
 *
 * Can be constructed from quaternion_t components (x, y, z, w) or Euler angles
 * (roll, pitch, yaw). Demonstrates rule of 5 and multiple constructors.
 */
class quaternion_t {
private:
    std::array<double, 4> components_;  // x, y, z, w

    // Helper: normalize the quaternion
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

    // Construct from quaternion_t components (x, y, z, w)
    quaternion_t(double const x, double const y, double const z, double const w)
        : components_{x, y, z, w} {
        normalize();
    }

    // Named constructor: from Euler angles (accepts ONLY radian_t type!)
    [[nodiscard]] static quaternion_t from_euler(literals::radian_t const roll,
                                               literals::radian_t const pitch,
                                               literals::radian_t const yaw) {
        // Conversion from Euler angles to quaternion
        double const cy = std::cos(yaw.value * 0.5);
        double const sy = std::sin(yaw.value * 0.5);
        double const cp = std::cos(pitch.value * 0.5);
        double const sp = std::sin(pitch.value * 0.5);
        double const cr = std::cos(roll.value * 0.5);
        double const sr = std::sin(roll.value * 0.5);

        return quaternion_t{
            sr * cp * cy - cr * sp * sy,  // x
            cr * sp * cy + sr * cp * sy,  // y
            cr * cp * sy - sr * sp * cy,  // z
            cr * cp * cy + sr * sp * sy   // w
        };
    }

    // Copy constructor
    quaternion_t(quaternion_t const& other) = default;

    // Move constructor
    quaternion_t(quaternion_t&& other) noexcept = default;

    // Copy assignment
    quaternion_t& operator=(quaternion_t const& other) = default;

    // Move assignment
    quaternion_t& operator=(quaternion_t&& other) noexcept = default;

    // Destructor
    ~quaternion_t() = default;

    // Accessors
    [[nodiscard]] double const& x() const { return components_[0]; }
    [[nodiscard]] double const& y() const { return components_[1]; }
    [[nodiscard]] double const& z() const { return components_[2]; }
    [[nodiscard]] double const& w() const { return components_[3]; }

    // Mutators
    [[nodiscard]] double& x() { return components_[0]; }
    [[nodiscard]] double& y() { return components_[1]; }
    [[nodiscard]] double& z() { return components_[2]; }
    [[nodiscard]] double& w() { return components_[3]; }

    // Equality comparison
    [[nodiscard]] bool operator==(quaternion_t const& other) const {
        return components_ == other.components_;
    }

    // Approximate equality
    [[nodiscard]] bool approx_equal(quaternion_t const& other, double const tolerance = 0.001) const {
        return std::abs(x() - other.x()) <= tolerance and
               std::abs(y() - other.y()) <= tolerance and
               std::abs(z() - other.z()) <= tolerance and
               std::abs(w() - other.w()) <= tolerance;
    }

    // Quaternion multiplication (composition of rotations)
    [[nodiscard]] quaternion_t operator*(quaternion_t const& other) const {
        return quaternion_t{
            w() * other.x() + x() * other.w() + y() * other.z() - z() * other.y(),
            w() * other.y() - x() * other.z() + y() * other.w() + z() * other.x(),
            w() * other.z() + x() * other.y() - y() * other.x() + z() * other.w(),
            w() * other.w() - x() * other.x() - y() * other.y() - z() * other.z()
        };
    }

    // Conjugate (inverse rotation for unit quaternions)
    [[nodiscard]] quaternion_t conjugate() const {
        return quaternion_t{-x(), -y(), -z(), w()};
    }
};

// ===== Transformation Class =====

/**
 * @brief Represents a 3D transformation_t (position_t + rotation).
 *
 * Stores a full 4x4 homogeneous transformation_t matrix as std::array<double, 16>.
 * Can only be constructed with both position_t and quaternion_t (no default constructor).
 * Matrix is stored in row-major order.
 */
class transformation_t {
private:
    std::array<double, 16> matrix_;  // 4x4 matrix in row-major order

    // Helper: Convert quaternion_t to 3x3 rotation matrix and build 4x4 transform
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

        // Build 4x4 matrix: [R | t]
        //                   [0 | 1]
        // Row-major order: [m00, m01, m02, m03, m10, m11, m12, m13, ...]
        return {
            1 - 2*(yy + zz), 2*(xy - wz),     2*(xz + wy),     pos.x().value,  // Row 0
            2*(xy + wz),     1 - 2*(xx + zz), 2*(yz - wx),     pos.y().value,  // Row 1
            2*(xz - wy),     2*(yz + wx),     1 - 2*(xx + yy), pos.z().value,  // Row 2
            0.0,             0.0,             0.0,             1.0             // Row 3
        };
    }

public:
    // Construct from position_t and quaternion
    transformation_t(position_t const& pos, quaternion_t const& rot)
        : matrix_{build_matrix(pos, rot)} {}

    // Copy constructor
    transformation_t(transformation_t const& other) = default;

    // Move constructor
    transformation_t(transformation_t&& other) noexcept = default;

    // Copy assignment
    transformation_t& operator=(transformation_t const& other) = default;

    // Move assignment
    transformation_t& operator=(transformation_t&& other) noexcept = default;

    // Destructor
    ~transformation_t() = default;

    // Extract position_t (from last column of matrix)
    [[nodiscard]] position_t get_position() const {
        return position_t{
            literals::meter_t{matrix_[3]},   // m03
            literals::meter_t{matrix_[7]},   // m13
            literals::meter_t{matrix_[11]}   // m23
        };
    }

    // Extract quaternion_t (from rotation part of matrix)
    [[nodiscard]] quaternion_t get_rotation() const {
        // Extract 3x3 rotation matrix elements
        double const m00 = matrix_[0], m01 = matrix_[1], m02 = matrix_[2];
        double const m10 = matrix_[4], m11 = matrix_[5], m12 = matrix_[6];
        double const m20 = matrix_[8], m21 = matrix_[9], m22 = matrix_[10];

        // Convert rotation matrix to quaternion_t (standard algorithm)
        double const trace = m00 + m11 + m22;

        if (trace > 0) {
            double const s = 0.5 / std::sqrt(trace + 1.0);
            return quaternion_t{
                (m21 - m12) * s,
                (m02 - m20) * s,
                (m10 - m01) * s,
                0.25 / s
            };
        } else if (m00 > m11 && m00 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m00 - m11 - m22);
            return quaternion_t{
                0.25 * s,
                (m01 + m10) / s,
                (m02 + m20) / s,
                (m21 - m12) / s
            };
        } else if (m11 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m11 - m00 - m22);
            return quaternion_t{
                (m01 + m10) / s,
                0.25 * s,
                (m12 + m21) / s,
                (m02 - m20) / s
            };
        } else {
            double const s = 2.0 * std::sqrt(1.0 + m22 - m00 - m11);
            return quaternion_t{
                (m02 + m20) / s,
                (m12 + m21) / s,
                0.25 * s,
                (m10 - m01) / s
            };
        }
    }

    // Equality comparison
    [[nodiscard]] bool operator==(transformation_t const& other) const {
        return matrix_ == other.matrix_;
    }

    // Approximate equality (uses strong type for position_t tolerance)
    [[nodiscard]] bool approx_equal(transformation_t const& other,
                                    literals::meter_t const pos_tolerance = literals::meter_t{0.001},
                                    double const rot_tolerance = 0.001) const {
        return near(get_position(), other.get_position(), pos_tolerance) and
               get_rotation().approx_equal(other.get_rotation(), rot_tolerance);
    }

    // Transform composition (4x4 matrix multiplication)
    [[nodiscard]] transformation_t operator*(transformation_t const& other) const {
        std::array<double, 16> result{};

        // Multiply: result = this->matrix_ * other.matrix_
        // Row-major order: result[row*4 + col]
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += matrix_[row * 4 + k] * other.matrix_[k * 4 + col];
                }
                result[row * 4 + col] = sum;
            }
        }

        // Create transformation_t from resulting matrix
        // Extract position_t and rotation from result matrix
        transformation_t temp{position_t{}, quaternion_t{}};
        temp.matrix_ = result;
        return temp;
    }

    // Transform a position_t (apply 4x4 matrix to homogeneous coordinate)
    [[nodiscard]] position_t transform_point(position_t const& point) const {
        // Treat point as [x, y, z, 1] in homogeneous coordinates
        double const x = point.x().value;
        double const y = point.y().value;
        double const z = point.z().value;

        // Matrix-vector multiplication: result = matrix_ * [x, y, z, 1]
        return position_t{
            literals::meter_t{matrix_[0]*x + matrix_[1]*y + matrix_[2]*z + matrix_[3]},
            literals::meter_t{matrix_[4]*x + matrix_[5]*y + matrix_[6]*z + matrix_[7]},
            literals::meter_t{matrix_[8]*x + matrix_[9]*y + matrix_[10]*z + matrix_[11]}
        };
    }

    // Inverse transformation_t (invert 4x4 matrix)
    [[nodiscard]] transformation_t inverse() const {
        // For SE(3) transformation_t, inverse is:
        // [R^T | -R^T*t]
        // [0   | 1     ]

        // Extract rotation (upper-left 3x3)
        double const r00 = matrix_[0], r01 = matrix_[1], r02 = matrix_[2];
        double const r10 = matrix_[4], r11 = matrix_[5], r12 = matrix_[6];
        double const r20 = matrix_[8], r21 = matrix_[9], r22 = matrix_[10];

        // Extract translation
        double const tx = matrix_[3];
        double const ty = matrix_[7];
        double const tz = matrix_[11];

        // Transpose rotation
        // Compute -R^T * t
        double const new_tx = -(r00*tx + r10*ty + r20*tz);
        double const new_ty = -(r01*tx + r11*ty + r21*tz);
        double const new_tz = -(r02*tx + r12*ty + r22*tz);

        std::array<double, 16> inv_matrix = {
            r00, r10, r20, new_tx,   // Row 0: R^T[0] | -R^T*t[0]
            r01, r11, r21, new_ty,   // Row 1: R^T[1] | -R^T*t[1]
            r02, r12, r22, new_tz,   // Row 2: R^T[2] | -R^T*t[2]
            0.0, 0.0, 0.0, 1.0       // Row 3
        };

        transformation_t result{position_t{}, quaternion_t{}};
        result.matrix_ = inv_matrix;
        return result;
    }
};

// ===== Tests =====

using namespace literals;

TEST(PositionTest, DefaultConstructor) {
    position_t const p;
    EXPECT_DOUBLE_EQ(p.x().value, 0.0);
    EXPECT_DOUBLE_EQ(p.y().value, 0.0);
    EXPECT_DOUBLE_EQ(p.z().value, 0.0);
}

TEST(PositionTest, ParameterizedConstructorWithStrongTypes) {
    position_t const p{1.0_m, 2.0_m, 3.0_m};
    EXPECT_DOUBLE_EQ(p.x().value, 1.0);
    EXPECT_DOUBLE_EQ(p.y().value, 2.0);
    EXPECT_DOUBLE_EQ(p.z().value, 3.0);
}

TEST(PositionTest, ConstAccessors) {
    position_t const p{1.0_m, 2.0_m, 3.0_m};
    // These should compile (const accessors return meter_t const&)
    [[maybe_unused]] literals::meter_t const& x = p.x();
    [[maybe_unused]] literals::meter_t const& y = p.y();
    [[maybe_unused]] literals::meter_t const& z = p.z();
}

TEST(PositionTest, Mutators) {
    position_t p{1.0_m, 2.0_m, 3.0_m};
    p.x() = 10.0_m;
    p.y() = 20.0_m;
    p.z() = 30.0_m;

    EXPECT_DOUBLE_EQ(p.x().value, 10.0);
    EXPECT_DOUBLE_EQ(p.y().value, 20.0);
    EXPECT_DOUBLE_EQ(p.z().value, 30.0);
}

TEST(PositionTest, UserDefinedLiterals) {
    position_t const p{1.0_m, 2.0_m, 3.0_m};
    EXPECT_DOUBLE_EQ(p.x().value, 1.0);
    EXPECT_DOUBLE_EQ(p.y().value, 2.0);
    EXPECT_DOUBLE_EQ(p.z().value, 3.0);
}

TEST(PositionTest, Equality) {
    position_t const p1{1.0_m, 2.0_m, 3.0_m};
    position_t const p2{1.0_m, 2.0_m, 3.0_m};
    position_t const p3{1.1_m, 2.0_m, 3.0_m};

    EXPECT_TRUE(p1 == p2);
    EXPECT_FALSE(p1 == p3);
}

TEST(PositionTest, ApproximateEquality) {
    position_t const p1{1.0_m, 2.0_m, 3.0_m};
    position_t const p2{1.0001_m, 2.0001_m, 3.0001_m};
    position_t const p3{1.1_m, 2.0_m, 3.0_m};

    EXPECT_TRUE(near(p1, p2, 0.001_m));
    EXPECT_FALSE(near(p1, p3, 0.001_m));
}

TEST(PositionTest, DistanceReturnsStrongType) {
    position_t const origin{0.0_m, 0.0_m, 0.0_m};
    position_t const p1{3.0_m, 4.0_m, 0.0_m};
    position_t const p2{1.0_m, 2.0_m, 2.0_m};

    auto const dist1 = distance(origin, p1);
    auto const dist2 = distance(origin, p2);

    EXPECT_NEAR(dist1.value, 5.0, 0.001);  // 3-4-5 triangle
    EXPECT_NEAR(dist2.value, 3.0, 0.001);  // sqrt(1+4+4)
}

TEST(PositionTest, VectorAddition) {
    position_t const p1{1.0_m, 2.0_m, 3.0_m};
    position_t const p2{4.0_m, 5.0_m, 6.0_m};
    auto const result = p1 + p2;

    EXPECT_DOUBLE_EQ(result.x().value, 5.0);
    EXPECT_DOUBLE_EQ(result.y().value, 7.0);
    EXPECT_DOUBLE_EQ(result.z().value, 9.0);
}

TEST(PositionTest, VectorSubtraction) {
    position_t const p1{5.0_m, 7.0_m, 9.0_m};
    position_t const p2{1.0_m, 2.0_m, 3.0_m};
    auto const result = p1 - p2;

    EXPECT_DOUBLE_EQ(result.x().value, 4.0);
    EXPECT_DOUBLE_EQ(result.y().value, 5.0);
    EXPECT_DOUBLE_EQ(result.z().value, 6.0);
}

TEST(PositionTest, ScalarMultiplication) {
    position_t const p{1.0_m, 2.0_m, 3.0_m};
    auto const result = p * 2.0;

    EXPECT_DOUBLE_EQ(result.x().value, 2.0);
    EXPECT_DOUBLE_EQ(result.y().value, 4.0);
    EXPECT_DOUBLE_EQ(result.z().value, 6.0);
}

TEST(QuaternionTest, DefaultConstructor) {
    quaternion_t const q;
    EXPECT_DOUBLE_EQ(q.x(), 0.0);
    EXPECT_DOUBLE_EQ(q.y(), 0.0);
    EXPECT_DOUBLE_EQ(q.z(), 0.0);
    EXPECT_DOUBLE_EQ(q.w(), 1.0);  // Identity
}

TEST(QuaternionTest, ComponentConstructor) {
    quaternion_t const q{0.0, 0.0, 0.0, 1.0};
    EXPECT_DOUBLE_EQ(q.x(), 0.0);
    EXPECT_DOUBLE_EQ(q.y(), 0.0);
    EXPECT_DOUBLE_EQ(q.z(), 0.0);
    EXPECT_DOUBLE_EQ(q.w(), 1.0);
}

TEST(QuaternionTest, NormalizationInConstructor) {
    quaternion_t const q{1.0, 1.0, 1.0, 1.0};
    // Should be normalized to unit quaternion
    double const magnitude = std::sqrt(q.x()*q.x() + q.y()*q.y() + q.z()*q.z() + q.w()*q.w());
    EXPECT_NEAR(magnitude, 1.0, 0.001);
}

TEST(QuaternionTest, FromEulerAnglesUsesRadianType) {
    auto const q = quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.0_rad);
    EXPECT_NEAR(q.x(), 0.0, 0.001);
    EXPECT_NEAR(q.y(), 0.0, 0.001);
    EXPECT_NEAR(q.z(), 0.0, 0.001);
    EXPECT_NEAR(q.w(), 1.0, 0.001);
}

TEST(QuaternionTest, FromEulerAnglesWithDegreeConversion) {
    // 90 degree_t rotation around Z axis (must convert to radians!)
    degree_t const deg = 90.0_deg;
    auto const q = quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(deg));
    EXPECT_NEAR(q.z(), 0.707, 0.01);  // sin(45°)
    EXPECT_NEAR(q.w(), 0.707, 0.01);  // cos(45°)
}

TEST(QuaternionTest, UserDefinedLiteralsRadians) {
    auto const angle = 1.57_rad;
    EXPECT_NEAR(angle.value, 1.57, 0.001);
}

TEST(QuaternionTest, UserDefinedLiteralsDegrees) {
    auto const angle = 90.0_deg;
    EXPECT_NEAR(angle.value, 90.0, 0.001);  // Stored as degrees
}

TEST(QuaternionTest, DegreeToRadianConversion) {
    auto const deg = 180.0_deg;
    auto const rad = to_radians(deg);
    EXPECT_NEAR(rad.value, std::numbers::pi, 0.001);
}

TEST(QuaternionTest, RadianToDegreeConversion) {
    auto const rad = 3.14159_rad;
    auto const deg = to_degrees(rad);
    EXPECT_NEAR(deg.value, 180.0, 0.1);  // Stored as degrees
}

TEST(QuaternionTest, ConversionFunctionsIdempotent) {
    // Converting same type returns same value
    auto const rad1 = 1.57_rad;
    auto const rad2 = to_radians(rad1);
    EXPECT_DOUBLE_EQ(rad1.value, rad2.value);

    auto const deg1 = 90.0_deg;
    auto const deg2 = to_degrees(deg1);
    EXPECT_DOUBLE_EQ(deg1.value, deg2.value);
}

TEST(QuaternionTest, CopyConstructor) {
    quaternion_t const q1{1.0, 2.0, 3.0, 4.0};
    quaternion_t const q2{q1};

    EXPECT_DOUBLE_EQ(q1.x(), q2.x());
    EXPECT_DOUBLE_EQ(q1.y(), q2.y());
    EXPECT_DOUBLE_EQ(q1.z(), q2.z());
    EXPECT_DOUBLE_EQ(q1.w(), q2.w());
}

TEST(QuaternionTest, MoveConstructor) {
    quaternion_t q1{1.0, 2.0, 3.0, 4.0};
    quaternion_t const q2{std::move(q1)};

    // q2 should have the values (normalized)
    EXPECT_GT(q2.w(), 0.0);
}

TEST(QuaternionTest, CopyAssignment) {
    quaternion_t const q1{1.0, 2.0, 3.0, 4.0};
    quaternion_t q2;
    q2 = q1;

    EXPECT_DOUBLE_EQ(q1.x(), q2.x());
    EXPECT_DOUBLE_EQ(q1.w(), q2.w());
}

TEST(QuaternionTest, Equality) {
    quaternion_t const q1{0.0, 0.0, 0.0, 1.0};
    quaternion_t const q2{0.0, 0.0, 0.0, 1.0};
    quaternion_t const q3{1.0, 0.0, 0.0, 1.0};

    EXPECT_TRUE(q1 == q2);
    EXPECT_FALSE(q1 == q3);
}

TEST(QuaternionTest, Multiplication) {
    quaternion_t const q1{0.0, 0.0, 0.0, 1.0};  // Identity
    quaternion_t const q2{1.0, 0.0, 0.0, 1.0};
    auto const result = q1 * q2;

    // Identity * q = q (after normalization)
    EXPECT_TRUE(result.approx_equal(q2, 0.01));
}

TEST(QuaternionTest, Conjugate) {
    quaternion_t const q{1.0, 2.0, 3.0, 4.0};
    auto const conj = q.conjugate();

    EXPECT_NEAR(conj.x(), -q.x(), 0.001);
    EXPECT_NEAR(conj.y(), -q.y(), 0.001);
    EXPECT_NEAR(conj.z(), -q.z(), 0.001);
    EXPECT_NEAR(conj.w(), q.w(), 0.001);
}

TEST(TransformationTest, ConstructorRequiresPositionAndRotation) {
    position_t const pos{1.0_m, 2.0_m, 3.0_m};
    quaternion_t const rot{0.0, 0.0, 0.0, 1.0};

    transformation_t const tf{pos, rot};

    EXPECT_TRUE(tf.get_position() == pos);
    EXPECT_TRUE(tf.get_rotation() == rot);
}

TEST(TransformationTest, CopyConstructor) {
    position_t const pos{1.0_m, 2.0_m, 3.0_m};
    quaternion_t const rot{0.0, 0.0, 0.0, 1.0};
    transformation_t const tf1{pos, rot};
    transformation_t const tf2{tf1};

    EXPECT_TRUE(tf1 == tf2);
}

TEST(TransformationTest, Equality) {
    position_t const pos{1.0_m, 2.0_m, 3.0_m};
    quaternion_t const rot{0.0, 0.0, 0.0, 1.0};
    transformation_t const tf1{pos, rot};
    transformation_t const tf2{pos, rot};

    EXPECT_TRUE(tf1 == tf2);
}

TEST(TransformationTest, ApproximateEquality) {
    transformation_t const tf1{position_t{1.0_m, 2.0_m, 3.0_m}, quaternion_t{0.0, 0.0, 0.0, 1.0}};
    transformation_t const tf2{position_t{1.0001_m, 2.0001_m, 3.0001_m}, quaternion_t{0.0, 0.0, 0.0, 1.0}};

    EXPECT_TRUE(tf1.approx_equal(tf2, 0.001_m));
}

TEST(TransformationTest, TransformComposition) {
    transformation_t const tf1{position_t{1.0_m, 0.0_m, 0.0_m}, quaternion_t{}};
    transformation_t const tf2{position_t{0.0_m, 1.0_m, 0.0_m}, quaternion_t{}};

    auto const composed = tf1 * tf2;

    // Result should have combined translation
    auto const result_pos = composed.get_position();
    EXPECT_NEAR(result_pos.x().value, 1.0, 0.1);
    EXPECT_NEAR(result_pos.y().value, 1.0, 0.1);
}

TEST(TransformationTest, TransformPoint) {
    transformation_t const tf{position_t{1.0_m, 2.0_m, 3.0_m}, quaternion_t{}};
    position_t const point{0.0_m, 0.0_m, 0.0_m};

    auto const transformed = tf.transform_point(point);

    EXPECT_NEAR(transformed.x().value, 1.0, 0.001);
    EXPECT_NEAR(transformed.y().value, 2.0, 0.001);
    EXPECT_NEAR(transformed.z().value, 3.0, 0.001);
}

TEST(TransformationTest, InverseTransformation) {
    transformation_t const tf{position_t{1.0_m, 2.0_m, 3.0_m}, quaternion_t{}};
    auto const tf_inv = tf.inverse();

    auto const inv_pos = tf_inv.get_position();
    EXPECT_NEAR(inv_pos.x().value, -1.0, 0.001);
    EXPECT_NEAR(inv_pos.y().value, -2.0, 0.001);
    EXPECT_NEAR(inv_pos.z().value, -3.0, 0.001);
}

TEST(IntegrationTest, RobotArmKinematicsWithStrongTypes) {
    using namespace literals;

    // Base of robot at origin
    transformation_t const base{position_t{0.0_m, 0.0_m, 0.0_m}, quaternion_t{}};

    // First joint: 1m up
    transformation_t const joint1{position_t{0.0_m, 0.0_m, 1.0_m}, quaternion_t{}};

    // Second joint: 0.5m forward, rotated 45 degrees
    // Note: must explicitly convert degrees to radians!
    degree_t const angle = 45.0_deg;
    transformation_t const joint2{
        position_t{0.5_m, 0.0_m, 0.0_m},
        quaternion_t::from_euler(0.0_rad, 0.0_rad, to_radians(angle))
    };

    // Compose transformations
    auto const end_effector = base * joint1 * joint2;

    // End effector should be approximately at (0.5, 0, 1)
    auto const pos = end_effector.get_position();
    EXPECT_NEAR(pos.z().value, 1.0, 0.1);
}

TEST(UserDefinedLiteralsTest, PositionAndQuaternionLiterals) {
    using namespace literals;

    // Test position_t literal (creates origin)
    auto const origin = 0.0_pos;
    EXPECT_DOUBLE_EQ(origin.x().value, 0.0);
    EXPECT_DOUBLE_EQ(origin.y().value, 0.0);
    EXPECT_DOUBLE_EQ(origin.z().value, 0.0);

    // Test quaternion_t literal (creates identity)
    auto const identity = 0.0_quat;
    EXPECT_DOUBLE_EQ(identity.x(), 0.0);
    EXPECT_DOUBLE_EQ(identity.y(), 0.0);
    EXPECT_DOUBLE_EQ(identity.z(), 0.0);
    EXPECT_DOUBLE_EQ(identity.w(), 1.0);

    // Can use in expressions
    auto const dist_from_origin = distance(origin, position_t{1.0_m, 0.0_m, 0.0_m});
    EXPECT_NEAR(dist_from_origin.value, 1.0, 0.001);
}

TEST(StrongTypesTest, PreventImplicitConversions) {
    // This test verifies compile-time type safety
    using namespace literals;

    // These should compile:
    [[maybe_unused]] position_t const p1{1.0_m, 2.0_m, 3.0_m};

    // Explicit degree_t to radian_t conversion required (both methods work)
    degree_t const deg_angle = 90.0_deg;
    radian_t const rad_angle = to_radians(deg_angle);
    [[maybe_unused]] auto const q1 = quaternion_t::from_euler(0.0_rad, 0.0_rad, rad_angle);

    // Can also convert radians to degrees
    radian_t const rad = 1.57_rad;
    [[maybe_unused]] degree_t const deg = to_degrees(rad);

    // These should NOT compile (uncomment to verify):
    // position_t const p2{1.0, 2.0, 3.0};  // ERROR: no conversion from double to meter
    // auto const q2 = quaternion_t::from_euler(0.0, 0.0, 1.57);  // ERROR: needs radian_t type
    // auto const q3 = quaternion_t::from_euler(0.0_deg, 0.0_deg, 90.0_deg);  // ERROR: needs radian, not degree

    SUCCEED();  // If we get here, the type system is working correctly
}
