// Exercise 2: Template Specialization for Rotation Averaging
// This exercise demonstrates template specialization with type-appropriate averaging algorithms

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <numbers>
#include <span>

// ============================================================================
// Strong Types and Classes from Exercise 1
// ============================================================================

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

constexpr meter_t operator""_m(long double v) {
    return meter_t{static_cast<double>(v)};
}
constexpr meter_t operator""_m(unsigned long long v) {
    return meter_t{static_cast<double>(v)};
}
constexpr radian_t operator""_rad(long double v) {
    return radian_t{static_cast<double>(v)};
}
constexpr radian_t operator""_rad(unsigned long long v) {
    return radian_t{static_cast<double>(v)};
}
}  // namespace literals

using literals::meter_t;
using literals::radian_t;

// position_t class
class position_t {
   private:
    std::array<meter_t, 3> coords_;

   public:
    position_t() : coords_{meter_t{0.0}, meter_t{0.0}, meter_t{0.0}} {}
    position_t(meter_t const x, meter_t const y, meter_t const z) : coords_{x, y, z} {}

    [[nodiscard]] meter_t const& x() const { return coords_[0]; }
    [[nodiscard]] meter_t const& y() const { return coords_[1]; }
    [[nodiscard]] meter_t const& z() const { return coords_[2]; }
    [[nodiscard]] meter_t& x() { return coords_[0]; }
    [[nodiscard]] meter_t& y() { return coords_[1]; }
    [[nodiscard]] meter_t& z() { return coords_[2]; }

    [[nodiscard]] position_t operator+(position_t const& other) const {
        return position_t{x() + other.x(), y() + other.y(), z() + other.z()};
    }
    [[nodiscard]] position_t operator-(position_t const& other) const {
        return position_t{x() - other.x(), y() - other.y(), z() - other.z()};
    }
    [[nodiscard]] position_t operator*(double const scalar) const {
        return position_t{x() * scalar, y() * scalar, z() * scalar};
    }
    [[nodiscard]] bool operator==(position_t const& other) const {
        return coords_ == other.coords_;
    }
};

// quaternion_t class
class quaternion_t {
   private:
    std::array<double, 4> components_;

    void normalize() {
        double const norm = std::sqrt(std::pow(components_[0], 2) + std::pow(components_[1], 2) +
                                      std::pow(components_[2], 2) + std::pow(components_[3], 2));
        if (norm > 1e-10) {
            for (auto& c : components_) {
                c /= norm;
            }
        }
    }

   public:
    quaternion_t() : components_{0.0, 0.0, 0.0, 1.0} {}
    quaternion_t(double const x, double const y, double const z, double const w)
        : components_{x, y, z, w} {
        normalize();
    }

    [[nodiscard]] static quaternion_t from_euler(radian_t const roll, radian_t const pitch,
                                                 radian_t const yaw) {
        double const cy = std::cos(yaw.value * 0.5);
        double const sy = std::sin(yaw.value * 0.5);
        double const cp = std::cos(pitch.value * 0.5);
        double const sp = std::sin(pitch.value * 0.5);
        double const cr = std::cos(roll.value * 0.5);
        double const sr = std::sin(roll.value * 0.5);

        return quaternion_t{sr * cp * cy - cr * sp * sy, cr * sp * cy + sr * cp * sy,
                            cr * cp * sy - sr * sp * cy, cr * cp * cy + sr * sp * sy};
    }

    [[nodiscard]] double const& x() const { return components_[0]; }
    [[nodiscard]] double const& y() const { return components_[1]; }
    [[nodiscard]] double const& z() const { return components_[2]; }
    [[nodiscard]] double const& w() const { return components_[3]; }
    [[nodiscard]] double& x() { return components_[0]; }
    [[nodiscard]] double& y() { return components_[1]; }
    [[nodiscard]] double& z() { return components_[2]; }
    [[nodiscard]] double& w() { return components_[3]; }

    [[nodiscard]] quaternion_t operator*(quaternion_t const& other) const {
        return quaternion_t{w() * other.x() + x() * other.w() + y() * other.z() - z() * other.y(),
                            w() * other.y() - x() * other.z() + y() * other.w() + z() * other.x(),
                            w() * other.z() + x() * other.y() - y() * other.x() + z() * other.w(),
                            w() * other.w() - x() * other.x() - y() * other.y() - z() * other.z()};
    }

    [[nodiscard]] quaternion_t operator+(quaternion_t const& other) const {
        return quaternion_t{x() + other.x(), y() + other.y(), z() + other.z(), w() + other.w()};
    }
    [[nodiscard]] quaternion_t operator-(quaternion_t const& other) const {
        return quaternion_t{x() - other.x(), y() - other.y(), z() - other.z(), w() - other.w()};
    }
    [[nodiscard]] quaternion_t operator/(double const scalar) const {
        return quaternion_t{x() / scalar, y() / scalar, z() / scalar, w() / scalar};
    }

    [[nodiscard]] quaternion_t conjugate() const { return quaternion_t{-x(), -y(), -z(), w()}; }
    [[nodiscard]] bool operator==(quaternion_t const& other) const {
        return components_ == other.components_;
    }
};

[[nodiscard]] bool near(quaternion_t const& q1, quaternion_t const& q2,
                        double const tolerance = 0.001) {
    // Quaternions q and -q represent the same rotation
    bool const same_sign =
        std::abs(q1.x() - q2.x()) <= tolerance && std::abs(q1.y() - q2.y()) <= tolerance &&
        std::abs(q1.z() - q2.z()) <= tolerance && std::abs(q1.w() - q2.w()) <= tolerance;
    bool const opposite_sign =
        std::abs(q1.x() + q2.x()) <= tolerance && std::abs(q1.y() + q2.y()) <= tolerance &&
        std::abs(q1.z() + q2.z()) <= tolerance && std::abs(q1.w() + q2.w()) <= tolerance;
    return same_sign || opposite_sign;
}

// transformation_t class
class transformation_t {
   private:
    std::array<double, 16> matrix_;

    explicit transformation_t(std::array<double, 16> const& mat) : matrix_{mat} {}

    static std::array<double, 16> build_matrix(position_t const& pos, quaternion_t const& rot) {
        double const x = rot.x();
        double const y = rot.y();
        double const z = rot.z();
        double const w = rot.w();

        double const xx = x * x, yy = y * y, zz = z * z;
        double const xy = x * y, xz = x * z, yz = y * z;
        double const wx = w * x, wy = w * y, wz = w * z;

        return {1 - 2 * (yy + zz),
                2 * (xy - wz),
                2 * (xz + wy),
                pos.x().value,
                2 * (xy + wz),
                1 - 2 * (xx + zz),
                2 * (yz - wx),
                pos.y().value,
                2 * (xz - wy),
                2 * (yz + wx),
                1 - 2 * (xx + yy),
                pos.z().value,
                0.0,
                0.0,
                0.0,
                1.0};
    }

   public:
    transformation_t(position_t const& pos, quaternion_t const& rot)
        : matrix_{build_matrix(pos, rot)} {}

    transformation_t(transformation_t const& other) = default;
    transformation_t(transformation_t&& other) noexcept = default;
    transformation_t& operator=(transformation_t const& other) = default;
    transformation_t& operator=(transformation_t&& other) noexcept = default;
    ~transformation_t() = default;

    [[nodiscard]] position_t position() const {
        return position_t{meter_t{matrix_[3]}, meter_t{matrix_[7]}, meter_t{matrix_[11]}};
    }

    [[nodiscard]] quaternion_t rotation() const {
        double const m00 = matrix_[0], m01 = matrix_[1], m02 = matrix_[2];
        double const m10 = matrix_[4], m11 = matrix_[5], m12 = matrix_[6];
        double const m20 = matrix_[8], m21 = matrix_[9], m22 = matrix_[10];

        double const trace = m00 + m11 + m22;

        if (trace > 0) {
            double const s = 0.5 / std::sqrt(trace + 1.0);
            return quaternion_t{(m21 - m12) * s, (m02 - m20) * s, (m10 - m01) * s, 0.25 / s};
        }

        if (m00 > m11 && m00 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m00 - m11 - m22);
            return quaternion_t{0.25 * s, (m01 + m10) / s, (m02 + m20) / s, (m21 - m12) / s};
        }

        if (m11 > m22) {
            double const s = 2.0 * std::sqrt(1.0 + m11 - m00 - m22);
            return quaternion_t{(m01 + m10) / s, 0.25 * s, (m12 + m21) / s, (m02 - m20) / s};
        }

        double const s = 2.0 * std::sqrt(1.0 + m22 - m00 - m11);
        return quaternion_t{(m02 + m20) / s, (m12 + m21) / s, 0.25 * s, (m10 - m01) / s};
    }

    [[nodiscard]] transformation_t operator*(transformation_t const& other) const {
        std::array<double, 16> result{};
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += matrix_[row * 4 + k] * other.matrix_[k * 4 + col];
                }
                result[row * 4 + col] = sum;
            }
        }
        return transformation_t{result};
    }

    [[nodiscard]] position_t operator*(position_t const& point) const {
        double const x = point.x().value;
        double const y = point.y().value;
        double const z = point.z().value;

        return position_t{meter_t{matrix_[0] * x + matrix_[1] * y + matrix_[2] * z + matrix_[3]},
                          meter_t{matrix_[4] * x + matrix_[5] * y + matrix_[6] * z + matrix_[7]},
                          meter_t{matrix_[8] * x + matrix_[9] * y + matrix_[10] * z + matrix_[11]}};
    }

    [[nodiscard]] bool operator==(transformation_t const& other) const {
        return matrix_ == other.matrix_;
    }

    // Allow access to matrix for averaging
    [[nodiscard]] std::array<double, 16> const& matrix() const { return matrix_; }
};

[[nodiscard]] bool near(transformation_t const& tf1, transformation_t const& tf2,
                        meter_t const pos_tolerance = meter_t{0.001},
                        double const rot_tolerance = 0.001) {
    position_t const p1 = tf1.position();
    position_t const p2 = tf2.position();
    bool const pos_near = std::abs(p1.x().value - p2.x().value) <= pos_tolerance.value &&
                          std::abs(p1.y().value - p2.y().value) <= pos_tolerance.value &&
                          std::abs(p1.z().value - p2.z().value) <= pos_tolerance.value;
    return pos_near && near(tf1.rotation(), tf2.rotation(), rot_tolerance);
}

// ============================================================================
// Primary Template: Linear Average
// ============================================================================

template <typename T>
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

// ============================================================================
// Quaternion Specialization: Eigenvector Method
// ============================================================================

template <>
quaternion_t average<quaternion_t>(std::span<quaternion_t const> quaternions) {
    if (quaternions.empty()) {
        return quaternion_t{};
    }

    if (quaternions.size() == 1) {
        return quaternions[0];
    }

    // Build 4x4 accumulator matrix M = Σ qᵢ ⊗ qᵢᵀ
    std::array<std::array<double, 4>, 4> M{};

    for (auto const& q : quaternions) {
        std::array<double, 4> const qv = {q.x(), q.y(), q.z(), q.w()};

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                M[i][j] += qv[i] * qv[j];
            }
        }
    }

    // Power iteration to find eigenvector with largest eigenvalue
    std::array<double, 4> v = {0.0, 0.0, 0.0, 1.0};

    constexpr int max_iterations = 20;
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::array<double, 4> v_new{};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                v_new[i] += M[i][j] * v[j];
            }
        }

        double const norm = std::sqrt(v_new[0] * v_new[0] + v_new[1] * v_new[1] +
                                      v_new[2] * v_new[2] + v_new[3] * v_new[3]);

        if (norm < 1e-10) {
            break;
        }

        for (size_t i = 0; i < 4; ++i) {
            v[i] = v_new[i] / norm;
        }
    }

    return quaternion_t{v[0], v[1], v[2], v[3]};
}

// ============================================================================
// Transformation Specialization: Matrix Projection Method
// ============================================================================

using Mat3 = std::array<std::array<double, 3>, 3>;

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

Mat3 mat3_transpose(Mat3 const& A) {
    Mat3 result{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

Mat3 mat3_inverse(Mat3 const& A) {
    double const det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                       A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    if (std::abs(det) < 1e-10) {
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

Mat3 mat3_identity() {
    return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
}

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

    return Z;
}

Mat3 extract_rotation_matrix(transformation_t const& tf) {
    auto const& m = tf.matrix();
    return {{{m[0], m[1], m[2]}, {m[4], m[5], m[6]}, {m[8], m[9], m[10]}}};
}

quaternion_t mat3_to_quaternion(Mat3 const& R) {
    double const trace = R[0][0] + R[1][1] + R[2][2];

    if (trace > 0) {
        double const s = 0.5 / std::sqrt(trace + 1.0);
        return quaternion_t{(R[2][1] - R[1][2]) * s, (R[0][2] - R[2][0]) * s,
                            (R[1][0] - R[0][1]) * s, 0.25 / s};
    }

    if (R[0][0] > R[1][1] && R[0][0] > R[2][2]) {
        double const s = 2.0 * std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]);
        return quaternion_t{0.25 * s, (R[0][1] + R[1][0]) / s, (R[0][2] + R[2][0]) / s,
                            (R[2][1] - R[1][2]) / s};
    }

    if (R[1][1] > R[2][2]) {
        double const s = 2.0 * std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]);
        return quaternion_t{(R[0][1] + R[1][0]) / s, 0.25 * s, (R[1][2] + R[2][1]) / s,
                            (R[0][2] - R[2][0]) / s};
    }

    double const s = 2.0 * std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]);
    return quaternion_t{(R[0][2] + R[2][0]) / s, (R[1][2] + R[2][1]) / s, 0.25 * s,
                        (R[1][0] - R[0][1]) / s};
}

template <>
transformation_t average<transformation_t>(std::span<transformation_t const> transforms) {
    if (transforms.empty()) {
        return transformation_t{position_t{meter_t{0}, meter_t{0}, meter_t{0}}, quaternion_t{}};
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
    position_t const avg_pos{meter_t{sum_x / n}, meter_t{sum_y / n}, meter_t{sum_z / n}};

    // Step 2: Average rotation matrices element-wise
    Mat3 R_mean{};
    for (auto const& tf : transforms) {
        Mat3 const R = extract_rotation_matrix(tf);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                R_mean[i][j] += R[i][j] / n;
            }
        }
    }

    // Step 3: Project R_mean to SO(3) via polar decomposition
    Mat3 const Rt = mat3_transpose(R_mean);
    Mat3 const RtR = mat3_multiply(Rt, R_mean);
    Mat3 const RtR_inv_sqrt = mat3_sqrt_inverse(RtR);
    Mat3 const R_avg = mat3_multiply(R_mean, RtR_inv_sqrt);

    quaternion_t const avg_rot = mat3_to_quaternion(R_avg);

    return transformation_t{avg_pos, avg_rot};
}

// ============================================================================
// Helper: Extract yaw from quaternion
// ============================================================================

double extract_yaw(quaternion_t const& q) {
    // yaw (z-axis rotation)
    double const siny_cosp = 2 * (q.w() * q.z() + q.x() * q.y());
    double const cosy_cosp = 1 - 2 * (q.y() * q.y() + q.z() * q.z());
    return std::atan2(siny_cosp, cosy_cosp);
}

// ============================================================================
// Tests
// ============================================================================

// ----- Primary Template Tests -----

TEST(RotationAveragingTest, LinearAverageDoubles) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_DOUBLE_EQ(average<double>(data), 3.0);
}

TEST(RotationAveragingTest, LinearAverageIntegers) {
    std::array<int, 4> const data = {10, 20, 30, 40};
    EXPECT_EQ(average<int>(data), 25);
}

TEST(RotationAveragingTest, LinearAverageEmpty) {
    std::vector<double> const empty;
    EXPECT_DOUBLE_EQ(average<double>(empty), 0.0);
}

TEST(RotationAveragingTest, LinearAverageSingleElement) {
    std::array<double, 1> const data = {42.0};
    EXPECT_DOUBLE_EQ(average<double>(data), 42.0);
}

// ----- Quaternion Specialization Tests -----

TEST(RotationAveragingTest, QuaternionAverageIdentity) {
    std::array<quaternion_t, 3> const quats = {quaternion_t{}, quaternion_t{}, quaternion_t{}};
    quaternion_t const avg = average<quaternion_t>(quats);

    // Average of identity quaternions should be identity
    EXPECT_TRUE(near(avg, quaternion_t{}, 0.001));
}

TEST(RotationAveragingTest, QuaternionAverageSingleElement) {
    using namespace literals;
    quaternion_t const q = quaternion_t::from_euler(0.1_rad, 0.2_rad, 0.3_rad);
    std::array<quaternion_t, 1> const quats = {q};

    quaternion_t const avg = average<quaternion_t>(quats);
    EXPECT_TRUE(near(avg, q, 0.001));
}

TEST(RotationAveragingTest, QuaternionAverageSymmetric) {
    using namespace literals;
    // Two rotations symmetric around identity
    quaternion_t const q1 = quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.1_rad);
    quaternion_t const q2 = quaternion_t::from_euler(0.0_rad, 0.0_rad, radian_t{-0.1});

    std::array<quaternion_t, 2> const quats = {q1, q2};
    quaternion_t const avg = average<quaternion_t>(quats);

    // Average should be close to identity
    EXPECT_TRUE(near(avg, quaternion_t{}, 0.01));
}

// ----- Transformation Specialization Tests -----

TEST(RotationAveragingTest, TransformAveragePositionOnly) {
    using namespace literals;
    position_t const p1{1.0_m, 0.0_m, 0.0_m};
    position_t const p2{2.0_m, 0.0_m, 0.0_m};
    position_t const p3{3.0_m, 0.0_m, 0.0_m};

    std::array<transformation_t, 3> const transforms = {transformation_t{p1, quaternion_t{}},
                                                        transformation_t{p2, quaternion_t{}},
                                                        transformation_t{p3, quaternion_t{}}};

    transformation_t const avg = average<transformation_t>(transforms);

    // Position should be linear average
    EXPECT_NEAR(avg.position().x().value, 2.0, 0.001);
    EXPECT_NEAR(avg.position().y().value, 0.0, 0.001);
    EXPECT_NEAR(avg.position().z().value, 0.0, 0.001);

    // Rotation should be identity
    EXPECT_TRUE(near(avg.rotation(), quaternion_t{}, 0.001));
}

TEST(RotationAveragingTest, TransformAverageRotationOnly) {
    using namespace literals;
    position_t const origin{0.0_m, 0.0_m, 0.0_m};

    quaternion_t const q1 = quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.1_rad);
    quaternion_t const q2 = quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.2_rad);
    quaternion_t const q3 = quaternion_t::from_euler(0.0_rad, 0.0_rad, 0.3_rad);

    std::array<transformation_t, 3> const transforms = {
        transformation_t{origin, q1}, transformation_t{origin, q2}, transformation_t{origin, q3}};

    transformation_t const avg = average<transformation_t>(transforms);

    // Position should stay at origin
    EXPECT_NEAR(avg.position().x().value, 0.0, 0.001);
    EXPECT_NEAR(avg.position().y().value, 0.0, 0.001);
    EXPECT_NEAR(avg.position().z().value, 0.0, 0.001);

    // Rotation should be around 0.2 rad (average of 0.1, 0.2, 0.3)
    double const yaw = extract_yaw(avg.rotation());
    EXPECT_NEAR(yaw, 0.2, 0.01);
}

// ----- Validation Test 1: Pure Yaw Rotations -----

TEST(RotationAveragingTest, PureYawRotationsMatch) {
    using namespace literals;

    // Create rotations with only yaw (rotation around Z axis)
    std::array<double, 3> const yaw_angles = {0.1, 0.2, 0.3};  // radians

    // Method 1: Average angles linearly
    double const avg_angle = average<double>(yaw_angles);  // = 0.2

    // Method 2: Convert to quaternions, average, extract yaw
    std::array<quaternion_t, 3> const quats = {
        quaternion_t::from_euler(0.0_rad, 0.0_rad, radian_t{yaw_angles[0]}),
        quaternion_t::from_euler(0.0_rad, 0.0_rad, radian_t{yaw_angles[1]}),
        quaternion_t::from_euler(0.0_rad, 0.0_rad, radian_t{yaw_angles[2]})};
    quaternion_t const avg_quat = average<quaternion_t>(quats);
    double const quat_yaw = extract_yaw(avg_quat);

    // Method 3: Convert to transforms, average, extract yaw
    position_t const origin{0.0_m, 0.0_m, 0.0_m};
    std::array<transformation_t, 3> const transforms = {transformation_t{origin, quats[0]},
                                                        transformation_t{origin, quats[1]},
                                                        transformation_t{origin, quats[2]}};
    transformation_t const avg_tf = average<transformation_t>(transforms);
    double const tf_yaw = extract_yaw(avg_tf.rotation());

    // All three should match!
    EXPECT_NEAR(avg_angle, quat_yaw, 0.001);
    EXPECT_NEAR(avg_angle, tf_yaw, 0.001);
    EXPECT_NEAR(quat_yaw, tf_yaw, 0.001);
}

// ----- Validation Test 2: Combined Pitch and Roll -----

TEST(RotationAveragingTest, CombinedPitchRollMatch) {
    using namespace literals;

    // Create rotations with both pitch and roll (no yaw)
    std::array<quaternion_t, 3> const quats = {
        quaternion_t::from_euler(0.1_rad, 0.2_rad, 0.0_rad),
        quaternion_t::from_euler(0.2_rad, 0.1_rad, 0.0_rad),
        quaternion_t::from_euler(0.15_rad, 0.15_rad, 0.0_rad)};

    position_t const origin{0.0_m, 0.0_m, 0.0_m};
    std::array<transformation_t, 3> const transforms = {transformation_t{origin, quats[0]},
                                                        transformation_t{origin, quats[1]},
                                                        transformation_t{origin, quats[2]}};

    quaternion_t const avg_quat = average<quaternion_t>(quats);
    transformation_t const avg_tf = average<transformation_t>(transforms);

    // Extract rotation from transform and compare to quaternion
    quaternion_t const tf_rotation = avg_tf.rotation();

    // Should match (within numerical tolerance)
    EXPECT_TRUE(near(avg_quat, tf_rotation, 0.01));
}

// ----- Additional Edge Cases -----

TEST(RotationAveragingTest, LargeRotationSet) {
    using namespace literals;

    // Average many small rotations around Z
    std::vector<quaternion_t> quats;
    std::vector<transformation_t> transforms;
    position_t const origin{0.0_m, 0.0_m, 0.0_m};

    for (int i = 0; i < 100; ++i) {
        double const angle = 0.1 + 0.002 * i;  // 0.1 to 0.298 radians
        quaternion_t const q = quaternion_t::from_euler(0.0_rad, 0.0_rad, radian_t{angle});
        quats.push_back(q);
        transforms.push_back(transformation_t{origin, q});
    }

    quaternion_t const avg_quat = average<quaternion_t>(quats);
    transformation_t const avg_tf = average<transformation_t>(transforms);

    // Both methods should agree
    EXPECT_TRUE(near(avg_quat, avg_tf.rotation(), 0.01));

    // Average should be around 0.199 radians (midpoint)
    double const yaw = extract_yaw(avg_quat);
    EXPECT_NEAR(yaw, 0.199, 0.02);
}

TEST(RotationAveragingTest, MixedPositionAndRotation) {
    using namespace literals;

    std::array<transformation_t, 3> const transforms = {
        transformation_t{position_t{1.0_m, 2.0_m, 3.0_m},
                         quaternion_t::from_euler(0.1_rad, 0.0_rad, 0.0_rad)},
        transformation_t{position_t{4.0_m, 5.0_m, 6.0_m},
                         quaternion_t::from_euler(0.2_rad, 0.0_rad, 0.0_rad)},
        transformation_t{position_t{7.0_m, 8.0_m, 9.0_m},
                         quaternion_t::from_euler(0.3_rad, 0.0_rad, 0.0_rad)}};

    transformation_t const avg = average<transformation_t>(transforms);

    // Position should be linear average
    EXPECT_NEAR(avg.position().x().value, 4.0, 0.001);
    EXPECT_NEAR(avg.position().y().value, 5.0, 0.001);
    EXPECT_NEAR(avg.position().z().value, 6.0, 0.001);
}
