// Exercise 2: Sensor Array Statistics with std::span
// This exercise demonstrates using std::span for generic, zero-copy operations

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

// Simplified quaternion_t for testing (supports arithmetic needed for averaging)
// In production, this would be shared from Exercise 1
class quaternion_t {
   private:
    std::array<double, 4> components_;

   public:
    quaternion_t() : components_{0.0, 0.0, 0.0, 1.0} {}
    quaternion_t(double const x, double const y, double const z, double const w)
        : components_{x, y, z, w} {}

    [[nodiscard]] double const& x() const { return components_[0]; }
    [[nodiscard]] double const& y() const { return components_[1]; }
    [[nodiscard]] double const& z() const { return components_[2]; }
    [[nodiscard]] double const& w() const { return components_[3]; }

    // Arithmetic operators for averaging
    [[nodiscard]] quaternion_t operator+(quaternion_t const& other) const {
        return quaternion_t{x() + other.x(), y() + other.y(), z() + other.z(), w() + other.w()};
    }

    [[nodiscard]] quaternion_t operator-(quaternion_t const& other) const {
        return quaternion_t{x() - other.x(), y() - other.y(), z() - other.z(), w() - other.w()};
    }

    [[nodiscard]] quaternion_t operator/(double const scalar) const {
        return quaternion_t{x() / scalar, y() / scalar, z() / scalar, w() / scalar};
    }
};

// Calculate the average of sensor readings - templated for any numeric type
template <typename T>
T average(std::span<T const> const readings) {
    if (readings.empty()) {
        return T{};
    }
    // Start with first element to avoid issues with T{} not being "zero" for all types
    // (e.g., quaternion_t{} is identity rotation {0,0,0,1}, not zero)
    T sum = readings[0];
    for (size_t i = 1; i < readings.size(); ++i) {
        sum = sum + readings[i];
    }
    return sum / static_cast<double>(readings.size());
}

// Find the minimum value in sensor readings - templated for any numeric type
template <typename T>
T min_value(std::span<T const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<T>::max();
    }
    return *std::min_element(readings.begin(), readings.end());
}

// Find the maximum value in sensor readings - templated for any numeric type
template <typename T>
T max_value(std::span<T const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<T>::lowest();
    }
    return *std::max_element(readings.begin(), readings.end());
}

// Count how many readings exceed a threshold - templated for any numeric type
template <typename T>
size_t count_above_threshold(std::span<T const> const readings, T const threshold) {
    return std::count_if(readings.begin(), readings.end(),
                         [threshold](T const value) { return value > threshold; });
}

// Normalize readings to [0, 1] range (modifies in-place) - templated for any numeric type
template <typename T>
void normalize(std::span<T> const readings) {
    if (readings.empty()) {
        return;
    }

    // Create const span for min/max computation
    std::span<T const> const const_readings{readings};
    T const min = min_value<T>(const_readings);
    T const max = max_value<T>(const_readings);
    T const range = max - min;

    if (range == T{}) {
        // All values are the same, set to 0.5
        std::fill(readings.begin(), readings.end(), static_cast<T>(0.5));
        return;
    }

    for (auto& value : readings) {
        value = (value - min) / range;
    }
}

// Sliding window average - templated for any numeric type
// Uses O(n) algorithm: compute first window sum, then slide by
// subtracting left element and adding right element
template <typename T>
std::vector<T> sliding_window_average(std::span<T const> const values, size_t const window_size) {
    if (values.size() < window_size || window_size == 0) {
        return {};
    }

    std::vector<T> result;
    result.reserve(values.size() - window_size + 1);

    // Compute sum of first window
    // Start with first value to avoid issues with T{} not being "zero" for all types
    T window_sum = values[0];
    for (size_t i = 1; i < window_size; ++i) {
        window_sum = window_sum + values[i];
    }
    result.push_back(window_sum / static_cast<double>(window_size));

    // Slide the window: remove leftmost, add new rightmost
    for (size_t i = window_size; i < values.size(); ++i) {
        window_sum = window_sum - values[i - window_size] + values[i];
        result.push_back(window_sum / static_cast<double>(window_size));
    }

    return result;
}

// ===== Tests =====

TEST(SensorArrayTest, AverageFromArray) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_DOUBLE_EQ(average<double>(data), 3.0);
}

TEST(SensorArrayTest, AverageFromVector) {
    std::vector<double> const data = {2.0, 4.0, 6.0};
    EXPECT_DOUBLE_EQ(average<double>(data), 4.0);
}

TEST(SensorArrayTest, AverageFromCArray) {
    double const data[] = {10.0, 20.0, 30.0};
    EXPECT_DOUBLE_EQ(average<double>(std::span{data}), 20.0);
}

TEST(SensorArrayTest, AverageEmpty) {
    std::vector<double> const empty;
    EXPECT_DOUBLE_EQ(average<double>(empty), 0.0);
}

TEST(SensorArrayTest, AverageWithIntegers) {
    std::array<int, 5> const data = {1, 2, 3, 4, 5};
    EXPECT_EQ(average<int>(data), 3);  // Integer division: 15/5 = 3
}

TEST(SensorArrayTest, AverageWithFloats) {
    std::array<float, 4> const data = {2.0f, 4.0f, 6.0f, 8.0f};
    EXPECT_FLOAT_EQ(average<float>(data), 5.0f);
}

TEST(SensorArrayTest, MinValueFromArray) {
    std::array<double, 5> const data = {5.0, 2.0, 8.0, 1.0, 4.0};
    EXPECT_DOUBLE_EQ(min_value<double>(data), 1.0);
}

TEST(SensorArrayTest, MinValueNegative) {
    std::vector<double> const data = {-5.0, -2.0, -10.0, -1.0};
    EXPECT_DOUBLE_EQ(min_value<double>(data), -10.0);
}

TEST(SensorArrayTest, MinValueWithIntegers) {
    std::array<int, 5> const data = {5, 2, 8, 1, 4};
    EXPECT_EQ(min_value<int>(data), 1);
}

TEST(SensorArrayTest, MaxValueFromArray) {
    std::array<double, 5> const data = {5.0, 2.0, 8.0, 1.0, 4.0};
    EXPECT_DOUBLE_EQ(max_value<double>(data), 8.0);
}

TEST(SensorArrayTest, MaxValueFromVector) {
    std::vector<double> const data = {1.5, 2.5, 3.5, 2.0};
    EXPECT_DOUBLE_EQ(max_value<double>(data), 3.5);
}

TEST(SensorArrayTest, MaxValueWithIntegers) {
    std::array<int, 4> const data = {10, 25, 15, 20};
    EXPECT_EQ(max_value<int>(data), 25);
}

TEST(SensorArrayTest, CountAboveThresholdSimple) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_EQ(count_above_threshold<double>(data, 3.0), 2);  // 4.0 and 5.0
}

TEST(SensorArrayTest, CountAboveThresholdNone) {
    std::array<double, 3> const data = {1.0, 2.0, 3.0};
    EXPECT_EQ(count_above_threshold<double>(data, 10.0), 0);
}

TEST(SensorArrayTest, CountAboveThresholdAll) {
    std::array<double, 3> const data = {5.0, 6.0, 7.0};
    EXPECT_EQ(count_above_threshold<double>(data, 1.0), 3);
}

TEST(SensorArrayTest, CountAboveThresholdExactMatch) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_EQ(count_above_threshold<double>(data, 3.0), 2);  // > not >=
}

TEST(SensorArrayTest, CountAboveThresholdWithIntegers) {
    std::array<int, 6> const data = {10, 20, 30, 40, 50, 60};
    EXPECT_EQ(count_above_threshold<int>(data, 35), 3);  // 40, 50, 60
}

TEST(SensorArrayTest, CountAboveThresholdWithFloats) {
    std::vector<float> const data = {1.5f, 2.5f, 3.5f, 4.5f};
    EXPECT_EQ(count_above_threshold<float>(data, 2.5f), 2);  // 3.5f, 4.5f
}

TEST(SensorArrayTest, NormalizeBasic) {
    std::vector<double> data = {2.0, 4.0, 6.0, 8.0};
    normalize<double>(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);  // min -> 0
    EXPECT_NEAR(data[1], 0.333, 0.01);
    EXPECT_NEAR(data[2], 0.666, 0.01);
    EXPECT_NEAR(data[3], 1.0, 0.001);  // max -> 1
}

TEST(SensorArrayTest, NormalizeAlreadyNormalized) {
    std::vector<double> data = {0.0, 0.5, 1.0};
    normalize<double>(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 1.0, 0.001);
}

TEST(SensorArrayTest, NormalizeNegativeValues) {
    std::vector<double> data = {-10.0, 0.0, 10.0};
    normalize<double>(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 1.0, 0.001);
}

TEST(SensorArrayTest, NormalizeSameValues) {
    std::vector<double> data = {5.0, 5.0, 5.0};
    normalize<double>(data);

    // All same -> set to 0.5
    EXPECT_NEAR(data[0], 0.5, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 0.5, 0.001);
}

TEST(SensorArrayTest, NormalizeEmpty) {
    std::vector<double> data;
    normalize<double>(data);  // Should not crash
    EXPECT_TRUE(data.empty());
}

TEST(SensorArrayTest, NormalizeWithFloats) {
    std::vector<float> data = {10.0f, 20.0f, 30.0f, 40.0f};
    normalize<float>(data);

    EXPECT_FLOAT_EQ(data[0], 0.0f);       // min -> 0
    EXPECT_NEAR(data[1], 0.333f, 0.01f);  // (20-10)/(40-10) ≈ 0.333
    EXPECT_NEAR(data[2], 0.666f, 0.01f);  // (30-10)/(40-10) ≈ 0.666
    EXPECT_FLOAT_EQ(data[3], 1.0f);       // max -> 1
}

TEST(SensorArrayTest, NormalizeWithIntegers) {
    std::vector<int> data = {0, 50, 100};
    normalize<int>(data);

    EXPECT_EQ(data[0], 0);  // min -> 0
    EXPECT_EQ(data[1], 0);  // (50-0)/(100-0) = 0 (integer division)
    EXPECT_EQ(data[2], 1);  // max -> 1
}

TEST(SensorArrayTest, WorksWithArraySubspan) {
    std::array<double, 10> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // Create a subspan (like a slice)
    auto const first_half = std::span{data}.subspan(0, 5);
    auto const second_half = std::span{data}.subspan(5, 5);

    EXPECT_DOUBLE_EQ(average<double>(first_half), 3.0);
    EXPECT_DOUBLE_EQ(average<double>(second_half), 8.0);
}

TEST(SensorArrayTest, ConstCorrectnessCompiles) {
    // This test verifies const-correctness at compile time
    std::array<double, 3> const const_data = {1.0, 2.0, 3.0};

    // These should compile (const span from const container)
    [[maybe_unused]] auto const avg = average<double>(const_data);
    [[maybe_unused]] auto const min = min_value<double>(const_data);
    [[maybe_unused]] auto const max = max_value<double>(const_data);
    [[maybe_unused]] auto const count = count_above_threshold<double>(const_data, 1.5);

    // normalize<double>(const_data);  // This should NOT compile (uncomment to verify)
}

TEST(SensorArrayTest, MutableSpanModifiesOriginal) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    std::span<double> const view{data};

    normalize<double>(view);

    // Original vector should be modified
    EXPECT_NEAR(data[0], 0.0, 0.001);
    EXPECT_NEAR(data[3], 1.0, 0.001);
}

TEST(SensorArrayTest, MultipleContainerTypes) {
    // Demonstrate that the same function works with different container types
    std::array<double, 4> arr = {1, 2, 3, 4};
    std::vector<double> vec = {1, 2, 3, 4};
    double c_arr[] = {1, 2, 3, 4};

    // All produce the same result
    EXPECT_DOUBLE_EQ(average<double>(arr), average<double>(vec));
    EXPECT_DOUBLE_EQ(average<double>(vec), average<double>(std::span{c_arr}));
}

// ===== Sliding Window Tests =====

TEST(SensorArrayTest, SlidingWindowAverageBasic) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto const windowed = sliding_window_average(std::span<double const>{data}, 3);

    ASSERT_EQ(windowed.size(), 3);         // 5 - 3 + 1 = 3
    EXPECT_NEAR(windowed[0], 2.0, 0.001);  // (1+2+3)/3
    EXPECT_NEAR(windowed[1], 3.0, 0.001);  // (2+3+4)/3
    EXPECT_NEAR(windowed[2], 4.0, 0.001);  // (3+4+5)/3
}

TEST(SensorArrayTest, SlidingWindowAverageWithIntegers) {
    std::vector<int> const int_data = {10, 20, 30, 40};
    auto const int_windowed = sliding_window_average(std::span<int const>{int_data}, 2);

    ASSERT_EQ(int_windowed.size(), 3);  // 4 - 2 + 1 = 3
    EXPECT_EQ(int_windowed[0], 15);     // (10+20)/2
    EXPECT_EQ(int_windowed[1], 25);     // (20+30)/2
    EXPECT_EQ(int_windowed[2], 35);     // (30+40)/2
}

TEST(SensorArrayTest, SlidingWindowAverageWindowSizeOne) {
    std::array<double, 4> const data = {1.0, 2.0, 3.0, 4.0};
    auto const windowed = sliding_window_average(std::span<double const>{data}, 1);

    ASSERT_EQ(windowed.size(), 4);
    EXPECT_DOUBLE_EQ(windowed[0], 1.0);
    EXPECT_DOUBLE_EQ(windowed[1], 2.0);
    EXPECT_DOUBLE_EQ(windowed[2], 3.0);
    EXPECT_DOUBLE_EQ(windowed[3], 4.0);
}

TEST(SensorArrayTest, SlidingWindowAverageFullWindow) {
    std::array<double, 5> const data = {2.0, 4.0, 6.0, 8.0, 10.0};
    auto const windowed = sliding_window_average(std::span<double const>{data}, 5);

    ASSERT_EQ(windowed.size(), 1);         // 5 - 5 + 1 = 1
    EXPECT_NEAR(windowed[0], 6.0, 0.001);  // (2+4+6+8+10)/5
}

TEST(SensorArrayTest, SlidingWindowAverageTooLarge) {
    std::array<double, 3> const data = {1.0, 2.0, 3.0};
    auto const windowed = sliding_window_average(std::span<double const>{data}, 5);

    EXPECT_TRUE(windowed.empty());
}

TEST(SensorArrayTest, SlidingWindowAverageEmpty) {
    std::vector<double> const empty;
    auto const windowed = sliding_window_average(std::span<double const>{empty}, 3);

    EXPECT_TRUE(windowed.empty());
}

TEST(SensorArrayTest, SlidingWindowAverageZeroWindow) {
    std::array<double, 3> const data = {1.0, 2.0, 3.0};
    auto const windowed = sliding_window_average(std::span<double const>{data}, 0);

    EXPECT_TRUE(windowed.empty());
}

// ===== Quaternion Tests =====

TEST(SensorArrayTest, AverageWithQuaternions) {
    // Test average with quaternion_t type
    // This demonstrates generic algorithms working with custom types that support +, -, /
    // Useful for averaging IMU orientation readings

    // Create 3 quaternion readings (component-wise for testing)
    std::array<quaternion_t, 3> const quats = {quaternion_t{0.1, 0.2, 0.3, 0.9},
                                               quaternion_t{0.2, 0.3, 0.4, 0.8},
                                               quaternion_t{0.3, 0.4, 0.5, 0.7}};

    // Compute average
    auto const avg = average<quaternion_t>(quats);

    // Expected: (0.1+0.2+0.3)/3, (0.2+0.3+0.4)/3, (0.3+0.4+0.5)/3, (0.9+0.8+0.7)/3
    // = (0.6/3, 0.9/3, 1.2/3, 2.4/3) = (0.2, 0.3, 0.4, 0.8)
    EXPECT_NEAR(avg.x(), 0.2, 0.001);
    EXPECT_NEAR(avg.y(), 0.3, 0.001);
    EXPECT_NEAR(avg.z(), 0.4, 0.001);
    EXPECT_NEAR(avg.w(), 0.8, 0.001);
}

TEST(SensorArrayTest, SlidingWindowAverageWithQuaternions) {
    // Test sliding window average with quaternion_t type
    // This demonstrates generic algorithms working with custom types that support +, -, /

    // Create 5 quaternion readings (component-wise for testing)
    std::array<quaternion_t, 5> const quats = {
        quaternion_t{0.1, 0.2, 0.3, 0.9}, quaternion_t{0.2, 0.3, 0.4, 0.8},
        quaternion_t{0.3, 0.4, 0.5, 0.7}, quaternion_t{0.4, 0.5, 0.6, 0.6},
        quaternion_t{0.5, 0.6, 0.7, 0.5}};

    // Compute 3-point sliding window average
    auto const windowed = sliding_window_average(std::span<quaternion_t const>{quats}, 3);

    ASSERT_EQ(windowed.size(), 3);  // 5 - 3 + 1 = 3

    // First window: average of quats[0], quats[1], quats[2]
    // Expected: (0.1+0.2+0.3)/3, (0.2+0.3+0.4)/3, (0.3+0.4+0.5)/3, (0.9+0.8+0.7)/3
    // = (0.6/3, 0.9/3, 1.2/3, 2.4/3) = (0.2, 0.3, 0.4, 0.8)
    EXPECT_NEAR(windowed[0].x(), 0.2, 0.001);
    EXPECT_NEAR(windowed[0].y(), 0.3, 0.001);
    EXPECT_NEAR(windowed[0].z(), 0.4, 0.001);
    EXPECT_NEAR(windowed[0].w(), 0.8, 0.001);

    // Second window: average of quats[1], quats[2], quats[3]
    // Expected: (0.2+0.3+0.4)/3, (0.3+0.4+0.5)/3, (0.4+0.5+0.6)/3, (0.8+0.7+0.6)/3
    // = (0.9/3, 1.2/3, 1.5/3, 2.1/3) = (0.3, 0.4, 0.5, 0.7)
    EXPECT_NEAR(windowed[1].x(), 0.3, 0.001);
    EXPECT_NEAR(windowed[1].y(), 0.4, 0.001);
    EXPECT_NEAR(windowed[1].z(), 0.5, 0.001);
    EXPECT_NEAR(windowed[1].w(), 0.7, 0.001);
}
