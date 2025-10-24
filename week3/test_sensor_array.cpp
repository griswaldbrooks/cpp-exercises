// Exercise 2: Sensor Array Statistics with std::span
// This exercise demonstrates using std::span for generic, zero-copy operations

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

// Calculate the average of sensor readings
double average(std::span<double const> const readings) {
    if (readings.empty()) {
        return 0.0;
    }
    double const sum = std::accumulate(readings.begin(), readings.end(), 0.0);
    return sum / static_cast<double>(readings.size());
}

// Find the minimum value in sensor readings
double min_value(std::span<double const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<double>::max();
    }
    return *std::min_element(readings.begin(), readings.end());
}

// Find the maximum value in sensor readings
double max_value(std::span<double const> const readings) {
    if (readings.empty()) {
        return std::numeric_limits<double>::lowest();
    }
    return *std::max_element(readings.begin(), readings.end());
}

// Count how many readings exceed a threshold
size_t count_above_threshold(std::span<double const> const readings, double const threshold) {
    return std::count_if(readings.begin(), readings.end(),
                         [threshold](double const value) { return value > threshold; });
}

// Normalize readings to [0, 1] range (modifies in-place)
void normalize(std::span<double> const readings) {
    if (readings.empty()) {
        return;
    }

    double const min = min_value(readings);
    double const max = max_value(readings);
    double const range = max - min;

    if (range == 0.0) {
        // All values are the same, set to 0.5
        std::fill(readings.begin(), readings.end(), 0.5);
        return;
    }

    for (auto& value : readings) {
        value = (value - min) / range;
    }
}

// ===== Tests =====

TEST(SensorArrayTest, AverageFromArray) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_DOUBLE_EQ(average(data), 3.0);
}

TEST(SensorArrayTest, AverageFromVector) {
    std::vector<double> const data = {2.0, 4.0, 6.0};
    EXPECT_DOUBLE_EQ(average(data), 4.0);
}

TEST(SensorArrayTest, AverageFromCArray) {
    double const data[] = {10.0, 20.0, 30.0};
    EXPECT_DOUBLE_EQ(average(std::span{data}), 20.0);
}

TEST(SensorArrayTest, AverageEmpty) {
    std::vector<double> const empty;
    EXPECT_DOUBLE_EQ(average(empty), 0.0);
}

TEST(SensorArrayTest, MinValueFromArray) {
    std::array<double, 5> const data = {5.0, 2.0, 8.0, 1.0, 4.0};
    EXPECT_DOUBLE_EQ(min_value(data), 1.0);
}

TEST(SensorArrayTest, MinValueNegative) {
    std::vector<double> const data = {-5.0, -2.0, -10.0, -1.0};
    EXPECT_DOUBLE_EQ(min_value(data), -10.0);
}

TEST(SensorArrayTest, MaxValueFromArray) {
    std::array<double, 5> const data = {5.0, 2.0, 8.0, 1.0, 4.0};
    EXPECT_DOUBLE_EQ(max_value(data), 8.0);
}

TEST(SensorArrayTest, MaxValueFromVector) {
    std::vector<double> const data = {1.5, 2.5, 3.5, 2.0};
    EXPECT_DOUBLE_EQ(max_value(data), 3.5);
}

TEST(SensorArrayTest, CountAboveThresholdSimple) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_EQ(count_above_threshold(data, 3.0), 2);  // 4.0 and 5.0
}

TEST(SensorArrayTest, CountAboveThresholdNone) {
    std::array<double, 3> const data = {1.0, 2.0, 3.0};
    EXPECT_EQ(count_above_threshold(data, 10.0), 0);
}

TEST(SensorArrayTest, CountAboveThresholdAll) {
    std::array<double, 3> const data = {5.0, 6.0, 7.0};
    EXPECT_EQ(count_above_threshold(data, 1.0), 3);
}

TEST(SensorArrayTest, CountAboveThresholdExactMatch) {
    std::array<double, 5> const data = {1.0, 2.0, 3.0, 4.0, 5.0};
    EXPECT_EQ(count_above_threshold(data, 3.0), 2);  // > not >=
}

TEST(SensorArrayTest, NormalizeBasic) {
    std::vector<double> data = {2.0, 4.0, 6.0, 8.0};
    normalize(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);  // min -> 0
    EXPECT_NEAR(data[1], 0.333, 0.01);
    EXPECT_NEAR(data[2], 0.666, 0.01);
    EXPECT_NEAR(data[3], 1.0, 0.001);  // max -> 1
}

TEST(SensorArrayTest, NormalizeAlreadyNormalized) {
    std::vector<double> data = {0.0, 0.5, 1.0};
    normalize(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 1.0, 0.001);
}

TEST(SensorArrayTest, NormalizeNegativeValues) {
    std::vector<double> data = {-10.0, 0.0, 10.0};
    normalize(data);

    EXPECT_NEAR(data[0], 0.0, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 1.0, 0.001);
}

TEST(SensorArrayTest, NormalizeSameValues) {
    std::vector<double> data = {5.0, 5.0, 5.0};
    normalize(data);

    // All same -> set to 0.5
    EXPECT_NEAR(data[0], 0.5, 0.001);
    EXPECT_NEAR(data[1], 0.5, 0.001);
    EXPECT_NEAR(data[2], 0.5, 0.001);
}

TEST(SensorArrayTest, NormalizeEmpty) {
    std::vector<double> data;
    normalize(data);  // Should not crash
    EXPECT_TRUE(data.empty());
}

TEST(SensorArrayTest, WorksWithArraySubspan) {
    std::array<double, 10> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // Create a subspan (like a slice)
    auto const first_half = std::span{data}.subspan(0, 5);
    auto const second_half = std::span{data}.subspan(5, 5);

    EXPECT_DOUBLE_EQ(average(first_half), 3.0);
    EXPECT_DOUBLE_EQ(average(second_half), 8.0);
}

TEST(SensorArrayTest, ConstCorrectnessCompiles) {
    // This test verifies const-correctness at compile time
    std::array<double, 3> const const_data = {1.0, 2.0, 3.0};

    // These should compile (const span from const container)
    [[maybe_unused]] auto const avg = average(const_data);
    [[maybe_unused]] auto const min = min_value(const_data);
    [[maybe_unused]] auto const max = max_value(const_data);
    [[maybe_unused]] auto const count = count_above_threshold(const_data, 1.5);

    // normalize(const_data);  // This should NOT compile (uncomment to verify)
}

TEST(SensorArrayTest, MutableSpanModifiesOriginal) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    std::span<double> const view{data};

    normalize(view);

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
    EXPECT_DOUBLE_EQ(average(arr), average(vec));
    EXPECT_DOUBLE_EQ(average(vec), average(std::span{c_arr}));
}
