cmake_minimum_required(VERSION 3.25)
project(week3_exercises CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find GTest
find_package(GTest REQUIRED)

# Enable testing
enable_testing()

# Exercise 1: Transform3D
add_executable(transform3d_test
    test_transform3d.cpp
)

target_link_libraries(transform3d_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME Transform3DTests COMMAND transform3d_test)

# Exercise 2: Sensor Array
add_executable(sensor_array_test
    test_sensor_array.cpp
)

target_link_libraries(sensor_array_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME SensorArrayTests COMMAND sensor_array_test)
