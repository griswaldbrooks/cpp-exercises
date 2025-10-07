cmake_minimum_required(VERSION 3.25)
project(week2_exercises CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find GTest
find_package(GTest REQUIRED)

# Enable testing
enable_testing()

# Exercise 1: Normalize Whitespace
add_executable(normalize_whitespace_test
    test_normalize_whitespace.cpp
)

target_link_libraries(normalize_whitespace_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME NormalizeWhitespaceTests COMMAND normalize_whitespace_test)

# Exercise 2: Split
add_executable(split_test
    test_split.cpp
)

target_link_libraries(split_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME SplitTests COMMAND split_test)

# Exercise 3: Serial Protocol Parser
add_executable(serial_protocol_test
    test_serial_protocol.cpp
)

target_link_libraries(serial_protocol_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME SerialProtocolTests COMMAND serial_protocol_test)

# Exercise 4: Command Joint State
add_executable(command_joint_state_test
    test_command_joint_state.cpp
)

target_link_libraries(command_joint_state_test
    PRIVATE
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME CommandJointStateTests COMMAND command_joint_state_test)
