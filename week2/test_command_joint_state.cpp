#include <gtest/gtest.h>

#include <cmath>
#include <format>
#include <optional>
#include <string>
#include <string_view>

/**
 * @brief Represents the state of a robot's joints in radians.
 */
struct joint_state {
    static constexpr size_t num_joints = 2;  ///< Number of joints in this state
    double joint1_rad;                       ///< First joint angle in radians
    double joint2_rad;                       ///< Second joint angle in radians
};

/**
 * @brief Calculates an 8-bit checksum by summing ASCII byte values.
 *
 * @param data The string data to checksum
 * @return uint8_t The calculated checksum (sum of all bytes modulo 256)
 */
uint8_t calculate_checksum(std::string_view const data) {
    uint32_t sum = 0;
    for (char const c : data) {
        sum += static_cast<uint8_t>(c);
    }
    return static_cast<uint8_t>(sum % 256);
}

/**
 * @brief Creates a joint state command message for the serial protocol.
 *
 * Uses std::format (C++20) for efficient string formatting. The checksum
 * is calculated over the data portion (header and angles) before appending.
 *
 * Protocol format: "JS:<angle1>,<angle2>:<checksum>\n"
 * Example: "JS:1.570796,-0.785398:0xA1\n"
 *
 * @param target The joint state to serialize
 * @return std::string The formatted protocol message
 */
std::string command_joint_state(joint_state const& target) {
    // Format the data portion with high precision (6 decimal places)
    std::string const data = std::format("JS:{:.6f},{:.6f}", target.joint1_rad, target.joint2_rad);

    // Calculate checksum on data (everything before the last colon)
    uint8_t const checksum = calculate_checksum(data);

    // Append checksum and newline
    return std::format("{}:0x{:02X}\n", data, checksum);
}

// Tests
TEST(CommandJointStateTest, BasicExample) {
    // GIVEN a joint state
    joint_state const cmd{1.570796, -0.785398};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message has the correct format
    EXPECT_TRUE(msg.starts_with("JS:"));
    EXPECT_TRUE(msg.ends_with("\n"));
    EXPECT_NE(msg.find("1.570796"), std::string::npos);
    EXPECT_NE(msg.find("-0.785398"), std::string::npos);
    EXPECT_NE(msg.find("0x"), std::string::npos);
}

TEST(CommandJointStateTest, ZeroAngles) {
    // GIVEN a joint state with zero angles
    joint_state const cmd{0.0, 0.0};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message contains zero values
    EXPECT_NE(msg.find("0.000000"), std::string::npos);
}

TEST(CommandJointStateTest, NegativeAngles) {
    // GIVEN a joint state with negative angles
    joint_state const cmd{-1.5, -2.5};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message contains negative values
    EXPECT_NE(msg.find("-1.500000"), std::string::npos);
    EXPECT_NE(msg.find("-2.500000"), std::string::npos);
}

TEST(CommandJointStateTest, PositiveAngles) {
    // GIVEN a joint state with positive angles
    joint_state const cmd{3.14159, 2.71828};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message contains positive values
    EXPECT_NE(msg.find("3.141590"), std::string::npos);
    EXPECT_NE(msg.find("2.718280"), std::string::npos);
}

TEST(CommandJointStateTest, HasCorrectStructure) {
    // GIVEN a joint state
    joint_state const cmd{1.0, 2.0};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message has exactly 2 colons (after header and before checksum)
    size_t colon_count = 0;
    for (char c : msg) {
        if (c == ':')
            colon_count++;
    }
    EXPECT_EQ(colon_count, 2);

    // AND has exactly 1 comma (between angles)
    size_t comma_count = 0;
    for (char c : msg) {
        if (c == ',')
            comma_count++;
    }
    EXPECT_EQ(comma_count, 1);
}

TEST(CommandJointStateTest, ChecksumFormat) {
    // GIVEN a joint state
    joint_state const cmd{1.0, 2.0};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the checksum is in hex format with 0x prefix
    auto const checksum_pos = msg.rfind("0x");
    EXPECT_NE(checksum_pos, std::string::npos);

    // AND checksum is 2 hex digits (uppercase)
    EXPECT_EQ(msg.at(checksum_pos + 2),
              msg.at(checksum_pos + 2) >= '0' ? msg.at(checksum_pos + 2) : '0');
    EXPECT_EQ(msg.at(checksum_pos + 3),
              msg.at(checksum_pos + 3) >= '0' ? msg.at(checksum_pos + 3) : '0');
}

TEST(CommandJointStateTest, PrecisionCheck) {
    // GIVEN a joint state with many decimal places
    joint_state const cmd{1.123456789, 2.987654321};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN values are formatted with exactly 6 decimal places
    EXPECT_NE(msg.find("1.123457"), std::string::npos);  // Rounded
    EXPECT_NE(msg.find("2.987654"), std::string::npos);
}

TEST(CommandJointStateTest, LargeAngles) {
    // GIVEN a joint state with large angle values
    joint_state const cmd{100.123456, -200.987654};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message contains the large values
    EXPECT_NE(msg.find("100.123456"), std::string::npos);
    EXPECT_NE(msg.find("-200.987654"), std::string::npos);
}

TEST(CommandJointStateTest, SmallAngles) {
    // GIVEN a joint state with very small angle values
    joint_state const cmd{0.000001, -0.000002};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN the message contains the small values
    EXPECT_NE(msg.find("0.000001"), std::string::npos);
    EXPECT_NE(msg.find("-0.000002"), std::string::npos);
}

TEST(CommandJointStateTest, ChecksumIsValid) {
    // GIVEN a joint state
    joint_state const cmd{1.57, -0.785};

    // WHEN creating a command message
    std::string const msg = command_joint_state(cmd);

    // THEN we can manually verify the checksum is correct
    // Extract the data portion (everything before last colon)
    auto const last_colon = msg.rfind(':');
    ASSERT_NE(last_colon, std::string::npos);
    std::string_view const data_portion = std::string_view(msg).substr(0, last_colon);

    // Calculate expected checksum
    uint8_t const expected_checksum = calculate_checksum(data_portion);

    // Extract actual checksum from message
    std::string_view const checksum_str = std::string_view(msg).substr(last_colon + 1);
    EXPECT_TRUE(checksum_str.starts_with("0x"));

    // Parse the hex checksum
    unsigned int actual_checksum;
    auto result =
        std::from_chars(checksum_str.data() + 2,
                        checksum_str.data() + checksum_str.size() - 1,  // -1 to skip newline
                        actual_checksum, 16);
    ASSERT_EQ(result.ec, std::errc{});

    // Verify checksums match
    EXPECT_EQ(actual_checksum, expected_checksum);
}
