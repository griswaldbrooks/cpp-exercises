#include <gtest/gtest.h>

#include <charconv>
#include <cmath>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Represents the state of a robot's joints in radians.
 */
struct joint_state {
    static constexpr size_t num_joints = 2;  ///< Number of joints in this state
    double joint1_rad;                       ///< First joint angle in radians
    double joint2_rad;                       ///< Second joint angle in radians
};

/**
 * @brief Configuration parameters for the serial protocol parser.
 */
struct protocol_config {
    char section_delimiter = ':';  ///< Delimiter between protocol sections (header, data, checksum)
    char angle_delimiter = ',';    ///< Delimiter between individual angle values
    std::string_view header = "JS";  ///< Protocol header identifier
    size_t expected_parts =
        3;  ///< Expected number of sections after splitting by section_delimiter
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
 * @brief Removes leading and trailing whitespace from a string view.
 *
 * @param str The input string view
 * @return std::string_view The trimmed string view (or empty if all whitespace)
 */
constexpr std::string_view trim(std::string_view const str) {
    constexpr std::string_view whitespace = " \t\n\r\f\v";
    auto const start = str.find_first_not_of(whitespace);
    if (start == std::string_view::npos)
        return "";  // All whitespace

    auto trimmed = str;
    trimmed.remove_prefix(start);
    auto const end = trimmed.find_last_not_of(whitespace);
    trimmed.remove_suffix(trimmed.size() - end - 1);
    return trimmed;
}

/**
 * @brief Splits a string view by a delimiter into a vector of string views.
 *
 * Uses C++20 ranges for lazy evaluation. The resulting string views reference
 * the original string data (zero-copy).
 *
 * @param str The input string view to split
 * @param delimiter The character delimiter to split on
 * @return std::vector<std::string_view> Vector of string view tokens
 */
std::vector<std::string_view> split_to_vector(std::string_view const str, char const delimiter) {
    // Create lazy view pipeline: split by delimiter, then convert subranges to string_views
    auto split_view = str | std::views::split(delimiter) | std::views::transform([](auto&& rng) {
                          return std::string_view(&*rng.begin(), std::ranges::distance(rng));
                      });

    // Materialize the lazy view into a vector
    return std::vector<std::string_view>(std::ranges::begin(split_view),
                                         std::ranges::end(split_view));
}

/**
 * @brief Parses a single floating-point angle from a string view.
 *
 * Uses std::from_chars for fast, exception-free parsing. Requires that the entire
 * string is consumed during parsing (no trailing garbage allowed).
 *
 * @param angle_str String view containing the angle value
 * @return std::optional<double> The parsed angle, or std::nullopt on parse failure
 */
std::optional<double> parse_angle(std::string_view const angle_str) {
    double angle;
    auto const [ptr, ec] =
        std::from_chars(angle_str.data(), angle_str.data() + angle_str.size(), angle);

    if (ec != std::errc{}) {
        return std::nullopt;
    }

    // Verify that the entire string was consumed (no trailing garbage)
    if (ptr != angle_str.data() + angle_str.size()) {
        return std::nullopt;
    }

    return angle;
}

/**
 * @brief Validates the checksum of a protocol message.
 *
 * The checksum is calculated over all characters from the beginning of the message
 * up to (but not including) the last section delimiter. The provided checksum must
 * be in hexadecimal format with "0x" prefix (e.g., "0xFD").
 *
 * @param message The full protocol message
 * @param checksum_str The checksum string from the message (hex format with 0x prefix)
 * @param section_delimiter The delimiter that separates the checksum from the data
 * @return bool True if checksum is valid, false otherwise
 */
bool validate_checksum(std::string_view const message, std::string_view const checksum_str,
                       char const section_delimiter) {
    // Verify checksum format (must start with "0x")
    if (not checksum_str.starts_with("0x")) {
        return false;
    }

    // Parse the provided checksum from hex string (skip "0x" prefix)
    unsigned int provided;
    auto const [ptr, ec] = std::from_chars(checksum_str.data() + 2,
                                           checksum_str.data() + checksum_str.size(), provided,
                                           16  // Base 16 for hexadecimal
    );

    if (ec != std::errc{}) {
        return false;
    }

    // Calculate expected checksum from message data (everything before last delimiter)
    auto const checksum_start = message.rfind(section_delimiter);
    uint8_t const calculated = calculate_checksum(message.substr(0, checksum_start));

    return provided == calculated;
}

/**
 * @brief Parses joint angles from a delimited string.
 *
 * Splits the input string by the angle delimiter and parses each angle value.
 * Expects exactly joint_state::num_joints angle values.
 * Uses .at() for bounds-checked access to prevent out-of-range errors.
 *
 * @param angles_str String containing the angle values separated by delimiters
 * @param angle_delimiter Character delimiter separating angle values
 * @return std::optional<joint_state> Parsed joint state, or std::nullopt on failure
 */
std::optional<joint_state> parse_joint_angles(std::string_view const angles_str,
                                              char const angle_delimiter) {
    // Split angles by angle delimiter
    auto const angles = split_to_vector(angles_str, angle_delimiter);
    if (angles.size() != joint_state::num_joints) {
        return std::nullopt;
    }

    // Parse each angle individually using .at() for bounds-checked access
    auto const joint1_opt = parse_angle(angles.at(0));
    if (not joint1_opt.has_value()) {
        return std::nullopt;
    }

    auto const joint2_opt = parse_angle(angles.at(1));
    if (not joint2_opt.has_value()) {
        return std::nullopt;
    }

    return joint_state{joint1_opt.value(), joint2_opt.value()};
}

/**
 * @brief Parses a joint state message from the serial protocol.
 *
 * Protocol format: "JS:<angle1>,<angle2>:<checksum>\n"
 * Example: "JS:1.57,-0.785:0xFD\n"
 *
 * Validation order:
 * 1. Trim whitespace and validate header
 * 2. Split into sections by section delimiter
 * 3. Validate checksum (early exit before expensive parsing)
 * 4. Parse joint angle values
 *
 * @param message The protocol message to parse
 * @param config Protocol configuration parameters (default: standard config)
 * @return std::optional<joint_state> Parsed joint state, or std::nullopt on any validation failure
 */
std::optional<joint_state> read_joint_state(std::string_view const message,
                                            protocol_config const config = {}) {
    // Trim whitespace and validate header
    auto const msg = trim(message);
    if (not msg.starts_with(config.header) or msg.empty()) {
        return std::nullopt;
    }

    // Split by section delimiter (expected: ["JS", "angles", "checksum"])
    auto const parts_vec = split_to_vector(msg, config.section_delimiter);
    if (parts_vec.size() != config.expected_parts) {
        return std::nullopt;
    }

    // Validate checksum early before expensive angle parsing
    if (not validate_checksum(msg, parts_vec[2], config.section_delimiter)) {
        return std::nullopt;
    }

    // Parse and return joint angles
    return parse_joint_angles(parts_vec[1], config.angle_delimiter);
}

// Tests for discovered bugs (these should fail until bugs are fixed)
TEST(BugTests, ChecksumWithOnlyPrefix) {
    // BUG: validate_checksum doesn't check if there are hex digits after "0x"
    // GIVEN a message with checksum that is just "0x" with no hex digits
    std::string_view const message = "JS:1.57,-0.785:0x\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse should fail (currently might succeed or crash)
    EXPECT_FALSE(state.has_value());
}

TEST(BugTests, MessageWithNoSectionDelimiterForChecksum) {
    // BUG: rfind can return npos, which causes invalid substr call
    // GIVEN a malformed message with no final section delimiter
    std::string_view const message = "JS1.57,-0.785";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse should fail gracefully (not crash)
    EXPECT_FALSE(state.has_value());
}

TEST(BugTests, EmptyTokenInSplit) {
    // BUG: split_to_vector may create invalid string_view from empty range
    // GIVEN a message with consecutive delimiters creating empty tokens
    std::string_view const message = "JS::0xFD\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse should fail gracefully (not crash with UB)
    EXPECT_FALSE(state.has_value());
}

TEST(BugTests, AngleWithTrailingGarbage) {
    // BUG: parse_angle doesn't verify entire string was consumed by from_chars
    // GIVEN angles with trailing non-numeric characters after valid numbers
    // Note: checksum 0x8E is correct for "JS:1.57xyz,-0.785abc"
    std::string_view const message = "JS:1.57xyz,-0.785abc:0x8E\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse should fail (currently succeeds by ignoring trailing chars)
    EXPECT_FALSE(state.has_value());
}

TEST(BugTests, ChecksumStringTooLarge) {
    // BUG: validate_checksum compares unsigned int to uint8_t
    // GIVEN a message with a checksum value > 255 (0xFF)
    std::string_view const message = "JS:1.0,2.0:0x1FF\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse should fail (0x1FF can't fit in 8-bit checksum)
    EXPECT_FALSE(state.has_value());
}

// Tests
TEST(ParseAngleTest, ValidPositiveAngle) {
    // GIVEN a valid positive angle string
    std::string_view const angle_str = "1.57";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse succeeds and returns the correct value
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(*result, 1.57, 0.001);
}

TEST(ParseAngleTest, ValidNegativeAngle) {
    // GIVEN a valid negative angle string
    std::string_view const angle_str = "-0.785";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse succeeds and returns the correct value
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(*result, -0.785, 0.001);
}

TEST(ParseAngleTest, ValidZeroAngle) {
    // GIVEN a zero angle string
    std::string_view const angle_str = "0.0";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse succeeds and returns zero
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(*result, 0.0, 0.001);
}

TEST(ParseAngleTest, ValidIntegerAngle) {
    // GIVEN an integer angle string (no decimal point)
    std::string_view const angle_str = "3";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse succeeds and returns the correct value
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(*result, 3.0, 0.001);
}

TEST(ParseAngleTest, InvalidNonNumericString) {
    // GIVEN a non-numeric string
    std::string_view const angle_str = "not_a_number";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse fails
    EXPECT_FALSE(result.has_value());
}

TEST(ParseAngleTest, InvalidEmptyString) {
    // GIVEN an empty string
    std::string_view const angle_str = "";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse fails
    EXPECT_FALSE(result.has_value());
}

TEST(ParseAngleTest, PartialNumberIsRejected) {
    // GIVEN a string with valid number followed by invalid characters
    // Note: We now require the entire string to be consumed
    std::string_view const angle_str = "1.23abc";

    // WHEN parsing the angle
    auto const result = parse_angle(angle_str);

    // THEN the parse fails due to trailing garbage
    EXPECT_FALSE(result.has_value());
}

TEST(SerialProtocolTest, ValidMessage) {
    // GIVEN a valid protocol message with correct header, angles, and checksum
    std::string_view const message = "JS:1.57,-0.785:0xFD\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse succeeds and returns correct joint angles
    ASSERT_TRUE(state.has_value());
    EXPECT_NEAR(state->joint1_rad, 1.57, 0.001);
    EXPECT_NEAR(state->joint2_rad, -0.785, 0.001);
}

TEST(SerialProtocolTest, ValidMessageNoNewline) {
    // GIVEN a valid protocol message without newline terminator
    std::string_view const message = "JS:0.0,3.14159:0xF6";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse succeeds and returns correct joint angles
    ASSERT_TRUE(state.has_value());
    EXPECT_NEAR(state->joint1_rad, 0.0, 0.001);
    EXPECT_NEAR(state->joint2_rad, 3.14159, 0.001);
}

TEST(SerialProtocolTest, InvalidHeader) {
    // GIVEN a message with incorrect header "XX" instead of "JS"
    std::string_view const message = "XX:1.0,2.0:0x50\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse fails
    EXPECT_FALSE(state.has_value());
}

TEST(SerialProtocolTest, InvalidChecksum) {
    // GIVEN a message with incorrect checksum
    std::string_view const message = "JS:1.0,2.0:0xFF\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse fails
    EXPECT_FALSE(state.has_value());
}

TEST(SerialProtocolTest, MalformedData) {
    // GIVEN a message with non-numeric angle data
    std::string_view const message = "JS:not_a_number,2.0:0x50\n";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse fails
    EXPECT_FALSE(state.has_value());
}

TEST(SerialProtocolTest, TrimWhitespace) {
    // GIVEN a valid message with leading and trailing whitespace
    std::string_view const message = "  JS:1.57,-0.785:0xFD\n  ";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse succeeds after trimming whitespace
    ASSERT_TRUE(state.has_value());
    EXPECT_NEAR(state->joint1_rad, 1.57, 0.001);
    EXPECT_NEAR(state->joint2_rad, -0.785, 0.001);
}

TEST(SerialProtocolTest, CustomDelimiters) {
    // GIVEN a message with custom delimiters and a custom protocol config
    std::string_view const message = "JS|1.57;-0.785|0x4E";
    protocol_config const custom{.section_delimiter = '|', .angle_delimiter = ';'};

    // WHEN parsing the message with custom config
    auto const state = read_joint_state(message, custom);

    // THEN the parse succeeds with custom delimiters
    ASSERT_TRUE(state.has_value());
    EXPECT_NEAR(state->joint1_rad, 1.57, 0.001);
    EXPECT_NEAR(state->joint2_rad, -0.785, 0.001);
}

TEST(SerialProtocolTest, WrongExpectedParts) {
    // GIVEN a valid message and a config expecting wrong number of parts
    std::string_view const message = "JS:1.57,-0.785:0xFD";
    protocol_config const wrong_config{
        .expected_parts = 4  // Expects 4 parts but message has 3
    };

    // WHEN parsing with misconfigured expected parts
    auto const state = read_joint_state(message, wrong_config);

    // THEN the parse fails
    EXPECT_FALSE(state.has_value());
}

TEST(SerialProtocolTest, WrongNumberOfAngles) {
    // GIVEN a message with 3 angles instead of the expected 2
    std::string_view const message = "JS:1.57,-0.785,2.1:0x5E";

    // WHEN parsing the message
    auto const state = read_joint_state(message);

    // THEN the parse fails because joint_state::num_joints is 2
    EXPECT_FALSE(state.has_value());
}

TEST(SerialProtocolTest, WrongHeader) {
    // GIVEN a message with "JS" header and a config expecting "XY" header
    std::string_view const message = "JS:1.57,-0.785:0xFD";
    protocol_config const wrong_config{
        .header = "XY"  // Expects "XY" header but message has "JS"
    };

    // WHEN parsing with misconfigured header
    auto const state = read_joint_state(message, wrong_config);

    // THEN the parse fails
    EXPECT_FALSE(state.has_value());
}
