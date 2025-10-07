# Week 2: `std::string` and `std::string_view` - Solutions

## Solution 1: Normalize Whitespace

### Simple Approach (Using loops and conditionals)

```cpp
std::string normalize_whitespace(std::string_view const input) {
    std::string result;
    bool in_whitespace = true;  // Start true to skip leading spaces

    for (char const c : input) {
        if (std::isspace(c)) {
            if (!in_whitespace) {
                result += ' ';
                in_whitespace = true;
            }
        } else {
            result += c;
            in_whitespace = false;
        }
    }

    // Remove trailing space if we added one
    if (!result.empty() && result.back() == ' ') {
        result.pop_back();
    }

    return result;
}
```

### Modern C++20 Approach

```cpp
#include <algorithm>
#include <ranges>

std::string normalize_whitespace(std::string_view const input) {
    // Compile-time constant for all whitespace characters
    constexpr std::string_view whitespace = " \t\n\r\f\v";

    // Trim leading/trailing whitespace using string_view methods (clearer than double-reverse)
    auto trimmed = input;
    auto const start = trimmed.find_first_not_of(whitespace);
    if (start == std::string_view::npos) return "";  // All whitespace

    trimmed.remove_prefix(start);
    auto const end = trimmed.find_last_not_of(whitespace);
    trimmed.remove_suffix(trimmed.size() - end - 1);

    // Create a lazy view that normalizes all whitespace to spaces (tabs, newlines, etc.)
    // This view doesn't allocate - it computes values on-the-fly during iteration
    auto const normalized = trimmed
        | std::views::transform([](char const c) { return std::isspace(c) ? ' ' : c; });

    // Use unique_copy to collapse consecutive spaces and write to result
    // This is where the actual allocation happens (via back_inserter)
    std::string result;
    std::ranges::unique_copy(normalized, std::back_inserter(result),
        [](char const a, char const b) { return a == ' ' and b == ' '; });

    return result;
}
```

---

## Solution 2: Split String

### Simple Approach (Using find in a loop)

```cpp
std::vector<std::string_view> split(std::string_view const str, char const delimiter) {
    std::vector<std::string_view> result;
    size_t start = 0;

    while (start < str.size()) {
        size_t const end = str.find(delimiter, start);

        if (end == std::string_view::npos) {
            // Last token
            result.push_back(str.substr(start));
            break;
        } else {
            result.push_back(str.substr(start, end - start));
            start = end + 1;
        }
    }

    // Handle trailing delimiter
    if (!str.empty() && str.back() == delimiter) {
        result.push_back(std::string_view{});
    }

    return result;
}
```

### Modern C++20 Approach (Using ranges)

```cpp
#include <ranges>

std::vector<std::string_view> split(std::string_view const str, char const delimiter) {
    // Create a lazy view pipeline: split by delimiter, then convert each subrange to string_view
    // views::split produces subranges, which we transform into string_views (no allocation yet)
    auto const split_view = str
        | std::views::split(delimiter)
        | std::views::transform([](auto&& rng) {
            // auto&& is a forwarding reference - accepts both lvalues and rvalues
            // rng.begin() returns an iterator to the first char
            // *rng.begin() dereferences to get the char
            // &*rng.begin() gets the address of that char (pointer to start of subrange)
            return std::string_view(&*rng.begin(), std::ranges::distance(rng));
        });

    // Materialize the view into a vector - this is where allocation happens
    // The vector constructor iterates the lazy view and copies the string_view objects
    return std::vector<std::string_view>(
        std::ranges::begin(split_view),
        std::ranges::end(split_view)
    );
}
```

---

## Solution 3: Read Joint State

### Simple Approach (Manual parsing)

```cpp
std::optional<joint_state> read_joint_state(std::string_view const message) {
    // Check header
    if (!message.starts_with("JS:")) {
        return std::nullopt;
    }

    // Find delimiters
    size_t const comma_pos = message.find(',', 3);
    if (comma_pos == std::string_view::npos) {
        return std::nullopt;
    }

    size_t const colon_pos = message.find(':', comma_pos);
    if (colon_pos == std::string_view::npos) {
        return std::nullopt;
    }

    // Extract substrings
    std::string_view const angle1_str = message.substr(3, comma_pos - 3);
    std::string_view const angle2_str = message.substr(comma_pos + 1, colon_pos - comma_pos - 1);
    std::string_view checksum_str = message.substr(colon_pos + 1);

    // Remove newline from checksum if present
    if (checksum_str.ends_with('\n')) {
        checksum_str.remove_suffix(1);
    }

    // Parse angles
    double joint1, joint2;
    try {
        joint1 = std::stod(std::string(angle1_str));
        joint2 = std::stod(std::string(angle2_str));
    } catch (...) {
        return std::nullopt;
    }

    // Validate checksum - checksum is calculated on everything before the last colon
    std::string_view const data_for_checksum = message.substr(0, colon_pos);
    uint8_t const calculated = calculate_checksum(data_for_checksum);

    // Parse provided checksum (hex format 0xXX)
    if (!checksum_str.starts_with("0x") || checksum_str.size() < 3) {
        return std::nullopt;
    }

    unsigned int provided_checksum;
    auto const result = std::from_chars(
        checksum_str.data() + 2,
        checksum_str.data() + checksum_str.size(),
        provided_checksum,
        16
    );

    if (result.ec != std::errc{} || provided_checksum != calculated) {
        return std::nullopt;
    }

    return joint_state{joint1, joint2};
}
```

### Modern C++20 Approach (Production-Ready)

```cpp
/**
 * @brief Represents the state of a robot's joints in radians.
 */
struct joint_state {
    static constexpr size_t num_joints = 2;  ///< Number of joints in this state
    double joint1_rad;  ///< First joint angle in radians
    double joint2_rad;  ///< Second joint angle in radians
};

/**
 * @brief Configuration parameters for the serial protocol parser.
 */
struct protocol_config {
    char section_delimiter = ':';   ///< Delimiter between protocol sections (header, data, checksum)
    char angle_delimiter = ',';     ///< Delimiter between individual angle values
    std::string_view header = "JS"; ///< Protocol header identifier
    size_t expected_parts = 3;      ///< Expected number of sections after splitting by section_delimiter
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
    if (start == std::string_view::npos) return "";  // All whitespace

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
    // Note: split_view must be mutable for std::ranges::begin/end
    auto split_view = str
        | std::views::split(delimiter)
        | std::views::transform([](auto&& rng) {
            return std::string_view(&*rng.begin(), std::ranges::distance(rng));
        });

    // Materialize the lazy view into a vector
    return std::vector<std::string_view>(
        std::ranges::begin(split_view),
        std::ranges::end(split_view)
    );
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
    auto const [ptr, ec] = std::from_chars(
        angle_str.data(),
        angle_str.data() + angle_str.size(),
        angle
    );

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
bool validate_checksum(std::string_view const message,
                       std::string_view const checksum_str,
                       char const section_delimiter) {
    // Verify checksum format (must start with "0x")
    if (not checksum_str.starts_with("0x")) {
        return false;
    }

    // Parse the provided checksum from hex string (skip "0x" prefix)
    unsigned int provided;
    auto const [ptr, ec] = std::from_chars(
        checksum_str.data() + 2,
        checksum_str.data() + checksum_str.size(),
        provided,
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
    if (not validate_checksum(msg, parts_vec.at(2), config.section_delimiter)) {
        return std::nullopt;
    }

    // Parse and return joint angles
    return parse_joint_angles(parts_vec.at(1), config.angle_delimiter);
}
```

---

## Solution 4: Command Joint State

### Simple Approach (Using stringstream)

```cpp
std::string command_joint_state(joint_state const& target) {
    std::stringstream ss;

    // Format the data portion
    ss << "JS:" << std::fixed << std::setprecision(6)
       << target.joint1_rad << "," << target.joint2_rad;

    std::string const data = ss.str();

    // Calculate checksum on the data portion (everything before the last colon)
    uint8_t const checksum = calculate_checksum(data);

    // Add checksum and newline
    ss << ":0x" << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
       << static_cast<int>(checksum) << "\n";

    return ss.str();
}
```

### Modern C++20 Approach (Using std::format)

```cpp
#include <format>

std::string command_joint_state(joint_state const& target) {
    // Format the data portion with high precision
    std::string const data = std::format("JS:{:.6f},{:.6f}",
                                         target.joint1_rad,
                                         target.joint2_rad);

    // Calculate checksum on data (everything before the last colon)
    uint8_t const checksum = calculate_checksum(data);

    // Append checksum and newline
    return std::format("{}:0x{:02X}\n", data, checksum);
}
```
