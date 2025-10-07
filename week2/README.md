# Week 2: `std::string` and `std::string_view` - Exercises

**Recommended:** Use [Godbolt.org](https://godbolt.org) to complete these exercises! Each exercise includes a pre-configured Godbolt link with starter code and tests.

Once you're ready to run the full test suite locally, see the [Local Testing Guide](TESTING.md).

---

## Exercise 1: String Manipulation Warm-up (Easy)

Write a function `std::string normalize_whitespace(std::string_view const input)` that:
- Removes leading and trailing whitespace
- Replaces all sequences of whitespace characters with a single space
- Returns the normalized string

**Example:**
```cpp
normalize_whitespace("  hello    world  \t\n")
// Returns: "hello world"
```

**Test in Godbolt:**
https://godbolt.org/z/zqEe6af43
```cpp
#include <cassert>
#include <iostream>
#include <string>
#include <string_view>

// Your solution here

int main() {
    assert(normalize_whitespace("  hello    world  ") == "hello world");
    assert(normalize_whitespace("\t\n  test  \n") == "test");
    assert(normalize_whitespace("already normalized") == "already normalized");
    assert(normalize_whitespace("   ") == "");

    std::cout << "All tests passed!\n";
    return 0;
}
```

---

## Exercise 2: String Parsing Warm-up (Medium)

Write a function `std::vector<std::string_view> split(std::string_view const str, char const delimiter)` that:
- Splits a string by a delimiter character
- Returns a vector of string views (no copying!)
- Handles consecutive delimiters (should produce empty views)

**Example:**
```cpp
split("one,two,three", ',')
// Returns: ["one", "two", "three"]

split("a,,b", ',')
// Returns: ["a", "", "b"]
```

**Test in Godbolt:**
https://godbolt.org/z/ccooEMxnd
```cpp
#include <string>
#include <string_view>
#include <vector>
#include <iostream>
#include <cassert>

// Your solution here

int main() {
    auto const result = split("one,two,three", ',');
    assert(result.size() == 3);
    assert(result[0] == "one");
    assert(result[1] == "two");
    assert(result[2] == "three");

    auto const result2 = split("a,,b", ',');
    assert(result2.size() == 3);
    assert(result2[1] == "");

    std::cout << "All tests passed!\n";
    return 0;
}
```

---

## Exercise 3: Serial Protocol Parser - Read Joint State (Hard)

### Protocol Specification

Your robot sends joint state data over a serial port using this protocol:

**Format:** `JS:<angle1>,<angle2>:<checksum>\n`

- `JS:` - Header (2 bytes + colon)
- `<angle1>` - First joint angle in radians (floating point, e.g., "1.57")
- `,` - Delimiter
- `<angle2>` - Second joint angle in radians (floating point, e.g., "-0.785")
- `:` - Checksum delimiter
- `<checksum>` - Simple 8-bit checksum (sum of all ASCII byte values from the start of the message, including the "JS:" header, up to but not including the last colon `:`, then modulo 256, formatted as hex with 0x prefix)
- `\n` - Newline terminator

**Example message:** `JS:1.57,-0.785:0xFD\n`

### Your Task

Write a function:
```cpp
struct joint_state {
    double joint1_rad;
    double joint2_rad;
};

std::optional<joint_state> read_joint_state(std::string_view const message);
```

The function should:
1. Validate the header starts with `"JS:"`
2. Parse the two floating-point joint angles
3. Validate the checksum
4. Return `std::nullopt` if any validation fails
5. Use modern C++20 features where appropriate

**Checksum calculation:** Sum all ASCII byte values from the start of the message (including the "JS:" header) up to but not including the last `:` character (the checksum delimiter), then take modulo 256. For example, in `JS:1.57,-0.785:0x4A`, you checksum the string `"JS:1.57,-0.785"`.

**Test in Godbolt:**
https://godbolt.org/z/o8hex3eqc
```cpp
#include <string>
#include <string_view>
#include <optional>
#include <iostream>
#include <cassert>
#include <charconv>
#include <cmath>

struct joint_state {
    double joint1_rad;
    double joint2_rad;
};

// Helper function to calculate checksum
uint8_t calculate_checksum(std::string_view const data) {
    uint32_t sum = 0;
    for (char const c : data) {
        sum += static_cast<uint8_t>(c);
    }
    return static_cast<uint8_t>(sum % 256);
}

// Your solution here
std::optional<joint_state> read_joint_state(std::string_view const message) {
    // TODO: Implement this
}

int main() {
    // Valid message
    auto const state1 = read_joint_state("JS:1.57,-0.785:0xFD\n");
    assert(state1.has_value());
    assert(std::abs(state1->joint1_rad - 1.57) < 0.001);
    assert(std::abs(state1->joint2_rad - (-0.785)) < 0.001);

    // Another valid message
    auto const state2 = read_joint_state("JS:0.0,3.14159:0xF6\n");
    assert(state2.has_value());
    assert(std::abs(state2->joint1_rad - 0.0) < 0.001);
    assert(std::abs(state2->joint2_rad - 3.14159) < 0.001);

    // Invalid header
    auto const state3 = read_joint_state("XX:1.0,2.0:0x50\n");
    assert(!state3.has_value());

    // Invalid checksum
    auto const state4 = read_joint_state("JS:1.0,2.0:0xFF\n");
    assert(!state4.has_value());

    // Malformed data
    auto const state5 = read_joint_state("JS:not_a_number,2.0:0x50\n");
    assert(!state5.has_value());

    std::cout << "All tests passed!\n";
    return 0;
}
```

**Hints:**
- Use `std::string_view::starts_with()` (C++20)
- Use `std::string_view::find()` to locate delimiters
- Use `std::from_chars()` or `std::stod()` for parsing numbers
- Consider using `std::string_view::substr()` to extract portions
- Think about how to cleanly handle each failure case

---

## Exercise 4: Serial Protocol Writer - Command Joint State (Hard)

### Your Task

Now implement the inverse operation - create a message to send joint commands to the robot.

Write a function:
```cpp
std::string command_joint_state(joint_state const& target);
```

The function should:
1. Format the joint angles with appropriate precision (6 decimal places)
2. Calculate and append the correct checksum in hex format (0xXX)
3. Add the newline terminator
4. Return the complete protocol string

**Test in Godbolt:**
https://godbolt.org/z/KjxePWPba
```cpp
#include <string>
#include <string_view>
#include <optional>
#include <iostream>
#include <cassert>
#include <format>  // C++20
#include <sstream>
#include <iomanip>

struct joint_state {
    double joint1_rad;
    double joint2_rad;
};

// Helper function to calculate checksum
uint8_t calculate_checksum(std::string_view const data) {
    uint32_t sum = 0;
    for (char const c : data) {
        sum += static_cast<uint8_t>(c);
    }
    return static_cast<uint8_t>(sum % 256);
}

// Your solution here
std::string command_joint_state(joint_state const& target) {
    // TODO: Implement this
}

int main() {
    joint_state const cmd1{1.570796, -0.785398};
    std::string const msg1 = command_joint_state(cmd1);

    std::cout << std::format("Generated message: {}", msg1);

    // Should start with header
    assert(msg1.starts_with("JS:"));

    // Should end with newline
    assert(msg1.ends_with("\n"));

    // Should contain both angles
    assert(msg1.find("1.570796") != std::string::npos);
    assert(msg1.find("-0.785398") != std::string::npos);

    // Should have checksum
    assert(msg1.find("0x") != std::string::npos);

    // Verify it round-trips with read_joint_state from Exercise 3
    // (Uncomment if you've completed Exercise 3)
    // auto const parsed = read_joint_state(msg1);
    // assert(parsed.has_value());
    // assert(std::abs(parsed->joint1_rad - cmd1.joint1_rad) < 0.0001);

    std::cout << "All tests passed!\n";
    return 0;
}
```

**Hints:**
- Use `std::format()` (C++20) or `std::stringstream` with `std::setprecision()`
- Calculate checksum on everything before the last `:` (including "JS:" header and the angle data)
- Format checksum as hex with `std::format("{:02X}", checksum)` or `std::hex`

---

## Learning Objectives

After completing these exercises, you should understand:

1. **`std::string_view` benefits**: Zero-copy string operations, perfect for parsing
2. **String searching**: `find()`, `starts_with()`, `ends_with()` methods
3. **String manipulation**: `substr()`, `remove_prefix()`, `remove_suffix()`
4. **Modern parsing**: `std::from_chars()` for safe, fast conversions
5. **Error handling**: Using `std::optional` for parse failures
6. **Formatting**: `std::format()` (C++20) vs older `iostream` methods
7. **Ranges**: C++20 ranges for more expressive data transformations
8. **Real-world protocol design**: Headers, checksums, delimiters, error handling

## Discussion Questions

1. Why use `std::string_view` for the parsing functions instead of `const std::string&`?
2. What are the lifetime considerations when returning `std::vector<std::string_view>`?
3. How would you modify the protocol to handle variable numbers of joints?
4. What are the trade-offs between `std::from_chars()` and `std::stod()`?
5. Why is `std::format()` preferred over `stringstream` in modern C++?
6. How does the checksum provide data integrity? What attacks does it NOT protect against?

---

## Solutions

✅ **[View Solutions with Detailed Explanations](week2_solutions.md)**
