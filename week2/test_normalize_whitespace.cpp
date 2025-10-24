#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <string>
#include <string_view>

/**
 * @brief Normalizes whitespace in a string.
 *
 * Removes leading and trailing whitespace and collapses all sequences
 * of whitespace characters into a single space.
 *
 * @param input The input string view
 * @return std::string The normalized string
 */
std::string normalize_whitespace(std::string_view const input) {
    // Compile-time constant for all whitespace characters
    constexpr std::string_view whitespace = " \t\n\r\f\v";

    // Trim leading/trailing whitespace using string_view methods
    auto trimmed = input;
    auto const start = trimmed.find_first_not_of(whitespace);
    if (start == std::string_view::npos)
        return "";  // All whitespace

    trimmed.remove_prefix(start);
    auto const end = trimmed.find_last_not_of(whitespace);
    trimmed.remove_suffix(trimmed.size() - end - 1);

    // Create a lazy view that normalizes all whitespace to spaces (tabs, newlines, etc.)
    // This view doesn't allocate - it computes values on-the-fly during iteration
    auto const normalized =
        trimmed | std::views::transform([](char const c) { return std::isspace(c) ? ' ' : c; });

    // Use unique_copy to collapse consecutive spaces and write to result
    // This is where the actual allocation happens (via back_inserter)
    std::string result;
    std::ranges::unique_copy(normalized, std::back_inserter(result),
                             [](char const a, char const b) { return a == ' ' and b == ' '; });

    return result;
}

// Tests
TEST(NormalizeWhitespaceTest, BasicExample) {
    // GIVEN a string with leading/trailing whitespace and multiple internal spaces
    std::string_view const input = "  hello    world  ";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN leading/trailing is removed and internal spaces collapsed
    EXPECT_EQ(result, "hello world");
}

TEST(NormalizeWhitespaceTest, TabsAndNewlines) {
    // GIVEN a string with tabs and newlines
    std::string_view const input = "\t\n  test  \n";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN all whitespace types are normalized
    EXPECT_EQ(result, "test");
}

TEST(NormalizeWhitespaceTest, AlreadyNormalized) {
    // GIVEN an already normalized string
    std::string_view const input = "already normalized";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN the string is unchanged
    EXPECT_EQ(result, "already normalized");
}

TEST(NormalizeWhitespaceTest, OnlyWhitespace) {
    // GIVEN a string with only whitespace
    std::string_view const input = "   ";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN the result is empty
    EXPECT_EQ(result, "");
}

TEST(NormalizeWhitespaceTest, EmptyString) {
    // GIVEN an empty string
    std::string_view const input = "";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN the result is empty
    EXPECT_EQ(result, "");
}

TEST(NormalizeWhitespaceTest, NoWhitespace) {
    // GIVEN a string with no whitespace
    std::string_view const input = "nowhitespace";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN the string is unchanged
    EXPECT_EQ(result, "nowhitespace");
}

TEST(NormalizeWhitespaceTest, MixedWhitespaceTypes) {
    // GIVEN a string with various whitespace types mixed together
    std::string_view const input = " \t hello \n\r world \f\v test ";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN all are normalized to single spaces
    EXPECT_EQ(result, "hello world test");
}

TEST(NormalizeWhitespaceTest, SingleWord) {
    // GIVEN a single word with whitespace around it
    std::string_view const input = "  word  ";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN only the word remains
    EXPECT_EQ(result, "word");
}

TEST(NormalizeWhitespaceTest, MultipleConsecutiveSpaces) {
    // GIVEN a string with many consecutive spaces
    std::string_view const input = "a          b";

    // WHEN normalizing whitespace
    auto const result = normalize_whitespace(input);

    // THEN collapsed to single space
    EXPECT_EQ(result, "a b");
}
