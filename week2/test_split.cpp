#include <gtest/gtest.h>

#include <ranges>
#include <string>
#include <string_view>
#include <vector>

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
std::vector<std::string_view> split(std::string_view const str, char const delimiter) {
    // Create a lazy view pipeline: split by delimiter, then convert each subrange to string_view
    // views::split produces subranges, which we transform into string_views (no allocation yet)
    // Note: split_view must be mutable for std::ranges::begin/end
    auto split_view = str | std::views::split(delimiter) | std::views::transform([](auto&& rng) {
                          // Handle empty ranges safely
                          if (std::ranges::empty(rng)) {
                              return std::string_view{};
                          }
                          return std::string_view(&*rng.begin(), std::ranges::distance(rng));
                      });

    // Materialize the view into a vector - this is where allocation happens
    // The vector constructor iterates the lazy view and copies the string_view objects
    return std::vector<std::string_view>(std::ranges::begin(split_view),
                                         std::ranges::end(split_view));
}

// Tests
TEST(SplitTest, BasicExample) {
    // GIVEN a comma-separated string
    std::string_view const input = "one,two,three";

    // WHEN splitting by comma
    auto const result = split(input, ',');

    // THEN we get three tokens
    ASSERT_EQ(result.size(), 3);
    EXPECT_EQ(result.at(0), "one");
    EXPECT_EQ(result.at(1), "two");
    EXPECT_EQ(result.at(2), "three");
}

TEST(SplitTest, ConsecutiveDelimiters) {
    // GIVEN a string with consecutive delimiters
    std::string_view const input = "a,,b";

    // WHEN splitting by comma
    auto const result = split(input, ',');

    // THEN we get empty token in the middle
    ASSERT_EQ(result.size(), 3);
    EXPECT_EQ(result.at(0), "a");
    EXPECT_EQ(result.at(1), "");
    EXPECT_EQ(result.at(2), "b");
}

TEST(SplitTest, EmptyString) {
    // GIVEN an empty string
    std::string_view const input = "";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN we get zero tokens (std::views::split behavior)
    ASSERT_EQ(result.size(), 0);
}

TEST(SplitTest, NoDelimiter) {
    // GIVEN a string with no delimiter
    std::string_view const input = "nodelmiter";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN we get the whole string as one token
    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result.at(0), "nodelmiter");
}

TEST(SplitTest, OnlyDelimiter) {
    // GIVEN a string with only delimiter
    std::string_view const input = ",";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN we get two empty tokens
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result.at(0), "");
    EXPECT_EQ(result.at(1), "");
}

TEST(SplitTest, LeadingDelimiter) {
    // GIVEN a string starting with delimiter
    std::string_view const input = ",abc";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN first token is empty
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result.at(0), "");
    EXPECT_EQ(result.at(1), "abc");
}

TEST(SplitTest, TrailingDelimiter) {
    // GIVEN a string ending with delimiter
    std::string_view const input = "abc,";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN last token is empty
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result.at(0), "abc");
    EXPECT_EQ(result.at(1), "");
}

TEST(SplitTest, DifferentDelimiter) {
    // GIVEN a pipe-separated string
    std::string_view const input = "a|b|c";

    // WHEN splitting by pipe
    auto const result = split(input, '|');

    // THEN we get three tokens
    ASSERT_EQ(result.size(), 3);
    EXPECT_EQ(result.at(0), "a");
    EXPECT_EQ(result.at(1), "b");
    EXPECT_EQ(result.at(2), "c");
}

TEST(SplitTest, SingleCharacter) {
    // GIVEN a single character
    std::string_view const input = "x";

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN we get one token
    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result.at(0), "x");
}

TEST(SplitTest, ZeroCopyBehavior) {
    // GIVEN a string
    std::string const original = "foo,bar,baz";
    std::string_view const input = original;

    // WHEN splitting
    auto const result = split(input, ',');

    // THEN the string_views point to the original data (zero-copy)
    ASSERT_EQ(result.size(), 3);
    EXPECT_EQ(result.at(0).data(), original.data());
    EXPECT_EQ(result.at(1).data(), original.data() + 4);  // After "foo,"
    EXPECT_EQ(result.at(2).data(), original.data() + 8);  // After "foo,bar,"
}
