# C++ STL Learning Plan - Weeks 2-8

A structured curriculum for learning the C++ Standard Template Library using C++20, designed for junior developers building on foundational knowledge of `std::vector`.

---

## Week 2: `std::string` and `std::string_view`

### Why Cover This
Strings are fundamental to almost every C++ program. Understanding string handling is critical for text processing, I/O, and API design. `string_view` is a C++17 addition that's essential for modern C++ performance optimization.

### Topics to Cover
- **Construction and Initialization**
  - Different constructor overloads
  - String literals and raw string literals (`R"(...)"`)
  - Move semantics with strings

- **Common Operations**
  - Concatenation (`+`, `+=`, `append()`)
  - Substring extraction (`substr()`)
  - Finding and searching (`find()`, `rfind()`, `find_first_of()`, etc.)
  - Replacement and modification (`replace()`, `insert()`, `erase()`)

- **C++20 Features**
  - `starts_with()` and `ends_with()` for prefix/suffix checking
  - Comparison with string views

- **`std::string_view`**
  - Non-owning view of string data
  - When to use: function parameters, avoiding copies
  - Lifetime considerations and dangling references (critical safety topic)
  - Performance benefits over `const std::string&`

### Exercises
- Implement a simple tokenizer/splitter
- Parse configuration key-value pairs
- Compare performance of `string` vs `string_view` parameters

---

## Week 3: `std::array` and `std::span`

### Why Cover This
`std::array` demonstrates how the STL wraps traditional C constructs with safer, more ergonomic interfaces. `std::span` (C++20) is crucial for modern C++ as it provides a safe, non-owning view of contiguous memory.

### Topics to Cover
- **`std::array` Fundamentals**
  - Compile-time fixed size
  - Stack allocation vs `vector`'s heap allocation
  - STL interface (iterators, `size()`, `at()`, etc.)
  - Structured bindings with arrays

- **When to Use `std::array`**
  - Known, fixed sizes at compile time
  - Performance-critical code (no heap allocation)
  - Returning arrays from functions (vs C-style arrays)

- **`std::span` (C++20)**
  - Non-owning view of contiguous sequences
  - Works with `std::array`, `std::vector`, C-style arrays
  - Dynamic vs fixed extent spans
  - Subspans and slicing

- **Comparison of Containers**
  - `std::array` vs C-style arrays vs `std::vector`
  - Memory layout and performance
  - Type safety and bounds checking

### Exercises
- Convert legacy code using C-style arrays to `std::array`
- Write functions accepting `std::span` that work with multiple container types
- Implement a fixed-size matrix class using `std::array`

---

## Week 4: `std::deque` and `std::list`

### Why Cover This
Understanding different container types and their performance characteristics is essential for writing efficient C++. This week introduces non-contiguous containers and teaches when `vector` isn't the right choice.

### Topics to Cover
- **`std::deque` (Double-Ended Queue)**
  - Internal structure (chunks of contiguous memory)
  - Efficient insertion/deletion at both ends
  - Random access capability
  - When to prefer over `vector`

- **`std::list` (Doubly-Linked List)**
  - Node-based structure
  - Iterator stability (iterators not invalidated on insertion/deletion)
  - Bidirectional but not random access
  - `splice()` for efficient list manipulation
  - When pointer/iterator stability matters

- **Performance Characteristics**
  - Big-O complexity for common operations
  - Cache locality and real-world performance
  - Memory overhead comparison
  - Why `vector` is usually still the default choice

- **Use Cases**
  - `deque`: Queue implementations, sliding window algorithms
  - `list`: When frequent middle insertions/deletions with iterator stability needed

### Exercises
- Implement a task queue using `std::deque`
- Benchmark insertion performance: `vector` vs `deque` vs `list`
- Write an LRU cache using `std::list` (preview of next week's `unordered_map`)

---

## Week 5: `std::map` and `std::unordered_map`

### Why Cover This
Associative containers are fundamental data structures for efficient key-value lookups. Understanding the trade-offs between ordered (tree-based) and unordered (hash-based) containers is crucial for performance-conscious development.

### Topics to Cover
- **`std::map` Fundamentals**
  - Red-black tree implementation
  - Automatically sorted by key
  - Logarithmic time complexity (O(log n))
  - Iterator ordering guarantees

- **`std::unordered_map` Fundamentals**
  - Hash table implementation
  - Constant average-time complexity (O(1))
  - No ordering guarantees
  - Load factor and bucket management

- **Common Operations**
  - Insertion: `insert()`, `emplace()`, `operator[]`, `try_emplace()`
  - Lookup: `find()`, `at()`, `operator[]`, C++20's `contains()`
  - Deletion: `erase()`
  - Iteration patterns

- **Custom Types**
  - Custom comparators for `std::map`
  - Custom hash functions and equality for `std::unordered_map`
  - Using `std::hash` specialization

- **C++20 Features**
  - `contains()` member function (cleaner than `find() != end()`)
  - Improved `erase_if()` support

- **Performance Considerations**
  - When to use `map` vs `unordered_map`
  - Ordered traversal requirements
  - Hash function quality impact

### Exercises
- Implement a word frequency counter
- Create a custom class with both hash function and comparator
- Benchmark lookup performance between `map` and `unordered_map`
- Build a simple symbol table or configuration registry

---

## Week 6: `std::set` and `std::unordered_set`

### Why Cover This
Sets are essential for maintaining unique collections and performing set operations. This week reinforces concepts from maps while introducing mathematical set semantics.

### Topics to Cover
- **`std::set` and `std::unordered_set` Fundamentals**
  - Relationship to `map`/`unordered_map` (just keys, no values)
  - Same tree/hash implementations
  - Automatic uniqueness enforcement
  - Similar performance characteristics to their map counterparts

- **Common Operations**
  - Insertion and uniqueness: `insert()`, `emplace()`
  - Lookup: `find()`, `count()`, `contains()` (C++20)
  - Set algorithms: `std::set_union()`, `std::set_intersection()`, `std::set_difference()`

- **Use Cases**
  - Maintaining unique collections
  - Fast membership testing
  - Mathematical set operations
  - Removing duplicates from sequences

- **`std::multiset` and `std::unordered_multiset`**
  - Allowing duplicate elements
  - When and why to use them
  - `equal_range()` for finding all equivalent elements

### Exercises
- Implement duplicate removal from a `vector`
- Build a spell-checker using a word set
- Solve set theory problems (intersection, union, difference)
- Track unique visitors/IDs in a system

---

## Week 7: Algorithm Basics (`<algorithm>`)

### Why Cover This
Algorithms are the heart of the STL. They promote generic, reusable code and are typically more efficient and correct than hand-written loops. This week introduces the most common algorithms and the concept of working with iterators.

### Topics to Cover
- **Non-Modifying Sequence Operations**
  - `std::find()` - finding elements
  - `std::find_if()` - finding with predicates
  - `std::count()` and `std::count_if()` - counting elements
  - `std::all_of()`, `std::any_of()`, `std::none_of()` - boolean checks

- **Sorting and Searching**
  - `std::sort()` - O(n log n) sorting
  - `std::stable_sort()` - maintaining relative order
  - `std::binary_search()` - O(log n) searching in sorted ranges
  - `std::lower_bound()`, `std::upper_bound()`, `std::equal_range()`
  - Custom comparators

- **Predicates and Lambdas**
  - Function objects and functors
  - Lambda expressions (capture lists, parameters, return types)
  - When to use lambdas vs named functions
  - `std::function` and callable objects

- **Iterator Concepts**
  - Input, output, forward, bidirectional, random-access iterators
  - `begin()` and `end()` idiom
  - Iterator invalidation rules (review from containers)

### Exercises
- Sort a vector of custom structs with multiple sort orders
- Find all even numbers in a container using `std::find_if`
- Implement custom predicates and compare with lambda equivalents
- Use `std::binary_search` on sorted data

---

## Week 8: More Algorithms and C++20 Ranges

### Why Cover This
This week covers transformative algorithms and introduces the powerful ranges library from C++20, which represents a paradigm shift in how we work with sequences and algorithms in modern C++.

### Topics to Cover
- **Modifying Sequence Operations**
  - `std::transform()` - applying functions to ranges
  - `std::copy()`, `std::copy_if()` - copying elements
  - `std::remove()`, `std::remove_if()` - the erase-remove idiom
  - `std::replace()`, `std::replace_if()` - replacing values
  - `std::reverse()`, `std::rotate()` - rearranging sequences

- **Numeric Algorithms (`<numeric>`)**
  - `std::accumulate()` - folding/reducing sequences
  - `std::reduce()` (C++17) - parallel-friendly reduction
  - `std::transform_reduce()` - map-reduce pattern
  - `std::iota()` - generating sequences

- **C++20 Ranges Introduction**
  - What ranges solve: composition and readability
  - Range adapters and views (lazy evaluation)
  - `std::ranges::views::filter`, `views::transform`
  - `std::ranges::views::take`, `views::drop`
  - Range-based algorithms (`std::ranges::sort`, `std::ranges::find`, etc.)

- **Ranges vs Traditional Algorithms**
  - Projection support (transforming during iteration)
  - Composability with `|` (pipe) operator
  - Safer (no iterator pair mismatches)
  - More expressive code

### Exercises
- Use `std::transform()` to convert a vector of strings to uppercase
- Implement the erase-remove idiom to filter a vector
- Calculate sum, product, and average using `std::accumulate`
- Rewrite traditional algorithm code using ranges and views
- Chain multiple range views to process data pipelines

---

## Teaching Notes

### General Approach
- Start each session with a brief review of the previous week
- Use compiler explorer (godbolt.org) to examine generated assembly for performance discussions
- Encourage writing benchmarks to verify performance claims
- Emphasize reading cppreference documentation directly
- Live-code examples and debug together

### Key Themes to Reinforce
- **RAII and ownership**: Who owns the memory?
- **Performance**: When does allocation happen? What's the complexity?
- **Iterator invalidation**: When do iterators become invalid?
- **Modern C++ idioms**: Prefer STL over hand-rolled solutions
- **Type safety**: Compile-time checks over runtime checks

### Progression Logic
1. **Weeks 2-3**: Sequential containers and views (building on `vector`)
2. **Weeks 4**: Alternative sequential containers (understanding trade-offs)
3. **Weeks 5-6**: Associative containers (new access patterns)
4. **Weeks 7-8**: Algorithms (bringing it all together)

This progression builds complexity gradually while reinforcing core concepts like iterators, templates, and performance characteristics throughout.

### Future Topics (Beyond Week 8)
- Smart pointers (`std::unique_ptr`, `std::shared_ptr`, `std::weak_ptr`)
- Function objects and `std::function`
- Iterators and custom iterator design
- More ranges and advanced views
- Concurrent containers and algorithms
- Advanced template metaprogramming with concepts (C++20)
