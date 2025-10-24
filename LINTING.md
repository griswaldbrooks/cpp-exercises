# Linting and Code Formatting

This repository uses `clang-format` for code formatting and `pre-commit` hooks to ensure code quality.

## Setup

Install the pre-commit hooks (first time only):

```bash
pixi run lint-install
```

## Usage

### Run linter on all files

```bash
pixi run lint
```

### Format all C++ files

```bash
pixi run format
```

### Check formatting without modifying files

```bash
pixi run format-check
```

## Pre-commit Hooks

Once installed, pre-commit hooks will automatically run on staged files when you commit:

- `clang-format`: Formats C++ code according to `.clang-format` rules
- `trailing-whitespace`: Removes trailing whitespace
- `end-of-file-fixer`: Ensures files end with a newline
- `check-yaml`: Validates YAML syntax
- `check-added-large-files`: Prevents committing large files
- `check-merge-conflict`: Detects merge conflict markers

## Configuration

- `.clang-format`: C++ code formatting rules (based on Google style with 4-space indent)
- `.pre-commit-config.yaml`: Pre-commit hook configuration

## Manual Run

To manually run hooks on specific files:

```bash
pre-commit run --files week3/test_transform3d.cpp
```

To skip hooks for a specific commit (not recommended):

```bash
git commit --no-verify
```
