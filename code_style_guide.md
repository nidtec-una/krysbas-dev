# Code Style Guide for KrySBAS

## Table of Contents
- [Motivation](#motivation)
- [Mandatory Rules](#mandatory-rules)
- [Highly Encouraged Rules](#highly-encouraged-rules)
- [Optional Rules and Good Practices](#optional-rules-and-good-practices)


## Motivation

Clean code is fundamental to develop human-readable, maintainable, and extendable code. This note is meant to be used as a guide for developing new code as part of the KrySBAS package. The adopted conventions and styles are largely based on the [MATLAB Style Guidelines 2.0](https://www.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0) (MATLAB Central File Exchange. Retrieved August 16, 2023.) with additional internal modifications/adaptions.

## Mandatory Rules

1. Variables must be named following the camelCase convention, e.g., `rowNumber`.
2. Constants must be named in uppercase with words separated by underscore, e.g., `MAX_ITERATIONS`.
3. Structures must be named following the CamelCase convention, e.g., `SolverProperties`.
4. Functions must be named following the snake_case convention, e.g., `restarted_gmres`.
5. Inline comments must be included when variables and constants are first introduced in the code, e.g., `MAX_ITERATIONS = 10; % maximum number of iterations`.
6. Code content must fit within the first 80 columns.
7. Indentation must occupy exactly 4 spaces.
8. White spaces must surround assignment and conditional operators, e.g., `=`, `&&`, `||`.
9. White spaces must surround binary operators, e.g., `+`, `-`, `/`, etc.
10. White spaces must surround commas, e.g., `foo(alpha, beta, gamma)`.
11. Comments and documentation must be written in English.
12. Header comments must be included in all functions and subfunctions (TODO: Decide on the specific style and syntax).
13. Built-in functions must be used whenever possible.
14. If a block of code appears in more than one `.m` file, it must be encapsulated into a function to avoid code duplication.
15. Subfunctions must be used only when they are associated with a single `.m` file. Subfunctions cannot be called externally.

## Highly Encouraged Rules

1. Number of objects should have meaningful names (e.g., `nRows` not `n`), unless context allows it.
2. Prefix iterators should be used (e.g., `iFile` not `i`) to avoid potential confusion, unless context allows it.
3. Loop result variables should be initialized before the loop starts.

## Optional Rules and Good Practices

1. Avoid using `i` or `j` as iterators (this overwrites the built-in variables for complex numbers in MATLAB).
2. Try writing comments while you code.
3. Try documenting how and why, not so much what the code does.
4. Partition/encapsulate code longer than two editor screens.
5. Instead of passing long lists of inputs or returning long lists of outputs, consider using structures.
6. Matrix names in uppercase and vector names in lowercase (unless they have a meaningful name, in which case they follow rule #1).
