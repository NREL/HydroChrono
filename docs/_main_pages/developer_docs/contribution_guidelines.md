---
layout: page
title: Contribution Guidelines
parent_section: Developer Documentation
---

# Contribution Guidelines

## Introduction

Contributions to the HydroChrono project are highly valued and appreciated. This section outlines the guidelines and best practices to follow when contributing. We follow the Google C++ Style Guide; below are some of its main points.

## Google C++ Style Guide Highlights

1. **General Naming Rules**:
   - Opt for clarity and completeness over brevity.
   - Use descriptive names, and avoid abbreviations unless they are more common than the full word.

2. **File Names**:
   - All lowercase, with words separated by underscores (`_`), e.g., `my_example_file.cpp`.

3. **Type Names**:
   - Camel case starting with an uppercase letter, e.g., `MyExampleClass`.

4. **Variable Names**:
   - Use snake case (all lowercase with underscores), e.g., `my_example_variable`.

5. **Class Member Variables**:
   - These have a trailing underscore, e.g., `class_member_variable_`.

6. **Function Names**:
   - Camel case starting with an uppercase letter, e.g., `DoSomethingInteresting()`.

7. **Constant Names**:
   - Use a `k` prefix followed by Camel case, e.g., `const int kDaysInAWeek = 7;`.

8. **Namespace Names**:
   - All lowercase, and based on project name or what makes sense for your project.
   - **Example**:

     ```cpp
     namespace hydrochrono { 
         // ... code ...
     }
     ```

9. **Enumerator Names**:
   - Use Pascal case for the enumeration type name. The individual enumerators within it should use Pascal case or Camel case depending on context.
   - **Example**:

     ```cpp
     enum class WaveType { 
         OceanWave, 
         RiverFlow, 
         Tsunami 
     };

     enum class colorOptions {
         redOption,
         blueOption,
         greenOption
     };
     ```

10. **Using `auto`**:
   - Use `auto` to avoid type names that are noisy, obvious, or unimportant.
   - Never use `auto` for defining a type, which isn't obvious from the context.

## Code Standards

- Adhere to the points mentioned above and the full Google C++ Style Guide.
- Provide comprehensive comments for any non-trivial code.
- Test your code to ensure existing functionality remains intact.

## Submitting Contributions

1. **Fork and Clone**:
   - Fork the HydroChrono repository on the platform where it's hosted (e.g., GitHub).
   - Clone your forked repository to your local machine.

2. **Create a New Branch**:
   - Create a new branch for each new feature or bugfix. Naming it descriptively can help, e.g., `feature-add-hydro-dynamics` or `bugfix-memory-leak`.

3. **Committing Changes**:
   - Write clear, concise commit messages that explain the change.
   - Split larger changes into multiple commits if possible.

4. **Syncing with Upstream**:
   - Regularly sync your fork and branch with the main repository (`upstream`) to keep up with changes.
   - Merge or rebase your branch with the latest from `upstream` before submitting a pull request.

5. **Code Reviews**:
   - Once you've pushed your branch to your fork, submit a pull request to the main HydroChrono repository.
   - Participate in the code review process. Respond to feedback, and make changes as requested.

6. **Testing**:
   - Ensure that your code passes all tests and doesn't introduce new issues.
   - Add new tests for your features to ensure future changes don't break your contribution.

7. **Documentation**:
   - Update relevant documentation pertaining to your changes.
   - Ensure examples, if provided, are clear and understandable.

## Contribution Best Practices

- **Collaboration**: Encourage collaboration with other contributors. Open discussions can lead to better solutions.
- **Stay Updated**: Regularly pull the latest changes from the `upstream` main branch and stay updated with the project's progress and changes.
- **Respect**: Respect the decisions of maintainers and the feedback of others. Every piece of feedback is to ensure the quality and consistency of the project.

## Conclusion

Whether you're adding new features, bug fixes, or simply improving documentation, your contributions are invaluable and greatly appreciated.
