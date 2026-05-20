# Contributing Guidelines

Thank you for your interest in contributing to the Bulk RNA-Seq Analysis Pipeline! This document provides guidelines for contributing code, documentation, and improvements.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Pull Request Process](#pull-request-process)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)

## Code of Conduct

This project adheres to the Contributor Covenant Code of Conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## How to Contribute

### Ways You Can Help

1. **Report Bugs** - Identify and report issues with the pipeline
2. **Suggest Features** - Propose new functionality or improvements
3. **Write Documentation** - Improve or expand existing documentation
4. **Improve Code** - Fix bugs, optimize performance, or add features
5. **Test** - Help identify edge cases and verify functionality
6. **Share Examples** - Contribute example datasets or use cases

### Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/Bulk-RNA-seq.git
   cd Bulk-RNA-seq
   ```
3. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```
4. **Make your changes** following the guidelines below
5. **Commit with clear messages**:
   ```bash
   git commit -m "Add feature: description of changes"
   ```
6. **Push to your fork** and **create a Pull Request**

## Reporting Bugs

### Before Submitting a Bug Report

- Check the existing issues to avoid duplicates
- Check the documentation and FAQ
- Collect version information and system details

### How to Submit a Bug Report

Open an issue with the following information:

```markdown
**Title:** Brief description of the bug

**Environment:**
- OS: [e.g., Linux, macOS]
- Python/R version: [if applicable]
- Tool versions: [HISAT2, DESeq2, etc.]

**Description:**
Clear description of what the bug is

**Steps to Reproduce:**
1. Step 1
2. Step 2
3. ...

**Expected Behavior:**
What should happen

**Actual Behavior:**
What actually happens

**Sample Data:**
[Attach small example FASTQ or provide SRR accession numbers]

**Relevant Logs:**
```
[Paste any error messages or logs]
```

**Additional Context:**
Any other relevant information
```

## Suggesting Enhancements

### Enhancement Request Format

```markdown
**Title:** Brief description of the enhancement

**Motivation:**
Why would this be useful?

**Proposed Solution:**
How would you implement this?

**Alternatives Considered:**
Other approaches you've thought of

**Additional Context:**
Any other relevant information
```

## Pull Request Process

### Before Creating a PR

1. Update documentation if needed
2. Test your changes thoroughly
3. Follow coding standards (see below)
4. Add comments for complex logic
5. Ensure backward compatibility

### PR Submission Checklist

- [ ] I've tested the changes locally
- [ ] I've updated relevant documentation
- [ ] I've added comments for complex code
- [ ] My code follows the project's style guidelines
- [ ] I haven't introduced new dependencies without discussion
- [ ] My changes don't break existing functionality
- [ ] I've signed off on the code

### PR Template

```markdown
**Description:**
Brief description of the changes

**Related Issue:**
Fixes #[issue number]

**Type of Change:**
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement

**Testing:**
Describe the tests you've run and results

**Documentation:**
- [ ] Updated README
- [ ] Updated code comments
- [ ] Updated configuration examples

**Breaking Changes:**
None / [describe if applicable]

**Checklist:**
- [ ] Code follows style guidelines
- [ ] No new compiler warnings
- [ ] Tests pass locally
- [ ] Changes are backward compatible
```

## Coding Standards

### Bash Scripts

```bash
#!/bin/bash
# Clear description of what the script does

set -e  # Exit on error

# Use meaningful variable names
VARIABLE_NAME="value"

# Add comments for complex logic
if [ condition ]; then
    # Do something
    command
fi

# Use functions for reusable code
function_name() {
    local input="$1"
    # Function body
}

# Error handling
if [ $? -ne 0 ]; then
    echo "Error message" >&2
    exit 1
fi
```

### R Scripts

```r
#!/usr/bin/env Rscript

# Header comment with script purpose
# Author: Your Name
# Date: YYYY-MM-DD

# Load required packages
library(package_name)

# Use meaningful variable names
variable_name <- value

# Comment complex operations
result <- operation(data,
                   param1 = value1,
                   param2 = value2)

# Use functions for repeated code
my_function <- function(param1, param2) {
    # Function body
    result <- computation(param1, param2)
    return(result)
}

# Error handling
tryCatch({
    result <- my_function(param1, param2)
}, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    quit(status = 1)
})
```

### General Guidelines

- **Naming:** Use descriptive, meaningful names
- **Comments:** Comment WHY, not WHAT (code should be clear)
- **Functions:** Break code into reusable functions
- **Error Handling:** Check for errors and provide helpful messages
- **Logging:** Add informative output during execution
- **Documentation:** Document inputs, outputs, and parameters

## Testing

### Before Submitting

Test with:
1. **Small dataset** - Verify basic functionality
2. **Full dataset** - Ensure scalability
3. **Different conditions** - Test various parameter combinations
4. **Edge cases** - Test boundary conditions

### Example Test Steps

```bash
# Set up test environment
bash setup.sh

# Download test data
fasterq-dump SRR11262284 --threads 4 --split-files -O $FASTQ

# Run pipeline
bash $SCRIPTS/fastqc_analysis.sh
bash $SCRIPTS/trimming.sh
bash $SCRIPTS/alignment.sh
bash $SCRIPTS/quantification.sh

# Verify outputs
ls -lh $RESULTS/
head $COUNTS/count_matrix.txt
```

## Documentation

### Documentation Standards

1. **Code Comments**
   - Explain complex logic
   - Document assumptions
   - Provide examples when helpful

2. **Function Documentation**
   ```bash
   # Function: process_data
   # Description: Processes input data and returns results
   # Arguments:
   #   $1 - Input file path
   #   $2 - Processing type
   # Returns:
   #   Output file path
   # Example:
   #   process_data "input.txt" "type1"
   ```

3. **README Updates**
   - Update if behavior changes
   - Add examples for new features
   - Document configuration options

4. **Changelog**
   - Document all significant changes
   - Note breaking changes
   - Provide upgrade instructions

## Code Review

All submissions will be reviewed for:

- ✓ Code quality and style
- ✓ Testing and validation
- ✓ Documentation completeness
- ✓ Backward compatibility
- ✓ Performance impact
- ✓ Security considerations

## Licensing

By contributing to this project, you agree that your contributions will be licensed under the same license as the project (MIT License).

## Questions?

- Check existing documentation
- Search previous issues
- Open a discussion issue
- Contact maintainers

---

**Thank you for contributing!** Your efforts help make this tool better for everyone.
