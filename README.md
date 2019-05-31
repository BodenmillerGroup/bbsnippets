[![Build Status](https://travis-ci.org/BodenmillerGroup/bbsnippets.svg?branch=master)](https://travis-ci.org/BodenmillerGroup/bbsnippets)

## Getting started

All snippets are separated by domain-specific folders. Each such folder can contain multiple subfolders for different
tasks/topics. It is recommended to split examples in different languages (R, Python, etc.) into a separate directory,
so programming tools can manage these snippets accordingly.

## Code conventions

We suggest authors of submitted examples to follow commonly accepted source code styling practices. There are
.editorconfig and .prettierrc config files in the root folder that can be used to apply such rules. For code formatting
these tools provide a strongly opinionated styling rules with good editor integration:

- Python: Black (https://black.readthedocs.io)
- R: styler (https://styler.r-lib.org/)
- Other languages: Prettier (https://prettier.io/)
