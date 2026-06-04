# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

CLANG_FORMAT ?= clang-format
JULIA ?= julia

.PHONY: format format-cpp format-jl docs tests

# Style/format
format: format-cpp format-jl

format-cpp:
	./scripts/format-source -cpp --clang-format $(CLANG_FORMAT)

format-jl:
	./scripts/format-source -jl --julia $(JULIA)

# Documentation
docs:
	$(RM) -r docs/build
	$(JULIA) --project=docs -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=docs --color=yes docs/make.jl

docs-generate-config:
	$(JULIA) --project=docs --color=yes docs/generate_config_docs.jl

# Tests
#
# Runs the [Regression] cases against an existing CMake build tree.
# Override BUILD_DIR if your build lives elsewhere; PALACE_REGRESSION_NUMPROC
# (default 2, set in test/unit/CMakeLists.txt) controls per-case rank count.
BUILD_DIR ?= build
tests:
	cd $(BUILD_DIR)/palace-build/test/unit && ctest -L "^regression$$" --output-on-failure
