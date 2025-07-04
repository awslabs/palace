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

# Tests
tests:
	$(JULIA) --project=test/examples -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/examples --color=yes test/examples/runtests.jl

build: ## Build python package
	rm -rf dist
	pip install build
	python -m build
