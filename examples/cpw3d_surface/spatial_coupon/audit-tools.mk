# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

AUDIT_SOURCE := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PALACE_BUILD ?= $(abspath $(AUDIT_SOURCE)/../../../build)
AUDIT_OUTPUT ?= /tmp/palace-coupon-audit
include $(PALACE_BUILD)/share/mfem/config.mk

.PHONY: all
all: $(AUDIT_OUTPUT)/audit_coupon_mesh $(AUDIT_OUTPUT)/audit_surface_resolution $(AUDIT_OUTPUT)/audit_mesh_measures $(AUDIT_OUTPUT)/audit_interface_ownership

$(AUDIT_OUTPUT)/%: $(AUDIT_SOURCE)/%.cpp
	@mkdir -p $(AUDIT_OUTPUT)
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS) -I$(AUDIT_SOURCE)/../../../palace $< $(MFEM_LIBS) -o $@
