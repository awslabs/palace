# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Helpers for patching external dependencies
#

if(__patch_helpers)
  return()
endif()
set(__patch_helpers YES)

# Download a patch file to PATCH_DIR, verify its hash, and append its path to the list
# variable PATCH_FILES in the caller's scope
function(download_patch PATCH_FILES PATCH_DIR FILENAME URL SHA256)
  set(PATCH_FILE "${PATCH_DIR}/${FILENAME}")
  file(DOWNLOAD
    "${URL}"
    "${PATCH_FILE}"
    EXPECTED_HASH "SHA256=${SHA256}"
    TLS_VERIFY ON
  )
  list(APPEND ${PATCH_FILES} "${PATCH_FILE}")
  set(${PATCH_FILES} "${${PATCH_FILES}}" PARENT_SCOPE)
endfunction()
