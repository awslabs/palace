// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_METADATA_FILE_HPP
#define PALACE_UTILS_METADATA_FILE_HPP

#include <ostream>
#include <nlohmann/json.hpp>

namespace palace
{

// Lane fix-response-memory round 2a R1 (decision 1475 (1)): the metadata file (palace.json)
// is streamed by nlohmann's serializer straight into the output stream (operator<< with
// width 2 = dump(2): the same serializer, indent character, ensure_ascii and error
// handler), followed by a newline, without building the whole text as a string first.
// Returns whether the stream accepted every byte (the caller fails loudly on a short write:
// a full disk must never leave a silently truncated or zero-filled palace.json).
inline bool StreamJsonMetadata(std::ostream &os, const nlohmann::json &meta)
{
  os.fill(' ');
  os.width(2);
  os << meta << '\n';
  os.flush();
  return static_cast<bool>(os);
}

}  // namespace palace

#endif  // PALACE_UTILS_METADATA_FILE_HPP
