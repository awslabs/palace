// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_METADATA_FILE_HPP
#define PALACE_UTILS_METADATA_FILE_HPP

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>

namespace palace
{

// Lane fix-response-memory round 2a R2 (decision 1475 (1)): an array of palace.json held in
// a scratch file instead of memory. Each entry is appended as the exact text dump(2) gives
// it at its nesting depth in palace.json (the array at `keys` from the top level), entries
// joined by ",\n"; the metadata write splices the file in place of the (empty) in-memory
// array (StreamJsonMetadata), so the file bytes are those of the in-memory form. The file
// is created on the first entry; every write is checked (the callers fail loudly).
class JsonArraySpill
{
public:
  JsonArraySpill(std::string path, std::vector<std::string> keys)
    : path(std::move(path)), keys(std::move(keys))
  {
  }

  // An owner destroyed before the final write (an abort, a unit fixture) leaves no scratch
  // file behind.
  ~JsonArraySpill()
  {
    if (os.is_open())
    {
      os.close();
      std::remove(path.c_str());
    }
  }
  JsonArraySpill(const JsonArraySpill &) = delete;
  JsonArraySpill &operator=(const JsonArraySpill &) = delete;

  const std::string &Path() const { return path; }
  const std::vector<std::string> &Keys() const { return keys; }
  std::size_t Size() const { return count; }
  std::uint64_t Bytes() const { return bytes; }

  // The indent of the array's entries in palace.json (dump(2): two per nesting level).
  unsigned int EntryIndent() const
  {
    return 2 * static_cast<unsigned int>(keys.size() + 1);
  }

  // Appends one entry; false on a failed (short) write.
  bool Append(const nlohmann::json &entry)
  {
    if (!os.is_open())
    {
      os.open(path, std::ios::binary | std::ios::trunc);
      if (!os.is_open())
      {
        return false;
      }
    }
    if (count > 0)
    {
      os.write(",\n", 2);
    }
    const std::string indent(EntryIndent(), ' ');
    os.write(indent.data(), static_cast<std::streamsize>(indent.size()));
    nlohmann::detail::serializer<nlohmann::json> serializer(
        nlohmann::detail::output_adapter<char>(os), ' ');
    serializer.dump(entry, true, false, 2, EntryIndent());
    count++;
    return static_cast<bool>(os);
  }

  // Flushes the appended entries to the file; false on a failed write. Records the size.
  bool Flush()
  {
    if (!os.is_open())
    {
      return count == 0;
    }
    os.flush();
    if (!os)
    {
      return false;
    }
    bytes = static_cast<std::uint64_t>(os.tellp());
    return true;
  }

  // Closes the file (checked) and removes it; false on a failed close or removal.
  bool Remove()
  {
    bool ok = true;
    if (os.is_open())
    {
      ok = Flush();
      os.close();
      ok = ok && !os.fail();
      ok = std::remove(path.c_str()) == 0 && ok;
    }
    count = 0;
    bytes = 0;
    return ok;
  }

private:
  std::string path;
  std::vector<std::string> keys;
  std::ofstream os;
  std::size_t count = 0;
  std::uint64_t bytes = 0;
};

namespace metadata_file
{

inline void WriteIndent(std::ostream &os, unsigned int indent)
{
  for (unsigned int i = 0; i < indent; i++)
  {
    os.put(' ');
  }
}

// The spilled array at `indent` as dump(2) writes an array: "[]" when empty, else "[\n",
// the entries (already at their indent, joined by ",\n"), "\n", the indent and "]".
inline bool StreamSpilledArray(std::ostream &os, const JsonArraySpill &spill,
                               unsigned int indent)
{
  if (spill.Size() == 0)
  {
    os.write("[]", 2);
    return static_cast<bool>(os);
  }
  std::ifstream in(spill.Path(), std::ios::binary);
  if (!in.is_open())
  {
    return false;
  }
  os.write("[\n", 2);
  std::vector<char> buffer(1 << 20);
  std::uint64_t copied = 0;
  while (in)
  {
    in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto n = in.gcount();
    if (n > 0)
    {
      os.write(buffer.data(), n);
      copied += static_cast<std::uint64_t>(n);
    }
  }
  if (!in.eof() || copied != spill.Bytes())
  {
    return false;
  }
  os.put('\n');
  WriteIndent(os, indent);
  os.put(']');
  return static_cast<bool>(os);
}

// `value` at `indent` as dump(2) writes it, the array at spill.Keys()[depth..] replaced by
// the spilled one: the objects on the path written member by member in their (sorted)
// order, every other member by nlohmann's own serializer at its indent.
inline bool StreamSpliced(std::ostream &os, const nlohmann::json &value,
                          const JsonArraySpill &spill, std::size_t depth,
                          unsigned int indent)
{
  const auto &keys = spill.Keys();
  if (depth == keys.size())
  {
    return StreamSpilledArray(os, spill, indent);
  }
  os.write("{\n", 2);
  nlohmann::detail::serializer<nlohmann::json> serializer(
      nlohmann::detail::output_adapter<char>(os), ' ');
  std::size_t i = 0;
  for (auto it = value.cbegin(); it != value.cend(); ++it, ++i)
  {
    WriteIndent(os, indent + 2);
    const std::string key = nlohmann::json(it.key()).dump();
    os.write(key.data(), static_cast<std::streamsize>(key.size()));
    os.write(": ", 2);
    if (it.key() == keys[depth])
    {
      if (!StreamSpliced(os, *it, spill, depth + 1, indent + 2))
      {
        return false;
      }
    }
    else
    {
      serializer.dump(*it, true, false, 2, indent + 2);
    }
    if (i + 1 < value.size())
    {
      os.write(",\n", 2);
    }
  }
  os.put('\n');
  WriteIndent(os, indent);
  os.put('}');
  return static_cast<bool>(os);
}

}  // namespace metadata_file

// Whether `meta` holds the spilled array's place: an EMPTY array at spill.Keys() (the
// in-memory placeholder; the entries live in the file).
inline bool HoldsSpillPlace(const nlohmann::json &meta, const JsonArraySpill &spill)
{
  const nlohmann::json *node = &meta;
  for (const auto &key : spill.Keys())
  {
    if (!node->is_object() || !node->contains(key))
    {
      return false;
    }
    node = &(*node)[key];
  }
  return node->is_array() && node->empty();
}

// Lane fix-response-memory round 2a R1 (decision 1475 (1)): the metadata file (palace.json)
// is streamed by nlohmann's serializer straight into the output stream (operator<< with
// width 2 = dump(2): the same serializer, indent character, ensure_ascii and error
// handler), followed by a newline, without building the whole text as a string first. R2:
// with a spill whose place `meta` holds (HoldsSpillPlace), the spilled entries are spliced
// in at that place; the bytes are those of `meta` holding the entries in memory. Returns
// whether the stream accepted every byte (and the spill file was read whole); the caller
// fails loudly on a short write: a full disk must never leave a silently incomplete
// palace.json.
inline bool StreamJsonMetadata(std::ostream &os, const nlohmann::json &meta,
                               const JsonArraySpill *spill = nullptr)
{
  if (spill && HoldsSpillPlace(meta, *spill))
  {
    if (!metadata_file::StreamSpliced(os, meta, *spill, 0, 0))
    {
      return false;
    }
    os.put('\n');
  }
  else
  {
    os.fill(' ');
    os.width(2);
    os << meta << '\n';
  }
  os.flush();
  return static_cast<bool>(os);
}

}  // namespace palace

#endif  // PALACE_UTILS_METADATA_FILE_HPP
