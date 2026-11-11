// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmshreader.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>
#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#define PALACE_GMSH_MMAP
#endif

namespace palace::mesh
{

namespace
{

// Read-only view of a whole file, memory-mapped when possible.
class FileView
{
private:
  const char *data = nullptr;
  std::size_t size = 0;
  void *map = nullptr;
  std::vector<char> buffer;

public:
  explicit FileView(const std::string &filename)
  {
#if defined(PALACE_GMSH_MMAP)
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0)
    {
      return;
    }
    struct stat st;
    if (fstat(fd, &st) == 0 && st.st_size > 0)
    {
      size = static_cast<std::size_t>(st.st_size);
      map = mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
      if (map == MAP_FAILED)
      {
        map = nullptr;
        buffer.resize(size);
        std::size_t off = 0;
        while (off < size)
        {
          ssize_t n = read(fd, buffer.data() + off, size - off);
          if (n <= 0)
          {
            break;
          }
          off += static_cast<std::size_t>(n);
        }
        if (off == size)
        {
          data = buffer.data();
        }
        else
        {
          size = 0;
        }
      }
      else
      {
        madvise(map, size, MADV_SEQUENTIAL);
        data = static_cast<const char *>(map);
      }
    }
    close(fd);
#else
    std::FILE *f = std::fopen(filename.c_str(), "rb");
    if (!f)
    {
      return;
    }
    std::fseek(f, 0, SEEK_END);
    long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    if (n > 0)
    {
      buffer.resize(static_cast<std::size_t>(n));
      if (std::fread(buffer.data(), 1, buffer.size(), f) == buffer.size())
      {
        data = buffer.data();
        size = buffer.size();
      }
    }
    std::fclose(f);
#endif
  }
  ~FileView()
  {
#if defined(PALACE_GMSH_MMAP)
    if (map)
    {
      munmap(map, size);
    }
#endif
  }
  FileView(const FileView &) = delete;
  FileView &operator=(const FileView &) = delete;
  const char *begin() const { return data; }
  const char *end() const { return data + size; }
  bool good() const { return data != nullptr && size > 0; }
};

inline bool IsSpace(char c)
{
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

inline bool IsDigit(char c)
{
  return c >= '0' && c <= '9';
}

// Whitespace-separated token scanner with the semantics of operator>> for well-formed
// numbers. Any token which operator>> would read differently (trailing characters,
// overflow) fails the scan, and the caller falls back to MFEM's reader.
class Scanner
{
private:
  const char *p, *end;

public:
  Scanner(const char *begin, const char *end) : p(begin), end(end) {}

  const char *Pos() const { return p; }

  void SkipSpace()
  {
    while (p < end && IsSpace(*p))
    {
      p++;
    }
  }

  bool AtTokenEnd() const { return p == end || IsSpace(*p); }

  bool Int(long long &v)
  {
    SkipSpace();
    bool neg = false;
    if (p < end && (*p == '-' || *p == '+'))
    {
      neg = (*p == '-');
      p++;
    }
    if (p == end || !IsDigit(*p))
    {
      return false;
    }
    long long x = 0;
    int digits = 0;
    while (p < end && IsDigit(*p))
    {
      x = 10 * x + (*p - '0');
      p++;
      if (++digits > 18)
      {
        return false;
      }
    }
    v = neg ? -x : x;
    return AtTokenEnd();
  }

  bool Int(int &v)
  {
    long long x;
    if (!Int(x) || x < std::numeric_limits<int>::min() ||
        x > std::numeric_limits<int>::max())
    {
      return false;
    }
    v = static_cast<int>(x);
    return true;
  }

  bool Double(double &v)
  {
    SkipSpace();
    if (p < end && *p == '+')
    {
      p++;
    }
    if (p == end || (*p != '-' && *p != '.' && !IsDigit(*p)))
    {
      return false;  // Also excludes inf/nan spellings
    }
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
    auto [ptr, ec] = std::from_chars(p, end, v);
    if (ec != std::errc())
    {
      return false;
    }
    p = ptr;
#else
    // strtod needs a terminated token: copy it.
    char tok[128];
    std::size_t n = 0;
    while (p + n < end && !IsSpace(p[n]) && n < sizeof(tok) - 1)
    {
      tok[n] = p[n];
      n++;
    }
    tok[n] = '\0';
    char *tend;
    v = std::strtod(tok, &tend);
    if (tend != tok + n || n == 0)
    {
      return false;
    }
    p += n;
#endif
    return AtTokenEnd();
  }

  // Native-endian binary value (bin_io::read of MFEM).
  template <typename T>
  bool Binary(T &v)
  {
    if (static_cast<std::size_t>(end - p) < sizeof(T))
    {
      return false;
    }
    std::memcpy(&v, p, sizeof(T));
    p += sizeof(T);
    return true;
  }

  // ChompNewline of the Gmsh reader: an optional '\r', then a required '\n'.
  bool Newline()
  {
    if (p < end && *p == '\r')
    {
      p++;
    }
    if (p == end || *p != '\n')
    {
      return false;
    }
    p++;
    return true;
  }

  // Next whitespace-separated token (as operator>> into a std::string).
  std::string_view Token()
  {
    SkipSpace();
    const char *b = p;
    while (p < end && !IsSpace(*p))
    {
      p++;
    }
    return std::string_view(b, p - b);
  }

  // Equivalent of the Gmsh reader's GoToNextSection: read lines (the first one is the rest
  // of the current line) until one starts with '$' but not with "$End", and return the
  // section name (the rest of that line, without a trailing '\r'), or "" at the end.
  std::string NextSection()
  {
    while (p < end)
    {
      const char *b = p;
      while (p < end && *p != '\n')
      {
        p++;
      }
      const char *e = p;
      if (p < end)
      {
        p++;  // Consume the newline
      }
      if (e > b && e[-1] == '\r')
      {
        e--;
      }
      std::string_view line(b, e - b);
      if (!line.empty() && line[0] == '$' &&
          (line.size() < 4 || line.compare(1, 3, "End") != 0))
      {
        return std::string(line.substr(1));
      }
    }
    return "";
  }

  // Equivalent of the Gmsh reader's ReadQuotedString.
  bool QuotedString(std::string &s)
  {
    while (p < end && *p != '"')
    {
      p++;
    }
    if (p == end)
    {
      return false;
    }
    p++;
    const char *b = p;
    while (p < end && *p != '"')
    {
      p++;
    }
    if (p == end)
    {
      return false;
    }
    s.assign(b, p - b);
    p++;
    return true;
  }
};

// First-order Gmsh element types: geometry and number of nodes.
bool GetGeometry(int el_type, mfem::Geometry::Type &geom, int &nv)
{
  switch (el_type)
  {
    case 15:
      geom = mfem::Geometry::POINT;
      nv = 1;
      return true;
    case 1:
      geom = mfem::Geometry::SEGMENT;
      nv = 2;
      return true;
    case 2:
      geom = mfem::Geometry::TRIANGLE;
      nv = 3;
      return true;
    case 3:
      geom = mfem::Geometry::SQUARE;
      nv = 4;
      return true;
    case 4:
      geom = mfem::Geometry::TETRAHEDRON;
      nv = 4;
      return true;
    case 5:
      geom = mfem::Geometry::CUBE;
      nv = 8;
      return true;
    case 6:
      geom = mfem::Geometry::PRISM;
      nv = 6;
      return true;
    case 7:
      geom = mfem::Geometry::PYRAMID;
      nv = 5;
      return true;
  }
  return false;  // High-order or unsupported
}

// Mesh built by the reader: gives access to the mesh type flags, which a mesh gets from
// FinalizeTopology otherwise.
class GmshMesh : public mfem::Mesh
{
public:
  GmshMesh(int dim, int nv, int ne, int nbe, int sdim) : mfem::Mesh(dim, nv, ne, nbe, sdim)
  {
  }
  void SetMeshGenAndAttributes()
  {
    SetMeshGen();
    SetAttributes();
  }
};

// Elements of one dimension, in file order.
struct ElementList
{
  std::vector<unsigned char> geom;
  std::vector<int> attr;
  std::vector<std::int64_t> offset = {0};
  std::vector<int> verts;
};

}  // namespace

std::unique_ptr<mfem::Mesh> ReadGmshMeshFast(const std::string &filename,
                                             bool finalize_topology)
{
  FileView file(filename);
  if (!file.good() || file.end() - file.begin() < 16 ||
      std::strncmp(file.begin(), "$MeshFormat", 11) != 0)
  {
    return nullptr;  // Not an uncompressed Gmsh file
  }
  Scanner in(file.begin(), file.end());

  // Header: same as mfem::Mesh::Loader (first line) and the GmshReader constructor.
  if (in.Token() != "$MeshFormat")
  {
    return nullptr;
  }
  int is_binary, data_size;
  if (in.Token() != "2.2" || !in.Int(is_binary) || (is_binary != 0 && is_binary != 1) ||
      !in.Int(data_size) || data_size != static_cast<int>(sizeof(double)))
  {
    return nullptr;  // Format 4.1: MFEM's reader
  }
  const bool binary = (is_binary == 1);
  if (binary)
  {
    int one;
    if (!in.Newline() || !in.Binary(one) || one != 1)
    {
      return nullptr;  // MFEM's reader reports the error
    }
  }
  // Integer and floating-point values of the node and element data.
  auto Int = [&in, binary](int &v) { return binary ? in.Binary(v) : in.Int(v); };
  auto Double = [&in, binary](double &v) { return binary ? in.Binary(v) : in.Double(v); };

  // Sections (ReadGmsh2Mesh).
  int num_vertices = -1;
  std::vector<double> coords;
  std::vector<int> vertex_map;  // Node tag to vertex index (-1 if none)
  double bb_min[3] = {mfem::infinity(), mfem::infinity(), mfem::infinity()};
  double bb_max[3] = {-mfem::infinity(), -mfem::infinity(), -mfem::infinity()};
  int space_dim = 0;
  ElementList elems[4];
  bool has_elements = false, has_positive_attrs = false, has_non_positive_attrs = false;
  std::unordered_map<int, std::unordered_map<int, std::string>> phys_names_by_dim;
  std::string section;
  do
  {
    section = in.NextSection();
    if (section == "Nodes")
    {
      if (!in.Int(num_vertices) || num_vertices < 0 || (binary && !in.Newline()))
      {
        return nullptr;
      }
      coords.resize(3 * static_cast<std::size_t>(num_vertices));
      std::vector<int> tags(num_vertices);
      int max_tag = 0;
      for (int v = 0; v < num_vertices; v++)
      {
        if (!Int(tags[v]) || tags[v] < 0)
        {
          return nullptr;
        }
        max_tag = std::max(max_tag, tags[v]);
        for (int d = 0; d < 3; d++)
        {
          double &c = coords[3 * static_cast<std::size_t>(v) + d];
          if (!Double(c))
          {
            return nullptr;
          }
          bb_min[d] = std::min(bb_min[d], c);
          bb_max[d] = std::max(bb_max[d], c);
        }
      }
      // Dense tag map (node tags are usually 1..N); MFEM's reader requires unique tags.
      if (static_cast<double>(max_tag) > 4.0 * num_vertices + 1024.0)
      {
        return nullptr;
      }
      vertex_map.assign(static_cast<std::size_t>(max_tag) + 1, -1);
      for (int v = 0; v < num_vertices; v++)
      {
        if (vertex_map[tags[v]] >= 0)
        {
          return nullptr;  // Non-unique node tags: MFEM's reader reports the error
        }
        vertex_map[tags[v]] = v;
      }
      // GetSpaceDimension of the Gmsh reader.
      constexpr double bb_tol = 1e-14;
      const double bb_size = std::max(
          bb_max[0] - bb_min[0], std::max(bb_max[1] - bb_min[1], bb_max[2] - bb_min[2]));
      space_dim = 1;
      if (bb_max[1] - bb_min[1] > bb_size * bb_tol)
      {
        space_dim++;
      }
      if (bb_max[2] - bb_min[2] > bb_size * bb_tol)
      {
        space_dim++;
      }
    }
    else if (section == "Elements")
    {
      int num_elements;
      if (num_vertices < 0 || !in.Int(num_elements) || num_elements < 0 ||
          (binary && !in.Newline()))
      {
        return nullptr;  // Elements before nodes: MFEM's reader
      }
      has_elements = true;
      // Add an element given its type, tags and nodes (NewElement/SetAttribute of the Gmsh
      // reader).
      auto AddElement = [&](int el_type, int n_tags, auto &&read_int) -> bool
      {
        int el_phys_tag = 0;
        for (int t = 0; t < n_tags; t++)
        {
          int tag;
          if (!read_int(tag))
          {
            return false;
          }
          if (t == 0)
          {
            el_phys_tag = tag;
          }
        }
        mfem::Geometry::Type geom;
        int nv;
        if (!GetGeometry(el_type, geom, nv))
        {
          return false;
        }
        auto &list = elems[mfem::Geometry::Dimension[geom]];
        list.geom.push_back(static_cast<unsigned char>(geom));
        if (el_phys_tag < 1)
        {
          has_non_positive_attrs = true;
          el_phys_tag = 1;
        }
        else
        {
          has_positive_attrs = true;
        }
        list.attr.push_back(el_phys_tag);
        for (int j = 0; j < nv; j++)
        {
          int tag;
          if (!read_int(tag) || tag < 0 || tag >= static_cast<int>(vertex_map.size()) ||
              vertex_map[tag] < 0)
          {
            return false;  // Unknown node tag
          }
          list.verts.push_back(vertex_map[tag]);
        }
        list.offset.push_back(static_cast<std::int64_t>(list.verts.size()));
        return true;
      };
      auto read_bin = [&in](int &v) { return in.Binary(v); };
      auto read_ascii = [&in](int &v) { return in.Int(v); };
      int num_el_read = 0;
      while (num_el_read < num_elements)
      {
        if (binary)
        {
          // Block header, then the elements of the block.
          int el_type, n_els, n_tags;
          if (!in.Binary(el_type) || !in.Binary(n_els) || !in.Binary(n_tags) || n_els < 0 ||
              n_tags < 0)
          {
            return nullptr;
          }
          for (int e = 0; e < n_els; e++)
          {
            int el_num;
            if (!in.Binary(el_num) || !AddElement(el_type, n_tags, read_bin))
            {
              return nullptr;
            }
            num_el_read++;
          }
          if (n_els == 0)
          {
            return nullptr;  // Would not terminate
          }
        }
        else
        {
          int el_num, el_type, n_tags;
          if (!in.Int(el_num) || !in.Int(el_type) || !in.Int(n_tags) || n_tags < 0 ||
              !AddElement(el_type, n_tags, read_ascii))
          {
            return nullptr;
          }
          num_el_read++;
        }
      }
    }
    else if (section == "PhysicalNames")
    {
      int num_names;
      if (!in.Int(num_names))
      {
        return nullptr;
      }
      for (int i = 0; i < num_names; i++)
      {
        int phys_dim, phys_tag;
        std::string name;
        if (!in.Int(phys_dim) || !in.Int(phys_tag) || !in.QuotedString(name))
        {
          return nullptr;
        }
        phys_names_by_dim[phys_dim][phys_tag] = name;
      }
    }
    else if (section == "Periodic")
    {
      return nullptr;  // Periodic meshes: MFEM's reader
    }
  } while (!section.empty());
  if (num_vertices < 0 || !has_elements)
  {
    return nullptr;
  }

  // CheckAttributes of the Gmsh reader.
  if (has_non_positive_attrs)
  {
    MFEM_VERIFY(!has_positive_attrs, "Non-positive element attribute in Gmsh mesh!\n"
                                     "By default Gmsh sets element tags (attributes)"
                                     " to '0' but MFEM requires that they be"
                                     " positive integers.\n"
                                     "Use \"Physical Curve\", \"Physical Surface\","
                                     " or \"Physical Volume\" to set tags/attributes"
                                     " for all curves, surfaces, or volumes in your"
                                     " Gmsh geometry to values which are >= 1.");
    MFEM_WARNING("Gmsh reader: all element attributes were zero.\n"
                 "MFEM only supports positive element attributes.\n"
                 "Setting all element attributes to 1.\n");
  }

  // AddElements of the Gmsh reader: elements of the highest dimension and boundary
  // elements of one dimension less.
  const int dim = !elems[3].geom.empty() ? 3 : (!elems[2].geom.empty() ? 2 : 1);
  auto &dom = elems[dim], &bdr = elems[dim - 1];

  // Vertices not used by the elements or boundary elements are removed as by
  // mfem::Mesh::RemoveUnusedVertices (keeping the order of the others), here before the
  // mesh is built since that function also regenerates the topology.
  std::vector<int> new_index(num_vertices, -1);
  for (const auto *list : {&dom, &bdr})
  {
    for (int v : list->verts)
    {
      new_index[v] = 0;
    }
  }
  int num_used = 0;
  for (int v = 0; v < num_vertices; v++)
  {
    if (new_index[v] == 0)
    {
      new_index[v] = num_used++;
    }
  }
  if (num_used < num_vertices)
  {
    for (auto *list : {&dom, &bdr})
    {
      for (int &v : list->verts)
      {
        v = new_index[v];
      }
    }
  }
  auto mesh = std::make_unique<GmshMesh>(dim, num_used, static_cast<int>(dom.geom.size()),
                                         static_cast<int>(bdr.geom.size()), space_dim);
  for (int v = 0; v < num_vertices; v++)
  {
    if (new_index[v] >= 0)
    {
      const double *x = coords.data() + 3 * static_cast<std::size_t>(v);
      mesh->AddVertex(x[0], x[1], x[2]);
    }
  }
  new_index = std::vector<int>();
  coords = std::vector<double>();
  vertex_map = std::vector<int>();
  auto AddElements = [&mesh](const ElementList &list, bool boundary)
  {
    for (std::size_t i = 0; i < list.geom.size(); i++)
    {
      mfem::Element *el = mesh->NewElement(list.geom[i]);
      el->SetVertices(list.verts.data() + list.offset[i]);
      el->SetAttribute(list.attr[i]);
      if (boundary)
      {
        mesh->AddBdrElement(el);
      }
      else
      {
        mesh->AddElement(el);
      }
    }
  };
  AddElements(dom, false);
  AddElements(bdr, true);
  for (auto &list : elems)
  {
    list = ElementList();
  }

  // AddPhysicalNames of the Gmsh reader.
  for (const auto &[attr, name] : phys_names_by_dim[dim - 1])
  {
    if (!mesh->bdr_attribute_sets.AttributeSetExists(name))
    {
      mesh->bdr_attribute_sets.CreateAttributeSet(name);
    }
    mesh->bdr_attribute_sets.AddToAttributeSet(name, attr);
  }
  for (const auto &[attr, name] : phys_names_by_dim[dim])
  {
    if (!mesh->attribute_sets.AttributeSetExists(name))
    {
      mesh->attribute_sets.CreateAttributeSet(name);
    }
    mesh->attribute_sets.AddToAttributeSet(name, attr);
  }
  mesh->RemoveUnusedVertices();  // Nothing to remove (see above)
  // The topology can be left to the finalization of the mesh when it updates the topology
  // anyway (meshes with simplices, and 2D meshes): see mfem::Mesh::Finalize.
  mesh->SetMeshGenAndAttributes();
  if (finalize_topology || !(dim == 2 || (dim == 3 && (mesh->MeshGenerator() & 1))))
  {
    mesh->FinalizeTopology();
  }
  return std::make_unique<mfem::Mesh>(std::move(*mesh));
}

}  // namespace palace::mesh
