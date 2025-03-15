// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "meshio.hpp"

#include <unordered_map>
#include <vector>
#include <mfem.hpp>

#define GMSH_BIN  // Use binary Gmsh format

namespace palace
{

namespace
{

inline int ElemTypeComsol(const std::string &type)
{
  if (!type.compare(0, 3, "tri"))  // 3-node triangle
  {
    return 2;
  }
  if (!type.compare(0, 4, "quad"))  // 4-node quadrangle
  {
    return 3;
  }
  if (!type.compare(0, 3, "tet"))  // 4-node tetrahedron
  {
    return 4;
  }
  if (!type.compare(0, 3, "hex"))  // 8-node hexahedron
  {
    return 5;
  }
  if (!type.compare(0, 5, "prism"))  // 6-node prism
  {
    return 6;
  }
  if (!type.compare(0, 3, "pyr"))  // 5-node pyramid
  {
    return 7;
  }
  if (!type.compare(0, 4, "tri2"))  // 6-node triangle
  {
    return 9;
  }
  if (!type.compare(0, 5, "quad2"))  // 9-node quadrangle
  {
    return 10;
  }
  if (!type.compare(0, 4, "tet2"))  // 10-node tetrahedron
  {
    return 11;
  }
  if (!type.compare(0, 4, "hex2"))  // 27-node hexahedron
  {
    return 12;
  }
  if (!type.compare(0, 6, "prism2"))  // 18-node prism
  {
    return 13;
  }
  if (!type.compare(0, 4, "pyr2"))  // 14-node pyramid
  {
    return 14;
  }
  return 0;  // Skip this element type
}

inline int ElemTypeNastran(const std::string &type)
{
  // Returns only the low-order type for a given keyword.
  if (!type.compare(0, 5, "CTRIA"))  // 3-node triangle
  {
    return 2;
  }
  if (!type.compare(0, 5, "CQUAD"))  // 4-node quadrangle
  {
    return 3;
  }
  if (!type.compare(0, 6, "CTETRA"))  // 4-node tetrahedron
  {
    return 4;
  }
  if (!type.compare(0, 5, "CHEXA"))  // 8-node hexahedron
  {
    return 5;
  }
  if (!type.compare(0, 6, "CPENTA"))  // 6-node prism
  {
    return 6;
  }
  if (!type.compare(0, 6, "CPYRAM"))  // 5-node pyramid
  {
    return 7;
  }
  return 0;  // Skip this element type
}

inline int HOElemTypeNastran(const int lo_type, const int num_nodes)
{
  // Get high-order element type for corresponding low-order type.
  if (lo_type == 2 && num_nodes > 3)  // 6-node triangle
  {
    MFEM_VERIFY(num_nodes == 6, "Invalid high-order Nastran element!");
    return 9;
  }
  if (lo_type == 3)
  {
    if (num_nodes == 9)  // 9-node quadrangle
    {
      return 10;
    }
    if (num_nodes == 8)  // 8-node quadrangle
    {
      return 16;
    }
    MFEM_VERIFY(num_nodes == 4, "Invalid high-order Nastran element!");
    return lo_type;
  }
  if (lo_type == 4 && num_nodes > 4)  // 10-node tetrahedron
  {
    MFEM_VERIFY(num_nodes == 10, "Invalid high-order Nastran element!");
    return 11;
  }
  if (lo_type == 5 && num_nodes > 8)  // 20-node hexahedron
  {
    MFEM_VERIFY(num_nodes == 20, "Invalid high-order Nastran element!");
    return 17;
  }
  if (lo_type == 6 && num_nodes > 6)  // 15-node prism
  {
    MFEM_VERIFY(num_nodes == 15, "Invalid high-order Nastran element!");
    return 18;
  }
  if (lo_type == 7 && num_nodes > 5)  // 13-node pyramid
  {
    MFEM_VERIFY(num_nodes == 13, "Invalid high-order Nastran element!");
    return 19;
  }
  return lo_type;
}

inline int LOElemTypeGmsh(int ho_type)
{
  if (ho_type == 9)  // 6-node triangle
  {
    return 2;
  }
  if (ho_type == 10 || ho_type == 16)  // 9- or 8-node quadrangle
  {
    return 3;
  }
  if (ho_type == 11)  // 10-node tetrahedron
  {
    return 4;
  }
  if (ho_type == 12 || ho_type == 17)  // 27- or 20-node hexahedron
  {
    return 5;
  }
  if (ho_type == 13 || ho_type == 18)  // 18- or 15-node prism
  {
    return 6;
  }
  if (ho_type == 14 || ho_type == 19)  // 14- or 13-node pyramid
  {
    return 7;
  }
  return ho_type;
}

constexpr int ElemNumNodes[] = {-1,  // 2-node edge
                                3,  4,  4,  8,  6,  5,
                                -1,  // 3-node edge
                                6,  9,  10, 27, 18, 14,
                                -1,  // 1-node node
                                8,  20, 15, 13};

// From COMSOL or Nastran to Gmsh ordering. See:
//   - https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
//   - https://tinyurl.com/yezswzfv
//   - https://tinyurl.com/4d32zxtn
constexpr int SkipElem[] = {-1};
constexpr int Msh3[] = {0, 1, 2};
constexpr int Msh4[] = {0, 1, 2, 3};
constexpr int Msh5[] = {0, 1, 2, 3, 4};
constexpr int Msh6[] = {0, 1, 2, 3, 4, 5};
constexpr int Msh8[] = {0, 1, 2, 3, 4, 5, 6, 7};
constexpr int Msh9[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

constexpr int MphQuad4[] = {0, 1, 3, 2};
constexpr int MphHex8[] = {0, 1, 3, 2, 4, 5, 7, 6};
constexpr int MphPyr5[] = {0, 1, 3, 2, 4};
constexpr int MphTri6[] = {0, 1, 2, 3, 5, 4};
constexpr int MphQuad9[] = {0, 1, 3, 2, 4, 7, 8, 5, 6};
constexpr int MphTet10[] = {0, 1, 2, 3, 4, 6, 5, 7, 9, 8};
constexpr int MphHex27[] = {0,  1,  3,  2,  4,  5,  7,  6,  8,  9,  20, 11, 13, 10,
                            21, 12, 22, 26, 23, 15, 24, 14, 16, 17, 25, 18, 19};
constexpr int MphWdg18[] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 15, 10, 16, 17, 11, 12, 13, 14};
constexpr int MphPyr14[] = {0, 1, 3, 2, 4, 5, 6, 13, 8, 10, 7, 9, 12, 11};

constexpr int NasTet10[] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8};
constexpr int NasHex20[] = {0,  1, 2,  3,  4,  5,  6,  7,  8,  11,
                            13, 9, 10, 12, 14, 15, 16, 18, 19, 17};
constexpr int NasWdg15[] = {0, 1, 2, 3, 4, 5, 6, 9, 7, 8, 10, 11, 12, 14, 13};
constexpr int NasPyr13[] = {0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12};

constexpr const int *ElemNodesComsol[] = {SkipElem, Msh3,     MphQuad4, Msh4,     MphHex8,
                                          Msh6,     MphPyr5,  SkipElem, MphTri6,  MphQuad9,
                                          MphTet10, MphHex27, MphWdg18, MphPyr14, SkipElem,
                                          SkipElem, SkipElem, SkipElem, SkipElem};
constexpr const int *ElemNodesNastran[] = {SkipElem, Msh3,     Msh4,     Msh4,     Msh8,
                                           Msh6,     Msh5,     SkipElem, Msh6,     Msh9,
                                           NasTet10, SkipElem, SkipElem, SkipElem, SkipElem,
                                           Msh8,     NasHex20, NasWdg15, NasPyr13};

// Get line, strip comments, leading/trailing whitespace. Should not be called if end of
// file is expected.
inline std::string GetLineComsol(std::ifstream &input)
{
  std::string str;
  std::getline(input, str);
  MFEM_VERIFY(input, "Unexpected read failure parsing mesh file!");
  const auto pos = str.find_first_of('#');
  if (pos != std::string::npos)
  {
    str.erase(pos);
  }
  const auto start = str.find_first_not_of(" \t");
  if (start == std::string::npos)
  {
    return "";
  }
  const auto stop = str.find_last_not_of(" \t");
  return str.substr(start, stop - start + 1);
}

inline std::string GetLineNastran(std::ifstream &input)
{
  std::string str;
  std::getline(input, str);
  MFEM_VERIFY(input.good(), "Unexpected read failure parsing mesh file!");
  str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
  return str[0] == '$' ? "" : str;
}

// COMSOL strings are parsed as an integer length followed by array of integers for the
// string characters.
inline std::string ReadStringComsol(std::istream &input)
{
  int n;
  std::string str;
  input >> n >> str;
  return str;
}

inline std::string ReadStringComsolBinary(std::istream &input)
{
  int n;
  input.read(reinterpret_cast<char *>(&n), sizeof(int));
  std::vector<int> vstr(n);
  input.read(reinterpret_cast<char *>(vstr.data()), (std::streamsize)(n * sizeof(int)));
  return std::string(vstr.begin(), vstr.end());
}

// Nastran has a special floating point format: "-7.-1" instead of "-7.E-01" or "2.3+2"
// instead of "2.3E+02".
inline double ConvertDoubleNastran(const std::string &str)
{
  double d;
  try
  {
    d = std::stod(str);
  }
  catch (const std::invalid_argument &ia)
  {
    const std::size_t start = str.find_first_not_of(' ');
    MFEM_VERIFY(start != std::string::npos,
                "Invalid number conversion parsing Nastran mesh!")
    std::string fstr = str.substr(start);
    std::size_t pos = fstr.find('+', 1);  // Skip leading +/- sign
    if (pos != std::string::npos)
    {
      fstr.replace(pos, 1, "E+");
    }
    else if ((pos = fstr.find('-', 1)) != std::string::npos)
    {
      fstr.replace(pos, 1, "E-");
    }
    d = std::stod(fstr);
  }
  return d;
}

inline void WriteNode(std::ostream &buffer, const int tag, const double *coord)
{
#if defined(GMSH_BIN)
  buffer.write(reinterpret_cast<const char *>(&tag), sizeof(int));
  buffer.write(reinterpret_cast<const char *>(coord), 3 * sizeof(double));
  // No newline for binary data
#else
  // Always 3D coordinates (user sets floating point format/precision on buffer).
  buffer << tag << ' ' << coord[0] << ' ' << coord[1] << ' ' << coord[2] << '\n';
#endif
}

inline void WriteElement(std::ostream &buffer, const int tag, const int type,
                         const int geom, const int nodes[])
{
#if defined(GMSH_BIN)
  const int data[3] = {tag, geom, geom};
  buffer.write(reinterpret_cast<const char *>(data), 3 * sizeof(int));
  buffer.write(reinterpret_cast<const char *>(nodes),
               (std::streamsize)(ElemNumNodes[type - 1] * sizeof(int)));
  // No newline for binary data
#else
  buffer << tag << ' ' << type << " 2 " << geom << ' ' << geom;
  for (int i = 0; i < ElemNumNodes[type - 1]; i++)
  {
    buffer << ' ' << nodes[i];
  }
  buffer << '\n';
#endif
}

void WriteGmsh(std::ostream &buffer, const std::vector<double> &node_coords,
               const std::vector<int> &node_tags,
               const std::unordered_map<int, std::vector<int>> &elem_nodes,
               const bool use_lo_type)
{
  // Write the Gmsh file header (version 2.2).
  buffer << "$MeshFormat\n2.2 "
         <<
#if defined(GMSH_BIN)
      "1 " <<
#else
      "0 " <<
#endif
      sizeof(double) << '\n';
#if defined(GMSH_BIN)
  const int one = 1;
  buffer.write(reinterpret_cast<const char *>(&one), sizeof(int));
  buffer << '\n';
#endif
  buffer << "$EndMeshFormat\n";

  // Write mesh nodes.
  const int num_nodes = (int)node_coords.size() / 3;
  MFEM_VERIFY(num_nodes > 0 && node_coords.size() % 3 == 0,
              "Gmsh nodes should always be in 3D space!");
  buffer << "$Nodes\n" << num_nodes << '\n';
  {
    if (!node_tags.empty())
    {
      // Use input node tags which should be positive but don't need to be contiguous.
      MFEM_VERIFY(node_tags.size() == (std::size_t)num_nodes,
                  "Invalid size for node tags!");
      for (int i = 0; i < num_nodes; i++)
      {
        WriteNode(buffer, node_tags[i], &node_coords[3 * i]);
      }
    }
    else
    {
      // Label nodes as contiguous starting at 1.
      for (int i = 0; i < num_nodes; i++)
      {
        WriteNode(buffer, i + 1, &node_coords[3 * i]);
      }
    }
  }
#if defined(GMSH_BIN)
  buffer << '\n';
#endif
  buffer << "$EndNodes\n";

  // Write mesh elements.
  int tot_num_elem = 0;
  for (const auto &[elem_type, nodes] : elem_nodes)
  {
    MFEM_VERIFY(elem_type > 0, "Invalid element type writing Gmsh elements!");
    const int &num_elem_nodes = ElemNumNodes[elem_type - 1];
    tot_num_elem += ((int)nodes.size()) / (num_elem_nodes + 1);
    MFEM_VERIFY(nodes.size() % (num_elem_nodes + 1) == 0,
                "Unexpected data size when writing elements!");
  }
  MFEM_VERIFY(tot_num_elem > 0, "No mesh elements parsed from COMSOL mesh file!");
  buffer << "$Elements\n" << tot_num_elem << '\n';
  {
    int tag = 1;  // Global element tag
    for (const auto &[elem_type, nodes] : elem_nodes)
    {
      const int elem_type_w = use_lo_type ? LOElemTypeGmsh(elem_type) : elem_type;
      const int &num_elem_nodes = ElemNumNodes[elem_type - 1];
      const int num_elem = (int)nodes.size() / (num_elem_nodes + 1);
#if defined(GMSH_BIN)
      // For binary output, write the element header for each type. Always have 2 tags
      // (physical + geometry)
      const int header[3] = {elem_type_w, num_elem, 2};
      buffer.write(reinterpret_cast<const char *>(header), 3 * sizeof(int));
#endif
      for (int i = 0; i < num_elem; i++)
      {
        WriteElement(buffer, tag++, elem_type_w,
                     nodes[i * (num_elem_nodes + 1)],        // Geometry tag
                     &nodes[i * (num_elem_nodes + 1) + 1]);  // Element nodes
      }
    }
  }
#if defined(GMSH_BIN)
  buffer << '\n';
#endif
  buffer << "$EndElements\n";
}

}  // namespace

namespace mesh
{

void ConvertMeshComsol(const std::string &filename, std::ostream &buffer,
                       bool remove_curvature)
{
  // Read a COMSOL format mesh.
  const int comsol_bin = !filename.compare(filename.length() - 7, 7, ".mphbin") ||
                         !filename.compare(filename.length() - 7, 7, ".MPHBIN");
  MFEM_VERIFY(!filename.compare(filename.length() - 7, 7, ".mphtxt") ||
                  !filename.compare(filename.length() - 7, 7, ".MPHTXT") || comsol_bin,
              "Invalid file extension for COMSOL mesh format conversion!");
  std::ifstream input(filename);
  if (!input.is_open())
  {
    MFEM_ABORT("Unable to open mesh file \"" << filename << "\"!");
  }

  // Parse COMSOL header. COMSOL encodes strings as integer-string pairs where the integer
  // is the string length. It also allows for blank lines and other whitespace wherever in
  // the file.
  {
    int version[2] = {-1, -1};
    int num_tags = -1;
    int num_types = -1;
    if (!comsol_bin)
    {
      while (num_types < 0)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (version[0] < 0)
          {
            sline >> version[0] >> version[1];
          }
          else if (num_tags < 0)
          {
            sline >> num_tags;
            int i = 0;
            while (i < num_tags)
            {
              if (!GetLineComsol(input).empty())
              {
                i++;
              }
            }
          }
          else if (num_types < 0)
          {
            sline >> num_types;
            int i = 0;
            while (i < num_types)
            {
              if (!GetLineComsol(input).empty())
              {
                i++;
              }
            }
          }
        }
      }
    }
    else
    {
      input.read(reinterpret_cast<char *>(version), 2 * sizeof(int));
      input.read(reinterpret_cast<char *>(&num_tags), sizeof(int));
      {
        int i = 0;
        while (i < num_tags)
        {
          ReadStringComsolBinary(input);
          i++;
        }
      }
      input.read(reinterpret_cast<char *>(&num_types), sizeof(int));
      {
        int i = 0;
        while (i < num_types)
        {
          ReadStringComsolBinary(input);
          i++;
        }
      }
    }
    MFEM_VERIFY(version[0] == 0 && version[1] == 1, "Invalid COMSOL file version!");
  }

  // Parse mesh objects until we get to the mesh. Currently only supports a single mesh
  // object in the file, and selections are ignored.
  while (true)
  {
    int object[3] = {-1, -1, -1};
    std::string object_class;
    if (!comsol_bin)
    {
      while (object_class.empty())
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (object[0] < 0)
          {
            sline >> object[0] >> object[1] >> object[2];
          }
          else if (object_class.empty())
          {
            object_class = ReadStringComsol(sline);
          }
        }
      }
    }
    else
    {
      input.read(reinterpret_cast<char *>(object), 3 * sizeof(int));
      object_class = ReadStringComsolBinary(input);
    }
    MFEM_VERIFY(object[0] == 0 && object[1] == 0 && object[2] == 1,
                "Invalid COMSOL object version!");

    // If yes, then ready to parse the mesh.
    if (!object_class.compare(0, 4, "Mesh"))
    {
      break;
    }

    // Otherwise, parse over the selection to the next object.
    MFEM_VERIFY(!object_class.compare(0, 9, "Selection"),
                "COMSOL mesh file only supports Mesh and Selection objects!");
    int version = -1;
    std::string label_str;
    std::string tag_str;
    int sdim = -1;
    int num_ent = -1;
    if (!comsol_bin)
    {
      while (num_ent < 0)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (version < 0)
          {
            sline >> version;
          }
          else if (label_str.empty())
          {
            label_str = ReadStringComsol(sline);
          }
          else if (tag_str.empty())
          {
            tag_str = ReadStringComsol(sline);
          }
          else if (sdim < 0)
          {
            sline >> sdim;
          }
          else if (num_ent < 0)
          {
            sline >> num_ent;
          }
        }
      }
    }
    else
    {
      input.read(reinterpret_cast<char *>(&version), sizeof(int));
      label_str = ReadStringComsolBinary(input);
      tag_str = ReadStringComsolBinary(input);
      input.read(reinterpret_cast<char *>(&sdim), sizeof(int));
      input.read(reinterpret_cast<char *>(&num_ent), sizeof(int));
    }

    // Parse over the entities in the selection.
    int i = 0;
    if (!comsol_bin)
    {
      while (i < num_ent)
      {
        if (!GetLineComsol(input).empty())
        {
          i++;
        }
      }
    }
    else
    {
      while (i < num_ent)
      {
        int dummy;
        input.read(reinterpret_cast<char *>(&dummy), sizeof(int));
        i++;
      }
    }
  }  // Repeat until Mesh is found

  // Parse the mesh object header.
  int sdim = -1;
  int num_nodes = -1;
  int nodes_start = -1;
  {
    int version = -1;
    if (!comsol_bin)
    {
      while (nodes_start < 0)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (version < 0)
          {
            sline >> version;
          }
          else if (sdim < 0)
          {
            sline >> sdim;
          }
          else if (num_nodes < 0)
          {
            sline >> num_nodes;
          }
          else if (nodes_start < 0)
          {
            sline >> nodes_start;
          }
        }
      }
    }
    else
    {
      input.read(reinterpret_cast<char *>(&version), sizeof(int));
      input.read(reinterpret_cast<char *>(&sdim), sizeof(int));
      input.read(reinterpret_cast<char *>(&num_nodes), sizeof(int));
      input.read(reinterpret_cast<char *>(&nodes_start), sizeof(int));
    }
    MFEM_VERIFY(version == 4, "Only COMSOL files with Mesh version 4 are supported!");
    MFEM_VERIFY(sdim == 2 || sdim == 3,
                "COMSOL mesh nodes are required to be in 2D or 3D space!");
    MFEM_VERIFY(num_nodes > 0, "COMSOL mesh file contains no nodes!");
    MFEM_VERIFY(nodes_start >= 0, "COMSOL mesh nodes have a negative starting tag!");
  }

  // Parse mesh nodes.
  std::vector<double> node_coords;
  {
    // Gmsh nodes are always 3D, so initialize to 0.0 in case z-coordinate isn't set.
    node_coords.resize(3 * num_nodes, 0.0);
    int i = 0;
    if (!comsol_bin)
    {
      while (i < num_nodes)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          for (int j = 0; j < sdim; j++)
          {
            sline >> node_coords[3 * i + j];
          }
          i++;
        }
      }
    }
    else
    {
      // Don't read as a single block in case sdim < 3.
      while (i < num_nodes)
      {
        input.read(reinterpret_cast<char *>(node_coords.data() + 3 * i),
                   (std::streamsize)(sdim * sizeof(double)));
        i++;
      }
    }
  }

  // Parse mesh elements. Store for each element of each type: [geometry tag, [node tags]].
  std::unordered_map<int, std::vector<int>> elem_nodes;
  {
    int num_elem_types = -1;
    if (!comsol_bin)
    {
      while (num_elem_types < 0)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (num_elem_types < 0)
          {
            sline >> num_elem_types;
          }
        }
      }
    }
    else
    {
      input.read(reinterpret_cast<char *>(&num_elem_types), sizeof(int));
    }
    MFEM_VERIFY(num_elem_types > 0, "COMSOL mesh file contains no elements!");

    int parsed_types = 0;  // COMSOL groups elements by type in file
    int elem_type = -1;
    int num_elem_nodes = -1;
    int num_elem = -1;
    int num_elem_geom = -1;
    bool skip_type = false;
    while (parsed_types < num_elem_types)
    {
      if (!comsol_bin)
      {
        auto line = GetLineComsol(input);
        if (!line.empty())
        {
          std::istringstream sline(line);
          if (elem_type < 0)
          {
            auto elem_str = ReadStringComsol(sline);
            MFEM_VERIFY(!elem_str.empty(),
                        "Unexpected empty element type found in COMSOL mesh file!");
            elem_type = ElemTypeComsol(elem_str);
            skip_type = (elem_type == 0);
            MFEM_VERIFY(skip_type || elem_nodes.find(elem_type) == elem_nodes.end(),
                        "Duplicate element types found in COMSOL mesh file!");
          }
          else if (num_elem_nodes < 0)
          {
            sline >> num_elem_nodes;
            MFEM_VERIFY(num_elem_nodes > 0,
                        "COMSOL element type " << elem_type << " has no nodes!");
            MFEM_VERIFY(skip_type || num_elem_nodes == ElemNumNodes[elem_type - 1],
                        "Mismatch between COMSOL and Gmsh element types!");
          }
          else if (num_elem < 0)
          {
            sline >> num_elem;
            MFEM_VERIFY(num_elem > 0,
                        "COMSOL mesh file has no elements of type " << elem_type << "!");
            std::vector<int> *data = nullptr;
            if (!skip_type)
            {
              data = &elem_nodes[elem_type];
              data->resize(num_elem * (num_elem_nodes + 1));  // Node tags + geometry tag
            }

            // Parse all element nodes.
            int i = 0;
            while (i < num_elem)
            {
              line = GetLineComsol(input);
              if (!line.empty())
              {
                if (!skip_type)
                {
                  std::istringstream isline(line);
                  for (int j = 0; j < num_elem_nodes; j++)
                  {
                    // Permute and reset to 1-based node tags.
                    const int &p = ElemNodesComsol[elem_type - 1][j];
                    isline >> (*data)[i * (num_elem_nodes + 1) + 1 + p];
                    (*data)[i * (num_elem_nodes + 1) + 1 + p] += (1 - nodes_start);
                  }
                }
                i++;
              }
            }
          }
          else if (num_elem_geom < 0)
          {
            sline >> num_elem_geom;
            MFEM_VERIFY(num_elem_geom == num_elem,
                        "COMSOL mesh file should have geometry tags for all elements!");
            std::vector<int> *data = nullptr;
            if (!skip_type)
            {
              MFEM_VERIFY(elem_nodes.find(elem_type) != elem_nodes.end(),
                          "Can't find expected element type!");
              data = &elem_nodes[elem_type];
              MFEM_VERIFY(data->size() == (std::size_t)num_elem * (num_elem_nodes + 1),
                          "Unexpected element data size!");
            }

            // Parse all element geometry tags (stored at beginning of element nodes). For
            // geometric entites in < 3D, the exported COMSOL tags are 0-based and need
            // correcting to 1-based for Gmsh.
            int i = 0;
            const int geom_start =
                (elem_type < 4 || (elem_type > 7 && elem_type < 11)) ? 1 : 0;
            while (i < num_elem)
            {
              line = GetLineComsol(input);
              if (!line.empty())
              {
                if (!skip_type)
                {
                  std::istringstream ssline(line);
                  ssline >> (*data)[i * (num_elem_nodes + 1)];
                  (*data)[i * (num_elem_nodes + 1)] += geom_start;
                }
                i++;
              }
            }

            // Debug
            // std::cout << "Finished parsing " << num_elem
            //           << " elements with type " << elem_type
            //           << " (parsed types " << parsed_types + 1 << ")\n";

            // Finished with this element type, on to the next.
            parsed_types++;
            elem_type = num_elem_nodes = num_elem = num_elem_geom = -1;
            skip_type = false;
          }
        }
      }
      else
      {
        auto elem_str = ReadStringComsolBinary(input);
        MFEM_VERIFY(!elem_str.empty(),
                    "Unexpected empty element type found in COMSOL mesh file!");
        elem_type = ElemTypeComsol(elem_str);
        skip_type = (elem_type == 0);
        MFEM_VERIFY(skip_type || elem_nodes.find(elem_type) == elem_nodes.end(),
                    "Duplicate element types found in COMSOL mesh file!");
        input.read(reinterpret_cast<char *>(&num_elem_nodes), sizeof(int));
        MFEM_VERIFY(num_elem_nodes > 0,
                    "COMSOL element type " << elem_type << " has no nodes!");
        MFEM_VERIFY(skip_type || num_elem_nodes == ElemNumNodes[elem_type - 1],
                    "Mismatch between COMSOL and Gmsh element types!");

        // Parse all element nodes.
        input.read(reinterpret_cast<char *>(&num_elem), sizeof(int));
        MFEM_VERIFY(num_elem > 0,
                    "COMSOL mesh file has no elements of type " << elem_type << "!");
        std::vector<int> *data = nullptr;
        if (!skip_type)
        {
          data = &elem_nodes[elem_type];
          data->resize(num_elem * (num_elem_nodes + 1));  // Node tags + geometry tag
        }
        int i = 0;
        std::vector<int> nodes(num_elem_nodes);
        while (i < num_elem)
        {
          input.read(reinterpret_cast<char *>(nodes.data()),
                     (std::streamsize)(num_elem_nodes * sizeof(int)));
          if (!skip_type)
          {
            for (int j = 0; j < num_elem_nodes; j++)
            {
              // Permute and reset to 1-based node tags.
              const int &p = ElemNodesComsol[elem_type - 1][j];
              (*data)[i * (num_elem_nodes + 1) + 1 + p] = nodes[j] + (1 - nodes_start);
            }
          }
          i++;
        }

        // Parse element geometry tags.
        input.read(reinterpret_cast<char *>(&num_elem_geom), sizeof(int));
        MFEM_VERIFY(num_elem_geom == num_elem,
                    "COMSOL mesh file should have geometry tags for all elements!");

        i = 0;
        const int geom_start = (elem_type < 4 || (elem_type > 7 && elem_type < 11)) ? 1 : 0;
        int geom_tag;
        while (i < num_elem)
        {
          input.read(reinterpret_cast<char *>(&geom_tag), sizeof(int));
          if (!skip_type)
          {
            (*data)[i * (num_elem_nodes + 1)] = geom_tag + geom_start;
          }
          i++;
        }

        // Debug
        // std::cout << "Finished parsing " << num_elem
        //           << " elements with type " << elem_type
        //           << " (parsed types " << parsed_types + 1 << ")\n";

        // Finished with this element type, on to the next.
        parsed_types++;
        elem_type = num_elem_nodes = num_elem = num_elem_geom = -1;
        skip_type = false;
      }
    }
  }

  // Finalize input, write the Gmsh mesh.
  input.close();
  std::vector<int> dummy;
  WriteGmsh(buffer, node_coords, dummy, elem_nodes, remove_curvature);
}

void ConvertMeshNastran(const std::string &filename, std::ostream &buffer,
                        bool remove_curvature)
{
  // Read a Nastran/BDF format mesh.
  MFEM_VERIFY(!filename.compare(filename.length() - 4, 4, ".nas") ||
                  !filename.compare(filename.length() - 4, 4, ".NAS") ||
                  !filename.compare(filename.length() - 4, 4, ".bdf") ||
                  !filename.compare(filename.length() - 4, 4, ".BDF"),
              "Invalid file extension for Nastran mesh format conversion!");
  std::ifstream input(filename);
  if (!input.is_open())
  {
    MFEM_ABORT("Unable to open mesh file \"" << filename << "\"!");
  }
  const int NASTRAN_CHUNK = 8;  // NASTRAN divides row into 10 columns of 8 spaces
  const int MAX_CHUNK = 9;      // Never read the 10-th chunk

  // Parse until bulk data starts.
  while (true)
  {
    auto line = GetLineNastran(input);
    if (line.length() > 0)
    {
      if (!line.compare(0, 10, "BEGIN BULK"))
      {
        break;
      }
    }
  }

  // Parse mesh nodes and elements. It is expected that node tags start at 1 and are
  // contiguous. Store for each element of each type: [geometry tag, [node tags]].
  std::vector<double> node_coords;
  std::vector<int> node_tags;
  std::unordered_map<int, std::vector<int>> elem_nodes;
  int elem_type;
  while (true)
  {
    auto line = GetLineNastran(input);
    if (line.length() > 0 && !input.eof())
    {
      if (!line.compare(0, 7, "ENDDATA"))
      {
        break;  // Done parsing file
      }
      else if (!line.compare(0, 5, "GRID*"))
      {
        // Coordinates in long field format (8 + 16 * 4 + 8).
        auto next = GetLineNastran(input);
        MFEM_VERIFY(!next.empty(), "Unexpected empty line parsing Nastran!");

        node_tags.push_back(std::stoi(line.substr(1 * NASTRAN_CHUNK, 2 * NASTRAN_CHUNK)));
        node_coords.insert(
            node_coords.end(),
            {ConvertDoubleNastran(line.substr(5 * NASTRAN_CHUNK, 2 * NASTRAN_CHUNK)),
             ConvertDoubleNastran(line.substr(7 * NASTRAN_CHUNK, 2 * NASTRAN_CHUNK)),
             ConvertDoubleNastran(next.substr(1 * NASTRAN_CHUNK, 2 * NASTRAN_CHUNK))});
      }
      else if (!line.compare(0, 4, "GRID"))
      {
        if (line.find_first_of(',') != std::string::npos)
        {
          // Free field format (comma separated).
          std::istringstream sline(line);

          std::string word;
          std::getline(sline, word, ',');  // Discard "GRID"

          std::getline(sline, word, ',');
          node_tags.push_back(std::stoi(word));

          std::getline(sline, word, ',');  // Discard coordinate system

          std::getline(sline, word, ',');
          double x = ConvertDoubleNastran(word);
          std::getline(sline, word, ',');
          double y = ConvertDoubleNastran(word);
          std::getline(sline, word, ',');
          double z = ConvertDoubleNastran(word);
          node_coords.insert(node_coords.end(), {x, y, z});
        }
        else
        {
          // Short format (10 * 8).
          node_tags.push_back(std::stoi(line.substr(1 * NASTRAN_CHUNK, NASTRAN_CHUNK)));
          node_coords.insert(
              node_coords.end(),
              {ConvertDoubleNastran(line.substr(3 * NASTRAN_CHUNK, NASTRAN_CHUNK)),
               ConvertDoubleNastran(line.substr(4 * NASTRAN_CHUNK, NASTRAN_CHUNK)),
               ConvertDoubleNastran(line.substr(5 * NASTRAN_CHUNK, NASTRAN_CHUNK))});
        }
      }
      else if ((elem_type = ElemTypeNastran(line)))
      {
        // Prepare to parse the element ID and nodes.
        const bool free = (line.find_first_of(',') != std::string::npos);

        // Get the element type, tag, and geometry attribute. Then get the element nodes on
        // this line.
        std::string elem_str;
        // int elem_tag;
        int geom_tag;
        std::vector<int> nodes;
        std::string word;
        if (!free)
        {
          elem_str = line.substr(0 * NASTRAN_CHUNK, NASTRAN_CHUNK);
          const std::size_t stop = elem_str.find_last_not_of(' ');
          MFEM_VERIFY(stop != std::string::npos, "Invalid element type string!");
          elem_str.resize(stop + 1);
          // elem_tag = std::stoi(line.substr(1*NASTRAN_CHUNK, NASTRAN_CHUNK));
          geom_tag = std::stoi(line.substr(2 * NASTRAN_CHUNK, NASTRAN_CHUNK));

          int i = 3;
          while (i < MAX_CHUNK)
          {
            word = line.substr((i++) * NASTRAN_CHUNK, NASTRAN_CHUNK);
            if (word.find_first_not_of(' ') == std::string::npos)
            {
              break;
            }
            nodes.push_back(std::stoi(word));
          }
        }
        else
        {
          std::istringstream sline(line);
          std::getline(sline, elem_str, ',');
          std::getline(sline, word, ',');
          // elem_tag = std::stoi(word);
          std::getline(sline, word, ',');
          geom_tag = std::stoi(word);

          int i = 3;
          while (i < MAX_CHUNK)
          {
            std::getline(sline, word, ',');
            if (word.find_first_not_of(' ') == std::string::npos)
            {
              break;
            }
            nodes.push_back(std::stoi(word));
            i++;
          }
        }

        // Handle line continuation.
        while (input.peek() == '+')
        {
          auto next = GetLineNastran(input);
          MFEM_VERIFY(!next.empty(), "Unexpected empty line parsing Nastran!");

          if (!free)
          {
            int i = 1;
            while (i < MAX_CHUNK)
            {
              word = next.substr((i++) * NASTRAN_CHUNK, NASTRAN_CHUNK);
              if (word.find_first_not_of(' ') == std::string::npos)
              {
                break;
              }
              nodes.push_back(std::stoi(word));
            }
          }
          else
          {
            std::istringstream snext(next);
            int i = 1;
            while (i < MAX_CHUNK)
            {
              std::getline(snext, word, ',');
              if (word.find_first_not_of(' ') == std::string::npos)
              {
                break;
              }
              nodes.push_back(std::stoi(word));
              i++;
            }
          }
        }

        // Save the element and its geometry tag.
        elem_type = HOElemTypeNastran(elem_type, (int)nodes.size());
        const int &num_elem_nodes = ElemNumNodes[elem_type - 1];
        MFEM_VERIFY((std::size_t)num_elem_nodes == nodes.size(),
                    "Mismatch between Nastran and Gmsh element types!");
        std::vector<int> &data = elem_nodes[elem_type];
        const int i = (int)data.size();
        data.resize(i + 1 + num_elem_nodes);
        data[i] = geom_tag;
        for (int j = 0; j < num_elem_nodes; j++)
        {
          // Permute back to Gmsh ordering.
          const int &p = ElemNodesNastran[elem_type - 1][j];
          data[i + 1 + p] = nodes[j];
        }
      }
    }
  }

  // Finalize input, write the Gmsh mesh.
  input.close();
  WriteGmsh(buffer, node_coords, node_tags, elem_nodes, remove_curvature);
}

}  // namespace mesh

}  // namespace palace
