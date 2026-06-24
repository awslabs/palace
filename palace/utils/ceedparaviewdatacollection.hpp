// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP
#define PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP

#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

// ParaView data collection with an extra registration path for point data that has
// already been evaluated in the exact order used by MFEM's VTU writer. This lets
// libCEED fill the VTU PointData buffer directly (domain or boundary) instead of
// materializing a GridFunction first or re-entering mfem::Coefficient evaluation during
// Save().
class CeedParaViewDataCollection : public mfem::ParaViewDataCollection
{
private:
  struct PointField
  {
    const Vector *values = nullptr;
    const std::vector<int> *bases = nullptr;
    int num_comp = 0;
  };

  std::fstream pvd_stream;
  std::map<std::string, PointField> domain_point_fields;
  std::map<std::string, PointField> boundary_point_fields;

  void WritePVTUHeader(std::ostream &os, bool appended_mesh);
  void SaveMeshVTU(std::ostream &os, int ref,
                   std::vector<std::vector<char>> *appended_blocks,
                   std::size_t *appended_offset);
  void SaveDataVTU(std::ostream &os, int ref);
  void SaveGFieldVTU(std::ostream &os, int ref, const FieldMapIterator &it,
                     std::vector<std::vector<char>> *appended_blocks,
                     std::size_t *appended_offset);
  void SavePointFieldVTU(std::ostream &os, int ref, const std::string &name,
                         const PointField &field, bool boundary,
                         std::vector<std::vector<char>> *appended_blocks,
                         std::size_t *appended_offset);

public:
  using mfem::ParaViewDataCollection::ParaViewDataCollection;

  void RegisterDomainPointField(const std::string &field_name, const Vector &values,
                                const std::vector<int> &bases, int num_comp);
  void DeregisterDomainPointField(const std::string &field_name);

  void RegisterBoundaryPointField(const std::string &field_name, const Vector &values,
                                  const std::vector<int> &bases, int num_comp);
  void DeregisterBoundaryPointField(const std::string &field_name);

  void Save() override;
};

}  // namespace palace

#endif  // PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP
