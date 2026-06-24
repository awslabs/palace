// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP
#define PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

// ParaView data collection with an extra registration path for point data that has
// already been evaluated in the exact order used by MFEM's boundary VTU writer. This
// avoids pretending precomputed libCEED output is an mfem::Coefficient and, more
// importantly, avoids recovering point identity from floating-point reference
// coordinates during Save().
class CeedParaViewDataCollection : public mfem::ParaViewDataCollection
{
private:
  struct BoundaryPointField
  {
    const Vector *values = nullptr;
    const std::vector<int> *bases = nullptr;
    int num_comp = 0;
  };

  std::fstream pvd_stream;
  std::map<std::string, BoundaryPointField> boundary_point_fields;

  void SaveDataVTU(std::ostream &os, int ref);
  void SaveBoundaryPointFieldVTU(std::ostream &os, int ref, const std::string &name,
                                 const BoundaryPointField &field);

public:
  using mfem::ParaViewDataCollection::ParaViewDataCollection;

  void RegisterBoundaryPointField(const std::string &field_name, const Vector &values,
                                  const std::vector<int> &bases, int num_comp);
  void DeregisterBoundaryPointField(const std::string &field_name);

  void Save() override;
};

}  // namespace palace

#endif  // PALACE_UTILS_CEED_PARAVIEW_DATA_COLLECTION_HPP
