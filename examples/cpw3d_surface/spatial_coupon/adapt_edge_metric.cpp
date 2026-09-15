// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Bounded native-MMG scout bridge. The caller owns geometry/metric validation,
// resource guards and publication. No partial/failed adaptation is accepted.
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <mmg/mmg3d/libmmg3d.h>

void Check(int result, const char *operation)
{
  if (result != 1)
  {
    throw std::runtime_error(operation);
  }
}

int main(int argc, char **argv)
{
  if (argc < 8 || argc > 11)
  {
    std::cerr << "input.meshb metric.f64 pins.txt output.meshb hmin hmax hgrad "
                 "[no-move|freeze-surface|freeze-matching|freeze-matching-no-move|freeze-"
                 "selected] "
                 "[fixed-triangles.txt [hausd]]\n";
    return 2;
  }
  const std::string mode = argc >= 9 ? argv[8] : "adapt";
  const bool freeze_matching = mode == "freeze-matching" ||
                               mode == "freeze-matching-no-move" ||
                               mode == "freeze-selected";
  if ((freeze_matching && argc != 10 && argc != 11) || (!freeze_matching && argc >= 10) ||
      (mode != "adapt" && mode != "no-move" && mode != "freeze-surface" &&
       !freeze_matching))
  {
    return 2;
  }
  if (std::filesystem::exists(argv[4]) ||
      std::filesystem::exists(std::string(argv[4]) + ".rejected.meshb"))
  {
    std::cerr << "Refusing to overwrite an adaptation attempt\n";
    return 2;
  }
  MMG5_pMesh mesh = nullptr;
  MMG5_pSol metric = nullptr;
  int status = 2;
  try
  {
    const double hmin = std::stod(argv[5]);
    const double hmax = std::stod(argv[6]);
    const double hgrad = std::stod(argv[7]);
    const double hausd = argc == 11 ? std::stod(argv[10]) : 1e-9;
    if (!std::isfinite(hmin) || !std::isfinite(hmax) || !std::isfinite(hgrad) ||
        !std::isfinite(hausd) || hausd <= 0 || hmin <= 0 || hmax < hmin ||
        (hgrad != -1 && hgrad <= 1))
    {
      throw std::runtime_error("Invalid metric size/gradation controls");
    }
    Check(MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric,
                          MMG5_ARG_end),
          "init");
    Check(MMG3D_loadMesh(mesh, argv[1]), "load mesh");
    int nv, ne, nprism, nt, nquad, nedge;
    Check(MMG3D_Get_meshSize(mesh, &nv, &ne, &nprism, &nt, &nquad, &nedge), "size");
    if (nv <= 0 || ne <= 0 || nt <= 0 || nprism || nquad)
    {
      throw std::runtime_error("Expected a tetrahedral multi-surface mesh");
    }
    std::vector<double> tensors(6ULL * nv);
    std::ifstream input(argv[2], std::ios::binary);
    input.read(reinterpret_cast<char *>(tensors.data()), tensors.size() * sizeof(double));
    if (!input || input.peek() != std::char_traits<char>::eof())
    {
      throw std::runtime_error("Metric byte count mismatch");
    }
    for (int i = 0; i < nv; i++)
    {
      const double *m = tensors.data() + 6ULL * i;
      const double minor = m[0] * m[3] - m[1] * m[1];
      const double det = m[0] * (m[3] * m[5] - m[4] * m[4]) -
                         m[1] * (m[1] * m[5] - m[4] * m[2]) +
                         m[2] * (m[1] * m[4] - m[3] * m[2]);
      for (int j = 0; j < 6; j++)
      {
        if (!std::isfinite(m[j]))
        {
          throw std::runtime_error("Nonfinite metric");
        }
      }
      if (!(m[0] > 0 && minor > 0 && det > 0))
      {
        throw std::runtime_error("Metric not SPD");
      }
    }
    Check(MMG3D_Set_solSize(mesh, metric, MMG5_Vertex, nv, MMG5_Tensor), "metric size");
    // Confirmed native API order: m11,m12,m13,m22,m23,m33, not Medit .sol order.
    Check(MMG3D_Set_tensorSols(metric, tensors.data()), "set tensors");
    std::vector<double> roundtrip(tensors.size());
    Check(MMG3D_Get_tensorSols(metric, roundtrip.data()), "get tensors");
    if (roundtrip != tensors)
    {
      throw std::runtime_error("Tensor API round trip failed");
    }
    for (int i = 1; i <= nedge; i++)
    {
      Check(MMG3D_Set_ridge(mesh, i), "ridge");
    }
    std::ifstream pins(argv[3]);
    int vertex, pin_count = 0;
    while (pins >> vertex)
    {
      if (vertex < 1 || vertex > nv)
      {
        throw std::runtime_error("Invalid pin");
      }
      Check(MMG3D_Set_corner(mesh, vertex), "corner");
      Check(MMG3D_Set_requiredVertex(mesh, vertex), "required vertex");
      pin_count++;
    }
    if (!pins.eof() || !pin_count)
    {
      throw std::runtime_error("Missing or malformed geometric pins");
    }
    if (mode == "no-move" || mode == "freeze-matching-no-move")
    {
      Check(MMG3D_Set_iparameter(mesh, metric, MMG3D_IPARAM_nomove, 1), "no move");
    }
    if (mode == "freeze-surface")
    {
      std::vector<int> triangles(3ULL * nt), refs(nt), required(nt);
      Check(MMG3D_Get_triangles(mesh, triangles.data(), refs.data(), required.data()),
            "get fixed surface");
      for (int i = 1; i <= nt; i++)
      {
        Check(MMG3D_Set_requiredTriangle(mesh, i), "required surface triangle");
      }
      for (const int node : triangles)
      {
        Check(MMG3D_Set_requiredVertex(mesh, node), "required surface vertex");
      }
      Check(MMG3D_Set_iparameter(mesh, metric, MMG3D_IPARAM_nosurf, 1), "fixed surface");
      std::cout << "Surface frozen; interior vertex movement enabled\n";
    }
    if (freeze_matching)
    {
      std::vector<int> triangles(3ULL * nt), refs(nt), required(nt);
      Check(MMG3D_Get_triangles(mesh, triangles.data(), refs.data(), required.data()),
            "get matching surface");
      std::ifstream fixed(argv[9]);
      int triangle, fixed_count = 0;
      while (fixed >> triangle)
      {
        if (triangle < 1 || triangle > nt)
        {
          throw std::runtime_error("Invalid fixed triangle index");
        }
        Check(MMG3D_Set_requiredTriangle(mesh, triangle), "fixed matching triangle");
        for (int j = 0; j < 3; j++)
        {
          Check(MMG3D_Set_requiredVertex(mesh, triangles[3ULL * (triangle - 1) + j]),
                "fixed matching vertex");
        }
        fixed_count++;
      }
      if (!fixed.eof() || !fixed_count)
      {
        throw std::runtime_error("Missing or malformed matching triangle list");
      }
      std::cout << "Selected surface triangles frozen: " << fixed_count << std::endl;
    }
    Check(MMG3D_Set_iparameter(mesh, metric, MMG3D_IPARAM_verbose, 4), "verbose");
    Check(MMG3D_Set_iparameter(mesh, metric, MMG3D_IPARAM_mem, 6000), "memory");
    Check(MMG3D_Set_iparameter(mesh, metric, MMG3D_IPARAM_nosizreq, 1), "required sizes");
    Check(MMG3D_Set_dparameter(mesh, metric, MMG3D_DPARAM_hmin, hmin), "hmin");
    Check(MMG3D_Set_dparameter(mesh, metric, MMG3D_DPARAM_hmax, hmax), "hmax");
    Check(MMG3D_Set_dparameter(mesh, metric, MMG3D_DPARAM_hgrad, hgrad), "hgrad");
    if (hgrad == -1)
    {
      Check(MMG3D_Set_dparameter(mesh, metric, MMG3D_DPARAM_hgradreq, -1),
            "disable required-entity gradation");
    }
    Check(MMG3D_Set_dparameter(mesh, metric, MMG3D_DPARAM_hausd, hausd), "hausd");
    std::cout << "Native tensor roundtrip passed; vertices=" << nv << " tets=" << ne
              << " triangles=" << nt << " feature edges=" << nedge << " pins=" << pin_count
              << std::endl;
    const int result = MMG3D_mmg3dlib(mesh, metric);
    std::cout << "MMG_RESULT=" << result << std::endl;
    if (result == MMG5_SUCCESS)
    {
      Check(MMG3D_saveMesh(mesh, argv[4]), "save mesh");
      Check(MMG3D_saveSol(mesh, metric, (std::string(argv[4]) + ".sol").c_str()),
            "save effective metric");
      status = 0;
    }
    else
    {
      const std::string rejected = std::string(argv[4]) + ".rejected.meshb";
      MMG3D_saveMesh(mesh, rejected.c_str());
      std::cerr << "Adaptation rejected; no accepted output\n";
    }
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << std::endl;
  }
  if (mesh || metric)
  {
    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric,
                   MMG5_ARG_end);
  }
  return status;
}
