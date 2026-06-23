// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_MODEL_HPP
#define PALACE_SCHEMA_TYPES_MODEL_HPP

// Mirrors palace::config::{Model, Refinement, BoxRefinement,
// SphereRefinement} (palace/utils/configfile.hpp). Descriptions match
// PR 716's scripts/schema/config/model.json.

#include <array>
#include <string>
#include <vector>

#include <rfl.hpp>

#include "common.hpp"
#include "schema/utils/annotations.hpp"

namespace palace::schema
{

struct Box
{
  PALACE_SCHEMA_DESC_REQUIRED(Levels,
                              "Levels of parallel mesh refinement inside this box region.",
                              palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC_REQUIRED(
      BoundingBoxMin,
      "Minimum coordinates `[x, y, z]` of the axis-aligned bounding box for "
      "this refinement region, in mesh length units.",
      Vector3) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC_REQUIRED(
      BoundingBoxMax,
      "Maximum coordinates `[x, y, z]` of the axis-aligned bounding box for "
      "this refinement region, in mesh length units.",
      Vector3) = {
    {0.0, 0.0, 0.0}
  };
};

struct Sphere
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Levels, "Levels of parallel mesh refinement inside this sphere region.",
      palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC_REQUIRED(Radius, "Radius of the sphere, in mesh length units.",
                              palace::schema::utils::XMin<double, 0>) = 1.0;

  PALACE_SCHEMA_DESC_REQUIRED(
      Center, "Center coordinates `[x, y, z]` of the sphere, in mesh length units.",
      Vector3) = {
    {0.0, 0.0, 0.0}
  };
};

struct Refinement
{
  PALACE_SCHEMA_DESC(Tol,
                     "Stop adaptive mesh refinement (AMR) when the norm of the estimated "
                     "error falls below this value. The error is reported in "
                     "`error-indicators.csv`.",
                     palace::schema::utils::XMin<double, 0>) = 1.0e-2;

  PALACE_SCHEMA_DESC(MaxIts, "Maximum number of AMR iterations to perform.",
                     palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC(MaxSize,
                     "The maximum allowable number of degrees of freedom for AMR. If an "
                     "adapted mesh exceeds this value no further adaptation will occur. A "
                     "value less than 1 means that no maximum size constraint will be "
                     "imposed.",
                     palace::schema::utils::Min<double, 0>) = 0.0;

  PALACE_SCHEMA_DESC(UpdateFraction,
                     "Dörfler marking fraction used to specify which elements to refine. "
                     "This marking strategy will mark the smallest number of elements that "
                     "make up \"UpdateFraction\" of the total error in the mesh. A larger "
                     "value will refine more elements per iteration, at the cost of the "
                     "final mesh being less efficient.",
                     palace::schema::utils::Open<double, 0, 1>) = 0.7;

  PALACE_SCHEMA_DESC(Nonconformal,
                     "Use nonconformal refinement in adaptation. Required for non-simplex "
                     "meshes.",
                     bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(
      MaxNCLevels, "Maximum number of nonconformal refinement levels. `0` means no limit.",
      palace::schema::utils::Min<int, 0>) = 1;

  PALACE_SCHEMA_DESC_ADVANCED(
      MaximumImbalance,
      "Maximum ratio of elements between the most- and least-loaded MPI "
      "ranks before repartitioning.",
      palace::schema::utils::Min<double, 1>) = 1.1;

  PALACE_SCHEMA_DESC_ADVANCED(SaveAdaptIterations,
                              "Save postprocessing results from each AMR iteration in a "
                              "subdirectory `iterationX`.",
                              bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(SaveAdaptMesh, "Save the final adapted mesh to disk.",
                              bool) = false;

  PALACE_SCHEMA_DESC(
      UniformLevels,
      "Levels of uniform parallel mesh refinement to be performed on the "
      "input mesh. If not performing AMR, these may be used as levels within "
      "a geometric multigrid scheme. If performing AMR the most refined mesh "
      "is used as the initial mesh and the coarser meshes cannot be used in "
      "a geometric multigrid scheme.",
      palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC_ADVANCED(
      SerialUniformLevels,
      "Levels of uniform serial mesh refinement applied before parallel "
      "distribution.",
      palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC(Boxes,
                     "Array of axis-aligned box refinement regions. All elements with a "
                     "node inside the box are marked for refinement.",
                     std::vector<Box>) = {};

  PALACE_SCHEMA_DESC(Spheres,
                     "Array of sphere refinement regions. All elements with a node inside "
                     "the sphere are marked for refinement.",
                     std::vector<Sphere>) = {};
};

struct Model
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Mesh,
      "Input mesh file path. An absolute path is recommended. If the "
      "provided mesh is nonconformal, it is assumed to come from a previous "
      "*Palace* AMR solve, and all mesh preprocessing checks and "
      "modifications (for example "
      "[CrackInternalBoundaryElements](@ref "
      "config-model-crackinternalboundaryelements)) are skipped.",
      std::string) = "";

  PALACE_SCHEMA_DESC(L0,
                     "Unit, relative to meters, for mesh vertex coordinates. For example, "
                     "a value of `1.0e-6` means the mesh coordinates are in μm.",
                     palace::schema::utils::XMin<double, 0>) = 1.0e-6;

  PALACE_SCHEMA_DESC(Lc,
                     "Characteristic length scale used for nondimensionalization, "
                     "specified in mesh length units. This keyword should typically not be "
                     "specified by the user. A value less than or equal to zero uses an "
                     "internally calculated length scale based on the bounding box of the "
                     "computational domain. A value of `1.0` will disable "
                     "nondimensionalization, so that all computations will take place in "
                     "the same units as the mesh.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC_ADVANCED(
      RemoveCurvature,
      "Project high-order nodes to the mesh surface, removing all curvature "
      "before the simulation.",
      bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(
      MakeSimplex, "Convert all mesh elements to simplices (tetrahedra/triangles).",
      bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(MakeHexahedral, "Convert all mesh elements to hexahedra.",
                              bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(ReorderElements,
                              "Reorder mesh elements to improve cache efficiency.",
                              bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(CleanUnusedElements,
                              "Remove elements not connected to any domain material.",
                              bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(
      CrackInternalBoundaryElements,
      "Duplicate nodes along internal boundary elements to create a crack.", bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(RefineCrackElements,
                              "Refine elements adjacent to cracked internal boundaries.",
                              bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(
      CrackDisplacementFactor,
      "Displacement factor applied to cracked nodes as a fraction of the "
      "local element size.",
      palace::schema::utils::Min<double, 0>) = 1.0e-12;

  PALACE_SCHEMA_DESC_ADVANCED(
      AddInterfaceBoundaryElements,
      "Add boundary elements at interfaces between domains that lack them.", bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(
      ExportPrerefinedMesh,
      "Export the mesh after preprocessing but before AMR refinement.", bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(ReorientTetMesh,
                              "Reorient tetrahedral elements to ensure positive Jacobians.",
                              bool) = false;

  PALACE_SCHEMA_DESC_ADVANCED(
      Partitioning,
      "Path to a mesh partitioning file. If empty, partitioning is computed "
      "automatically.",
      std::string) = "";

  PALACE_SCHEMA_DESC(Refinement, "Configuration for adaptive and uniform mesh refinement.",
                     Refinement) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_MODEL_HPP
