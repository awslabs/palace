// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"

#include <array>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <mfem.hpp>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Floating point precision for mesh IO. This precision is important, make sure nothing is
// lost!
const auto MSH_FLT_PRECISION = std::numeric_limits<double>::max_digits10;

// Load the serial mesh from disk.
mfem::Mesh LoadMesh(const std::string &);

// Optionally reorder mesh elements based on MFEM's internal reordeing tools for improved
// cache usage.
void ReorderMesh(mfem::Mesh &);

// Generate element-based mesh partitioning, using either a provided file or METIS.
std::unique_ptr<int[]> GetMeshPartitioning(mfem::Mesh &, int, const std::string &);

// Cleanup the provided serial mesh by removing unnecessary domain and elements, adding
// boundary elements for material interfaces and exterior boundaries, and adding boundary
// elements for subdomain interfaces.
std::map<int, std::array<int, 2>> CheckMesh(mfem::Mesh &, const std::unique_ptr<int[]> &,
                                            const IoData &, bool, bool, bool);

// Given a serial mesh on the root processor and element partitioning, create a parallel
// mesh over the given communicator.
mfem::ParMesh DistributeMesh(MPI_Comm, mfem::Mesh &, const std::unique_ptr<int[]> &);

// Get list of domain and boundary attribute markers used in configuration file for mesh
// cleaning.
void GetUsedAttributeMarkers(const IoData &, int, int, mfem::Array<int> &,
                             mfem::Array<int> &);

}  // namespace

namespace mesh
{

mfem::ParMesh ReadMesh(MPI_Comm comm, const IoData &iodata, bool reorder, bool clean,
                       bool add_bdr, bool unassembled, Timer &timer)
{
  // On root, read the serial mesh (converting format if necessary), and do all necessary
  // serial preprocessing. When finished, distribute the mesh to all processes. Count disk
  // I/O time separately for the mesh read from file.
  auto t0 = timer.Now();

  // Optionally reorder elements (and vertices) based on spatial location after loading
  // the serial mesh.
  auto mesh = LoadMesh(iodata.model.mesh);
  if (reorder)
  {
    ReorderMesh(mesh);
  }

  Mpi::Barrier(comm);
  timer.io_time += timer.Now() - t0;

  // Generate the parallel mesh partitioning on the root process.
  const auto partitioning =
      GetMeshPartitioning(mesh, Mpi::Size(comm), iodata.model.partition);

  // Clean up unused domain elements from the mesh, add new boundary elements for material
  // interfaces if not present, and optionally (when running unassembled) add subdomain
  // interface boundary elements.
  static_cast<void>(CheckMesh(mesh, partitioning, iodata, clean, add_bdr, unassembled));

  // Construct the parallel mesh data structure by distributing the serial mesh from the
  // root process. The serial mesh and partitioning are deleted inside.
  return DistributeMesh(comm, mesh, partitioning);
}

void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
{
  // Prepare for uniform and region-based refinement.
  MFEM_VERIFY(mesh.size() == 1,
              "Input mesh vector before refinement has more than a single mesh!");
  int uniform_ref_levels = iodata.model.refinement.uniform_ref_levels;
  int max_region_ref_levels = 0;
  for (const auto &box : iodata.model.refinement.GetBoxes())
  {
    if (max_region_ref_levels < box.ref_levels)
    {
      max_region_ref_levels = box.ref_levels;
    }
  }
  for (const auto &sphere : iodata.model.refinement.GetSpheres())
  {
    if (max_region_ref_levels < sphere.ref_levels)
    {
      max_region_ref_levels = sphere.ref_levels;
    }
  }
  if (iodata.solver.linear.mat_gmg)
  {
    mesh.reserve(1 + uniform_ref_levels + max_region_ref_levels);
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the coarse mesh so that the refinements are still true refinements of
  // the original mesh (required for geometric multigrid). Otherwise, it happens after
  // refinement.
  if (iodata.model.reorient_tet && mesh.capacity() > 1)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    mesh[0]->ReorientTetMesh();
#pragma GCC diagnostic pop
  }

  // Uniformly refine the mesh further in parallel, saving the level meshes for geometric
  // coarsening later on if desired.
  for (int l = 0; l < uniform_ref_levels; l++)
  {
    if (mesh.capacity() > 1)
    {
      mesh.push_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
    }
    mesh.back()->UniformRefinement();
  }

  // Proceed with region-based refinement, level-by-level for all regions. Currently support
  // box and sphere region shapes. Any overlap between regions is ignored (take the union,
  // don't double-refine).
  if (max_region_ref_levels > 0 &&
      (mesh[0]->MeshGenerator() & 2 || mesh[0]->MeshGenerator() & 4 ||
       mesh[0]->MeshGenerator() & 8))
  {
    // XX TODO: Region-based refinement won't work if the ParMesh has been constructed from
    //          a conforming mesh, but non-conforming refinement is needed. Unclear if the
    //          current mesh distribution scheme will work even for a conforming serial mesh
    //          which is a NCMesh after Mesh::EnsureNCMesh is called.
    MFEM_ABORT("Region-based refinement is currently only supported for simplex meshes!");
  }
  int region_ref_level = 0;
  bool use_nodes = (mesh.back()->GetNodes() != nullptr);
  int ref = use_nodes ? mesh.back()->GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  int dim = mesh.back()->SpaceDimension();
  while (region_ref_level < max_region_ref_levels)
  {
    // Mark elements for refinement in all regions. An element is marked for refinement if
    // any of its vertices are inside any refinement region for the given level.
    mfem::Array<mfem::Refinement> refs;
    for (int i = 0; i < mesh.back()->GetNE(); i++)
    {
      bool refine = false;
      mfem::DenseMatrix pointmat;
      if (use_nodes)
      {
        mfem::ElementTransformation *T = mesh.back()->GetElementTransformation(i);
        mfem::Geometry::Type geo = mesh.back()->GetElementGeometry(i);
        mfem::RefinedGeometry *RefG = mfem::GlobGeometryRefiner.Refine(geo, ref);
        T->Transform(RefG->RefPts, pointmat);
      }
      else
      {
        mfem::Array<int> verts;
        mesh.back()->GetElementVertices(i, verts);
        pointmat.SetSize(dim, verts.Size());
        for (int j = 0; j < verts.Size(); j++)
        {
          const double *coord = mesh.back()->GetVertex(verts[j]);
          for (int d = 0; d < dim; d++)
          {
            pointmat(d, j) = coord[d];
          }
        }
      }
      for (const auto &box : iodata.model.refinement.GetBoxes())
      {
        if (region_ref_level < box.ref_levels)
        {
          for (int j = 0; j < pointmat.Width(); j++)
          {
            // Check if the point is inside the box.
            int d = 0;
            for (; d < pointmat.Height(); d++)
            {
              if (pointmat(d, j) < box.bbmin[d] || pointmat(d, j) > box.bbmax[d])
              {
                break;
              }
            }
            if (d == dim)
            {
              refine = true;
              break;
            }
          }
          if (refine)
          {
            break;
          }
        }
      }
      if (refine)
      {
        refs.Append(mfem::Refinement(i));
        continue;
      }
      for (const auto &sphere : iodata.model.refinement.GetSpheres())
      {
        if (region_ref_level < sphere.ref_levels)
        {
          for (int j = 0; j < pointmat.Width(); j++)
          {
            // Check if the point is inside the sphere.
            double dist = 0.0;
            for (int d = 0; d < pointmat.Height(); d++)
            {
              double s = pointmat(d, j) - sphere.center[d];
              dist += s * s;
            }
            if (dist <= sphere.r * sphere.r)
            {
              refine = true;
              break;
            }
          }
          if (refine)
          {
            break;
          }
        }
      }
      if (refine)
      {
        refs.Append(mfem::Refinement(i));
      }
    }

    // Do the refinement. For tensor element meshes, this may make the mesh nonconforming
    // (adds hanging nodes).
    if (mesh.capacity() > 1)
    {
      mesh.push_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
    }
    mesh.back()->GeneralRefinement(refs, -1);
    region_ref_level++;
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the mesh after refinement if there is a single mesh (doesn't work with
  // h-refinement geometric multigrid).
  if (iodata.model.reorient_tet && mesh.size() == 1)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    mesh[0]->ReorientTetMesh();
#pragma GCC diagnostic pop
  }

  // Print some mesh information.
  mfem::Vector bbmin, bbmax;
  mesh[0]->GetBoundingBox(bbmin, bbmax);
  const double Lc = iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0);
  Mpi::Print(mesh[0]->GetComm(), "\nMesh curvature order: {:d}\nMesh bounding box:\n",
             mesh[0]->GetNodes() ? mesh[0]->GetNodes()->FESpace()->GetMaxElementOrder()
                                 : 1);
  if (mesh[0]->SpaceDimension() == 3)
  {
    Mpi::Print(mesh[0]->GetComm(),
               " (Xmin, Ymin, Zmin) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n"
               " (Xmax, Ymax, Zmax) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n",
               bbmin[0] * Lc, bbmin[1] * Lc, bbmin[2] * Lc, bbmax[0] * Lc, bbmax[1] * Lc,
               bbmax[2] * Lc);
  }
  else
  {
    Mpi::Print(mesh[0]->GetComm(),
               " (Xmin, Ymin) = ({:+.3e}, {:+.3e}) m\n"
               " (Xmax, Ymax) = ({:+.3e}, {:+.3e}) m\n",
               bbmin[0] * Lc, bbmin[1] * Lc, bbmax[0] * Lc, bbmax[1] * Lc);
  }
  Mpi::Print(mesh[0]->GetComm(), "\n{}", (mesh.size() > 1) ? "Coarse " : "");
  mesh[0]->PrintInfo();
  if (mesh.size() > 1)
  {
    Mpi::Print(mesh[0]->GetComm(), "\nRefined ");
    mesh.back()->PrintInfo();
  }
}

void AttrToMarker(int max_attr, const mfem::Array<int> &attrs, mfem::Array<int> &marker)
{
  MFEM_VERIFY(attrs.Size() == 0 || attrs.Max() <= max_attr,
              "Invalid attribute number present (" << attrs.Max() << ")!");
  marker.SetSize(max_attr);
  if (attrs.Size() == 1 && attrs[0] == -1)
  {
    marker = 1;
  }
  else
  {
    marker = 0;
    for (auto attr : attrs)
    {
      MFEM_VERIFY(attr > 0, "Attribute number less than one!");
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeate attribute in attribute list!");
      marker[attr - 1] = 1;
    }
  }
}

void AttrToMarker(int max_attr, const std::vector<int> &attrs, mfem::Array<int> &marker)
{
  MFEM_VERIFY(attrs.empty() || *std::max_element(attrs.begin(), attrs.end()) <= max_attr,
              "Invalid attribute number present ("
                  << *std::max_element(attrs.begin(), attrs.end()) << ")!");
  marker.SetSize(max_attr);
  if (attrs.size() == 1 && attrs[0] == -1)
  {
    marker = 1;
  }
  else
  {
    marker = 0;
    for (auto attr : attrs)
    {
      MFEM_VERIFY(attr > 0, "Attribute number less than one!");
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeate attribute in attribute list!");
      marker[attr - 1] = 1;
    }
  }
}

void GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                    mfem::Vector &max)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetBoundingBox(mesh, marker, bdr, min, max);
}

void GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                    mfem::Vector &min, mfem::Vector &max)
{
  int dim = mesh.SpaceDimension();
  min.SetSize(dim);
  max.SetSize(dim);
  for (int d = 0; d < dim; d++)
  {
    min(d) = mfem::infinity();
    max(d) = -mfem::infinity();
  }
  if (mesh.GetNodes() == nullptr)
  {
    auto BBUpdate = [&mesh, &dim, &min, &max](mfem::Array<int> &verts) -> void
    {
      for (int j = 0; j < verts.Size(); j++)
      {
        const double *coord = mesh.GetVertex(verts[j]);
        for (int d = 0; d < dim; d++)
        {
          if (coord[d] < min(d))
          {
            min(d) = coord[d];
          }
          if (coord[d] > max(d))
          {
            max(d) = coord[d];
          }
        }
      }
    };
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mfem::Array<int> verts;
        mesh.GetBdrElementVertices(i, verts);
        BBUpdate(verts);
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mfem::Array<int> verts;
        mesh.GetElementVertices(i, verts);
        BBUpdate(verts);
      }
    }
  }
  else
  {
    int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    auto BBUpdate = [&ref, &min, &max](mfem::ElementTransformation *T,
                                       mfem::Geometry::Type &geo) -> void
    {
      mfem::DenseMatrix pointmat;
      mfem::RefinedGeometry *RefG = mfem::GlobGeometryRefiner.Refine(geo, ref);
      T->Transform(RefG->RefPts, pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        for (int d = 0; d < pointmat.Height(); d++)
        {
          if (pointmat(d, j) < min(d))
          {
            min(d) = pointmat(d, j);
          }
          if (pointmat(d, j) > max(d))
          {
            max(d) = pointmat(d, j);
          }
        }
      }
    };
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
        mfem::Geometry::Type geo = mesh.GetBdrElementGeometry(i);
        BBUpdate(T, geo);
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
        mfem::Geometry::Type geo = mesh.GetElementGeometry(i);
        BBUpdate(T, geo);
      }
    }
  }
  auto *Min = min.HostReadWrite();
  auto *Max = max.HostReadWrite();
  Mpi::GlobalMin(dim, Min, mesh.GetComm());
  Mpi::GlobalMax(dim, Max, mesh.GetComm());
}

void GetSurfaceNormal(mfem::ParMesh &mesh, int attr, mfem::Vector &normal)
{
  mfem::Array<int> marker(mesh.bdr_attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetSurfaceNormal(mesh, marker, normal);
}

void GetSurfaceNormal(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                      mfem::Vector &normal)
{
  int dim = mesh.SpaceDimension();
  mfem::Vector nor(dim);
  normal.SetSize(dim);
  normal = 0.0;
  bool init = false;
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    if (!marker[mesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
    const mfem::IntegrationPoint &ip =
        mfem::Geometries.GetCenter(mesh.GetBdrElementGeometry(i));
    T->SetIntPoint(&ip);
    mfem::CalcOrtho(T->Jacobian(), nor);
    if (!init)
    {
      normal = nor;
      init = true;
    }
    else
    {
      // Check orientation and make sure consistent on this process. If a boundary has
      // conflicting normal definitions, use the first value.
      if (nor * normal < 0.0)
      {
        normal -= nor;
      }
      else
      {
        normal += nor;
      }
    }
  }
  // If different processors have different normal orientations, take that from the lowest
  // rank processor.
  MPI_Comm comm = mesh.GetComm();
  int rank = Mpi::Size(comm);
  mfem::Vector glob_normal(dim);
  if (init)
  {
    rank = Mpi::Rank(comm);
  }
  Mpi::GlobalMin(1, &rank, comm);
  if (rank == Mpi::Size(comm))
  {
    // No boundary elements of attribute attr.
    normal = 0.0;
    return;
  }
  if (rank == Mpi::Rank(comm))
  {
    glob_normal = normal;
  }
  {
    auto *GlobNormal = glob_normal.HostReadWrite();
    Mpi::Broadcast(dim, GlobNormal, rank, comm);
  }
  if (init && normal * glob_normal < 0.0)
  {
    normal.Neg();
  }
  {
    auto *Normal = normal.HostReadWrite();
    Mpi::GlobalSum(dim, Normal, comm);
  }
  normal /= normal.Norml2();
  // if (dim == 3)
  // {
  //   Mpi::Print(comm, " Surface normal {:d} = ({:+.3e}, {:+.3e}, {:+.3e})", attr,
  //   normal(0),
  //              normal(1), normal(2));
  // }
  // else
  // {
  //   Mpi::Print(comm, " Surface normal {:d} = ({:+.3e}, {:+.3e})", attr, normal(0),
  //              normal(1));
  // }
}

}  // namespace mesh

namespace
{

mfem::Mesh LoadMesh(const std::string &path)
{
  // Read the (serial) mesh from the given mesh file. Handle preparation for refinement and
  // orientations here to avoid possible reorientations and reordering later on. MFEM
  // supports a native mesh format (.mesh), VTK/VTU, Gmsh, as well as some others. We use
  // built-in converters for the types we know, otherwise rely on MFEM to do the conversion
  // or error out if not supported.
  std::filesystem::path mfile(path);
  if (mfile.extension() == ".mphtxt" || mfile.extension() == ".mphbin" ||
      mfile.extension() == ".nas" || mfile.extension() == ".bdf")
  {
    // Put translated mesh in temporary string buffer.
    std::stringstream fi(std::stringstream::in | std::stringstream::out);
    // fi << std::fixed;
    fi << std::scientific;
    fi.precision(MSH_FLT_PRECISION);

#if 0
    // Put translated mesh in temporary storage (directory is created and destroyed in
    // calling function).
    std::string tmp = iodata.problem.output;
    if (tmp.back() != '/')
    {
      tmp += '/';
    }
    tmp += "tmp/serial.msh";
    std::ofstream fo(tmp);
    // mfem::ofgzstream fo(tmp, true);  // Use zlib compression if available
    // fo << std::fixed;
    fo << std::scientific;
    fo.precision(MSH_FLT_PRECISION);
#endif

    if (mfile.extension() == ".mphtxt" || mfile.extension() == ".mphbin")
    {
      mesh::ConvertMeshComsol(path, fi);
      // mesh::ConvertMeshComsol(path, fo);
    }
    else
    {
      mesh::ConvertMeshNastran(path, fi);
      // mesh::ConvertMeshNastran(path, fo);
    }

#if 0
    std::ifstream fi(tmp);
    // mfem::ifgzstream fi(tmp);
    if (!fi.good())
    {
      MFEM_ABORT("Unable to open translated mesh file \"" << tmp << "\"!");
    }
#endif
    return mfem::Mesh(fi, 1, 1, true);
  }
  // Otherwise, just rely on MFEM load the mesh.
  std::ifstream fi(path);
  if (!fi.good())
  {
    MFEM_ABORT("Unable to open mesh file \"" << path << "\"!");
  }
  auto mesh = mfem::Mesh(fi, 1, 1, true);
  mesh.EnsureNodes();
  return mesh;
}

void ReorderMesh(mfem::Mesh &mesh)
{
  mfem::Array<int> ordering;

  // (Faster) Hilbert reordering.
  mesh.GetHilbertElementOrdering(ordering);
  mesh.ReorderElements(ordering);
}

std::unique_ptr<int[]> GetMeshPartitioning(mfem::Mesh &mesh, int size,
                                           const std::string &partition)
{
  MFEM_VERIFY(size <= mesh.GetNE(), "Mesh partitioning must have parts <= mesh elements ("
                                        << size << " vs. " << mesh.GetNE() << ")!");
  if (partition.length() == 0)
  {
    const int part_method = 1;
    std::unique_ptr<int[]> partitioning(mesh.GeneratePartitioning(size, part_method));
    Mpi::Print("Finished partitioning mesh into {:d} subdomain{}\n", size,
               (size > 1) ? "s" : "");
    return partitioning;
  }
  // User can optionally specify a mesh partitioning file as generated from the MFEM
  // mesh-explorer miniapp, for example. It has the format:
  //
  //   number_of_elements <NE>
  //   number_of_processors <NPART>
  //   <part[0]>
  //     ...
  //   <part[NE-1]>
  //
  int nel, np;
  std::ifstream part_ifs(partition);
  part_ifs.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
  part_ifs >> nel;
  if (nel != mesh.GetNE())
  {
    MFEM_ABORT("Invalid partitioning file (number of elements)!");
  }
  part_ifs.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
  part_ifs >> np;
  if (np != size)
  {
    MFEM_ABORT("Invalid partitioning file (number of processors)!");
  }
  auto partitioning = std::make_unique<int[]>(mesh.GetNE());
  int i = 0;
  while (i < mesh.GetNE())
  {
    part_ifs >> partitioning[i++];
  }
  Mpi::Print("Read mesh partitioning into {:d} subdomain{} from disk\n", size,
             (size > 1) ? "s" : "");
  return partitioning;
}

std::map<int, std::array<int, 2>> CheckMesh(mfem::Mesh &orig_mesh,
                                            const std::unique_ptr<int[]> &partitioning,
                                            const IoData &iodata, bool clean_elem,
                                            bool add_bdr, bool add_subdomain)
{
  // - Check that all external boundaries of the mesh have a corresponding boundary
  //   condition.
  // - If desired, create a new mesh which has added boundary elements for all material
  //   interfaces if these elements do not yet exist.
  // - If desired, create a new mesh which has removed all domain elements which do not have
  //   an associated material property specified in the input file.
  MFEM_VERIFY(orig_mesh.Dimension() == 3 && !orig_mesh.Nonconforming(),
              "Nonconforming or 2D meshes have not been tested yet!");
  mfem::Array<int> mat_marker, bdr_marker;
  GetUsedAttributeMarkers(iodata, orig_mesh.attributes.Max(),
                          orig_mesh.bdr_attributes.Max(), mat_marker, bdr_marker);
  bool warn = false;
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    int attr = orig_mesh.GetBdrAttribute(be);
    if (!bdr_marker[attr - 1])
    {
      int f, o, e1, e2;
      orig_mesh.GetBdrElementFace(be, &f, &o);
      orig_mesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)  // Internal boundary elements are allowed to have no BC
      {
        warn = true;
        break;
      }
    }
  }
  if (warn)
  {
    Mpi::Warning("One or more external boundary attributes has no associated boundary "
                 "condition!\n\"PMC\"/\"ZeroCharge\" condition is assumed!\n");
  }

  // Mapping from new interface boundary attribute tags to vector of neighboring domain
  // attributes (when adding new boundary elements).
  std::map<int, std::array<int, 2>> new_attr_map;
  if (!clean_elem && !add_bdr && !add_subdomain)
  {
    return new_attr_map;
  }

  // Count deleted or added domain and boundary elements.
  int new_ne = orig_mesh.GetNE();
  int new_nbdr = orig_mesh.GetNBE();
  mfem::Array<bool> elem_delete, bdr_delete;
  mfem::Array<int> orig_bdr_faces, add_bdr_faces;
  elem_delete.SetSize(orig_mesh.GetNE(), false);
  bdr_delete.SetSize(orig_mesh.GetNBE(), false);
  orig_bdr_faces.SetSize(orig_mesh.GetNumFaces(), -1);
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    int f, o;
    orig_mesh.GetBdrElementFace(be, &f, &o);
    MFEM_VERIFY(orig_bdr_faces[f] < 0,
                "Mesh should not define boundary elements multiple times!");
    orig_bdr_faces[f] = be;
  }
  if (add_bdr || add_subdomain)
  {
    add_bdr_faces.SetSize(orig_mesh.GetNumFaces(), -1);
  }

  if (clean_elem)
  {
    // Delete domain and boundary elements which have no associated material or BC attribute
    // from the mesh.
    for (int e = 0; e < orig_mesh.GetNE(); e++)
    {
      int attr = orig_mesh.GetAttribute(e);
      if (!mat_marker[attr - 1])
      {
        elem_delete[e] = true;
        new_ne--;
      }
    }

    // Make sure to remove any boundary elements which are no longer attached to elements in
    // the domain.
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be >= 0)
      {
        int e1, e2;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        if ((e1 < 0 || elem_delete[e1]) && (e2 < 0 || elem_delete[e2]))
        {
          // Mpi::Print("Deleting an unattached boundary element!\n");
          bdr_delete[be] = true;
          new_nbdr--;
        }
      }
    }
    if (new_ne < orig_mesh.GetNE())
    {
      Mpi::Print("Removed {:d} unmarked domain elements from the mesh\n",
                 orig_mesh.GetNE() - new_ne);
    }
    if (new_nbdr < orig_mesh.GetNBE())
    {
      Mpi::Print("Removed {:d} unattached boundary elements from the mesh\n",
                 orig_mesh.GetNBE() - new_nbdr);
    }
  }
  int new_ne_step1 = new_ne;
  int new_nbdr_step1 = new_nbdr;

  if (add_bdr)
  {
    // Add new boundary elements at material interfaces or on the exterior boundary of the
    // simulation domain, if there is not already a boundary element present.
    MFEM_VERIFY(!orig_mesh.Nonconforming(), "Adding material interface boundary elements "
                                            "is not supported for nonconforming meshes!");
    int add_bdr_ext = 0, add_bdr_int = 0;
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        if (e1 < 0 || elem_delete[e1] || e2 < 0 || elem_delete[e2])
        {
          // Mpi::Print("Adding exterior boundary element!\n");
          add_bdr_faces[f] = 1;
          add_bdr_ext++;
        }
        else if (orig_mesh.GetAttribute(e1) != orig_mesh.GetAttribute(e2))
        {
          // Add new boundary element at material interface between two domains.
          // Mpi::Print("Adding material interface boundary element!\n");
          add_bdr_faces[f] = 1;
          add_bdr_int++;
        }
      }
    }
    new_nbdr += (add_bdr_ext + add_bdr_int);
    if (add_bdr_ext > 0)
    {
      Mpi::Print("Added {:d} boundary elements for exterior boundaries to the mesh\n",
                 add_bdr_ext);
    }
    if (add_bdr_int > 0)
    {
      Mpi::Print("Added {:d} boundary elements for material interfaces to the mesh\n",
                 add_bdr_int);
    }
  }
  int new_ne_step2 = new_ne;
  int new_nbdr_step2 = new_nbdr;

  if (add_subdomain)
  {
    // Add new boundary elements at interfaces between elements beloning to different
    // subdomains. This uses similar code to mfem::Mesh::PrintWithPartitioning.
    MFEM_VERIFY(partitioning, "Cannot add subdomain interface boundary elements without "
                              "supplied mesh partitioning!");
    MFEM_VERIFY(!orig_mesh.Nonconforming(), "Adding subdomain interface boundary elements "
                                            "is not supported for nonconforming meshes!");
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        if (e1 >= 0 && !elem_delete[e1] && e2 >= 0 && !elem_delete[e2] &&
            partitioning[e1] != partitioning[e2])
        {
          // Internal face is connected to two elements belonging to different subdomains
          // (this works for conforming meshes).
          add_bdr_faces[f] = 2;
          new_nbdr += 2;
        }
      }
      // else
      // {
      //   // This face is attached to a boundary element. We could define a new boundary
      //   // element with opposite orientation to ensure both subdomains in the distributed
      //   // ParMesh have the boundary element.
      // }
    }
    if (new_nbdr > new_nbdr_step2)
    {
      Mpi::Print("Added boundary elements for subdomain interfaces to the mesh\n",
                 new_nbdr - new_nbdr_step2);
    }
  }

  // Create the new mesh.
  if (new_ne == new_ne_step1 && new_ne_step1 == new_ne_step2 &&
      new_ne_step2 == orig_mesh.GetNE() && new_nbdr == new_nbdr_step1 &&
      new_nbdr_step1 == new_nbdr_step2 && new_nbdr_step2 == orig_mesh.GetNBE())
  {
    return new_attr_map;
  }
  auto new_mesh = mfem::Mesh(orig_mesh.Dimension(), orig_mesh.GetNV(), new_ne, new_nbdr,
                             orig_mesh.SpaceDimension());

  // Copy vertices and non-deleted domain and boundary elements.
  for (int v = 0; v < orig_mesh.GetNV(); v++)
  {
    new_mesh.AddVertex(orig_mesh.GetVertex(v));
  }
  for (int e = 0; e < orig_mesh.GetNE(); e++)
  {
    if (!elem_delete[e])
    {
      mfem::Element *ne = orig_mesh.GetElement(e)->Duplicate(&new_mesh);
      new_mesh.AddElement(ne);
    }
  }
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    if (!bdr_delete[be])
    {
      mfem::Element *ne = orig_mesh.GetBdrElement(be)->Duplicate(&new_mesh);
      new_mesh.AddBdrElement(ne);
    }
  }

  // Add new boundary elements.
  if (add_bdr || add_subdomain)
  {
    auto FlipVertices = [](mfem::Element *e)
    {
      mfem::Array<int> v;
      e->GetVertices(v);
      int start = 0, end = v.Size() - 1;
      while (start < end)
      {
        int t = v[start];
        v[start] = v[end];
        v[end] = t;
        start++;
        end--;
      }
      e->SetVertices(v.HostRead());
    };

    // 1-based, some boundary attributes may be empty since they were removed from the
    // original mesh, but to keep indices the same as config file we don't compact the
    // list.
    int max_bdr_attr = orig_mesh.bdr_attributes.Max();
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      if (add_bdr_faces[f] > 0)
      {
        // Assign new unique attribute based on attached elements (we want the material
        // properties on the face to average those on the elements). This is used later on
        // when integrating the transmission condition on the subdomain interface. Save the
        // inverse so that the attributes of e1 and e2 can be easily referenced using the
        // new attribute. Since attributes are in 1-based indexing, a, b > 0.
        int e1, e2, a = 0, b = 0;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        if (e1 >= 0 && !elem_delete[e1] && e2 >= 0 && !elem_delete[e2])
        {
          a = std::max(orig_mesh.GetAttribute(e1), orig_mesh.GetAttribute(e2));
          b = (a == orig_mesh.GetAttribute(e1)) ? orig_mesh.GetAttribute(e2)
                                                : orig_mesh.GetAttribute(e1);
        }
        else if (e1 >= 0 && !elem_delete[e1])
        {
          a = orig_mesh.GetAttribute(e1);
          b = 0;
        }
        else if (e2 >= 0 && !elem_delete[e2])
        {
          a = orig_mesh.GetAttribute(e2);
          b = 0;
        }
        MFEM_VERIFY(a + b > 0, "Invalid new boundary element attribute!");
        int new_attr = max_bdr_attr + (a * (a - 1)) / 2 + b;  // At least max_bdr_attr+1
        if (new_attr_map.find(new_attr) == new_attr_map.end())
        {
          new_attr_map.emplace(new_attr, std::array<int, 2>{a, b});
        }

        // Add the boundary elements with the new boundary attribute.
        mfem::Element *ne = orig_mesh.GetFace(f)->Duplicate(&new_mesh);
        ne->SetAttribute(new_attr);
        new_mesh.AddBdrElement(ne);
        if (add_bdr_faces[f] > 1)
        {
          // Flip order of vertices to reverse normal direction of second added element.
          ne = orig_mesh.GetFace(f)->Duplicate(&new_mesh);
          FlipVertices(ne);
          ne->SetAttribute(new_attr);
          new_mesh.AddBdrElement(ne);
          // Mpi::Print("Adding two BE with attr {:d} from elements {:d} and {:d}\n",
          //            new_attr, a, b);
        }
      }
    }
  }

  // Finalize new mesh and replace the old one. If a curved mesh, set up the new mesh by
  // projecting nodes onto the new mesh for the non-trimmed vdofs (accounts for new
  // boundary elements too since no new dofs are added). See the MFEM trimmer miniapp for
  // reference.
  new_mesh.FinalizeTopology();
  new_mesh.Finalize();
  new_mesh.RemoveUnusedVertices();
  if (orig_mesh.GetNodes())
  {
    const mfem::GridFunction *nodes = orig_mesh.GetNodes();
    const mfem::FiniteElementSpace *fes = nodes->FESpace();

    mfem::Ordering::Type ordering = fes->GetOrdering();
    int order = fes->GetMaxElementOrder();
    int sdim = orig_mesh.SpaceDimension();
    bool discont = dynamic_cast<const mfem::L2_FECollection *>(fes->FEColl()) != nullptr;

    new_mesh.SetCurvature(order, discont, sdim, ordering);
    mfem::GridFunction *new_nodes = new_mesh.GetNodes();
    const mfem::FiniteElementSpace *new_fes = new_nodes->FESpace();

    // The element loop works because we know the mapping from old_mesh to new_mesh element
    // indices from the insertion order.
    mfem::Array<int> vdofs, new_vdofs;
    mfem::Vector loc_vec;
    int te = 0;
    for (int e = 0; e < orig_mesh.GetNE(); e++)
    {
      if (!elem_delete[e])
      {
        fes->GetElementVDofs(e, vdofs);
        nodes->GetSubVector(vdofs, loc_vec);
        new_fes->GetElementVDofs(te, new_vdofs);
        new_nodes->SetSubVector(new_vdofs, loc_vec);
        te++;
      }
    }
  }
  orig_mesh = std::move(new_mesh);
  return new_attr_map;
}

mfem::ParMesh DistributeMesh(MPI_Comm comm, mfem::Mesh &mesh,
                             const std::unique_ptr<int[]> &partitioning)
{
  return mfem::ParMesh(comm, mesh, partitioning.get());
}

void GetUsedAttributeMarkers(const IoData &iodata, int n_mat, int n_bdr,
                             mfem::Array<int> &mat_marker, mfem::Array<int> &bdr_marker)
{
  mfem::Array<int> mat_attr, bdr_attr;
  mat_attr.Reserve(static_cast<int>(iodata.domains.attributes.size()));
  for (auto attr : iodata.domains.attributes)
  {
    mat_attr.Append(attr);
  }
  bdr_attr.Reserve(static_cast<int>(iodata.boundaries.attributes.size()));
  for (auto attr : iodata.boundaries.attributes)
  {
    bdr_attr.Append(attr);
  }
  mesh::AttrToMarker(n_mat, mat_attr, mat_marker);
  mesh::AttrToMarker(n_bdr, bdr_attr, bdr_marker);
}

}  // namespace

}  // namespace palace
