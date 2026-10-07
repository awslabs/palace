// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_MIRROR_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_MIRROR_HPP

#include <array>
#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "models/surfaceresponseidentification.hpp"

namespace mfem
{

class ParMesh;

}  // namespace mfem

namespace palace
{

//
// Mirror-extended identification across planar NATURAL truncation planes (boundary-cut
// DESIGN 2.2, decisions 442 / 454 / 455; rule F-DB-c). On a planar natural (homogeneous
// Neumann) face the half-domain solution is the restriction of the mirror-symmetric full
// problem, so the device's truth near a window cut, a symmetry plane or a natural outline
// is the MIRRORED geometry's: an edge meeting the plane at theta is half of a 2 theta
// corner, an edge parallel at d < R half of a 2 d strip / gap, a straight perpendicular
// meeting forms no new feature. The band of the metal perimeter within kMirrorBandOverR x R
// of each such plane is reflected and appended to the identification input as IMAGE
// segments (joined to the real chain at the truncation vertex on the plane), the
// identification runs unchanged on the extended input, and the result is MERGED onto the
// unextended run so that the real features keep their Ids, portions and order wherever the
// band formed no new feature (the bitwise requirements 2.2.2 (a)-(c)).
//

// A truncation plane: n . x = c with n the OUTWARD unit normal of the domain, fitted from
// the boundary faces of one truncation attribute (an exterior, nonmetal, non-interface
// attribute under no boundary condition: metaledge's simulation-cut surfaces). Status:
// "Natural" (planar, vertical: mirrors), "Unsupported" (not parallel to the process normal:
// a z cut would flip the metal's process side), "NonPlanar" (the attribute's faces are
// chords of a curved surface: two of its plane groups meet at a dihedral turn below the
// arc rule's joint cap, kArcMaxJointTurnDegrees).
struct MirrorPlane
{
  int attribute = 0;
  std::array<double, 3> normal{};
  double offset = 0.0;
  int faces = 0;
  double max_deviation = 0.0;
  std::string status;
  // The bounding box of the plane's faces: a plane acts only where a point's PROJECTION
  // onto it falls within its faces (a non-convex domain — a window with a re-entrant
  // outline, the symmetry fixture's chevron — has planes whose infinite extension passes
  // through the domain elsewhere).
  std::array<double, 3> box_min{}, box_max{};
  bool Mirrors() const { return status == "Natural"; }
  bool NearFaces(const std::array<double, 3> &p, double margin) const
  {
    const double d = Inside(p);
    for (int k = 0; k < 3; k++)
    {
      const double projected = p[k] + d * normal[k];
      if (projected < box_min[k] - margin || projected > box_max[k] + margin)
      {
        return false;
      }
    }
    return true;
  }
  // Signed distance of a point into the domain (>= 0 inside).
  double Inside(const std::array<double, 3> &p) const
  {
    return offset - (normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2]);
  }
  std::array<double, 3> Reflect(const std::array<double, 3> &p) const
  {
    const double d = Inside(p);
    return {p[0] + 2.0 * d * normal[0], p[1] + 2.0 * d * normal[1],
            p[2] + 2.0 * d * normal[2]};
  }
  // The linear part of the reflection on a direction.
  std::array<double, 3> ReflectDirection(const std::array<double, 3> &v) const
  {
    const double dot = normal[0] * v[0] + normal[1] * v[1] + normal[2] * v[2];
    return {v[0] - 2.0 * dot * normal[0], v[1] - 2.0 * dot * normal[1],
            v[2] - 2.0 * dot * normal[2]};
  }
};

// Width of the mirror band in matching radii: the 2R interaction reach of every
// identification rule plus R of margin, so that no real feature farther than 2R from a
// plane can change and every image artefact at the band's end lies beyond 2R of real
// perimeter (DESIGN 2.2.2).
constexpr double kMirrorBandOverRadius = 3.0;

// Fit the truncation planes of the given boundary attributes (collective over the mesh's
// communicator; identical on every rank: the faces are gathered and the planes sorted by
// (attribute, quantised normal, quantised offset)). `process_normal` decides Vertical.
std::vector<MirrorPlane> FitTruncationPlanes(const mfem::ParMesh &mesh,
                                             const std::set<int> &truncation_attributes,
                                             const std::array<double, 3> &process_normal,
                                             double matching_radius);

// Reflect a point outside the domain through the mirror planes it lies beyond, in canonical
// order (plane index ascending; the composition for a box edge / corner), as long as the
// point lies within `band` of every plane it crosses. Returns the image and the planes
// used, or nullopt when the point is beyond the band of a plane or lies beyond a plane that
// does not mirror; a point inside every plane is returned unchanged with no reflection. Any
// point STRICTLY beyond a plane (within its faces' projection, `tolerance`) is reflected: a
// sample a roundoff beyond a cut face is outside the mesh for the locator.
struct ReflectedPoint
{
  std::array<double, 3> point{};
  std::vector<int> planes;
};
std::optional<ReflectedPoint> ReflectIntoDomain(const std::array<double, 3> &p,
                                                const std::vector<MirrorPlane> &planes,
                                                double band, double tolerance);

// The mirror-band extension of an identification input (a pure function of (input,
// planes)): every physical segment with a point within band x R of a mirroring plane is
// reflected through every non-empty subset of the planes it is near (the group generated
// by the reflections: one image per plane, the double image at a box edge, ...) and
// appended with image_of / mirror_planes; image vertices are created per (real vertex,
// plane subset) except that a real vertex ON a plane is its own image (the joint); a real
// truncation vertex joined to its image is typed by the geometric joint noise rule
// (REGULAR: the image chain continues the real chain's id; CORNER: a new chain) and loses
// its truncation flag; an image vertex whose real counterpart has further physical segments
// not imaged is an ImageBandCut (ENDPOINT on the truncation boundary: never a feature).
// Image faces (for the CrossLayer rule) are the reflected faces within the band.
struct MirrorExtension
{
  IdentificationInput input;
  std::size_t real_segments = 0;
  std::size_t real_vertices = 0;
  int image_segments = 0;
  int image_vertices = 0;
  int image_faces = 0;
  std::vector<std::size_t> joined_vertices;  // real truncation vertices joined (ascending)
  std::vector<std::size_t> band_cut_vertices;
  // Per joined vertex: the planes it lies on and whether the joint is straight (REGULAR).
  std::vector<std::vector<int>> joined_planes;
  std::vector<bool> joined_straight;
};
MirrorExtension ExtendIdentificationInputAcrossMirrorPlanes(
    const IdentificationInput &input, const std::vector<MirrorPlane> &planes,
    double band_over_radius = kMirrorBandOverRadius);

// The start of the REAL arm's translational cells of a virtual (mirror-formed) corner of
// the given angle: s_half = (R + s) / 2 with s = R / max(|cos theta|, |sin theta|) the
// arm's exit from the coupon's matching square (DESIGN 2.2.3: the half-weight coupon
// integrates E(real [0, R]) + E(real [0, s]) over two, which equals E(real [0, s_half]) up
// to the second-order term (E[R, s_half] - E[s_half, s]) / 2; exactly R at 90 degrees).
double HalfCornerArmStart(double angle_degrees, double matching_radius);

// Merge the identification of the extended input onto the unextended run (DESIGN 2.2.2
// (a)-(c)): the real features, Ids, portions and order are the unextended run's wherever
// the images formed no new feature (a feature of the extended run whose real portions are
// exactly one real feature's portions with the same type and key: Mirror Status
// "Continued"); a feature of the extended run touching an image whose real portions are not
// such a feature is MIRROR-FORMED (the 2 theta corner, the 2 d strip / gap, a bent stack,
// a two-vertex cluster): the real features' portions it overlaps are clipped to it, a real
// feature it replaces wholly lends it its Id, otherwise it is numbered after every real
// feature; its portions keep the image portions (segment index >= real_segments, framed by
// the caller for the pair / stack placement, never placed, never counted); image-only
// features are counted and dropped. The merged result's segments are the real ones followed
// by the image ones; the joint vertices take the extended run's reading. `diagnostics`
// receives the MirrorBand record.
struct MirrorMergeSummary
{
  int continued_features = 0;
  int mirror_formed_features = 0;
  int image_only_features = 0;
  std::vector<int> mirror_formed_ids;  // feature ids in the merged result
  // Mirror-formed features of a topology without a mirror placement (stacks, clusters,
  // curved pairs): not merged, recorded {Feature, Type, Key, RealLength, ImageLength,
  // Planes, Status "Unmerged"} for the discovery.
  nlohmann::json unmerged_features = nlohmann::json::array();
};
MirrorMergeSummary MergeMirrorIdentification(const IdentificationResult &real,
                                             const IdentificationResult &extended,
                                             const MirrorExtension &extension,
                                             const std::vector<MirrorPlane> &planes,
                                             IdentificationResult &merged);

// The Identification.Diagnostics.MirrorBand record.
nlohmann::json DescribeMirrorBand(const std::vector<MirrorPlane> &planes,
                                  const MirrorExtension &extension,
                                  const MirrorMergeSummary &summary,
                                  const IdentificationResult &merged,
                                  double band_over_radius, double coordinate_scale);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_MIRROR_HPP
