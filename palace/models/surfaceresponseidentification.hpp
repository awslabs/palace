// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_IDENTIFICATION_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_IDENTIFICATION_HPP

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "utils/labels.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

//
// Geometry identification for the surface-response correction: a pure function of the metal
// perimeter, the per-segment process frames and the matching radius R which produces the
// canonical per-segment contract of SURFACE-RESPONSE-IDENTIFICATION.md: every perimeter
// segment is assigned to feature portions or to exactly one exclusion record, every corner
// / endpoint / junction vertex to exactly one vertex feature or cluster, and every feature
// carries a translation-, rotation- and mirror-invariant signature with a stable hash. The
// library is not an input; matching is a separate lookup by signature.
//

struct IdentificationSegment
{
  std::array<double, 3> p0{};
  std::array<double, 3> p1{};
  std::array<std::size_t, 2> vertices{};
  int chain = -1;
  bool truncation = false;
  int conductor = 0;
  std::map<InterfaceDielectric, int> targets;
  std::array<double, 3> gap_direction{};
  std::array<double, 3> process_normal{};
  std::string boundary_law;

  // A segment excluded before identification (no target interface, ...) carries its class
  // and reason and takes no part in the geometry.
  std::optional<std::pair<std::string, std::string>> exclusion;
};

struct IdentificationVertex
{
  std::array<double, 3> coordinate{};
  std::vector<std::size_t> segments;
  std::optional<MetalEdgeVertexType> physical_type;
  bool on_truncation_boundary = false;
};

struct IdentificationInput
{
  double radius = 0.0;
  std::vector<IdentificationSegment> segments;
  std::vector<IdentificationVertex> vertices;
};

struct IdentifiedPortion
{
  std::size_t segment = 0;
  double s0 = 0.0;
  double s1 = 0.0;
};

struct IdentifiedFeature
{
  int id = 0;
  std::string type;
  nlohmann::json signature;
  std::string signature_key;
  std::string hash;
  int chirality = 1;
  double length = 0.0;
  std::vector<IdentifiedPortion> portions;
  std::vector<std::size_t> vertices;
  std::array<double, 3> origin{};
  std::array<std::array<double, 3>, 3> axes{};

  // Filled by the matching pass.
  std::optional<std::string> matched_model;
};

struct IdentifiedSegment
{
  // Canonical key: the two endpoints in lexicographic order (portions run from key[0]).
  std::array<std::array<double, 3>, 2> key{};
  double length = 0.0;
  int chain = -1;
  std::vector<std::array<double, 3>> portions;  // {s0, s1, feature id}
  std::optional<std::pair<std::string, std::string>> exclusion;
};

struct IdentifiedVertex
{
  std::size_t vertex = 0;
  std::string type;  // Corner | Endpoint | Junction | TruncationCut
  double turn_degrees = 0.0;
  int feature = -1;
};

struct IdentificationExclusion
{
  std::string cls;
  std::string reason;
  int count = 0;
  double length = 0.0;
};

struct IdentificationResult
{
  double radius = 0.0;
  std::array<double, 3> reference_process_normal{};
  std::vector<IdentifiedFeature> features;
  std::vector<IdentifiedSegment> segments;
  std::vector<IdentifiedVertex> vertices;
  std::vector<IdentificationExclusion> exclusions;
  double perimeter_length = 0.0;
  double assigned_length = 0.0;
  double excluded_length = 0.0;
  std::string geometry_digest;

  // Manifest "Identification" object; the length scale converts mesh units for output.
  nlohmann::json ToJson(double length_scale) const;
};

IdentificationResult IdentifyMetalPerimeter(const IdentificationInput &input);

// Canonical signature of a set of straight edge portions and vertices in a common frame,
// shared by the device features and the library models so that both sides are hashed by
// the same function. Every portion is {p0, p1, gap direction, process normal, conductor,
// interface types, law}; the result is the minimal serialisation over the candidate frames.
struct SignaturePortion
{
  std::array<double, 3> p0{};
  std::array<double, 3> p1{};
  std::array<double, 3> gap_direction{};
  int conductor = 0;
  std::vector<std::string> interfaces;
  std::string boundary_law;
};

struct SignatureVertex
{
  std::array<double, 3> point{};
  std::string type;
  double turn_degrees = 0.0;
};

struct CanonicalSignature
{
  nlohmann::json signature;
  std::string key;
  std::string hash;
  int chirality = 1;
  std::array<double, 3> origin{};
  std::array<std::array<double, 3>, 3> axes{};
};

CanonicalSignature CanonicalClusterSignature(const std::vector<SignaturePortion> &portions,
                                             const std::vector<SignatureVertex> &vertices,
                                             const std::array<double, 3> &process_normal,
                                             double radius);

// Canonical signature of parallel edges over a common longitudinal interval: offsets / R
// from the lowest edge, gap side (+1 toward increasing offset), conductor labels by first
// appearance, interface types and law; minimal over the two lateral orientations (mirror).
struct TranslationalEdge
{
  double offset = 0.0;
  int gap_sign = 1;
  int conductor = 0;
  std::vector<std::string> interfaces;
  std::string boundary_law;
};

struct TranslationalSignature
{
  nlohmann::json signature;
  int chirality = 1;
};

TranslationalSignature CanonicalTranslationalSignature(std::vector<TranslationalEdge> edges,
                                                       double radius);

// Signature grids. Signature coordinates (portion endpoints / R, offsets / R, corner
// radii / R) come from bisections at analytic region boundaries and from chip-scale mesh
// coordinates whose roundoff is ~ulp(|p|) (1e-12 um at 10 mm); the grid must be far above
// that roundoff so that translated / rotated copies of one feature hash identically, and far
// below any resolution the response can depend on (the response varies on the scale R).
// 1e-6 R (2 pm at R = 2 um) satisfies both; angles use the same relative grid in degrees.
constexpr double kSignatureLengthQuantumOverRadius = 1.0e-6;
constexpr double kSignatureAngleQuantumDegrees = 1.0e-6;

// Corner / junction signatures (without "Type"), shared by device features and library
// models so that one canonicalisation produces both keys.
nlohmann::json CanonicalCornerSignature(const std::vector<std::string> &interfaces,
                                        const std::string &boundary_law,
                                        double angle_degrees, double corner_radius_over_R);
nlohmann::json CanonicalJunctionSignature(const std::vector<std::string> &interfaces,
                                          const std::string &boundary_law,
                                          std::vector<double> arm_angles_degrees);

// Feature signature key and hash shared by device features and library models.
std::pair<std::string, std::string> SignatureKeyAndHash(nlohmann::json signature,
                                                        const std::string &type);

std::string Sha256Hex(const std::string &text);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_IDENTIFICATION_HPP
