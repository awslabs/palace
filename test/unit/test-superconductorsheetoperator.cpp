// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/superconductorsheetoperator.hpp"

using namespace palace;

// The finite-thickness London kinetic sheet inductance L_ksq = lambda * coth(d/lambda) is
// the experimentally validated form (Tolpygo et al., Supercond. Sci. Technol. 34, 085005
// (2021), for Nb films down to 120 nm). It reduces to the thin-film Pearl limit lambda^2/d
// for d << lambda and saturates at lambda for d >> lambda; the old sheet used the thin-film
// expression, which undercounts the kinetic inductance by (d/lambda)*coth(d/lambda) once
// d >~ lambda (e.g. ~2.3x at the MIT-LL SFQ5ee geometry d=200 nm, lambda=90 nm).
TEST_CASE("Superconductor sheet kinetic inductance", "[superconductor][Serial]")
{
  auto Lksq = [](double lam, double d)
  { return SuperconductorSheetOperator::KineticSheetInductance(lam, d); };

  SECTION("Matches lambda*coth(d/lambda)")
  {
    for (double lam : {0.05, 0.09, 0.4, 1.0})
    {
      for (double d : {0.02, 0.1, 0.2, 0.5, 2.0})
      {
        CHECK_THAT(Lksq(lam, d),
                   Catch::Matchers::WithinRel(lam / std::tanh(d / lam), 1.0e-12));
      }
    }
  }

  SECTION("Thin-film limit d << lambda -> lambda^2/d")
  {
    // coth(x) -> 1/x + x/3, so lambda*coth(d/lambda) -> lambda^2/d to O((d/lambda)^2).
    const double lam = 0.4, d = 1.0e-3;  // d/lambda = 2.5e-3
    CHECK_THAT(Lksq(lam, d), Catch::Matchers::WithinRel(lam * lam / d, 1.0e-5));
  }

  SECTION("Thick-film limit d >> lambda -> lambda")
  {
    const double lam = 0.09, d = 10.0;  // d/lambda ~ 111
    CHECK_THAT(Lksq(lam, d), Catch::Matchers::WithinRel(lam, 1.0e-6));
  }

  SECTION("Finite-thickness enhancement over thin-film at MIT-LL geometry")
  {
    // d = 200 nm, lambda = 90 nm (d/lambda = 2.222): the finite-thickness value exceeds the
    // thin-film lambda^2/d by the factor (d/lambda)*coth(d/lambda), ~2.28x here.
    const double lam = 0.09, d = 0.2;
    const double thin = lam * lam / d;
    const double factor = Lksq(lam, d) / thin;
    CHECK_THAT(factor, Catch::Matchers::WithinRel((d / lam) / std::tanh(d / lam), 1.0e-12));
    CHECK(factor > 2.2);
    CHECK(factor < 2.4);
  }
}
