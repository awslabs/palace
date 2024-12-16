// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_EXCITATIONS_HPP
#define PALACE_UTILS_EXCITATIONS_HPP

#include <cmath>

namespace palace::excitations
{

//
// Define temporal excitation functions for transient simulations.
//

inline mfem::real_t pulse_sinusoidal(mfem::real_t t, mfem::real_t omega, mfem::real_t t0)
{
  // g(t) = sin(ω*(t-t0))
  return std::sin(omega * (t - t0));
}

inline mfem::real_t dpulse_sinusoidal(mfem::real_t t, mfem::real_t omega, mfem::real_t t0)
{
  // g(t) = sin(ω*(t-t0))
  return omega * std::cos(omega * (t - t0));
}

inline mfem::real_t pulse_gaussian(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ts = t - t0;
  return std::exp(-0.5 * ts * ts / (tau * tau));
}

inline mfem::real_t dpulse_gaussian(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ootau2 = 1.0 / (tau * tau);
  mfem::real_t ts = t - t0;
  return -ts * ootau2 * std::exp(-0.5 * ts * ts * ootau2);
}

inline mfem::real_t pulse_gaussian_diff(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = -(t-t0)/τ²*exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ootau2 = 1.0 / (tau * tau);
  mfem::real_t ts = t - t0;
  return -ts * ootau2 * std::exp(-0.5 * ts * ts * ootau2);
}

inline mfem::real_t dpulse_gaussian_diff(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = -(t-t0)/τ²*exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ootau2 = 1.0 / (tau * tau);
  mfem::real_t ts = t - t0;
  mfem::real_t ts2 = ts * ts;
  return -ootau2 * (1.0 - ts2 * ootau2) * std::exp(-0.5 * ts2 * ootau2);
}

inline mfem::real_t pulse_gaussian_mod(mfem::real_t t, mfem::real_t omega, mfem::real_t tau,
                                       mfem::real_t t0)
{
  // g(t) = sin(ω*(t-t0))*exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ts = t - t0;
  return std::sin(omega * ts) * std::exp(-0.5 * ts * ts / (tau * tau));
}

inline mfem::real_t dpulse_gaussian_mod(mfem::real_t t, mfem::real_t omega,
                                        mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = sin(ω*(t-t0))*exp(-0.5*(t-t0)²/τ²)
  mfem::real_t ootau2 = 1.0 / (tau * tau);
  mfem::real_t ts = t - t0;
  return (-ts * ootau2 * std::sin(omega * ts) + omega * std::cos(omega * ts)) *
         std::exp(-0.5 * ts * ts * ootau2);
}

inline mfem::real_t pulse_ramp(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = 0, t <= t0
  //        (t-t0)/τ, t0 < t <= τ
  //        1, t > τ+t0
  return (t <= t0) ? 0.0 : ((t - t0 >= tau) ? 1.0 : (t - t0) / tau);
}

inline mfem::real_t dpulse_ramp(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = 0, t <= t0
  //        (t-t0)/τ, t0 < t <= τ
  //        1, t > τ
  return (t <= t0) ? 0.0 : ((t - t0 >= tau) ? 0.0 : 1.0 / tau);
}

inline mfem::real_t pulse_smootherstep(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = 0, t <= t0
  //        6*((t-t0)/τ)⁵-15*((t-t0)/τ)⁴+10*((t-t0)/τ)³, t0 < t <= τ+t0
  //        1, t > τ+t0
  mfem::real_t ts = (t <= t0) ? 0.0 : ((t - t0 >= tau) ? 1.0 : (t - t0) / tau);
  mfem::real_t ts2 = ts * ts;
  return ts * ts2 * (6.0 * ts2 - 15.0 * ts + 10.0);
}

inline mfem::real_t dpulse_smootherstep(mfem::real_t t, mfem::real_t tau, mfem::real_t t0)
{
  // g(t) = 0, t <= t0
  //        6*((t-t0)/τ)⁵-15*((t-t0)/τ)⁴+10*((t-t0)/τ)³, t0 < t <= τ
  //        1, t > τ
  mfem::real_t ts = (t <= t0) ? 0.0 : ((t - t0 >= tau) ? 1.0 : (t - t0) / tau);
  mfem::real_t ts2 = ts * ts;
  return ts2 / tau * (30.0 * ts2 - 60.0 * ts + 30.0);
}

}  // namespace palace::excitations

#endif  // PALACE_UTILS_EXCITATIONS_HPP
