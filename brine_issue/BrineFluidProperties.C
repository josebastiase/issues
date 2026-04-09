//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BrineFluidProperties.h"

#include <algorithm>
#include <cmath>

registerMooseObject("FluidPropertiesApp", BrineFluidProperties);

namespace
{
  // SimpleFluidProperties calls inRegion() and rejects p < 0.
  // In unsaturated (Richards) flow the liquid pressure can be negative (suction).
  // Clamp to 1 Pa — with a large bulk modulus density is essentially constant anyway.
  inline Real safePressure(Real p) { return p < 1.0 ? 1.0 : p; }

  constexpr Real T_FALLBACK = 293.15; // K — used when T is NaN
  inline Real safeT(Real T) { return std::isfinite(T) ? T : T_FALLBACK; }
} // namespace

InputParameters
BrineFluidProperties::validParams()
{
  InputParameters params = MultiComponentFluidProperties::validParams();
  params.addParam<UserObjectName>("water_fp",
                                  "The name of the FluidProperties UserObject for water");
  params.addClassDescription("Fluid properties for brine");
  return params;
}

BrineFluidProperties::BrineFluidProperties(const InputParameters & parameters)
  : MultiComponentFluidProperties(parameters), _water_fp_derivs(true)
{
  // We maintain a Water97FluidProperties object ONLY for Henry's constant formulation.
  // All other property calculations are done using _water_fp (which we set to SimpleFluidProperties
  // by default, unless user supplies a water_fp).
  const std::string water97_name = name() + ":water97";
  {
    const std::string class_name = "Water97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    if (_tid == 0)
      _fe_problem.addUserObject(class_name, water97_name, params);
  }
  _water97_fp = &_fe_problem.getUserObject<Water97FluidProperties>(water97_name);

  // Choose the water backend used for rho/mu/h/cp/k/etc.
  if (parameters.isParamSetByUser("water_fp"))
  {
    // User-supplied SinglePhaseFluidProperties for the "water part"
    _water_fp = &getUserObject<SinglePhaseFluidProperties>("water_fp");

    // NOTE:
    // We intentionally do NOT enforce _water_fp->fluidName() == "water" anymore,
    // because SimpleFluidProperties (and others) may use different fluidName() strings.
  }
  else
  {
    // Default: construct a SimpleFluidProperties UserObject for water-like behavior.
    const std::string water_name = name() + ":water";
    {
      const std::string class_name = "SimpleFluidProperties";
      InputParameters params = _app.getFactory().getValidParams(class_name);
      // Defaults from docs are generally reasonable; user can override by supplying water_fp.
      if (_tid == 0)
        _fe_problem.addUserObject(class_name, water_name, params);
    }
    _water_fp = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(water_name);
  }

  // SinglePhaseFluidProperties UserObject for NaCl (used mainly for molar mass)
  const std::string nacl_name = name() + ":nacl";
  {
    const std::string class_name = "NaClFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    if (_tid == 0)
      _fe_problem.addUserObject(class_name, nacl_name, params);
  }
  _nacl_fp = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(nacl_name);

  // Molar mass of NaCl and H2O
  _Mnacl = _nacl_fp->molarMass();
  _Mh2o = _water_fp->molarMass();
}

BrineFluidProperties::~BrineFluidProperties() {}

const SinglePhaseFluidProperties &
BrineFluidProperties::getComponent(unsigned int component) const
{
  switch (component)
  {
    case WATER:
      return *_water_fp;

    case NACL:
      return *_nacl_fp;

    default:
      mooseError("BrineFluidProperties::getComponent has been provided an incorrect component");
  }
}

std::string
BrineFluidProperties::fluidName() const
{
  return "brine";
}

FPADReal
BrineFluidProperties::molarMass(const FPADReal & xnacl) const
{
  return 1.0 / (xnacl / _Mnacl + (1.0 - xnacl) / _Mh2o);
}

Real
BrineFluidProperties::molarMass(Real xnacl) const
{
  return 1.0 / (xnacl / _Mnacl + (1.0 - xnacl) / _Mh2o);
}

Real
BrineFluidProperties::molarMassNaCl() const
{
  return _Mnacl;
}

Real
BrineFluidProperties::molarMassH2O() const
{
  return _Mh2o;
}

FPADReal
BrineFluidProperties::rho_from_p_T_X(const FPADReal & pressure,
                                     const FPADReal & temperature,
                                     const FPADReal & xnacl) const
{
  // NOTE: The Driesner virtual-temperature (Tv) trick is intentionally bypassed.
  // Tv maps brine state to an equivalent pure-water state tuned for Water97's
  // reference — it gives wrong results with SimpleFluidProperties.
  // We evaluate at the real (p,T) and scale by the NaCl molar mass ratio.
  const Real p_safe = safePressure(pressure.value());
  const Real T_safe = safeT(temperature.value());

  FPADReal water_density;
  if (_water_fp_derivs)
  {
    Real rho, drho_dp, drho_dT;
    _water_fp->rho_from_p_T(p_safe, T_safe, rho, drho_dp, drho_dT);
    water_density = rho;
    water_density.derivatives() = pressure.derivatives() * drho_dp +
                                   temperature.derivatives() * drho_dT;
  }
  else
    water_density = _water_fp->rho_from_p_T(p_safe, T_safe);

  return water_density * molarMass(xnacl) / _Mh2o;
}


Real
BrineFluidProperties::rho_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  const Real rho_w = _water_fp->rho_from_p_T(safePressure(pressure), safeT(temperature));
  return rho_w * molarMass(xnacl) / _Mh2o;
}

void
BrineFluidProperties::rho_from_p_T_X(Real pressure,
                                     Real temperature,
                                     Real xnacl,
                                     Real & rho,
                                     Real & drho_dp,
                                     Real & drho_dT,
                                     Real & drho_dx) const
{
  FPADReal p = safePressure(pressure);
  Moose::derivInsert(p.derivatives(), 0, pressure >= 1.0 ? 1.0 : 0.0);
  FPADReal T = safeT(temperature);
  Moose::derivInsert(T.derivatives(), 1, std::isfinite(temperature) ? 1.0 : 0.0);
  FPADReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPADReal ad_rho = this->rho_from_p_T_X(p, T, x);

  rho = ad_rho.value();
  drho_dp = ad_rho.derivatives()[0];
  drho_dT = ad_rho.derivatives()[1];
  drho_dx = ad_rho.derivatives()[2];
}

Real
BrineFluidProperties::mu_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  const Real p_eff = safePressure(pressure);
  const Real T_eff = safeT(temperature);

  const Real mol = massFractionToMolalConc(xnacl);
  const Real mol2 = mol * mol;
  const Real mol3 = mol2 * mol;

  const Real Tc = T_eff - _T_c2k;

  const Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 +
                 0.629e-3 * Tc * (1.0 - std::exp(-0.7 * mol));

  const Real water_viscosity = _water_fp->mu_from_p_T(p_eff, T_eff);

  return a * water_viscosity;
}

void
BrineFluidProperties::mu_from_p_T_X(Real pressure,
                                    Real temperature,
                                    Real xnacl,
                                    Real & mu,
                                    Real & dmu_dp,
                                    Real & dmu_dT,
                                    Real & dmu_dx) const
{
  const bool t_bad = !std::isfinite(temperature);
  const Real p_eff = safePressure(pressure);
  const Real T_eff = safeT(temperature);

  Real muw, dmuw_dp, dmuw_dT;
  _water_fp->mu_from_p_T(p_eff, T_eff, muw, dmuw_dp, dmuw_dT);

  Real mol = massFractionToMolalConc(xnacl);
  Real dmol_dx = 1.0 / ((1.0 - xnacl) * (1.0 - xnacl) * _Mnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  Real Tc = T_eff - _T_c2k;

  Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 +
           0.629e-3 * Tc * (1.0 - std::exp(-0.7 * mol));
  Real da_dx =
      (0.0816 + 0.0244 * mol + 3.84e-4 * mol2 + 4.403e-4 * Tc * std::exp(-0.7 * mol)) * dmol_dx;
  Real da_dT = 0.629e-3 * (1.0 - std::exp(-0.7 * mol));

  mu = a * muw;
  dmu_dp = pressure < 1.0 ? 0.0 : (a * dmuw_dp);
  dmu_dT = t_bad ? 0.0 : (da_dT * muw + a * dmuw_dT);
  dmu_dx = da_dx * muw;
}

FPADReal
BrineFluidProperties::h_from_p_T_X(const FPADReal & pressure,
                                   const FPADReal & temperature,
                                   const FPADReal & xnacl) const
{
  // Driesner Th virtual-temperature mapping bypassed — evaluate at real (p,T).
  const Real p_safe = safePressure(pressure.value());
  const Real T_safe = safeT(temperature.value());

  FPADReal enthalpy;
  if (_water_fp_derivs)
  {
    Real h, dh_dp, dh_dT;
    _water_fp->h_from_p_T(p_safe, T_safe, h, dh_dp, dh_dT);
    enthalpy = h;
    enthalpy.derivatives() = pressure.derivatives() * dh_dp +
                              temperature.derivatives() * dh_dT;
  }
  else
    enthalpy = _water_fp->h_from_p_T(p_safe, T_safe);

  return enthalpy;
}


Real
BrineFluidProperties::h_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  return _water_fp->h_from_p_T(safePressure(pressure), safeT(temperature));
}

void
BrineFluidProperties::h_from_p_T_X(Real pressure,
                                   Real temperature,
                                   Real xnacl,
                                   Real & h,
                                   Real & dh_dp,
                                   Real & dh_dT,
                                   Real & dh_dx) const
{
  const bool t_bad = !std::isfinite(temperature);
  FPADReal p = safePressure(pressure);
  Moose::derivInsert(p.derivatives(), 0, pressure >= 1.0 ? 1.0 : 0.0);
  FPADReal T = safeT(temperature);
  Moose::derivInsert(T.derivatives(), 1, t_bad ? 0.0 : 1.0);
  FPADReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPADReal ad_h = h_from_p_T_X(p, T, x);

  h = ad_h.value();
  dh_dp = ad_h.derivatives()[0];
  dh_dT = ad_h.derivatives()[1];
  dh_dx = ad_h.derivatives()[2];
}

Real
BrineFluidProperties::cp_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  // Driesner Th virtual-temperature mapping bypassed — evaluate at real (p,T).
  return _water_fp->cp_from_p_T(safePressure(pressure), safeT(temperature));
}


FPADReal
BrineFluidProperties::e_from_p_T_X(const FPADReal & pressure,
                                   const FPADReal & temperature,
                                   const FPADReal & xnacl) const
{
  const FPADReal p_eff = pressure;
  const FPADReal T_eff = temperature;

  FPADReal enthalpy = h_from_p_T_X(p_eff, T_eff, xnacl);
  FPADReal density = rho_from_p_T_X(p_eff, T_eff, xnacl);

  return enthalpy - p_eff / density;
}

Real
BrineFluidProperties::e_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  const Real p_eff = safePressure(pressure);
  const Real T_eff = safeT(temperature);

  Real enthalpy = h_from_p_T_X(p_eff, T_eff, xnacl);
  Real density = rho_from_p_T_X(p_eff, T_eff, xnacl);

  return enthalpy - p_eff / density;
}

void
BrineFluidProperties::e_from_p_T_X(Real pressure,
                                   Real temperature,
                                   Real xnacl,
                                   Real & e,
                                   Real & de_dp,
                                   Real & de_dT,
                                   Real & de_dx) const
{
  const bool t_bad = !std::isfinite(temperature);
  FPADReal p = safePressure(pressure);
  Moose::derivInsert(p.derivatives(), 0, pressure >= 1.0 ? 1.0 : 0.0);
  FPADReal T = safeT(temperature);
  Moose::derivInsert(T.derivatives(), 1, t_bad ? 0.0 : 1.0);
  FPADReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPADReal ad_e = e_from_p_T_X(p, T, x);

  e = ad_e.value();
  de_dp = ad_e.derivatives()[0];
  de_dT = ad_e.derivatives()[1];
  de_dx = ad_e.derivatives()[2];
}

Real
BrineFluidProperties::k_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  const Real p_eff = safePressure(pressure);
  const Real T_eff = safeT(temperature);

  Real mol = massFractionToMolalConc(xnacl);
  Real Tc = T_eff - _T_c2k;

  Real S = 100.0 * _Mnacl * mol / (1.0 + _Mnacl * mol);
  Real lambdaw = _water_fp->k_from_p_T(p_eff, T_eff);
  Real lambda = 1.0 - (2.3434e-3 - 7.924e-6 * Tc + 3.924e-8 * Tc * Tc) * S +
                (1.06e-5 - 2.0e-8 * Tc - 1.2e-10 * Tc * Tc) * S * S;

  return lambda * lambdaw;
}

Real
BrineFluidProperties::vaporPressure(Real temperature, Real xnacl) const
{
  Real T_eff = safeT(temperature);

  Real mol = massFractionToMolalConc(xnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  Real a = 1.0 + 5.93582e-6 * mol - 5.19386e-5 * mol2 + 1.23156e-5 * mol3;
  Real b = 1.1542e-6 * mol + 1.41254e-7 * mol2 - 1.92476e-8 * mol3 - 1.70717e-9 * mol * mol3 +
           1.0539e-10 * mol2 * mol3;

  Real th20 = std::exp(std::log(T_eff) / (a + b * T_eff));

  // SimpleFluidProperties has no vaporPressure(); delegate to Water97.
  return _water97_fp->vaporPressure(th20);
}

Real
BrineFluidProperties::haliteSolubility(Real temperature) const
{
  Real T_eff = safeT(temperature);
  Real Tc = T_eff - _T_c2k;

  return (26.18 + 7.2e-3 * Tc + 1.06e-4 * Tc * Tc) / 100.0;
}

Real
BrineFluidProperties::massFractionToMolalConc(Real xnacl) const
{
  return xnacl / ((1.0 - xnacl) * _Mnacl);
}

Real
BrineFluidProperties::massFractionToMoleFraction(Real xnacl) const
{
  Real Mbrine = molarMass(xnacl);
  return xnacl * Mbrine / _Mnacl;
}

FPADReal
BrineFluidProperties::massFractionToMoleFraction(const FPADReal & xnacl) const
{
  FPADReal Mbrine = molarMass(xnacl);
  return xnacl * Mbrine / _Mnacl;
}

Real
BrineFluidProperties::henryConstant(Real temperature, const std::vector<Real> & coeffs) const
{
  // Henry remains Water97-based by design.
  return _water97_fp->henryConstant(safeT(temperature), coeffs);
}

void
BrineFluidProperties::henryConstant(Real temperature,
                                    const std::vector<Real> & coeffs,
                                    Real & Kh,
                                    Real & dKh_dT) const
{
  _water97_fp->henryConstant(safeT(temperature), coeffs, Kh, dKh_dT);
}

ADReal
BrineFluidProperties::henryConstant(const ADReal & temperature,
                                    const std::vector<Real> & coeffs) const
{
  Real Kh, dKh_dT;
  henryConstant(temperature.value(), coeffs, Kh, dKh_dT);

  // If T is NaN, safeT() used fallback and this should be treated constant
  if (!std::isfinite(temperature.value()))
    return ADReal(Kh);

  ADReal henry = Kh;
  henry.derivatives() = temperature.derivatives() * dKh_dT;
  return henry;
}
