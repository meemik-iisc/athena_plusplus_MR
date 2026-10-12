//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file runaway_bh.cpp
//! \brief Runaway black hole moving through a uniform CGM (Bondi-Hoyle-Lyttleton flow).
//!
//! Simulated in the BH rest frame using cylindrical coordinates (R, phi, z) with the BH
//! at the origin, on the symmetry axis. The CGM enters at the inner-x3 boundary moving
//! in +z. The flow is axisymmetric, so phi is covered by only a few cells with periodic
//! boundaries. Gravity is a point mass, a = -GM r_vec/(r^2+eps^2)^{3/2} with
//! r = sqrt(R^2+z^2), added as a user source term. A sink region r < r_sink is reset
//! every step and the removed mass is recorded as the accretion rate.
//!
//! The z-grid is a sinh mapping that concentrates cells around z = 0 (set x3rat < 0).

// C headers

// C++ headers
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
Real gm, eps2, r_sink2, d_sink, p_sink;
Real d_amb, p_amb, v_amb;
Real zgrid_b, zgrid_x0, zgrid_a;

void PointMassGravity(MeshBlock *pmb, const Real time, const Real dt,
                      const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                      const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                      AthenaArray<Real> &cons_scalar);
void InflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);
Real SinhMeshX3(Real x, RegionSize rs);
Real MdotHst(MeshBlock *pmb, int iout);
Real MaccHst(MeshBlock *pmb, int iout);
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Read parameters and enroll gravity, inflow boundary, z-grid and history outputs
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in runaway_bh.cpp ProblemGenerator" << std::endl
        << "This problem requires --coord=cylindrical" << std::endl;
    ATHENA_ERROR(msg);
  }

  gm    = pin->GetReal("problem", "GM_bh");
  d_amb = pin->GetReal("problem", "d_amb");
  p_amb = pin->GetReal("problem", "p_amb");
  v_amb = pin->GetReal("problem", "v_amb");
  Real r_sink = pin->GetReal("problem", "r_sink");
  Real eps    = pin->GetOrAddReal("problem", "softening", 0.5*r_sink);
  r_sink2 = SQR(r_sink);
  eps2 = SQR(eps);
  // sink gas: low density, same temperature as the ambient CGM
  d_sink = pin->GetOrAddReal("problem", "d_sink", 1.0e-2*d_amb);
  p_sink = d_sink*p_amb/d_amb;

  EnrollUserExplicitSourceFunction(PointMassGravity);
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, InflowInnerX3);
  }

  // sinh z-grid: z(x) = A sinh(b (x - x0)), x in [0,1], solved so that
  // z(0) = x3min and z(1) = x3max. Larger b concentrates more cells at z = 0.
  if (mesh_size.x3rat < 0.0) {
    zgrid_b = pin->GetReal("problem", "zgrid_b");
    Real z1 = mesh_size.x3min, z2 = mesh_size.x3max;
    if (z1 >= 0.0 || z2 <= 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in runaway_bh.cpp" << std::endl
          << "sinh z-grid requires x3min < 0 < x3max" << std::endl;
      ATHENA_ERROR(msg);
    }
    // bisection for x0: sinh(b x0)/sinh(b (1-x0)) = -z1/z2 (monotonic in x0)
    Real lo = 0.0, hi = 1.0;
    for (int n=0; n<200; ++n) {
      Real mid = 0.5*(lo + hi);
      Real f = std::sinh(zgrid_b*mid)/std::sinh(zgrid_b*(1.0 - mid)) + z1/z2;
      if (f > 0.0) hi = mid; else lo = mid;
    }
    zgrid_x0 = 0.5*(lo + hi);
    zgrid_a = z2/std::sinh(zgrid_b*(1.0 - zgrid_x0));
    EnrollUserMeshGenerator(X3DIR, SinhMeshX3);
  }

  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, MdotHst, "mdot");
  EnrollUserHistoryOutput(1, MaccHst, "m_acc");
}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//! \brief ruser_meshblock_data[0](0): mass removed by the sink in the last step
//!        ruser_meshblock_data[0](1): cumulative mass removed
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(2);
  ruser_meshblock_data[0](0) = 0.0;
  ruser_meshblock_data[0](1) = 0.0;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Uniform CGM moving in +z everywhere
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gm1 = peos->GetGamma() - 1.0;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = d_amb;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = d_amb*v_amb;
        phydro->u(IEN,k,j,i) = p_amb/gm1 + 0.5*d_amb*SQR(v_amb);
      }
    }
  }
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//! \brief Sink: reset cells inside r_sink to low-density gas at rest, record removed mass
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  Real gm1 = peos->GetGamma() - 1.0;
  Real dm = 0.0;
  for (int k=ks; k<=ke; ++k) {
    Real z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real r2 = SQR(pcoord->x1v(i)) + SQR(z);
        if (r2 < r_sink2) {
          Real d = phydro->u(IDN,k,j,i);
          if (d > d_sink) dm += (d - d_sink)*pcoord->GetCellVolume(k,j,i);
          phydro->u(IDN,k,j,i) = d_sink;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = p_sink/gm1;
          phydro->w(IDN,k,j,i) = d_sink;
          phydro->w(IVX,k,j,i) = 0.0;
          phydro->w(IVY,k,j,i) = 0.0;
          phydro->w(IVZ,k,j,i) = 0.0;
          phydro->w(IPR,k,j,i) = p_sink;
        }
      }
    }
  }
  ruser_meshblock_data[0](0) = dm;
  ruser_meshblock_data[0](1) += dm;
}

namespace {
//----------------------------------------------------------------------------------------
//! \brief Point-mass gravity of the BH at the origin, spherical 1/r^2 force

void PointMassGravity(MeshBlock *pmb, const Real time, const Real dt,
                      const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                      const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                      AthenaArray<Real> &cons_scalar) {
  Coordinates *pco = pmb->pcoord;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real z = pco->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rcyl = pco->x1v(i);
        Real r2 = SQR(rcyl) + SQR(z) + eps2;
        Real fac = -gm/(r2*std::sqrt(r2));
        Real a_r = fac*rcyl;
        Real a_z = fac*z;
        Real den = prim(IDN,k,j,i);
        cons(IM1,k,j,i) += dt*den*a_r;
        cons(IM3,k,j,i) += dt*den*a_z;
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) += dt*den*(prim(IVX,k,j,i)*a_r + prim(IVZ,k,j,i)*a_z);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Fixed CGM inflow at the inner-x3 (upstream) boundary

void InflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,kl-k,j,i) = d_amb;
        prim(IVX,kl-k,j,i) = 0.0;
        prim(IVY,kl-k,j,i) = 0.0;
        prim(IVZ,kl-k,j,i) = v_amb;
        prim(IPR,kl-k,j,i) = p_amb;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief sinh z-grid concentrating cells around the BH at z = 0

Real SinhMeshX3(Real x, RegionSize rs) {
  return zgrid_a*std::sinh(zgrid_b*(x - zgrid_x0));
}

//----------------------------------------------------------------------------------------
//! \brief History outputs: accretion rate in the last step, and cumulative accreted mass

Real MdotHst(MeshBlock *pmb, int iout) {
  Real dt = pmb->pmy_mesh->dt;
  return (dt > 0.0) ? pmb->ruser_meshblock_data[0](0)/dt : 0.0;
}

Real MaccHst(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[0](1);
}
} // namespace
