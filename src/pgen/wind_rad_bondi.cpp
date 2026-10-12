//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
    Real M_bh, GMbh, v_bh, epsilon, n_cgm, T_cgm, rho_cgm, pres_cgm, gamma_gas, mu, gm1;

void AddBHGrav(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                AthenaArray<Real> &cons_scalar);
void WindZBC(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
            FaceField &b, Real time, Real dt,
            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
}

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

  M_bh     = pin->GetReal("problem", "M_bh");
  v_bh     = pin->GetReal("problem", "v_bh");
  epsilon  = pin->GetReal("problem", "epsilon");
  n_cgm    = pin->GetReal("problem", "n_cgm");
  T_cgm    = pin->GetReal("problem", "T_cgm");
  gamma_gas= pin->GetReal("hydro", "gamma"); 
  mu       = pin->GetOrAddReal("problem", "mu",1.0); 
  
  gm1      = gamma_gas - 1.0;
  
  rho_cgm  = n_cgm*mu;
  pres_cgm = (rho_cgm*T_cgm)/(mu*punit->code_temperature_mu_cgs);
  Real G_code   = punit->grav_const_code;
  GMbh     = G_code*M_bh;  
  punit->PrintCodeUnits();
  EnrollUserExplicitSourceFunction(AddBHGrav);
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, WindZBC);
  }

  // sinh z-grid: z(x) = A sinh(b (x - x0)), x in [0,1], solved so that
  // z(0) = x3min and z(1) = x3max. Larger b concentrates more cells at z = 0.
//   if (mesh_size.x3rat < 0.0) {
//     zgrid_b = pin->GetReal("problem", "zgrid_b");
//     Real z1 = mesh_size.x3min, z2 = mesh_size.x3max;
//     if (z1 >= 0.0 || z2 <= 0.0) {
//       std::stringstream msg;
//       msg << "### FATAL ERROR in runaway_bh.cpp" << std::endl
//           << "sinh z-grid requires x3min < 0 < x3max" << std::endl;
//       ATHENA_ERROR(msg);
//     }
//     // bisection for x0: sinh(b x0)/sinh(b (1-x0)) = -z1/z2 (monotonic in x0)
//     Real lo = 0.0, hi = 1.0;
//     for (int n=0; n<200; ++n) {
//       Real mid = 0.5*(lo + hi);
//       Real f = std::sinh(zgrid_b*mid)/std::sinh(zgrid_b*(1.0 - mid)) + z1/z2;
//       if (f > 0.0) hi = mid; else lo = mid;
//     }
//     zgrid_x0 = 0.5*(lo + hi);
//     zgrid_a = z2/std::sinh(zgrid_b*(1.0 - zgrid_x0));
//     EnrollUserMeshGenerator(X3DIR, SinhMeshX3);
//   }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Uniform CGM moving in -z everywhere (wind enters from outer-x3)
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho_cgm;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = -rho_cgm*v_bh;
        phydro->u(IEN,k,j,i) = pres_cgm/gm1 + 0.5*rho_cgm*SQR(v_bh);
      }
    }
  }
}

namespace {
//----------------------------------------------------------------------------------------
//! \brief Point-mass gravity of the BH at the origin, spherical 1/r^2 force

void AddBHGrav(MeshBlock *pmb, const Real time, const Real dt,
                      const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                      const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                      AthenaArray<Real> &cons_scalar) {
  Coordinates *pco = pmb->pcoord;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real z = pco->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rcyl = pco->x1v(i);
        Real r2 = SQR(rcyl) + SQR(z) + SQR(epsilon);

        Real vr     = prim(IVX,k,j,i);
        Real vz     = prim(IVZ,k,j,i);
        Real rho    = prim(IDN,k,j,i);
        Real grad_phi_by_r  = (GMbh)/pow(r2,1.5);
            
        cons(IM1,k,j,i) -= rho*grad_phi_by_r*dt*rcyl;
        cons(IM3,k,j,i) -= rho*grad_phi_by_r*dt*z;
        Real v_dot_r    = (vr*rcyl)+(vz*z);
        cons(IEN,k,j,i) -= rho*v_dot_r*grad_phi_by_r*dt;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Fixed CGM inflow at the outer-x3 (upstream) boundary, moving in -z

void WindZBC(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,ku+k,j,i) = rho_cgm;
        prim(IVX,ku+k,j,i) = 0.0;
        prim(IVY,ku+k,j,i) = 0.0;
        prim(IVZ,ku+k,j,i) = -v_bh;
        prim(IPR,ku+k,j,i) = pres_cgm;
      }
    }
  }
}
} // namespace
