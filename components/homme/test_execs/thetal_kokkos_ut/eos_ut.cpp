#include <catch2/catch.hpp>

#include <iostream>
#include <random>

#include "EquationOfState.hpp"
#include "Types.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"

using namespace Homme;

// ============= EQUATION OF STATE ================ //

extern "C" {
void init_f90 (const Real* hyai_ptr, const Real& ps0, const bool& hydrostatic);
void pnh_and_exner_from_eos_f90(const int& num_elems,
                                const Real*& vtheta_dp,
                                const Real*& dp,
                                const Real*& phi_i,
                                      Real*& pnh,
                                      Real*& exner,
                                      Real*& dpnh_dp_i);
void phi_from_eos_f90(const int& num_elems,
                      const Real*& phis,
                      const Real*& vtheta_dp,
                      const Real*& dp,
                      Real*& phi_i);
} // extern "C"

TEST_CASE("eos", "eos") {

  constexpr auto LAST_LEV_P = ColInfo<NUM_INTERFACE_LEV>::LastPack;
  constexpr auto LAST_INTERFACE_VEC_IDX = ColInfo<NUM_INTERFACE_LEV>::LastVecEnd;

  constexpr int num_elems = 2;

  // F90 and CXX views (and host mirrors of cxx views)
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  vtheta_dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]> phi_i_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  pnh_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  exner_f90("",num_elems);
  HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]> dpnh_dp_i_f90("",num_elems);
  HostViewManaged<Real*[NP][NP]> phis_f90("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   vtheta_dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> phi_i_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   pnh_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   exner_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> dpnh_dp_i_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> dp_i_cxx("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   p_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> p_i_cxx("",num_elems);

  ExecViewManaged<Real*[NP][NP]> phis_cxx("",num_elems);

  auto h_vtheta_dp = Kokkos::create_mirror_view(vtheta_dp_cxx);
  auto h_phi_i = Kokkos::create_mirror_view(phi_i_cxx);
  auto h_dp = Kokkos::create_mirror_view(dp_cxx);
  auto h_pnh = Kokkos::create_mirror_view(pnh_cxx);
  auto h_exner = Kokkos::create_mirror_view(exner_cxx);
  auto h_dpnh_dp_i = Kokkos::create_mirror_view(dpnh_dp_i_cxx);

  std::random_device rd;
  using rngAlg = std::mt19937_64;
  rngAlg engine(rd());
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);

  HybridVCoord hvcoord;
  hvcoord.random_init(rd());

  decltype(hvcoord.hybrid_ai)::HostMirror hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  const Real* hyai_ptr = hyai.data();

  ColumnOps col_ops;
  EquationOfState eos;

  SECTION ("pnh_and_exner") {
    for (bool hydrostatic : {false,true}) {
      // Init f90
      init_f90(hyai_ptr,hvcoord.ps0,hydrostatic);

      eos.init(hydrostatic,hvcoord);

      // Create random inputs
      genRandArray(vtheta_dp_cxx,engine,pdf);
      genRandArray(dp_cxx,engine,pdf);

      sync_to_host(vtheta_dp_cxx,vtheta_dp_f90);
      sync_to_host(dp_cxx,dp_f90);

      // Note: to avoid errors in pnh_and_exner_from_eos, we need phi to be increasing.
      //       Rather than using a constraint (which may call the function many times,
      //       we simply ask that there are no duplicates, then we sort it later.
      auto sort_and_chek = [](const ExecViewManaged<Scalar[NUM_LEV_P]>::HostMirror v)->bool {
        Real* start = reinterpret_cast<Real*>(v.data());
        Real* end   = reinterpret_cast<Real*>(v.data()) + NUM_LEV_P*VECTOR_SIZE;
        std::sort(start,end);
        std::reverse(start,end);
        auto it = std::unique(start,end);
        return it==end;
      };
      for (int ie=0; ie<num_elems; ++ie) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++ jgp) {
            genRandArray(Homme::subview(phi_i_cxx,ie,igp,jgp),engine,pdf,sort_and_chek);
          }
        }
      }
      sync_to_host(phi_i_cxx,phi_i_f90);

      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                           KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                             [&](const int idx) {
          const int igp = idx / NP;
          const int jgp = idx % NP;

          auto pnh = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
          auto exner = Homme::subview(exner_cxx,kv.ie,igp,jgp);
          if (hydrostatic) {
            auto dp  = Homme::subview(dp_cxx,kv.ie,igp,jgp);
            auto p_i = Homme::subview(p_i_cxx,kv.ie,igp,jgp);

            eos.compute_hydrostatic_p(kv,dp,p_i,pnh);
            eos.compute_exner(kv,pnh,exner);
          } else {
            // Compute pnh and exner
            auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
            auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);

            eos.compute_pnh_and_exner(kv,vtheta_dp,phi_i,pnh,exner);
          }
        });
        kv.team_barrier();
      });
      Kokkos::fence();

      // Compute dpnh_dp_i
      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                           KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                             [&](const int idx) {
          const int igp = idx / NP;
          const int jgp = idx % NP;

          auto dp   = Homme::subview(dp_cxx,kv.ie,igp,jgp);
          auto dp_i = Homme::subview(dp_i_cxx,kv.ie,igp,jgp);

          col_ops.compute_interface_values(kv,dp,dp_i);

          auto pnh = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
          auto dpnh_dp_i = Homme::subview(dpnh_dp_i_cxx,kv.ie,igp,jgp);

          eos.compute_dpnh_dp_i(kv,pnh,dp_i,dpnh_dp_i);
        });
        kv.team_barrier();
      });
      Kokkos::fence();

      // Run f90
      const Real* vtheta_dp_ptr = vtheta_dp_f90.data();
      const Real* dp_ptr = dp_f90.data();
      const Real* phi_i_ptr = phi_i_f90.data();
      Real* pnh_ptr = pnh_f90.data();
      Real* exner_ptr = exner_f90.data();
      Real* dpnh_dp_i_ptr = dpnh_dp_i_f90.data();
      pnh_and_exner_from_eos_f90(num_elems,vtheta_dp_ptr,dp_ptr,phi_i_ptr,pnh_ptr,exner_ptr,dpnh_dp_i_ptr);

      // Now, compare results
      Kokkos::deep_copy(h_exner,exner_cxx);
      Kokkos::deep_copy(h_pnh,pnh_cxx);
      Kokkos::deep_copy(h_dpnh_dp_i,dpnh_dp_i_cxx);

      for (int ie=0; ie<num_elems; ++ie) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++jgp) {
            for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
              const int ilev = k / VECTOR_SIZE;
              const int ivec = k % VECTOR_SIZE;

              REQUIRE (h_exner(ie,igp,jgp,ilev)[ivec] == exner_f90(ie,k,igp,jgp));
              REQUIRE (h_pnh(ie,igp,jgp,ilev)[ivec] == pnh_f90(ie,k,igp,jgp));
              REQUIRE (h_dpnh_dp_i(ie,igp,jgp,ilev)[ivec] == dpnh_dp_i_f90(ie,k,igp,jgp));
            }
            REQUIRE (h_dpnh_dp_i(ie,igp,jgp,LAST_LEV_P)[LAST_INTERFACE_VEC_IDX] == dpnh_dp_i_f90(ie,NUM_INTERFACE_LEV-1,igp,jgp));
          }
        }
      }
    }
  }

  SECTION ("phi") {
    // Generate random inputs
    genRandArray(vtheta_dp_cxx,engine,pdf);
    genRandArray(dp_cxx,engine,pdf);
    genRandArray(phis_cxx,engine,pdf);

    Homme::sync_to_host(vtheta_dp_cxx,vtheta_dp_f90);
    Homme::sync_to_host(dp_cxx,dp_f90);
    Homme::sync_to_host(phis_cxx,phis_f90);

    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        auto dp  = Homme::subview(dp_cxx,kv.ie,igp,jgp);
        auto p_i = Homme::subview(p_i_cxx,kv.ie,igp,jgp);
        auto p   = Homme::subview(p_cxx,kv.ie,igp,jgp);

        // Compute p_i from dp with an exclusive, forward scan sum
        p_i(0)[0] = hvcoord.hybrid_ai(0)*hvcoord.ps0;
        col_ops.column_scan_mid_to_int<true>(kv,dp,p_i);

        // Compute p from p_i
        col_ops.compute_midpoint_values(kv,p_i,p);

        // Now compute the geopotential
        auto phis = phis_cxx(kv.ie,igp,jgp);
        auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
        auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);

        eos.compute_phi_i(kv,phis,vtheta_dp,p,phi_i);
      });
      kv.team_barrier();
    });
    Kokkos::fence();

    // Run f90
    const Real* phis_ptr = phis_f90.data();
    const Real* vtheta_dp_ptr = vtheta_dp_f90.data();
    const Real* dp_ptr = dp_f90.data();
          Real* phi_i_ptr = phi_i_f90.data();
    phi_from_eos_f90(num_elems,phis_ptr,vtheta_dp_ptr,dp_ptr,phi_i_ptr);

    // Now, compare results
    Kokkos::deep_copy(h_phi_i,phi_i_cxx);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            REQUIRE (h_phi_i(ie,igp,jgp,ilev)[ivec] == phi_i_f90(ie,k,igp,jgp));
          }
        }
      }
    }
  }

  // Note: this test cannot be a BFB test, meaning that you cannot ask that
  //       phi_in == phi(pnh(phi_in)) bit-for-bit, due to roundoffs in the
  //       calculations. Therefore, we enforce |phi(pnh(phi_in))-phi_in|<tol.
  SECTION ("phi_pnh_composition") {
    // This test checks that phi = phi(pnh(phi))

    // The two eos routines are reciprocally dual only in non-hydrostatic mode
    eos.init(false,hvcoord);

    // Create random inputs
    genRandArray(vtheta_dp_cxx,engine,pdf);
    genRandArray(dp_cxx,engine,pdf);
    // Note: to avoid errors in pnh_and_exner_from_eos, we need phi to be increasing.
    //       Rather than using a constraint (which may call the function many times,
    //       we simply ask that there are no duplicates, then we sort it later.
    // auto sort_and_chek = [](const ExecViewManaged<Scalar[NUM_LEV_P]>::HostMirror v)->bool {
    //   Real* start = reinterpret_cast<Real*>(v.data());
    //   Real* end   = reinterpret_cast<Real*>(v.data()) + NUM_LEV_P*VECTOR_SIZE;
    //   std::sort(start,end);
    //   std::reverse(start,end);
    //   auto it = std::unique(start,end);
    //   return it==end;
    // };
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++ jgp) {
          // genRandArray(Homme::subview(phi_i_cxx,ie,igp,jgp),engine,pdf,sort_and_chek);
          auto phi = Homme::subview(phi_i_cxx,ie,igp,jgp);
          auto hphi = Kokkos::create_mirror_view(phi);
          auto rhphi = viewAsReal(hphi);
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            rhphi(k) = NUM_INTERFACE_LEV-k;
          }
          Kokkos::deep_copy(phi,hphi);
        }
      }
    }

    // Store initial phi.
    decltype(phi_i_cxx)::HostMirror h_phi_i_in("",num_elems);
    Kokkos::deep_copy(h_phi_i_in,phi_i_cxx);

    // phis should be phi at the bottom. Make it happen
    auto h_phis = Kokkos::create_mirror_view(phis_cxx);
    using Helper = ColInfo<NUM_INTERFACE_LEV>;
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++ jgp) {
          h_phis(ie,igp,jgp) = h_phi_i_in(ie,igp,jgp,int(Helper::LastPack))[Helper::LastVecEnd];
        }
      }
    }
    Kokkos::deep_copy(phis_cxx,h_phis);

    // Compute pnh(phi)
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        // Compute pnh and exner
        auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
        auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);
        auto pnh = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
        auto exner = Homme::subview(exner_cxx,kv.ie,igp,jgp);

        eos.compute_pnh_and_exner(kv,vtheta_dp,phi_i,pnh,exner);
      });
      kv.team_barrier();
    });
    Kokkos::fence();

    // Compute phi(pnh) using exponentials
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        auto p = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
        auto phis = phis_cxx(kv.ie,igp,jgp);
        auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
        auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);

        eos.compute_phi_i(kv,phis,vtheta_dp,p,phi_i);
      });
      kv.team_barrier();
    });
    Kokkos::fence();

    // Now, compare results
    Kokkos::deep_copy(h_phi_i,phi_i_cxx);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            REQUIRE (compare_answers(h_phi_i_in(ie,igp,jgp,ilev)[ivec],h_phi_i(ie,igp,jgp,ilev)[ivec])<1e-14);
          }
        }
      }
    }

    // Compute phi(pnh,exner)
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        auto p = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
        auto exner = Homme::subview(exner_cxx,kv.ie,igp,jgp);
        auto phis = phis_cxx(kv.ie,igp,jgp);
        auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
        auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);

        eos.compute_phi_i(kv,phis,vtheta_dp,p,exner,phi_i);
      });
      kv.team_barrier();
    });
    Kokkos::fence();

    // Now, compare results
    Kokkos::deep_copy(h_phi_i,phi_i_cxx);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            REQUIRE (compare_answers(h_phi_i_in(ie,igp,jgp,ilev)[ivec],h_phi_i(ie,igp,jgp,ilev)[ivec])<1e-14);
          }
        }
      }
    }
  }
}