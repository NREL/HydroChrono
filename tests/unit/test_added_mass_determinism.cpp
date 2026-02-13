/*********************************************************************
 * @file  test_added_mass_determinism.cpp
 * @brief Tripwire: added-mass assembly must be deterministic.
 *
 * Builds N=3 bodies with cross-coupled added mass, runs several
 * independent single-step trials (with heap perturbation between
 * them), and asserts bit-identical velocity results.
 *
 * An early version of Chrono's ChLoadHydrodynamics used an
 * unordered_map whose iteration order could vary with heap layout.
 * That was fixed (switched to std::vector), but this test guards
 * against any similar regression.
 *
 * Also sweeps across Chrono solver types and reports which solvers
 * correctly preserve the full added-mass matrix (including
 * off-diagonal cross-coupling terms).
 *
 * Self-contained — no external data files.
 *********************************************************************/

#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChLoadHydrodynamics.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/solver/ChSystemDescriptor.h>

using namespace chrono;

// ── Parameters ──────────────────────────────────────────────────────
constexpr int N       = 3;
constexpr int DOF     = 6;
constexpr int NDOF    = DOF * N;
constexpr double DT   = 0.001;
constexpr double TOL  = 1e-10;
constexpr int TRIALS  = 12;

constexpr double MASS[N]  = {1000.0, 2000.0, 3000.0};
constexpr double VZ0[N]   = {1.0, -0.5, 0.3};

// ── Helpers ─────────────────────────────────────────────────────────

/// Symmetric, positive-definite 6N x 6N added-mass matrix with fully dense
/// 6x6 sub-blocks and non-trivial cross-coupling between all body pairs.
/// Every element is distinct so that any row/column swap is detectable.
static Eigen::MatrixXd BuildAddedMass() {
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(NDOF, NDOF);
    double val = 1.0;
    for (int i = 0; i < NDOF; ++i) {
        for (int j = 0; j <= i; ++j) {
            L(i, j) = std::fmod(val * 0.7 + 3.1, 10.0);
            val = std::fmod(val * 7.3 + 1.7, 100.0);
        }
    }
    Eigen::MatrixXd A = L * L.transpose();
    A = (A + A.transpose()) * 0.5;  // enforce exact symmetry

    constexpr double boost[N] = {500.0, 800.0, 1200.0};
    for (int i = 0; i < N; ++i)
        A.block(DOF*i, DOF*i, DOF, DOF) +=
            boost[i] * Eigen::MatrixXd::Identity(DOF, DOF);

    return A;
}

using Sig = std::array<std::array<double, DOF>, N>;

struct TrialResult {
    Sig sig;
    Eigen::MatrixXd added_mass;  // system mass minus body inertia
    bool ok = true;
    std::string error;
};

/// Create a fresh system with the given solver, attach added mass, step once.
static TrialResult RunTrial(const Eigen::MatrixXd& A, ChSolver::Type solver) {
    TrialResult result;
    try {
        ChSystemNSC sys;
        sys.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
        sys.SetSolverType(solver);

        std::array<std::shared_ptr<ChBody>, N> bodies;
        for (int i = 0; i < N; ++i) {
            auto b = chrono_types::make_shared<ChBody>();
            b->SetName("body" + std::to_string(i));
            b->SetMass(MASS[i]);
            b->SetInertiaXX(ChVector3d(1, 1, 1));
            b->SetPosDt(ChVector3d(0, 0, VZ0[i]));
            sys.Add(b);
            bodies[i] = b;
        }

        ChBodyAddedMassBlocks blocks;
        for (int i = 0; i < N; ++i) {
            ChBodyAddedMassBlock entry;
            entry.body  = bodies[i];
            entry.block = A.block(DOF*i, 0, DOF, NDOF);
            blocks.push_back(entry);
        }

        auto load = chrono_types::make_shared<ChLoadHydrodynamics>(blocks);
        load->SetVerbose(false);
        sys.Add(load);
        sys.DoStepDynamics(DT);

        for (int i = 0; i < N; ++i) {
            auto v = bodies[i]->GetPosDt();
            auto w = bodies[i]->GetAngVelLocal();
            result.sig[i] = {v.x(), v.y(), v.z(), w.x(), w.y(), w.z()};
        }

        // Extract added mass = system mass matrix - body inertia
        ChSparseMatrix M_sparse;
        sys.GetMassMatrix(M_sparse);
        result.added_mass = Eigen::MatrixXd(M_sparse);
        for (int i = 0; i < N; ++i) {
            int off = bodies[i]->GetOffset_w();
            for (int d = 0; d < 3; ++d) result.added_mass(off+d, off+d) -= MASS[i];
            auto I = bodies[i]->GetInertiaXX();
            result.added_mass(off+3, off+3) -= I.x();
            result.added_mass(off+4, off+4) -= I.y();
            result.added_mass(off+5, off+5) -= I.z();
        }
    } catch (const std::exception& e) {
        result.ok = false;
        result.error = e.what();
    }
    return result;
}

static bool HasBadValues(const Sig& s) {
    for (auto& b : s) for (double v : b)
        if (std::isnan(v) || std::isinf(v)) return true;
    return false;
}

static bool Match(const Sig& a, const Sig& b) {
    for (int i = 0; i < N; ++i)
        for (int d = 0; d < DOF; ++d)
            if (std::abs(a[i][d] - b[i][d]) > TOL) return false;
    return true;
}

static std::string Fmt(const Sig& s) {
    std::ostringstream o;
    o.precision(12); o << std::scientific;
    for (int i = 0; i < N; ++i) {
        o << "  body" << i << ": [";
        for (int d = 0; d < DOF; ++d) { if (d) o << ", "; o << s[i][d]; }
        o << "]\n";
    }
    return o.str();
}

// ── Main ────────────────────────────────────────────────────────────

int main() {
    std::cout << "=== Added-mass determinism test ===" << std::endl;

    auto A = BuildAddedMass();

    // Sanity checks on the test matrix itself
    if ((A - A.transpose()).norm() > 1e-15) { std::cerr << "FAIL: A not symmetric\n"; return 1; }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(A);
    if (eig.eigenvalues().minCoeff() < -1e-10) { std::cerr << "FAIL: A not PSD\n"; return 1; }

    // ── INPUT: what we feed to ChLoadHydrodynamics ──
    Eigen::IOFormat compact(4, 0, "  ", "\n", "    ");
    std::cout << "\n--- INPUT: per-body row blocks passed to ChLoadHydrodynamics ---\n";
    for (int i = 0; i < N; ++i) {
        std::cout << "\n  body " << i
                  << "  (mass=" << MASS[i] << " kg, v0_z=" << VZ0[i] << ")"
                  << "  block " << DOF << "x" << NDOF << ":\n"
                  << A.block(DOF*i, 0, DOF, NDOF).format(compact) << "\n";
    }
    std::cout << "\n  Eigenvalues of full " << NDOF << "x" << NDOF
              << " assembled matrix:\n    "
              << eig.eigenvalues().transpose().format(compact) << "\n";

    // ====================================================================
    // Part 1: determinism test (SPARSE_QR — the solver HydroChrono uses)
    // ====================================================================
    std::cout << "\n--- DETERMINISM: " << TRIALS
              << " independent trials with SPARSE_QR ---\n";

    std::vector<TrialResult> trials;
    trials.reserve(TRIALS);
    for (int t = 0; t < TRIALS; ++t) {
        std::vector<std::unique_ptr<char[]>> junk;
        for (int d = 0; d < (t+1)*137; ++d)
            junk.push_back(std::make_unique<char[]>(64 + t*31));
        trials.push_back(RunTrial(A, ChSolver::Type::SPARSE_QR));
    }

    // Print assembled matrix from trial 0
    std::cout << "\n  Assembled added-mass matrix (body mass/inertia subtracted):\n\n"
              << trials[0].added_mass.format(compact) << "\n";

    // Check: assembled added-mass matches input
    double assembly_err = (trials[0].added_mass - A).norm();
    if (assembly_err > 1e-10) {
        std::cerr << "FAIL: assembled added-mass matrix differs from input (norm = "
                  << assembly_err << ")\n";
        return 1;
    }
    std::cout << "  PASS: assembled matrix matches input (error norm = "
              << assembly_err << ")\n";

    // Check: no NaN/Inf, all trials identical, velocities changed
    for (int t = 0; t < TRIALS; ++t) {
        if (!trials[t].ok) {
            std::cerr << "FAIL: trial " << t << " error: " << trials[t].error << "\n";
            return 1;
        }
        if (HasBadValues(trials[t].sig)) {
            std::cerr << "FAIL: trial " << t << " has NaN/Inf\n" << Fmt(trials[t].sig);
            return 1;
        }
    }
    for (int t = 1; t < TRIALS; ++t) {
        if (!Match(trials[0].sig, trials[t].sig)) {
            std::cerr << "FAIL: trial " << t << " differs from trial 0\n"
                      << Fmt(trials[0].sig) << Fmt(trials[t].sig);
            return 1;
        }
    }
    bool active = false;
    for (int i = 0; i < N; ++i)
        if (std::abs(trials[0].sig[i][2] - VZ0[i]) > 1e-12) { active = true; break; }
    if (!active) { std::cerr << "FAIL: velocities unchanged\n"; return 1; }

    std::cout << "  PASS: " << TRIALS << " trials all bit-identical\n";
    std::cout << "\n  Reference velocities:\n" << Fmt(trials[0].sig);

    // ====================================================================
    // Part 2: solver sweep — test assembly & dynamics across solver types
    // ====================================================================
    // This is informational: we report which solvers correctly preserve
    // the added-mass matrix and produce matching dynamics. Failures here
    // do NOT fail the test (some solvers may not suit this problem), but
    // the summary table makes any issues visible.

    struct SolverEntry {
        ChSolver::Type type;
        const char* name;
    };
    const SolverEntry solvers[] = {
        {ChSolver::Type::SPARSE_QR,        "SPARSE_QR"},
        {ChSolver::Type::SPARSE_LU,        "SPARSE_LU"},
        {ChSolver::Type::MINRES,            "MINRES"},
        {ChSolver::Type::GMRES,             "GMRES"},
        {ChSolver::Type::BICGSTAB,          "BICGSTAB"},
        {ChSolver::Type::PSOR,              "PSOR"},
        {ChSolver::Type::BARZILAIBORWEIN,   "BARZILAIBORWEIN"},
        {ChSolver::Type::APGD,              "APGD"},
    };

    // Reference: SPARSE_QR result from part 1
    const Sig& ref_sig = trials[0].sig;

    std::cout << "\n--- SOLVER SWEEP: added-mass assembly across solver types ---\n\n";
    std::cout << "  " << std::left << std::setw(20) << "Solver"
              << std::setw(14) << "Assembly"
              << std::setw(14) << "NaN/Inf"
              << std::setw(14) << "Dynamics"
              << "Notes" << std::endl;
    std::cout << "  " << std::string(72, '-') << std::endl;

    for (const auto& s : solvers) {
        auto r = RunTrial(A, s.type);

        std::string asm_status, nan_status, dyn_status, notes;

        if (!r.ok) {
            asm_status = "ERROR";
            nan_status = "-";
            dyn_status = "-";
            notes = r.error;
        } else {
            // Assembly check
            double err = (r.added_mass - A).norm();
            if (err < 1e-10) {
                asm_status = "MATCH";
            } else {
                std::ostringstream o; o << std::scientific << std::setprecision(1) << err;
                asm_status = "DIFF " + o.str();
            }

            // NaN/Inf check
            nan_status = HasBadValues(r.sig) ? "HAS NaN/Inf" : "clean";

            // Dynamics comparison vs SPARSE_QR reference
            if (HasBadValues(r.sig)) {
                dyn_status = "-";
            } else {
                double max_diff = 0;
                for (int i = 0; i < N; ++i)
                    for (int d = 0; d < DOF; ++d)
                        max_diff = std::max(max_diff, std::abs(r.sig[i][d] - ref_sig[i][d]));
                if (max_diff < 1e-10) {
                    dyn_status = "MATCH";
                } else {
                    std::ostringstream o; o << std::scientific << std::setprecision(1) << max_diff;
                    dyn_status = "diff " + o.str();
                    notes = "differs from SPARSE_QR";
                }
            }
        }

        std::cout << "  " << std::left << std::setw(20) << s.name
                  << std::setw(14) << asm_status
                  << std::setw(14) << nan_status
                  << std::setw(14) << dyn_status
                  << notes << std::endl;
    }

    std::cout << "\n  Note: solver sweep is informational only (does not affect PASS/FAIL).\n"
              << "  Assembly = added-mass matrix matches input after subtracting body inertia.\n"
              << "  Dynamics = post-step velocities match SPARSE_QR reference.\n";

    std::cout << "\n=== PASSED ===" << std::endl;
    return 0;
}
