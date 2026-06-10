// Order-resolved singularities pipeline (2026-06-07): CLI vs C ABI
// snapshot + collector pin for the `emit_sings` extension of
// find_lr_orders.
//
// `emit_sings` is the opt-in per-face kinematic-divisor collector
// (notes/ristretto/order_resolved_singularities.md, Implementation
// step 2): when set, the find_lr_orders response additionally carries
//   "sings"       — the canonical irreducible KINEMATIC divisors (free
//                   of all integration variables, ANY degree) collected
//                   across the whole LR-search DP, deduped by the
//                   engine's proportionality canonical form;
//   "sings_total" — their count.
// No chi-vetting in C++ (the Mathematica side owns vetting); the
// collector observes factors BEFORE/ASIDE from the deg-2 letter cap that
// bounds the LR VERDICT, and never alters control flow.
//
// Cases:
//   A* — ABI sanity (non-NULL, single envelope) on the emit_sings req
//   N1 — byte-identity WITHOUT the option, CLI vs ABI (modulo timing).
//        This is the local transport-neutral form of gate #1; the
//        old-vs-new byte-identity gate is run out-of-band in the build
//        protocol (a HEAD baseline binary cannot be linked here).
//   C1 — byte-identity WITH emit_sings, CLI vs ABI (modulo timing):
//        the collector is transport-neutral.
//   P* — collector pin: on the massive-bubble fixture, the sings array
//        is EXACTLY {"mm","s","s - 4*mm"} (sings_total = 3), in
//        first-encounter order.  Each entry is hand-verified in the
//        report as a genuine irreducible kinematic divisor of the walk
//        (lc(F,x) = mm; disc(F,x) = x2^2 * s * (s - 4 mm)).
//   G1 — the no-option response carries NO "sings" / "sings_total"
//        field (the field is strictly opt-in).
//
// Mirrors test_c_abi_find_lr_orders_scan_cli_snapshot.cpp; the
// strip_timing / run_cli helpers are intentionally duplicated, per the
// "snapshot tests are independent" lineage expectation.
//
// Usage:
//   test_c_abi_find_lr_orders_emit_sings_cli_snapshot <path-to-hyperflint-cli>

#include "hyperflint/c_abi.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>  // close, unlink, mkstemp

namespace {

std::string strip_timing(const std::string& s) {
    const std::string marker = ",\"timing_compute_s\":";
    size_t a = s.find(marker);
    if (a == std::string::npos) return s;
    size_t b = a + marker.size();
    auto is_num_char = [](char c) {
        return (c >= '0' && c <= '9') || c == '.' || c == '-' || c == '+'
            || c == 'e' || c == 'E'
            || c == 'n' || c == 'N' || c == 'a' || c == 'A'
            || c == 'i' || c == 'I' || c == 'f' || c == 'F'
            || c == '(' || c == ')' || c == '_';
    };
    while (b < s.size() && is_num_char(s[b])) ++b;
    return s.substr(0, a) + s.substr(b);
}

std::string run_cli(const std::string& cli_path,
                    const std::string& request) {
    char req_path[] = "/tmp/hf_cli_emitsings_snapshot_req_XXXXXX";
    char out_path[] = "/tmp/hf_cli_emitsings_snapshot_out_XXXXXX";
    int rfd = ::mkstemp(req_path);
    int ofd = ::mkstemp(out_path);
    if (rfd < 0 || ofd < 0) {
        std::cerr << "[FAIL] mkstemp failed\n";
        return "";
    }
    ::close(rfd);
    ::close(ofd);
    {
        std::ofstream f(req_path);
        f << request;
    }
    std::string cmd = cli_path + " eval-json < " + req_path
                    + " > " + out_path + " 2>/dev/null";
    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        std::cerr << "[WARN] CLI returned exit code " << rc << "\n";
    }
    std::ifstream f(out_path);
    std::stringstream ss;
    ss << f.rdbuf();
    ::unlink(req_path);
    ::unlink(out_path);
    std::string out = ss.str();
    if (!out.empty() && out.back() == '\n') out.pop_back();
    return out;
}

size_t count_substr(const std::string& hay, const std::string& need) {
    size_t n = 0, pos = 0;
    while ((pos = hay.find(need, pos)) != std::string::npos) {
        ++n;
        pos += need.size();
    }
    return n;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <path-to-hyperflint-cli>\n";
        return 2;
    }
    const std::string cli_path = argv[1];

    int n_pass = 0;
    int n_fail = 0;
    auto check = [&](bool cond, const char* tag) {
        if (cond) { std::printf("[PASS] %s\n", tag); ++n_pass; }
        else      { std::printf("[FAIL] %s\n", tag); ++n_fail; }
    };

    // Massive-bubble fixture: U = x1 + x2, F = mm*(x1+x2)^2 - s*x1*x2,
    // integration vars {x1, x2}, kinematics {s, mm}.  The LR walk
    // produces lc(F,x) = mm, disc(F,x) = x2^2 * s * (s - 4 mm), and
    // res(U,F,x) ~ s*x2^2 — so the kinematic divisors are exactly
    // {mm, s, s - 4*mm}, the last being the threshold class.
    const std::string req_base =
        "\"polys\":[\"x1 + x2\",\"mm*(x1 + x2)^2 - s*x1*x2\"],"
        "\"xvars\":[\"x1\",\"x2\"],\"coeff_vars\":[\"s\",\"mm\"]";
    const std::string req_noopt =
        "{\"op\":\"find_lr_orders\"," + req_base + "}";
    const std::string req_emit =
        "{\"op\":\"find_lr_orders\"," + req_base + ",\"emit_sings\":true}";

    // ---- N1: byte-identity WITHOUT the option (CLI vs ABI) ----
    char* noopt_abi_raw = hf_find_lr_orders(req_noopt.c_str());
    std::string noopt_abi = noopt_abi_raw ? std::string(noopt_abi_raw)
                                          : std::string();
    hf_free_string(noopt_abi_raw);
    std::string noopt_cli = run_cli(cli_path, req_noopt);
    check(strip_timing(noopt_abi) == strip_timing(noopt_cli),
          "N1 no-option byte-identical CLI vs ABI (timing stripped)");

    // ---- G1: no-option response carries NO sings field ----
    check(noopt_abi.find("\"sings\":") == std::string::npos,
          "G1a no-option: no \"sings\" field");
    check(noopt_abi.find("\"sings_total\":") == std::string::npos,
          "G1b no-option: no \"sings_total\" field");

    // ---- A*/C1: emit_sings ABI sanity + CLI-vs-ABI byte identity ----
    char* emit_abi_raw = hf_find_lr_orders(req_emit.c_str());
    check(emit_abi_raw != nullptr, "A1 hf_find_lr_orders(emit): non-NULL");
    std::string emit_abi = emit_abi_raw ? std::string(emit_abi_raw)
                                        : std::string();
    hf_free_string(emit_abi_raw);
    check(count_substr(emit_abi, "\"schema_version\":") == 1,
          "A2 ABI: schema_version exactly once (no double-wrap)");
    check(count_substr(emit_abi, "\"hf_version\":") == 1,
          "A3 ABI: hf_version exactly once (no double-wrap)");

    std::string emit_cli = run_cli(cli_path, req_emit);
    check(!emit_cli.empty(), "B1 CLI(emit): non-empty stdout");
    check(emit_cli.rfind("{\"op\":\"find_lr_orders\"", 0) == 0,
          "B2 CLI(emit): starts with {\"op\":\"find_lr_orders\"");
    bool emit_equal = (strip_timing(emit_abi) == strip_timing(emit_cli));
    check(emit_equal,
          "C1 emit_sings byte-identical CLI vs ABI (timing stripped)");
    if (!emit_equal) {
        std::cerr << "    CLI stripped: " << strip_timing(emit_cli) << "\n";
        std::cerr << "    ABI stripped: " << strip_timing(emit_abi) << "\n";
    }

    // ---- P*: collector pin ----
    // Exact serialized fragment for the bubble fixture.  Each entry is a
    // genuine irreducible kinematic divisor of the walk (see header):
    //   mm        = lc(F, x)               (leading coefficient)
    //   s         = kinematic factor of disc(F, x) and res(U, F, x)
    //   s - 4*mm  = the bubble threshold (the other disc factor)
    // First-encounter order in the DP is deterministic for fixed input
    // + memo configuration (verified stable across repeated runs at
    // implementation time).  Canonical form: proportionality
    // representative => leading sign normalized (s - 4*mm, not
    // 4*mm - s nor -(s - 4*mm)).
    const std::string pin =
        "\"sings\":[\"mm\",\"s\",\"s - 4*mm\"],\"sings_total\":3";
    check(emit_abi.find(pin) != std::string::npos,
          "P1 ABI: sings == {mm, s, s - 4*mm}, total 3 (exact pin)");
    check(emit_cli.find(pin) != std::string::npos,
          "P2 CLI: sings == {mm, s, s - 4*mm}, total 3 (exact pin)");

    std::printf("[SUMMARY] PASS=%d FAIL=%d\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
